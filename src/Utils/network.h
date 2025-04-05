#pragma once

#include <map>
#include <unordered_map>
#include <vector>
#include<queue>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <memory>
#include <string>
#include <limits>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <ctime>
#include "utility_v.h"
#include "common.h"
#include "codon.h"

using namespace std;


namespace EnsembleDesign {

// Lattice class
class Lattice {
public:
    unordered_map<IndexType, vector<NodeType>> nodes;
    unordered_map<NodeType, vector<NodeNucWType>, hash_pair> left_edges;
    unordered_map<NodeType, vector<NodeNucWType>, hash_pair> right_edges;

    IndexType node_count; // Counter for nodes
    IndexType edge_count; // Counter for edges

    Lattice(): nodes(), left_edges(), right_edges(), node_count(0), edge_count(0) {};

    void add_edge(NodeType n1, NodeType n2, NucType nuc, ScoreType weight = 0.0f){
        right_edges[n1].push_back(make_tuple(n2, nuc, weight));
        left_edges[n2].push_back(make_tuple(n1, nuc, weight));
        edge_count++;
    }

    void add_node(NodeType n1){
        IndexType pos = get<0>(n1);
        nodes[pos].push_back(n1);
        node_count++;
    }

    int get_node_count() const {
        return node_count;
    }

    int get_edge_count() const {
        return edge_count;
    }
};

class WDFA {
public:
    vector<vector<NodeType>> nodes;
    vector<EdgeType> edges;
    vector<Parameter> parameters;

    unordered_map<NodeType, vector<EdgeType*>, hash_pair> left_edges;
    unordered_map<NodeType, vector<EdgeType*>, hash_pair> right_edges;
    unordered_map<NodeType, unordered_map<NodeType, vector<EdgeType*>, hash_pair>, hash_pair> auxiliary_left_edges;
    unordered_map<NodeType, unordered_map<NodeType, vector<EdgeType*>, hash_pair>, hash_pair> auxiliary_right_edges;

    WDFA() : nodes(), edges(), left_edges(), right_edges(), auxiliary_left_edges(), auxiliary_right_edges() {};

    void reserve_capacity(size_t node_count, size_t edge_count) {
        nodes.resize(node_count);
        edges.reserve(edge_count);
        parameters.reserve(edge_count);
    }

    void add_node(const NodeType& node){
        IndexType pos = get<0>(node);
        assert(nodes.size() > pos);
        nodes[pos].push_back(node);
    }

    void add_edge(const NodeType& n1, const NodeType& n2, const NucType& nuc, ScoreType weight, ScoreType grad, ScoreType cai){
        parameters.emplace_back(weight, grad, cai);
        edges.emplace_back(n1, n2, nuc, &parameters.back());

        assert(edges.capacity() > edges.size());
        EdgeType* edge_ptr = &edges.back();
        right_edges[n1].push_back(edge_ptr);
        left_edges[n2].push_back(edge_ptr);
        auxiliary_right_edges[n1][n2].push_back(edge_ptr);
        auxiliary_left_edges[n2][n1].push_back(edge_ptr);
    }


    EdgeType* findRightEdgePtr(const NodeType& node, const NucType& nuc) {
        for (const auto &edge_ptr : right_edges[node]) {
            if (std::get<2>(*edge_ptr) == nuc) {
                return edge_ptr;
            }
        }
        return nullptr;
    }

    EdgeType* findLeftEdgePtr(const NodeType& node, const NucType& nuc) {
        for (const auto &edge_ptr : left_edges[node]) {
            if (std::get<2>(*edge_ptr) == nuc) {
                return edge_ptr;
            }
        }
        return nullptr;
    }

    vector<EdgeType*> findPath(const string& nucs_str, const NodeType& start, const NodeType& end) {
        std::vector<EdgeType*> path;
        NodeType current_node = start;

        for (size_t i = 0; i < nucs_str.length(); ++i) {
            NucType nuc = GET_ACGU_NUC(nucs_str[i]);
            EdgeType* edge_ptr = findRightEdgePtr(current_node, nuc);
            if (edge_ptr == nullptr) {
                // No valid edge found, path cannot be completed
                break;
            }

            // Add the edge to the path
            path.push_back(edge_ptr);


            // Check if the next node is the end node and we are at the last nucleotide
            if (i == nucs_str.length() - 1 && current_node == end) {
                return path;
            }

            // Move to the next node
            current_node = std::get<1>(*edge_ptr);
        }

        // If the correct path is not found
        return vector<EdgeType*>();
    }

    void set_all_weights_to_neg_inf() {
        for (auto& param : parameters) {
            param.weight = util::value_min<ScoreType>();
        }
    }

    void apply_momentum(ScoreType t) {
        for (auto& param : parameters) {
            auto weight_ = param.weight;
            param.weight = std::log(std::max((1 + t) * Fast_Exp(param.weight) - t * Fast_Exp(param.prev_weight), 1e-10));
            param.prev_weight = weight_;
        }
    }

    void set_path_probability(const string& nucs, const NodeType& start_node, ScoreType epsilon) {
        set_all_weights_to_neg_inf();

        NodeType current_node = start_node;

        for (char nuc_char : nucs) {
            NucType nuc = GET_ACGU_NUC(nuc_char);

            auto edge_ptr = findRightEdgePtr(current_node, nuc);
            if (edge_ptr == nullptr) {
                throw std::runtime_error("Edge not found for nucleotide: " + std::string(1, nuc_char));
            }

            // Set the probability of the found edge
            std::get<3>(*edge_ptr)->weight = std::log(1 - epsilon);

            // Update probabilities for other edges
            size_t num_other_edges = right_edges[current_node].size() - 1;
            for (const auto &other_edge_ptr : right_edges[current_node]) {
                if (other_edge_ptr != edge_ptr) {
                    if (epsilon < 1e-20) {
                        std::get<3>(*other_edge_ptr)->weight = std::log(1e-20); // log(0) = -inf
                    } else {
                        std::get<3>(*other_edge_ptr)->weight = std::log(epsilon / num_other_edges);
                    }
                }
            }

            // Move to the next node
            current_node = std::get<1>(*edge_ptr);
        }
    }

    void random_init_with_mfe_path(const string& nucs, const NodeType& start_node, ScoreType epsilon, unsigned int seed) {
        randomize_weights(seed);
        validate_and_project_probabilities();

        std::mt19937 rng(seed + 2333); // Random number generator
        NodeType current_node = start_node;
        for (char nuc_char : nucs) {
            NucType nuc = GET_ACGU_NUC(nuc_char);


            auto edge_ptr = findRightEdgePtr(current_node, nuc);
            if (edge_ptr == nullptr) {
                throw std::runtime_error("Edge not found for nucleotide: " + std::string(1, nuc_char));
            }

            // Get the other edges
            std::vector<EdgeType*> other_edges;
            for (auto& e_ptr : right_edges[current_node]) {
                if (e_ptr != edge_ptr) {
                    other_edges.push_back(e_ptr);
                }
            }

            // Special case: only one right edge
            if (other_edges.empty()) {
                std::get<3>(*edge_ptr)->weight = 0; // log(1) = 0
            }
            else {

                // Set the probability of the found edge
                std::get<3>(*edge_ptr)->weight = std::log(1 - epsilon);

                // Randomly distribute epsilon among the other edges
                std::vector<ScoreType> random_fractions(other_edges.size(), 0.0);
                std::uniform_real_distribution<ScoreType> dist(0.0, 1.0);
                ScoreType total_fraction = 0.0;

                for (auto &fraction: random_fractions) {
                    fraction = dist(rng); // Assign random fraction
                    total_fraction += fraction;
                }

                // Normalize and assign weights
                for (size_t i = 0; i < other_edges.size(); ++i) {
                    random_fractions[i] /= total_fraction; // Normalize
                    random_fractions[i] *= epsilon; // Scale to epsilon
                    auto w = random_fractions[i];
                    if (w < 1e-10) w = 1e-10;
                    std::get<3>(*other_edges[i])->weight = std::log(w);
                }
            }


            // Move to the next node
            current_node = std::get<1>(*edge_ptr);
        }
    }

    void randomize_weights(unsigned int seed) {
        std::mt19937 rng(seed); // Random number generator
        std::uniform_real_distribution<ScoreType> dist(0.0, 1.0); // Distribution in the range (0, 1]

        for (auto& param : parameters) {
            param.weight = log(dist(rng)); // Assign a random weight
        }
    }

    void zero_gradients(){
        for(auto& param: parameters){
            param.zero_gradient();
        }
    }
    void gradient_decent(double learning_rate, ScoreType Z){
        for(auto& param : parameters){
            Fast_LogPlusEquals(param.weight, log(learning_rate) + param.gradient - Z);
        }
    }
    void validate_and_project_probabilities(){
        for(auto& i_nodes : nodes) {
            for (auto &node: i_nodes) {
                validate_and_project_probabilities_onenode_(node);
            }
        }
    }

    pair<string, ScoreType> get_best_nuc_sequence(int seq_length) {
        string best_sequence;
        ScoreType best_logprob = util::value_zero<ScoreType>();
        for (int i = 0; i < seq_length; i += 3) {
            if (i + 3 > seq_length) break;
            auto best_codon = find_best_codon(i);
            best_sequence += best_codon.first;
            best_logprob += best_codon.second;

        }
        return make_pair(best_sequence, best_logprob);
    }

private:

    pair<string, ScoreType> find_best_codon(int start) {
        NodeType start_node(start, 0);
        NodeType end_node(start + 3, 0);

        struct Path {
            NodeType node;
            std::string nuc_sequence;
            ScoreType log_prob;
        };

        std::queue<Path> q;
        q.push({start_node, "", 0.0}); // Initialize queue with start node

        ScoreType max_logprob = util::value_min<ScoreType>();
        string max_prob_sequence;

        while (!q.empty()) {
            Path current = q.front();
            q.pop();

            if (current.node == end_node) {
                if (current.log_prob > max_logprob) {
                    max_logprob = current.log_prob;
                    max_prob_sequence = current.nuc_sequence;
                }
                continue;
            }

            auto it = right_edges.find(current.node);
            if (it != right_edges.end()) {
                for (const auto& edge_ptr : it->second) {
                    NodeType next_node = std::get<1>(*edge_ptr);
                    NucType nuc = std::get<2>(*edge_ptr);
                    ScoreType edge_weight = std::get<3>(*edge_ptr)->weight;

                    // Assuming GET_ACGU converts NucType to a nucleotide character
                    char nuc_char = GET_ACGU(nuc);

                    Path next_path = {next_node, current.nuc_sequence + nuc_char, current.log_prob + edge_weight};
                    q.push(next_path);
                }
            }
        }


        return make_pair(max_prob_sequence, max_logprob);
    }

    void validate_and_project_probabilities_onenode_(const NodeType& node) {
        auto it = right_edges.find(node);
        if (it == right_edges.end()) {
            return; // No right edges for this node
        }

        // Convert log probabilities to actual probabilities
        std::vector<double> probabilities;
        double sum = 0.0;
        for (auto& edge_ptr : it->second) {
            double prob = std::exp(std::get<3>(*edge_ptr)->weight);
            probabilities.push_back(prob);
            sum += prob;
        }

        // Check if sum is close to 1
        const double EPSILON = 1e-6;
        if (std::abs(sum - 1.0) > EPSILON) {
            // Project probabilities back to the simplex
            auto v_proj = ProjectToSimplex(probabilities);

            // Update the weights with projected probabilities
            for (size_t i = 0; i < v_proj.size(); ++i) {
                std::get<3>(*it->second[i])->weight = std::log(v_proj[i]);
            }
        }
    }

};


unordered_map<string, Lattice> read_wheel(string const &filename) {
    unordered_map<string, Lattice> aa_graphs;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) {
        printf("Unable to open coding_wheel file\n");
        exit(1);   // call system to stop
    }

    vector<string> stuff;
    vector<string> option_splited;
    string aa;
    IndexType i;

    for (string line; getline(inFile, line);) {
        stuff = util::split(line, '\t');

        aa = stuff[0];
        Lattice graph = Lattice();
        graph.add_node(make_pair(0,0)); // always func1 with node (0,0)

        char last_first = 0;
        vector<string>::iterator iter = stuff.begin();
        ++iter; // position 0 is aa name
        i = 0;
        while(iter != stuff.end()){
            string option = *iter;
            option_splited = util::split(option, ' ');
            char first = option_splited[0][0];
            char second = option_splited[1][0];
            string thirds = option_splited[2];
            NodeType n2 = make_pair(2, i);
            graph.add_node(n2);
            NodeType n1;
            if (first != last_first) {
                n1 = make_pair(1, i);
                graph.add_node(n1);
                graph.add_edge(make_pair(0, 0), n1, GET_ACGU_NUC(first));
            }
            else {
                n1 = make_pair(1, i-1);
            }
            last_first = first;
            graph.add_edge(n1, n2, GET_ACGU_NUC(second));
            for (auto& third : thirds) {
                graph.add_edge(n2, make_pair(0,0), GET_ACGU_NUC(third));
            }
            i++; iter++;
        }
        aa_graphs[aa] = graph;
    }
    inFile.close();
    return aa_graphs;
}


unordered_map<string, Lattice> read_wheel_with_weights(const std::string& filename,
        std::unordered_map<std::string, std::unordered_map<NodeType, double, hash_pair>>& nodes_with_best_weight,
        std::unordered_map<std::string, std::unordered_map<NodeNucNodeType, double, std::hash<NodeNucNodeType>>>& edges_with_best_weight,
        const Codon& codon) {
    unordered_map<string, Lattice> aa_graphs;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) 
        throw std::runtime_error("Unable to open coding_wheel file\n");

    vector<string> stuff;
    vector<string> option_splited;
    string aa;
    IndexType i;

    for (string line; getline(inFile, line);) {
        stuff = util::split(line, '\t');
        aa = stuff[0];
        Lattice graph = Lattice();
        graph.add_node(make_pair(0,0)); // always func1 with node (0,0)

        char last_first = 0;
        vector<string>::iterator iter = stuff.begin();
        ++iter; // position 0 is aa name
        i = 0;
        while(iter != stuff.end()){
            string option = *iter;
            option_splited = util::split(option, ' ');
            char first = option_splited[0][0];
            char second = option_splited[1][0];
            string thirds = option_splited[2];
            NodeType n2 = make_pair(2, i);
            graph.add_node(n2);
            NodeType n1;
            if (first != last_first) {
                n1 = make_pair(1, i);
                graph.add_node(n1);
                auto first_num = GET_ACGU_NUC(first);

                double weight = 0.0f;
                if (nodes_with_best_weight[aa].count(make_pair(0, 0))) {
                    weight = edges_with_best_weight[aa][make_tuple(make_pair(0, 0), first_num, n1)] / nodes_with_best_weight[aa][make_pair(0, 0)];
                }

                graph.add_edge(make_pair(0, 0), n1, first_num, weight);
            }
            else {
                n1 = make_pair(1, i-1);
            }
            
            last_first = first;

            auto second_num = GET_ACGU_NUC(second);

            double weight = 0.0f;
            if (nodes_with_best_weight[aa].count(n1)) {
                weight = edges_with_best_weight[aa][make_tuple(n1, second_num, n2)] / nodes_with_best_weight[aa][n1];
            }

            graph.add_edge(n1, n2, second_num, weight);

            for (auto& third : thirds) {

                std::string three_nums = std::string(1, first) + std::string(1, second) + std::string(1, third);

                double weight = 0.0f;
                if (nodes_with_best_weight[aa].count(n2)) {
                    weight = codon.get_weight(aa, three_nums) / nodes_with_best_weight[aa][n2];
                } else {
                    weight = codon.get_weight(aa, three_nums);
                }

                graph.add_edge(n2, make_pair(0,0), GET_ACGU_NUC(third), weight);
            }
            i++; iter++;
        }
        aa_graphs[aa] = graph;
    }

    inFile.close();
    return aa_graphs;
}


unordered_map<string, Lattice> read_wheel_with_weights_log(const std::string& filename,
        std::unordered_map<std::string, std::unordered_map<NodeType, double, hash_pair>>& nodes_with_best_weight,
        std::unordered_map<std::string, std::unordered_map<NodeNucNodeType, double, std::hash<NodeNucNodeType>>>& edges_with_best_weight,
        const Codon& codon, double lambda_) {
    unordered_map<string, Lattice> aa_graphs;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) 
        throw std::runtime_error("Unable to open coding_wheel file\n");

    vector<string> stuff;
    vector<string> option_splited;
    string aa;
    IndexType i;

    for (string line; getline(inFile, line);) {
        stuff = util::split(line, '\t');
        aa = stuff[0];
        Lattice graph = Lattice();
        graph.add_node(make_pair(0,0)); // always func1 with node (0,0)

        char last_first = 0;
        vector<string>::iterator iter = stuff.begin();
        ++iter; // position 0 is aa name
        i = 0;
        while(iter != stuff.end()){
            string option = *iter;
            option_splited = util::split(option, ' ');
            char first = option_splited[0][0];
            char second = option_splited[1][0];
            string thirds = option_splited[2];
            NodeType n2 = make_pair(2, i);
            graph.add_node(n2);
            NodeType n1;
            if (first != last_first) {
                n1 = make_pair(1, i);
                graph.add_node(n1);
                auto first_num = GET_ACGU_NUC(first);

                double weight = 1.0f;
                if (nodes_with_best_weight[aa].count(make_pair(0, 0))) {
                    weight = lambda_ * log(edges_with_best_weight[aa][make_tuple(make_pair(0, 0), first_num, n1)] / nodes_with_best_weight[aa][make_pair(0, 0)]);
                }

                graph.add_edge(make_pair(0, 0), n1, first_num, weight);
            }
            else {
                n1 = make_pair(1, i-1);
            }
            
            last_first = first;

            auto second_num = GET_ACGU_NUC(second);

            double weight = 1.0f;
            if (nodes_with_best_weight[aa].count(n1)) {
                weight = lambda_ * log(edges_with_best_weight[aa][make_tuple(n1, second_num, n2)] / nodes_with_best_weight[aa][n1]);
            }

            graph.add_edge(n1, n2, second_num, weight);

            for (auto& third : thirds) {

                std::string three_nums = std::string(1, first) + std::string(1, second) + std::string(1, third);

                double weight = 1.0f;
                if (nodes_with_best_weight[aa].count(n2)) {
                    weight = lambda_ *  log(codon.get_weight(aa, three_nums) / nodes_with_best_weight[aa][n2]);
                } else {
                    weight = lambda_ *  log(codon.get_weight(aa, three_nums));
                }

                graph.add_edge(n2, make_pair(0,0), GET_ACGU_NUC(third), weight);
            }
            i++; iter++;
        }
        aa_graphs[aa] = graph;
    }

    inFile.close();
    return aa_graphs;
}

void prepare_codon_unit_lattice(const std::string& wheel_path, const Codon& codon,
        std::unordered_map<string, Lattice>& aa_graphs_with_ln_weights_ret,
        std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<double, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>>&
                best_path_in_one_codon_unit_ret,
        std::unordered_map<std::string, std::string>& aa_best_path_in_a_whole_codon_ret, double lambda_) {

    std::unordered_map<std::string, std::unordered_map<NodeType, ScoreType, hash_pair>> nodes_with_best_weight;
    std::unordered_map<std::string, std::unordered_map<NodeNucNodeType, ScoreType, std::hash<NodeNucNodeType>>> edges_with_best_weight;

    unordered_map<string, Lattice> aa_graphs_with_ln_weights;
    unordered_map<string, Lattice> aa_graphs_with_weights = read_wheel_with_weights(wheel_path, nodes_with_best_weight, edges_with_best_weight, codon);

    for (auto& aa_aa_elem : aa_graphs_with_weights) {
        auto& aa = aa_aa_elem.first;
        auto& aa_elem = aa_aa_elem.second;
        for (auto& node_at_2 : aa_elem.nodes[2]) {
            for (auto& node_at_3_nuc_weight : aa_elem.right_edges[node_at_2]) {
                auto node_at_3 = std::get<0>(node_at_3_nuc_weight);
                auto nuc = std::get<1>(node_at_3_nuc_weight);
                auto weight = std::get<2>(node_at_3_nuc_weight);
                nodes_with_best_weight[aa][node_at_2] = max(nodes_with_best_weight[aa][node_at_2], weight);
                edges_with_best_weight[aa][make_tuple(node_at_2,nuc,node_at_3)] = weight;
            }
        }

        for (auto& node_at_1 : aa_elem.nodes[1]) {
            for (auto& node_at_2_nuc_weight : aa_elem.right_edges[node_at_1]) {
                auto node_at_2 = std::get<0>(node_at_2_nuc_weight);
                auto nuc = std::get<1>(node_at_2_nuc_weight);
                nodes_with_best_weight[aa][node_at_1] = max(nodes_with_best_weight[aa][node_at_1], nodes_with_best_weight[aa][node_at_2]);
                edges_with_best_weight[aa][make_tuple(node_at_1,nuc,node_at_2)] = nodes_with_best_weight[aa][node_at_2];
            }
        }

        for (auto& node_at_0 : aa_elem.nodes[0]) {
            for (auto& node_at_1_nuc_weight : aa_elem.right_edges[node_at_0]) {
                auto node_at_1 = std::get<0>(node_at_1_nuc_weight);
                auto nuc = std::get<1>(node_at_1_nuc_weight);
                nodes_with_best_weight[aa][node_at_0] = max(nodes_with_best_weight[aa][node_at_0], nodes_with_best_weight[aa][node_at_1]);
                edges_with_best_weight[aa][make_tuple(node_at_0,nuc,node_at_1)] = nodes_with_best_weight[aa][node_at_1];
            }
        }
    }

    aa_graphs_with_ln_weights = read_wheel_with_weights_log(wheel_path,  nodes_with_best_weight, edges_with_best_weight, codon, lambda_);

    std::unordered_map<std::string, 
                       std::unordered_map<std::tuple<NodeType, NodeType>, 
                       std::tuple<double, NucType, NucType>,
                       std::hash<std::tuple<NodeType, NodeType>>>>
                       best_path_in_one_codon_unit;


    for (auto& aa_graph : aa_graphs_with_ln_weights) {
        auto& aa = aa_graph.first;
        auto& graph = aa_graph.second;
        for (auto& node_0 : graph.nodes[0]) {
            for (auto& node_1_nuc_log_w : graph.right_edges[node_0]) {
                auto node_1 = std::get<0>(node_1_nuc_log_w);
                auto nuc = std::get<1>(node_1_nuc_log_w);
                auto log_weight = std::get<2>(node_1_nuc_log_w);

                if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_0,node_1)))
                    best_path_in_one_codon_unit[aa][make_tuple(node_0,node_1)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc);

                double current_log_weight = std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_0,node_1)]);
                if (current_log_weight < log_weight) {
                    best_path_in_one_codon_unit[aa][make_tuple(node_0,node_1)] = make_tuple(log_weight,nuc,k_void_nuc);
                }
            }
        }

        for (auto& node_1 : graph.nodes[1]) {
            for (auto& node_2_nuc_log_w : graph.right_edges[node_1]) {
                auto node_2 = std::get<0>(node_2_nuc_log_w);
                auto nuc = std::get<1>(node_2_nuc_log_w);
                auto log_weight = std::get<2>(node_2_nuc_log_w);

                if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_1,node_2)))
                    best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc);

                double current_log_weight = std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)]);
                if (current_log_weight < log_weight) {
                    best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)] = make_tuple(log_weight,nuc,k_void_nuc);
                }

                auto temp = best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)];
            }
        }

        for (auto& node_2 : graph.nodes[2]) {
            for (auto& node_3_nuc_log_w : graph.right_edges[node_2]) {
                auto node_3 = std::get<0>(node_3_nuc_log_w);
                auto nuc = std::get<1>(node_3_nuc_log_w);
                auto log_weight = std::get<2>(node_3_nuc_log_w);

                if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_2,node_3)))
                    best_path_in_one_codon_unit[aa][make_tuple(node_2,node_3)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc);

                double current_log_weight = std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_2,node_3)]);
                if (current_log_weight < log_weight) {
                    best_path_in_one_codon_unit[aa][make_tuple(node_2,node_3)] = make_tuple(log_weight,nuc,k_void_nuc);
                }
            }
        }

        for (auto& node_0 : graph.nodes[0]) {
            for (auto& node_1_nuc_0_log_weight_0 : graph.right_edges[node_0]) {
                auto& node_1 = std::get<0>(node_1_nuc_0_log_weight_0);
                auto& nuc_0 = std::get<1>(node_1_nuc_0_log_weight_0);
                auto log_weight_0 = std::get<2>(node_1_nuc_0_log_weight_0);
                for (auto& node_2_nuc_1_log_weight_1 : graph.right_edges[node_1]) {
                    auto& node_2 = std::get<0>(node_2_nuc_1_log_weight_1);
                    auto& nuc_1 = std::get<1>(node_2_nuc_1_log_weight_1);
                    auto log_weight_1 = std::get<2>(node_2_nuc_1_log_weight_1);

                    if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_0,node_2)))
                        best_path_in_one_codon_unit[aa][make_tuple(node_0,node_2)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc);

                    if (std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_0,node_2)]) < log_weight_0 + log_weight_1)
                        best_path_in_one_codon_unit[aa][make_tuple(node_0,node_2)] = make_tuple(log_weight_0 + log_weight_1, nuc_0, nuc_1);
                }
            }
        }

        for (auto& node_1 : graph.nodes[1]) {
            for (auto& node_2_nuc_1_log_weight_1 : graph.right_edges[node_1]) {
                auto& node_2 = std::get<0>(node_2_nuc_1_log_weight_1);
                auto& nuc_1 = std::get<1>(node_2_nuc_1_log_weight_1);
                auto log_weight_1 = std::get<2>(node_2_nuc_1_log_weight_1);
                for (auto& node_3_nuc_2_log_weight_2 : graph.right_edges[node_2]) {
                    auto& node_3 = std::get<0>(node_3_nuc_2_log_weight_2);
                    auto& nuc_2 = std::get<1>(node_3_nuc_2_log_weight_2);
                    auto log_weight_2 = std::get<2>(node_3_nuc_2_log_weight_2);

                    if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_1,node_3)))
                        best_path_in_one_codon_unit[aa][make_tuple(node_1,node_3)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc);

                    if (std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_1,node_3)]) < log_weight_1 + log_weight_2)
                        best_path_in_one_codon_unit[aa][make_tuple(node_1,node_3)] = make_tuple(log_weight_1 + log_weight_2, nuc_1, nuc_2);
                }
            }
        }
    }

    std::unordered_map<std::string, double> max_path;
    std::unordered_map<std::string, std::string> aa_best_path_in_a_whole_codon;

    for (auto& aa_path_weight : codon.aa_table_) {
        auto& aa = aa_path_weight.first; // char
        for (auto& path_weight : aa_path_weight.second) {
            if (max_path[aa] < path_weight.second) {
                max_path[aa] = path_weight.second;
                aa_best_path_in_a_whole_codon[aa] = path_weight.first;
            }
        }
    }

    aa_graphs_with_ln_weights_ret = aa_graphs_with_ln_weights;
    best_path_in_one_codon_unit_ret = best_path_in_one_codon_unit;
    aa_best_path_in_a_whole_codon_ret = aa_best_path_in_a_whole_codon;
}



WDFA get_dfa(unordered_map<string, Lattice> aa_graphs, vector<string> aa_seq) {
    IndexType total_nodes = 0;
    IndexType total_edges = 0;

    for(const auto& aa : aa_seq) {
        Lattice graph = aa_graphs[aa];
        total_nodes += graph.get_node_count();
        total_edges += graph.get_edge_count();
    }

    WDFA dfa = WDFA();
    dfa.reserve_capacity(total_nodes + 10, total_edges + 10);

    NodeType start_node = make_pair(3 * static_cast<IndexType>(aa_seq.size()), 0);
    dfa.add_node(start_node);

    IndexType i = 0;
    IndexType i3;

    for(auto& aa : aa_seq) {
        i3 = i * 3;
        Lattice graph = aa_graphs[aa];
        for (IndexType pos = 0; pos <= 2; pos++) {
            for(auto& node : graph.nodes[pos]) {
                NumType num1 = node.second;
                NodeType node1 = make_pair(i3 + pos, num1);
                dfa.add_node(node1);
                for (auto& edge : graph.right_edges[node]) {
                    NodeType n2 = get<0>(edge);
                    NucType nuc = get<1>(edge);
                    NumType num2 = get<1>(n2);
                    NodeType node2 = make_pair(i3 + pos + 1, num2);
                    auto weight = util::value_zero<ScoreType>();
                    auto grad = util::value_min<ScoreType>();
                    auto cai = get<2>(edge);
                    dfa.add_edge(node1, node2, nuc, weight, grad, cai);
                }
            }
        }
        i++;
    }

    return dfa;
}

}
