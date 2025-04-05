#pragma once

#include <memory>
#include <string>
#include <cstring>
#include <limits>
#include <vector>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <stdexcept>
#include <unordered_set>

#include "Utils/network.h"
#include "Utils/codon.h"
#include "Utils/common.h"
#include "Utils/utility_v.h"

using namespace std;

namespace EnsembleDesign {

string get_nuc_from_dfa_cai(WDFA& dfa, const NodeType& start_node, const NodeType& end_node,
        const std::vector<std::string>& protein, std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, 
        std::tuple<double, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>>&
        best_path_in_one_codon_unit, std::unordered_map<std::string, std::string>& aa_best_path_in_a_whole_codon) {

    IndexType s_index = start_node.first;
    IndexType t_index = end_node.first;

    if (s_index >= t_index)
        return "";

    auto aa_left = protein[s_index / 3]; // tri letter
    auto aa_right = protein[t_index / 3];
    auto start_node_re_index = make_pair(s_index % 3, start_node.second);
    auto end_node_re_index = make_pair(t_index % 3, end_node.second);
    if (t_index - s_index < 3) {
        if (s_index / 3 == t_index / 3) {
            std::string temp_seq = "";
            auto& nucs = best_path_in_one_codon_unit[aa_left][make_tuple(start_node_re_index, end_node_re_index)];
            temp_seq.append(1, GET_ACGU(std::get<1>(nucs)));
            if (std::get<2>(nucs) != k_void_nuc) 
                temp_seq.append(1, GET_ACGU(std::get<2>(nucs)));

            if (temp_seq.length() != end_node.first - start_node.first) {
                assert(false);
            }
            return temp_seq;
        } else {
            std::string temp_left = "";
            std::string temp_right = "";
            if (s_index % 3 != 0) {
                auto& nucs = best_path_in_one_codon_unit[aa_left][make_tuple(start_node_re_index, make_pair(0, 0))];
                temp_left.append(1, GET_ACGU(std::get<1>(nucs)));
                if (std::get<2>(nucs) != k_void_nuc) 
                    temp_left.append(1, GET_ACGU(std::get<2>(nucs)));
            }

            if (t_index % 3 != 0) {
                auto& nucs = best_path_in_one_codon_unit[aa_right][make_tuple(make_pair(0, 0), end_node_re_index)];
                temp_right.append(1, GET_ACGU(std::get<1>(nucs)));
                if (std::get<2>(nucs) != k_void_nuc) 
                    temp_right.append(1, GET_ACGU(std::get<2>(nucs)));
            }

            assert((temp_left + temp_right).length() == end_node.first - start_node.first);

            return temp_left + temp_right;
        }

    } else {

        std::string temp_left = "";
        std::string temp_mid = "";
        std::string temp_right = "";

        if (s_index % 3 != 0) {
            auto& nucs = best_path_in_one_codon_unit[aa_left][make_tuple(start_node_re_index, make_pair(0, 0))];
            temp_left.append(1, GET_ACGU(std::get<1>(nucs)));
            if (std::get<2>(nucs) != k_void_nuc) 
                temp_left.append(1, GET_ACGU(std::get<2>(nucs)));
        }

        IndexType protein_start_index = s_index / 3;
        if (s_index % 3 != 0)
            protein_start_index++;

        IndexType protein_end_index = t_index / 3;

        if (protein_start_index != protein_end_index) {
            for (IndexType protein_index = protein_start_index; protein_index < protein_end_index; ++protein_index) {
                
                std::string nucs;
                auto aa_tri = protein[protein_index];
                if (k_map_3_1.count(aa_tri)) {
                    nucs = aa_best_path_in_a_whole_codon[std::string(1, k_map_3_1[aa_tri])];
                } else if (aa_best_path_in_a_whole_codon.count(aa_tri)) {
                    nucs = aa_best_path_in_a_whole_codon[aa_tri];
                } else {
                    assert(false);
                }

                for (auto nuc : nucs) {
                    temp_mid.append(1, nuc);
                }
            }
        }

        if (t_index % 3 != 0) {
            auto& nucs = best_path_in_one_codon_unit[aa_right][make_tuple(make_pair(0, 0), end_node_re_index)];
            temp_right.append(1, GET_ACGU(std::get<1>(nucs)));
            if (std::get<2>(nucs) != k_void_nuc) 
                temp_right.append(1, GET_ACGU(std::get<2>(nucs)));
        }

        assert((temp_left + temp_mid + temp_right).length() == end_node.first - start_node.first);

        return temp_left + temp_mid + temp_right;
    }
}

class Optimizer {
public:

    int beam_size=0;
    int num_epochs;
    double learning_rate;
    double epsilon;
    bool enable_cube_pruning=false;
    unsigned int rand_seed=0;
    std::string init_solution;

    using State_t = State<ScoreType>;
    using DFA_t = WDFA;
    using ScoreInnerDate_t = ScoreInnerDate;
    using NextPair_t = vector<unordered_map<NodeType, vector<NextPairType>, hash_pair>>;
    using NextPairSet_t = vector<unordered_map<NodeType, set<NextPairType>, hash_pair>>;
    using BestC_t = unordered_map<NodeType, State_t, hash_pair>;
    using BestX_t = unordered_map<NodeType, unordered_map<NodeType, unordered_map<NucPairType, State_t>, hash_pair>, hash_pair>;
    using BestM_t = unordered_map<NodeType, unordered_map<NodeType, State_t, hash_pair>, hash_pair>;
    using sorted_BestM_t = unordered_map<NodeType, vector<ScoreInnerDate_t>, hash_pair>;
    using PrefixScore_t = unordered_map<NodeType, ScoreType, hash_pair>;

    Optimizer(int beam_size_, int num_epochs_, double learning_rate_, double epsilon_, std::string init_solution_, bool is_verbose_, unsigned int rand_seed_);

    void optimize(DFA_t& dfa,
        Codon& codon, 
        std::string& aa_seq, 
        std::vector<std::string>& p, 
        std::unordered_map<std::string, std::string>& aa_best_in_codon,
        std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>,
        std::tuple<double, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>>& best_path_in_one_codon,
        std::unordered_map<string, Lattice>& aa_graphs_with_ln_weights);
    void clear_parser_states();
    void MFE(DFA_t& dfa);
    void forward(DFA_t& dfa);
    void backward(DFA_t& dfa);
    ScoreType beam_test(DFA_t& dfa, int beam);

private:
    
    template <PhaseOption phase>
    void hairpin_beam(IndexType j, IndexType j_num, DFA_t& dfa);
    
    template <PhaseOption phase>
    void Multi_beam(IndexType j, IndexType j_num, DFA_t& dfa);
    
    template <PhaseOption phase>
    void P_beam(IndexType j, IndexType j_num, DFA_t& dfa);
    
    template <PhaseOption phase>
    void M2_beam(IndexType j, IndexType j_num, DFA_t& dfa);
    
    template <PhaseOption phase>
    void M_beam(IndexType j, IndexType j_num, DFA_t& dfa);
    
    template <PhaseOption phase>
    void C_beam(IndexType j, IndexType j_num, DFA_t& dfa);



    template <PhaseOption phase>
    void update_state(State_t &u, ScoreType w = util::value_zero<ScoreType>(), const Parameters &params = Parameters()){
        if (phase == PhaseOption::Inside){
            update_inside(u, w, params);
        }
        else if (phase == PhaseOption::Outside){
            update_outside(u, w, params);
        }
        else{
            throw std::logic_error("Function not yet implemented");
        }
    }
    void update_inside(State_t &u, const ScoreType w, const Parameters &params) {
        ScoreType logprob = util::value_zero<ScoreType>();
        for (const auto& param : params) {
            logprob += param->weight;
        }

        Fast_LogPlusEquals(u.inside, w + logprob);
    }
    void update_outside(State_t &u, const ScoreType w, const Parameters &params) {
        ScoreType logprob = util::value_zero<ScoreType>();
        for (const auto& param : params) {
            logprob += param->weight;
        }

        for (const auto& param : params) {
            param->update_gradient(w + logprob - param->weight + u.outside);
        }
    }


    template <PhaseOption phase>
    void update_state(State_t &u, State_t &v, ScoreType w = util::value_zero<ScoreType>(), const Parameters &params = Parameters(), NodeType pre_node = NodeType()) {
        if (phase == PhaseOption::Inside){
            update_inside(u, v, w, params, pre_node);
        }
        else if (phase == PhaseOption::Outside){
            update_outside(u, v, w, params);
        }
        else{
            throw std::logic_error("Function not yet implemented");
        }
    }
    void update_inside(State_t &u, const State_t &v, const ScoreType w, const Parameters &params, NodeType pre_node) {
        ScoreType logprob = util::value_zero<ScoreType>();
        for (const auto& param : params) {
            logprob += param->weight;
        }

        Fast_LogPlusEquals(u.inside, v.inside + w + logprob);
    }
    void update_outside(const State_t &u, State_t &v, const ScoreType w, const Parameters &params) {
        ScoreType logprob = util::value_zero<ScoreType>();
        for (const auto& param : params) {
            logprob += param->weight;
        }

        for (const auto& param : params) {
            param->update_gradient(v.inside + w + logprob - param->weight + u.outside);
        }

        Fast_LogPlusEquals(v.outside, u.outside + w + logprob);
    }


    template <PhaseOption phase>
    void update_state(State_t &u, State_t &v1, State_t &v2, ScoreType w = util::value_zero<ScoreType>(), const Parameters &params = Parameters()) {
        if (phase == PhaseOption::Inside){
            update_inside(u, v1, v2, w, params);
        }
        else if (phase == PhaseOption::Outside){
            update_outside(u, v1, v2, w, params);
        }
        else{
            throw std::logic_error("Function not yet implemented");
        }
    }
    void update_inside(State_t &u, const State_t &v1, const State_t &v2, const ScoreType w, const Parameters &params) {
        ScoreType logprob = util::value_zero<ScoreType>();
        for (const auto& param : params) {
            logprob += param->weight;
        }

        Fast_LogPlusEquals(u.inside, v1.inside + v2.inside + w + logprob);
    }
    void update_outside(const State_t &u, State_t &v1, State_t &v2, const ScoreType w, const Parameters &params) {
        ScoreType logprob = util::value_zero<ScoreType>();
        for (const auto& param : params) {
            logprob += param->weight;
        }

        for (const auto& param : params) {
            param->update_gradient(v1.inside + v2.inside + w + logprob - param->weight + u.outside);
        }

        Fast_LogPlusEquals(v1.outside, u.outside + v2.inside + w + logprob);
        Fast_LogPlusEquals(v2.outside, v1.inside + u.outside + w + logprob);
    }


    void get_next_pair(DFA_t& dfa);
    void get_next_pair_set();

    void get_prev_pair(DFA_t& dfa);
    void get_prev_pair_set();

    void preprocess(DFA_t& dfa);

    ScoreType quickselect_partition(std::vector<ScoreInnerDate_t>& scores,
        ScoreType lower, ScoreType upper);

    ScoreType quickselect(std::vector<ScoreInnerDate_t>& scores,
        const ScoreType lower, const ScoreType upper, const IndexType k);

    template <PhaseOption phase>
    void beam_prune(DFA_t& dfa, BestX_t& bestX, const IndexType j);

    template <PhaseOption phase>
    double beam_prune(DFA_t& dfa, BestM_t& bestX, const IndexType j);

    bool is_verbose;
    double EPS = 1e-8;

    IndexType seq_length; 

    BestX_t bestH, bestP, bestMulti;
    BestM_t bestM2, bestM;
    BestC_t bestC;

    NextPair_t next_pair;
    NextPairSet_t next_pair_set;

    NextPair_t prev_pair;
    NextPairSet_t prev_pair_set;

    NextPair_t next_list;
    NextPair_t prev_list;

    vector<vector<vector<ScoreType>>> bulge_score;
    vector<vector<ScoreType>> stacking_score;

    std::unordered_map<string, Lattice> aa_graphs_with_ln_w;

    std::vector<std::string> protein;
    std::unordered_map<std::string, std::string> aa_best_path_in_a_whole_codon;
    std::unordered_map<std::string, 
                       std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<FinalScoreType, NucType, NucType>, 
                       std::hash<std::tuple<NodeType, NodeType>>>> best_path_in_one_codon_unit;


    vector<ScoreInnerDate_t> reserved_scores;

#ifdef SPECIAL_HP
    unordered_map<NodeType, unordered_map<NodeType, unordered_map<int8_t, vector<tuple<string, ScoreType, FinalScoreType>>>, hash_pair>, hash_pair> hairpin_seq_score_cai;
    void special_hp(DFA_t& dfa, int8_t hairpin_length);
#endif
};



}
