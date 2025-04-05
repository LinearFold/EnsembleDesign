#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <cassert>
#include <unordered_map>
#include <map>
#include <set>
#include <unordered_set>

#include "optimizer.h"

using namespace std;

#define tetra_hex_tri -1
#define is_candidate_list
#define MIN_CUBE_PRUNING_SIZE 20

namespace EnsembleDesign {

template <PhaseOption phase>
void Optimizer::hairpin_beam(IndexType j, IndexType j_num, DFA_t& dfa) {

    if (beam_size > 0 and j_num == 0) {
        beam_prune<phase>(dfa, bestH, j);
    }

    auto j_node = make_pair(j, j_num);

    for (const auto& j_redge_ptr : dfa.right_edges[j_node]) {

        // Dereference the edge pointer and access tuple elements
        const auto& j1_node = std::get<1>(*j_redge_ptr); // Right node of the edge
        const auto& j_nuc   = std::get<2>(*j_redge_ptr); // Nucleotide/Weight
        const auto& j_param = std::get<3>(*j_redge_ptr); // Accessing the weight from the Parameter in the tuple

        if (j + 4 >= dfa.nodes.size()) continue;

        for (const auto& j4_node : dfa.nodes[j + 4]) {
            const auto& jnext_list = next_pair[j_nuc][j4_node];

            if (jnext_list.empty()) continue;

            for (const auto& jnext_node_nuc_param : jnext_list) { // Item Type: tuple<NodeType, NucType, Parameter*>
                const auto& jnext_node  = std::get<0>(jnext_node_nuc_param);
                const auto& jnext_nuc   = std::get<1>(jnext_node_nuc_param);
                const auto& jnext_param = std::get<2>(jnext_node_nuc_param);
                IndexType jnext = jnext_node.first;
                IndexType hairpin_length = jnext + 1 - j; //special hairpin

                NucPairType pair_nuc = NTP(j_nuc, jnext_nuc);

#ifdef SPECIAL_HP // TODO
                /**/
                if (hairpin_length == 5 or hairpin_length == 6 or hairpin_length == 8){
                    for(auto & seq_score_weight : hairpin_seq_score_cai[j_node][jnext_node][NTP(j_nuc, jnext_nuc)]){
                            auto sp_hairpin_seq =  get<0>(seq_score_weight);
                            auto sp_hairpin_path = dfa.findPath(sp_hairpin_seq, j_node, jnext_node);
                            assert(!sp_hairpin_path.empty());
                            ScoreType score_sp_hairpin =  get<1>(seq_score_weight) / kT;
                            Parameters params;
                            for (const auto& edge_ptr : sp_hairpin_path) {
                                const auto& param_ptr = std::get<3>(*edge_ptr);
                                params.push_back(param_ptr);
                            }

                            update_state<phase>(bestH[jnext_node][j_node][pair_nuc], score_sp_hairpin, params);
                    }

                    continue;
                }
#endif

                for (const auto &j1_redge_ptr : dfa.right_edges[j1_node]) {
                    const auto& j2_node  = std::get<1>(*j1_redge_ptr); // Right node of the edge
                    const auto& j2_num   = j2_node.second;
                    const auto& j1_nuc   = std::get<2>(*j1_redge_ptr); // Nucleotide/Weight
                    const auto& j1_param = std::get<3>(*j1_redge_ptr); // Accessing the weight from the Parameter in the tuple

                    for (const auto& jnext_1_node2ledges : dfa.auxiliary_left_edges[jnext_node]){

                        NodeType jnext_1_node = jnext_1_node2ledges.first;
                        NumType jnext_1_num = jnext_1_node.second;
                        if (jnext - j == 4 and (jnext_1_num != j2_num and dfa.nodes[j+2].size() == dfa.nodes[jnext-1].size())) continue;

                        for (const auto& jnext_ledge_ptr : jnext_1_node2ledges.second){

                            NucType jnext_1_nuc = std::get<2>(*jnext_ledge_ptr);
                            const auto& jnext_1_param = std::get<3>(*jnext_ledge_ptr);

                            ScoreType score_hairpin = - func4(j, jnext, j_nuc, j1_nuc, jnext_1_nuc, jnext_nuc, tetra_hex_tri) / kT;
                            Parameters params = {j_param, jnext_param, j1_param, jnext_1_param};

                            update_state<phase>(bestH[jnext_node][j_node][pair_nuc], score_hairpin, params);

                        }
                    }
                }
            }
        }
    }

    // for every state h in H[j]
    //    1. extend h(i, j) to h(i, jnext)
    //    2. generate p(i, j)

    for (auto &i_node_elem : bestH[j_node]) {
        auto i_node = i_node_elem.first;
        auto i = i_node.first;
        auto i_num = i_node.second;

        for (auto &pair_elem: bestH[j_node][i_node]) {
            NucPairType pair_nuc = pair_elem.first;

            auto i_nuc = PTLN(pair_nuc);
            auto j_nuc = PTRN(pair_nuc);


            for (const auto &j1_node2redges: dfa.auxiliary_right_edges[j_node]) { // Item Type: unordered_map<NodeType, vector<EdgeType*>>
                const auto &j1_node = j1_node2redges.first;
                const auto &jnext_list = next_pair[i_nuc][j1_node];

                if (jnext_list.empty()) continue;

                for (const auto &jnext_node_nuc_param: jnext_list) { // Item Type: tuple<NodeType, NucType, Parameter*>
                    const auto &jnext_node = std::get<0>(jnext_node_nuc_param);
                    const auto &jnext_nuc = std::get<1>(jnext_node_nuc_param);
                    const auto &jnext_param = std::get<2>(jnext_node_nuc_param);
                    IndexType jnext = jnext_node.first;
                    IndexType hairpin_length = jnext + 1 - i;

                    NucPairType  pair_nuc_i_jnext = NTP(i_nuc, jnext_nuc);

#ifdef SPECIAL_HP // TODO
                    /**/

                    if (hairpin_length == 5 or hairpin_length == 6 or hairpin_length == 8){
                        for(auto & seq_score_weight : hairpin_seq_score_cai[i_node][jnext_node][NTP(i_nuc, jnext_nuc)]){
                                auto sp_hairpin_seq =  get<0>(seq_score_weight);
                                auto sp_hairpin_path = dfa.findPath(sp_hairpin_seq, i_node, jnext_node);
                                assert(!sp_hairpin_path.empty());
                                ScoreType score_sp_hairpin =  get<1>(seq_score_weight) / kT;
                                Parameters params;
                                for (const auto& edge_ptr : sp_hairpin_path) {
                                    const auto& param = std::get<3>(*edge_ptr);
                                    params.push_back(param);
                                }

                                update_state<phase>(bestH[jnext_node][i_node][pair_nuc_i_jnext], score_sp_hairpin, params);
                        }

                        continue;
                    }

#endif

                    for (const auto &i_redge_ptr: dfa.right_edges[i_node]) {
                        NucType i_nuc_ = std::get<2>(*i_redge_ptr);
                        if (i_nuc != i_nuc_) continue;
                        NodeType i1_node = std::get<1>(*i_redge_ptr); // right node on the edge
                        const auto &i_param = std::get<3>(*i_redge_ptr);


                        for (const auto &i1_redge_ptr: dfa.right_edges[i1_node]) {
                            const auto &i2_node = std::get<1>(*i1_redge_ptr); // Right node of the edge
                            const auto &i1_nuc = std::get<2>(*i1_redge_ptr); // Nucleotide/Weight
                            const auto &i1_param = std::get<3>(
                                    *i1_redge_ptr); // Accessing the weight from the Parameter in the tuple

                            for (const auto &jnext_ledge_ptr: dfa.left_edges[jnext_node]) {
                                const auto &jnext_1_node = std::get<0>(*jnext_ledge_ptr); // Light node of the edge
                                const auto &jnext_1_nuc = std::get<2>(*jnext_ledge_ptr);
                                const auto &jnext_1_param = std::get<3>(*jnext_ledge_ptr);

                                ScoreType score_hairpin =
                                        -func4(i, jnext, i_nuc, i1_nuc, jnext_1_nuc, jnext_nuc,
                                                         tetra_hex_tri) / kT;
                                Parameters params = {i_param, jnext_param, i1_param, jnext_1_param};

                                update_state<phase>(bestH[jnext_node][i_node][pair_nuc_i_jnext], score_hairpin, params);

                            }
                        }
                    }
                }
            }

            auto& state_H = pair_elem.second;

            for (const auto &j_redge_ptr: dfa.right_edges[j_node]) {
                NucType j_nuc_ = std::get<2>(*j_redge_ptr);
                if (j_nuc != j_nuc_) continue;
                NodeType j1_node = get<1>(*j_redge_ptr); // right node on the edge

                update_state<phase>(bestP[j1_node][i_node][pair_nuc], state_H);
            }
        }
    }

}


template <PhaseOption phase>
void Optimizer::Multi_beam(IndexType j, IndexType j_num, DFA_t& dfa){
    if (beam_size > 0 and j_num == 0) {
        beam_prune<phase>(dfa, bestMulti, j);
    }

    NodeType j_node = make_pair(j, j_num);


    for (auto& i_node_elem : bestMulti[j_node]) {
        auto i_node = i_node_elem.first;
        auto i = i_node.first;

        for (auto &pair_elem: bestMulti[j_node][i_node]) {
            auto pair_nuc = pair_elem.first;
            auto i_nuc = PTLN(pair_nuc);
            auto j_1_nuc = PTRN(pair_nuc);

            auto &state_Multi = pair_elem.second;

            const auto &jnext_list = next_pair[i_nuc][j_node];

            if (!jnext_list.empty()) {

                for (const auto &jnext_node_nuc_param: jnext_list) { // Item Type: tuple<NodeType, NucType, Parameter*>
                    const auto &jnext_node = std::get<0>(jnext_node_nuc_param);
                    const auto &jnext_nuc = std::get<1>(jnext_node_nuc_param);
                    const auto &jnext_param = std::get<2>(jnext_node_nuc_param);
                    IndexType jnext = jnext_node.first;

                    for (const auto &jnext_redge_ptr: dfa.right_edges[jnext_node]) {
                        const auto &jnext1_node = std::get<1>(*jnext_redge_ptr); // right node of the edge
                        const auto &jnext_nuc_ = std::get<2>(*jnext_redge_ptr);
                        if (jnext_nuc != jnext_nuc_) continue;

                        Parameters params = {jnext_param};

                        NucPairType pair_nuc_i_jnext = NTP(i_nuc, jnext_nuc);

                        update_state<phase>(bestMulti[jnext1_node][i_node][pair_nuc_i_jnext], state_Multi, 0, params, state_Multi.pre_node);

                    }
                }
            }
            //  2. generate multi(i, j) -> p(i, j)

            ScoreType score_P_eq_Multi = -func12(i, j, i_nuc, -1, -1, j_1_nuc, seq_length) / kT;
            update_state<phase>(bestP[j_node][i_node][pair_nuc], state_Multi, score_P_eq_Multi);
        }
    }
}


template <PhaseOption phase>
void Optimizer::P_beam(IndexType j, IndexType j_num, DFA_t& dfa){
    if (beam_size > 0 and j_num == 0) {
        beam_prune<phase>(dfa, bestP, j); //ZL debug
    }
    auto j_node = make_pair(j, j_num);

    if (j < seq_length){

        for (auto &i_node_elem : bestP[j_node]){
            auto i_node = i_node_elem.first;
            auto i = i_node.first;

            if (i <= 0) continue;

            for (auto &pair_elem : bestP[j_node][i_node]) {
                auto pair_nuc = pair_elem.first;
                auto i_nuc = PTLN(pair_nuc);
                auto j_1_nuc = PTRN(pair_nuc);

                auto &state_P = pair_elem.second;

                // stacking
                for (const auto &j_redge_ptr: dfa.right_edges[j_node]) {
                    const auto &j1_node = std::get<1>(*j_redge_ptr);
                    const auto &j_nuc = std::get<2>(*j_redge_ptr);
                    const auto &j_param = std::get<3>(*j_redge_ptr);

                    for (const auto &i_ledge_ptr: dfa.left_edges[i_node]) {
                        const auto &i_1_node = std::get<0>(*i_ledge_ptr);
                        const auto &i_1_nuc = std::get<2>(*i_ledge_ptr);
                        const auto &i_1_param = std::get<3>(*i_ledge_ptr);
                        auto outer_pair = NTP(i_1_nuc, j_nuc);
                        if (_allowed_pairs[i_1_nuc][j_nuc]) {
                            ScoreType score_stacking = stacking_score[outer_pair - 1][pair_nuc - 1] / kT;
                            Parameters params = {i_1_param, j_param};

                            update_state<phase>(bestP[j1_node][i_1_node][outer_pair], state_P, score_stacking, params);
                        }
                    }
                }

                // right bulge: ((...)..)
                for (const auto &j1_node2redges: dfa.auxiliary_right_edges[j_node]) {
                    auto j1_node = j1_node2redges.first;

                    for (const auto &i_ledge_ptr: dfa.left_edges[i_node]) {
                        const auto &i_1_node = std::get<0>(*i_ledge_ptr);
                        const auto &i_1_nuc = std::get<2>(*i_ledge_ptr);
                        const auto &i_1_param = std::get<3>(*i_ledge_ptr);

                        const auto &q_list = next_list[i_1_nuc][j1_node];

                        for (const auto &q_node_nuc_param: q_list) { // Item Type: tuple<NodeType, NucType, Parameter*>

                            auto q_node = std::get<0>(q_node_nuc_param);

                            auto q = q_node.first;
                            auto q_num = q_node.second;

                            if (q - j > SINGLE_MAX_LEN) break;

                            auto q_nuc = std::get<1>(q_node_nuc_param);
                            const auto &q_param = std::get<2>(q_node_nuc_param);
                            auto outer_pair = NTP(i_1_nuc, q_nuc);

                            for (const auto &q1_node2redges: dfa.auxiliary_right_edges[q_node]) {
                                auto q1_node = q1_node2redges.first;
                                auto q_first_redge_ptr = q1_node2redges.second[0];
                                auto q_nuc_ = std::get<2>(*q_first_redge_ptr);
                                if (dfa.nodes[q].size() == 1 and dfa.nodes[q + 1].size() == 2 and
                                    q_nuc != q_nuc_)
                                    continue;

                                ScoreType score_bulge = bulge_score[outer_pair - 1][pair_nuc - 1][q - j - 1] / kT;
                                Parameters params = {i_1_param, q_param};

                                update_state<phase>(bestP[q1_node][i_1_node][outer_pair], state_P, score_bulge, params);

                                break;
                            }
                        }
                    }
                }

                // left bulge: (..(...))
                for (const auto &j_redge_ptr: dfa.right_edges[j_node]) {
                    const auto &j1_node = std::get<1>(*j_redge_ptr);
                    const auto &j_nuc = std::get<2>(*j_redge_ptr);
                    const auto &j_param = std::get<3>(*j_redge_ptr);

                    for (const auto &i_1_node2ledges: dfa.auxiliary_left_edges[i_node]) {
                        auto i_1_node = i_1_node2ledges.first;
                        const auto &p_list = prev_list[j_nuc][i_1_node];

                        for (const auto &p_node_nuc_param: p_list) { // Item Type: tuple<NodeType, NucType, Parameter*>

                            auto p_node = std::get<0>(p_node_nuc_param);

                            auto p = p_node.first;
                            auto p_num = p_node.second;

                            if (i - p > SINGLE_MAX_LEN) break;

                            auto p_1_nuc = std::get<1>(p_node_nuc_param);
                            auto outer_pair = NTP(p_1_nuc, j_nuc);

                            for (const auto &p_ledge_ptr: dfa.left_edges[p_node]) {

                                const auto &p_1_nuc_ = std::get<2>(*p_ledge_ptr);
                                if (p_1_nuc != p_1_nuc_) continue;

                                const auto &p_1_node = std::get<0>(*p_ledge_ptr);
                                const auto &p_1_param = std::get<3>(*p_ledge_ptr);

                                ScoreType score_bulge = bulge_score[outer_pair - 1][pair_nuc - 1][i - p - 1] / kT;
                                Parameters params = {p_1_param, j_param};

                                update_state<phase>(bestP[j1_node][p_1_node][outer_pair], state_P, score_bulge, params);
                            }
                        }
                    }
                }

                // internal loop (...(...)...)
                for (const auto &j1_node2redges: dfa.auxiliary_right_edges[j_node]) {
                    auto j1_node = j1_node2redges.first;
                    auto j1_num = j1_node.second;

                    for (const auto &i_ledge_ptr: dfa.left_edges[i_node]) {
                        const auto &i_1_node = std::get<0>(*i_ledge_ptr);
                        const auto &i_1_num = i_1_node.second;
                        const auto &i_1_nuc = std::get<2>(*i_ledge_ptr);
                        const auto &i_1_param = std::get<3>(*i_ledge_ptr);

                        for (IndexType p = i - 1;
                             p > max(i - SINGLE_MAX_LEN, 0); --p) {//ZL, i-(p-1)<=len => i - len < p
                            vector <NodeType> p_node_list;

                            if (p == i - 1)
                                p_node_list.push_back(i_1_node);
                            else if (p == i - 2) // hzhang: N.B. add this p, i-1, i o--o--o
                                for (const auto &p_node_dict: dfa.auxiliary_left_edges[i_1_node])
                                    p_node_list.push_back(p_node_dict.first);
                            else
                                p_node_list = dfa.nodes[p];

                            for (const auto &p_node: p_node_list) {
                                for (const auto &p_redge_ptr: dfa.right_edges[p_node]) {

                                    const auto &p1_node = std::get<1>(*p_redge_ptr);
                                    const auto &p1_num = p1_node.second;
                                    const auto &p_nuc = std::get<2>(*p_redge_ptr);
                                    const auto &p_param = std::get<3>(*p_redge_ptr);

                                    if (p == i - 1 and p_nuc != i_1_nuc) continue;
                                    else if (p == i - 2 and p1_num != i_1_num) continue;
                                    else if (p == i - 3 and p1_num != i_1_num and
                                             dfa.nodes[p + 1].size() == dfa.nodes[i - 1].size())
                                        continue;

                                    for (const auto &p_ledge_ptr: dfa.left_edges[p_node]) {
                                        const auto &p_1_node = std::get<0>(*p_ledge_ptr);
                                        const auto &p_1_nuc = std::get<2>(*p_ledge_ptr);
                                        const auto &p_1_param = std::get<3>(*p_ledge_ptr);

                                        auto q_list = next_list[p_1_nuc][j1_node];

                                        for (const auto &q_node_nuc_param: q_list) { // Item Type: tuple<NodeType, NucType, Parameter*>

                                            const auto &q_node = std::get<0>(q_node_nuc_param);
                                            const auto &q = q_node.first;
                                            const auto &q_num = q_node.second;

                                            if (i - p + q - j >
                                                SINGLE_MAX_LEN) //check if q is still in the internal loop limit boundary.
                                                break;

                                            const auto &q_nuc = std::get<1>(q_node_nuc_param);
                                            const auto &q_param = std::get<2>(q_node_nuc_param);

                                            for (const auto &q1_node2redges: dfa.auxiliary_right_edges[q_node]) {
                                                const auto &q1_node = q1_node2redges.first;
                                                const auto &q_first_redge_ptr = q1_node2redges.second[0];
                                                const auto &q_nuc_ = std::get<2>(*q_first_redge_ptr);

                                                if (dfa.nodes[q].size() == 1 and dfa.nodes[q + 1].size() == 2 and
                                                    q_nuc != q_nuc_)
                                                    continue;

                                                auto &new_state_P = bestP[q1_node][p_1_node][NTP(p_1_nuc, q_nuc)];

                                                // p_1 p ... i_1 i ... j_1 j ... q_1 q

                                                for (const auto &j_redge_ptr: j1_node2redges.second) {
                                                    const auto &j_nuc = std::get<2>(
                                                            *j_redge_ptr);   //nucj_weightj.first;
                                                    const auto &j_param = std::get<3>(
                                                            *j_redge_ptr); //nucj_weightj.second;

                                                    if (q == j + 1) { // case: p_1 p ... i_1 i ... j_1 j q
                                                        ScoreType score_internal =
                                                                -func5(p - 1, q, i, j - 1, p_1_nuc, p_nuc,
                                                                                j_nuc, q_nuc, i_1_nuc, i_nuc, j_1_nuc,
                                                                                j_nuc) / kT;
                                                        Parameters params = {p_param, p_1_param, j_param, q_param};

                                                        if (p == i - 1) { // case: p_1 p i ... j_1 j q
                                                        } else {
                                                            params.push_back(i_1_param);
                                                        }

                                                        update_state<phase>(new_state_P, state_P, score_internal, params);
                                                    } else if (q == j + 2) { // case: p_1 p ... i_1 i ... j_1 j q_1 q
                                                        for (const auto &q_1_node2ledges: dfa.auxiliary_left_edges[q_node]) {
                                                            const auto &q_1_node = q_1_node2ledges.first;
                                                            const auto &q_1_num = q_1_node.second;
                                                            if (q_1_num != j1_num) continue;

                                                            for (const auto &q_ledge_ptr: q_1_node2ledges.second) {
                                                                const auto &q_1_nuc = std::get<2>(*q_ledge_ptr);
                                                                const auto &q_1_param = std::get<3>(*q_ledge_ptr);
                                                                ScoreType score_internal =
                                                                        -func5(p - 1, q, i, j - 1, p_1_nuc,
                                                                                        p_nuc, q_1_nuc, q_nuc, i_1_nuc,
                                                                                        i_nuc, j_1_nuc, j_nuc) / kT;
                                                                Parameters params = {p_param, p_1_param, j_param,
                                                                                     q_1_param, q_param};

                                                                 if (p == i - 1) { // case: p_1 p i ... j_1 j q_1 q
                                                                } else {
                                                                    params.push_back(i_1_param);
                                                                }

                                                                update_state<phase>(new_state_P, state_P, score_internal, params);
                                                            }
                                                            if (dfa.nodes[q - 1].size() == 2) break;
                                                        }
                                                    } else if (q == j + 3) { // case: p_1 p ... i_1 i ... j_1 j . q_1 q
                                                        for (const auto &q_1_node2ledges: dfa.auxiliary_left_edges[q_node]) {
                                                            const auto &q_1_node = q_1_node2ledges.first;
                                                            const auto &q_1_num = q_1_node.second;
                                                            if (q_1_num != j1_num and
                                                                dfa.nodes[q - 1].size() == dfa.nodes[j + 1].size())
                                                                continue;
                                                            for (const auto &q_ledge_ptr: q_1_node2ledges.second) {
                                                                const auto &q_1_nuc = std::get<2>(*q_ledge_ptr);
                                                                const auto &q_1_param = std::get<3>(*q_ledge_ptr);
                                                                ScoreType score_internal =
                                                                        -func5(p - 1, q, i, j - 1, p_1_nuc,
                                                                                        p_nuc, q_1_nuc, q_nuc, i_1_nuc,
                                                                                        i_nuc, j_1_nuc, j_nuc) / kT;
                                                                Parameters params = {p_param, p_1_param, j_param,
                                                                                     q_1_param, q_param};

                                                                if (p == i - 1) { // case: p_1 p i ... j_1 j . q_1 q
                                                                } else {
                                                                    params.push_back(i_1_param);
                                                                }


                                                                update_state<phase>(new_state_P, state_P, score_internal, params);
                                                            }
                                                            if (dfa.nodes[q - 1].size() == 2) break;
                                                        }
                                                    } else {// case: p_1 p ... i_1 i ... j_1 j ... q_1 q
                                                        for (const auto &q_1_node2ledges: dfa.auxiliary_left_edges[q_node]) {
                                                            const auto &q_1_node = q_1_node2ledges.first;
                                                            for (const auto &q_ledge_ptr: q_1_node2ledges.second) {
                                                                const auto &q_1_nuc = std::get<2>(*q_ledge_ptr);
                                                                const auto &q_1_param = std::get<3>(*q_ledge_ptr);
                                                                ScoreType score_internal =
                                                                        -func5(p - 1, q, i, j - 1, p_1_nuc,
                                                                                        p_nuc, q_1_nuc, q_nuc, i_1_nuc,
                                                                                        i_nuc, j_1_nuc, j_nuc) / kT;
                                                                Parameters params = {p_param, p_1_param, j_param, q_1_param, q_param};

                                                                //double weight_left;
                                                                if (p == i - 1) {
                                                                } else {
                                                                    params.push_back(i_1_param);
                                                                }

                                                                update_state<phase>(new_state_P, state_P, score_internal, params);

                                                            }
                                                            if (dfa.nodes[q - 1].size() == 2) break;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // M = P and M_P = P

    bool use_cube_pruning = enable_cube_pruning && beam_size > MIN_CUBE_PRUNING_SIZE && bestP[j_node].size() > MIN_CUBE_PRUNING_SIZE;

    // M = P and M2 = M + P
    for (auto &i_node_elem : bestP[j_node]) {
        auto i_node = i_node_elem.first;
        auto i = i_node.first;

        for (auto &pair_elem: bestP[j_node][i_node]) {
            auto pair_nuc = pair_elem.first;
            auto i_nuc = PTLN(pair_nuc);
            auto j_1_nuc = PTRN(pair_nuc);
            auto &state_P = pair_elem.second;

            if (i > 0 and j < seq_length) {

                ScoreType score_M_eq_P =
                        -func9(i, j - 1, j - 1, -1, i_nuc, j_1_nuc, -1, seq_length) / kT;


                update_state<phase>(bestM[j_node][i_node], state_P, score_M_eq_P);

                // M2 = M + P
                if(!bestM[i_node].empty() && !use_cube_pruning){
                    //cout<<"not using cube pruning "<<endl;
                    for (auto &m_node_m_state : bestM[i_node]){
                        auto m_node = m_node_m_state.first;
                        auto& state_M = m_node_m_state.second;

                        update_state<phase>(bestM2[j_node][m_node], state_M, state_P, score_M_eq_P);
                    }
                }
            }

            // C = C + P
            ScoreType score_C_eq_C_plus_P = - func14(i, j-1, i_nuc, j_1_nuc, seq_length) / kT;

            if (i > 0){
                auto& state_C = bestC[i_node];
                if (state_C.inside == util::value_min<ScoreType>()) continue;
                update_state<phase>(bestC[j_node], state_C, state_P, score_C_eq_C_plus_P);
            }
            else{
                update_state<phase>(bestC[j_node], state_P, score_C_eq_C_plus_P);
            }
        }
    }

}

template <PhaseOption phase>
void Optimizer::M2_beam(IndexType j, IndexType j_num, DFA_t& dfa){
    if (beam_size > 0 and j_num == 0) {
        beam_prune<phase>(dfa, bestM2, j); //ZL debug
    }

    auto j_node = make_pair(j, j_num);

    for (auto &i_node_state : bestM2[j_node]){
        auto i_node = i_node_state.first;
        auto& state_M2 = i_node_state.second;
        auto i = i_node.first;

        // 1. multi-loop
        for (IndexType p = i-1; p >= max(i - SINGLE_MAX_LEN, 0); --p){
            vector<NodeType> p_node_list;
            if (p == i - 1)
                for(const auto& p_node_dict : dfa.auxiliary_left_edges[i_node])
                    p_node_list.push_back(p_node_dict.first);
            else p_node_list = dfa.nodes[p];

            for (const auto &p_node : p_node_list){
                for (const auto &p_redge_ptr : dfa.right_edges[p_node]){
                    const auto& p1_node = std::get<1>(*p_redge_ptr); // right node of the edge
                    const auto& p_nuc   = std::get<2>(*p_redge_ptr);
                    const auto& p_param = std::get<3>(*p_redge_ptr);

                    if(p == i - 1 and p1_node != i_node) continue;
                    if(p == i - 2 and dfa.nodes[p+1].size() == dfa.nodes[i].size() and p1_node.second != i_node.second) continue;

                    const auto& q_list = next_pair[p_nuc][j_node];

                    for (const auto& q_node_nuc_param : q_list){ // Item Type: tuple<NodeType, NucType, Parameter*>
                        const auto& q_node  = std::get<0>(q_node_nuc_param);
                        const auto& q_nuc   = std::get<1>(q_node_nuc_param);
                        const auto& q_param = std::get<2>(q_node_nuc_param);
                        IndexType q = q_node.first;

                        if (i - p + q - j - 1 > SINGLE_MAX_LEN) continue; //ZL, i-p-1+q-j
                        auto outer_pair = NTP(p_nuc, q_nuc);
                        for (const auto& q_redge_ptr : dfa.right_edges[q_node]){
                            const auto& q_nuc_ = std::get<2>(*q_redge_ptr);
                            if (q_nuc != q_nuc_) continue;
                            const auto& q1_node = std::get<1>(*q_redge_ptr);
                            Parameters params = {p_param, q_param};

                            update_state<phase>(bestMulti[q1_node][p_node][outer_pair], state_M2, 0, params, j_node);
                            break;
                        }
                    }
                }
            }
        }
        //  2. M = M2
        update_state<phase>(bestM[j_node][i_node], state_M2);
    }
}


template <PhaseOption phase>
void Optimizer::M_beam(IndexType j, IndexType j_num, DFA_t& dfa)
{
    ScoreType threshold = util::value_min<ScoreType>();
    if (beam_size > 0 and j_num == 0) {
        threshold = beam_prune<phase>(dfa, bestM, j);
    }

    auto j_node = make_pair(j, j_num);
    // M = M + unpaired
    for (auto &i_node_state : bestM[j_node]){

        auto i_node = i_node_state.first;
        auto& state_M = i_node_state.second;

        for (const auto& j_redge_ptr : dfa.right_edges[j_node]){
            const auto& j1_node = std::get<1>(*j_redge_ptr);
            const auto& j_nuc   = std::get<2>(*j_redge_ptr);
            const auto& j_param = std::get<3>(*j_redge_ptr);

            Parameters params = {j_param};
            update_state<phase>(bestM[j1_node][i_node], state_M, 0, params);
        }
    }
}

template <PhaseOption phase>
void Optimizer::C_beam(IndexType j, IndexType j_num, DFA_t& dfa)
{
    auto j_node = make_pair(j, j_num);

    auto& state_C = bestC[j_node];
    //  C = C + unpaired
    for (const auto& j_redge_ptr : dfa.right_edges[j_node]){
        const auto& j1_node = std::get<1>(*j_redge_ptr);
        const auto& j_nuc   = std::get<2>(*j_redge_ptr);
        const auto& j_param = std::get<3>(*j_redge_ptr);

        Parameters params = {j_param};
        update_state<phase>(bestC[j1_node], state_C, 0, params);
    }
}

void Optimizer::get_next_pair(DFA_t& dfa) {
    //vector<tuple<NodeType, NucType, double>> temp_vector;
    vector<NextPairType> temp_vector;
    for (NucType i_nuc = 0; i_nuc < NOTON; i_nuc++) {
        for (IndexType j = seq_length; j > 0; j--) {
            for (const auto& j_node : dfa.nodes[j]) {
                for (const auto& j_1_node2ledges : dfa.auxiliary_left_edges[j_node]) {
                    temp_vector.clear();
                    const auto& j_1_node = j_1_node2ledges.first;
                    for (const auto& j_ledge_ptr : j_1_node2ledges.second){
                        const auto& j_1_nuc = std::get<2>(*j_ledge_ptr);
                        const auto& j_1_param = std::get<3>(*j_ledge_ptr);
                        if (_allowed_pairs[i_nuc][j_1_nuc])
                            temp_vector.push_back(make_tuple(j_1_node, j_1_nuc, j_1_param));
                    }
                    if(temp_vector.size() == 0){
                        if (next_pair[i_nuc][j_1_node].size() > 0 and next_pair[i_nuc][j_node].size() > 0) {
                            // merge
                            IndexType index1 = std::get<0>(next_pair[i_nuc][j_1_node][0]).first;
                            IndexType index2 = std::get<0>(next_pair[i_nuc][j_node][0]).first;
                            if(index1/3 == index2/3)
                                next_pair[i_nuc][j_1_node].insert(next_pair[i_nuc][j_1_node].end(),
                                                                     next_pair[i_nuc][j_node].begin(),
                                                                     next_pair[i_nuc][j_node].end());
                            else if(index1 > index2){
                                next_pair[i_nuc][j_1_node].clear();
                                next_pair[i_nuc][j_1_node].insert(next_pair[i_nuc][j_1_node].end(),
                                                                     next_pair[i_nuc][j_node].begin(),
                                                                     next_pair[i_nuc][j_node].end());
                            }
                        }else if (next_pair[i_nuc][j_node].size() > 0)
                            next_pair[i_nuc][j_1_node].insert(next_pair[i_nuc][j_1_node].end(),
                                                                     next_pair[i_nuc][j_node].begin(),
                                                                     next_pair[i_nuc][j_node].end());
                    }
                    else
                        next_pair[i_nuc][j_1_node].insert(next_pair[i_nuc][j_1_node].end(),
                                                                     temp_vector.begin(),
                                                                     temp_vector.end());
                }
            }
        }
    }
}

void Optimizer::get_next_pair_set() {

    for(NucType i_nuc=0; i_nuc < NOTON; i_nuc++){
        for (const auto& j_node_vnuc : next_pair[i_nuc]) {
            NodeType j_node = j_node_vnuc.first;
            next_pair_set[i_nuc][j_node] =
                set<NextPairType>(j_node_vnuc.second.begin(), j_node_vnuc.second.end());
        }
    }
    for(NucType i_nuc=0; i_nuc < NOTON; i_nuc++){
        for (const auto& j_node_vnuc : next_pair_set[i_nuc]) {
            NodeType j_node = j_node_vnuc.first;
            next_pair[i_nuc][j_node].clear();
            for(const auto& item : next_pair_set[i_nuc][j_node]){
                next_pair[i_nuc][j_node].push_back(item);
            }
        }
    }
}

void Optimizer::get_prev_pair(DFA_t& dfa) {
    vector<NextPairType> temp_vector;
    for (NucType i_nuc = 0; i_nuc < NOTON; i_nuc++) {
        for (IndexType j = 0; j < seq_length; j++) {
            for (const auto& j_node : dfa.nodes[j]) {
                for (const auto& j1_node2redges : dfa.auxiliary_right_edges[j_node]) {
                    temp_vector.clear();
                    const auto& j1_node = j1_node2redges.first;
                    for (const auto& j_redge_ptr : j1_node2redges.second){
                        const auto& j_nuc = std::get<2>(*j_redge_ptr);
                        const auto& j_param = std::get<3>(*j_redge_ptr);
                        if (_allowed_pairs[i_nuc][j_nuc])
                            temp_vector.push_back(make_tuple(j1_node, j_nuc, j_param));
                    }
                    if(temp_vector.size() == 0){
                        if (prev_pair[i_nuc][j1_node].size() > 0 and prev_pair[i_nuc][j_node].size() > 0) {
                            // merge
                            IndexType index1 = std::get<0>(prev_pair[i_nuc][j1_node][0]).first-1;
                            IndexType index2 = std::get<0>(prev_pair[i_nuc][j_node][0]).first-1;
                            if(index1/3 == index2/3)
                                prev_pair[i_nuc][j1_node].insert(prev_pair[i_nuc][j1_node].end(),
                                                                     prev_pair[i_nuc][j_node].begin(),
                                                                     prev_pair[i_nuc][j_node].end());
                            else if(index1 < index2){
                                prev_pair[i_nuc][j1_node].clear();
                                prev_pair[i_nuc][j1_node].insert(prev_pair[i_nuc][j1_node].end(),
                                                                     prev_pair[i_nuc][j_node].begin(),
                                                                     prev_pair[i_nuc][j_node].end());
                            }
                        }else if (prev_pair[i_nuc][j_node].size() > 0)
                            prev_pair[i_nuc][j1_node].insert(prev_pair[i_nuc][j1_node].end(),
                                                                     prev_pair[i_nuc][j_node].begin(),
                                                                     prev_pair[i_nuc][j_node].end());
                    }
                    else
                        prev_pair[i_nuc][j1_node].insert(prev_pair[i_nuc][j1_node].end(),
                                                                     temp_vector.begin(),
                                                                     temp_vector.end());
                }
            }
        }
    }
}


void Optimizer::get_prev_pair_set() {

    for(NucType i_nuc=0; i_nuc < NOTON; i_nuc++){
        for (const auto& j_node_vnuc : prev_pair[i_nuc]) {
            NodeType j_node = j_node_vnuc.first;
            prev_pair_set[i_nuc][j_node] =
                set<NextPairType>(j_node_vnuc.second.begin(), j_node_vnuc.second.end());
        }
    }
    for(NucType i_nuc=0; i_nuc < NOTON; i_nuc++){
        for (const auto& j_node_vnuc : prev_pair_set[i_nuc]) {
            NodeType j_node = j_node_vnuc.first;
            prev_pair[i_nuc][j_node].clear();
            for(const auto& item : prev_pair_set[i_nuc][j_node]){
                prev_pair[i_nuc][j_node].push_back(item);
            }
        }
    }
}

#ifdef SPECIAL_HP
void Optimizer::special_hp(DFA_t& dfa, int8_t hairpin_length) {
    int8_t hairpin_type = HAIRPINTYPE(hairpin_length);
    vector<tuple<NodeType, string, double, NodeType>> queue;
    vector<tuple<NodeType, string, double, NodeType>> frontier; 
    // vector
    for(IndexType i=0; i<=seq_length - hairpin_length; i++){
        for(NodeType i_node : dfa.nodes[i]){
            int count = hairpin_length;
            queue.clear();
            queue.push_back(make_tuple(i_node, "", double(0.), i_node));
            while(count > 0){
                count --;
                frontier.clear();
                for(const auto& node_str : queue){
                    NodeType cur_node = std::get<0>(node_str);
                    string cur_str = std::get<1>(node_str);
                    double cur_lncai = std::get<2>(node_str);
                    for(const auto& redge_ptr : dfa.right_edges[cur_node]){
                        NodeType new_node = std::get<1>(*redge_ptr);
                        string new_str = cur_str + GET_ACGU(std::get<2>(*redge_ptr));
                        double new_total_lncai = cur_lncai + std::get<3>(*redge_ptr)->cai_score;
                        frontier.push_back(make_tuple(new_node, new_str, new_total_lncai, cur_node));
                    }
                }
                queue.swap(frontier);
            }
            for(auto node_str : queue){
                auto j_node = std::get<3>(node_str);
                auto temp_seq = std::get<1>(node_str);
                auto cai_score = std::get<2>(node_str);
                auto hairpin_length = temp_seq.size();
                int8_t hairpin_type = HAIRPINTYPE(hairpin_length);
                NucType i_nuc = GET_ACGU_NUC(temp_seq[0]);
                NucType j_nuc = GET_ACGU_NUC(temp_seq[temp_seq.size() - 1]);
                auto temp_nucpair = NTP(i_nuc, j_nuc);

                ScoreType special_hairpin_score = func2(temp_seq, hairpin_type);
                if(special_hairpin_score == SPECIAL_HAIRPIN_SCORE_BASELINE){

                    auto newscore = - func4(0, hairpin_length - 1, GET_ACGU_NUC(temp_seq[0]), GET_ACGU_NUC(temp_seq[1]), GET_ACGU_NUC(temp_seq[hairpin_length-2]), GET_ACGU_NUC(temp_seq[hairpin_length-1]), tetra_hex_tri);
                    hairpin_seq_score_cai[i_node][j_node][temp_nucpair].push_back(make_tuple(temp_seq, newscore, cai_score));
                }

                else{
                    hairpin_seq_score_cai[i_node][j_node][temp_nucpair].push_back(make_tuple(temp_seq, special_hairpin_score, cai_score));
                }
            }
        }
    }
}
#endif

void Optimizer::preprocess(DFA_t& dfa) {

    vector<NextPairType> new_q_list, new_p_list;
    set<NodeType> visited;

    // next_list
    NodeType init_node = make_pair(0, 0);
    for (NucType i_nuc=1; i_nuc < NOTON; i_nuc++) {
        visited.clear();
        auto q_list = next_pair[i_nuc][init_node];

        while (!q_list.empty()){
            new_q_list.clear();

            for (const auto& q_node_nuc_param: q_list){

                auto q_node = std::get<0>(q_node_nuc_param);
                auto q = q_node.first;
                auto q_num = q_node.second;
                auto q_nuc = std::get<1>(q_node_nuc_param);

                // q_node
                next_list[i_nuc][q_node].push_back(q_node_nuc_param);
                // q-1 is special
                for(const auto& q_1_node2ledges : dfa.auxiliary_left_edges[q_node]){
                    NodeType q_1_node = q_1_node2ledges.first;
                    next_list[i_nuc][q_1_node].push_back(q_node_nuc_param);
                }
                for(IndexType j=q-2; j>=max(0, q-SINGLE_MAX_LEN-1); j--)
                    for(NodeType j_node : dfa.nodes[j])
                        next_list[i_nuc][j_node].push_back(q_node_nuc_param);

                for(const auto& q1_node2redges : dfa.auxiliary_right_edges[q_node]){
                    const auto& q1_node = q1_node2redges.first;
                    const auto& q_first_redge_ptr = q1_node2redges.second[0];
                    const auto& q_nuc_  = std::get<2>(*q_first_redge_ptr);

                    if(dfa.nodes[q].size() == 1 and dfa.nodes[q+1].size() == 2 and q_nuc != q_nuc_) continue;

                
                    if(visited.find(q1_node) == visited.end()){
                        visited.insert(q1_node);
                        new_q_list.insert(new_q_list.end(), next_pair[i_nuc][q1_node].cbegin(), next_pair[i_nuc][q1_node].cend());
                    }
                    break;
                }
            }
            q_list.swap(new_q_list);
        }
    }

    // prev_list
    init_node = make_pair(seq_length, 0);
    for (NucType j_nuc=1; j_nuc < NOTON; j_nuc++) {
        visited.clear();
        auto p_list = prev_pair[j_nuc][init_node];

        while (!p_list.empty()){
            new_p_list.clear();

            for (const auto& p_node_nuc_param : p_list){
                auto p_node = std::get<0>(p_node_nuc_param);
                auto p = p_node.first;
                auto p_num = p_node.second;
                auto p_1_nuc = std::get<1>(p_node_nuc_param);

                // p_node
                prev_list[j_nuc][p_node].push_back(p_node_nuc_param);
                // p+1 is special
                for(const auto& p1_node2redges : dfa.auxiliary_right_edges[p_node]){
                    NodeType p1_node = p1_node2redges.first;
                    prev_list[j_nuc][p1_node].push_back(p_node_nuc_param);
                }
                for(IndexType i=p+2; i<=min(seq_length, p+SINGLE_MAX_LEN+1); i++)
                    for(NodeType i_node : dfa.nodes[i])
                        prev_list[j_nuc][i_node].push_back(p_node_nuc_param);

                for(const auto& p_ledge_ptr : dfa.left_edges[p_node]){
                                    
                    NucType p_1_nuc_ = std::get<2>(*p_ledge_ptr);
                    if(p_1_nuc != p_1_nuc_) continue;
                    NodeType p_1_node = std::get<0>(*p_ledge_ptr);
                
                    if(visited.find(p_1_node) == visited.end()){
                        visited.insert(p_1_node);
                        new_p_list.insert(new_p_list.end(), prev_pair[j_nuc][p_1_node].cbegin(), prev_pair[j_nuc][p_1_node].cend());
                    }
                }
            }
            p_list.swap(new_p_list);
        }
    }

    // stacking energy computation
    ScoreType newscore;
    for(int8_t outer_pair=1; outer_pair<=6; outer_pair++){
        auto i_1_nuc = PTLN(outer_pair);
        auto j_nuc = PTRN(outer_pair);
        for(int8_t inner_pair=1; inner_pair<=6; inner_pair++){
            auto i_nuc = PTLN(inner_pair);
            auto j_1_nuc = PTRN(inner_pair);
            newscore = - func5(0, 1, 1, 0,
                                i_1_nuc, i_nuc, j_1_nuc, j_nuc,
                                i_1_nuc, i_nuc, j_1_nuc, j_nuc);
            stacking_score[outer_pair-1][inner_pair-1] = newscore;

            for (IndexType l=0; l<=SINGLE_MAX_LEN; l++){
                newscore = - func5(0, l+2, 1, 0,
                              i_1_nuc, i_nuc, j_1_nuc, j_nuc,
                              i_1_nuc, i_nuc, j_1_nuc, j_nuc);

                bulge_score[outer_pair-1][inner_pair-1][l] = newscore;
            }
        }   
    }
#ifdef SPECIAL_HP
    // Triloops
    special_hp(dfa, 5);
    // Tetraloop37
    special_hp(dfa, 6);
    // Hexaloops
    special_hp(dfa, 8);
#endif
}


void Optimizer::clear_parser_states() {
    bestC.clear();
    bestH.clear();
    bestP.clear();
    bestMulti.clear();
    bestM2.clear();
    bestM.clear();

}

void Optimizer::forward(DFA_t& dfa){
    auto start_node = make_pair(0, 0);
    auto end_node = make_pair(seq_length, 0);

    bestC[start_node].inside = 0;
    C_beam<PhaseOption::Inside>(0, 0, dfa);
    C_beam<PhaseOption::Inside>(0, 1, dfa);
    hairpin_beam<PhaseOption::Inside>(0, 0, dfa);
    hairpin_beam<PhaseOption::Inside>(0, 1, dfa);

    for (IndexType j = 1; j <= seq_length; ++j) {

        hairpin_beam<PhaseOption::Inside>(j, 0, dfa);
        hairpin_beam<PhaseOption::Inside>(j, 1, dfa);
        Multi_beam<PhaseOption::Inside>(j, 0, dfa);
        Multi_beam<PhaseOption::Inside>(j, 1, dfa);
        P_beam<PhaseOption::Inside>(j, 0, dfa);
        P_beam<PhaseOption::Inside>(j, 1, dfa);
        M2_beam<PhaseOption::Inside>(j, 0, dfa);
        M2_beam<PhaseOption::Inside>(j, 1, dfa);

        if (j < seq_length) {
            M_beam<PhaseOption::Inside>(j, 0, dfa);
            M_beam<PhaseOption::Inside>(j, 1, dfa);
            C_beam<PhaseOption::Inside>(j, 0, dfa);
            C_beam<PhaseOption::Inside>(j, 1, dfa);
        }

    }

}

void Optimizer::backward(DFA_t& dfa){
    auto start_node = make_pair(0, 0);
    auto end_node = make_pair(seq_length, 0);
    bestC[end_node].outside = 0;


    for (IndexType j = seq_length; j > 0 ; --j) {

        if (j < seq_length) {
            C_beam<PhaseOption::Outside>(j, 0, dfa);
            C_beam<PhaseOption::Outside>(j, 1, dfa);
            M_beam<PhaseOption::Outside>(j, 0, dfa);
            M_beam<PhaseOption::Outside>(j, 1, dfa);
        }

        M2_beam<PhaseOption::Outside>(j, 0, dfa);
        M2_beam<PhaseOption::Outside>(j, 1, dfa);
        P_beam<PhaseOption::Outside>(j, 0, dfa);
        P_beam<PhaseOption::Outside>(j, 1, dfa);
        Multi_beam<PhaseOption::Outside>(j, 0, dfa);
        Multi_beam<PhaseOption::Outside>(j, 1, dfa);
        hairpin_beam<PhaseOption::Outside>(j, 0, dfa);
        hairpin_beam<PhaseOption::Outside>(j, 1, dfa);

    }

    hairpin_beam<PhaseOption::Outside>(0, 0, dfa);
    hairpin_beam<PhaseOption::Outside>(0, 1, dfa);
    C_beam<PhaseOption::Outside>(0, 0, dfa);
    C_beam<PhaseOption::Outside>(0, 1, dfa);

}


ScoreType Nesterov_coef(ScoreType a){
    return (1 + std::sqrt(4 * a * a + 1)) / 2;
}

void Optimizer::optimize(
        DFA_t& dfa, Codon& codon, std::string& aa_seq, std::vector<std::string>& p, 
        std::unordered_map<std::string, std::string>& aa_best_in_codon,
        std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<double, NucType, NucType>,
            std::hash<std::tuple<NodeType, NodeType>>>>& best_path_in_one_codon,
            std::unordered_map<string, Lattice>& aa_graphs_with_ln_weights) {
    

    protein = p;
    aa_graphs_with_ln_w = aa_graphs_with_ln_weights;
    aa_best_path_in_a_whole_codon = aa_best_in_codon;
    best_path_in_one_codon_unit = best_path_in_one_codon;

    seq_length = 3 * static_cast<IndexType>(aa_seq.size());
    next_pair.resize(NOTON);
    next_pair_set.resize(NOTON);
    get_next_pair(dfa);
    get_next_pair_set();

    prev_pair.resize(NOTON);
    prev_pair_set.resize(NOTON);
    get_prev_pair(dfa);
    get_prev_pair_set();

    next_list.resize(NOTON);
    prev_list.resize(NOTON);
    stacking_score.resize(6, vector<ScoreType>(6));
    bulge_score.resize(6, vector<vector<ScoreType>>(6, vector<ScoreType>(SINGLE_MAX_LEN+1)));

    preprocess(dfa);


    auto start_node = make_pair(0, 0);
    auto end_node = make_pair(seq_length, 0);


    if (epsilon < 1.0){
        dfa.random_init_with_mfe_path(init_solution, start_node, epsilon, rand_seed);
    }
    else {
        dfa.randomize_weights(rand_seed);
    }
    dfa.validate_and_project_probabilities();

    ScoreType best_obj = 0;
    auto best_nuc_seq = dfa.get_best_nuc_sequence(seq_length);

    ScoreType prev_coef(0), curr_coef(Nesterov_coef(0));
    for(int k=1; k<=num_epochs; ++k){
        prev_coef = curr_coef;
        curr_coef = Nesterov_coef(curr_coef);
        dfa.apply_momentum((prev_coef - 1) / curr_coef);
        dfa.validate_and_project_probabilities();

        dfa.zero_gradients();
        clear_parser_states();
        forward(dfa);
        backward(dfa);
        auto end_state = bestC[end_node];
        dfa.gradient_decent(learning_rate, end_state.inside);
        dfa.validate_and_project_probabilities();

        int width = int(log10(num_epochs)) + 1;
        cout<<"Iteration ["<< setw(width) << setfill(' ') << k <<"]: obj: "<<std::setprecision(4) << fixed<< - end_state.inside <<endl;

        if(end_state.inside > best_obj){
            best_obj = end_state.inside;
            best_nuc_seq = dfa.get_best_nuc_sequence(seq_length);
        }
    }

    cout << string(seq_length, '-') << endl;
    cout<<"Final mRNA sequence: "<<best_nuc_seq.first<<endl;
}

Optimizer::Optimizer(int beam_size_, int num_epochs_, double learning_rate_, double epsilon_, std::string init_solution_, bool is_verbose_, unsigned int rand_seed_)
        : beam_size(beam_size_), num_epochs(num_epochs_), learning_rate(learning_rate_), epsilon(epsilon_),
          init_solution(init_solution_), is_verbose(is_verbose_), rand_seed(rand_seed_) {
        func1();
}

ScoreType Optimizer::quickselect_partition(vector<ScoreInnerDate_t>& scores, ScoreType lower, ScoreType upper) {
    double pivot = scores[upper].newscore;
    while (lower < upper) {
        while (scores[lower].newscore < pivot)
            lower += 1;
        while (scores[upper].newscore > pivot)
            upper -= 1;
        if (abs(scores[lower].newscore - scores[upper].newscore) < EPS) // hzhang: do not use == for double type
            lower += 1;
        else if (lower < upper)
            swap(scores[lower], scores[upper]);
    }
    return upper;
}

ScoreType Optimizer::quickselect(vector<ScoreInnerDate_t>& scores, const ScoreType lower,
        const ScoreType upper, const IndexType k) {
    if (lower == upper)
        return scores[lower].newscore;
    ScoreType split = quickselect_partition(scores, lower, upper);
    IndexType length = split - lower + 1;
    if (length == k)
        return scores[split].newscore;
    else if (k < length)
        return quickselect(scores, lower, split-1, k);
    else
        return quickselect(scores, split+1, upper, k-length);
}

template <PhaseOption phase>
void Optimizer::beam_prune(DFA_t& dfa, BestX_t& bestX, const IndexType j) {
    if (phase != PhaseOption::Inside) return;
    vector<ScoreInnerDate_t> scores;
    IndexType size = 0;
    for (auto& j_node : {NodeType{j, 0}, NodeType{j, 1}}) {
        for (auto& i_node_elem : bestX[j_node]) {
            auto i_node = i_node_elem.first;
            for (auto& nuc_pair_elem : bestX[j_node][i_node]) {

                ScoreType newscore = 0.;


                if (phase == PhaseOption::Inside){
                    if (bestC[i_node].inside == util::value_min<ScoreType>()){
                        newscore = nuc_pair_elem.second.inside;
                    }
                    else{
                        newscore = nuc_pair_elem.second.inside + bestC[i_node].inside;
                    }
                }
                else{
                    throw std::logic_error("Function not yet implemented");
                }

                scores.push_back({newscore, j_node, i_node, nuc_pair_elem.first});
                size += 1;
            }
        }
    }
    if (size > beam_size) {
        ScoreType threshold = quickselect(scores, 0, size-1, size-beam_size+1);
        for (auto& scores_inner_data : scores) {
            if (scores_inner_data.newscore < threshold) {
                bestX[scores_inner_data.j_node][scores_inner_data.i_node].erase(scores_inner_data.nuc_pair);
            }
        }
    }
}

template <PhaseOption phase>
ScoreType Optimizer::beam_prune(DFA_t& dfa, BestM_t& bestX, const IndexType j) {
    if (phase != PhaseOption::Inside) return util::value_min<ScoreType>();
    reserved_scores.clear();
    IndexType size = 0;
    for (auto& j_node : {NodeType{j, 0}, NodeType{j, 1}}) {
        for (auto& i_node_elem : bestX[j_node]) {
            auto i_node = i_node_elem.first;

            ScoreType newscore = 0.;


            if (phase == PhaseOption::Inside){
                if (bestC[i_node].inside == util::value_min<ScoreType>()){
                    newscore = i_node_elem.second.inside;
                }
                else{
                    newscore = i_node_elem.second.inside + bestC[i_node].inside;
                }
            }
            else{
                throw std::logic_error("Function not yet implemented");
            }

            reserved_scores.push_back({newscore, j_node, i_node, 0});
            size += 1;
        }
    }

    ScoreType threshold = util::value_min<ScoreType>();

    if (size > beam_size) {
        threshold = quickselect(reserved_scores, 0, size-1, size-beam_size+1);
        for (auto& scores_inner_data : reserved_scores) {
            if (scores_inner_data.newscore < threshold) {
                bestX[scores_inner_data.j_node].erase(scores_inner_data.i_node);
            }
        }
    }

    return threshold;
}


}