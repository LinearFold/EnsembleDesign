#pragma once

#ifndef common_h
#define common_h

#include <utility>
#include <functional>
#include <array>
#include <string>
#include <set>
#include <map>
#include <tuple>
#include <exception>
#include <list>
#include "base.h"


#define kT 61.63207755
#define NEG_INF -2e20

namespace EnsembleDesign {


using SizeType               = size_t;
using ScoreType              = double; //int32_t;
using IndexType              = int32_t; //if less than 10000, only int16_t is needed here
using NucType                = int8_t;
using NumType                = int32_t;
using NucPairType            = int8_t;
using FinalScoreType         = double;

inline void Fast_LogPlusEquals (ScoreType &x, ScoreType y);
// Define the Parameter class
class Parameter {
public:
    ScoreType weight, prev_weight;
    ScoreType gradient;
    ScoreType cai_score;

    Parameter(ScoreType w = util::value_min<ScoreType>(), ScoreType g = util::value_min<ScoreType>(), ScoreType c = util::value_min<ScoreType>())
            : weight(w), prev_weight(w), gradient(g), cai_score(c) {};

    // Reset the gradient to 0
    void zero_gradient() {
        gradient = util::value_min<ScoreType>();
    }

    // Update gradient using log-plus
    void update_gradient(ScoreType delta) {
        Fast_LogPlusEquals(gradient, delta);
    }
};

using Parameters             = std::vector<Parameter*>;
using NodeType               = std::pair<IndexType, NumType>;
using EdgeType               = std::tuple<NodeType, NodeType, NucType, Parameter*>;
using IndexWType             = std::pair<IndexType, ScoreType>;
using NodeNucType            = std::pair<NodeType, NucType>;
using NodeNucWType           = std::tuple<NodeType, NucType, ScoreType>;
using NodeNucNodeType        = std::tuple<NodeType, NucType, NodeType>;
using NextPairType           = std::tuple<NodeType, NucType, Parameter*>;

enum class PhaseOption : std::uint8_t {
    Inside = 0,
    Outside
};

template <typename ScoreType>
struct State {
    ScoreType inside = util::value_min<ScoreType>();
    ScoreType outside = util::value_min<ScoreType>();
    NodeType pre_node = NodeType();
};

struct ScoreInnerDate {
    ScoreType newscore;
    NodeType j_node;
    NodeType i_node;
    int nuc_pair;
};


inline ScoreType LogExpPlusOne(ScoreType x) {
    return std::log1p(std::exp(x)); // log1p(exp(x)) is more accurate for small x
}

inline void LogPlusEquals(ScoreType &x, ScoreType y) {
    if (x < y) std::swap(x, y);
    if (y > -std::numeric_limits<ScoreType>::infinity()) {
        x = std::log1p(std::exp(y - x)) + x; // log1p(exp(y-x)) + x is numerically stable
    }
}

inline ScoreType Fast_LogExpPlusOne(ScoreType x){

    // Bounds for tolerance of 7.05e-06: (0, 11.8625)
    // Approximating interval: (0, 0.661537) --> ((T(-0.0065591595)*x+T(0.1276442762))*x+T(0.4996554598))*x+T(0.6931542306);
    // Approximating interval: (0.661537, 1.63202) --> ((T(-0.0155157557)*x+T(0.1446775699))*x+T(0.4882939746))*x+T(0.6958092989);
    // Approximating interval: (1.63202, 2.49126) --> ((T(-0.0128909247)*x+T(0.1301028251))*x+T(0.5150398748))*x+T(0.6795585882);
    // Approximating interval: (2.49126, 3.37925) --> ((T(-0.0072142647)*x+T(0.0877540853))*x+T(0.6208708362))*x+T(0.5909675829);
    // Approximating interval: (3.37925, 4.42617) --> ((T(-0.0031455354)*x+T(0.0467229449))*x+T(0.7592532310))*x+T(0.4348794399);
    // Approximating interval: (4.42617, 5.78907) --> ((T(-0.0010110698)*x+T(0.0185943421))*x+T(0.8831730747))*x+T(0.2523695427);
    // Approximating interval: (5.78907, 7.81627) --> ((T(-0.0001962780)*x+T(0.0046084408))*x+T(0.9634431978))*x+T(0.0983148903);
    // Approximating interval: (7.81627, 11.8625) --> ((T(-0.0000113994)*x+T(0.0003734731))*x+T(0.9959107193))*x+T(0.0149855051);
    // 8 polynomials needed.

    assert(ScoreType(0.0000000000) <= x && x <= ScoreType(11.8624794162) && "Argument out-of-range.");
    if (x < ScoreType(3.3792499610))
    {
        if (x < ScoreType(1.6320158198))
        {
            if (x < ScoreType(0.6615367791))
                return ((ScoreType(-0.0065591595)*x+ScoreType(0.1276442762))*x+ScoreType(0.4996554598))*x+ScoreType(0.6931542306);
            return ((ScoreType(-0.0155157557)*x+ScoreType(0.1446775699))*x+ScoreType(0.4882939746))*x+ScoreType(0.6958092989);
        }
        if (x < ScoreType(2.4912588184))
            return ((ScoreType(-0.0128909247)*x+ScoreType(0.1301028251))*x+ScoreType(0.5150398748))*x+ScoreType(0.6795585882);
        return ((ScoreType(-0.0072142647)*x+ScoreType(0.0877540853))*x+ScoreType(0.6208708362))*x+ScoreType(0.5909675829);
    }
    if (x < ScoreType(5.7890710412))
    {
        if (x < ScoreType(4.4261691294))
            return ((ScoreType(-0.0031455354)*x+ScoreType(0.0467229449))*x+ScoreType(0.7592532310))*x+ScoreType(0.4348794399);
        return ((ScoreType(-0.0010110698)*x+ScoreType(0.0185943421))*x+ScoreType(0.8831730747))*x+ScoreType(0.2523695427);
    }
    if (x < ScoreType(7.8162726752))
        return ((ScoreType(-0.0001962780)*x+ScoreType(0.0046084408))*x+ScoreType(0.9634431978))*x+ScoreType(0.0983148903);
    return ((ScoreType(-0.0000113994)*x+ScoreType(0.0003734731))*x+ScoreType(0.9959107193))*x+ScoreType(0.0149855051);
}

inline void Fast_LogPlusEquals (ScoreType &x, ScoreType y)
{
    if (x < y) std::swap (x, y);
    if (y > ScoreType(NEG_INF/2) && x-y < ScoreType(11.8624794162))
        x = Fast_LogExpPlusOne(x-y) + y;
}

inline ScoreType Fast_Exp(ScoreType x)
{
    // Bounds for tolerance of 4.96e-05: (-9.91152, 0)
    // Approximating interval: (-9.91152, -5.86228) --> ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014);
    // Approximating interval: (-5.86228, -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740);
    // Approximating interval: (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221);
    // Approximating interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
    // Approximating interval: (-1.48054, -0.672505) --> ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728);
    // Approximating interval: (-0.672505, -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077);
    // 6 polynomials needed.

    if (x < ScoreType(-2.4915033807))
    {
        if (x < ScoreType(-5.8622823336))
        {
            if (x < ScoreType(-9.91152))
                return ScoreType(0);
            return ((ScoreType(0.0000803850)*x+ScoreType(0.0021627428))*x+ScoreType(0.0194708555))*x+ScoreType(0.0588080014);
        }
        if (x < ScoreType(-3.8396630909))
            return ((ScoreType(0.0013889414)*x+ScoreType(0.0244676474))*x+ScoreType(0.1471290604))*x+ScoreType(0.3042757740);
        return ((ScoreType(0.0072335607)*x+ScoreType(0.0906002677))*x+ScoreType(0.3983111356))*x+ScoreType(0.6245959221);
    }
    if (x < ScoreType(-0.6725053211))
    {
        if (x < ScoreType(-1.4805375919))
            return ((ScoreType(0.0232410351)*x+ScoreType(0.2085645908))*x+ScoreType(0.6906367911))*x+ScoreType(0.8682322329);
        return ((ScoreType(0.0573782771)*x+ScoreType(0.3580258429))*x+ScoreType(0.9121133217))*x+ScoreType(0.9793091728);
    }
    if (x < ScoreType(0))
        return ((ScoreType(0.1199175927)*x+ScoreType(0.4815668234))*x+ScoreType(0.9975991939))*x+ScoreType(0.9999505077);
    return (x > ScoreType(46.052) ? ScoreType(1e20) : expf(x));
}

std::vector<ScoreType> ProjectToSimplex(const std::vector<ScoreType>& v) {
    int n_dims = v.size();
    std::vector<ScoreType> u = v;

    // Sort u in descending order
    std::sort(u.begin(), u.end(), std::greater<ScoreType>());

    std::vector<ScoreType> lambdas(n_dims);
    ScoreType cum_sum = 0.0;
    for (int i = 0; i < n_dims; ++i) {
        cum_sum += u[i];
        lambdas[i] = 1.0 - cum_sum;
    }

    std::vector<int> indexes(n_dims);
    std::iota(indexes.begin(), indexes.end(), 1);  // Fill with 1, 2, ..., n_dims

    int rho = 0;
    for (int i = 0; i < n_dims; ++i) {
        if (u[i] + lambdas[i] / indexes[i] > 0) {
            rho = indexes[i];
        }
    }

    std::vector<ScoreType> v_proj(n_dims);
    for (int i = 0; i < n_dims; ++i) {
        v_proj[i] = std::max(v[i] + lambdas[rho - 1] / rho, 1e-32);
    }

    return v_proj;
}

}

#endif /* common_h */
