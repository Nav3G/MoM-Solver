#include "quadrature.hpp"
#include <cmath>
#include <cassert>

namespace mom {
    GLRule gauss_legendre(int order) {
        switch (order) {
            case 3: {
                GLRule rule;
                rule.nodes   = { -std::sqrt(3.0/5.0), 0.0, +std::sqrt(3.0/5.0) };
                rule.weights = { 5.0/9.0, 8.0/9.0, 5.0/9.0 };
                return rule;
            }
            default:
                assert(false && "unsupported GL order");
                return {};   // unreachable, satisfies the compiler
        }
    }
}