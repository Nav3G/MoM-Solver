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

    // Returns the self-term double integral:
    // integral over [0,L] x [0,L] of Lambda_obs(l) * Lambda_src(l') * G(R) dl' dl
    // using analytic inner integration (small-kr approximation) + GL outer.
    std::complex<double> singular_self_LambdaLambda(const Geometry& geom,
                                                    int s,                    // segment index (s_obs == s_src)
                                                    Ramp ramp_obs,
                                                    Ramp ramp_src,
                                                    double k,
                                                    int order) {
        double L_obs = geom.lengths[s];
        double a = geom.a;
        GLRule rule = gauss_legendre(order);

        std::complex<double> result = 0.0;
        for (int i = 0; i < order; ++i) {
            double x_obs = rule.nodes[i];
            double w_i = rule.weights[i];

            double l_obs = (L_obs / 2) * (1 + x_obs);
            double Lambda_obs = (ramp_obs == Ramp::Up) ? (l_obs / L_obs)
                                                      : (1.0 - l_obs / L_obs);

            std::complex<double> inner = (ramp_src == Ramp::Up) ? S_up(l_obs, L_obs, a, k)
                                                               : S_down(l_obs, L_obs, a, k);

            result += w_i * Lambda_obs * inner;
        }

        result *= (L_obs / 2.0);
        return result;
    }

    // Returns the self-term double integral:
    // integral over [0,L] x [0,L] of G(R) dl' dl
    // using analytic inner integration (small-kr approximation) + GL outer.
    std::complex<double> singular_self_G(const Geometry& geom,
                                         int s,
                                         double k,
                                         int order) {
        double L_obs = geom.lengths[s];
        double a = geom.a;
        GLRule rule = gauss_legendre(order);

        std::complex<double> result = 0.0;
        for (int i = 0; i < order; ++i) {
            double x_obs = rule.nodes[i];
            double w_i = rule.weights[i];

            double l_obs = (L_obs / 2) * (1 + x_obs);

            std::complex<double> inner = S_bare(l_obs, L_obs, a, k);

            result += w_i * inner;
        }

        result *= (L_obs / 2.0);
        return result;                                    
    }
}