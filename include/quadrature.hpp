#pragma once
#include "constants.hpp"
#include "geometry.hpp"
#include "basis.hpp"
#include <vector>
#include <complex>
#include <cmath>
#include <Eigen/Dense>

namespace mom {
    // Gauss-Legendre rule.
    struct GLRule {
        std::vector<double> nodes;    // on [-1, 1]
        std::vector<double> weights;
    };

    GLRule gauss_legendre(int order);   // returns the rule for n = 3, 5, 7

    // --------------------------- Integral driver ---------------------------------
    // Inner integral over [0, L] of e^{-jkr}/(4 pi r) dx', small-kr approximation.
    // Used for the Z_phi self-term.
    inline std::complex<double> S_bare(double x, double L, double a, double k) {
        double log = std::asinh(x / a) - std::asinh((x - L) / a);

        return (1.0 / (4.0 * constants::pi)) * (log - constants::j * k * L);
    }

    // Inner integral of (x'/L) * e^{-jkr}/(4 pi r) dx', small-kr approximation.
    // Used for Z_A self-term with an up-ramp source.
    inline std::complex<double> S_up(double x, double L, double a, double k) {
        double root_0 = std::sqrt(a * a + x * x);
        double root_L = std::sqrt(a * a + (x - L) * (x - L));
        
        double log = std::asinh(x / a) - std::asinh((x - L) / a);

        double static_part = (1.0 / (4.0 * constants::pi * L)) 
                            * (root_L - root_0 + x * log);
        std::complex<double> correction = constants::j * k * L / (8.0 * constants::pi);

        return static_part - correction;
    }

    // Inner integral of (1 - x'/L) * e^{-jkr}/(4 pi r) dx', via S_bare - S_up.
    inline std::complex<double> S_down(double x, double L, double a, double k) {
        return S_bare(x, L, a, k) - S_up(x, L, a, k);
    }

    // Returns the self-term double integral:
    // integral over [0,L] x [0,L] of Lambda_obs(l) * Lambda_src(l') * G(R) dl' dl
    // using analytic inner integration (small-kr approximation) + GL outer.
    std::complex<double> singular_self_LambdaLambda(const Geometry& geom,
                                                    int s,                    // segment index (s_obs == s_src)
                                                    Ramp ramp_obs,
                                                    Ramp ramp_src,
                                                    double k,
                                                    int order);
    
    // Returns the self-term double integral:
    // integral over [0,L] x [0,L] of G(R) dl' dl
    // using analytic inner integration (small-kr approximation) + GL outer.
    std::complex<double> singular_self_G(const Geometry& geom,
                                         int s,
                                         double k,
                                         int order);
    
    template <typename F>
    std::complex<double> double_integral(const Geometry& geom,
                                         int s_obs, int s_src,
                                         int order,
                                         F&& f) {
        GLRule rule = gauss_legendre(order); 
        
        double L_obs = geom.lengths[s_obs];
        double L_src = geom.lengths[s_src];

        Eigen::Vector3d n_obs_left = geom.nodes[s_obs];
        Eigen::Vector3d n_obs_right = geom.nodes[s_obs + 1];
        Eigen::Vector3d n_src_left = geom.nodes[s_src];
        Eigen::Vector3d n_src_right = geom.nodes[s_src + 1];

        std::complex<double> result = 0.0;
        for (int i = 0; i < order; ++i) {
            double x_obs = rule.nodes[i];
            Eigen::Vector3d r = 0.5 * (1.0 - x_obs) * n_obs_left
                                + 0.5 * (1.0 + x_obs) * n_obs_right;

            for (int j = 0; j < order; ++j) {
                double x_src = rule.nodes[j];
                Eigen::Vector3d r_prime = 0.5 * (1.0 - x_src) * n_src_left
                                         + 0.5 * (1.0 + x_src) * n_src_right;

                result += rule.weights[i] * rule.weights[j] * f(r, r_prime);
            }
        }

        result *= (L_obs / 2.0) * (L_src / 2.0);
        return result;
    }
}