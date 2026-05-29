#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include "constants.hpp"
#include "geometry.hpp"

namespace mom {
    // Gauss-Legendre rule.
    struct GLRule {
        std::vector<double> nodes;    // on [-1, 1]
        std::vector<double> weights;
    };

    GLRule gauss_legendre(int order);   // returns the rule for n = 3, 5, 7

    // Integral driver.
    inline double static_inner_integral(double l, 
                                        double l_a, 
                                        double l_b, 
                                        double a) {
        return (std::asinh((l - l_a) / a) - std::asinh((l - l_b) / a)) / (4.0 * constants::pi);                                     
    }
    
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

    template <typename F_full>
    std::complex<double> singular_double_integral(const Geometry& geom,
                                                int s_obs, int s_src,
                                                int order,
                                                F_full&& f_full) {
        assert(s_obs == s_src && "singular_double_integral only handles self-term");                                            

        GLRule rule = gauss_legendre(order); 
        
        double L_obs = geom.lengths[s_obs];
        double L_src = geom.lengths[s_src];
        double a = geom.a;
        
        Eigen::Vector3d n_obs_left = geom.nodes[s_obs];
        Eigen::Vector3d n_obs_right = geom.nodes[s_obs + 1];
        Eigen::Vector3d n_src_left = geom.nodes[s_src];
        Eigen::Vector3d n_src_right = geom.nodes[s_src + 1];
        
        std::complex<double> smooth_part = 0.0;
        std::complex<double> singular_part = 0.0;
        for (int i = 0; i < order; ++i) {
            double x_obs = rule.nodes[i];
            Eigen::Vector3d r = 0.5 * (1.0 - x_obs) * n_obs_left
                                + 0.5 * (1.0 + x_obs) * n_obs_right;
            // Singular term (analytic).
            double l_i = (L_obs / 2.0) * (1.0 + x_obs);
            singular_part += rule.weights[i] * static_inner_integral(l_i, 0.0, L_src, a);

            for (int j = 0; j < order; ++j) {
                double x_src = rule.nodes[j];
                Eigen::Vector3d r_prime = 0.5 * (1.0 - x_src) * n_src_left
                                         + 0.5 * (1.0 + x_src) * n_src_right;
                // Smooth term.
                double R_eff = std::sqrt((r - r_prime).squaredNorm() + a * a);
                smooth_part += rule.weights[i] * rule.weights[j] 
                                               * (f_full(r, r_prime) - 1.0 / (4 * constants::pi * R_eff));
            }
        }

        smooth_part *= (L_obs / 2.0) * (L_src / 2.0);
        singular_part *= (L_obs / 2.0);
        return smooth_part + singular_part;
    }
}