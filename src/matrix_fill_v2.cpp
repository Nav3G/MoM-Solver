#include "matrix_fill_v2.hpp"
#include "quadrature.hpp"
#include "kernel.hpp"
#include "basis.hpp"
#include <complex>
#include <array>

namespace mom {
    struct SubPair { int s; Ramp sigma; };

    std::complex<double> integral_G(const Geometry& geom,
                                    const SolverParams& params,
                                    int s_obs, int s_src,
                                    int order) {
        double k = wavenumber(params);
        double a = geom.a;  
        
        auto kernel_lam = [k, a](const Eigen::Vector3d& r, 
                                 const Eigen::Vector3d& r_prime) {
            return green_3d(r, r_prime, k, a);
        };

        if (s_obs == s_src) { 
            return singular_double_integral(geom, s_obs, s_src, order, kernel_lam);
        } else {
            return double_integral(geom, s_obs, s_src, order, kernel_lam);
        }
    }

    std::complex<double> integral_LambdaLambdaG(const Geometry& geom,
                                                const SolverParams& params,
                                                int s_obs, int s_src,
                                                Ramp ramp_obs, Ramp ramp_src,
                                                int order) {
        double k = wavenumber(params);
        double a = geom.a; 
        double L_obs = geom.lengths[s_obs];
        double L_src = geom.lengths[s_src];
        Eigen::Vector3d n_obs_left = geom.nodes[s_obs];
        Eigen::Vector3d n_src_left = geom.nodes[s_src];

        auto integrand = [k, a, L_obs, L_src, n_obs_left, n_src_left, ramp_obs, ramp_src]
                        (const Eigen::Vector3d& r, const Eigen::Vector3d& r_prime) {
            double l_obs = (r - n_obs_left).norm();
            double l_src = (r_prime - n_src_left).norm();
            
            double Lambda_obs = (ramp_obs == Ramp::Up) ? (l_obs / L_obs)
                                           : (1.0 - l_obs / L_obs);
            double Lambda_src = (ramp_src == Ramp::Up) ? (l_src / L_src)
                                           : (1.0 - l_src / L_src);        
            
            return Lambda_obs * Lambda_src * green_3d(r, r_prime, k, a);
        };

        return double_integral(geom, s_obs, s_src, order, integrand);
    }

    Eigen::MatrixXcd fill_impedance_matrix_v2(const Geometry& geom,
                                              const SolverParams& params,
                                              int order) {
        int N_u = num_rooftops(geom);
        Eigen::MatrixXcd Z = Eigen::MatrixXcd::Zero(N_u, N_u);      
        
        double omega = 2.0 * constants::pi * params.frequency;
        std::complex<double> prefactor_A = constants::j * omega * constants::mu_0;
        std::complex<double> prefactor_phi = 1.0 / (constants::j * omega * constants::eps_0);

        for (int m = 0; m < N_u; ++m) {
            std::pair<int, int> segments_m = rooftop_segments(m);
            int s_m_left = segments_m.first;
            int s_m_right = segments_m.second;

            for (int n = 0; n < N_u; ++n) {
                std::pair<int, int> segments_n = rooftop_segments(n);
                int s_n_left = segments_n.first;
                int s_n_right = segments_n.second;

                std::complex<double> Z_A   = 0.0;
                std::complex<double> Z_phi = 0.0;

                std::array<SubPair, 2> m_sides = {{ {s_m_left, Ramp::Up}, {s_m_right, Ramp::Down} }};
                std::array<SubPair, 2> n_sides = {{ {s_n_left, Ramp::Up}, {s_n_right, Ramp::Down} }};

                for (auto [s, sigma_m] : m_sides) {
                    for (auto [s_prime, sigma_n] : n_sides) {
                        double t_dot_t = geom.tangents[s].dot(geom.tangents[s_prime]);
                        auto I_A   = integral_LambdaLambdaG(geom, params, s, s_prime, sigma_m, sigma_n, order);
                        Z_A += t_dot_t * I_A;

                        double dLm = rooftop_derivative(geom, m, s);
                        double dLn = rooftop_derivative(geom, n, s_prime);
                        auto I_phi = integral_G(geom, params, s, s_prime, order);
                        Z_phi += dLm * dLn * I_phi;
                    }
                }

                Z(m, n) = prefactor_A * Z_A + prefactor_phi * Z_phi;
            }
        }

        return Z;
    }
}