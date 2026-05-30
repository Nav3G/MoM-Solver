#include "matrix_fill_v2.hpp"
#include "quadrature.hpp"
#include "kernel.hpp"
#include "basis.hpp"
#include <complex>
#include <array>

namespace mom {
    // Bundles "which segment + which ramp direction" into one value.
    // Used to represent one side (left or right) of a rooftop.
    struct SubPair { int s; Ramp sigma; };

    // ---------------------------------------------------------------------
    // Integrates the bare Green's function G(r, r') over segments s_obs x s_src.
    // Routes to the singularity-subtracted quadrature when s_obs == s_src.
    // ---------------------------------------------------------------------
    std::complex<double> integral_G(const Geometry& geom,
                                    const SolverParams& params,
                                    int s_obs, int s_src,
                                    int order) {
    double k = wavenumber(params);

    if (s_obs == s_src) {
        return singular_self_G(geom, s_obs, k, order);
    } else {
        double a = geom.a;
        auto kernel_lam = [k, a](const Eigen::Vector3d& r,
                                 const Eigen::Vector3d& r_prime) {
            return green_3d(r, r_prime, k, a);
        };
        
        return double_integral(geom, s_obs, s_src, order, kernel_lam);
    }
}

    // ---------------------------------------------------------------------
    // Integrates Lambda_obs(l) * Lambda_src(l') * G(r, r') over s_obs x s_src.
    // Lambda is the scalar ramp (Up: l/L, Down: 1 - l/L) on each segment.
    // ---------------------------------------------------------------------
    std::complex<double> integral_LambdaLambdaG(const Geometry& geom,
                                                const SolverParams& params,
                                                int s_obs, int s_src,
                                                Ramp ramp_obs, Ramp ramp_src,
                                                int order) {
        double k = wavenumber(params);

        if (s_obs == s_src) { 
            return singular_self_LambdaLambda(geom, s_obs, ramp_obs, ramp_src, k, order);
        } else {
            double a = geom.a;
            double L_obs = geom.lengths[s_obs];
            double L_src = geom.lengths[s_src]; 
            Eigen::Vector3d n_obs_left = geom.nodes[s_obs];
            Eigen::Vector3d n_src_left = geom.nodes[s_src];     

            // Capture everything the integrand needs to compute Lambda on the fly:
            // segment length and left endpoint give us l = |r - n_left|, and
            // ramp direction picks the formula.
            auto integrand = [k, a, L_obs, L_src, n_obs_left, n_src_left, ramp_obs, ramp_src]
                            (const Eigen::Vector3d& r, const Eigen::Vector3d& r_prime) {
                double l_obs = (r - n_obs_left).norm();
                double l_src = (r_prime - n_src_left).norm();
                
                // Ternary: condition ? value_if_true : value_if_false.
                double Lambda_obs = (ramp_obs == Ramp::Up) ? (l_obs / L_obs)
                                                        : (1.0 - l_obs / L_obs);
                double Lambda_src = (ramp_src == Ramp::Up) ? (l_src / L_src)
                                                        : (1.0 - l_src / L_src);        
                
                return Lambda_obs * Lambda_src * green_3d(r, r_prime, k, a);
            };

            return double_integral(geom, s_obs, s_src, order, integrand);
        }                                            
    }

    // ---------------------------------------------------------------------
    // Builds the N_u x N_u impedance matrix Z.
    // For each test rooftop m and source rooftop n, we sum over the 4 sub-segment
    // pairs (m_left/m_right) x (n_left/n_right) and accumulate Z_A and Z_phi.
    // ---------------------------------------------------------------------
    Eigen::MatrixXcd fill_impedance_matrix_v2(const Geometry& geom,
                                              const SolverParams& params,
                                              int order) {
        int N_u = num_rooftops(geom);
        Eigen::MatrixXcd Z = Eigen::MatrixXcd::Zero(N_u, N_u);      
        
        double omega = 2.0 * constants::pi * params.frequency;
        std::complex<double> prefactor_A = constants::j * omega * constants::mu_0;
        std::complex<double> prefactor_phi = 1.0 / (constants::j * omega * constants::eps_0);

        // Outer two loops: pick test rooftop m and source rooftop n.
        for (int m = 0; m < N_u; ++m) {
            // rooftop_segments(m) returns std::pair<int,int> = {left_seg, right_seg}.
            auto [s_m_left, s_m_right] = rooftop_segments(m);

            for (int n = 0; n < N_u; ++n) {
                std::pair<int, int> segments_n = rooftop_segments(n);
                int s_n_left  = segments_n.first;
                int s_n_right = segments_n.second;

                std::complex<double> Z_A   = 0.0;
                std::complex<double> Z_phi = 0.0;

                // -- std::array with aggregate initialization. --
                // std::array<T, N> is a fixed-size array known at compile time.
                // The OUTER {{ }} initializes the std::array (which internally
                // wraps a C-array, hence the doubled braces). Each INNER {...}
                // builds a SubPair by aggregate-initializing its two fields.
                // So m_sides is: [{s_m_left, Up}, {s_m_right, Down}].
                std::array<SubPair, 2> m_sides = {{ {s_m_left,  Ramp::Up},
                                                    {s_m_right, Ramp::Down} }};
                std::array<SubPair, 2> n_sides = {{ {s_n_left,  Ramp::Up},
                                                    {s_n_right, Ramp::Down} }};

                // -- range-for with structured bindings. --
                // 'for (auto X : container)' iterates over container, binding X
                // to each element in turn. Writing 'auto [s, sigma_m]' unpacks
                // the SubPair's two fields into local variables on each iteration.
                // Equivalent longhand:
                //     for (const auto& sub : m_sides) {
                //         int s = sub.s;
                //         Ramp sigma_m = sub.sigma;
                //         ...
                //     }
                // The two nested loops together iterate over all 4 combinations
                // of (m-side) x (n-side) -- the 2x2 sub-segment sum.
                for (auto [s, sigma_m] : m_sides) {
                    for (auto [s_prime, sigma_n] : n_sides) {
                        // --- Vector-potential contribution ---
                        // The rooftop's vector form is Lambda * t_hat, so the dot
                        // product of the two vector basis functions contributes
                        // (t_s . t_s') as a scalar factor.
                        double t_dot_t = geom.tangents[s].dot(geom.tangents[s_prime]);
                        auto I_A = integral_LambdaLambdaG(geom, params,
                                                          s, s_prime,
                                                          sigma_m, sigma_n, order);
                        Z_A += t_dot_t * I_A;

                        // --- Scalar-potential contribution ---
                        // The divergence of rooftop m on segment s is +/- 1/L_s
                        // (sign depends on whether s is the left or right segment
                        // of rooftop m). rooftop_derivative returns that value.
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