#include "geometry.hpp"
#include "solver_params.hpp"
#include "basis.hpp"
#include "quadrature.hpp"
#include "kernel.hpp"
#include "matrix_fill_v2.hpp"
#include "excitation_v2.hpp"
#include "io.hpp"
#include "constants.hpp"

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>

int main() {
    using namespace mom;

    // Antenna parameters
    const double f      = 300e6;
    const double lambda = constants::c / f;
    const double h      = 0.235 * lambda;
    const double a      = 0.005 * lambda;
    const int    N_s    = 50;                          // even, so a node sits at the bend
    const double alpha  = constants::pi / 3.0;         // 45 degree half-angle

    // Build V-dipole
    Geometry v_dipole = make_v_dipole(h, a, N_s, alpha);

    // Geometry sanity prints
    std::cout << std::setprecision(6);
    std::cout << "--- V-dipole geometry, alpha = " << alpha << " rad ---\n";
    std::cout << "Lower tip:   " << v_dipole.nodes[0].transpose()       << "\n";
    std::cout << "Bend node:   " << v_dipole.nodes[N_s / 2].transpose() << "\n";
    std::cout << "Upper tip:   " << v_dipole.nodes[N_s].transpose()     << "\n";
    std::cout << "Tangent dot across bend: "
              << v_dipole.tangents[N_s / 2 - 1].dot(v_dipole.tangents[N_s / 2])
              << "  (expected cos(2*alpha) = " << std::cos(2.0 * alpha) << ")\n\n";

    // Solver setup
    auto params = SolverParams{f, N_s / 2};            // feed_node at the bend
    const int feed_rooftop = N_s / 2 - 1;              // rooftop anchored at the bend node
    const int gl_order = 3;

    std::complex<double> V0 = {1.0, 0.0};

    // Assemble + solve
    auto Z = fill_impedance_matrix_v2(v_dipole, params, gl_order);
    auto v = delta_gap_v2(v_dipole, params, V0, feed_rooftop);
    Eigen::VectorXcd alpha_vec = Z.colPivHouseholderQr().solve(v);

    auto Z_in = V0 / alpha_vec(feed_rooftop);
    std::cout << "Z_in = " << Z_in.real() << " + j" << Z_in.imag() << " ohms\n";

    // Outputs (note: write_current_csv writes z-coordinates only, which loses
    // x-information for bent wires. For now you'll see a folded projection.)
    write_matrix_csv("output/Z_v.csv",    Z);
    write_current_csv("output/I_v.csv",   v_dipole, alpha_vec);

    return 0;
}