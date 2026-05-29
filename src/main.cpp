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
#include <fstream>
#include <iomanip>
#include <complex>

int main() {
    using namespace mom;

    // 1. Configure dipole (half-wave).
    const double f = 300e6;
    const double lambda = constants::c / f;
    const double h = 0.235 * lambda;
    const double a = 0.005 * lambda;
    const int N = 50;

    Geometry dipole = make_straight_dipole(h, a, N);
    auto params = SolverParams{f, (N - 1) / 2};

    std::complex<double> V0 = {1.0, 0.0};

    // -------------------- v2 path --------------------
    const int feed_rooftop = 25; 
    const int gl_order = 3;

    auto Z_v2 = fill_impedance_matrix_v2(dipole, params, gl_order);
    auto v_v2 = delta_gap_v2(dipole, params, V0, feed_rooftop);
    Eigen::VectorXcd alpha_v2 = Z_v2.colPivHouseholderQr().solve(v_v2);

    auto Z_in_v2 = V0 / alpha_v2(feed_rooftop);
    std::cout << "Z_in (v2) = " << Z_in_v2.real() << " + j" << Z_in_v2.imag() << " ohms\n";

    write_matrix_csv("output/Z_v2.csv", Z_v2);
    write_current_csv("output/I_v2.csv", dipole, alpha_v2);

    return 0;
}