#include "dipole.hpp"
#include "matrix_fill.hpp"
#include "excitation.hpp"
#include "io.hpp"
#include "constants.hpp"

#include <Eigen/Dense>
#include <iostream>
#include <complex>

int main() {
    using namespace mom;

    // 1. Configure dipole (half-wave).
    const double f = 300e6;
    const double lambda = constants::c / f;
    const double h = lambda / 4.0;
    const double a = lambda / 200.0;
    const int N = 41;
    Dipole d{h, a, f, N};
    
    // 2. Build Z.
    auto Z = fill_impedance_matrix(d);

    // 3. Build v.
    std::complex<double> V0 = {1.0, 0.0};
    auto v = magnetic_frill(d, V0, 2.3);

    // 4. Solve.
    Eigen::VectorXcd alpha = Z.colPivHouseholderQr().solve(v);

    auto Z_in = V0 / alpha(feed_segment_index(d));
    std::cout << "Z_in = " << Z_in.real() << " + j" << Z_in.imag() << " ohms\n";

    // 5. IO.
    write_current_csv("output/I.csv", d, alpha);
    write_matrix_csv("output/Z.csv", Z);

    return 0;
}