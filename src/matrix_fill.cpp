#include "matrix_fill.hpp"
#include "kernel.hpp"

namespace mom {
    std::complex<double> integrate_segment(double z_m, double z_n_center, double k, double a, double delta) {
        return delta * kernel(z_m, z_n_center, k, a, delta);
    }

    Eigen::MatrixXcd fill_impedance_matrix(const Dipole& dipole) {
        int N = dipole.N; 
        Eigen::MatrixXcd Z = Eigen::MatrixXcd::Zero(N, N);

        const double k = wavenumber(dipole);
        const double delta = segment_width_delta(dipole);
        for (int m = 0; m < N; ++m) {
            const double z_m = segment_center(dipole, m);
            for (int n = 0; n < N; ++n) {
                Z(m, n) = integrate_segment(z_m, segment_center(dipole, n), k, dipole.a, delta);
            }
        }

        return Z;
    }
}
