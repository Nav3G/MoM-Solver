#include "excitation.hpp"
#include <cmath>

namespace mom {
    Eigen::VectorXcd delta_gap(const Dipole& dipole, 
                               std::complex<double> V0) {
        Eigen::VectorXcd v = Eigen::VectorXcd::Zero(dipole.N);
        v(feed_segment_index(dipole)) = V0;
        return v;
    }

    Eigen::VectorXcd magnetic_frill(const Dipole& dipole,
                                    std::complex<double> V0,
                                    double b_over_a) {
        Eigen::VectorXcd v = Eigen::VectorXcd::Zero(dipole.N);

        double k = wavenumber(dipole);
        double delta = segment_width_delta(dipole);
        double a = dipole.a;
        double b = b_over_a * a;
        std::complex<double> norm = -V0 / (2 * std::log(b_over_a));

        for (int m = 0; m < dipole.N; ++m) {
            double z_m = segment_center(dipole, m);
            double R1 = std::sqrt(z_m * z_m + a * a);
            double R2 = std::sqrt(z_m * z_m + b * b);
            std::complex<double> E_inc = norm * (std::exp(-constants::j * k * R1) / R1 - 
                                                 std::exp(-constants::j * k * R2) / R2);
            v(m) = -delta * E_inc;
        }
        return v;
    }
}