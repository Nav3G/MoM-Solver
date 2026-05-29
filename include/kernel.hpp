#pragma once

#include "constants.hpp"
#include <Eigen/Dense>
#include <complex>
#include <cmath>

namespace mom {
    // Thin-wire distance (3D)
    inline double R_3d(const Eigen::Vector3d& r, const Eigen::Vector3d& r_prime, double a) {
        return std::sqrt((r - r_prime).squaredNorm() + a * a);
    }

    // Vector Green's funciton G(z, z')
    inline std::complex<double> green_3d(const Eigen::Vector3d& r,
                                         const Eigen::Vector3d& r_prime,
                                         double k, double a) {
        auto R_eff = R_3d(r, r_prime, a);
        return std::exp(-constants::j * k * R_eff) / (4.0 * constants::pi * R_eff);                                    
    }
}