#pragma once

#include "constants.hpp"
#include <complex>
#include <cmath>

namespace mom {
    // Thin-wire distance 
    inline double R(double z, double z_prime, double a) {
        return std::sqrt((z - z_prime) * (z - z_prime) + a * a);
    }

    // Scalar Green's funciton G(z, z')
    inline std::complex<double> green(double z, double z_prime, double k, double a) {
        auto r = R(z, z_prime, a);
        return std::exp(-constants::j * k * r) / (4.0 * constants::pi * r);
    }

    // Pocklington kernel
    inline std::complex<double> kernel(double z, 
                                       double z_prime, 
                                       double k, 
                                       double a, 
                                       double delta) {
        auto g = green(z, z_prime, k, a);
        auto g_right = green(z + delta, z_prime, k, a);
        auto g_left = green(z - delta, z_prime, k, a);
        return (g_right - 2.0 * g + g_left) / (delta * delta) + k * k * g;
    }
}