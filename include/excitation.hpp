#pragma once

#include "dipole.hpp"
#include <Eigen/Dense>
#include <complex>

namespace mom {
    // Delta-gap voltage source at center segment. Returns the RHS vector v of length N, all zero except
    // v(feed_idx) = V0.
    Eigen::VectorXcd delta_gap(const Dipole& dipole, 
                               std::complex<double> V0);

    // Magnetic frill voltage source at the center, modeling a coaxial feed
    // of inner radius a (=wire radius) and outer radius b.
    // Standard choice: b/a = 2.3 (corresponds to 50-ohm coax).
    Eigen::VectorXcd magnetic_frill(const Dipole& dipole,
                                    std::complex<double> V0,
                                    double b_over_a);
}