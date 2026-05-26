#pragma once

#include "dipole.hpp"
#include <Eigen/Dense>
#include <complex>

namespace mom {
    // Integrate the kernel over one source segment using midpoint rule (v1).
    // Returns delta * K(z_m, z_n_center) for now.
    std::complex<double> integrate_segment(double z_m, 
                                           double z_n_center, 
                                           double k, 
                                           double a, 
                                           double delta);

    // Build the N x N complex impedance matrix Z by point-matching:
    // testing at segment centers z_m, integrating over source segments z_n.
    Eigen::MatrixXcd fill_impedance_matrix(const Dipole& dipole);
}