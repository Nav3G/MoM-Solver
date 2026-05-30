#pragma once

#include <Eigen/Dense>
#include <vector>

namespace mom {
    struct Geometry {
        std::vector<Eigen::Vector3d> nodes;      // size N_s + 1
        std::vector<Eigen::Vector3d> positions;  // size N_s (segment centers)
        std::vector<Eigen::Vector3d> tangents;   // size N_s (unit)
        std::vector<double> lengths;             // size N_s
        double a;                                // wire radius

        int num_segments() const { return static_cast<int>(positions.size()); }
        int num_nodes() const { return static_cast<int>(nodes.size()); }
    };

    Geometry make_straight_dipole(double h, double a, int N_s);

    Geometry make_v_dipole(double h, double a, int N_s, double half_angle_rad);
}