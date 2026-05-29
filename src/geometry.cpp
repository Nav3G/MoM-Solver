#include "geometry.hpp"
#include <cassert>

namespace mom {
    Geometry make_straight_dipole(double h, double a, int N_s) {
        assert(N_s % 2 == 0);
        assert(h > 0);
        assert(a > 0);

        Geometry g;
        g.a = a;

        double delta = 2.0 * h / N_s;

        // Node fill.
        g.nodes.reserve(N_s + 1);
        for (int i = 0; i <= N_s; ++i) {
            double z_i = -h + i * delta;
            g.nodes.emplace_back(0.0, 0.0, z_i);
        }

        // Per-segment data, derived from nodes.
        g.positions.reserve(N_s);
        g.tangents.reserve(N_s);
        g.lengths.reserve(N_s);
        for (int s = 0; s < N_s; ++s) {
            Eigen::Vector3d edge = g.nodes[s + 1] - g.nodes[s];
            double L = edge.norm();
            g.lengths.emplace_back(L);

            Eigen::Vector3d mid = 0.5 * (g.nodes[s] + g.nodes[s + 1]);
            g.positions.emplace_back(mid);

            g.tangents.emplace_back(edge / L);
        }

        return g;
    }
}