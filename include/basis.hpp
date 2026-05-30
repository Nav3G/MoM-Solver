#pragma once
#include "geometry.hpp"

namespace mom {
    enum class Ramp { Up, Down };

    // Number of rooftops on this geometry.
    inline int num_rooftops(const Geometry& geom) {
        return geom.num_segments() - 1;
    }

    // The two segments touched by rooftop m. Returns {left_seg, right_seg}.
    inline std::pair<int, int> rooftop_segments(int m) {
        return {m, m + 1};
    }

    // Derivative of rooftop m on segment s: +1/L_s on left, -1/L_s on right, 0 otherwise.
    inline double rooftop_derivative(const Geometry& geom, int m, int s) {
        if (s == m) { return 1.0 / geom.lengths[s] ; }
        else if (s == m + 1) { return - 1.0 / geom.lengths[s]; } 
        else { return 0; }
    }

    inline int rooftop_at_node(int node_index) {
        return node_index - 1;
    }
}