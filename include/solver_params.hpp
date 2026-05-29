#pragma once
#include "constants.hpp" 

namespace mom {
    struct SolverParams {
        double frequency;
        int feed_node;
    };

    inline double wavenumber(const SolverParams& p) {
        return 2 * constants::pi * p.frequency / constants::c;
    }
}