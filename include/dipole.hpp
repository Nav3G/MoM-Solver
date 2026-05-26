#pragma once

#include "constants.hpp"
#include <cassert>

namespace mom {
    struct Dipole {
        double h;       // half-length [m]
        double a;       // radius [m]
        double f;       // frequency [Hz]
        int N;          // number of segments (odd)
    };

    inline double wavenumber(const Dipole& dipole) {
        return 2.0 * constants::pi * dipole.f / constants::c;
    }

    inline double segment_width_delta(const Dipole& dipole) {
        return 2.0 * dipole.h / static_cast<double>(dipole.N);
    }

    inline double segment_center(const Dipole& dipole, int n) {
        assert(n >= 0 && n < dipole.N);

        return -dipole.h + (n + 0.5) * segment_width_delta(dipole);
    }

    inline int feed_segment_index(const Dipole& dipole) {
        assert(dipole.N % 2 == 1);

        return (dipole.N - 1) / 2;
    }
}