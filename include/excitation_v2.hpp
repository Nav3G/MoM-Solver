#pragma once
#include "geometry.hpp"
#include "solver_params.hpp"
#include <Eigen/Dense>
#include <complex>

namespace mom {
    // Delta-gap voltage source at the rooftop indexed by feed_rooftop.
    // Returns v of length N_u with v[feed_rooftop] = V0, all others zero.
    Eigen::VectorXcd delta_gap_v2(const Geometry& geom,
                                  const SolverParams& params,
                                  std::complex<double> V0,
                                  int feed_rooftop);
}