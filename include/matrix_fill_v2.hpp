#pragma once
#include "geometry.hpp"
#include "solver_params.hpp"
#include <complex>
#include <Eigen/Dense>

namespace mom {
    enum class Ramp { Up, Down };

    std::complex<double> integral_G(const Geometry& geom,
                                    const SolverParams& params,
                                    int s_obs, int s_src,
                                    int order);
                                    
    std::complex<double> integral_LambdaLambdaG(const Geometry& geom,
                                                const SolverParams& params,
                                                int s_obs, int s_src,
                                                Ramp ramp_obs, Ramp ramp_src,
                                                int order);

    Eigen::MatrixXcd fill_impedance_matrix_v2(const Geometry& geom,
                                              const SolverParams& params,
                                              int order);                                       
}