#pragma once

#include "geometry.hpp"
#include "basis.hpp"
#include <Eigen/Dense>
#include <string>

namespace mom {
    void write_current_csv(const std::string& filename,
                           const Geometry& geom,
                           const Eigen::VectorXcd& alpha);

    void write_matrix_csv(const std::string& filename,
                          const Eigen::MatrixXcd& M);
}