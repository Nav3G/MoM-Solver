#include "excitation_v2.hpp"
#include "basis.hpp"

namespace mom {
    Eigen::VectorXcd delta_gap_v2(const Geometry& geom,
                                  const SolverParams& params,
                                  std::complex<double> V0,
                                  int feed_rooftop) {
        int N_u = num_rooftops(geom);
        Eigen::VectorXcd v = Eigen::VectorXcd::Zero(N_u);
        v(feed_rooftop) = V0;
        
        return v;
    }
}