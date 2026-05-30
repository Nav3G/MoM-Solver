#include "io.hpp"
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <stdexcept>

namespace mom {
    void write_current_csv(const std::string& filename,
                           const Geometry& geom,
                           const Eigen::VectorXcd& alpha) {
        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Failed to open " + filename);
        }
        out << std::setprecision(15);
        out << "s,Re_I,Im_I,abs_I\n";

        // Cumulative arc length at each node.
        int N_nodes = geom.num_nodes();
        std::vector<double> arc(N_nodes, 0.0);
        for (int s = 0; s < N_nodes - 1; ++s) {
            arc[s + 1] = arc[s] + geom.lengths[s];
        }

        double total = arc[N_nodes - 1];

        for (int m = 0; m < num_rooftops(geom); ++m) {
            auto a = alpha(m);
            double s_anchor = arc[m + 1] - total / 2.0;
            out << s_anchor << ","
                << a.real() << ","
                << a.imag() << ","
                << std::abs(a) << "\n";
        }        
    }

    void write_matrix_csv(const std::string& filename,
                          const Eigen::MatrixXcd& M) {
        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Failed to open " + filename);
        }
        out << std::setprecision(15);
        out << "m,n,Re_Z,Im_Z\n";
        for (int m = 0; m < M.rows(); ++m) {
            for (int n = 0; n < M.cols(); ++n) {
                auto entry = M(m, n);
                out << m << "," 
                    << n << "," 
                    << entry.real() << "," 
                    << entry.imag() 
                    << "\n";
            }
        }
    }
}