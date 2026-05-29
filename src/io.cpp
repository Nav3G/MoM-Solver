#include "io.hpp"
#include <fstream>
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
        out << "z,Re_I,Im_I,abs_I\n";
        for (int m = 0; m < num_rooftops(geom); ++m) {
            auto a = alpha(m);
            double z_anchor = geom.nodes[m + 1].z();
            out << z_anchor << ","
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