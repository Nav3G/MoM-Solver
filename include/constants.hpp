#pragma once
#include <complex>

namespace mom::constants {

    inline constexpr double pi = 3.14159265358979323846;
    inline constexpr double c = 299792458.0; 
    inline constexpr double mu_0 = 4.0 * pi * 1e-7;
    inline constexpr double eps_0 = 1.0 / (mu_0 * c * c);
    inline constexpr double eta_0 = mu_0 * c;
    inline const std::complex<double> j(0.0, 1.0);
    
}