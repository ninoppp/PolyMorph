#include "const.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

void validate_parameters() {
  assert(D_mu.size() == NUM_SPECIES && D_CV.size() == NUM_SPECIES);
  assert(p_mu.size() == NUM_SPECIES && p_CV.size() == NUM_SPECIES);
  assert(k_mu.size() == k_CV.size());
}

template <typename T>
std::string to_string(const std::vector<T>& v) {
    std::string s = "[";
    for (const T& x : v) {
        s += std::to_string(x) + ", ";
    }
    s.pop_back();
    s.pop_back();
    s += "]";
    return s;
}

void write_config() {
    std::ofstream config("simulation.cfg");
    config
        << "Date=" << __DATE__ << std::endl
        << "Time=" << __TIME__ << std::endl
        << "D0=" << D0 << std::endl
        << "k0=" << k0 << std::endl
        << "p0=" << p0 << std::endl
        << "D_mu=" << to_string(D_mu) << std::endl
        << "k_mu=" << to_string(k_mu) << std::endl
        << "p_mu=" << to_string(p_mu) << std::endl
        << "threshold_mu=" << to_string(threshold_mu) << std::endl 
        << "D_CV=" << to_string(D_CV) << std::endl
        << "k_CV=" << to_string(k_CV) << std::endl
        << "p_CV=" << to_string(p_CV) << std::endl
        << "threshold_CV=" << to_string(threshold_CV) << std::endl
        << "dx=" << dx << std::endl
        << "dt=" << dt << std::endl
        << "Nf=" << Nf << std::endl
        << "Ns=" << Ns << std::endl
        << "Nr=" << Nr << std::endl;
    config.close();
}
