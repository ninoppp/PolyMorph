#ifndef UTILS_H
#define UTILS_H

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <random>

#include "const.h"

// easier to read than a pair
struct Index {
  int i; 
  int j;
  Index(int i, int j): i(i), j(j) {}
};

void validate_parameters() {
    assert(D_mu.size() == NUM_SPECIES && D_CV.size() == NUM_SPECIES);
    assert(p_mu.size() == NUM_SPECIES && p_CV.size() == NUM_SPECIES);
    assert(k_mu.size() == k_CV.size());
    assert(threshold_mu.size() == threshold_CV.size());
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

// returns the means for a vector of lognormal distributions
std::vector<double> get_means(const std::vector<std::lognormal_distribution<>>& dists) {
    std::vector<double> means;
    for (const auto& dist : dists) {
        double mu_prime = dist.m(); // Log-mean of the underlying normal distribution
        double sigma = dist.s(); // Standard deviation of the underlying normal distribution
        double mean = std::exp(mu_prime + 0.5 * sigma * sigma);
        means.push_back(mean);
    }
    return means;
}

// creates a vector of lognormal distributions from vectors of means and CVs
std::vector<std::lognormal_distribution<>> create_lognormal(const std::vector<double>& mu, const std::vector<double>& CV) {
    std::vector<std::lognormal_distribution<>> dists;
    for (int i = 0; i < mu.size(); i++) {
        double sigma = std::sqrt(std::log(1 + CV[i]*CV[i]));  // Standard deviation of the underlying normal
        double mu_prime = std::log(mu[i]) - 0.5 * sigma * sigma;  // Adjust mean of the underlying normal
        dists.push_back(std::lognormal_distribution<>(mu_prime, sigma));
    }   
    return dists;
}

// returns a vector of samples from a vector of lognormal distributions
std::vector<double> sample(std::vector<std::lognormal_distribution<>>& dists, std::mt19937& rng) {
    std::vector<double> samples;
    auto mean = get_means(dists);
    for (int i = 0; i < dists.size(); i++) {
        double x = dists[i](rng);
        int max_tries = 50;
        while (x > cutoff_factor*mean[i]) {
            x = dists[i](rng);
            if (max_tries-- == 0) {
                x = mean[i];
                break;
            } 
        }
        samples.push_back(x);
    }
    return samples;
}

// writes a nice welcome message to console to lift my mood
void welcome() {
    std::cout << "--------------------------" << std::endl
            << "|  Welcome to PolyMorph  |" << std::endl
            << "--------------------------" << std::endl;
}

void write_config() {
    std::ofstream config("simulation.cfg");
    config
        << "Date=" << __DATE__ << std::endl
        << "Time=" << __TIME__ << std::endl
        << "NUM_SPECIES=" << NUM_SPECIES << std::endl
        << "NUM_KIN=" << NUM_KIN << std::endl
        << "ADVECTION_DILUTION=" << ADVECTION_DILUTION << std::endl
        << "D0=" << to_string(D0) << std::endl
        << "k0=" << to_string(k0) << std::endl
        << "p0=" << to_string(p0) << std::endl
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
        << "Nr=" << Nr << std::endl
        << "cutoff_factor=" << cutoff_factor << std::endl
        << "Output::u=" << Output::u << std::endl
        << "Output::D=" << Output::D << std::endl
        << "Output::p=" << Output::p << std::endl
        << "Output::k=" << Output::k << std::endl
        << "Output::parent_idx=" << Output::parent_idx << std::endl
        << "Output::threshold=" << Output::threshold << std::endl
        << "Output::flag=" << Output::flag << std::endl
        << "Output::velocity=" << Output::velocity << std::endl;
        // expand if more parameters become relevant
    config.close();
}

#endif