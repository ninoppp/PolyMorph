#ifndef UTILS_H
#define UTILS_H

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include <sys/time.h>
#include <omp.h>

#include "const.h"

// easier to read than a pair
struct Index {
  int i; 
  int j;
  Index(int i, int j): i(i), j(j) {}
};

std::vector<Index> neighbors(Index idx) {
    return {Index(idx.i-1, idx.j), Index(idx.i+1, idx.j), Index(idx.i, idx.j-1), Index(idx.i, idx.j+1)};
}

double walltime(){
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        std::cout << "Unexpected error in get_wall_time" << std::endl;
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void validate_parameters() {
    assert(D_mu.size() == NUM_SPECIES && D_CV.size() == NUM_SPECIES);
    assert(p_mu.size() == NUM_SPECIES && p_CV.size() == NUM_SPECIES);
    assert(k_mu.size() == NUM_KIN && k_CV.size() == NUM_KIN);
    assert(threshold_mu.size() == threshold_CV.size());
    assert(NUM_SPECIES >= threshold_mu.size());
    assert(D0.size() == NUM_SPECIES && p0.size() == NUM_SPECIES);
    assert(k0.size() == NUM_KIN);
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

double mean(const std::vector<double>& v) {
    double sum = 0;
    for (const double& x : v) {
        sum += x;
    }
    return sum / v.size();
}

double variance(const std::vector<double>& v) {
    double m = mean(v);
    double sum = 0;
    for (const double& x : v) {
        sum += (x - m) * (x - m);
    }
    return sum / v.size();
}

double stddev(const std::vector<double>& v) {
    return std::sqrt(variance(v));
}

// get the mean of a lognormal distribution
double get_mean(const std::lognormal_distribution<>& dist) {
    double mu_prime = dist.m(); // Log-mean of the underlying normal distribution
    double sigma = dist.s(); // Standard deviation of the underlying normal distribution
    return std::exp(mu_prime + 0.5 * sigma * sigma);
}

// get the means of a vector of lognormal distributions
std::vector<double> get_means(const std::vector<std::lognormal_distribution<>>& dists) {
    std::vector<double> means;
    for (const auto& dist : dists) {
        means.push_back(get_mean(dist));
    }
    return means;
}

// creates a vector of lognormal distributions from vectors of means and CVs
std::vector<std::lognormal_distribution<>> create_lognormal(const std::vector<double>& mu, const std::vector<double>& CV) {
    assert(mu.size() == CV.size());
    std::vector<std::lognormal_distribution<>> dists;
    for (int i = 0; i < mu.size(); i++) {
        double sigma = std::sqrt(std::log(1 + CV[i]*CV[i]));  // Standard deviation of the underlying normal
        double mu_prime = std::log(mu[i]) - 0.5 * sigma * sigma;  // Adjust mean of the underlying normal
        dists.push_back(std::lognormal_distribution<>(mu_prime, sigma));
    }   
    return dists;
}

// sample from a lognormal distribution, with a cutoff factor
double sample(std::lognormal_distribution<>& dist, std::mt19937& rng) {
    double mean = get_mean(dist);
    double x = dist(rng);   
    int max_tries = 50;
    while (x > cutoff_factor*mean) {
        x = dist(rng);
        if (max_tries-- == 0) {
            x = mean;
            break;
        } 
    }
    return x;
}

// returns a vector of samples from a vector of lognormal distributions
std::vector<double> sample(std::vector<std::lognormal_distribution<>>& dists, std::mt19937& rng) {
    std::vector<double> samples;
    auto mean = get_means(dists);
    for (int i = 0; i < dists.size(); i++) {
        double x = sample(dists[i], rng);
        samples.push_back(x);
    }
    return samples;
}

// writes a nice welcome message to console to lift my mood
void welcome() {
    std::cout << "--------------------------" << std::endl
            << "|  Welcome to PolyMorph  |" << std::endl
            << "--------------------------" << std::endl;
    std::cout << "simulation started at " << __DATE__ << ", " << __TIME__ << std::endl;
    std::cout << "max threads = " << omp_get_max_threads() << std::endl;
}

// saving simluations parameters
void write_config() {
    std::ofstream config("simulation.cfg");
    config
        << "Date=" << __DATE__ << std::endl
        << "Time=" << __TIME__ << std::endl
        << "h=" << h << std::endl
        << "lmin=" << lmin << std::endl
        << "lmax=" << lmax << std::endl
        << "Q=" << Q << std::endl
        << "alpha_mu=" << alpha_mu << std::endl
        << "alpha_CV=" << alpha_CV << std::endl
        << "beta=" << beta << std::endl
        << "Amin=" << Amin << std::endl
        << "Amax_mu=" << Amax_mu << std::endl
        << "Amax_CV=" << Amax_CV << std::endl
        << "gam=" << gam << std::endl
        << "ka=" << ka << std::endl
        << "kl=" << kl << std::endl
        << "kb=" << kb << std::endl
        << "kr=" << kr << std::endl
        << "kh=" << kh << std::endl
        << "sh=" << sh << std::endl
        << "ss=" << ss << std::endl
        << "theta=" << theta << std::endl
        << "mu=" << mu << std::endl
        << "rho=" << rho << std::endl
        << "g=" << g << std::endl
        << "gl=" << gl << std::endl
        << "cv=" << cv << std::endl
        << "cd=" << cd << std::endl
        << "cc=" << cc << std::endl
        << "dt=" << dt << std::endl
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
        << "chemotaixs_strength=" << to_string(chemotaxis_strength) << std::endl
        << "RNG_SEED=" << RNG_SEED << std::endl
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