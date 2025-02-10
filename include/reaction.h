#ifndef REACTION_H
#define REACTION_H

#include <vector>
#include "const.h"

/**
 * @brief Reaction models for the reaction-diffusion equations
 * 
 * Each reaction takes two vectors: 
 * - c: the local concentrations of each species
 * - k: all local kinetic parameters for the reaction
 * 
 * And returns a vector
 * - r: the local reaction contributions for each species
 * 
 * Local meaning at a specific grid point. 
 */

using Reaction = std::function<std::vector<double>(const std::vector<double>& c, const std::vector<double>& k, double t)>;

namespace Reactions {

Reaction linearDegradation = [](const std::vector<double>& c, const std::vector<double>& k, double t) {
    std::vector<double> r(NUM_SPECIES);
    for (int i = 0; i < NUM_SPECIES; i++) {
        r[i] = -k[i]*c[i];
    }
    return r;
};

Reaction inhibition = [](const std::vector<double>& c, const std::vector<double>& k, double t) {
    std::vector<double> r(NUM_SPECIES);
    r = {
        -k[0]*c[0] - k[2]*c[0]*c[1], // species 2 inhibits species 1
        -k[1]*c[1]
    };
    return r;
};

Reaction turing = [](const std::vector<double>& c, const std::vector<double>& k, double t) {
    std::vector<double> r(NUM_SPECIES);
    const double a = k[0];
    const double b = k[1];
    const double gamma = k[2];
    r = {
        gamma * (a - c[0] + c[0]*c[0]*c[1]),
        gamma * (b -        c[0]*c[0]*c[1])
    };
    return r;
};

}

#endif