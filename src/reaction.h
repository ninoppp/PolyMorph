#ifndef REACTION_H
#define REACTION_H

#include <vector>
#include "const.h"

using Reaction = std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&)>;

namespace Reactions {

    Reaction linearDegradation = [](const std::vector<double>& u, const std::vector<double>& k) {
        std::vector<double> r(NUM_SPECIES);
        for (int i = 0; i < NUM_SPECIES; i++) {
            r[i] = -k[i]*u[i];
        }
        return r;
    };

    Reaction inhibition = [](const std::vector<double>& u, const std::vector<double>& k) {
        std::vector<double> r(NUM_SPECIES);
        r = {
            -k[0]*u[0] - k[2]*u[0]*u[1], // species 2 inhibits species 1
            -k[1]*u[1]
        };
        return r;
    };

    Reaction turing = [](const std::vector<double>& u, const std::vector<double>& k) {
        std::vector<double> r(NUM_SPECIES);
        const double a = k[0];
        const double b = k[1];
        const double gamma = k[2];
        r = {
            gamma * (a - u[0] + u[0]*u[0]*u[1]),
            gamma * (b -        u[0]*u[0]*u[1])
        };
        return r;
    };
}

#endif