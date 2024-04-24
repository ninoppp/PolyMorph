#include <vector>
#include "const.h"

struct Reaction {
    virtual std::vector<double> operator()(std::vector<double> u, std::vector<double> k) {
        std::vector<double> r(num_species);
        return r;
    }
};

struct LinearDegradation : Reaction {
    std::vector<double> operator()(std::vector<double> u, std::vector<double> k) {
        std::vector<double> r(num_species);
        for (int i = 0; i < num_species; i++) {
            r[i] = -k[i]*u[i];
        }
        return r;
    }
};