#include <vector>
#include "const.h"

// ToDo: replace functors with lambdas

struct Reaction {
    std::vector<double> operator()(const std::vector<double>& u, const std::vector<double>& k) {
        std::vector<double> r(NUM_SPECIES);
        return r;
    }
};

struct LinearDegradation : Reaction {
    std::vector<double> operator()(const std::vector<double>& u, const std::vector<double>& k) {
        std::vector<double> r(NUM_SPECIES);
        for (int i = 0; i < NUM_SPECIES; i++) {
            r[i] = -k[i]*u[i];
        }
        return r;
    }
};

struct Inhibition : Reaction {
    std::vector<double> operator()(const std::vector<double>& u, const std::vector<double>& k) {
        std::vector<double> r(NUM_SPECIES);
        r = {
            -k[0]*u[0] - k[2]*u[0]*u[1], // species 2 inhibits species 1
            -k[1]*u[1]
        };
        return r;
    }
};

struct Turing : Reaction { // ToDo: get a,b,gamma as kinetic coefficients
    double a, b, gamma;
    Turing(double a, double b, double gamma) : a(a), b(b), gamma(gamma) {}
    std::vector<double> operator()(const std::vector<double>& u, const std::vector<double>& k) {
        std::vector<double> r(NUM_SPECIES);
        r = {
            gamma * (a - u[0] + u[0]*u[0]*u[1]),
            gamma * (b -        u[0]*u[0]*u[1])
        };
        return r;
    }
};  