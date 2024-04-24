#include <vector>
#include "const.h"

struct Reaction {
    virtual double* operator()(double u[num_species], double k[num_k]) {
        double r[num_species];
        return r;
    }
};

struct LinearDegradation : Reaction {
    double* operator()(double u[num_species], double k[num_k]) {
        double r[num_species];
        for (int i = 0; i < num_species; i++) {
            r[i] = -k[0] * u[i];
        }
        return r;
    }
};