#ifndef DOMAIN_H
#define DOMAIN_H

#include "geometry.h"
#include "const.h"

struct Domain {
    double x0, y0, x1, y1;
    // Growth normal to the domain.
    // E, N, W, S. Positive values mean growth, negative mean shrinkage. 
    double growth_rate[4]; // [L/T]
    double stiffness = kr; // [1/T^2] repulsion stiffness per vertex mass

    Domain() : x0(0), y0(0), x1(1), y1(1) {}
    Domain(double x0, double y0, double x1, double y1) : 
            x0(x0), y0(y0), x1(x1), y1(y1) {}

    double width() const { return x1 - x0; }
    double height() const { return y1 - y0; }

    bool contains(const Point& r) const {
        return r.x >= x0 && r.x <= x1 && r.y >= y0 && r.y <= y1;
    }

    void step(double dt) {
        x1 += growth_rate[0] * dt;
        y1 += growth_rate[1] * dt;
        x0 -= growth_rate[2] * dt;
        y0 -= growth_rate[3] * dt;
    }

    // uniform growth in all directions
    void set_growth_rate(double rate) {
        growth_rate[0] = growth_rate[1] = growth_rate[2] = growth_rate[3] = rate;
    }

    // different growth rates in x and y directions
    void set_growth_rate(double rate_x, double rate_y) {
        growth_rate[0] = rate_x;
        growth_rate[1] = rate_y;
        growth_rate[2] = rate_x;
        growth_rate[3] = rate_y;
    }

    // different growth rates in all directions
    void set_growth_rate(double rate_E, double rate_N, double rate_W, double rate_S) {
        growth_rate[0] = rate_E;
        growth_rate[1] = rate_N;
        growth_rate[2] = rate_W;
        growth_rate[3] = rate_S;
    }
};

#endif // DOMAIN_H