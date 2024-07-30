#include "ensembleController.h"
#include "utils.h"

int main() {
    welcome();
    validate_parameters();
    write_config("velocity");
    assert(ADVECTION_DILUTION_EN && "Advection-dilution must be enabled");
    double L = 25;
    Domain domain(-L/2, -L/2, L/2, L/2);
    domain.set_growth_rate(3.0);
    Ensemble ensemble("ensemble/tissue_127.off", domain); 
    Solver solver(domain, dx, Reactions::linearDegradation); // init solver
    //solver.boundary.west = {BoundaryCondition::Type::Dirichlet, 0};
    Interpolator interpolator(ensemble, solver);
    ensemble.is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p == Nr}; };
    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step(); 
            interpolator.scatter();
            solver.step(dt);
            interpolator.gather();
            domain.step(dt); // only needed if domain is growing/shrinking
        } 
        ensemble.output(f);
        solver.output(f);
    }
}
