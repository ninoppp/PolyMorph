#include "../include/ensembleController.h"

int main() {
    std::cout << "DEBUGGING NOW" << std::endl;
    double L = 10;
    Domain domain(-L/2, -L/2, L/2, L/2);
    Ensemble ensemble("ensemble/default.off", domain); 
    Solver solver(domain, dx, Reactions::linearDegradation); // init solver
    //solver.boundary.west = {BoundaryCondition::Type::Dirichlet, 0};

    Interpolator interpolator(ensemble, solver);
    ensemble.is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p == Nr}; };

    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    std::cout << "alpha of starting cell: " << ensemble.polygons[0].alpha << std::endl;
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
    return 0;
}