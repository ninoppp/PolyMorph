#include "ensembleController.h"

void chemotaxis_experiment() {
    assert(NUM_SPECIES == 1 && "Chemotaxis assumes one species");
    assert(chemotaxis_strength[0] != 0 && "Chemotaxis strength must be nonzero");
    Domain domain(-30, -15, 30, 15);
    Ensemble ensemble("ensemble/tissues_varwidth/30_0.off", domain);
    Solver solver(domain, dx, Reactions::linearDegradation); 
    solver.boundary.north = {BoundaryCondition::Type::Dirichlet, 0};
    solver.boundary.south = {BoundaryCondition::Type::Dirichlet, 0};
    solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};
    solver.boundary.west = {BoundaryCondition::Type::Dirichlet, 1};
    Interpolator interpolator(ensemble, solver);
    //ensemble.is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p % 30 == 0}; };
    ensemble.set_flag = [](const Polygon& p) { return p.vertices[0].p % 3 == 0; };
    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    EnsembleController::stop_growth(ensemble);
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step(); 
            interpolator.scatter();
            solver.step(dt);
            interpolator.gather();
        } 
        ensemble.output(f);
        solver.output(f);
    }
}

int main() {
    welcome();
    validate_parameters();
    write_config();
    chemotaxis_experiment();
}