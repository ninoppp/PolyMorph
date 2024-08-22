#include "ensembleController.h"

/** [Usage example]
*
* @brief Standard simulation example. Exponential growth, starting from a single cell
*/

int main(int argc, char* argv[]) {
    welcome();
    validate_parameters(); // checks that correct number of kinetic parameters are set
    write_config();
    double L = 60;
    // set up core components
    Domain domain(-L/2, -L/2, L/2, L/2); // initialize square domain
    Ensemble ensemble("ensemble/default.off", domain); // init ensemble
    Reaction reaction = Reactions::linearDegradation; // define reaction model
    Solver solver(domain, dx, reaction); // init solver
    Interpolator interpolator(ensemble, solver); // init interpolator
    // optional: set boundary conditions (default zero-flux)
    solver.boundary.west = {BoundaryCondition::Type::Dirichlet, 0}; 
    // optional: set production lambda (default no cell-based production)
    ensemble.is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p == 0}; }; // starting cell (index 0) produces. Vector for multiple species
    // optional: when are cells "flagged"
    ensemble.set_flag = [](const Polygon& p) { return p.u[0] < p.threshold[0]; }; // cells are flagged when local concentration drops below threshold

    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step(); // advance mechanical ensemble
            interpolator.scatter(); // interpolate data to grid
            solver.step(dt); // advance chemical solver
            interpolator.gather(); // interpolate data back to ensemble
            EnsembleController::stop_growth_if_flagged(ensemble); // flagged cells stop growing
            //domain.step(dt); // only needed if domain is growing/shrinking
        } 
        ensemble.output(f); // print frame
        solver.output(f); // print frame
    }
    ensemble.writeOFF("final_state.off"); // save final state, can reuse as input later
}