#include "experiments.h"
#include "ensembleController.h"
#include "utils.h"

void default_testrun() {
    double L = 50;
    Domain domain(-L/2, -L/2, L/2, L/2);
    Ensemble ensemble("ensemble/default.off", domain); 
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
            //domain.step(dt); // only needed if domain is growing/shrinking
        } 
        ensemble.output(f);
        solver.output(f);
    }
    ensemble.writeOFF("ensemble/typical.off");
    // relax
    EnsembleController::stop_growth(ensemble);
    for (std::size_t f = Nf + 1; f <= Nf + 20; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step();
            interpolator.scatter();
            solver.step(dt);
            interpolator.gather();
        }
        ensemble.output(f);
        solver.output(f);
    }
    ensemble.writeOFF("ensemble/typical_relaxed.off");
}

int main(int argc, char* argv[]) {
  welcome();
  validate_parameters();
  write_config();
  //default_testrun();
  Domain domain(-50, -50, 50, 50);
  Ensemble tmp ("ensemble/typical.off", domain);
  std::cout << "Number of polygons: " << tmp.polygons.size() << std::endl;
}