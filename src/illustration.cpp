#include "experiments.h"
#include "ensembleController.h"
#include "utils.h"

void illustration() {
    double L = 50;
    Domain domain(-L/2, -L/2, L/2, L/2);
    Ensemble ensemble("ensemble/tissue_127.off", domain); 
    Solver solver(domain, 0.26, Reactions::linearDegradation); // init solver
    //solver.boundary.west = {BoundaryCondition::Type::Dirichlet, 0};
    Interpolator interpolator(ensemble, solver);
    ensemble.is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p == Nr}; };
    double grad_cv = 10;
    ensemble.D_dist = create_lognormal(D_mu, {grad_cv}); // FIX! has to be done before constructin  
    ensemble.k_dist = create_lognormal(k_mu, {grad_cv});
    ensemble.p_dist = create_lognormal(p_mu, {grad_cv});
    ensemble.is_producing = [domain](const Polygon& p) { // left side is producing
        return std::vector<bool> {p.midpoint().x < domain.x0 + 10}; 
    };
    EnsembleController::redraw_params_from_dists(ensemble); // apply the new distributions
    ensemble.step();
    interpolator.scatter();
    solver.step(dt);
    interpolator.gather();
    //ensemble.output(0); // print the initial state
    //solver.output(0); // print the initial state
    // Measure actual mu of D
    double Dmu = 0;
    for (auto& p : ensemble.polygons) {
        Dmu += p.D[0];
    }
    Dmu /= ensemble.polygons.size();
    std::cout << "Dmu=" << Dmu << std::endl;
}

int main(int argc, char* argv[]) {
  welcome();
  validate_parameters();
  //write_config();
  illustration();
}