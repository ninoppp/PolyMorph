#include "ensembleController.h"

constexpr double C0 = 1.0;
const double lambda = std::sqrt(D0[0] / k0[0]); // using background values for this because no cells anyways

// standard model
double analytic_solution(double x) {
    return C0 * std::exp(-x/lambda); 
}

int main(int argc, char* argv[]) {
    validate_parameters();
    write_config("convergence");

    std::ofstream file("convergence.csv");
    file << "dx,abs_err" << std::endl;
    
    double L = 100;
    double W = 3;
    Domain domain(-L/2, -W/2, L/2, W/2);
    int num_steps = 1e5;
    double dt = 1e-4; // maybe make smaller

    double DX[] = {4.0, 3.0, 2.0, 1.0, 0.75, 0.5, 0.25, 0.125, 0.0625};
    
    for (double dx : DX) {
        
        Solver solver(domain, dx, Reactions::linearDegradation); // init solver
        solver.boundary.west = {BoundaryCondition::Type::Dirichlet, C0};
        solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};

        solver.output(0); // print the initial state
        for (int s = 0; s < num_steps; ++s) {
            solver.step(dt);
        }
    }

    file.close();
    return 0;
}

