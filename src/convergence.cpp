#include "ensembleController.h"

constexpr double C0 = 1.0;
const double lambda = std::sqrt(D0[0] / k0[0]); // using background values for this because no cells anyways

// standard model
double analytic_solution(double x) {
    return C0 * std::exp(-x/lambda); 
}

double RMSE(std::vector<double> &solution, double dx) {
    double error = 0;
    for (int i = 0; i < solution.size(); i++) {
        double x = i * dx;
        error += std::pow(analytic_solution(x) - solution[i], 2);
    }
    return std::sqrt(error / solution.size());
}

double infinity_norm(std::vector<double> &solution, double dx) {
    double error = 0;
    for (int i = 0; i < solution.size(); i++) {
        double x = i * dx;
        error = std::max(error, std::abs(analytic_solution(x) - solution[i]));
    }
    return error;
}

std::vector<double> grid_to_vec(Grid<std::vector<double>>& grid) {
    std::vector<double> vec;
    for (int i = 0; i < grid.sizeX(); i++) {
        vec.push_back(grid(i, 1)[0]);
    }
    return vec;
}

int main(int argc, char* argv[]) {
    validate_parameters();
    write_config("convergence");

    std::ofstream file("convergence.csv");
    file << "dx,abs_err" << std::endl;

    std::ofstream file_solution("solution.csv");
    file_solution << "t,x,analytic,simulated" << std::endl;
    
    double L = 100;
    double W = 3;
    Domain domain(-L/2, -W/2, L/2, W/2);
    int num_steps = 1e5;
    double dt = 1e-4; // maybe make smaller

    double DX[] = {0.3};// {4.0, 3.0, 2.0, 1.0, 0.75, 0.5, 0.3, 0.25, 0.125, 0.0625};
    
    for (double dx : DX) {
        
        Solver solver(domain, dx, Reactions::linearDegradation); // init solver
        solver.boundary.west = {BoundaryCondition::Type::Dirichlet, C0};
        solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};

        for (int s = 0; s < num_steps; ++s) {
            solver.step(dt);
            if (s % 1000 == 0) {
                auto sol = grid_to_vec(solver.u);
                for (int i = 0; i < solver.Nx; i++) {
                    file_solution << s * dt << "," << i * dx << "," << analytic_solution(i * dx) << "," << sol[i] << std::endl;
                }
                double error = RMSE(sol, dx);
                std::cout << "dx=" << dx << " step=" << s << " RMSE=" << error << std::endl;
                file << dx << "," << error << std::endl;
            }
        }
    }

    file.close();
    return 0;
}

