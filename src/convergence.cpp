#include "ensembleController.h"


/** [Usage example]
*
* @brief Convergence analysis of the finite difference solver
*/

constexpr double C0 = 1.0;
const double lambda = std::sqrt(D0[0] / k0[0]); // using background values for this because we don't use cells here anyways

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

// takes the middle row of the long grid and stores it in a vector for comparison
std::vector<double> grid_to_vec(Grid<std::vector<double>>& grid) {
    std::vector<double> vec;
    for (int i = 0; i < grid.sizeX(); i++) {
        vec.push_back(grid(i, 1)[0]);
    }
    return vec;
}

int main(int argc, char* argv[]) {
    welcome();
    validate_parameters();
    write_config("convergence");

    std::ofstream file("convergence.csv");
    file << "dx,rmse,inf_norm" << std::endl;

    std::ofstream file_solution("solution.csv"); // solution values to animate the convergence
    file_solution << "dx,t,x,analytic,simulated" << std::endl;
    
    double dt = 1e-4; // maybe make smaller
    double T = 10; // final time
    int num_steps = T / dt;

    double DX[] = {5.0, 4.0, 3.0, 2.0, 1.5, 1.0, 0.75, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.125, 0.1};
    
    for (double dx : DX) {
        double L = 300; // long enough to model x -> infinity
        double W = 2*dx; // 3 nodes in y direction
        Domain domain (-L/2, -W/2, L/2, W/2);
        Solver solver(domain, dx, Reactions::linearDegradation); // init solver
        solver.boundary.west = {BoundaryCondition::Type::Dirichlet, C0};
        solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};

        for (int s = 0; s < num_steps; ++s) {
            solver.step(dt);
            if (s % 1000 == 0) {
                
                auto sol = grid_to_vec(solver.u);
                double rmse = RMSE(sol, dx);
                double inf_norm = infinity_norm(sol, dx);
                std::stringstream ss;
                ss << "dx=" << dx << " t=" << s * dt << " step=" << s << " RMSE=" << rmse << " inf_norm=" << inf_norm;
                LOG(ss.str());
                for (int i = 0; i < solver.Nx; i++) {
                    file_solution << dx << "," << s * dt << "," << i * dx << "," << analytic_solution(i * dx) << "," << sol[i] << std::endl;
                }
            }   
        }
        auto sol = grid_to_vec(solver.u);
        double rmse = RMSE(sol, dx);
        double inf_norm = infinity_norm(sol, dx);
        file << dx << "," << rmse << "," << inf_norm << std::endl;
    }

    file.close();
    return 0;
}

