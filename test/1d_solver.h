#include <vector>
#include "../solver.h"

struct Solver {  
    std::vector<double> u; // concentration of the morphogen
    double D; // diffusion coefficient. Later also a grid datastructure
    double dx; // grid spacing
    double dt; // time step
    Reaction R; // reaction term

    // initialize the grid with a given initial condition
    Solver(std::vector<double> u0, const double D = 1.0, const double dx = 0.1, 
            double dt = 1e-4, Reaction R = LinearDegradation(0.1)) {
        this->u = u0;
        this->D = D;
        this->dx = dx;
        this->dt = dt;
        this->R = R;
    }

    void step() { 
        std::vector<double> unew;
        // Forward Euler with central differences ToDo: adapt for variable diffusion coefficient
        // maybe separate inner nodes and boundary to efficiently parallelize and vectorize inner nodes
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {   
                // mirror past-boundary nodes
                double n = (j == Ny-1) ? u(i, j-1) : u(i, j+1);
                double s = (j == 0)    ? 1         : u(i, j-1); // dirichlet boundary condition
                double e = (i == Nx-1) ? u(i-1, j) : u(i+1, j);
                double w = (i == 0)    ? u(i+1, j) : u(i-1, j);
                unew(i, j) = u(i, j) + dt * (
                    D / (dx*dx) * (n + s + e + w - 4 * u(i, j))
                    + R(u(i, j), i, j)
                ); 
            }
        }      
        u = unew; 
    }
};