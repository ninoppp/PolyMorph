#include <vector>
#include <functional>
#include <fstream>
#include <cmath>

constexpr int N = 10; // size of the FD grid. Question: does N stay constant or does dx stay constant?

struct Grid {
    double G[N][N]; // ToDo: change datastructure to std::vector
    double& operator()(int i, int j) { return G[i][j]; }
};

struct Reaction {
    double operator()(double u, int i, int j) {
        return 0.0; // ToDo: implement reaction function
    }
};

struct LinearDegradation : Reaction {
    double k;
    LinearDegradation(double k): k(k) {}
    double operator()(double u, int i, int j) {
        return -k * u;
    }
};

struct Solver {    
    Grid u; // concentration of the morphogen
    double D; // diffusion coefficient. Later also a grid datastructure
    double dx; // grid spacing
    // dt comes from polyhoop main program

    // initialize the grid with a given initial condition
    Solver(const Grid u0, const double D = 1.0, const double dx = 0.1) {
        this->u = u0;
        this->D = D;
        this->dx = dx;
    }

    void step(Grid &u, const double dt, double dx, double D, Reaction R) { 
        Grid unew;
        // Forward Euler with central differences ToDo: adapt for variable diffusion coefficient
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                unew(i, j) = u(i, j) + dt * (
                    D / (dx*dx) * (u(i + 1, j) + u(i - 1, j) + u(i, j + 1) + u(i, j - 1) - 4 * u(i, j))
                    + R(u(i, j), i, j)
                );
            }
        }
        // Zero-flux boundary conditions
        for (int i = 0; i < N; i++) {
            unew(i, 0) = unew(i, 1);
            unew(i, N - 1) = unew(i, N - 2);
            unew(0, i) = unew(1, i);
            unew(N - 1, i) = unew(N - 2, i);
        }
        u = unew; 
    }

    void output(const size_t frame) {    // f: frame number
        char filename [16]; 
        snprintf(filename, 16, "rd_frame%06zu.vts", frame);        
        std::ofstream file(filename);
        // TODO write VTK file
    }
};

// helper function for IC
Grid create_gaussian() {
    Grid u;
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            u(i, j) = std::exp(-((i - N/2)*(i - N/2) + (j - N/2)*(j - N/2)) / 10.0);
    u(N/2, N/2) = 1.0;
    return u;
}