#include <vector>
#include <functional>
#include <fstream>

constexpr int N = 100; // size of the FD grid

struct Grid {
    double G[N][N]; // ToDo: change datastructure to make N more flexible. 
    double& operator()(int i, int j) { return G[i][j]; } // access element
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
    double D; // diffusion coefficient

    void step(Grid &u, double dt, double dx, double D, Reaction R) { // D later changes to a grid too
        Grid unew;
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

    void output(int frame) {    // f: frame number
        std::string filename = "frame" + std::to_string(frame) + ".vts";
        std::ofstream file(filename);
        // TODO write VTK file
    }
};
