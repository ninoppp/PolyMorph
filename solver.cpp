#include <vector>
#include <functional>

constexpr int N = 100; // size of the FD grid

struct Solver {
    struct Grid {
        double G[N][N]; // ToDo: change datastructure to make N more flexible. 
        double& operator()(int i, int j) { return G[i][j]; } // access element
    };

    struct Reaction {
        double operator()(double u, int i, int j) {
            return 0.0; // ToDo: implement reaction function
        }
    };

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
        u = unew; 
    }
};

