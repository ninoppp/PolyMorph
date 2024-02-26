#include "../solver.h"
#include <iostream>
#include <iomanip> // Add this line

// #include <gtest/gtest.h>

void print(std::string s) {
    std::cout << s << std::endl;
}

void print_grid(Grid u, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << std::setw(4) << std::fixed << std::setprecision(4) << u(i, j) << " ";            
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    print("2x2 test grid");
    Grid u;
    u(0, 0) = 1.0;
    u(0, 1) = 2.0;
    u(1, 0) = 3.0;
    u(1, 1) = 4.0;
    print_grid(u, 2);

    print("Gaussian test grid");
    Grid u_gauss = create_gaussian();
    print_grid(u_gauss, N);

    // ToDo: create plot of u_gauss

    print("Stepping previous grid once");
    Solver s(u_gauss, 1.0, 0.1, 1e-4);
    s.step(u_gauss, s.dt, s.dx, s.D, LinearDegradation(0.1));
    print_grid(u_gauss, N);

    return 0;
}

