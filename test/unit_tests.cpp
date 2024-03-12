#include "../solver.h"
#include <iostream>
#include <iomanip> // Add this line
#include <fstream>

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

// only run solver same as in polyhoop
void run_solver() {
    Grid u;
    Solver solver(u, 0.03, 0.01, 1e-4, LinearDegradation(0.01));
    for (int f = 0; f < 100; f++) {
        for (int s = 0; s < 3000; s++) {
            solver.step();
        }
        //solver.output(f);
        std::cout << "frame " << f << std::endl;
    }
    // final (hopefully steady) state
    std::ofstream file("test/final_frame.txt");
    for (int j = 0; j < Ny; j++) {
        file << solver.u(1,j) << std::endl;
    }
    file.close();
}
    
int main() {    
    run_solver();
    return 0;
}

