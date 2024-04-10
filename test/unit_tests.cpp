#include <iostream>
#include <iomanip> // Add this line
#include <fstream>
#include <cmath>

//#include "../solver.h"
#include "1d_solver.h"
// #include <gtest/gtest.h>

constexpr double c0 = 1;
constexpr double D = 0.03;
constexpr double j = 0.01; //influx, negative gradient at left boundary
constexpr double dx = 0.01;
constexpr double k = 0.01;
constexpr double dt = 1e-3;
constexpr double L = 100.0;
constexpr int N = L/dx;

constexpr unsigned Ns = 10000;
constexpr unsigned Nf = 100;

const double lambda = std::sqrt(D/k);
const double tau = (1+L/10/lambda ) / (2*k); // time to stead state from Roman
//constexpr double tau = L*L / (D * N*N * M_PI*M_PI); // from script for double dirichlet

void print(std::string s) {
    std::cout << s << std::endl;
}
/*/
void print_grid(Grid<double> u, int N) {
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
    Grid<double> u;
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
}*/

void run_1d_solver() {
    std::cout << "running 1d solver" << std::endl;
    std::vector<double> u0(N, 0.0);
    u0[0] = 1.0; // dirichlet boundary condition
    Solver1D solver(u0, D, dx, dt, k);
    for (int f = 0; f < Nf; f++) {
        for (int s = 0; s < Ns; s++) {
            solver.step();
        }
        //solver.output(f);
        double time = f*Ns*dt;
        std::cout << 
            "frame " << f << 
            " norm: " << solver.norm << 
            " time: " << time << 
            " tau: " << tau << std::endl;
        /*if (time > tau) {
            break;
        }*/
    }
    // final (hopefully steady) state
    std::ofstream file("test/final_frame.txt");
    for (int i = 0; i < N; i++) {
        file << solver.u[i] << std::endl;
    }
    file.close();
}

int main() {    
    run_1d_solver();
    return 0;
}
