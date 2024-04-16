#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "../solver.h"
#include "1d_solver.h"

constexpr double c0 = 1;
constexpr double D = 0.03;
constexpr double j = 0.03; //influx, negative gradient at left boundary
constexpr double dx = 0.01;
constexpr double k = 0.05;
constexpr double dt = 1e-3;
constexpr double L = 50.0;
constexpr int N = L/dx;

constexpr unsigned Ns = 1000;
constexpr unsigned Nf = 100;

const double lambda = std::sqrt(D/k);
const double tau_standard = (1+L/10/lambda ) / (2*k); // time to stead state from Roman
//constexpr double tau = L*L / (D * N*N * M_PI*M_PI); // from script for double dirichlet
const double tau_neumann = L*L / (D*M_PI*M_PI + k*L*L); // not working somehow
const double tau = tau_standard;

void print(std::string s) {
    std::cout << s << std::endl;
}

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
    std::cout << "running 2d solver" << std::endl;
    Grid<double> u(N, 3);
    Solver solver(u, D, dx, dt, LinearDegradation(k));
    for (int f = 0; f < Nf; f++) {
        for (int s = 0; s < Ns; s++) {
            solver.step();
        }
        //solver.output(f);
        double time = f*Ns*dt;
        std::cout << 
            "frame " << f <<
            " time: " << time << 
            " tau: " << tau << std::endl;
    }
    // final (hopefully steady) state
    std::ofstream file("test/final_frame.txt");
    for (int i = 0; i < solver.Nx; i++) {
        file << solver.u(i,1) << std::endl;
    }
    file.close();
}

void run_1d_solver() {
    std::cout << "running 1d solver" << std::endl;
    std::vector<double> u0(N, 0.0);
    u0[0] = 1.0; // dirichlet boundary condition
    Solver1D solver(u0, D, dx, dt, k, "2xNeumann", j);
    for (int f = 0; f < Nf; f++) {
        for (int s = 0; s < Ns; s++) {
            solver.step();
        }
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
    //run_1d_solver();
    run_solver();
    return 0;
}
