#include "const.h"
#include "polyhoop.h"
#include "solver.h"
#include "interpolator.h"
#include "chemistry.h"
#include "walltime.h"
#include <iostream>
#include <omp.h>

int main() {
    std::cout << "Running benchmark" << std::endl;
    rng.seed(90178009);
    Ensemble ensemble("ensemble/default.off"); // read the input file
    
    size_t L = 100;
    size_t N = L/dx; 
    size_t frames = 100;
    size_t timesteps = 1000;
    int num_threads = omp_get_num_threads();

    Grid<double> u0(L, L); // initial condition, just zeros
    Solver solver(u0, D0, dx, dt, k0); // init solver
    Interpolator interpolator(ensemble, solver);
    Chemistry chemistry(ensemble);

    std::ofstream file("benchmark.txt");
    file << "num_threads frame total_vertices time" << std::endl;

    for (std::size_t f = 1; f <= frames; ++f) {
        double start = walltime();
        for (std::size_t s = 0; s < timesteps; ++s) {
            ensemble.step();
        } 
        double end = walltime();
        size_t total_vertices = 0;
        for (auto& cell : ensemble.polygons) {
            total_vertices += cell.vertices.size();
        }
        file << num_threads << " " << f << " " << total_vertices << " " << end - start << std::endl;
    }
}