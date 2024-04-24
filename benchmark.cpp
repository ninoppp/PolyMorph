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
    Ensemble ensemble("ensemble/box_100x50.off"); // read the input file
    
    size_t N = 256;
    size_t frames = 100;
    size_t timesteps = 2000;

    Grid<double> u0(N, N); // initial condition, just zeros
    Solver solver(u0, D0, dx, dt, k0); // init solver
    Interpolator interpolator(ensemble, solver);
    Chemistry chemistry(ensemble);
    chemistry.is_producing = [solver](const Polygon& p) { return p.midpoint().x < solver.box_position_x + 10; }; // heavyside

    std::ofstream file("benchmark_ensemble_box.txt", std::ios::app);
    file << "num_threads frame num_gridpoints num_vertices time" << std::endl;

    int num_threads = omp_get_max_threads();
    for (std::size_t f = 1; f <= frames; ++f) {
        double start = walltime();
        for (std::size_t s = 0; s < timesteps; ++s) {
            ensemble.step();
        } 
        double end = walltime();
        size_t num_vertices = 0;
        for (const Polygon& p : ensemble.polygons) {
            num_vertices += p.vertices.size();
        }
        file << num_threads << " " << f << " " << N*N << " " << num_vertices << " " << end - start << std::endl;
    }
    file.close();
    return 0;
}