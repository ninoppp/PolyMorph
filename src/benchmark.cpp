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
    
    size_t frames = 10;
    size_t timesteps = 1000;

    std::ofstream file("benchmark_solver.csv", std::ios::app);
    file << "num_threads frame num_gridpoints num_vertices time" << std::endl;

    int num_threads = omp_get_max_threads();
    
    for (int N = 2048; N <= 8192; N *= 2) {

        Grid<std::vector<double>> u0(N, N); // initial condition, just zeros
        Reaction R = linearDegradation();
        Solver solver(u0, dx, R); // init solver
        Interpolator interpolator(ensemble, solver);
        Chemistry chemistry(ensemble);
        //chemistry.is_producing = [solver](const Polygon& p) { return p.midpoint().x < solver.x0 + 10; }; // heavyside

        for (std::size_t f = 1; f <= frames; ++f) {
            double start = walltime();
            for (std::size_t s = 0; s < timesteps; ++s) {
                solver.step(dt);
            } 
            double end = walltime();
            size_t num_vertices = 0;
            for (const Polygon& p : ensemble.polygons) {
                num_vertices += p.vertices.size();
            }
            file << num_threads << " " << f << " " << N*N << " " << num_vertices << " " << end - start << std::endl;
        }
    }
    file.close();
    return 0;
}