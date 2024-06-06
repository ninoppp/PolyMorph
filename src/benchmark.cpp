#include "const.h"
#include "ensemble.h"
#include "solver.h"
#include "interpolator.h"
#include "ensembleController.h"
#include <iostream>
#include <omp.h>

/*
void solver_scaling() {
    std::cout << "Running benchmark" << std::endl;
    Domain domain(-20, -20, 20, 20); // domain
    Ensemble ensemble("ensemble/default.off"); // read the input file
    ensemble.rng.seed(90178009);

    size_t frames = 10;
    size_t timesteps = 1000;

    std::ofstream file("benchmark_solver.csv", std::ios::app);
    file << "num_threads frame num_gridpoints num_vertices time" << std::endl;

    int num_threads = omp_get_max_threads();
    
    for (int N = 2048; N <= 8192; N *= 2) {

        Grid<std::vector<double>> u0(N, N); // initial condition, just zeros
        Reaction R = Reactions::linearDegradation;
        Solver solver(, dx, R); // init solver
        Interpolator interpolator(ensemble, solver);
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
}*/

void bench() {

    std::ofstream file("benchmark.csv", std::ios::app);
    file << "step,gridpoints,polygons,ensemble,solver,scatter,gather,all" << std::endl; 
    
    // TODO: measure with different gridpoint and polygon counts
    
    Domain domain(-30, -30, 30, 30); // ToDo: adjust to tissue size
    Ensemble ensemble("ensemble/default.off", domain);
    ensemble.rng.seed(90178009);
    Solver solver(domain, dx, Reactions::linearDegradation);
    solver.boundary.west = {BoundaryCondition::Type::Dirichlet, 1};
    solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};
    Interpolator interpolator(ensemble, solver);

    for (int i = 0; i < 1000; i++) {
        double all = walltime();
        double ens = walltime();
        ensemble.step();
        ens = walltime() - ens;
        double sol = walltime();
        solver.step(dt);
        sol = walltime() - sol;
        double scat = walltime();
        interpolator.scatter();
        scat = walltime() - scat;
        double gath = walltime();
        interpolator.gather();
        gath = walltime() - gath;
        all = walltime() - all;
        file << i << "," << solver.Nx * solver.Ny << "," << ensemble.polygons.size() << "," << ens << "," << sol << "," << scat << "," << gath << "," << all << std::endl;
    }
}

int main() {
    std::cout << "Running benchmark. Nothing here yet" << std::endl;
}
