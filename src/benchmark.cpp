#include "const.h"
#include "ensemble.h"
#include "solver.h"
#include "interpolator.h"
#include "ensembleController.h"
#include <iostream>
#include <omp.h>

void bench() {
    std::cout << "Benchmarking..." << std::endl;

    std::ofstream file("benchmark.csv", std::ios::app);
    file << "step,gridpoints,polygons,vertices,num_species,ensemble,solver,scatter,gather,all,threads" << std::endl; 
    
    int num_threads = omp_get_max_threads();

    for (int num_polygons = 1; num_polygons <= 1e5; num_polygons *= 10) {
        Domain domain(-30, -30, 30, 30);
        std::string off_filename = "ensemble/tissues_bynum/" + std::to_string(num_polygons) + ".off";
        Ensemble ensemble(off_filename.c_str(), domain);
        Solver solver(domain, dx, Reactions::linearDegradation);
        solver.boundary.west = {BoundaryCondition::Type::Dirichlet, 1};
        solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};
        Interpolator interpolator(ensemble, solver);
        EnsembleController::stop_growth(ensemble);

        for (int i = 0; i < 10000; i++) {
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
            int num_vertices = 0;
            for (const Polygon& p : ensemble.polygons) {
                num_vertices += p.vertices.size();
            }
            file << i << "," << solver.Nx * solver.Ny << "," << ensemble.polygons.size() << "," << num_vertices << "," << NUM_SPECIES << "," 
                 << ens << "," << sol << "," << scat << "," << gath << "," << all << num_threads << std::endl;
        }
    }
    file.close();
    write_config();
}


void grow() {
    std::cout << "Benchmarking..." << std::endl;

    std::ofstream file("benchmark_grow_ens.csv", std::ios::app);
    file << "gridpoints,polygons,vertices,num_species,time,threads" << std::endl; 
    
    int num_threads = omp_get_max_threads();

    Domain domain(-30, -30, 30, 30);
    Ensemble ensemble("ensemble/default.off", domain);
    Solver solver(domain, dx, Reactions::linearDegradation);
    solver.boundary.west = {BoundaryCondition::Type::Dirichlet, 1};
    solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};
    Interpolator interpolator(ensemble, solver);

    int num_polygons = 1;
    double start = walltime();
    while (true) {
        ensemble.step();
        //interpolator.scatter();
        //solver.step(dt);
        //interpolator.gather();
        // record every doubling of polygons
        if (ensemble.polygons.size() >= num_polygons * 2) {
            num_polygons *= 2;
            double timestamp = walltime();
            int num_vertices = 0;
            for (const Polygon& p : ensemble.polygons) {
                num_vertices += p.vertices.size();
            }
            file << solver.Nx * solver.Ny << "," << ensemble.polygons.size() << "," << num_vertices << "," << NUM_SPECIES << "," << timestamp - start << "," << num_threads << std::endl;
            std::cout << "Polygons: " << ensemble.polygons.size() << " Time: " << timestamp - start << std::endl;
        }
        // stop when corners of domain are filled
        if (solver.parent_idx(1, 1) >= Nr && solver.parent_idx(solver.Nx - 2, solver.Ny - 2) >= Nr 
            && solver.parent_idx(1, solver.Ny - 2) >= Nr && solver.parent_idx(solver.Nx - 2, 1) >= Nr) {
            break;
        }
    }
    
    file.close();
}

int main(int argc, char* argv[]) {
    write_config("bench");
    grow();
    return 0;
}