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

        for (int i = 0; i < 200; i++) {
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
