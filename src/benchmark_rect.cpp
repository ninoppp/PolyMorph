#include <fstream>
#include <omp.h>
#include <cassert>

#include "utils.h"
#include "ensembleController.h"

/*
This measures runtime of all components respective to the number of gridpoints and polygons.
No scaling.
*/
int main(int argc, char* argv[]) {
    const char* nodeID_str = getenv("SLURM_NODEID");
    if (!nodeID_str) {
        std::cerr << "SLURM_NODEID is not set." << std::endl;
        return 1;
    }
    int nodeID = std::atoi(nodeID_str);
    assert(Amax_mu == 0);
    write_config();

    double length = 60;
    const double widths[] = {2, 4, 6, 8, 10, 13, 17, 20, 25, 30, 35, 40}; // 12 widths
    const double DX[] = {2.0, 1.0, 0.75, 0.5, 0.25, 0.125}; // 6 dx vallues

    std::string csv_filename = "benchmark_rect_" + std::to_string(nodeID) + ".csv";
    std::ofstream file(csv_filename, std::ios::app);
    file << "width,gridpoints,polygons,vertices,num_species,ensemble,solver,scatter,gather,all,threads" << std::endl;

    #pragma omp parallel for collapse(2) 
    for (int j = 0; j < 6; j++) {
        for (int i = 0; i < 12; i++) {
            double dx = DX[j];
            int width = widths[i];
            Domain domain(-length/2, -width/2, length/2, width/2);

            std::string off_file = "ensemble/tissues_varwidth/" + std::to_string(width) + "_0" + ".off";
            Ensemble ensemble(off_file.c_str(), domain);
            ensemble.is_producing = [domain](const Polygon& p) { // left side is producing
                return std::vector<bool> {p.midpoint().x < domain.x0 + 10}; 
            };
            int num_vertices = 0;
            for (const Polygon& p : ensemble.polygons) {
                num_vertices += p.vertices.size();
            }
            Solver solver(domain, dx, Reactions::linearDegradation);
            Interpolator interpolator(ensemble, solver);
            solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};
            #pragma omp critical
            std::cout << "Benchmarking with " << ensemble.polygons.size() << " polygons and " << solver.Nx * solver.Ny << " gridpoints" << std::endl;

            double all = walltime();
            double ens = walltime();
            for (int step = 0; step < Ns; ++step) {
                ensemble.step();
            }
            ens = walltime() - ens;
            double sol = walltime();
            for (int step = 0; step < Ns; ++step) {
                solver.step(dt);
            }
            sol = walltime() - sol;
            double scat = walltime();
            for (int step = 0; step < Ns; ++step) {
                interpolator.scatter();
            }
            scat = walltime() - scat;
            double gath = walltime();
            for (int step = 0; step < Ns; ++step) {
                interpolator.gather();
            }
            gath = walltime() - gath;
            all = walltime() - all;
            
            #pragma omp critical
            file << width << "," << solver.Nx * solver.Ny << "," << ensemble.polygons.size() << "," << num_vertices << "," << NUM_SPECIES << "," 
                    << ens << "," << sol << "," << scat << "," << gath << "," << all << "," << omp_get_num_threads() << std::endl;
        }
    }
    file.close();
    return 0;
}