#include <fstream>
#include <omp.h>
#include <cassert>

#include "utils.h"
#include "ensembleController.h"

/** [NOT intended as usage example]
*
* @brief This program measures runtime of all components (ensemble, solver, scatter, gather) respective to the number of gridpoints and polygons.
* Designed to be run on a cluster node with 128 threads.
*/

int main(int argc, char* argv[]) {
    const char* nodeID_str = getenv("SLURM_NODEID");
    if (!nodeID_str) {
        std::cerr << "SLURM_NODEID is not set." << std::endl;
        return 1;
    }
    int nodeID = std::atoi(nodeID_str);
    write_config();

    double length = 60;
    const double widths[] = {2, 4, 6, 8, 10, 13, 17, 20, 25, 30, 35, 40}; // 12 widths
    const double DX[] = {2.0, 1.0, 0.75, 0.5, 0.25, 0.125}; // 6 dx vallues

    std::string csv_filename = "benchmark_rect.csv";
    std::ofstream file(csv_filename, std::ios::app);
    file << "width,dx,thread_num,gridpoints,polygons,vertices,num_species,ensemble,solver,scatter,gather,all" << std::endl;

    omp_set_dynamic(0);
    omp_set_nested(0);

    #pragma omp parallel for num_threads(128)
    for (int t = 0; t < 128; t++) {
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
                file << width << "," << dx << "," << omp_get_thread_num() << "," << solver.Nx * solver.Ny << "," << ensemble.polygons.size() << "," << num_vertices << "," << NUM_SPECIES << "," 
                        << ens << "," << sol << "," << scat << "," << gath << "," << all << std::endl;
            }
        }
    }
    file.close();
    return 0;
}