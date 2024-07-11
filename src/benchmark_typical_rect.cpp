#include "ensembleController.h"

/*!
    * @brief Benchmark a typical simulation of rectangular tissue. Run with 8 threads
*/

int main (int argc, char *argv[]) {
    const char* nodeID_str = getenv("SLURM_NODEID");
    if (!nodeID_str) {
        std::cerr << "SLURM_NODEID is not set." << std::endl;
        return 1;
    }
    assert(alpha_mu == 1);
    const int nodeID = std::atoi(nodeID_str);
    const double L = 60;
    const double W = 30;
    const double DX[] = {0.2}; //{1.47, 1.04, 0.738, 0.522, 0.369, 0.301, 0.261}; // specifc number due to nodes per polygon estimation
    const int repetitions = 30;
    write_config();

    // setup output
    std::ofstream file("benchmark_typical_rect_02_" + std::to_string(nodeID) + ".csv", std::ios::app);
    file << "node_id,dx,gridpoints,time_ensemble,time_all,threads" << std::endl;
    for (int r = 0; r < repetitions; ++r) {
        for (double dx : DX) {
            std::cout << "Node ID: " << nodeID << " dx: " << dx << " threads: " << omp_get_max_threads() << std::endl;
            double time_ensemble, time_all;
            double num_polygons, num_gridpoints;

            { // ensemble only
                Domain domain(-L/2, -W/2, L/2, W/2);
                Ensemble ensemble("ensemble/tissues_varwidth/30_0.off", domain); 
                Solver solver(domain, dx, Reactions::linearDegradation); // init solver
                Interpolator interpolator(ensemble, solver);
                ensemble.is_producing = [domain](const Polygon& p) { return std::vector<bool> {p.midpoint().x < domain.x0 + 10}; };

                double start = walltime();
                for (std::size_t f = 1; f <= Nf; ++f) {
                    for (std::size_t s = 0; s < Ns; ++s) {
                        ensemble.step(); 
                    } 
                }
                double end  = walltime();
                time_ensemble = end - start;
                std::cout << "Node ID: " << nodeID << " finished ens-only in " << time_ensemble << " seconds. Total polygons: " << ensemble.polygons.size() << std::endl;
            }

            { // ensemble + solver + interpolation
                Domain domain(-L/2, -W/2, L/2, W/2);
                Ensemble ensemble("ensemble/tissues_varwidth/30_0.off", domain); 
                Solver solver(domain, dx, Reactions::linearDegradation); // init solver
                Interpolator interpolator(ensemble, solver);
                ensemble.is_producing = [domain](const Polygon& p) { return std::vector<bool> {p.midpoint().x < domain.x0 + 10}; };

                double start = walltime();
                for (std::size_t f = 1; f <= Nf; ++f) {
                    for (std::size_t s = 0; s < Ns; ++s) {
                        ensemble.step(); 
                        interpolator.scatter();
                        solver.step(dt);
                        interpolator.gather();
                    } 
                }
                double end  = walltime();
                time_all = end - start;
                num_gridpoints = solver.Nx * solver.Ny;
                num_polygons = ensemble.polygons.size();
                std::cout << "Node ID: " << nodeID << " finished full in " << time_all << " seconds. Total polygons: " << num_polygons << std::endl;
            }
            file << nodeID << "," << dx << "," << num_gridpoints << "," << time_ensemble << "," << time_all << "," << omp_get_max_threads() << std::endl;
        }
    }
}