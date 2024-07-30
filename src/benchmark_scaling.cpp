#include "ensembleController.h"

/** [NOT intended as usage example]
*
* @brief Parallel scaling analysis
*
* Benchmarking done on a packed tissue without growth, similar to the ones used for positional error analysis.
* Using 3 different grid spacings.
* Meant to be run on multiple nodes simultaneously.
*/

int main (int argc, char *argv[]) {
    const char* nodeID_str = getenv("SLURM_NODEID");
    if (!nodeID_str) {
        std::cerr << "SLURM_NODEID is not set." << std::endl;
        return 1;
    }
    assert(alpha_mu == 0);
    const int nodeID = std::atoi(nodeID_str);
    const double L = 60;
    const double W = 40;
    const double DX[] = {0.5, 0.3, 0.2};
    const int repetitions = 48;
    write_config("scaling");

    // setup output
    std::ofstream file("scaling.csv", std::ios::app);
    file << "node_id,dx,gridpoints,polygons,time_ensemble,time_all,threads" << std::endl;

    for (int i = 0; i < repetitions; i++) {
        for (double dx : DX) {
            std::cout << "Node ID: " << nodeID << " dx: " << dx << " threads: " << omp_get_max_threads() << std::endl;
            
            Domain domain(-L/2, -W/2, L/2, W/2);
            Ensemble ensemble("ensemble/tissues_varwidth/40_0.off", domain); 
            Solver solver(domain, dx, Reactions::linearDegradation); // init solver
            Interpolator interpolator(ensemble, solver);
            ensemble.is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p == Nr}; };

            double start = walltime();
            for (int i = 0; i < 1000; i++) {
                ensemble.step(); 
            }
            double time_ensemble  = walltime() - start;
            std::cout << "Node ID: " << nodeID << " finished ens-only in " << time_ensemble << " seconds." << std::endl;

            start = walltime();
            for (int i = 0; i < 1000; i++) {
                ensemble.step(); 
                interpolator.scatter();
                solver.step(dt);
                interpolator.gather();
            } 
            double time_all  = walltime() - start;
            std::cout << "Node ID: " << nodeID << " finished full in " << time_all << " seconds." << std::endl;
            
            file << nodeID << "," << dx << "," << solver.Nx * solver.Ny << "," << ensemble.polygons.size() << "," << time_ensemble << "," << time_all << "," << omp_get_max_threads() << std::endl;
        }
    }
}