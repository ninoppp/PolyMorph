#include "ensembleController.h"

/** [NOT intended as usage example]
*
* @brief Benchmark different number of species. Not implemented yet (TODO)
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
    const double N[] = {1, 2, 4, 8};
    write_config("scaling");

    // setup output
    std::ofstream file("scaling.csv", std::ios::app);
    file << "node_id,dx,gridpoints,polygons,species,time,threads" << std::endl;

    for (int NUM_SPECIES : N) {
        std::cout << "Node ID: " << nodeID << " dx: " << dx << " threads: " << omp_get_max_threads() << std::endl;
        Domain domain(-L/2, -W/2, L/2, W/2);
        Ensemble ensemble("ensemble/tissues_varwidth/40_0.off", domain); 
        Solver solver(domain, dx, Reactions::linearDegradation); // init solver
        Interpolator interpolator(ensemble, solver);
        ensemble.is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p == Nr}; };

        double start = walltime();
        for (int i = 0; i < 1000; i++) {
            ensemble.step(); 
            interpolator.scatter();
            solver.step(dt);
            interpolator.gather();
        } 
        double time  = walltime() - start;
        std::cout << "Node ID: " << nodeID << " finished full in " << time << " seconds." << std::endl;
        
        file << nodeID << "," << dx << "," << solver.Nx * solver.Ny << "," << ensemble.polygons.size() << "," << NUM_SPECIES << "," << time << "," << omp_get_max_threads() << std::endl;
    }
}