#include "ensembleController.h"
#include "utils.h"
#include <fstream>
#include <omp.h>

int main(int argc, char* argv[]) {
    const char* nodeID_str = getenv("SLURM_NODEID");
    if (!nodeID_str) {
        std::cerr << "SLURM_NODEID is not set." << std::endl;
        return 1;
    }
    int nodeID = std::atoi(nodeID_str);
    std::cout << "Node ID: " << nodeID << std::endl;
    write_config();

    std::string csv_filename = "positional_error_cv_" + std::to_string(nodeID) + ".csv";
    std::ofstream file(csv_filename, std::ios::app);
    file << "threshold,grad_cv,seed,readout_pos,prec_zone_width,time,num_threads" << std::endl;
    
    Domain domain(-30, -15, 30, 15);
    double cv[] = {1.0, 3.0, 5.0, 10}; /*{0.01, 0.03, 0.07,
                    0.1, 0.3, 0.5, 0.7, 
                    1.0, 3.0, 5.0, 7.0, 
                    10};*/

    //omp_set_nested(1);
    omp_set_dynamic(0);

    #pragma omp parallel for collapse(2) num_threads(120) // 120 calculations per node. 1200 total
    for (int i = 0; i < 4; i++) { // 12 cv values TODO revert back
        for (int seed = nodeID*30; seed < (nodeID+1)*30; seed++) { // 10 seeds
            
            double grad_cv = cv[i];
            
            std::string off_file = "ensemble/tissues_varwidth/30_" + std::to_string(seed) + ".off";
            Ensemble ensemble(off_file.c_str(), domain, seed);

            #pragma omp critical
            std::cout << "Core " << omp_get_thread_num() << " calculating with gcv=" << grad_cv << " seed=" << seed << " and off file " << off_file << std::endl;
            
            // create distributions
            ensemble.D_dist = create_lognormal(D_mu, {grad_cv}); // FIX! has to be done before constructing
            ensemble.k_dist = create_lognormal(k_mu, {grad_cv});
            ensemble.p_dist = create_lognormal(p_mu, {grad_cv});
            ensemble.is_producing = [domain](const Polygon& p) { // left side is producing
                return std::vector<bool> {p.midpoint().x < domain.x0 + 10}; 
            };
            EnsembleController::redraw_params_from_dists(ensemble); // apply the new distributions

            Solver solver(domain, dx, Reactions::linearDegradation);
            Interpolator interpolator(ensemble, solver);
            solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};

            double start = walltime();
            ensemble.step(); // update boxes
            interpolator.scatter();
            int num_steps = 200000;
            for (int step = 0; step < num_steps; step++) {
                solver.step(dt);
            } 
            double end = walltime();

            interpolator.gather();
            EnsembleController::apply_flag(ensemble);
            double readout_pos = EnsembleController::mean_readout_position(ensemble, solver);
            double prec_zone_width = EnsembleController::get_precision_zone_width(ensemble);

            #pragma omp critical
            {
                std::cout << "Core " << omp_get_thread_num() << " finished calculating with gcv=" << grad_cv << " seed=" << seed << " in " << end - start << " seconds" << std::endl;
                file << threshold_mu[0] << "," << grad_cv << "," << seed << "," << readout_pos << "," << prec_zone_width << "," << end - start << "," << omp_get_num_threads() << std::endl;
                ensemble.output(i * 1000 + seed); // to inspect the final state (1.2k frames)
            }
        }
    }
    file.close();
    return 0;
}