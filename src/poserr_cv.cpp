#include "ensembleController.h"
#include "utils.h"
#include <fstream>
#include <omp.h>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Node ID not provided." << std::endl;
        return 1;
    }
    int nodeID = std::stoi(argv[1]);

    std::string csv_filename = "positional_error_cv_" + std::to_string(nodeID) + ".csv";
    std::ofstream file(csv_filename, std::ios::app);
    file << "threshold,grad_cv,seed,readout_pos,prec_zone_width,time,num_threads" << std::endl;
    
    Domain domain(-30, -15, 30, 15);
    double cv[] = {0.01, 0.03, 0.07,
                    0.1, 0.3, 0.7, 
                    1.0, 3.0, 7.0, 
                    10};

    //omp_set_nested(1);
    omp_set_dynamic(0);

    for (double grad_cv : cv) { // 12 iterations
        double lnCV = std::log(1 + grad_cv*grad_cv);

        #pragma omp parallel for num_threads(128) // 128 seeds parallel for this specific cv. Same seed was also used for tissue generation
        for (int seed = nodeID*128; seed < (nodeID+1)*128; seed++) { // ToDo: maybe add more nodes to seeds

            #pragma omp critical
            std::cout << "Core " << omp_get_thread_num() << " calculating with gcv=" << grad_cv << " seed=" << seed << std::endl;

            std::string off_file = "ensemble/tissues_60x30/" + std::to_string(seed) + ".off";
            Ensemble ensemble(off_file.c_str(), domain, seed);
            
            // create distributions
            ensemble.D_dist = create_lognormal(D_mu, {lnCV});
            ensemble.k_dist = create_lognormal(k_mu, {lnCV});
            ensemble.p_dist = create_lognormal(p_mu, {lnCV});
            ensemble.is_producing = [domain](const Polygon& p) { // left side is producing
                return std::vector<bool> {p.midpoint().x < domain.x0 + 10}; 
            };

            Solver solver(domain, dx, Reactions::linearDegradation);
            Interpolator interpolator(ensemble, solver);
            solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};

            double start = walltime();
            ensemble.step(); // update boxes
            interpolator.scatter();
            int num_steps = 100000; // TODO verify
            for (int step = 0; step < num_steps; step++) {
                solver.step(dt);
            } 
            interpolator.gather();
            double end = walltime();

            EnsembleController::apply_flag(ensemble);
            double readout_pos = EnsembleController::mean_readout_position(ensemble, solver);
            double prec_zone_width = EnsembleController::get_precision_zone_width(ensemble);

            #pragma omp critical
            file << threshold_mu[0] << "," << grad_cv << "," << seed << "," << readout_pos << "," << prec_zone_width << "," << end - start << "," << omp_get_num_threads() << std::endl;
        }
    }
    file.close();
}