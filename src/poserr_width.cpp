#include <fstream>
#include <omp.h>
#include <cassert>

#include "utils.h"
#include "ensembleController.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Node ID not provided." << std::endl;
        return 1;
    }
    int nodeID = std::stoi(argv[1]);

    int length = 60;
    const int widths[] = {2, 4, 6, 8, 10, 13, 17, 20, 25, 30, 35, 40}; // 12 widths
    assert(D_CV[0] == 0.3 && k_CV[0] == 0.3 && p_CV[0] == 0.3); // use CV=0.3 for noisy gradient

    std::string csv_filename = "positional_error_width_" + std::to_string(nodeID) + ".csv";
    std::ofstream file(csv_filename, std::ios::app);
    file << "threshold,width,seed,readout_pos,prec_zone_width,time,num_threads" << std::endl;

    //omp_set_nested(1);
    omp_set_dynamic(0);

    #pragma omp parallel for collapse(2) num_threads(120)
    for (int seed = nodeID*10; seed < (nodeID+1)*10; seed++) {
        for (int i = 0; i < 12; i++) {
            int width = widths[i];
            Domain domain(-length/2, -width/2, length/2, width/2);

            #pragma omp critical
            std::cout << "Core " << omp_get_thread_num() << " calculating with width=" << width << " seed=" << seed << std::endl;

            std::string off_file = "ensemble/tissues_varwidth/" + std::to_string(width) + "_" + std::to_string(seed) + ".off";
            Ensemble ensemble(off_file.c_str(), domain, seed);
            ensemble.is_producing = [domain](const Polygon& p) { // left side is producing
                return std::vector<bool> {p.midpoint().x < domain.x0 + 10}; 
            };
            // CONTINUE HERE

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
            file << threshold_mu[0] << "," << width << "," << seed << "," << readout_pos << "," << prec_zone_width << "," << end - start << "," << omp_get_num_threads() << std::endl;
        }
    }
    file.close();
}