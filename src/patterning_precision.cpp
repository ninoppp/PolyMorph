#include "ensembleController.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <omp.h>

void generate_OFF(int nodeID) { // 100 60x30 tissues per node

  Domain domain(-30, -15, 30, 15);

  #pragma omp parallel for num_threads(100)
  for (int seed = nodeID*100; seed < (nodeID+1)*100; seed++) {
    #pragma omp critical
    std::cout << "Core " << omp_get_thread_num() << " generating OFF with seed=" << seed << std::endl;
    Ensemble ensemble("ensemble/default.off", domain, seed);
    EnsembleController::grow_tissue(ensemble);
    std::cout << "Core " << omp_get_thread_num() << " finished generating" << std::endl;

    ensemble.output(seed);
    ensemble.writeOFF("ensemble/tissues_60x30/" + std::to_string(seed) + ".off");
  }  
}

void generate_width_OFF(int nodeID) { // 120 tissues per node (12 widths with 10 different seeds each)
  int length = 60;
  const int width[] = {2, 4, 6, 8, 10, 13, 17, 20, 25, 30, 35, 40}; // 12 widths
  #pragma omp parallel for collapse(2) num_threads(120)
  for (int seed = nodeID*10; seed < (nodeID+1)*10; seed++) {
    for (int i = 0; i < 12; i++) {
      const int w = width[i];
      Domain domain(-length/2, -w/2, length/2, w/2);
      Ensemble ensemble("ensemble/default.off", domain, seed);
      EnsembleController::grow_tissue(ensemble);
      ensemble.output(seed);
      ensemble.writeOFF("ensemble/tissues_varwidth/" + std::to_string(w) + "_" + std::to_string(seed));
    }
  }
}

void positional_error_experiment() {
  Domain domain(-30, -15, 30, 15);

  std::ofstream file("positional_error.csv");
  file << "thresh_cv,grad_cv,seed,readout_pos,prec_zone_width,time,num_threads" << std::endl;

  double cv[] = {0.01, 0.03, 0.07,
                 0.1, 0.3, 0.7, 
                 1.0, 3.0, 7.0, 
                 10};

  //omp_set_nested(1);
  omp_set_dynamic(0);


  for (double grad_cv : cv) {
    
    double lnCV = std::log(1 + grad_cv*grad_cv);
    #pragma omp parallel for num_threads(32)
    for (int seed = 0; seed < 32; seed++) {
      #pragma omp critical
      std::cout << "Calculating with tcv=" << threshold_CV[0] << " gcv=" << grad_cv << " seed=" << seed << std::endl;

      double start = walltime();
      Ensemble ensemble("ensemble/rect_60x30_nobox.off", domain, seed);
      ensemble.D_dist = create_lognormal(D_mu, {lnCV});
      ensemble.k_dist = create_lognormal(k_mu, {lnCV});
      ensemble.p_dist = create_lognormal(p_mu, {lnCV});
      ensemble.is_producing = [domain](const Polygon& p) { 
        return std::vector<bool> {p.midpoint().x < domain.x0 + 10}; 
      };

      Solver solver(domain, dx, Reactions::linearDegradation);
      Interpolator interpolator(ensemble, solver);
      solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};

      ensemble.step(); // update boxes
      interpolator.scatter();
      int num_steps = 100000; 
      for (int step = 0; step < num_steps; step++) {
        solver.step(dt);
      } 
      interpolator.gather();
      double end = walltime();

      EnsembleController::apply_flag(ensemble);
      double readout_pos = EnsembleController::mean_readout_position(ensemble, solver);
      double prec_zone_width = EnsembleController::get_precision_zone_width(ensemble);
      #pragma omp critical
      file << threshold_CV[0] << "," << grad_cv << "," << seed << "," << readout_pos << "," << prec_zone_width << "," << end - start << "," << omp_get_num_threads() << std::endl;
    }
  }

  file.close();
}



int main() {
  const char* nodeID_str = getenv("SLURM_NODEID");
  if (!nodeID_str) {
      std::cerr << "SLURM_NODEID is not set." << std::endl;
      return 1;
  }
  int nodeID = std::atoi(nodeID_str);
  write_config();
  generate_width_OFF(nodeID);
  return 0;
}