#include "ensembleController.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <omp.h>

/**
 * Generate tissues with different cell area variability Amax_CV
 */

int main(int argc, char* argv[]) {
  const char* nodeID_str = getenv("SLURM_NODEID");
  if (!nodeID_str) {
      std::cerr << "SLURM_NODEID is not set." << std::endl;
      return 1;
  }
  const int nodeID = 0;//std::atoi(nodeID_str);

  assert(kh == 0); // adhesion stiffness should zero
  assert(dt == 0.5 * 1e-4); // improve stability
 
  const double cv[] = {0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0};
  const int num_cv = 11;

  write_config();
  const double W = 20; // only for test. check size of lattice sim
  const double L = 20; 

  //#pragma omp parallel for collapse(2) num_threads(120)
  //for (int seed = nodeID*10; seed < (nodeID+1)*10; seed++) {
    for (int i = 0; i < 11; i++) {
      int seed = 0; // tmp
      double Amax_CV = cv[i];

      #pragma omp critical
      std::cout << "Core " << omp_get_thread_num() << " generating OFF with Amax_CV=" << Amax_CV << " and seed=" << seed << std::endl;

      const double start = walltime();
      Domain domain(-L/2, -W/2, L/2, W/2);      
      Ensemble ensemble("ensemble/default.off", domain, seed);
      ensemble.Amax_dist = create_lognormal({Amax_mu}, {Amax_CV})[0]; // change area distribution
      EnsembleController::redraw_params_from_dists(ensemble); // update all params to account for updated distribution
      EnsembleController::grow_tissue(ensemble); // grow tissue until domain is filled
      const double end = walltime();

      #pragma omp critical
      std::cout << "Core " << omp_get_thread_num() << " finished generating in " << end - start << " seconds. Measuring real CV now" << std::endl;

      // measure real CV
      std::vector<double> areas;
      for (int p = 0; p < ensemble.polygons.size(); p++) {
        areas.push_back(ensemble.polygons[p].A);
      }
      double real_CV = stddev(areas) / mean(areas); 
      int real_CV_trunc = int(real_CV * 1e4); // 6 digits assuming cv <= 99
      
      #pragma omp critical
      {
        std::cout << "Core " << omp_get_thread_num() << " used CV=" << Amax_CV << " and measured real CV=" << real_CV << std::endl;
        ensemble.output(real_CV_trunc);
        ensemble.write_OFF("ensemble/tissues_cell_area_cv/" + std::to_string(real_CV_trunc) + ".off");
      }
    }
  //}
  return 0;
}