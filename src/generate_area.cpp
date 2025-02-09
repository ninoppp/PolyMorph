#include "ensembleController.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <omp.h>

/**
 * Generate tissues with different cell area Amax_mu
 * 1 to 40 micrometer diameter
 */

int main(int argc, char* argv[]) {
  const char* nodeID_str = getenv("SLURM_NODEID");
  if (!nodeID_str) {
      std::cerr << "SLURM_NODEID is not set." << std::endl;
      //return 1;
  }
  const int nodeID = 0;//std::atoi(nodeID_str);

  assert(kh == 0); // adhesion stiffness should zero
  // maybe add empirical stability condition for dt here

  const double a_max[] = {M_PI * 2}; // M_PI / 2, M_PI, 

  write_config();
  const double W = 10; // only for test. check size of lattice sim
  const double L = 10; 
  
  //#pragma omp parallel for collapse(2) num_threads(120)
  //for (int seed = nodeID*10; seed < (nodeID+1)*10; seed++) {
    for (int i = 0; i < 1; i++) {
      int seed = 0; // tmp
      double Amax_mu = a_max[i];

      #pragma omp critical
      std::cout << "Core " << omp_get_thread_num() << " generating OFF with Amax_mu=" << Amax_mu << " and seed=" << seed << std::endl;

      const double start = walltime();
      Domain domain(-L/2, -W/2, L/2, W/2);      
      Ensemble ensemble("ensemble/default.off", domain, seed);
      ensemble.Amax_dist = create_lognormal({Amax_mu}, {Amax_CV})[0]; // change area distribution
      EnsembleController::redraw_params_from_dists(ensemble); // update all params to account for updated distribution
      EnsembleController::grow_tissue(ensemble, true); // grow tissue until domain is filled
      const double end = walltime();

      #pragma omp critical
      std::cout << "Core " << omp_get_thread_num() << " finished generating in " << end - start << " seconds. Measuring real area now" << std::endl;

      // measure real CV
      std::vector<double> areas;
      for (int p = 0; p < ensemble.polygons.size(); p++) {
        areas.push_back(ensemble.polygons[p].A);
      }
      double real_area = mean(areas); 
      int real_area_trunc = int(real_area * 1e4); // 6 digits assuming area <= 99
      
      #pragma omp critical
      {
        std::cout << "Core " << omp_get_thread_num() << " used Amax_mu=" << Amax_mu << " and measured real area=" << real_area << std::endl;
        ensemble.output(real_area_trunc);
        ensemble.write_OFF("ensemble/tissues_cell_area/" + std::to_string(real_area_trunc) + ".off");
      }
    }
  //}
  return 0;
}