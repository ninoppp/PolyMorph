#include "ensembleController.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <omp.h>

int main(int argc, char* argv[]) {
  const char* nodeID_str = getenv("SLURM_NODEID");
  if (!nodeID_str) {
      std::cerr << "SLURM_NODEID is not set." << std::endl;
      //return 1;
  }
  int nodeID = 0;//std::atoi(nodeID_str);

  //assert(kh == 0); // adhesion stiffness should zero

  double cv[] = {0.01, 0.03, 0.07,
                    0.1, 0.3, 0.5, 0.7, 
                    1.0, 3.0, 5.0, 7.0, 10};

  write_config();
  Domain domain(-10, -10, 10, 10);

  //#pragma omp parallel for collapse(2) num_threads(100)
  //for (int seed = nodeID*10; seed < (nodeID+1)*10; seed++) {
    //for (int i = 0; i < 10; i++) {
      int seed = 0; // tmp
      double Amax_CV = 10;//cv[i];

      #pragma omp critical
      std::cout << "Core " << omp_get_thread_num() << " generating OFF with alpha_cv=" << Amax_CV << " and seed=" << seed << std::endl;

      double start = walltime();
      Domain domain(-30, -15, 30, 15);
      Ensemble ensemble("ensemble/default.off", domain, seed);
      ensemble.Amax_dist = create_lognormal({Amax_mu}, {Amax_CV})[0];
      EnsembleController::redraw_params_from_dists(ensemble);
      EnsembleController::grow_tissue(ensemble);
      double end = walltime();

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
        ensemble.writeOFF("ensemble/tissues_vararea/" + std::to_string(real_CV_trunc) + ".off");
      }
    //}
  //}
  return 0;
}