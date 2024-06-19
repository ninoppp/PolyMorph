#include "ensembleController.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <omp.h>

int main(int argc, char* argv[]) {
  /*if (argc < 2) {
      std::cerr << "Node ID not provided." << std::endl;
      return 1;
  }
  int nodeID = std::stoi(argv[1]);*/
  int nodeID = 0;

  double cv[] = {0.01, 0.03, 0.07,
                    0.1, 0.3, 0.5, 0.7, 
                    1.0, 3.0, 5.0, 7.0, 10};

  write_config();

  #pragma omp parallel for collapse(2) num_threads(120)
  for (int seed = nodeID*10; seed < (nodeID+1)*10; seed++) {
    for (int i = 0; i < 8; i++) {
      double area_cv = cv[0]; // TODO remove hotfix

      #pragma omp critical
      std::cout << "Core " << omp_get_thread_num() << " generating OFF with alpha_cv=" << area_cv << " and seed=" << seed << std::endl;

      double start = walltime();
      Domain domain(-30, -15, 30, 15);
      Ensemble ensemble("ensemble/default.off", domain, seed);
      //ensemble.Amax_dist = create_lognormal(p_mu, {area_cv})[0];
      //EnsembleController::redraw_params_from_dists(ensemble); // not necessarily needed since starting with 1 cell only
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
        ensemble.output(real_CV_trunc);
        ensemble.writeOFF("ensemble/tissues_vararea/" + std::to_string(real_CV_trunc) + ".off");
      }
    }
  }
  return 0;
}