#include "ensembleController.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <omp.h>

int main(int argc, char* argv[]) {
  assert(Amax_CV != 0 && alpha_CV != 0);
  assert(kh == 0); // adhesion stiffness should zero

  int nodeID = 0;
  if (argc < 2) {
    std::cerr << "Node ID not provided. Assuming 0" << std::endl;
  } else {
    nodeID = std::stoi(argv[1]);
    std::cout << "Node ID: " << nodeID << std::endl;
  }

  write_config();
  int length = 60;
  const int width[] = {2, 4, 6, 8, 10, 13, 17, 20, 25, 30, 35, 40}; // 12 widths

  #pragma omp parallel for collapse(2) num_threads(120)
  for (int seed = nodeID*10; seed < (nodeID+1)*10; seed++) {
    for (int i = 0; i < 12; i++) {

      const int w = width[i];
      Domain domain(-length/2, -w/2, length/2, w/2);
      #pragma omp critical
      std::cout << "Core " << omp_get_thread_num() << " generating OFF with width = " << w << " and seed=" << seed << std::endl;
      
      double start = walltime();
      Ensemble ensemble("ensemble/default.off", domain, seed);
      EnsembleController::grow_tissue(ensemble);
      double end = walltime();

      #pragma omp critical
      std::cout << "Core " << omp_get_thread_num() << " finished generating in " << end - start << " seconds" << std::endl;

      ensemble.output(w * 1000 + seed);
      ensemble.writeOFF("ensemble/tissues_varwidth/" + std::to_string(w) + "_" + std::to_string(seed) + ".off");
    }
  }
  return 0;
}