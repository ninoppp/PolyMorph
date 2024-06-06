#include "ensembleController.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <omp.h>

int main(int argc, char* argv[]) {
  if (argc < 2) {
      std::cerr << "Node ID not provided." << std::endl;
      return 1;
  }
  int nodeID = std::stoi(argv[1]);

  write_config();
  Domain domain(-30, -15, 30, 15);

  #pragma omp parallel for num_threads(128)
  for (int seed = nodeID*128; seed < (nodeID+1)*128; seed++) {
    
    #pragma omp critical
    std::cout << "Core " << omp_get_thread_num() << " generating OFF with seed=" << seed << std::endl;

    double start = walltime();
    Ensemble ensemble("ensemble/default.off", domain, seed);
    EnsembleController::grow_tissue(ensemble);
    double end = walltime();

    #pragma omp critical
    std::cout << "Core " << omp_get_thread_num() << " finished generating in " << end - start << " seconds" << std::endl;

    ensemble.output(seed);
    ensemble.writeOFF("ensemble/tissues_60x30/" + std::to_string(seed) + ".off");
  }  
  return 0;
}
