#include "experiments.h"
#include "ensembleController.h"
#include "utils.h"

int main() {
  welcome();
  validate_parameters();
  write_config();
  //default_testrun();
  //positional_error_experiment();
  //chemotaxis_experiment();
  //turing_patterns_experiment();
  //relax_tissue();
  Domain domain(-30, -15, 30, 15);
  #pragma omp parallel for
  for (int seed = 0; seed < 10; seed++) {
    #pragma omp critical
    std::cout << "Running with seed=" << seed << std::endl;
    Ensemble ensemble("ensemble/default.off", domain, seed);
    EnsembleController::grow_tissue(ensemble);
    ensemble.output(seed);
    char filename [16];
    snprintf(filename, 16, "poserr-%05zu.off", seed);
    ensemble.writeOFF(filename);
  }
}

