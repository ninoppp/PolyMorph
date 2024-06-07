#include "experiments.h"
#include "ensembleController.h"
#include "utils.h"

int main() {
  welcome();
  validate_parameters();
  write_config();
  //default_testrun();

  #pragma omp parallel for
  for (int seed = 0; seed < 10; seed++) {
    std::mt19937 rng; // random number generator
    rng.seed(seed);
    std::lognormal_distribution<> dist = create_lognormal({Amax_mu}, {Amax_CV})[0];
    #pragma omp critical
    std::cout << "seed=" << seed << " val=" << dist(rng) << std::endl;
  }
  return 0;
}

