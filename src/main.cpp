#include "const.h"
#include "experiments.h"
#include "utils.h"

int main() {
  welcome();
  validate_parameters();
  rng.seed(RNG_SEED);
  write_config();
  //default_testrun();
  //positional_error_experiment();
  //chemotaxis_experiment();
  turing_patterns_experiment();
}

