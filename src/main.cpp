#include "solver.h"
#include <iostream>
#include "const.h"
#include "polyhoop.h"
#include "interpolator.h"
#include "chemistry.h"
#include "experiments.h"
#include "utils.h"
#include <cassert>

int main() {
  welcome();
  validate_parameters();
  rng.seed(RNG_SEED);
  write_config();
  //default_testrun();
  //positional_error_experiment();
  chemotaxis_experiment();
}

