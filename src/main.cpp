#include "experiments.h"
#include "utils.h"

int main() {
  welcome();
  validate_parameters();
  write_config();
  default_testrun();
  //positional_error_experiment();
  //chemotaxis_experiment();
  //turing_patterns_experiment();
  //relax_tissue();
  Domain domain(-5, -5, 30, 20);
  Ensemble ensemble = EnsembleController::grow_tissue(42);
  ensemble.writeOFF("new.off");
}

