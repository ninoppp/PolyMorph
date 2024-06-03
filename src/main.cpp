#include "const.h"
#include "experiments.h"
#include "utils.h"
#include "domain.h"

int main() {
  welcome();
  validate_parameters();
  write_config();
  //default_testrun();
  //positional_error_experiment();
  //chemotaxis_experiment();
  //turing_patterns_experiment();
  //relax_tissue();
  Domain domain(-5, -5, 30, 20);
  Ensemble ensemble = grow_tissue(domain);
  ensemble.writeOFF("new.off");
}

