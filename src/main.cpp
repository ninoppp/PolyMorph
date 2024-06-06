#include "experiments.h"
#include "ensembleController.h"
#include "utils.h"

int main() {
  welcome();
  validate_parameters();
  write_config();
  std::string s = "ensemble/tissues_60x30/" + std::to_string(2) + ".off";
  std::string off_file = "ensemble/tissues_60x30/" + std::to_string(69) + ".off";
  std::cout << off_file.c_str() << std::endl;
  //default_testrun();
  //positional_error_experiment();
  //chemotaxis_experiment();
  //turing_patterns_experiment();
  //relax_tissue();
}

