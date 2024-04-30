#include "solver.h"
#include <iostream>
#include "const.h"
#include "polyhoop.h"
#include "interpolator.h"
#include "chemistry.h"
#include "experiments.h"
#include "utils.h"
#include <cassert>

void welcome() {
  std::cout << "--------------------------" << std::endl
            << "|  Welcome to PolyMorph  |" << std::endl
            << "--------------------------" << std::endl;
}


int main() {
  welcome();
  validate_parameters();
  rng.seed(90178009);
  write_config();
  default_testrun();
}

