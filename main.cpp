#include "const.h"
#include "polyhoop.h"
#include "solver.h"
#include "interpolator.h"
#include "chemistry.h"
#include "experiments.h"
#include <iostream>

void welcome() {
  std::cout << "--------------------------" << std::endl
            << "|  Welcome to PolyMorph  |" << std::endl
            << "--------------------------" << std::endl;
}

int main()
{
  welcome();
  rng.seed(90178009);

  sharpness_experiment();
  write_config();
}

