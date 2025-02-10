#include "ensembleController.h"

/** [Usage example]
*
* @brief Not working yet. ToDo: Have to fine-tune parameters and initial condition
*/

int main() {
  assert(NUM_SPECIES == 2 && "Turing requires two species");
  assert(NUM_KIN == 3 && "Turing requires 3 kinetic coefficients (a, b, gamma)");
  assert(p_mu[0] == 0 && p_mu[1] == 0 && "Production rates must be zero");
  double L = 60;
  Domain domain(-L/2, -L/2, L/2, L/2);
  Ensemble ensemble("ensemble/tissue_127.off", domain);
  Solver solver(domain, dx, Reactions::turing);
  //solver.c = solver.noisy_ic(1, 0.1);
  Interpolator interpolator(ensemble, solver);
  ensemble.output(0); // print the initial state
  solver.output(0); // print the initial state
  for (std::size_t f = 1; f <= Nf; ++f) {
      for (std::size_t s = 0; s < Ns; ++s) {
          ensemble.step(); 
          interpolator.scatter();
          solver.step(dt);
          interpolator.gather();
      } 
      ensemble.output(f);
      solver.output(f);
  }
}