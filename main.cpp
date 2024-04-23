#include "const.h"
#include "polyhoop.h"
#include "solver.h"
#include "interpolator.h"
#include "chemistry.h"
#include <iostream>

void welcome() {
  std::cout << "--------------------------" << std::endl
            << "|  Welcome to PolyMorph  |" << std::endl
            << "--------------------------" << std::endl;
  // ToDo: log configuration for reproducibility
}

int main()
{
  welcome();
  rng.seed(90178009);
  Ensemble ensemble("ensemble/default.off"); // read the input file
  
  unsigned L = 50;
  unsigned N = L/dx; 
  Grid<double> u0(N, N); // initial condition, just zeros
  Solver solver(u0, D0, dx, dt, k0); // init solver

  Interpolator interpolator(ensemble, solver);
  Chemistry chemistry(ensemble);
  chemistry.is_producing = [](const Polygon& p) { return p.vertices[0].p == Nr; }; // mother cell
  
  ensemble.output(0); // print the initial state
  solver.output(0); // print the initial state
  for (std::size_t f = 1; f <= Nf; ++f)
  {
    for (std::size_t s = 0; s < Ns; ++s) 
    {
      ensemble.step();
      //chemistry.update();
      interpolator.scatter(); 
      solver.step();
      interpolator.gather();
    } 
    ensemble.output(f); // print a frame
    solver.output(f); // print a frame
  }
  //ensemble.writeOFF("new.off");
}