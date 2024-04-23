#include "const.h"
#include "polyhoop.h"
#include "solver.h"
#include "interpolator.h"
#include "chemistry.h"
#include <iostream>

void sharpness_experiment() { // TODO: modify distributions
  double CV[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0};
  std::ofstream file("sharpness.txt");
  for (double cv : CV) {
    Ensemble ensemble("ensemble/rect_100x50.off");
    Chemistry chemistry(ensemble);
    Grid<double> u0(int(100/dx), int(50/dx)); // initial condition, just zeros
    Solver solver(u0, D0, dx, dt, k0); // init solver
    Interpolator interpolator(ensemble, solver);
    chemistry.is_producing = [solver](const Polygon& p) { return p.midpoint().x < solver.box_position_x + 10; }; // heavyside
    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    ensemble.step(); // update boxes and everything
    chemistry.update();
    interpolator.scatter(); 
    for (std::size_t f = 1; f <= Nf; ++f)
    {
      for (std::size_t s = 0; s < Ns; ++s) 
      {
        solver.step();
      } 
    }
    interpolator.gather(); // get u
    chemistry.update(); // flag
    ensemble.output(int(cv*100)); // print a frame
    solver.output(int(cv*100)); // print a frame
    file << "CV=" << cv << " sharpness=" << chemistry.get_border_sharpness() << std::endl;
  }
  file.close();
}

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
  
  /*Ensemble ensemble("ensemble/rect_100x50.off"); // read the input file
  
  unsigned L = 100;
  unsigned N = L/dx; 
  Grid<double> u0(int(100/dx), int(50/dx)); // initial condition, just zeros
  Solver solver(u0, D0, dx, dt, k0); // init solver

  Interpolator interpolator(ensemble, solver);
  Chemistry chemistry(ensemble);
  //chemistry.is_producing = [](const Polygon& p) { return p.vertices[0].p == Nr; }; // mother cell
  chemistry.is_producing = [](const Polygon& p) { return p.midpoint().x < -40; }; // heavyside

  ensemble.output(0); // print the initial state
  solver.output(0); // print the initial state
  ensemble.step(); // update boxes and everything
  chemistry.update();
  interpolator.scatter(); // scatter initial conditions
  for (std::size_t f = 1; f <= Nf; ++f)
  {
    for (std::size_t s = 0; s < Ns; ++s) 
    {
      solver.step();
    } 
    interpolator.gather();
    chemistry.update();
    ensemble.output(f); // print a frame
    solver.output(f); // print a frame
  }*/
  //ensemble.writeOFF("new.off");
}

