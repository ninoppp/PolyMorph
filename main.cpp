#include "const.h"
#include "polyhoop.h"
#include "solver.h"
#include "interpolator.h"
#include "chemistry.h"
#include <iostream>

void sharpness_experiment() {
  double CV[] = {0.5, 0.7, 1.0};
  std::ofstream file("sharpness.txt");
  for (double cv : CV) {
    // adjust distributions
    double lnCV = std::log(1 + cv*cv);
    D_dist = std::lognormal_distribution<double>(std::log(D_mu) - lnCV/2, std::sqrt(lnCV));
    k_dist = std::lognormal_distribution<double>(std::log(k_mu) - lnCV/2, std::sqrt(lnCV));
    p_dist = std::lognormal_distribution<double>(std::log(p_mu) - lnCV/2, std::sqrt(lnCV));
    // ToDo: independent
    threshold_dist = std::lognormal_distribution<double>(std::log(threshold_mu) - lnCV/2, std::sqrt(lnCV));

    Ensemble ensemble("ensemble/rect_100x50.off");
    Chemistry chemistry(ensemble);
    Grid<double> u0(int(102/dx), int(52/dx)); // initial condition, just zeros
    Solver solver(u0, D0, dx, dt, k0); // init solver
    Interpolator interpolator(ensemble, solver);
    chemistry.is_producing = [solver](const Polygon& p) { return p.midpoint().x < solver.box_position_x + 10; }; // heavyside
    ensemble.output(cv*10000); // print the initial state
    solver.output(cv*10000); // print the initial state
    ensemble.step(); // update boxes and everything
    chemistry.update();
    interpolator.scatter(); 
    for (std::size_t f = 1; f <= Nf; ++f)
    {
      for (std::size_t s = 0; s < Ns; ++s) 
      {
        solver.step();
      } 
      interpolator.gather(); // get u
      chemistry.update(); // flag
      ensemble.output(int(cv*10000 + f)); // print a frame
      solver.output(int(cv*10000 + f)); // print a frame
    }
    file << "CV=" << cv << " sharpness=" << chemistry.get_border_sharpness() << std::endl;
  }
  file << "config:" << std::endl
      << "D0=" << D0 << std::endl
      << "D_mu=" << D_mu << std::endl
      << "k_mu=" << k_mu << std::endl
      << "p_mu=" << p_mu << std::endl
      << "threshold_mu=" << threshold_mu << std::endl
      << "dx=" << dx << std::endl
      << "dt=" << dt << std::endl
      << "Nf=" << Nf << std::endl
      << "Ns=" << Ns << std::endl;
      
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

