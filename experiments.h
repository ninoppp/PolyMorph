#include "const.h"
#include "polyhoop.h"
#include "solver.h"
#include "interpolator.h"
#include "chemistry.h"
#include <iostream>

void write_config() {
    std::ofstream config("simulation.cfg");
    config
        << "Date=" << __DATE__ << std::endl
        << "Time=" << __TIME__ << std::endl
        << "D0=" << D0 << std::endl
        << "k0=" << k0 << std::endl
        << "p0=" << p0 << std::endl
        << "D_mu=" << D_mu << std::endl
        << "k_mu=" << k_mu << std::endl
        << "p_mu=" << p_mu << std::endl
        << "threshold_mu=" << threshold_mu << std::endl 
        << "D_CV=" << D_CV << std::endl
        << "k_CV=" << k_CV << std::endl
        << "p_CV=" << p_CV << std::endl
        << "threshold_CV=" << threshold_CV << std::endl
        << "dx=" << dx << std::endl
        << "dt=" << dt << std::endl
        << "Nf=" << Nf << std::endl
        << "Ns=" << Ns << std::endl
        << "Nr=" << Nr << std::endl;
    config.close();
}

void default_testrun() {
    Ensemble ensemble("ensemble/default.off"); // read the input file
    unsigned L = 100;
    unsigned N = L/dx; 
    Grid<double> u0(N, N); // initial condition, just zeros
    Solver solver(u0, D0, dx, dt, k0); // init solver
    Interpolator interpolator(ensemble, solver);
    Chemistry chemistry(ensemble);
    chemistry.is_producing = [](const Polygon& p) { return p.vertices[0].p == Nr; }; // mother cell
    
    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step(); 
            chemistry.update();
            interpolator.scatter();
            solver.step();
            interpolator.gather();
        } 
        ensemble.output(f);
        solver.output(f);
    }
}

void sharpness_experiment() {
  std::ofstream file("sharpness.csv");
  file << "ThreshCV GradCV sharpness" << std::endl;
  double CV[] = {0.1, 0.2, 0.3, 0.5, 0.7, 1.0};

  for (double thresh_cv : CV) {
    double tlnCV = std::log(1 + thresh_cv*thresh_cv);
    threshold_dist = std::lognormal_distribution<double>(std::log(threshold_mu) - tlnCV/2, std::sqrt(tlnCV));
    
    for (double grad_cv : CV) {
      for (int rep = 0; rep < 5; rep++) { // repeat each experiment 5 times

        double lnCV = std::log(1 + grad_cv*grad_cv);
        D_dist = std::lognormal_distribution<double>(std::log(D_mu) - lnCV/2, std::sqrt(lnCV));
        k_dist = std::lognormal_distribution<double>(std::log(k_mu) - lnCV/2, std::sqrt(lnCV));
        p_dist = std::lognormal_distribution<double>(std::log(p_mu) - lnCV/2, std::sqrt(lnCV));

        Ensemble ensemble("ensemble/rect_100x50.off");
        Chemistry chemistry(ensemble);
        Grid<double> u0(int(102/dx), int(52/dx)); // initial condition, just zeros
        Solver solver(u0, D0, dx, dt, k0); // init solver
        Interpolator interpolator(ensemble, solver);
        chemistry.is_producing = [solver](const Polygon& p) { return p.midpoint().x < solver.box_position_x + 10; }; // heavyside
        
        //ensemble.output(cv*10000); // print the initial state
        //solver.output(cv*10000); // print the initial state
        ensemble.step(); // update boxes and everything
        chemistry.update();
        interpolator.scatter(); 
        for (std::size_t f = 1; f <= Nf; ++f)
        {
          for (std::size_t s = 0; s < Ns; ++s) 
          {
            solver.step();
          } 
          //ensemble.output(int(cv*10000 + f));
          //solver.output(int(cv*10000 + f));
        }
        interpolator.gather(); // get u
        chemistry.update(); // flag
        file << thresh_cv << " " << grad_cv << " " << chemistry.get_border_sharpness() << std::endl;
      }
    } 
  }
  file.close();
}

void differentiation_experiment() {
    Ensemble ensemble("ensemble/singlecell.off"); // read the input file
    unsigned L = 50;
    unsigned N = L/dx; 
    Grid<double> u0(N, N); // initial condition, just zeros
    Solver solver(u0, D0, dx, dt, k0); // init solver
    Interpolator interpolator(ensemble, solver);
    Chemistry chemistry(ensemble);
    chemistry.growth_control = true; // stop growth if flagged
    chemistry.is_producing = [](const Polygon& p) { return p.vertices[0].p == Nr; }; // mother cell
    
    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step(); 
            chemistry.update();
            interpolator.scatter();
            solver.step();
            interpolator.gather();
        } 
        ensemble.output(f);
        solver.output(f);
        if (f == Nf/2) { // stop production after half-time
            chemistry.is_producing = [](const Polygon& p) { return false; }; 
        }
    }
}