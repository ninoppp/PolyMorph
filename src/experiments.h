#include "const.h"
#include "polyhoop.h"
#include "solver.h"
#include "interpolator.h"
#include "chemistry.h"
#include "utils.h"
#include "domain.h"
#include <iostream>
#include <cmath>

/*! \brief This file contains multiple "main" functions for different experiments
 *
 *  The main routines are called from the main() function in main.cpp
 *  and contain the logic for different simulation runs.
 * 
 * ToDo: split into separate cpp files in an experiments folder
 */

void default_testrun() {
    double L = 40;
    Domain domain(-L/2, -L/2, L/2, L/2);
    Ensemble ensemble("ensemble/default.off", domain); 
    Solver solver(domain, dx, Reactions::linearDegradation); // init solver
    //solver.boundary.west = {BoundaryCondition::Type::Dirichlet, 0};

    Interpolator interpolator(ensemble, solver);
    Chemistry chemistry(ensemble, solver);
    chemistry.is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p == Nr}; };

    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step(); 
            chemistry.update();
            interpolator.scatter();
            solver.step(dt);
            interpolator.gather();
            domain.step(dt); // only needed if domain is growing/shrinking
        } 
        ensemble.output(f);
        solver.output(f);
    }
}

void chemotaxis_experiment() {
    assert(NUM_SPECIES == 1 && "Chemotaxis assumes one species");
    double L = 40;
    Domain domain(-L/2, -L/2, L/2, L/2);
    Ensemble ensemble("ensemble/tissue_127.off", domain);
    unsigned N = L/dx; 
    Solver solver(domain, dx, Reactions::linearDegradation); 
    Interpolator interpolator(ensemble, solver);
    Chemistry chemistry(ensemble, solver);
    chemistry.is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p == Nr}; };
    chemistry.flag = [](const Polygon& p) { return p.vertices[0].p % 2 == 1; }; // flag every 2nd cell
    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step(); 
            chemistry.update();
            interpolator.scatter();
            solver.step(dt);
            interpolator.gather();
        } 
        ensemble.output(f);
        solver.output(f);
    }
}

void turing_patterns_experiment() {
  assert(NUM_SPECIES == 2 && "Turing requires two species");
  assert(NUM_KIN == 3 && "Turing requires 3 kinetic coefficients (a, b, gamma)");
  assert(p_mu[0] == 0 && p_mu[1] == 0 && "Production rates must be zero");
  double L = 60;
  Domain domain(-L/2, -L/2, L/2, L/2);
  Ensemble ensemble("ensemble/tissue_127.off", domain);
  Solver solver(domain, dx, Reactions::turing);
  solver.u = solver.noisy_ic(1, 0.1);
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

void two_opposing() {
  Domain domain(-50, -25, 50, 25);
  Ensemble ensemble("ensemble/rect_100x50.off", domain); // read the input file
  assert(Nr == 1 && "Nr must be 1 for running in box");
  Reaction R = Reactions::inhibition;
  Solver solver(domain, dx, R); // init solver
  Interpolator interpolator(ensemble, solver);
  Chemistry chemistry(ensemble, solver);
  chemistry.is_producing = [domain](const Polygon& p) { 
      return std::vector<bool> {p.midpoint().x < domain.x0 + 10, 
                                p.midpoint().x > domain.x0 + 90}; 
    };
  
  ensemble.output(0); // print the initial state
  solver.output(0); // print the initial state
  ensemble.step(); 
  chemistry.update();
  interpolator.scatter();
  for (std::size_t f = 1; f <= Nf; ++f) {
      for (std::size_t s = 0; s < Ns; ++s) {
          solver.step(dt);
      } 
      chemistry.update();
      interpolator.gather();
      ensemble.output(f);
      solver.output(f);
  }
}


void grow_tissue() {
  Domain domain(-30, -15, 30, 15);
  Ensemble ensemble("ensemble/default.off", domain); // read the input file
  assert(Nr == 0 && "Nr must be 0");
  double eps = 0.1;
  
  for (int f = 0; f < Nf; f++) {
    for (int s = 0; s < Ns; s++) {
      ensemble.step();
    }
    std::cout << "Frame " << f << " num_polygons " << ensemble.polygons.size() << std::endl;
  }
  ensemble.writeOFF("rect_60x30_nobox.off");
}


void positional_error_experiment() { // ToDo: fix
    std::ofstream file("sharpness.csv");
    file << "frame border_x pos_err" << std::endl;
    Domain domain(-30, -15, 30, 15);
    Ensemble ensemble("ensemble/rect_60x30_nobox.off", domain); // read the input file
    assert(Nr == 0 && "Nr must be 0");
    Reaction reaction = Reactions::linearDegradation;
    Solver solver(domain, dx, reaction); // init solver
    solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};
    Interpolator interpolator(ensemble, solver);
    Chemistry chemistry(ensemble, solver);
    chemistry.is_producing = [domain](const Polygon& p) { 
      return std::vector<bool> {p.midpoint().x < domain.x0 + 10}; 
    };
    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    chemistry.update();
    interpolator.scatter();
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            solver.step(dt);
        } 
        interpolator.gather();
        chemistry.update();
        ensemble.output(f);
        solver.output(f);
        double readout_pos = chemistry.mean_readout_position();
        //file << f << "," << mean << std::endl;
    }
    file.close();
}

/*void sharpness_experiment() {
  std::ofstream file("sharpness.csv");
  file << "ThreshCV GradCV sharpness" << std::endl;
  double CV[] = {0.1, 0.2, 0.3, 0.5, 0.7, 1.0};

  for (double thresh_cv : CV) {
    double tlnCV = std::log(1 + thresh_cv*thresh_cv);
    std::lognormal_distribution threshold_dist = std::lognormal_distribution<double>(std::log(threshold_mu[0]) - tlnCV/2, std::sqrt(tlnCV));
    
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
        chemistry.is_producing = [solver](const Polygon& p) { return p.midpoint().x < solver.x0 + 10; }; // heavyside
        
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
        file << thresh_cv << " " << grad_cv << " " << chemistry.get_precision_zone_width() << std::endl;
      }
    } 
  }
  file.close();
}*/

/*
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
}*/