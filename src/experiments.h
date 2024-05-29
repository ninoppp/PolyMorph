#include "const.h"
#include "polyhoop.h"
#include "solver.h"
#include "interpolator.h"
#include "measurements.h"
#include "utils.h"
#include "domain.h"
#include <omp.h>
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
    is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p == Nr}; };

    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step(); 
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
    is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p == Nr}; };
    set_flag = [](const Polygon& p) { return p.vertices[0].p % 2 == 1; }; // flag every 2nd cell
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

// generate tissue filling out the given domain
void grow_tissue(Domain& domain, const char* filename) {
  Ensemble ensemble("ensemble/default.off", domain); // read the input file
  assert(Nr == 0 && "Nr must be 0");
  double eps = 0.1;
  for (int f = 0; f < Nf; f++) {
    for (int s = 0; s < Ns; s++) {
      ensemble.step();
    }
    std::cout << "Frame " << f << " num_polygons " << ensemble.polygons.size() << std::endl;
  }
  ensemble.writeOFF(filename);
}


void positional_error_experiment() {
  Domain domain(-30, -15, 30, 15);
  is_producing = [domain](const Polygon& p) { 
    return std::vector<bool> {p.midpoint().x < domain.x0 + 10}; 
  };

  std::ofstream file("positional_error.csv");
  file << "thresh_cv,grad_cv,seed,readout_pos,time" << std::endl;

  double cv[] = {0.01, 0.05, 0.1, 0.3, 0.7, 1.0};
  double cv_small[] = {0.3};

  omp_set_nested(1);
  omp_set_dynamic(0);

  for (double thresh_cv : cv_small) {
    for (double grad_cv : cv_small) {
      
      double lnCV = std::log(1 + grad_cv*grad_cv);
      double tlnCV = std::log(1 + thresh_cv*thresh_cv);
      D_dist = create_lognormal(D_mu, {lnCV});
      k_dist = create_lognormal(k_mu, {lnCV});
      p_dist = create_lognormal(p_mu, {lnCV});
      threshold_dist = create_lognormal(threshold_mu, {tlnCV});

      #pragma omp parallel for num_threads(16)
      for (int seed = 0; seed < 16; seed++) {

        double start = walltime();
        Ensemble ensemble("ensemble/rect_60x30_nobox.off", domain);
        ensemble.rng.seed(seed);
        Solver solver(domain, dx, Reactions::linearDegradation);
        Interpolator interpolator(ensemble, solver);
        solver.boundary.east = {BoundaryCondition::Type::Dirichlet, 0};

        ensemble.step(); // update boxes
        interpolator.scatter();
        int num_steps = 100000; 
        for (int step = 0; step < num_steps; step++) {
          solver.step(dt);
        } 
        interpolator.gather();
        double end = walltime();

        Measurements::apply_flag(ensemble);
        double readout_pos = Measurements::mean_readout_position(ensemble, solver);
        #pragma omp critical
        file << thresh_cv << "," << grad_cv << "," << seed << "," << readout_pos << "," << end - start << std::endl;
      }
    }
  }
  file.close();
}

void relax_tissue() {
  Domain domain(-30, -15, 30, 15);
  Ensemble ensemble("ensemble/rect_60x30_nobox.off", domain);
  for (int f = 0; f < Nf; f++) {
    for (int s = 0; s < Ns; s++) {
      ensemble.step();
    }
    ensemble.output(f);
  }
  ensemble.writeOFF("rect_60x30_nobox_relaxed.off");
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