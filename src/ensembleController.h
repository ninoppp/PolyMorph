#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include <numeric>
#include "ensemble.h"
#include "utils.h"
#include "solver.h"
#include "interpolator.h"

/*! \file ensembleController.h
 *  \brief Facade-like routines to interact with ensemble during experiments
 */

namespace EnsembleController {

  // currently assumes left to right gradient
  double get_precision_zone_width(Ensemble& ensemble) {
    double xmin = ensemble.x1;
    double xmax = ensemble.x0;
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      auto& cell = ensemble.polygons[p];
      for (auto& vertex : cell.vertices) {
        if (cell.flag) {
          xmin = std::min(xmin, vertex.r.x); // leftmost flagged point
        } else {
          xmax = std::max(xmax, vertex.r.x); // rightmost unflagged
        }
      }
    }
    std::cout << "xmin=" << xmin << " xmax=" << xmax << std::endl;
    return xmax - xmin;
  }

  double mean_readout_position(Ensemble& ensemble, Solver& solver) {
    std::vector<double> border_x;
    // find first flagged grid point
    for (int j = 0; j < solver.Ny; j++) {
      for (int i = 0; i < solver.Nx; i++) {
        if (solver.parent_idx(i, j) >= Nr && ensemble.polygons[solver.parent_idx(i, j)].flag) {
          border_x.push_back(solver.domain.x0 + i * solver.dx);
          break;
        }
      }
    }
    const double mean_border_x = std::accumulate(border_x.begin(), border_x.end(), 0.0) / border_x.size();
    return mean_border_x;
  }

  // assign all flags without need for solver step
  void apply_flag(Ensemble& ensemble) {
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      auto& cell = ensemble.polygons[p];
      cell.flag = ensemble.set_flag(cell);
    }
  }

  // generate tissue filling out the given domain
  // domain should center at (0,0) to produce an uniform tissue
  Ensemble grow_tissue(Domain& domain) {
    if (beta != 0.8) {
      std::cout << "Warning: beta=" << beta << ". It's recommended to grow with beta=0.8" << std::endl;
    }
    if (kh != 0) {
      std::cout << "Warning: kh=" << kh << ". It's recommended to grow without adhesion" << std::endl;
    }
    Ensemble ensemble("ensemble/default.off", domain);
    Solver solver(domain, dx, Reactions::linearDegradation);
    Interpolator interpolator(ensemble, solver);
    int max = INT32_MAX;
    int f = 1;
    while (f < max) {
      for (int i = 0; i < Ns; i++) {
        ensemble.step();
      }
      interpolator.scatter();
      //ensemble.output(f);
      // check if all 4 corners are occupied by a polygon
      if (max == INT32_MAX &&
          solver.parent_idx(1, 1) >= Nr &&
          solver.parent_idx(solver.Nx - 2, solver.Ny - 2) >= Nr &&
          solver.parent_idx(1, solver.Ny - 2) >= Nr &&
          solver.parent_idx(solver.Nx - 2, 1) >= Nr) {
        std::cout << "All corners filled. Relaxing now" << std::endl;
        max = f + f/3; // run again for 1/4 of the time to relax the tissue
      }
      f++;
    }
    return ensemble;
  }

  // generate tissue with a given number of polygons
  Ensemble grow_tissue(int num_polygons) {
    double polygon_diameter = 2 * std::sqrt(Amax_mu / M_PI);
    double L = 2 * std::sqrt(num_polygons) * polygon_diameter;
    Domain domain(-L/2, -L/2, L/2, L/2);
    Ensemble ensemble("ensemble/default.off", domain);
    while (ensemble.polygons.size() < num_polygons) {
      ensemble.step();
    }
    ensemble.output(0);
    return ensemble;
  }
};

#endif