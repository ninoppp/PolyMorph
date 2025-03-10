#ifndef ENSEMBLE_CONTROLLER_H
#define ENSEMBLE_CONTROLLER_H

#include <numeric>

#include "interpolator.h"
#include "ensemble.h"
#include "solver.h"

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
        if (cell.cell_type) {
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
        if (solver.parent_idx(i, j) >= Nr && ensemble.polygons[solver.parent_idx(i, j)].cell_type) {
          border_x.push_back(solver.domain.x0 + i * solver.dx);
          break;
        }
      }
    }
    const double mean_border_x = std::accumulate(border_x.begin(), border_x.end(), 0.0) / border_x.size();
    return mean_border_x;
  }

  // NOT TESTED YET. 
  // readout positions for each cell_type (corresponding to a different threshold)
  std::vector<double> mean_readout_positions(Ensemble& ensemble, Solver& solver) {
    int num_readouts = threshold_mu.size();
    std::vector<double> mean_border_x(num_readouts);
    for (int cell_type = 1; cell_type <= num_readouts; cell_type++) {
      std::vector<double> border_x;
      // find first flagged grid point
      for (int j = 0; j < solver.Ny; j++) {
        for (int i = 0; i < solver.Nx; i++) {
          if (solver.parent_idx(i, j) >= Nr && ensemble.polygons[solver.parent_idx(i, j)].cell_type == cell_type) {
            border_x.push_back(solver.domain.x0 + i * solver.dx);
            break;
          }
        }
      }
      const double mean = std::accumulate(border_x.begin(), border_x.end(), 0.0) / border_x.size();
      mean_border_x[cell_type - 1] = mean;
    }
    return mean_border_x;
  }

  // assign all flags without need for solver step
  void apply_flag(Ensemble& ensemble) {
    #pragma omp parallel for
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      auto& cell = ensemble.polygons[p];
      cell.cell_type = ensemble.cellTypeEffect(cell, cell.c, cell.grad_c, 0);
    }
  }

  // stops proliferation completely, both growth and cell division
  void stop_growth(Ensemble& ensemble) {
    #pragma omp parallel for
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      ensemble.polygons[p].alpha = 0;
      ensemble.polygons[p].Amax = MAXFLOAT; 
    }
  }

  // void stop_growth_if_flagged(Ensemble& ensemble) {
  //   #pragma omp parallel for
  //   for (int p = Nr; p < ensemble.polygons.size(); p++) {
  //     if (ensemble.polygons[p].cell_type) {
  //       ensemble.polygons[p].alpha = 0;
  //       ensemble.polygons[p].Amax = MAXFLOAT; 
  //     } else if (!ensemble.polygons[p].cell_type && ensemble.polygons[p].alpha == 0) {
  //       ensemble.polygons[p].alpha = ensemble.polygons[p].alpha;
  //       ensemble.polygons[p].Amax = sample(ensemble.Amax_dist, ensemble.rng);
  //     }
  //   }
  // }

  // only stops cell divisions (keeps growing) ToDo: fix instability with exploding cells
  void stop_division(Ensemble& ensemble) {
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      // ToDo: test if this is sufficient or if additional condition has to be added in ensemble code
      ensemble.polygons[p].Amax = MAXFLOAT; 
    }
  }

  // generate tissue filling out the given domain
  // domain should center at (0,0) to produce a uniform tissue
  void grow_tissue(Ensemble& ensemble, bool output_frames = false) {
    if (beta != 0.8) {
      std::cout << "Warning: beta=" << beta << ". It's recommended to grow with beta=0.8" << std::endl;
    }
    if (kh != 0) {
      std::cout << "Warning: kh=" << kh << ". It's recommended to grow without adhesion" << std::endl;
    }
    Solver solver(ensemble.domain, dx, Reactions::linearDegradation);
    Interpolator interpolator(ensemble, solver);
    int max = INT32_MAX;
    int f = 1;
    while (f < max) {
      for (int i = 0; i < Ns; i++) {
        ensemble.step();
      }
      interpolator.scatter();
      if (output_frames) {
        ensemble.output(f);
      }
      // check if all 4 corners are occupied by a polygon
      if (max == INT32_MAX &&
          solver.parent_idx(1, 1) >= Nr &&
          solver.parent_idx(solver.Nx - 2, solver.Ny - 2) >= Nr &&
          solver.parent_idx(1, solver.Ny - 2) >= Nr &&
          solver.parent_idx(solver.Nx - 2, 1) >= Nr) {
        std::cout << "All corners filled. Relaxing now" << std::endl;
        stop_growth(ensemble);// change to stop cell division
        max = f + f/3; // run again for 1/3 of the time to relax the tissue
      }
      f++;
    }
  }

  // generate tissue with a given number of polygons
  Ensemble grow_tissue_by_num(int num_polygons) {
    double polygon_diameter = 2 * std::sqrt(Amax_mu / M_PI);
    double L = 4 * std::sqrt(num_polygons) * polygon_diameter;
    Domain domain(-L/2, -L/2, L/2, L/2);
    Ensemble ensemble("ensemble/default.off", domain);
    size_t steps = 0;
    while (ensemble.polygons.size() < num_polygons) {
      ensemble.step();
      steps++;
    }
    // relaxing: run again for 1/3 of the time to relax the tissue
    stop_growth(ensemble);
    for (int i = 0; i < steps/3; i++) {
      ensemble.step();
    }
    ensemble.output(0);
    return ensemble;
  }  

  // regenerate all distribution-drawn parameters
  // needed when manipulating dists after ensemble initialization
  void redraw_params_from_dists(Ensemble& ensemble) {
    #pragma omp parallel for
    for (int p = 0; p < ensemble.polygons.size(); p++) {
      auto& cell = ensemble.polygons[p];
      // vectors
      cell.D = sample(ensemble.D_dist, ensemble.rng, true);
      cell.p = sample(ensemble.p_dist, ensemble.rng);
      cell.k = sample(ensemble.k_dist, ensemble.rng);
      cell.threshold = sample(ensemble.threshold_dist, ensemble.rng);
      // scalars
      cell.Amax = ensemble.Amax_dist(ensemble.rng); // TODO: use sample function
      cell.alpha = ensemble.alpha_dist(ensemble.rng);
    }
  }

  // avg number of children. 
  // only works when calling "gather" before this
  double get_avg_gridpoints_per_polygon(Ensemble& ensemble) {
    int sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      sum += ensemble.polygons[p].children.size();
    }
    return sum / (ensemble.polygons.size() - Nr);
  }
};

#endif