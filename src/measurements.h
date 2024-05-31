#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include <numeric>

#include "ensemble.h"
#include "utils.h"

namespace Measurements {
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

  void apply_flag(Ensemble& ensemble) {
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      auto& cell = ensemble.polygons[p];
      cell.flag = set_flag(cell);
    }
  }
};

#endif