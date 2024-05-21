#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include <numeric>

#include "polyhoop.h"
#include "utils.h"

// Handles polygon modification due to signaling effects (growth control, chemotaxis, etc.)
struct Chemistry {
  Ensemble& ensemble; // ToDo: don't store ens and solver but pass to functions. move some functions to ensemble and solver
  Solver& solver;
  bool growth_control = false; // set true to stop growth when below threshold
  std::function<std::vector<bool>(Polygon&)> is_producing = [](Polygon& p) { return std::vector(NUM_SPECIES, false); };
  std::function<bool(Polygon&)> flag = [](Polygon& p) { return p.u[0] < p.threshold[0]; };
  Chemistry(Ensemble& ensemble, Solver& solver) : ensemble(ensemble), solver(solver) {}

  void update() {
    #pragma omp parallel for
    for (int i = Nr; i < ensemble.polygons.size(); i++) {
      auto& cell = ensemble.polygons[i];
      // set production. ToDo: move to ensemble
      std::vector<bool> producing = is_producing(cell); // which species are produced by the cell
      if (producing.size() != NUM_SPECIES) {
        std::cerr << "is_producing function must return a vector of size NUM_SPECIES" << std::endl;
        exit(1);
      }
      std::vector<double> production_rate = sample(p_dist, rng);
      for (int i = 0; i < NUM_SPECIES; i++) {
        if (cell.p[i] == 0 && producing[i]) {
          cell.p[i] = production_rate[i];
        } else if (cell.p[i] > 0 && !producing[i]) {
          cell.p[i] = 0;
        }
      }

      cell.flag = flag(cell); // set flag ToDo: move to interpolator or ensemble   
      if (growth_control && cell.flag) cell.alpha = 0;
    }
  }  

  // currently assumes left to right gradient
  double get_precision_zone_width() {
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

  // returns the standard deviation of the x-coordinates of the border vertices
  std::pair<double, double> get_positional_error() { // PROBLEM: This doesn't work. Also counts right end cells plus cells beyond the border. FIX? use fact that adjecent vertex very close by
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

    if (border_x.size() == 0) return {0, 0};
    const double mean_border_x = std::accumulate(border_x.begin(), border_x.end(), 0.0) / border_x.size();
    double variance = 0;
    for (auto x : border_x) {
      variance += (x - mean_border_x) * (x - mean_border_x);
    }
    double std_dev = std::sqrt(variance / border_x.size());
    return {mean_border_x, std_dev};
  }
};

#endif