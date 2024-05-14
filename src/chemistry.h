#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include <numeric>

#include "polyhoop.h"
#include "utils.h"

// Handles polygon modification due to signaling effects (growth control, chemotaxis, etc.)
struct Chemistry {
  Ensemble& ensemble;
  bool growth_control = false; // set true to stop growth when below threshold
  std::function<std::vector<bool>(Polygon&)> is_producing = [](Polygon& p) { return std::vector(NUM_SPECIES, false); };

  Chemistry(Ensemble& ensemble) : ensemble(ensemble) {}

  void update() { // ToDo: split into flag() and produce()
    #pragma omp parallel for
    for (int i = Nr; i < ensemble.polygons.size(); i++) {
      auto& cell = ensemble.polygons[i];
      // set production
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

      // flag below threshold
      if (!cell.flag && cell.u[0] < cell.threshold[0]) {  // atm only consider first species
        cell.flag = true;
        if (growth_control) cell.alpha = 0;
      } else if (cell.flag && cell.u[0] > cell.threshold[0]){
        cell.flag = false;
        if (growth_control) cell.alpha = cell.alpha0; 
      }
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
  double get_positional_error() { // PROBLEM: This doesn't work. Also counts right end cells plus cells beyond the border. FIX? use fact that adjecent vertex very close by
    double mean_border_x; 
    std::vector<double> border_x;
    // find all vertices that lie at the border
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      const auto& cell = ensemble.polygons[p];
      if (cell.flag) {
        for (auto& vertex : cell.vertices) {
          const std::size_t bxi = (vertex.r.x - ensemble.x0) / ensemble.bs + 1; // box index in x direction
          const std::size_t byi = (vertex.r.y - ensemble.y0) / ensemble.bs + 1; // box index in y direction
          for (size_t bxj = bxi - 1; bxj <= bxi; bxj++) { // check this box and one to the left for neighbouring vertices
            for (Vertex* v = ensemble.first[bxi * ensemble.Ny + byi]; v; v = v->next) { // loop over vertices in box
              if (vertex.r < v->r && !ensemble.polygons[v->p].flag) { // if vertex is to the left and in an unflagged cell
                border_x.push_back(vertex.r.x);
              }
            }
          }
        }
      }
    }
    if (border_x.size() == 0) return 0;
    mean_border_x = std::accumulate(border_x.begin(), border_x.end(), 0.0) / border_x.size();
    double variance = 0;
    for (auto x : border_x) {
      variance += (x - mean_border_x) * (x - mean_border_x);
    }
    double std_dev = std::sqrt(variance / border_x.size());
    return std_dev;
  }

  void chemotaxis() {
    #pragma omp parallel for
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      auto& cell = ensemble.polygons[p];
      for (auto& vertex : cell.vertices) {
        // ToDo
        // move towards higher concentration
        // manipulate vertex acceleration
        // Q: how to get morphogen gradient cheaply?
        ; 
      }
    }
  }
};

#endif