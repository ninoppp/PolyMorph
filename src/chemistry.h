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

  double get_positional_error() {
    return 0;
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