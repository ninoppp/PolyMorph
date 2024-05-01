#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "solver.h"
#include "polyhoop.h"

// takes care of the data scattering and gathering between ensemble and solver
struct Interpolator {
  Ensemble& ensemble;
  Solver& solver;
  int istart, jstart, iend, jend; // limits of the grid points to be updated (within ensemble bounds)
  Interpolator(Ensemble& ensemble, Solver& solver) : ensemble(ensemble), solver(solver) {}
  
  // Search algorithm to find parent polygon for a grid point.
  // Complexity O(vertices per polygon)
  std::size_t find_parent(Point grid_point) {
    std::size_t bxi = (grid_point.x - ensemble.x0) / ensemble.bs + 1; // box index in x direction
    const std::size_t byi = (grid_point.y - ensemble.y0) / ensemble.bs + 1; // box index in y direction
    std::unordered_set<std::size_t> checked_polygons; // store checked polygons to avoid checking them again
    bool last_iteration = false; // abort search as soon as we can be sure it's a background node
    while (bxi < ensemble.Nx) {
      for (Vertex* v = ensemble.first[bxi * ensemble.Ny + byi]; v; v = v->next) { // loop over vertices in box
        if (checked_polygons.find(v->p) == checked_polygons.end()) {  // new polygon encountered
          if (v->p >= Nr && ensemble.polygons[v->p].contains(grid_point)) { // don't want to check rigid polygons
            return v->p;
          } else {
            checked_polygons.insert(v->p);
          }
        }
      }
      if (last_iteration) return -2; // background node
      if (!checked_polygons.empty()) last_iteration = true; // only go 1 more layer
      ++bxi;
    }
    return -2; // background node (reached boundary of ensemble box)
  }

  // scatter coefficients D, k from polygons to grid points
  void scatter() { // ToDo: make this function prettier
    Grid<int>& prev_idx = solver.parent_idx; // stores the polygon index of the cell in which a grid point lies (its parent)
    Grid<int> new_idx(solver.Nx, solver.Ny, -1); // negative indices indicate a background node. ToDo: could make this in place
    
    istart = std::max(int((ensemble.x0 - solver.x0) / solver.dx) + 1, 0);
    jstart = std::max(int((ensemble.y0 - solver.y0) / solver.dx) + 1, 0);
    iend = std::min(size_t((ensemble.x1 - solver.x0) / solver.dx), solver.Nx);
    jend = std::min(size_t((ensemble.y1 - solver.y0) / solver.dx), solver.Ny); // plus 1 maybe problem
    /*assert(solver.x0 + istart * solver.dx >= ensemble.x0);
    assert(solver.y0 + jstart * solver.dx >= ensemble.y0);
    assert(solver.x0 + iend * solver.dx <= ensemble.x1);
    assert(solver.y0 + jend * solver.dx <= ensemble.y1);*/

    #pragma omp parallel for collapse(2)
    for (int i = istart; i < iend; i++) {
      for (int j = jstart; j < jend; j++) { 
        // spatial coordinates of grid point
        const double x = solver.x0 + i * solver.dx;
        const double y = solver.y0 + j * solver.dx;
        const Point grid_point(x, y);
        // check if still the same parent (ignore rigid)
        if (prev_idx(i, j) >= int(Nr) && ensemble.polygons[prev_idx(i, j)].contains(grid_point)) { // ToDo: check if in bounds
          new_idx(i, j) = prev_idx(i, j);
        } 
        else {
          new_idx(i, j) = find_parent(grid_point);
        }
        // scatter values
        if (new_idx(i, j) < int(Nr)) { // is background node (treat rigid as BG)
          solver.D(i, j) = D0; // background diffusion
          solver.k(i, j) = k0; // background degradation
          solver.p(i, j) = p0; // background production (should be zero)
        } else { 
          solver.D(i, j) = ensemble.polygons[new_idx(i, j)].D;
          solver.k(i, j) = ensemble.polygons[new_idx(i, j)].k;
          solver.p(i, j) = ensemble.polygons[new_idx(i, j)].p;
        }
      }
    }
    solver.parent_idx = new_idx;  // ToDo: update in place
  }

  // gather concentration u from grid points to polygons
  // important: depends on scatter being called every iteration to build the parent_idx
  void gather() {
    // get all children from parent idx built during scatter()
    for (auto& cell : ensemble.polygons) {
      cell.children.clear();
    }
    for (int i = istart; i < iend; i++) {
      for (int j = jstart; j < jend; j++) {
        if (solver.parent_idx(i, j) >= int(Nr)) { // skip background nodes
          auto& cell = ensemble.polygons[solver.parent_idx(i, j)]; 
          cell.children.push_back(Index(i, j)); // cannot parallelize this part
        }
      }
    }
    // accumulate data from children (average concentration)
    #pragma omp parallel for
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      auto& cell = ensemble.polygons[p];
      cell.u = std::vector<double>(NUM_SPECIES, 0.0);
      for (const Index& idx : cell.children) {
        for (int i = 0; i < NUM_SPECIES; i++){
          cell.u[i] += solver.u(idx)[i];
        }
      }
      if (cell.children.size() > 0) { // avoid division by zero if cells exceed RD box
        for (int i = 0; i < NUM_SPECIES; i++) {
          cell.u[i] /= cell.children.size();
        }
      }
    }
  }
};

#endif