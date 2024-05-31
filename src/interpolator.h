#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "solver.h"
#include "ensemble.h"
#include "grid.h"
#include "utils.h"

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

  // returns the interpolated velocity at a grid point inside a parent polygon
  Point interior_velocity(const Point& grid_point, const int parent_idx) {
    const Polygon& parent = ensemble.polygons[parent_idx];
    std::vector<double> weights;
    double total_weight = 0;
    for (int i = 0; i < parent.vertices.size(); i++) {
      const Vertex& vertex = parent.vertices[i];
      double distance = (vertex.r - grid_point).length();
      double weight = 1.0 / (distance + 1e-6 * h); // avoid division by zero TODO: scale epsilon with length scale
      weights.push_back(weight);
      total_weight += weight;
    }    
    Point vel = Point(0, 0);
    for (int i = 0; i < parent.vertices.size(); i++) {
      const Vertex& vertex = parent.vertices[i];
      vel = vel + weights[i] / total_weight * vertex.v;
    }
    return vel;
  }

  // ToDo: try with 8 directions
  Point bilinear_vel_interpolation(int i, int j) {
    // find first non-background node or boundary node in each direction
    int i_left = i - 1;
    while (i_left > 0 && solver.parent_idx(i_left, j) < Nr) {
      i_left--;
    }
    int i_right = i + 1;
    while (i_right < solver.Nx-1 && solver.parent_idx(i_right, j) < Nr) {
      i_right++;
    }
    int j_down = j - 1;
    while (j_down > 0 && solver.parent_idx(i, j_down) < Nr) {
      j_down--;
    }
    int j_up = j + 1;
    while (j_up < solver.Ny-1 && solver.parent_idx(i, j_up) < Nr) {
      j_up++;
    }
    // calculate distances
    const double dx_left = i - i_left;
    const double dx_right = i_right - i;
    const double dy_down = j - j_down;
    const double dy_up = j_up - j;
    // calculate weights
    const double w_left = dx_right / (dx_left + dx_right);
    const double w_right = dx_left / (dx_left + dx_right);
    const double w_down = dy_up / (dy_down + dy_up);
    const double w_up = dy_down / (dy_down + dy_up);
    // interpolate velocity
    Point vel = w_left * solver.velocity(i_left, j) 
              + w_right * solver.velocity(i_right, j) 
              + w_down * solver.velocity(i, j_down) 
              + w_up * solver.velocity(i, j_up);
    return vel;
  }

  // scatter coefficients D, k from polygons to grid points
  void scatter() { // ToDo: make this function prettier
    Grid<int>& prev_idx = solver.parent_idx; // stores the polygon index of the cell in which a grid point lies (its parent)
    Grid<int> new_idx(solver.Nx, solver.Ny, -1); // negative indices indicate a background node. ToDo: could make this in place
    
    istart = std::max(int((ensemble.x0 - solver.domain.x0) / solver.dx) + 1, 0);
    jstart = std::max(int((ensemble.y0 - solver.domain.y0) / solver.dx) + 1, 0);
    iend = std::min(int((ensemble.x1 - solver.domain.x0) / solver.dx), solver.Nx);
    jend = std::min(int((ensemble.y1 - solver.domain.y0) / solver.dx), solver.Ny); // plus 1 maybe problem
    #if DEBUG 
    assert(solver.x0 + istart * solver.dx >= ensemble.x0);
    assert(solver.y0 + jstart * solver.dx >= ensemble.y0);
    assert(solver.x0 + iend * solver.dx <= ensemble.x1);
    assert(solver.y0 + jend * solver.dx <= ensemble.y1);
    #endif
    // ToDo: print warning if ensemble box greater than solver box (only once)

    #pragma omp parallel for collapse(2)
    for (int i = istart; i < iend; i++) {
      for (int j = jstart; j < jend; j++) { 
        // spatial coordinates of grid point
        const double x = solver.domain.x0 + i * solver.dx;
        const double y = solver.domain.y0 + j * solver.dx;
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
          // set velocity to zero ToDo: remove duplication (done below in vel-field interpolation)
          if (ADVECTION_DILUTION) {
            //solver.velocity(i, j) = Point(0, 0);
          }
        } else { 
          solver.D(i, j) = ensemble.polygons[new_idx(i, j)].D;
          solver.k(i, j) = ensemble.polygons[new_idx(i, j)].k;
          solver.p(i, j) = ensemble.polygons[new_idx(i, j)].p;
          // interpolate velocity
          if (ADVECTION_DILUTION) {
            solver.velocity(i, j) = interior_velocity(grid_point, new_idx(i, j));
          }
        }
      }
    }
    solver.parent_idx = new_idx;  // ToDo: update in place

    // interpolate remaining velocity field
    if (ADVECTION_DILUTION) {
      // set boundary to domain velocity
      for (int i = 0; i < solver.Nx; i++) {
        solver.velocity(i, 0) = solver.domain.growth_rate[3] * Point(0, -1); // south
        solver.velocity(i, solver.Ny - 1) = solver.domain.growth_rate[1] * Point(0, 1); // north
      }
      for (int j = 0; j < solver.Ny; j++) {
        solver.velocity(0, j) = solver.domain.growth_rate[2] * Point(-1, 0); // west
        solver.velocity(solver.Nx - 1, j) = solver.domain.growth_rate[0] * Point(1, 0); // east
      }
      // interior points
      #pragma omp parallel for collapse(2)
      for (int i = 1; i < solver.Nx - 1; i++) {
        for (int j = 1; j < solver.Ny - 1; j++) {
          if (solver.parent_idx(i, j) < int(Nr)) { // only treat background nodes
            solver.velocity(i, j) = bilinear_vel_interpolation(i, j);
          } 
        }
      }
    }
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
      // store gradient at vertices
      for (auto& vertex : cell.vertices) {
        const int i = std::round((vertex.r.x - solver.domain.x0) / solver.dx);
        const int j = std::round((vertex.r.y - solver.domain.y0) / solver.dx);
        if (i >= 0 && i < solver.Nx && j >= 0 && j < solver.Ny) {
          vertex.grad_u = solver.grad_u(i, j);
        }
      }
    }
  }
};

#endif