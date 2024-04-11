// PolyHoop
// Copyright (c) 2023 Roman Vetter
// vetterro@ethz.ch
// ETH Zurich

#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <algorithm>
#include <cassert>
#include <unordered_set>

#include "solver.h"

constexpr double h = 0.01; // [L] edge thickness
constexpr double lmin = 0.02; // [L] minimum edge length
constexpr double lmax = 0.2; // [L] maximum edge length
constexpr double Q = 1; // [-] isoparametric ratio
constexpr double alpha_mu = 1; // [L^2/T] mean area growth rate
constexpr double alpha_CV = 0; // [-] coefficient of variation of area growth rate
constexpr double Amin = 0; // [L^2] minimum area
constexpr double Amax_mu = M_PI; // [L^2] mean maximum area
constexpr double Amax_CV = 0; // [-] coefficient of variation of maximum area
constexpr double gam = 1e3; // [L/T^2] line tension per vertex mass
constexpr double ka = 1e5; // [1/(L^2*T^2)] area stiffness per vertex mass
constexpr double kl = 1e4; // [L/T^2] edge contractility stiffness per vertex mass
constexpr double kb = 0; // [L^3/T^2] bending stiffness per vertex mass
constexpr double kr = 1e7; // [1/T^2] repulsion stiffness per vertex mass
constexpr double kh = 1e6; // [1/T^2] adhesion stiffness per vertex mass
constexpr double sh = 0.01; // [L] adhesion hardening zone size
constexpr double ss = 0.01; // [L] adhesion softening zone size
constexpr double theta = 0; // [-] fusion threshold
constexpr double mu = 0; // [-] dynamic friction coefficient
constexpr double rho = 0; // [1/L^2] fluid mass density per vertex mass
constexpr double g = 0; // [L/T^2] gravitational acceleration
constexpr double gl = 0; // [L/T^2] edge gravitational acceleration
constexpr double cv = 10; // [1/T] viscous damping rate
constexpr double cd = 0; // [-] drag coefficient
constexpr double cc = 30; // [1/T] collision damping rate
constexpr double dt = 1e-4; // [T] time step

constexpr std::size_t Nf = 100; // number of output frames
constexpr std::size_t Ns = 1000; // number of time steps between frames // default 1000
constexpr std::size_t Nr = 1; // number of rigid polygons

constexpr double drmax = h + sh + ss; // maximum interaction distance

// PolyMorh extension
constexpr double dx = 0.1; // [L] grid spacing for solver
constexpr double D_mu = 3.0; // [L^2/T] diffusion coefficient mean
constexpr double k_mu = 1.0; // [1/T] degradation rate mean 
constexpr double p_mu = 6.0; // [1/T] production rate mean
constexpr double D_CV = 0.1; // [-] coefficient of variation of diffusion
constexpr double k_CV = 0.1; // [-] coefficient of variation of degradation rate
constexpr double p_CV = 0.1; // [-] coefficient of variation of production rate
constexpr double D0 = 3.0; // [L^2/T] diffusion coefficient background
constexpr double k0 = 1.0; // [1/T] reaction rate background
constexpr double p0 = 0.0; // [1/T] reaction rate background

std::mt19937 rng; // random number generator
const double Amax_lnCV = std::log(1 + Amax_CV*Amax_CV);
const double alpha_lnCV = std::log(1 + alpha_CV*alpha_CV);
std::lognormal_distribution<> Amax_dist(std::log(Amax_mu) - Amax_lnCV/2, std::sqrt(Amax_lnCV)); // division area distribution
std::lognormal_distribution<> alpha_dist(std::log(alpha_mu) - alpha_lnCV/2, std::sqrt(alpha_lnCV)); // area growth rate distribution
std::uniform_real_distribution<> uni_dist;

// PolyMorph extension
const double D_lnCV = std::log(1 + D_CV*D_CV);
const double k_lnCV = std::log(1 + k_CV*k_CV);
const double p_lnCV = std::log(1 + p_CV*p_CV);
std::lognormal_distribution<> D_dist(std::log(D_mu) - D_lnCV/2, std::sqrt(D_lnCV)); // diffusion coefficient distribution ToDo: "magic numer sigma"
std::lognormal_distribution<> k_dist(std::log(k_mu) - k_lnCV/2, std::sqrt(k_lnCV)); // reaction rate distribution
std::lognormal_distribution<> p_dist(std::log(p_mu) - p_lnCV/2, std::sqrt(p_lnCV)); // production rate distribution

struct Point  // basically a 2D vector
{
  double x, y; // coordinates
  Point () : x(0), y(0) {} // default constructor
  Point (const double x, const double y) : x(x), y(y) {} // constructor
  Point operator+(const Point& r) const { return {x + r.x, y + r.y}; } // vector sum
  Point operator-(const Point& r) const { return {x - r.x, y - r.y}; } // vector difference
  Point cross() const { return {-y, x}; } // perpendicular vector
  double operator*(const Point& r) const { return x * r.x + y * r.y; } // dot product
  double operator^(const Point& r) const { return x * r.y - y * r.x; } // wedge product
  double length2() const { return x * x + y * y; } // squared norm
  double length() const { return std::sqrt(length2()); } // norm
  void add(const double a, const Point& r)
  {
    #pragma omp atomic
    x += a * r.x;
    #pragma omp atomic
    y += a * r.y;
  }
};
Point operator*(const double a, const Point& r) { return {a * r.x, a * r.y}; }

double point_edge_dist2(const Point& r0, const Point& r1, const Point& r2, double& xi, double& xit, Point& dr)
{
  const Point r12 = r2 - r1;
  const Point r10 = r0 - r1;
  xi = r10 * r12 / (r12 * r12);
  xit = std::min(std::max(xi, 0.), 1.);
  dr = r10 - xit * r12;
  return dr * dr;
}

struct Vertex
{
  Point r, v, a; // position, velocity, acceleration
  std::size_t p; // polygon index
  Vertex* next; // pointer to next vertex in same box
  double l0; // rest length of edge to the right
};

struct Polygon
{
  std::vector<Vertex> vertices; // vertex list in counter-clockwise orientation
  bool phase; // phase of the enclosed medium // NM: can extend this to label source cells
  double A0, A, Amax, alpha; // target, actual & division area, area growth rate
  // Polymorph extension. 
  double D, k, p; // diffusion coefficient, degradation rate, production rate
  double u = 0; // local morphogen concentration
  std::vector<Index> children; // stores the indices of the grid points within the polygon
  // END Polymorph extension
  double area()
  {
    A = 0;
    for (std::size_t i = vertices.size() - 1, j = 0; j < vertices.size(); i = j++)
      A += vertices[i].r ^ vertices[j].r;
    return A /= 2;
  }
  double perimeter() const
  {
    double L = 0;
    for (std::size_t i = vertices.size() - 1, j = 0; j < vertices.size(); i = j++)
      L += (vertices[j].r - vertices[i].r).length();
    return L;
  }
  double perimeter0() const
  {
    double L0 = 0;
    for (auto& v : vertices)
      L0 += v.l0;
    return L0;
  }
  bool contains(const Point& r) const
  {
    bool in = false;
    for (std::size_t i = vertices.size() - 1, j = 0; j < vertices.size(); i = j++)
      if ((vertices[i].r.y > r.y) != (vertices[j].r.y > r.y) &&
          r.x < (vertices[i].r.x - vertices[j].r.x) * (r.y - vertices[j].r.y) / (vertices[i].r.y - vertices[j].r.y) + vertices[j].r.x)
        in = !in;
    return in; // true if the point lies inside this polygon
  }
  // Polymorph extension. Calculate (geometric?) midpoint of polygon
  Point midpoint() const 
  {
    Point midp;
    for (const Vertex& v : vertices) {
      midp = midp + v.r;
    }
    return 1/vertices.size() * midp;
  }
};

// methods for selecting producing cells ToDo: move this to a better place
auto heavyside = [](const Polygon& p) { return p.midpoint().x < 0.5; }; // ToDo: remove magic-number 0.5
auto mother_cell = [](const Polygon& p) { return p.vertices[0].p == Nr; }; // workaround to get polygon index
// choose method
auto is_producing = mother_cell;

struct Ensemble
{
  std::vector<Polygon> polygons; // list of polygons
  std::vector<Vertex*> first; // pointer to first vertex in each box
  std::size_t Nx, Ny; // number of boxes in x,y direction
  double x0, y0, x1, y1, bs; // box offset, box end, and size
  double t; // time
  
  Ensemble(const char* name) : t(0)
  {
    // read OFF file header
    std::ifstream file(name);
    file.ignore(3); // ignore "OFF"
    std::size_t Nv, Np, Ne; // number of vertices, polygons, edges
    file >> Nv >> Np >> Ne; // Ne unused
    
    // read all vertices
    std::vector<Point> points(Nv);
    std::vector<int> z(Nv);
    for (std::size_t i = 0; i < Nv; ++i)
      file >> points[i].x >> points[i].y >> z[i];
    
    // read all polygons
    polygons.resize(Np);
    for (std::size_t p = 0, j; p < Np; ++p)
    {
      file >> Nv;
      for (std::size_t i = 0; i < Nv; ++i)
      {
        file >> j;
        polygons[p].vertices.push_back({points[j], {0, 0}, {0, 0}, p});
      }
      polygons[p].phase = std::abs(z[j]) % 2; // the z coordinate is used as the phase 
      polygons[p].A0 = polygons[p].area(); // use the initial area as target area
      polygons[p].Amax = Amax_dist(rng);
      polygons[p].alpha = alpha_dist(rng);
      // PolyMorph extension
      polygons[p].D = D_dist(rng);
      polygons[p].k = k_dist(rng);
      polygons[p].p = is_producing(polygons[p]) ? p_dist(rng) : 0;
      // end PolyMorph extension
      for (std::size_t i = Nv - 1, j = 0; j < Nv; i = j++)
        polygons[p].vertices[i].l0 = (polygons[p].vertices[j].r - polygons[p].vertices[i].r).length(); // edge rest length
    }
  }
  
  void remove(std::size_t p)
  {
    polygons[p] = polygons.back(); // copy the last polygon to index p
    for (auto& v : polygons[p].vertices)
      v.p = p; // let the vertices know about their new position
    polygons.pop_back(); // remove the last polygon
  }
  
  void step()
  {
    // polygon fusion
    while (theta > 0) // this is either if(false) or while(true)
    {
      volatile bool fuse = false;
      std::size_t prem; // id of the polygon to be removed
      Polygon pnew [2]; // new polygons to be added (either one or two)
      boxes(); // place all vertices into boxes
      #pragma omp parallel for
      for (std::size_t p1 = Nr; p1 < polygons.size(); ++p1)
      {
        auto& v1 = polygons[p1].vertices;
        for (std::size_t i = v1.size() - 1, k = i - 1, j = 0; !fuse && j < v1.size(); k = i, i = j++)
        {
          const std::size_t bxi = (v1[i].r.x - x0) / bs + 1; // box index in x direction
          const std::size_t byi = (v1[i].r.y - y0) / bs + 1; // box index in y direction
          for (std::size_t bxj = bxi - 1; bxj <= bxi + 1; ++bxj)
            for (std::size_t byj = byi - 1; byj <= byi + 1; ++byj)
              for (Vertex* v = first[bxj * Ny + byj]; v; v = v->next)
                if (v->p >= Nr && v != &v1[k] && v != &v1[i] && v != &v1[j])
                {
                  // select the closer of the two edges
                  Vertex* vn [2] = {&v1[k], &v1[j]}; // pointers to neighbors of vertex i
                  double xi [2], xit, dr2 [2];
                  Point dr;
                  for (unsigned int a = 0; a < 2; ++a)
                    dr2[a] = point_edge_dist2(v->r, v1[i].r, vn[a]->r, xi[a], xit, dr);
                  const unsigned int a = dr2[0] > dr2[1];
                  
                  if (xi[a] <= 1 && dr2[a] < theta * h * theta * h)
                  {
                    // determine if one of the polygons is inside the other
                    std::size_t p2 = v->p;
                    bool inside [2] = {p1 != p2 && polygons[p2].contains(v1[i].r), p1 != p2 && polygons[p1].contains(v->r)};
                    const bool in = inside[0] || inside[1];
                    
                    // walk along both polygons in both directions to find end vertices that are separated enough
                    std::size_t p [2] = {p1, p2};
                    std::size_t vend [4] = {i, i, 0, 0}; // end vertices number 0 to 3
                    auto& v2 = polygons[p2].vertices;
                    vend[2] = vend[3] = v - v2.data();
                    bool valid = true, close = true;
                    for (unsigned int b = 0; valid && close; b = 1 - b)
                      for (unsigned int d = (a + b * !in) % 2, c = 0; valid && close && c < 2; d = 1 - d, ++c)
                      {
                        vend[2*b+d] = (vend[2*b+d] + polygons[p[b]].vertices.size() + 2*d-1) % polygons[p[b]].vertices.size();
                        close = point_edge_dist2(v1[vend[0]].r, v1[vend[1]].r, v2[vend[2+in]].r, xi[0], xit, dr) < h * h ||
                                point_edge_dist2(v1[vend[1]].r, v1[vend[0]].r, v2[vend[3-in]].r, xi[0], xit, dr) < h * h ||
                                point_edge_dist2(v2[vend[2]].r, v1[vend[  in]].r, v2[vend[3]].r, xi[0], xit, dr) < h * h ||
                                point_edge_dist2(v2[vend[3]].r, v1[vend[1-in]].r, v2[vend[2]].r, xi[0], xit, dr) < h * h;
                        if (vend[2*b+d] == vend[2*b+!d] || (p1 == p2 && (vend[0] == vend[3] || vend[1] == vend[2])))
                          valid = false; // if any of the end vertices are equal, don't fuse
                      }
                    if (!valid) continue;

                    const double ls1 = polygons[p1].perimeter() / polygons[p1].perimeter0(); // average stretch ratio of polygon 1
                    
                    if (p1 != p2) // fusion type I: merge two polygons into one
                    {
                      # pragma omp critical
                      if (!fuse)
                      {
                        fuse = true;
                        prem = p2; // polygon 2 will be removed
                        
                        // merge the vertex lists, inverting the vertex order for polygons enclosed by others
                        const double ls2 = polygons[p2].perimeter() / polygons[p2].perimeter0(); // average stretch ratio of polygon 2
                        for (unsigned int b = 0; b < 2; ++b)
                        {
                          const bool d = inside[b];
                          for (std::size_t c = vend[2*b+!d]; ; c = (c + polygons[p[b]].vertices.size() + 1-2*d) % polygons[p[b]].vertices.size())
                          {
                            pnew[0].vertices.push_back(polygons[p[b]].vertices[c]);
                            pnew[0].vertices.back().p = p1;
                            if (d && pnew[0].vertices.size() > 1)
                              pnew[0].vertices[pnew[0].vertices.size()-2].l0 = pnew[0].vertices.back().l0; // shift l0 if vertex order is inverted
                            if (c == vend[2*b+d]) break;
                          }
                          // for external polygons, connect vertices 0-3 & 1-2, for polygons inside each other connect vertices 0-2 & 1-3
                          pnew[0].vertices.back().l0 = (v1[vend[b]].r - v2[vend[3-(b!=in)]].r).length() / (ls1 + ls2) * 2;
                        }
                        
                        // adjust polygon properties
                        pnew[0].phase = polygons[p[inside[0]]].phase; // use the phase of the outer polygon
                        pnew[0].A0 = (1-2*inside[0]) * polygons[p1].A0 + (1-2*inside[1]) * polygons[p2].A0;
                        pnew[0].Amax  = !inside[0] * polygons[p1].Amax  + !inside[1] * polygons[p2].Amax;
                        pnew[0].alpha = !inside[0] * polygons[p1].alpha + !inside[1] * polygons[p2].alpha;
                        pnew[0].area(); // compute the new area
                      }
                    }
                    else // fusion type II: split a polygon into two
                    {
                      // check if arclength distance between end vertices is large enough for splitting
                      for (unsigned int d = 0; valid && d < 2; ++d)
                      {
                        double l = 0;
                        for (std::size_t a = vend[2*d+1], b = (a + 1) % v1.size(); l <= h && b != vend[2-2*d]; b = ((a = b) + 1) % v1.size())
                          l += (v1[b].r - v1[a].r).length();
                        if (l <= h) valid = false;
                      }
                      if (!valid) continue;
                      
                      // distribute the vertices to the two polygons
                      p[1] = p2 = polygons.size();
                      Polygon pn [2];
                      for (unsigned int b = 0; valid && b < 2; ++b)
                      {
                        for (std::size_t c = vend[3-2*b]; ; c = (c + 1) % v1.size())
                        {
                          pn[b].vertices.push_back(v1[c]);
                          pn[b].vertices.back().p = p[b];
                          if (c == vend[2*b]) break;
                        }
                        if (pn[b].vertices.size() < std::max(5., (M_PI - 1) * h / lmax + 1))
                          valid = false; // new polygon must have enough vertices
                      }
                      
                      // if one polygon ends up inside the other, reverse its vertex order
                      for (unsigned int b = 0; valid && b < 2; ++b)
                      {
                        if ((inside[b] = pn[1-b].contains(pn[b].vertices[0].r))) // assignment, not an equality test
                        {
                          std::reverse(pn[b].vertices.begin(), pn[b].vertices.end());
                          for (std::size_t c = 0; c < pn[b].vertices.size() - 1; ++c)
                            pn[b].vertices[c].l0 = pn[b].vertices[c+1].l0;
                        }
                        pn[b].phase = polygons[p1].phase != inside[b]; // flip the phase if needed
                        pn[b].vertices.back().l0 = (pn[b].vertices[0].r - pn[b].vertices.back().r).length() / ls1; // connect the end vertices
                        pn[b].area(); // compute the new area
                        if (pn[b].A < std::max(Amin, (h + lmax) * (h + lmax)))
                          valid = false; // new polygon must have large enough area
                      }
                      if (!valid) continue;
                      
                      # pragma omp critical
                      if (!fuse)
                      {
                        fuse = true;
                        pnew[0] = pn[0];
                        pnew[1] = pn[1];
                        
                        // adjust polygon properties
                        pnew[0].Amax = pnew[1].Amax = polygons[p1].Amax;
                        if (inside[0] || inside[1])
                        {
                          pnew[inside[0]].A0 = pn[inside[1]].A + polygons[p1].A0;
                          pnew[inside[0]].alpha = polygons[p1].alpha;
                          pnew[inside[1]].A0 = pn[inside[1]].A;
                          pnew[inside[1]].alpha = 0;
                        }
                        else
                        {
                          const double f = pn[0].A / (pn[0].A + pn[1].A);
                          pnew[0].A0 = polygons[p1].A0 * f;
                          pnew[0].alpha = polygons[p1].alpha * f;
                          pnew[1].A0 = polygons[p1].A0 * (1 - f);
                          pnew[1].alpha = polygons[p1].alpha * (1 - f);
                        }
                      }
                    }
                  }
                }
        }
      }
      if (!fuse) break; // break out of the while(true) loop if no fusion was prepared
      
      // add/remove the affected polygon(s), then repeat from the beginning, recomputing the boxes
      polygons[pnew[0].vertices[0].p] = pnew[0];
      if (pnew[1].vertices.size() == 0) // fusion type I: remove the second old polygon
        remove(prem);
      else // fusion type II: add the second new polygon
        polygons.push_back(pnew[1]);
    }
    
    // polygon removal and division
    for (std::size_t p = Nr; p < polygons.size(); ++p)
    {
      auto& v = polygons[p].vertices;
      
      if (polygons[p].A < Amin || v.size() < 3)
        remove(p--);
      else if (polygons[p].A > polygons[p].Amax && v.size() > 5)
      {
        // compute polygon centroid, inertia tensor and division axis
        Point c{0, 0};
        double Ixx = 0, Iyy = 0, Ixy = 0;
        for (std::size_t i = v.size() - 1, j = 0; j < v.size(); i = j++)
        {
          const Point rsum = v[i].r + v[j].r;
          const double w = v[i].r ^ v[j].r;
          c.add(w, rsum);
          Ixx += w * (rsum.y * rsum.y - v[i].r.y * v[j].r.y);
          Iyy += w * (rsum.x * rsum.x - v[i].r.x * v[j].r.x);
          Ixy += w * (rsum.x * rsum.y + v[i].r.x * v[i].r.y + v[j].r.x * v[j].r.y);
        }
        c = 1 / (6 * polygons[p].A) * c;
        Ixx =  Ixx / 12 - c.y * c.y * polygons[p].A;
        Iyy =  Iyy / 12 - c.x * c.x * polygons[p].A;
        Ixy = -Ixy / 24 + c.x * c.y * polygons[p].A;
        const double dI = (Ixx - Iyy) / 2;
        const double eig = dI + Iyy + std::sqrt(dI * dI + Ixy * Ixy);
        Point axis = (Ixx < Iyy ? Point{Ixy, eig - Ixx} : Point{eig - Iyy, Ixy});
        
        // use a random cell division axis if the inertia tensor is proportional to the identity matrix
        const double l = axis.length();
        if (l == 0)
        {
          const double angle = 2 * M_PI * uni_dist(rng);
          axis = {std::cos(angle), std::sin(angle)};
        }
        else
          axis = 1 / l * axis; // otherwise, normalize the axis
        
        // divide the polygon in half, ensuring the end vertices are separated enough
        std::size_t vend [4];
        for (unsigned int a = 0, i = v.size() - 1, j = 0; a < 4; i = j++)
        {
          const double d0 = axis ^ (v[i].r - c);
          const double d1 = axis ^ (v[j].r - c);
          if (d0 * d1 <= 0)
          {
            vend[a] = i;
            vend[a+1] = j;
            for (int d = 0; (v[vend[a+1]].r - v[vend[a]].r).length2() < h * h; d = 1 - d)
              vend[a+d] = (vend[a+d] + v.size() + 2*d-1) % v.size();
            j = (vend[a+1] + 1) % v.size();
            a += 2;
          }
        }
        
        // make sure the new edges are stretched just as much as the rest of the polygon
        const double ls = polygons[p].perimeter() / polygons[p].perimeter0(); // average stretch ratio
        v[vend[0]].l0 = (v[vend[3]].r - v[vend[0]].r).length() / ls;
        v[vend[2]].l0 = (v[vend[1]].r - v[vend[2]].r).length() / ls;
        
        // distribute the vertices to the two new polygons
        const double A0 = polygons[p].A0;
        std::vector<Vertex> vold;
        vold.swap(v);
        polygons[p] = {{vold[vend[2]]}, polygons[p].phase, 0, 0, Amax_dist(rng), alpha_dist(rng), D_dist(rng), k_dist(rng)}; // new polygon 1
        polygons.push_back({{vold[vend[0]]}, polygons[p].phase, 0, 0, Amax_dist(rng), alpha_dist(rng), D_dist(rng), k_dist(rng)}); // new polygon 2
        for (std::size_t i = vend[1]; i != vend[2]; i = (i + 1) % vold.size())
          polygons[p].vertices.push_back(vold[i]);
        for (std::size_t i = vend[3]; i != vend[0]; i = (i + 1) % vold.size())
          polygons.back().vertices.push_back(vold[i]);
        for (auto& v2 : polygons.back().vertices)
          v2.p = polygons.size() - 1;
        
        // update polygon areas, splitting the target area in proportion to the actual area
        polygons[p].area();
        polygons.back().area();
        polygons[p].A0 = A0 * polygons[p].A / (polygons[p].A + polygons.back().A);
        polygons.back().A0 = A0 - polygons[p].A0;
        // Polymorph extension: update polygon production rate (note: has to happen after vertices are assigned)
        polygons[p].p = is_producing(polygons[p]) ? p_dist(rng) : 0;
        polygons.back().p = is_producing(polygons.back()) ? p_dist(rng) : 0;
      }
    }
    
    #pragma omp parallel for
    for (std::size_t p = 0; p < polygons.size(); ++p)
    {
      auto& v = polygons[p].vertices;
      
      // refine or coarsen the polygon
      for (std::size_t i = v.size() - 1, k = i - 1, j = 0; j < v.size(); )
      {
        const double l2 = (v[j].r - v[i].r).length2();
        if (l2 > lmax * lmax) // refine long edges, including rigid polygons
        {
          // insert a new vertex in the middle
          const Point rnew = 0.5 * (v[i].r + v[j].r);
          const Point vnew = 0.5 * (v[i].v + v[j].v);
          v[i].l0 /= 2;
          v.insert(v.begin() + j, {rnew, vnew, {0, 0}, p, 0, v[i].l0});
          if (k > j) ++k;
          if (i > j) ++i;
        }
        else if (l2 < lmin * lmin && p >= Nr && v.size() > 3) // coarsen only non-rigid polygons with more than 3 vertices
        {
          // merge the two vertices
          v[k].l0 += v[i].l0 / 2;
          v[i].l0 = v[i].l0 / 2 + v[j].l0;
          v[i].r = 0.5 * (v[i].r + v[j].r);
          v[i].v = 0.5 * (v[i].v + v[j].v);
          v.erase(v.begin() + j);
          if (k > j) --k;
          if (i > j) --i;
        }
        else // continue checking the next edge
          k = i, i = j++;
      }
      
      // compute vertex accelerations
      polygons[p].area(); // update the polygon area
      const double ls = polygons[p].perimeter0() / std::sqrt(Q * 4 * M_PI * polygons[p].A0); // inverse stretch ratio
      for (std::size_t i = v.size() - 1, k = i - 1, j = 0; j < v.size(); k = i, i = j++)
      {
        const Point e1 = v[i].r - v[k].r;
        const Point e2 = v[j].r - v[i].r;
        const double l1 = e1.length();
        const double l2 = e2.length();
        const Point n = (e1 + e2).cross(); // unnormalized inward normal vector
        v[i].a = {0, -gl}; // edge gravitational acceleration
        v[i].a.add(ka / 2 * (polygons[p].A - polygons[p].A0), n); // area compressibility
        v[i].a.add(-kl, (ls / v[k].l0 - 1 / l1) * e1 - (ls / v[i].l0 - 1 / l2) * e2); // edge contractility
        v[i].a.add(-gam, (1 / l1) * e1 - (1 / l2) * e2); // line tension
        v[i].a.add(rho * g / 6 * (2*polygons[p].phase-1), (v[k].r.y + v[i].r.y + v[j].r.y) * n - Point{0, e1 ^ e2}); // hydrostatic pressure
        v[i].a.add(-cv - rho * cd / 4 * std::abs(v[i].v * n), v[i].v); // viscous damping and drag

        // bending
        const Point e0 = v[k].r - v[(k + v.size() - 1) % v.size()].r;
        const Point e3 = v[(j + 1) % v.size()].r - v[j].r;
        const double l0 = e0.length();
        const double l3 = e3.length();
        const double b0 = l0 * l1 + e0 * e1;
        const double b1 = l1 * l2 + e1 * e2;
        const double b2 = l2 * l3 + e2 * e3;
        const double a0 = (e0 ^ e1) / b0;
        const double a1 = (e1 ^ e2) / b1;
        const double a2 = (e2 ^ e3) / b2;
        v[i].a.add(-8 * kb, (a0 / ((l0 + l1) * b0)) * (e0.cross() - a0 * e0)
                          - (a1 / ((l1 + l2) * b1)) * (n - a1 * (e1 - e2))
                          + (a2 / ((l2 + l3) * b2)) * (e3.cross() + a2 * e3));
      }
    }
    
    // polygon-polygon interaction
    boxes(); // place all vertices into boxes
    #pragma omp parallel for
    for (std::size_t p = Nr; p < polygons.size(); ++p)
      for (std::size_t i = polygons[p].vertices.size() - 1, k = i - 1, j = 0; j < polygons[p].vertices.size(); k = i, i = j++)
      {
        const std::size_t bxi = (polygons[p].vertices[i].r.x - x0) / bs + 1; // box index in x direction // NM: is the +1 because of extra boxes?
        const std::size_t byi = (polygons[p].vertices[i].r.y - y0) / bs + 1; // box index in y direction
        for (std::size_t bxj = bxi - 1; bxj <= bxi + 1; ++bxj)  //NM: interaction within 1 box in each direction
          for (std::size_t byj = byi - 1; byj <= byi + 1; ++byj)
            for (Vertex* v = first[bxj * Ny + byj]; v; v = v->next) //NM: in first we get the starting vertex from each box
              if (v != &polygons[p].vertices[k] && v != &polygons[p].vertices[i] && v != &polygons[p].vertices[j])
                interaction(v, &polygons[p].vertices[i], &polygons[p].vertices[k], &polygons[p].vertices[j]);
      }
    
    // time integration (semi-implicit Euler method)
    #pragma omp parallel for
    for (std::size_t p = Nr; p < polygons.size(); ++p)
    {
      for (auto& v : polygons[p].vertices)
      {
        v.v.add(dt, v.a); // update vertex velocity
        v.r.add(dt, v.v); // update vertex position
      }
      polygons[p].A0 += polygons[p].alpha * dt; // apply the area growth rate
      polygons[p].area(); // compute the new polygon area
    }
    t += dt; // advance the time
  } // NM: end step()
  
  void boxes()
  {
    // compute the global bounding box and maximum squared edge length
    double xmin = polygons[0].vertices[0].r.x, xmax = xmin;
    double ymin = polygons[0].vertices[0].r.y, ymax = ymin;
    double l2max = 0;
    #pragma omp parallel for reduction(min:xmin,ymin) reduction(max:xmax,ymax,l2max)
    for (std::size_t p = 0; p < polygons.size(); ++p)
    {
      for (std::size_t i = polygons[p].vertices.size() - 1, j = 0; j < polygons[p].vertices.size(); i = j++)
      {
        if      (polygons[p].vertices[i].r.x < xmin) xmin = polygons[p].vertices[i].r.x;
        else if (polygons[p].vertices[i].r.x > xmax) xmax = polygons[p].vertices[i].r.x;
        if      (polygons[p].vertices[i].r.y < ymin) ymin = polygons[p].vertices[i].r.y;
        else if (polygons[p].vertices[i].r.y > ymax) ymax = polygons[p].vertices[i].r.y;
        l2max = std::max(l2max, (polygons[p].vertices[j].r - polygons[p].vertices[i].r).length2());
      }
    }
    x0 = xmin;
    y0 = ymin;
    x1 = xmax;
    y1 = ymax;
    
    // place vertices in boxes
    bs = std::sqrt(l2max) + drmax; // box size
    Nx = (xmax - x0) / bs + 3; // number of boxes in x direction (with an extra column on both ends)
    Ny = (ymax - y0) / bs + 3; // number of boxes in y direction (with an extra row on both ends)
    first.assign(Nx * Ny, 0); // clear the boxes
    for (auto& p : polygons)
    {
      for (auto& v : p.vertices)
      {
        const std::size_t bx = (v.r.x - x0) / bs + 1; // box index in x direction
        const std::size_t by = (v.r.y - y0) / bs + 1; // box index in y direction
        const std::size_t b = bx * Ny + by; // global box index
        v.next = first[b];
        first[b] = &v; // place this vertex in front of the list in this box
      }
    }
  }
  
  void interaction(Vertex* v0, Vertex* v1, Vertex* v1n0, Vertex* v1n1)
  {
    // select the closer of the two edges
    Vertex* v1n [2] = {v1n0, v1n1}; // pointers to neighbors of vertex 1
    double xi [2], xit [2], dr2 [2];
    Point dr [2];
    for (unsigned int j = 0; j < 2; ++j)
      dr2[j] = point_edge_dist2(v0->r, v1->r, v1n[j]->r, xi[j], xit[j], dr[j]);
    const unsigned int j = dr2[0] > dr2[1];
    
    if (xi[j] <= 1 && dr2[j] < drmax * drmax)
    {
      Vertex* v2 = v1n[j]; // pointer to vertex 2
      
      // check if arclength distance is large enough to allow for self-contact
      if (lmin < h && v0->p == v1->p)
      {
        auto& v = polygons[v1->p].vertices;
        const std::size_t j1 = v1 - v.data(), j2 = v2 - v.data();
        for (unsigned int d = 0; d < 2; ++d)
        {
          double l = 0;
          for (std::size_t i = v0 - v.data(), k; l <= h && i != j1 && i != j2; i = k)
          {
            k = (i + v.size() + 2*d-1) % v.size();
            l += (v[i].r - v[k].r).length();
          }
          if (l <= h) return;
        }
      }
      
      // trilinear traction-separation law including viscous damping
      const double dr_abs = std::sqrt(dr2[j]);
      const Point dv = v0->v - v1->v - xit[j] * (v2->v - v1->v);
      const double dndv_dr = dr[j] * dv / dr2[j];
      double da_dr = -cc * dndv_dr; // viscous dashpot
      if (dr_abs < h)
        da_dr += kr * (h / dr_abs - 1); // linear repulsion
      else if (v0->p != v1->p)
      {
        if (dr_abs < h + sh)
          da_dr += kh * (h / dr_abs - 1); // linear adhesion
        else
          da_dr += kh * sh / ss * (1 - drmax / dr_abs); // softening adhesion
      }
      Point a = da_dr * dr[j]; // acceleration vector
      
      // add dynamic Coulomb friction
      const Point dvt = dv - dndv_dr * dr[j];
      const double dvt_abs = dvt.length();
      if (dvt_abs > 0 && v0->p != v1->p)
        a.add(-mu * std::abs(da_dr) * dr_abs / dvt_abs, dvt);
      
      // distribute the acceleration to the involved vertices
      v0->a.add(1, a);
      v1->a.add(xit[j]-1, a);
      v2->a.add(-xit[j], a);
    }
  }
  
  void output(const std::size_t f)
  {
    // print the frame number, simulated time, number of polygons, and total number of vertices
    std::size_t Nv = 0;
    for (auto& p : polygons)
      Nv += p.vertices.size();
    std::cout << "frame " << f << ", t=" << t << ", " << polygons.size() << " polygons, " << Nv << " vertices\n";
    
    // write a VTK polygon file containing the areas, perimeters, and coordination numbers
    char name [16];
    snprintf(name, 16, "frame%06zu.vtp", f);
    std::ofstream file(name);
    file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <PolyData>\n";
    file << "    <Piece NumberOfPoints=\"" << Nv << "\" NumberOfPolys=\"" << polygons.size() << "\">\n";
    file << "      <CellData>\n";
    file << "        <DataArray type=\"Float64\" Name=\"area\" format=\"ascii\">\n";
    for (auto& p : polygons)
      file << p.A << " ";
    file << "\n";
    file << "        </DataArray>\n";
    // polymorph extension
    // u
    file << "        <DataArray type=\"Float64\" Name=\"u\" format=\"ascii\">\n";
    for (auto& p : polygons)
      file << p.u << " ";
    file << "\n";
    file << "        </DataArray>\n";
    // D
    file << "        <DataArray type=\"Float64\" Name=\"D\" format=\"ascii\">\n";
    for (auto& p : polygons)
      file << p.D << " ";
    file << "\n";
    file << "        </DataArray>\n";
    // k
    file << "        <DataArray type=\"Float64\" Name=\"k\" format=\"ascii\">\n";
    for (auto& p : polygons)
      file << p.k << " ";
    file << "\n";
    file << "        </DataArray>\n";
    // p
    file << "        <DataArray type=\"Float64\" Name=\"p\" format=\"ascii\">\n";
    for (auto& p : polygons)
      file << p.p << " ";
    file << "\n";
    file << "        </DataArray>\n";
    // end polymorph extension
    file << "        <DataArray type=\"Float64\" Name=\"perimeter\" format=\"ascii\">\n";
    for (auto& p : polygons)
      file << p.perimeter() << " ";
    file << "\n";
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Float64\" Name=\"neighbors\" format=\"ascii\">\n";
    boxes(); // place all vertices into boxes
    #pragma omp parallel for ordered schedule(static,1)
    for (std::size_t p = 0; p < polygons.size(); ++p)
    {
      double xi, xit;
      Point dr;
      std::set<std::size_t> n; // set of neighboring polygons
      for (std::size_t i = polygons[p].vertices.size() - 1, k = i - 1, j = 0; j < polygons[p].vertices.size(); k = i, i = j++)
      {
        const std::size_t bxi = (polygons[p].vertices[i].r.x - x0) / bs + 1; // box index in x direction
        const std::size_t byi = (polygons[p].vertices[i].r.y - y0) / bs + 1; // box index in y direction
        for (std::size_t bxj = bxi - 1; bxj <= bxi + 1; ++bxj)
          for (std::size_t byj = byi - 1; byj <= byi + 1; ++byj)
            for (Vertex* v = first[bxj * Ny + byj]; v; v = v->next)
              if (p != v->p) // do not count contact of a polygon with itself
                if (point_edge_dist2(v->r, polygons[p].vertices[i].r, polygons[p].vertices[k].r, xi, xit, dr) < drmax * drmax ||
                    point_edge_dist2(v->r, polygons[p].vertices[i].r, polygons[p].vertices[j].r, xi, xit, dr) < drmax * drmax)
                  n.insert(v->p);
      }
      #pragma omp ordered
      file << n.size() << " ";
    }
    file << "\n";
    file << "        </DataArray>\n";
    file << "      </CellData>\n";
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (auto& p : polygons)
      for (auto& v : p.vertices)
        file << v.r.x << " " << v.r.y << " " << p.phase << " ";
    file << "\n";
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    file << "      <Polys>\n";
    file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < Nv; ++i)
      file << i << " ";
    file << "\n";
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    std::size_t offset = 0;
    for (auto& p : polygons)
      file << (offset += p.vertices.size()) << " ";
    file << "\n";
    file << "        </DataArray>\n";
    file << "      </Polys>\n";
    file << "    </Piece>\n";
    file << "  </PolyData>\n";
    file << "</VTKFile>\n";
  }
};

// takes care of the data scattering and gathering between ensemble and solver
struct Interpolator {
  Ensemble& ensemble;
  Solver& solver;
  Interpolator(Ensemble& ensemble, Solver& solver) : ensemble(ensemble), solver(solver) {}
  
  // Search algorithm to find parent polygon for a grid point.
  // Complexity almost O(1). More precicely O(vertices per polygon) which is almost constant
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
      if (last_iteration) return -1; // background node
      if (!checked_polygons.empty()) last_iteration = true; // only go 1 more layer
      ++bxi;
    }
    return -1; // background node (reached boundary of ensemble box)
  }

  // scatter coefficients D, k from polygons to grid points
  void scatter() { // ToDo: make this function prettier
    Grid<int>& prev_idx = solver.parent_idx; // stores the polygon index of the cell in which a grid point lies (its parent)
    Grid<int> new_idx(solver.Nx, solver.Ny, -1); // negative indices indicate a background node. ToDo: could make this in place
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < solver.Nx; i++) {
      for (int j = 0; j < solver.Ny; j++) { 
        // spatial coordinates of grid point
        const double x = solver.box_position_x + i * solver.dx;
        const double y = solver.box_position_y + j * solver.dx;
        const Point grid_point(x, y);
        // check if outside the ensemble box ToDo: Limit loop to those boundaries
        if (x < ensemble.x0 || x > ensemble.x1 || y < ensemble.y0 || y > ensemble.y1) {
          new_idx(i, j) = -2; // external background node
        } 
        // check if still the same parent
        else if (prev_idx(i, j) >= 0 && ensemble.polygons[prev_idx(i, j)].contains(grid_point)) { 
          new_idx(i, j) = prev_idx(i, j);
        } 
        // onion layer search algorithm
        else {
          new_idx(i, j) = find_parent(grid_point);
        }
        // scatter values
        if (new_idx(i, j) < 0) { // is background node
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
    solver.parent_idx = new_idx;
  }

  // gather concentration u from grid points to polygons
  // important: depends on scatter being called every iteration to build the parent_idx
  void gather() {
    // get all children from parent idx built during scatter()
    for (auto& cell : ensemble.polygons) {
      cell.children.clear();
    }
    for (int i = 0; i < solver.Nx; i++) {
      for (int j = 0; j < solver.Ny; j++) {
        if (solver.parent_idx(i, j) >= 0) { // skip background nodes
          auto& cell = ensemble.polygons[solver.parent_idx(i, j)]; 
          cell.children.push_back(Index(i, j)); 
        }
      }
    }
    // accumulate data from children
    for (auto& cell : ensemble.polygons) {
      // average concentration
      cell.u = 0;
      for (const Index& idx : cell.children) {
        cell.u += solver.u(idx);
      }
      if (cell.children.size() > 0) { // avoid division by zero if cells exceed RD box
        cell.u /= cell.children.size();
      }
    }
  }
};


int main()
{
  // ToDo: fix seed of rng

  Ensemble ensemble("ensemble.off"); // read the input file
  ensemble.output(0); // print the initial state
  
  Grid<double> u0(100, 100); // initial condition, just zeros
  Solver solver(u0, D_mu, dx, dt, LinearDegradation(k_mu));
  solver.output(0); // print the initial state
  
  Interpolator interpolator(ensemble, solver);

  for (std::size_t f = 1; f <= Nf; ++f)
  {
    for (std::size_t s = 0; s < Ns; ++s) 
    {
      ensemble.step();
      interpolator.scatter(); 
      solver.step();
      interpolator.gather();
    }
    // for testing purposes only every frame. no interaction yet. 
    //interpolator.scatter(); 
    //interpolator.gather();
    
    ensemble.output(f); // print a frame
    solver.output(f); // print a frame
  }
}
