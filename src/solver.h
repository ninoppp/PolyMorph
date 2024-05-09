#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <functional>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iostream>
#include <cassert>

#include "reaction.h"
#include "grid.h"
#include "geometry.h"

enum class BoundaryCondition {
    Dirichlet,
    Neumann,
    Mixed, // 1 at west boundary, 0 at east boundary, zero-flux at north and south
};

struct Solver { 
    double x0, y0; // offset/position of the grid
    size_t Nx, Ny; // number of grid points
    double dx; // grid spacing
    Reaction R; // reaction term
    Grid<int> parent_idx; // polygon idx
    Grid<std::vector<double>> u; // concentrations
    Grid<std::vector<double>> D; // diffusion coefficients
    Grid<std::vector<double>> p; // production rates
    Grid<std::vector<double>> k; // kinetic coefficients
    Grid<Point> velocity; // velocity field
 
    // initialize the grid with a given initial condition
    Solver(const Grid<std::vector<double>> u0, const double dx, Reaction R) {
        this->u = u0;
        this->dx = dx;
        this->R = R;
        this->Nx = u.sizeX();
        this->Ny = u.sizeY();
        std::cout << "solver dimensions Nx=" << Nx << " Ny=" << Ny << std::endl;
        // ToDo: make x0, y0 optional constructor arguments
        x0 = -0.5 * Nx * dx; // midpoint at 0
        y0 = -0.5 * Ny * dx;
        std::cout << "dx=" << dx << std::endl;
        std::cout << "x0=" << x0 << " y0=" << y0 << std::endl;
        // initialize with background values
        parent_idx = Grid<int>(Nx, Ny, -2);
        D = Grid<std::vector<double>>(Nx, Ny, D0);
        p = Grid<std::vector<double>>(Nx, Ny, p0);
        k = Grid<std::vector<double>>(Nx, Ny, k0); 
        if (ADVECTION_DILUTION) {
            velocity = Grid<Point>(Nx, Ny, Point(0, 0));
        }
    }

    void step(double dt) {
        Grid<std::vector<double>> unew(Nx, Ny, std::vector<double>(NUM_SPECIES, 0.0));
        // Forward Euler with central differences
        // inner nodes
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < Nx-1; i++) {
            for (int j = 1; j < Ny-1; j++) {   
                const std::vector<double> reaction = R(u(i, j), k(i, j));
                // mirror past-boundary nodes
                for (int sp = 0; sp < NUM_SPECIES; sp++) { // don't parallelize inner loop. likely to be vectorized
                    // calculate diffusion term
                    const double n = u(i, j+1)[sp];
                    const double s = u(i, j-1)[sp];
                    const double e = u(i+1, j)[sp]; 
                    const double w = u(i-1, j)[sp];
                    const double diffusion = D(i, j)[sp] / (dx*dx) * (e + w + anisotropy * (n + s) - 2 * (1 + anisotropy) * u(i, j)[sp]);
                    // calculate advection term
                    const Point grad_u = Point((e - w) / (2 * dx), (n - s) / (2 * dx));
                    const double advection = velocity(i, j) * grad_u;
                    // calculate dilution term
                    const double dvdx = (velocity(i+1, j).x - velocity(i-1, j).x) / (2 * dx);
                    const double dvdy = (velocity(i, j+1).y - velocity(i, j-1).y) / (2 * dx);
                    const double dilution = u(i, j)[sp] * (dvdx + dvdy);
                    // update grid point
                    unew(i, j)[sp] = u(i, j)[sp] + dt * (diffusion + reaction[sp] + p(i, j)[sp]); 
                    if(ADVECTION_DILUTION) {
                        unew(i, j)[sp] -= dt * (advection + dilution);
                    }
                }
            }
        }
        // north boundary // ToDo: maybe parallelize, maybe not. Maybe merge back into loop above (add Neumann coeff and Dirichlet option)
        for (int i = 1; i < Nx-1; i++) {
            int j = Ny-1;
            const std::vector<double> reaction = R(u(i, j), k(i, j));
            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                const double n = u(i, j-1)[sp]; // mirror past-boundary nodes
                const double s = u(i, j-1)[sp];
                const double e = u(i+1, j)[sp]; 
                const double w = u(i-1, j)[sp];
                const double diffusion = D(i, j)[sp] / (dx*dx) * (e + w + anisotropy * (n + s) - 2 * (1 + anisotropy) * u(i, j)[sp]);
                const Point grad_u = Point((e - w) / (2 * dx), (n - s) / (2 * dx));
                const double advection = velocity(i, j) * grad_u;
                const double dvdx = (velocity(i+1, j).x - velocity(i-1, j).x) / (2 * dx);
                const double dvdy = 0; // zero-gradient velocity at boundary
                const double dilution = u(i, j)[sp] * (dvdx + dvdy);
                unew(i, j)[sp] = u(i, j)[sp] + dt * (diffusion + reaction[sp] + p(i, j)[sp]); 
                if(ADVECTION_DILUTION) {
                    unew(i, j)[sp] -= dt * (advection + dilution);
                } 
            }
        }
        // south boundary
        for (int i = 1; i < Nx-1; i++) {
            int j = 0;
            const std::vector<double> reaction = R(u(i, j), k(i, j));
            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                const double n = u(i, j+1)[sp];
                const double s = u(i, j+1)[sp]; // mirror past-boundary nodes
                const double e = u(i+1, j)[sp]; 
                const double w = u(i-1, j)[sp];
                const double diffusion = D(i, j)[sp] / (dx*dx) * (e + w + anisotropy * (n + s) - 2 * (1 + anisotropy) * u(i, j)[sp]);
                const Point grad_u = Point((e - w) / (2 * dx), (n - s) / (2 * dx));
                const double advection = velocity(i, j) * grad_u;
                const double dvdx = (velocity(i+1, j).x - velocity(i-1, j).x) / (2 * dx);
                const double dvdy = 0; // zero-gradient velocity at boundary
                const double dilution = u(i, j)[sp] * (dvdx + dvdy);
                unew(i, j)[sp] = u(i, j)[sp] + dt * (diffusion + reaction[sp] + p(i, j)[sp]); 
                if(ADVECTION_DILUTION) {
                    unew(i, j)[sp] -= dt * (advection + dilution); 
                }
            }
        }
        // east boundary
        for (int j = 1; j < Ny-1; j++) {
            int i = Nx-1;
            const std::vector<double> reaction = R(u(i, j), k(i, j));
            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                const double n = u(i, j+1)[sp];
                const double s = u(i, j-1)[sp];
                const double e = u(i-1, j)[sp]; // mirror past-boundary nodes
                const double w = u(i-1, j)[sp];
                const double diffusion = D(i, j)[sp] / (dx*dx) * (e + w + anisotropy * (n + s) - 2 * (1 + anisotropy) * u(i, j)[sp]);
                const Point grad_u = Point((e - w) / (2 * dx), (n - s) / (2 * dx));
                const double advection = velocity(i, j) * grad_u;
                const double dvdx = 0; // zero-gradient velocity at boundary
                const double dvdy = (velocity(i, j+1).y - velocity(i, j-1).y) / (2 * dx);
                const double dilution = u(i, j)[sp] * (dvdx + dvdy);
                unew(i, j)[sp] = u(i, j)[sp] + dt * (diffusion + reaction[sp] + p(i, j)[sp]); 
                if(ADVECTION_DILUTION) {
                    unew(i, j)[sp] -= dt * (advection + dilution); 
                }
            }
        }
        // west boundary
        for (int j = 1; j < Ny-1; j++) {
            int i = 0;
            const std::vector<double> reaction = R(u(i, j), k(i, j));
            for (int sp = 0; sp < NUM_SPECIES; sp++) {
                const double n = u(i, j+1)[sp];
                const double s = u(i, j-1)[sp];
                const double e = u(i+1, j)[sp]; 
                const double w = u(i+1, j)[sp]; // mirror past-boundary nodes
                const double diffusion = D(i, j)[sp] / (dx*dx) * (e + w + anisotropy * (n + s) - 2 * (1 + anisotropy) * u(i, j)[sp]);
                const Point grad_u = Point((e - w) / (2 * dx), (n - s) / (2 * dx));
                const double advection = velocity(i, j) * grad_u;
                const double dvdx = 0; // zero-gradient velocity at boundary
                const double dvdy = (velocity(i, j+1).y - velocity(i, j-1).y) / (2 * dx);
                const double dilution = u(i, j)[sp] * (dvdx + dvdy);
                unew(i, j)[sp] = u(i, j)[sp] + dt * (diffusion + reaction[sp] + p(i, j)[sp]); 
                if(ADVECTION_DILUTION) {
                    unew(i, j)[sp] -= dt * (advection + dilution); 
                }
            }
        }
        // ---- update state ----
        u = unew;
    }

    void output(const std::size_t frame) {
        char filename [19]; 
        snprintf(filename, 19, "rd_frame%06zu.vts", frame);        
        std::ofstream file(filename);

        file << "<?xml version=\"1.0\"?>" << std::endl;
        file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">" << std::endl;
        file << "<StructuredGrid WholeExtent=\"0 " << Nx-1 << " 0 " << Ny-1 << " 0 0\">" << std::endl;
        file << "<Piece Extent=\"0 " << Nx-1 << " 0 " << Ny-1 << " 0 0\">" << std::endl;
        // define points
        file << "<Points>" << std::endl;
        file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double x = x0 + i * dx;
                double y = y0 + j * dx;
                file << x << " " << y << " 0" << std::endl;
            }
        }
        file << "</DataArray>" << std::endl;
        file << "</Points>" << std::endl;
        file << "<PointData Scalars=\"scalars\">" << std::endl; // start point data
        
        if (Output::u) file << u.to_vtk("u");
        if (Output::parent_idx) file << parent_idx.to_vtk("parent_idx");
        if (Output::D) file << D.to_vtk("D");
        if (Output::k) file << k.to_vtk("k");
        if (Output::p) file << p.to_vtk("p");
        if (Output::velocity && ADVECTION_DILUTION) file << velocity.to_vtk("velocity");
        
        file << "</PointData>" << std::endl;    // end of point data
        file << "</Piece>" << std::endl;
        file << "</StructuredGrid>" << std::endl;
        file << "</VTKFile>" << std::endl;
        file.close();         
    }
};

#endif