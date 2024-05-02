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
    }

    void step(double dt) { 
        Grid<std::vector<double>> unew(Nx, Ny, std::vector<double>(NUM_SPECIES, 0.0));
        // Forward Euler with central differences
        // ToDo: separate inner nodes and boundary to efficiently parallelize and vectorize inner nodes
        // while allowing different boundary conditions
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {   
                const std::vector<double> r = R(u(i, j), k(i, j));
                assert(r.size() == NUM_SPECIES);
                for (int sp = 0; sp < NUM_SPECIES; sp++) { // don't parallelize inner loop, access to same array 
                    // mirror past-boundary nodes
                    const double n = (j == Ny-1) ? u(i, j-1)[sp] : u(i, j+1)[sp]; // ToDo: option for different BDC
                    const double s = (j == 0)    ? u(i, j+1)[sp] : u(i, j-1)[sp];
                    const double e = (i == Nx-1) ? u(i-1, j)[sp] : u(i+1, j)[sp]; 
                    const double w = (i == 0)    ? u(i+1, j)[sp] : u(i-1, j)[sp]; 
                    unew(i, j)[sp] = u(i, j)[sp] + dt * (
                        D(i, j)[sp] / (dx*dx) * (n + s + e + w - 4 * u(i, j)[sp])
                        + r[sp]
                        + p(i, j)[sp]
                    ); 
                }
            }
        }
        // temporary dirichlet 1 bdc for 2nd morphogen
        /*for (int j = 0; j < Ny; j++) {
            unew(0, j)[1] = 1.0; 
        }*/
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
        
        file << "</PointData>" << std::endl;    // end of point data
        file << "</Piece>" << std::endl;
        file << "</StructuredGrid>" << std::endl;
        file << "</VTKFile>" << std::endl;
        file.close();         
    }
};

#endif