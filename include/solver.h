#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <functional>
#include <fstream>
#include <iostream>

#include "reaction.h"
#include "grid.h"
#include "geometry.h"
#include "domain.h"

struct BoundaryCondition {
    enum class Type {
        Dirichlet,
        Neumann
    };

    Type type;
    double value; // value at boundary (Dirichlet case) or derivative at boundary (Neumann case)
};

struct Boundary { // ToDo: allow multiple species
    BoundaryCondition north, south, east, west;
    
    static Boundary zeroFlux() {
        return {
            {BoundaryCondition::Type::Neumann, 0},
            {BoundaryCondition::Type::Neumann, 0},
            {BoundaryCondition::Type::Neumann, 0},
            {BoundaryCondition::Type::Neumann, 0}
        };
    }
};

struct Solver { 
    Domain& domain;
    int Nx, Ny; // number of grid points
    double dx; // grid spacing
    Boundary boundary; // boundary conditions
    Reaction R = [](const std::vector<double>& u, const std::vector<double>& k) { return std::vector<double>(NUM_SPECIES, 0.0); }; // reaction term
    Grid<int> parent_idx; // polygon idx
    Grid<std::vector<double>> u; // concentrations
    Grid<std::vector<double>> D; // diffusion coefficients
    Grid<std::vector<double>> p; // production rates
    Grid<std::vector<double>> k; // kinetic coefficients
    Grid<Point> velocity; // velocity field
    Grid<std::vector<Point>> grad_u; // concentration gradient

    Solver(Domain& domain, const double dx, Reaction R) : domain(domain), R(R), dx(dx) {
        this->boundary = Boundary::zeroFlux();
        this->Nx = domain.width() / dx + 1;
        this->Ny = domain.height() / dx + 1;
        this->u = Grid<std::vector<double>>(Nx, Ny, std::vector<double>(NUM_SPECIES, 0.0));

        std::cout << "solver dimensions Nx=" << Nx << " Ny=" << Ny << std::endl;
        std::cout << "domain: [" << domain.x0 << ", " << domain.x1 << "] x [" << domain.y0 << ", " << domain.y1 << "]" << std::endl;
        std::cout << "dx=" << dx << std::endl;

        // initialize with background values
        parent_idx = Grid<int>(Nx, Ny, -2);
        D = Grid<std::vector<double>>(Nx, Ny, D0); // ToDo: Benchmark and maybe remove D,p,k grids (just use values from parent polygon)
        p = Grid<std::vector<double>>(Nx, Ny, p0);
        k = Grid<std::vector<double>>(Nx, Ny, k0);
        grad_u = Grid<std::vector<Point>>(Nx, Ny, std::vector<Point>(NUM_SPECIES, Point(0, 0)));
        if (ADVECTION_DILUTION) {
            velocity = Grid<Point>(Nx, Ny, Point(0, 0));
        }
    }

    // rescale all grids
    void rescale(size_t Nx_new, size_t Ny_new, int offset_x, int offset_y) {
        u.rescale(Nx_new, Ny_new, offset_x, offset_y, std::vector<double>(NUM_SPECIES, 0.0));
        D.rescale(Nx_new, Ny_new, offset_x, offset_y, D0);
        p.rescale(Nx_new, Ny_new, offset_x, offset_y, p0);
        k.rescale(Nx_new, Ny_new, offset_x, offset_y, k0);
        parent_idx.rescale(Nx_new, Ny_new, offset_x, offset_y, -2);
        if (ADVECTION_DILUTION) {
            velocity.rescale(Nx_new, Ny_new, offset_x, offset_y, Point(0, 0));
        }
        this->Nx = Nx_new;
        this->Ny = Ny_new;
        std::cout << "rescaled solver to Nx=" << Nx << " Ny=" << Ny << std::endl;
    }

    void step(double dt) {
        // resize grids if necessary TODO reenable
        if (RESIZE_GRID) {
            if (domain.width() >= (Nx + 1) * dx || domain.height() >= (Ny + 1) * dx
                || domain.width() <= (Nx - 1) * dx || domain.height() <= (Ny - 1) * dx) {
                int Nx_new = floor(domain.width() / dx);
                int Ny_new = floor(domain.height() / dx);
                rescale(Nx_new, Ny_new, 0, 0);
            }
        }

        Grid<std::vector<double>> unew(Nx, Ny, std::vector<double>(NUM_SPECIES, 0.0));

        // Forward Euler with central differences
        #pragma omp parallel for
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {   
                const std::vector<double> reaction = R(u(i, j), k(i, j));
                for (int sp = 0; sp < NUM_SPECIES; sp++) {
                    // dirichlet boundary conditions
                    if      (i == 0     && boundary.west.type  == BoundaryCondition::Type::Dirichlet) { unew(i, j)[sp] = boundary.west.value; continue; } 
                    else if (i == Nx-1  && boundary.east.type  == BoundaryCondition::Type::Dirichlet) { unew(i, j)[sp] = boundary.east.value; continue; } 
                    else if (j == 0     && boundary.south.type == BoundaryCondition::Type::Dirichlet) { unew(i, j)[sp] = boundary.south.value; continue; } 
                    else if (j == Ny-1  && boundary.north.type == BoundaryCondition::Type::Dirichlet) { unew(i, j)[sp] = boundary.north.value; continue; } 
                    else {
                        // account for neumann BDC
                        const double n = (j == Ny-1) ? u(i, j-1)[sp] + 2*dx*boundary.north.value : u(i, j+1)[sp]; 
                        const double s = (j == 0)    ? u(i, j+1)[sp] - 2*dx*boundary.south.value : u(i, j-1)[sp];
                        const double e = (i == Nx-1) ? u(i-1, j)[sp] + 2*dx*boundary.east.value  : u(i+1, j)[sp];
                        const double w = (i == 0)    ? u(i+1, j)[sp] - 2*dx*boundary.west.value  : u(i-1, j)[sp];
                        // calculate diffusion term
                        const double diffusion = D(i, j)[sp] / (dx*dx) * (e + w + anisotropy * (n + s) - 2 * (1 + anisotropy) * u(i, j)[sp]);
                        // update grid point
                        unew(i, j)[sp] = u(i, j)[sp] + dt * (diffusion + reaction[sp] + p(i, j)[sp]); 
                        grad_u(i, j)[sp] = {(e - w) / (2 * dx), (n - s) / (2 * dx)};
                        
                        if(ADVECTION_DILUTION) {
                            const double advection = velocity(i, j) * grad_u(i, j)[sp]; // dot product
                            const double dvdx = (j == Ny-1 || j == 0) ? 0 : (velocity(i, j+1).y - velocity(i, j-1).y) / (2 * dx);
                            const double dvdy = (i == Nx-1 || i == 0) ? 0 : (velocity(i+1, j).x - velocity(i-1, j).x) / (2 * dx);
                            const double dilution = u(i, j)[sp] * (dvdx + dvdy);
                            // update grid point
                            unew(i, j)[sp] -= dt * (advection + dilution);
                        }                    
                    }
                }
            }
        }
        // update state
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
                double x = domain.x0 + i * dx;
                double y = domain.y0 + j * dx;
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

        //std::cout << "Solver frame with: (Nx,Ny) =" << Nx << "," << Ny << std::endl
        //    << "Nx*dx=" << Nx*dx << " domain width=" << domain.width() << std::endl;  

        //if (frame == Nf) output_pvd(); // output pvd file on last frame. Necessary for varying grid sizes          
    }

    void output_pvd() {  // TODO remove or fix visualization  
        std::ofstream file("rd_frame.pvd");
        file << "<?xml version=\"1.0\"?>" << std::endl;
        file << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl;
        file << "<Collection>" << std::endl;
        for (std::size_t f = 0; f <= Nf; ++f) {
            char filename [19]; 
            snprintf(filename, 19, "rd_frame%06zu.vts", f);
            file << "<DataSet timestep=\"" << f << "\" file=\"" << filename << "\"/>" << std::endl;
        }
        file << "</Collection>" << std::endl;
        file << "</VTKFile>" << std::endl;
        file << std::endl;
        file.close();
    }

    // generate noisy initial condition
    Grid<std::vector<double>> noisy_ic(double mean, double stddev) {
        std::default_random_engine generator;
        generator.seed(RNG_SEED);
        std::normal_distribution dist(mean, stddev);
        Grid<std::vector<double>> u0(Nx, Ny, std::vector<double>(NUM_SPECIES, 0.0));
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                for (int sp = 0; sp < NUM_SPECIES; sp++) {
                    u0(i, j)[sp] = dist(generator);
                }
            }
        }
        return u0;
    }
};

#endif