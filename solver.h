#include <vector>
#include <functional>
#include <fstream>
#include <cmath>

// size of the FD grid. 
// x always corresponds to i and y to j.
// ToDo: make members of grid
constexpr int Nx = 100; 
constexpr int Ny = 100;

struct Index {
  const int i; 
  const int j;
  Index(int i, int j): i(i), j(j) {}
};

template<typename T>
struct Grid {
    std::vector<std::vector<T>> data;
    Grid () : data(Nx, std::vector<T>(Ny, T(0))) {} // constructor always creates Nx x Ny grid
    Grid (double value) : data(Nx, std::vector<T>(Ny, value)) {}
    T& operator()(int i, int j) { return data[i][j]; }
    T& operator()(Index idx) { return data[idx.i][idx.j]; }
};

struct Reaction {
    double operator()(double u, int i, int j) {
        return 0.0;
    }
};

struct LinearDegradation : Reaction {
    double k;
    LinearDegradation(double k): k(k) {}
    double operator()(double u, int i, int j) {
        return -k * u;
    }
};

struct Solver { 
    double box_position_x;  // ToDo: change this ugly datastruct. Maybe a dimensions struct
    double box_position_y;   
    double D; // diffusion coefficient. Later also a grid datastructure
    double dx; // grid spacing
    double dt; // time step
    Reaction R; // reaction term
    Grid<double> u; // concentration of the morphogen
    Grid<int> parent_idx; // polygon idx

    // initialize the grid with a given initial condition
    Solver(const Grid<double> u0, const double D = 1.0, const double dx = 0.1, 
            double dt = 1e-4, Reaction R = LinearDegradation(0.1)) {
        this->u = u0;
        this->D = D;
        this->dx = dx;
        this->dt = dt;
        this->R = R;
        box_position_x = - Nx/2 * dx; // midpoint at 0
        box_position_y = - Ny/2 * dx;
        parent_idx = Grid<int>(-1); // initialize as all background nodes. Maybe change to std::map?
    }

    void step() { 
        Grid<double> unew;
        // Forward Euler with central differences ToDo: adapt for variable diffusion coefficient
        // maybe separate inner nodes and boundary to efficiently parallelize and vectorize inner nodes
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {   
                // mirror past-boundary nodes
                const double n = (j == Ny-1) ? u(i, j-1) : u(i, j+1);
                const double s = (j == 0)    ? 1         : u(i, j-1); // dirichlet boundary condition
                const double e = (i == Nx-1) ? u(i-1, j) : u(i+1, j);
                const double w = (i == 0)    ? u(i+1, j) : u(i-1, j);
                unew(i, j) = u(i, j) + dt * (
                    D / (dx*dx) * (n + s + e + w - 4 * u(i, j))
                    + R(u(i, j), i, j)
                ); 
            }
        }      
        u = unew; 
    }


    void output(const std::size_t frame) {    // f: frame number
        char filename [19]; 
        snprintf(filename, 19, "rd_frame%06zu.vts", frame);        
        std::ofstream file(filename);

        file << "<?xml version=\"1.0\"?>" << std::endl;
        file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">" << std::endl;
        file << "<StructuredGrid WholeExtent=\"0 " << Nx-1 << " 0 " << Ny-1 << " 0 0\">" << std::endl;
        file << "<Piece Extent=\"0 " << Nx-1 << " 0 " << Ny-1 << " 0 0\">" << std::endl;
        file << "<Points>" << std::endl;
        file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double x = box_position_x + i * dx;
                double y = box_position_y + j * dx;
                file << x << " " << y << " 0" << std::endl;
            }
        }
        file << "</DataArray>" << std::endl;
        file << "</Points>" << std::endl;
        file << "<PointData Scalars=\"scalars\">" << std::endl;
        // u
        file << "<DataArray type=\"Float64\" Name=\"u\" format=\"ascii\">" << std::endl;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                file << u(i, j) << " ";
            }
            file << std::endl;
        }
        file << "</DataArray>" << std::endl;
        // parent_idx
        file << "<DataArray type=\"Int32\" Name=\"parent_idx\" format=\"ascii\">" << std::endl;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                file << parent_idx(i, j) << " ";
            }
            file << std::endl;
        }
        file << "</DataArray>" << std::endl;
        file << "</PointData>" << std::endl;

        file << "</Piece>" << std::endl;
        file << "</StructuredGrid>" << std::endl;
        file << "</VTKFile>" << std::endl;
        file.close();         
    }
};

// helper function for IC

Grid<double> create_gaussian() {
    Grid<double> u;
    for(int i = 0; i < Nx; i++)
        for(int j = 0; j < Ny; j++)
            u(i, j) = std::exp(-((i - Nx/2)*(i - Nx/2) + (j - Ny/2)*(j - Ny/2)) / 100.0);
    u(Nx/2, Ny/2) = 1.0;
    return u;
}