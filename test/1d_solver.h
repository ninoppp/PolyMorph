#include <vector>

struct Reaction {
    double operator()(double u, int i, int j) {
        return 0.0; // ToDo: implement reaction function
    }
};

struct LinearDegradation : Reaction {
    double k;
    LinearDegradation(double k): k(k) {}
    double operator()(double u, int i, int j) {
        return -k * u;
    }
};

struct Solver1D {  
    std::vector<double> u; // concentration of the morphogen
    double D; // diffusion coefficient. Later also a grid datastructure
    double dx; // grid spacing
    double dt; // time step
    Reaction R; // reaction term

    // initialize the grid with a given initial condition
    Solver1D(std::vector<double> u0, const double D = 1.0, const double dx = 0.1, 
            double dt = 1e-4, Reaction R = LinearDegradation(0.1)) {
        this->u = u0;
        this->D = D;
        this->dx = dx;
        this->dt = dt;
        this->R = R;
    }

    void step() { 
        std::vector<double> unew (u.size(), 0.0);
        // Forward Euler with central differences ToDo: adapt for variable diffusion coefficient
        // maybe separate inner nodes and boundary to efficiently parallelize and vectorize inner nodes
        #pragma omp parallel for
        for (int i = 1; i < u.size() - 1; i++) {
            unew[i] = u[i] + dt * (
                D / (dx*dx) * (u[i+1] + u[i-1] - 2 * u[i])
                + R(u[i], 0, 0)
            ); 
        }
        unew[0] = 1;
        unew[u.size() - 1] = unew[u.size() - 2]; // consider last cell as ghost cell
        u = unew; 
    }
};