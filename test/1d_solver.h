#include <vector>
#include <cmath>

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
    double D0; // diffusion coefficient. Later also a grid datastructure
    double j; // influx
    double dx; // grid spacing
    double dt; // time step
    double k;
    double norm = 0; // keep track of change in u
    std::string mode = "dirichlet"; // boundary condition mode

    // initialize the grid with a given initial condition
    Solver1D(std::vector<double> u0, const double D0, const double dx, 
            double dt, double k, std::string mode = "dirichlet", const double j = 0) {
        this->u = u0;
        this->D0 = D0;
        this->dx = dx;
        this->dt = dt;
        this->k = k;
        this->mode = mode;
        this->j = j;
    }
    //tau(x) = (1+x/lambda ) / (2k) time to stead state

    void step() { 
        std::vector<double> unew (u.size(), 0.0);
        // Forward Euler with central differences ToDo: adapt for variable diffusion coefficient
        // maybe separate inner nodes and boundary to efficiently parallelize and vectorize inner nodes
        #pragma omp parallel for
        for (int i = 1; i < u.size() - 1; i++) {
            unew[i] = u[i] + dt * (
                D0 / (dx*dx) * (u[i+1] + u[i-1] - 2 * u[i])
                - k*u[i]
            ); 
        }
        // BDC
        unew[0] = 1;
        unew[u.size() - 1] = 0;

        if (mode == "1xNeumann" || mode == "2xNeumann") { // zero flux right side
            int i = u.size() - 1;
            unew[i] = u[i] + dt * (
                D0 / (dx*dx) * (2*u[i-1] - 2*u[i])
                - k*u[i]
            );
        }
        if (mode == "2xNeumann") { // influx left side
            double u_left = u[1] + 2*dx*j;
            unew[0] = u[0] + dt * (
                D0 / (dx*dx) * (u[1] - 2*u[0] + u_left)
                - k*u[0]
            );
        }

        norm = get_norm(u, unew);
        u = unew; 
    }

    double get_norm(std::vector<double> u1, std::vector<double> u2) {
        double sum = 0.0;
        for (int i = 0; i < u1.size(); i++) {
            sum += (u1[i] - u2[i]) * (u1[i] - u2[i]);
        }
        return std::sqrt(sum);
    }
};