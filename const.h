#ifndef CONST_H
#define CONST_H

#include <cmath>

constexpr double h = 0.01; // [L] edge thickness
constexpr double lmin = 0.02; // [L] minimum edge length
constexpr double lmax = 0.2; // [L] maximum edge length
constexpr double Q = 1; // [-] isoparametric ratio
constexpr double alpha_mu = 1; // [L^2/T] mean area growth rate
constexpr double alpha_CV = 0; // [-] coefficient of variation of area growth rate
constexpr double beta = 0.9; // [-] minimum area fraction for growth
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
constexpr double dt = 1e-4; // [T] time step // default 1e-4

constexpr std::size_t Nf = 100; // number of output frames
constexpr std::size_t Ns = 3000; // number of time steps between frames // default 1000
constexpr std::size_t Nr = 1; // number of rigid polygons

constexpr double drmax = h + sh + ss; // maximum interaction distance

// PolyMorph extension
constexpr double dx = 0.4; // [L] grid spacing for solver
constexpr double D_mu = 128.0; // [L^2/T] diffusion coefficient mean
constexpr double k_mu = 1.0; // [1/T] degradation rate mean 
constexpr double p_mu = 24.0; // [1/T] production rate mean
constexpr double threshold_mu = 0.1; // [-] threshold for morphogen concentration mean
constexpr double D_CV = 0.3; // [-] coefficient of variation of diffusion
constexpr double k_CV = 0.3; // [-] coefficient of variation of degradation rate
constexpr double p_CV = 0.3; // [-] coefficient of variation of production rate
constexpr double threshold_CV = 0.3; // [-] coefficient of variation of threshold
constexpr double D0 = 64.0; // [L^2/T] diffusion coefficient background
constexpr double k0 = 0.0; // [1/T] reaction rate background
constexpr double p0 = 0.0; // [1/T] reaction rate background

#endif