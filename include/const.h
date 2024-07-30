#ifndef CONST_H
#define CONST_H

#include <cmath>
#include <vector>

#define DEBUG false

// -- Original Polyhoop parameters --
constexpr double h = 0.01; // [L] edge thickness
constexpr double lmin = 0.02; // [L] minimum edge length
constexpr double lmax = 0.2; // [L] maximum edge length
constexpr double Q = 1; // [-] isoparametric ratio
constexpr double alpha_mu = 1; // [L^2/T] mean area growth rate
constexpr double alpha_CV = 0.03; // [-] coefficient of variation of area growth rate
constexpr double beta = 0.8; // [-] minimum area fraction for growth
constexpr double Amin = 0; // [L^2] minimum area
constexpr double Amax_mu = M_PI; // [L^2] mean maximum area
constexpr double Amax_CV = 0.03; // [-] coefficient of variation of maximum area // default 0
constexpr double gam = 1e4; // [L/T^2] line tension per vertex mass
constexpr double ka = 1e5; // [1/(L^2*T^2)] area stiffness per vertex mass // default 1e5
constexpr double kl = 1e5; // [L/T^2] edge contractility stiffness per vertex mass
constexpr double kb = 0; // [L^3/T^2] bending stiffness per vertex mass
constexpr double kr = 1e7; // [1/T^2] repulsion stiffness per vertex mass // default 1e7
constexpr double kh = 1e6; // [1/T^2] adhesion stiffness per vertex mass // default 1e6
constexpr double sh = 0.01; // [L] adhesion hardening zone size
constexpr double ss = 0.01; // [L] adhesion softening zone size
constexpr double theta = 0; // [-] fusion threshold. Not fully supported in PolyMorph, keep 0 to be safe.
constexpr double mu = 0; // [-] dynamic friction coefficient
constexpr double rho = 0; // [1/L^2] fluid mass density per vertex mass
constexpr double g = 0; // [L/T^2] gravitational acceleration
constexpr double gl = 0; // [L/T^2] edge gravitational acceleration
constexpr double cv = 20; // [1/T] viscous damping rate
constexpr double cd = 0; // [-] drag coefficient
constexpr double cc = 30; // [1/T] collision damping rate
constexpr double dt = 1e-4; // [T] time step

constexpr std::size_t Nf = 200; // number of output frames
constexpr std::size_t Ns = 1000; // number of time steps between frames // default 1000
constexpr int Nr = 0; // number of rigid polygons (keep as int or things will break)

constexpr double drmax = h + sh + ss; // maximum interaction distance

// -- New PolyMorph parameters --
constexpr bool ADVECTION_DILUTION_EN = false; // enable advection-dilution and calculate velocity field
constexpr bool CHEMOTAXIS_EN = false; // enable chemotaxis. slightly performance enhancing if disabled
constexpr bool RESIZE_GRID_EN = false; // enable resizing of grid if domain changes/grows (not optimized. can be computationally expensive)
constexpr int NUM_SPECIES = 1; // [-] number of diffusable species (defines size of D,p vectors)
constexpr int NUM_KIN = 1; // [-] number of kinetic coefficients (defines size of k vectors)

const std::vector<double> k0 =   {1}; // [?] reaction coefficients background
const std::vector<double> k_mu = {1}; // [?] reaction coefficients mean
const std::vector<double> k_CV = {0}; // [-] reaction coefficients CV
const std::vector<double> D0 =   {32}; // [L^2/T] diffusivity background (recommended to not be zero)
const std::vector<double> D_mu = {32}; // [L^2/T] diffusivity mean
const std::vector<double> D_CV = {0}; // [-] diffusivity CV
const std::vector<double> p0 =   {0}; // [1/(L^2*T)] production rate background (usually zero)
const std::vector<double> p_mu = {1}; // [1/(L^2*T)] production rate mean
const std::vector<double> p_CV = {0}; // [-] production rate CV
const std::vector<double> threshold_mu = {0.001}; // [1/L^2] concentration threshold mean (can use any number of thresholds, but define your ensemble.set_flag accordingly!)
const std::vector<double> threshold_CV = {0.0}; // [-] threshold CV

constexpr double dx = 0.3; // [L] grid spacing for solver
constexpr double anisotropy = 1.0; // [-] diffusion anisotropy (default 1)
constexpr double dist_cutoff_factor = 2.0; // [-] used to cut off lognormal dists at mu*factor to maintain stability. only used for diffusivity D atm.
constexpr double velocity_cutoff_radius = 3.0 * Amax_mu; // [L] radius for IDW interpolation of velocity field. Coupled to Amax to account for length scale
constexpr double domain_bd_stiffness = kr / 2; // [1/T^2] domain boundary stiffness (too high can cause instabilities)

const std::vector<double> chemotaxis_strength = {0}; // force proportional to concentration gradient
const std::vector<int> chem_affect_flag = {1}; // which cell types (flag) to affect with each species

constexpr int RNG_SEED = 90178009; // random number generator seed

// define what to write to vtk files for visualization/debugging. disable things to save space 
namespace Output { 
    constexpr bool u = true; // concentration
    constexpr bool D = true; // diffusion coefficient
    constexpr bool p = true; // production rate
    constexpr bool k = true; // kinetic coefficients
    constexpr bool parent_idx = true; // polygon idx grid
    constexpr bool threshold = true; // threshold
    constexpr bool flag = true; // boolean polygon flag 
    constexpr bool velocity = true; // velocity field
    constexpr bool grad_u = true; // concentration gradient
}; 

#endif