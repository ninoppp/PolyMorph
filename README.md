# PolyMorph
Coupling [PolyHoop](https://www.sciencedirect.com/science/article/pii/S0010465524000511) with a finite difference solver for morphogenetic problems.  

[Bachelor's thesis](https://www.overleaf.com/read/ytngwfdzfdrx#3782c2) project of Nicolas Müller, ETH Zurich. 

## Abstract

PolyMorph is an extension of PolyHoop, a 2D mechanical tissue simulation software for soft particle dynamics developed by Vetter et al. PolyMorph adds a numerical solver for reaction-diffusion equations on top of it which allows, among other things, the modeling of morphogenetic problems with realistic cell shapes. 

The simulation software is programmed in C++. It makes use of finite-difference approximations, explicit Euler time integration, scatter-gather interpolation between Eulerian and Lagrangian reference frames and shared memory parallelism with OpenMP. 

The solver's correctness and convergence was verified and its computational performance quantified through parallel scaling analysis and runtime comparison to the purely mechanical simulation. The coupled system is then used to study patterning precision with noisy morphogen gradients and the results are compared to previous findings.

## Usage Instructions: Quick Start 
```shell
$ git clone https://github.com/MuellerNico/PolyMorph.git
$ cd PolyMorph/
$ sh run.sh
```

Alternatively, a standard makefile is provided if you prefer that over the shell script:

```shell
$ cd PolyMorph/
$ make
$ OMP_NUM_THREADS=8 ./polymorph
```

The script contains a few additional functionalities:

```shell
$ sh run.sh c # cleanup leftover files (e.g. after aborting)
$ sh run.sh m # move output files
$ sh run.sh {filename} # runs {filename}.cpp from the src folder. If omitted main.cpp is used
$ sh run.sh chemotaxis # example: compile and run src/chemotaxis.cpp
```
Remember to set the correct number of threads (default: 8) and desired output folder (default: ``PolyMorph/out``) inside ``run.sh``.

## Features
- Most PolyHoop capabilities (see limitations)
- Finite difference solver for reaction-diffusion equations
- Advection and dilution terms (Can be switched on and off. Computationally more expensive due to velocity field interpolation)
- Coupling via scatter-gather method
- Supports any number of diffusable species
- Supports general, customizable reaction terms (involving any number of interacting species and any number of kinetic coefficients)
- Dirichlet and Neumann boundary conditions
- Rectangular domains (static or expanding with constant speed)
- Chemotaxis mechanism 
- Helper functions to extract measurement data (readout position, precision-zone width, etc.)

## Output
A standard usage produces 3 output types:
- A series of ``.vtp`` frames containing the polygons to be visualized in paraview.
- A series of ``.vts`` frames containing the grid of the finite difference solver. Can be loaded into the same paraview session. 
- ``simulation.cfg`` saves all parameters (set in const.h) used for a simulation run for reproducibility. 

Additionally, depending on the experiment:
- A ``.off`` file for saving the polygon ensemble at a certain time point (usually at the end) to later load as the input/starting point of another simulation. 
- A ``.csv`` file containing measurement data

## Folder Structure
`root`: makefile, run.sh, [Euler](https://scicomp.ethz.ch/wiki/Euler) job scripts, binary executable  
`/include`: Contains all header files which make up the core of this software  
`/src`: Contains cpp files with main() functions; a default testrun (main.cpp) plus a few example experiments.  
`/out`: Default output folder for files  
`/ensemble`: Input files with pre-grown tissues for various tasks  
`/results`: Data and plotting scripts for my personal experiments  
`/bench`: Benchmarking data and plots 

## Source Code Documentation
In the following the most important components of the software are briefly explained. Also refer to my [my thesis report](https://www.overleaf.com/read/ytngwfdzfdrx#3782c2) for some additional implementation details. 

- ``const.h``: Contains all parameters and some settings (like enabling advection-dilution or chemotaxis). Treat this like a configuration-file. 

- ``domain.h``: Defines a very simple rectangular domain within which all calculations will happen. The boundary conditions of the solver will be applied at the boundaries of this domain. A constant growth factor can be set for each direction to turn the static domain into an expanding one.  

- ``ensemble.h``: The ensemble struct contains the core of the original PolyHoop software. It is responsible for the mechanical part of the simulation. It represents biological cells as polygons whose dynamics are governed by their potential energy.  

- ``geometry.h``: Contains the geometric structs used by the ensemble: Point, Vertex, Polygon. 

- ``grid.h``: Defines a simple matrix-like grid data structure used for the solver fields. Under the hood it is a simple nested std::vector. 

- ``solver.h``: The solver handles the reaction-diffusion part of the simulation. It uses central finite difference approximations and explicit euler time stepping. 
Additionally this file defines the available boundary conditions (Dirichlet or Neumann) which can be set individually for each side of the rectangular domain. By default all boundaries are treated as zero-flux.

- ``reaction.h``: Defines a reaction term (e.g. degradation) as a std::function taking two vectors of concentrations and kinetic coefficients and returning a vector with the values of the reaction term. Some example reactions are provided but the user should define their own reactions for their respective experiments/applications. 

- ``interpolator.h``: The interpolator takes care of moving data back and forth between the ensemble (polygons) and the solver (grid). It does this with a scatter-gather process similar to the particle-in-cell (PIC) method otherwise commonly used in plasma physics. Additionally, a velocity field is interpolated using the movement of the cells, more precisely their vertices.     

- ``utils.h``: Various quality-of-life functions.

- ``ensembleController.h``: A collection of functions which interact with the ensemble grouped in this namespace. They provide functionalities used in different experiments (e.g. determine the mean readout position or stop growth cells). These functions should make it easier for the user to build their experiment without having to change the source code of the ensemble. 

## Usage Instructions: Detailed

### Examples
In the ``/src`` folder there are multiple examples to help you get started with PolyMorph. Read the description at the beginning of the file. 

### Adapt to your needs

When simulating your own scenario, the three files you want to primarily consider are:

- ``const.h``: Basically the input file. Here you can set all parameters for your scenario
- ``ensembleController.h``: As described above this file shows you how to interact with the simulation in a very simple way without having to modify any core functionalities of the software. Take the already existing functions in this file as inspiration on how to implement your additions. 
- ``reaction.h``: Definitely have a quick look at how reactions are defined so you understand how the kinetic parameters set in const.h are being used. The file is only 50 sparse lines long. 

### Production lambda function
The function ``std::function<std::vector<bool>(const Polygon&)> Ensemble::is_producing`` defines which cells produce a diffusable species. Producing cells draw a production constant ``p`` from the respective lognormal distribution. Non-producing cells have ``p=0``. The function takes the polygon/cell as an argument and returns a boolean vector of length `NUM_SPECIES` to define which species are produced by this cell. Some examples:

```C++  
// Default behavior: always false (no cell based production)
ensemble.is_producing = [](const Polygon& p) { 
    return std::vector(NUM_SPECIES, false); 
};

// System with only one species produced by the starting cell (index 0)
ensemble.is_producing = [](const Polygon& p) { 
    return std::vector<bool> {
        p.vertices[0].p == 0
    };
};

// System with two species produced by cells at two ends of the domain: 
ensemble.is_producing = [domain](const Polygon& p) { 
    return std::vector<bool>{ 
        p.midpoint().x < domain.x0 + 10, 
        p.midpoint().x > domain.x1 - 10
    }; 
};``
```

### Flagging lambda function
Similarly you can apply an integer flag to certain cells with ``std::function<int(const Polygon&)> Ensemble::set_flag``. This flag can e.g. be used to represent different cell types and model domain boundaries. Example: 

```C++
// System with 1 species. 
// Cells with concentration below their threshold: 1. Otherwise: 0
ensemble.set_flag = [](const Polygon& p) { 
    return p.u[0] < p.threshold[0]; 
};
```


Good luck, have fun!

## Limitations
- Does not support polygon fusion.
- Nested polygons (one inside the other) should technically work fine but might hold unwanted behavior or even crash. Not tested thoroughly. 
- Use of rigid polygons was also not tested thoroughly enough.  
