# PolyMorph
Coupling PolyHoop with an FDM-solver for morphogenetic problems. 

## Usage instructions
Quick start: 
```shell
cd PolyMorph
sh run.sh
```

Alternatively a standard makefile is provided:

```shell
cd PolyMorph
make
OMP_NUM_THREADS=8 ./polymorph
```

The script contains a few additional functionalities:

```shell
sh run.sh c # cleanup leftover files (e.g. after aborting)
sh run.sh m # move output files
sh run.sh {filename} # compiles {filename}.cpp from the src folder. If omitted main.cpp is used
```

## Output
A standard usage produces 3 output types:
- A series of .vtp files containing the polygons to be visualized in paraview.
- A series of .vts files containing the grid of the finite difference solver. Can be loaded into the same paraview session. 
- simulation.cfg contains all parameters (from const.h) used for this simulation to allow reproducibility. 
- Optionally .off for saving the polygon ensemble at a certain time point (usually at the end) to later load as the input/starting point of another simulation. 

## Structure
- include: Contains all header files and the core of this software
- src: Contains cpp files with main() functions. A default testrun (main.cpp) plus a few example experiments. 
- out: Default output folder for output files (.csv, .vtp, .vts, .off, .cfg)
- makefile, run.sh, euler job scripts

## Documentation
In the following the most important components of the code are briefly explained

### const.h
Contains all parameters and some settings (like enabling advection-dilution or chemotaxis). Treat this like a configuration-file. 

### domain.h
Defines a very simple rectangular domain within which all calculations will happen. The boundary conditions of the solver will be applied at the boundaries of this domain. The size of the domain can be changed by setting a growth factor for each direction or  (not advised) by directly changing its size and/or coordinates.  

### ensemble.h
The ensemble struct contains the core of the original PolyHoop software. It is responsible for the mechanical part of the simulation. 

### geometry.h
Contains the geometric structs used in the ensemble: Point, Vertex, Polygon. 

### grid.h
Defines a simple matrix-like grid data structure used for the solver fields. Under the hood it is a simple nested std::vector. 

### solver.h
The solver handles the reaction-diffusion part of the simulation. It uses central finite difference approximations and explicit euler timestepping. 
Additionally this file defines the available boundary conditions (Dirichlet or Neumann) which can be set individually for each side of the rectangular domain. By default all boundaries are zero-flux.

### reaction.h
Defines a reaction term (e.g. degradation) as a std::function taking two vectors of concentrations and kinetic coefficients and returning a vector with the values of the reaction term. Some example reactions are provided but the user should define their own reactions for their respective experiments/applications. 

### interpolator.h
The interpolator takes care of moving data back and forth between the ensemble (polygons) and the solver (grid). It does this with a scatter-gather process similar to the particle-in-cell (PIC) method for N-body problems.  

### utils.h
Various quality-of-life functions.

### ensembleController.h
A collection of functions which interact with the ensemble. They provide functionalities used in different experiments (e.g. determine the mean readout position or stop growth of all cells). These functions should make it easier for the user to build their experiment without having to change the source code of the ensemble. 




