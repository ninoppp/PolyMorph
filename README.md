# PolyMorph
Coupling PolyHoop with an FDM-solver for morphogenetic problems.  
Bachelor's thesis project of Nicolas Müller, ETH Zurich. 

## Usage Instructions: Quick Start 
```shell
$ cd PolyMorph
$ sh run.sh
```

Alternatively a standard makefile is provided:

```shell
$ cd PolyMorph
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
Remember to set the correct number of threads (default: 8) and desired output folder (default: ``./out``) inside ``run.sh``.

## Features
- Most PolyHoop capabilities (see limitations)
- Coupled FDM-solver for reaction-diffusion equations
- Advection-dilution terms (Can be switched on and off. Computationally more expensive due to velocity field interpolation)
- Supports any number of diffusable species
- Supports general, customizable reaction terms (involving any number of interacting species and any number of kinetic coefficients)
- Dirichlet and Neumann boundary conditions
- Chemotaxis mechanism
- Helper functions to extract measurement data (e.g. readout position, precision-zone width)

## Output
A standard usage produces 3 output types:
- A series of ``.vtp`` frames containing the polygons to be visualized in paraview.
- A series of ``.vts`` frames containing the grid of the finite difference solver. Can be loaded into the same paraview session. 
- ``simulation.cfg`` saves all parameters (set in const.h) used for this simulation run for reproducibility. 

Additionally, depending on the experiment:
- A ``.off`` file for saving the polygon ensemble at a certain time point (usually at the end) to later load as the input/starting point of another simulation. 
- A ``.csv`` file containing measurement data

## Folder Structure
`root`: makefile, run.sh, euler job scripts, binary executable  
`/include`: Contains all header files which make up the core of this software  
`/src`: Contains cpp files with main() functions; a default testrun (main.cpp) plus a few example experiments.  
`/out`: Default output folder for files 

## Source Code Documentation
In the following the most important components of the code are briefly explained. For a more detailed explanation see the /documentation (ToDo: doxygen) and the report (ToDo: link or pdf). 

- ``const.h``: Contains all parameters and some settings (like enabling advection-dilution or chemotaxis). Treat this like a configuration-file. 

- ``domain.h``: Defines a very simple rectangular domain within which all calculations will happen. The boundary conditions of the solver will be applied at the boundaries of this domain. The size of the domain can be changed by setting a growth factor for each direction or  (not advised) by directly changing its size and/or coordinates.  

- ``ensemble.h``: The ensemble struct contains the core of the original PolyHoop software. It is responsible for the mechanical part of the simulation. 

- ``geometry.h``: Contains the geometric structs used in the ensemble: Point, Vertex, Polygon. 

- ``grid.h``: Defines a simple matrix-like grid data structure used for the solver fields. Under the hood it is a simple nested std::vector. 

- ``solver.h``: The solver handles the reaction-diffusion part of the simulation. It uses central finite difference approximations and explicit euler timestepping. 
Additionally this file defines the available boundary conditions (Dirichlet or Neumann) which can be set individually for each side of the rectangular domain. By default all boundaries are treated as zero-flux.

- ``reaction.h``: Defines a reaction term (e.g. degradation) as a std::function taking two vectors of concentrations and kinetic coefficients and returning a vector with the values of the reaction term. Some example reactions are provided but the user should define their own reactions for their respective experiments/applications. 

- ``interpolator.h``: The interpolator takes care of moving data back and forth between the ensemble (polygons) and the solver (grid). It does this with a scatter-gather process similar to the particle-in-cell (PIC) method for N-body problems.  

- ``utils.h``: Various quality-of-life functions.

- ``ensembleController.h``: A collection of functions which interact with the ensemble. They provide functionalities used in different experiments (e.g. determine the mean readout position or stop growth of all cells). These functions should make it easier for the user to build their experiment without having to change the source code of the ensemble. 

## Usage Instructions: Detailed

ToDo


## Limitations
- Does not support polygon fusion.
- Nested polygons (one inside the other) should technically work fine but might hold unwanted behavior or even crash. Not tested thoroughly. 
- Use of rigid polygons was also not tested thoroughly enough. 
