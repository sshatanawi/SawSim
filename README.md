Sawsan’s Simulation (SawSim): A three-dimensional transient groundwater flow model 
Auther: Sawsan M. Shatanawi
Washington State University

Overview
SawSim is a Fortran-based simulation code developed by Sawsan Shatanawi for modeling three-dimensional transient groundwater flow in a two-layer aquifer system using finite difference approach. The simulation is designed to capture the behavior of an unconfined aquifer (top layer) and a confined aquifer (bottom layer) separated by a non-storage confining unit (3D Quasi system). SawSim adopts MODFLOW principles and uses PETSc’s iterative solvers (PICARD ITERATION + GMRES SOLVER) for solving the large system of equations.. The main objectives of Sawsan are:
* Simulate Transient Flow: Model time-dependent groundwater flow over a specified number of time steps.
* Handle Non-linearity: Use Picard iteration to manage the non-linear aspects of the groundwater flow equations.
* Compute Inter-Cell Conductance: Calculate conductance between grid cells using harmonic means, ensuring stable and accurate representation of spatial flow.
* Support Multiple Sources/Sinks: Incorporate boundary conditions (fixed head or no-flow), wells, river interactions, and upper soil interactions. 
* Enable Coupling with Land Surface Models: Provide a framework where recharge and surface elevation (currently hard-coded for testing) can later be dynamically coupled to outputs from the Noah MP land surface model.
Features
* Domain Description:
o Grid discretization is defined by NCOL × NROW × NLAY cells.
o Two groundwater layers: a top unconfined aquifer and a bottom confined aquifer.
o User-Specified domain size, time stepping (steady or transient), boundary conditions, and source arrays (e.g. wells, rivers, etc.).
o Heterogeneous hydraulic conductivity (K), storage (Sy, Ss) read from data arrays
o Uses harmonic means to compute lateral and vertical conductance between adjacent grid cells.
* Numerical Methods:
o Time-Stepping: Uses a user-specified time step size (DELTAT) over NTSP steps.
o Picard Iteration: Iteratively solves non-linear groundwater flow equations until the convergence tolerance is met.
* PETSc Solver:
o The code creates PETSc matrix and vector objects for assembling and solving the system.
o Use GMRES as Krylov solver with a Jacobi preconditioner to solve the large system of linear equations and Picard iteration to handle the nonlinearity  
* Boundary Conditions and Sources:
o Supports Dirichlet (fixed head) and no-flow boundary conditions.
o No-flow boundaries remain part of the computational domain but do not exchange water with neighboring cells. However, they do not correctly reflect groundwater head since they are not part of the physical flow system. They should be considered outside the domain of interest and are used solely to define model boundaries.  
o External sources include recharge, wells, and river exchanges().
o Does not simulate dewatering (dry) layers; cells remain active with computed heads even if storage approaches zero.
* Integration with Noah MP:
o Recharge and surface elevation routines are currently implemented for testing
o When integrated with Noah MP, the routines GW_recharge() and the surface elevation, NCOL and NROW assignment will be replaced by direct calls to Noah MP outputs to obtain dynamic recharge rates and surface elevations, surface discretization. 
 Sketch
***************************************************************************

(k=1 is top/first layer in z)
(y direction goes downward: j=1 is up, j=2 is below)
(top)
i=1        i=2        i=3
<---- left to right ---->
+----------+----------+----------+
j=1   |(1,1,1)   |(2,1,1)   |(3,1,1)   |
+----------+----------+----------+
j=2   |(1,2,1)   |(2,2,1)   |(3,2,1)   |
+----------+----------+----------+

(k=2 is deeper/second layer in z)
(y direction goes downward: j=1 is up, j=2 is below)
(bottom)
i=1        i=2        i=3
<---- left to right ---->
+----------+----------+----------+
j=1   |(1,1,2)   |(2,1,2)   |(3,1,2)   |
+----------+----------+----------+
j=2   |(1,2,2)   |(2,2,2)   |(3,2,2)   |
+----------+----------+----------+

***************************************************************************
Dependencies
* PETSc: SawSim relies on PETSc (the Portable, Extensible Toolkit for Scientific Computation) for its linear algebra operations and iterative solvers.
o Download and Installation: PETSc can be obtained from the PETSc website. Follow the PETSc installation instructions appropriate for your platform.
Input Files
SawSim reads several input files that provide model parameters and initial conditions:
* GW_config.txt: Contains scalar parameters for the domain (NCOL, NROW, NLAY), time-stepping parameters (NTSP, DELTAT, MAX_ITER, picard_tol), physical properties (soil thickness, vadose thickness, aquifer thicknesses), and file paths. Some of parameters can be use directly from LSM when integrated. 
* Hydraulic Conductivity File (K_FILE): Lists the hydraulic conductivity components (Kx, Ky, Kz) for each grid cell.
* Storage File (STORAGE_FILE): Provides specific yield (SY) for the unconfined aquifer and specific storage (SS) for the confined aquifer.
* Boundary Condition File (BC_FILE): Specifies the boundary condition type and values for each grid cell.
* Source File (SOURCE_FILE): Defines the type of external source (e.g., well, river) for each grid cell.
* Initial Head File (HEAD0_FILE): Contains the initial groundwater head values.
* River File (RIVER_FILE, Optional): Provides river geometry and stage information if river interaction is included.
* Well File (WELL_FILE, Optional): Contains pumping (well) information when applicable.
Solver and Time-Dependent Components
* Time Dependency:
The simulation is fully transient. The state of the system is updated at each time step (controlled by DELTAT and NTSP). Storage terms and external sources are also updated as functions of time.
the steady state can be simulated by setting the specific yield and specific storage to zero or setting NTSP to 1.0
* Iterative Solver:
PETSc’s Krylov subspace solver is employed; GMRES solver with a Jacobi preconditioner. Solver type and parameters (such as tolerance and maximum iterations) can be changed via PETSc options.
* Picard Iteration:
To resolve the non-linear relationship between hydraulic head and conductance (especially in the unconfined aquifer), Picard iteration is used. The iteration continues until the infinity norm of the difference between successive head vectors falls below picard_tol.
Integration with Noah MP
For full coupling with the Noah MP land surface model:
* Recharge:
The subroutine GW_recharge() is designed to eventually accept recharge rates computed by Noah MP rather than using the current hard coded test value.
* Surface Elevation:
The current flat surface elevation is for test purposes only. In a coupled simulation, surface elevation should be read from or computed by Noah MP or it can becalled as input file
* Future Modifications:
Modify the corresponding routines to call Noah MP modules for dynamic water table and recharge computations.
Code Structure
* Main Program (SawSim_main.f90):
Initializes PETSc, calls initialization routines (GW_initialize()), runs the solver (GW_solver()), and finalizes PETSc.
* Module SawSim_declaration:
Contains all global variable declarations (domain geometry, aquifer properties, arrays for boundary conditions, sources, storage, etc.).
* Module SawSim_subroutines:
Implements all subroutines for reading configuration and data files, constructing the grid and harmonic conductance arrays, assembling the PETSc matrix and right-hand side vector, solving the system, and handling input/output.

* Before running any testcase, delete GW_results.csv 
