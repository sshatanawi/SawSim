# Sawsan's Simulation (SawSim): A Three-Dimensional Transient Groundwater Flow Model

**Author:** Sawsan M. Shatanawi  
**Institution:** Washington State University

---

## Overview

**SawSim** is a groundwater simulation model developed in Fortran for simulating three-dimensional transient groundwater flow in a two-layer aquifer system. It is based on a finite difference formulation and follows principles from MODFLOW while leveraging PETSc's powerful linear solver capabilities.

The upper layer represents an **unconfined aquifer**, and the lower layer is a **confined aquifer**. These are separated by a non-storage confining unit (3D-quasi system).

SawSim is designed for researchers and students working in hydrogeology, groundwater-surface water interactions, and land surface model integration, especially those interested in enhancing land-atmosphere interactions by accounting for groundwater dynamics.

---

## Key Objectives

- **Simulate Transient Flow:** Time-evolving groundwater head and storage.
- **Handle Non-linearity:** Using Picard iteration for unconfined flow.
- **Capture Inter-Cell Conductance:** Calculate conductance between grid cells using harmonic means, ensuring stable and accurate representation of lateral and vertical flow.
- **Support Multiple Inputs:** Incorporate boundary conditions (fixed head or no-flow), wells, river interactions, and upper soil interactions.
- **Enable Coupling with LSMs:** Designed for future integration with Noah-MP.

---

## Model Features

### Domain Setup

- Structured 3D finite difference grid (NCOL x NROW x NLAY)
- Layer 1: Unconfined aquifer
- Layer 2: Confined aquifer
- Heterogeneous hydraulic conductivity (K) and storage parameters (Sy, Ss)
- Fully transient solver with flexible time stepping (NTSP, DELTAT)
- Uses harmonic means to compute lateral and vertical conductance between adjacent grid cells.

### Numerical Methods

- **Time-Stepping:** Explicit control over time step and simulation duration
- **Non-linear Solver:** Picard iteration for water table dynamics
- **Linear Solver:** PETSc GMRES with Jacobi preconditioning

### Boundary Conditions

- **Dirichlet:** Fixed-head boundaries
- **No-flow:** Internalized but do not exchange water (used to set domain edges).No flow boundaries remain part of the computational domain but do not exchange water with neighboring cells. However, they do not correctly reflect groundwater head since they are not part of the physical flow system. They should be considered outside the domain of interest and are used solely to define model boundaries.  
- **External Sources:** Recharge, wells, river-aquifer interaction

### Special Notes

- Cells do not go dry but may approach zero storage: Does not simulate dewatering (dry) layers; cells remain active with computed heads even if storage approaches zero.
- Recharge and surface elevation currently hardcoded (for testing)
- Ready for dynamic coupling with LSMs

---

## Dependencies

- **PETSc** ([https://petsc.org](https://petsc.org)): Used for all matrix and solver operations
- Requires MPI and compatible Fortran compiler (e.g., `ifort`)

---

## Input Files

- `GW_config.txt`: Scalar model parameters (domain size, time steps, etc.)
- `K_FILE`: Hydraulic conductivity values (Kx, Ky, Kz)
- `STORAGE_FILE`: Specific yield/storage values (Sy, Ss)
- `BC_FILE`: Boundary condition types and values
- `SOURCE_FILE`: Source/sink types (wells, rivers)
- `HEAD0_FILE`: Initial groundwater head
- `RIVER_FILE` (optional): River geometry and stage
- `WELL_FILE` (optional): Pumping rates and well locations

---

## Solver Behavior

- Full transient simulations with time-dependent storage, sources, and boundary conditions
- Convergence via Picard iteration using PETSc GMRES solver
- Infinity norm-based stopping criteria (`picard_tol`)
- Optionally run steady-state (set `NTSP=1` or zero storage terms)

---

## Integration with Noah-MP (In Progress)

SawSim is designed to support future coupling with the Noah-MP Land Surface Model:

- `GW_recharge()` subroutine will accept Noah-MP recharge values instead of fixed input
- Surface elevation can be passed dynamically from Noah-MP grid outputs
- NCOL/NROW/grid resolution can be synchronized with Noah-MP domain settings

---

## Sketch of Grid System

```
Layer 1 (k=1) - Unconfined
      i=1        i=2        i=3
      <---- left to right ---->
      +----------+----------+----------+
j=1   |(1,1,1)   |(2,1,1)   |(3,1,1)   |
      +----------+----------+----------+
j=2   |(1,2,1)   |(2,2,1)   |(3,2,1)   |
      +----------+----------+----------+

Layer 2 (k=2) - Confined
      i=1        i=2        i=3
      <---- left to right ---->
      +----------+----------+----------+
j=1   |(1,1,2)   |(2,1,2)   |(3,1,2)   |
      +----------+----------+----------+
j=2   |(1,2,2)   |(2,2,2)   |(3,2,2)   |
      +----------+----------+----------+
```

---

## Code Structure

- `SawSim_main.F90`: Initializes and runs the model
- `SawSim_declaration.F90`: Module with variable declarations
- `SawSim_subroutines.F90`: All routines for reading inputs, solving, and I/O

---

## Before Running

- Always delete or clear `GW_results.csv` before running a new test

- Export PETSc paths in your terminal:
  ```bash
  export PETSC_DIR=/path/to/petsc
  export PETSC_ARCH=arch-linux-c-opt
  ```
- Then compile:
  ```bash
  cd hrldas
  make gw_sawsim
  ```

---

## Who Can Use SawSim?

- **Hydrologists & Hydrogeologists** exploring confined/unconfined dynamics
- **Model developers** extending land surface models with groundwater components
- **Students & Researchers** learning groundwater flow and solver coupling

Noah-MP developers can also use SawSim to evaluate confined aquifer dynamics or enhance TWS closure through groundwater feedback.

---

For feedback, support, or collaboration inquiries, please contact:  
**Sawsan M. Shatanawi**  
Washington State University
