!************************************************************************
! Sawsan's Simulation (SawSim): A 3D TRANSIENT GROUNDWATER FLOW (UNCONF/CONF)
! Author: Sawsan Shatanawi
!***********************************************************************
! PICARD ITERATION + PCG SOLVER (PETSc)
!
! Features:
!   - NROW x NCOL x NLAY domain
!   - 2 layers (layer1=unconfined, layer2=confined)
!   - Seperated by confining unit does not store water
!   - Transient: time-stepping
!   - Picard iteration for non-linearity 
!   - Conductance between cells uses harmonic mean
!   - Boundary Condition (BC): user-specified fixed head or no-flow
!   - Well or river infiltration: user sets arrays for source/sink
!   - PCG solver from PETSc
!
! Sketch
!(k=1 is top/first layer in z)
!(y direction goes downward: j=1 is up, j=2 is below)
!        i=1        i=2        i=3
!          <---- left to right ---->
!       +----------+----------+----------+
! j=1   |(1,1,1)   |(2,1,1)   |(3,1,1)   |
!(top)  +----------+----------+----------+
! j=2   |(1,2,1)   |(2,2,1)   |(3,2,1)   |
!       +----------+----------+----------+
!
!(k=2 is deeper/second layer in z)
!(y direction goes downward: j=1 is up, j=2 is below)
!        i=1        i=2        i=3
!          <---- left to right ---->
!       +----------+----------+----------+
! j=1   |(1,1,2)   |(2,1,2)   |(3,1,2)   |
!       +----------+----------+----------+
! j=2   |(1,2,2)   |(2,2,2)   |(3,2,2)   |
!       +----------+----------+----------+
!************************************************************************
program SawSim_main
#include <petsc/finclude/petsc.h>
! The main program:
!   1) Initialize PETSc
!   2) gw_initialize => read config & data, build conduction
!   3) SawSim_solver => create local PETSc, time loop, etc.
!   4) finalize PETSc

   use petsc
   use SawSim_declaration
   use SawSim_subroutines
   implicit none

   integer :: iErr

   call PetscInitialize(iErr)
   if (iErr/=0) stop "PetscInitialize failed"

   ! read config, read K/stor/BC/source 
   call GW_initialize()
   
   ! solve
   call GW_solver()

   call GW_cleanup ()

   call PetscFinalize(iErr)
   if (iErr/=0) stop "PetscFinalize failed"
end program SawSim_main
