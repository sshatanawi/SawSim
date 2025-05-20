#!/bin/bash
module load intel/20.2
module load petsc/3.20.0
module load miniconda3
source activate myenv
module list
#rm *genmod* *.o GW.exe
make -f Makefile
#./SawSim_GW.exe
gdb ./SawSim_GW.exe
#srun -n 1 ./GW.exe -start_in_debugger
