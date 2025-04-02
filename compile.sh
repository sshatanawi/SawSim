#!/bin/bash
module load intel/20.2
module load petsc/3.20.0
module load miniconda3
source activate myenv
module list
#rm *genmod* *.o GW.exe
make -f makefile
#./GW.exe
gdb ./GW.exe
#srun -n 1 ./GW.exe -start_in_debugger
