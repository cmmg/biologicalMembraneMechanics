#!/bin/sh
#This file is called submit-script.sh
#SBATCH --partition=pre # default "univ", if not specified
#SBATCH --time=0-23:00:00 # run time in days-hh:mm:ss
#SBATCH --nodes=2# require 2 nodes
#SBATCH --ntasks-per-node=16            # (by default, "ntasks"="cpus")
#SBATCH --mem-per-cpu=4000# RAM per CPU core, in MB (default 4 GB/core)
#SBATCH --job-name="r320h160C1"
#SBATCH --error=error.err
#SBATCH --output=out.out
#SBATCH --mail-user=krudraraju@wisc.edu
#SBATCH --mail-type=FAIL
#Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output

#Now list your executable command (or a string of them).
# Example for non-SLURM-compiled code:
module load boost-1.49 cmake-2.8.11.2 gcc-4.9.0 mpi/gcc/openmpi-2.0.1
export PETSC_DIR=/home/krudraraju/software/petsc/petsc-3.8.4; export PETSC_ARCH=shared-optimized
export PETIGA_DIR=/home/krudraraju/software/PetIGA
export TRILINOS_DIR=/home/krudraraju/software/trilinos/trilinos-12.12.1-Source/install
export LD_LIBRARY_PATH=/home/krudraraju/software/trilinos/trilinos-12.12.1-Source/install/lib

cd /home/krudraraju/workspace/repos/shellpinching/ResultTubeMeshr320h160C1

mpiexec -np 32 ./main -iga_view -ts_monitor -snes_monitor -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -ts_adapt_type none -ts_max_snes_failures 500 -snes_max_it 100 -snes_max_funcs 50000 -snes_type newtontr


