#!/bin/bash
#SBATCH --partition=cn
#SBATCH --job-name=Binarybh_mod
#SBATCH --time=96:00:00

#SBATCH --nodes=16
#SBATCH --ntasks-per-node=128
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1

#SBATCH --account=bg-phys-02
#SBATCH --qos=bg-phys-02

#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

#SBATCH --no-requeue
#SBATCH --mem=251G

# Clean environment and load modules
module purge
module load gmp/6/latest
module load mpfr/latest
module load mpc/latest
module load isl/latest
module load gcc/9/9.2.0
module load openmpi/4/gcc/4.1.6
module load openblas/0/0.3.23-gcc
module load hdf5/1/1.12/latest-gcc-openmpi
module load compiler-rt/2024.1.0
module load tbb/2021.12
module load mkl/latest
module load gsl/2/2.7-gcc

# Export GRChombo paths
export CHOMBO_HOME=/valhalla/projects/bg-phys-02/Marta/Chombo/lib
export GRCHOMBO_SOURCE=/valhalla/projects/bg-phys-02/Marta/GRChombo/Source

export LD_LIBRARY_PATH=/valhalla/projects/bg-phys-02/Marta/software/petsc/petsc-3.14.4/arch-linux-c-opt/lib:$LD_LIBRARY_PATH

# --- OpenMP and Networking (From Working Job) ---
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=master
export UCX_NET_DEVICES=mlx5_0:1

# --- Job Execution ---
cd ${SLURM_SUBMIT_DIR}
mpirun ./Main_UnequalMassBinaryBH4dST3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex params_two.txt



