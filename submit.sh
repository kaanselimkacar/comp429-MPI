#!/usr/bin/env bash
#
# You should only work under the /scratch/users/<username> directory.
#
# Example job submission script
#
# -= Resources =-
#
#SBATCH --job-name=cardiacsim
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --cores-per-socket=16
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --output=output.out

#SBATCH --qos=users
#SBATCH --account=users
module load openmpi/4.0.1
export PATH=/kuacc/apps/openmpi/4.0.1/bin/:$PATH
export LD_LIBRARY_PATH=/kuacc/apps/openmpi/4.0.1/lib/:$LD_LIBRARY_PATH
export PATH=/kuacc/apps/openmpi/4.0.1/include:$PATH

################################################################################
##################### !!! DO NOT EDIT ABOVE THIS LINE !!! ######################
################################################################################

echo "Running Job...!"
echo "==============================================================================="
echo "Running compiled binary..."

export KMP_AFFINITY=verbose,compact
export OMP_NUM_THREADS=4
export SRUN_CPUS_PER_TASK=4
echo "OMP_NUM_THREADS = 4, 4 MPI"

#x=16
#y=1

echo "Parallel version with 1 processes"
mpirun -np 4 -cpus-per-proc 4 --mca orte_base_help_aggregate 0 --mca btl_openib_allow_ib 1 ./cardiacsim -n 1024 -t 100 -y 4 -x 1 -o 4

#echo "Parallel version with 16 processes"
#mpirun -np 16 --mca orte_base_help_aggregate 0 --mca btl_openib_allow_ib 1 ./cardiacsim -n 1024 -t 100 -y $y -x $x

#echo "Parallel version with 4 processes"
#mpirun -np 4 --mca orte_base_help_aggregate 0 --mca btl_openib_allow_ib 1 ./cardiacsim -n 1024 -t 100 -y 4 -x 1



#echo "Parallel version with 8 processes"
#mpirun -np 8 ./cardiacsim -n 1024 -t 100 -y 8 -x 1

#echo "Parallel version with 8 processes"
#mpirun -np 8 ./cardiacsim -n 1024 -t 100 -y 4 -x 2


#echo "Parallel version with 8 processes"
#mpirun -np 8 ./cardiacsim -n 1024 -t 100 -y 2 -x 4

#echo "Parallel version with 8 processes"
#mpirun -np 8 ./cardiacsim -n 1024 -t 100 -y 1 -x 8

#echo "Parallel version with 15 processes"
#mpirun -np 15 ./cardiacsim -n 1024 -t 100 -y 1 -x 15

#echo "Parallel version with 14 processes"
#mpirun -np 14 ./cardiacsim -n 1024 -t 100 -y 1 -x 14

#echo "Parallel version with 13 processes"
#mpirun -np 13 ./cardiacsim -n 1024 -t 100 -y 1 -x 13

#echo "Parallel version with 12 processes"
#mpirun -np 12 ./cardiacsim -n 1024 -t 100 -y 4 -x 3

#echo "Parallel version with 11 processes"
#mpirun -np 11 ./cardiacsim -n 1024 -t 100 -y 11 -x 1

#echo "Parallel version with 11 processes"
#mpirun -np 11 ./cardiacsim -n 1024 -t 100 -y 1 -x 11

#echo "Parallel version with 9 processes"
#mpirun -np 9 ./cardiacsim -n 1024 -t 100 -y 3 -x 3

#echo "Parallel version with 9 processes"
#mpirun -np 9 ./cardiacsim -n 1024 -t 100 -y 9 -x 1

#echo "Parallel version with 9 processes"
#mpirun -np 9 ./cardiacsim -n 1024 -t 100 -y 1 -x 9

#echo "Parallel version with 8 processes"
#mpirun -np 8 ./cardiacsim -n 1024 -t 100 -y 1 -x 8


#echo "Parallel version with 8 processes"
#mpirun -np 8 ./cardiacsim -n 1024 -t 100 -y 2 -x 4
#echo "Parallel version with 9 processes"
#mpirun -np 8 ./cardiacsim -n 1024 -t 100 -y 3 -x 3

#echo "Parallel version with 10 processes"
#mpirun -np 10 ./cardiacsim -n 1024 -t 100 -y 10 -x 1

#echo "Parallel version with 16 processes"
#mpirun -np 16 ./cardiacsim -n 1024 -t 100 -y 16 -x 1


#Different configuration of MPI+OpenMP
#[1 + 16] [2 + 8] [4 + 4] [8 + 2] [ 1 + 16]

#export KMP_AFFINITY=verbose,compact

#echo "MPI1 + OMP16"
#export OMP_NUM_THREADS=16
#mpirun -np 1 ./cardiacsim -n 1024 -t 100 -y 1 -x 1

#echo "MPI2 + OMP8"
#export OMP_NUM_THREADS=8
#export SRUN_CPUS_PER_TASK=8
#mpirun -np 2 -cpus-per-proc 8 ./cardiacsim -n 1024 -t 100 -y 2 -x 1
