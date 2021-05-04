#!/bin/bash
#SBATCH -J CR61-test                
#SBATCH -o polymerisation.out
#SBATCH -e polymerisation.e%j            
#SBATCH -N 1                   
#SBATCH -n 8 
#SBATCH -p short,comp
#SBATCH -t 10:00            

module load vmd
module load openmpi

polymatic_autogenerate.sh 
cp -r ../polymatic/scripts .

export LAMMPS_EXEC="mpirun -np $SLURM_NTASKS ~/p2015120004/apps/clammps/build/lmp_mpi"
python3 ../polymatic/polym_loop.py --controlled
