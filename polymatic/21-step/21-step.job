#!/bin/sh
#PBS -P k96
#PBS -l storage=scratch/k96+gdata/k96
#PBS -l mem=192gb
#PBS -l ncpus=48
#PBS -l jobfs=300gb
#PBS -l walltime=2:00:00
#PBS -l wd

module load openmpi/4.0.2
mpirun -np $PBS_NCPUS /g/data/k96/apps/lammps-3Mar20/bin/lmp_nci -in 21-step.in -v tk 300 >& lammps.out
