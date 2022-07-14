#!/bin/bash
#SBATCH -J 21-step
#SBATCH -e run.%j.err
#SBATCH -o run.%j.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p skx-dev
#SBATCH -t 2:00:00

export OMP_NUM_THREADS=48

module load gromacs/2019.6

# Step 1
gmx grompp -f 1.mdp -c ../npt.gro -p ../topol.top -o 1.tpr -maxwarn 2
ibrun mdrun_mpi -s 1.tpr -deffnm 1 -g 1.log

# Steps 2-21
for step in {2..21}
do
  prev=$((step-1))
  gmx grompp -f $step.mdp -c $prev.gro -t $prev.cpt -p ../topol.top -o $step.tpr -maxwarn 2
  ibrun mdrun_mpi -s $step.tpr -deffnm $step -g $step.log
done

# save a structure for futher use
echo 0 | gmx trjconv -f 21.cpt -s 21.tpr -o last_step.gro
echo 0 | gmx trjconv -f 21.cpt -s 21.tpr -o last_step_unwrapped.gro -pbc mol
