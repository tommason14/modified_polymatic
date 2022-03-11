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
# Step 2
gmx grompp -f 2.mdp -c 1.gro -t 1.cpt -p ../topol.top -o 2.tpr -maxwarn 2
ibrun mdrun_mpi -s 2.tpr -deffnm 2 -g 2.log
# Step 3
gmx grompp -f 3.mdp -c 2.gro -t 2.cpt -p ../topol.top -o 3.tpr -maxwarn 2
ibrun mdrun_mpi -s 3.tpr -deffnm 3 -g 3.log
# Step 4
gmx grompp -f 4.mdp -c 3.gro -t 3.cpt -p ../topol.top -o 4.tpr -maxwarn 2
ibrun mdrun_mpi -s 4.tpr -deffnm 4 -g 4.log
# Step 5
gmx grompp -f 5.mdp -c 4.gro -t 4.cpt -p ../topol.top -o 5.tpr -maxwarn 2
ibrun mdrun_mpi -s 5.tpr -deffnm 5 -g 5.log
# Step 6
gmx grompp -f 6.mdp -c 5.gro -t 5.cpt -p ../topol.top -o 6.tpr -maxwarn 2
ibrun mdrun_mpi -s 6.tpr -deffnm 6 -g 6.log
# Step 7
gmx grompp -f 7.mdp -c 6.gro -t 6.cpt -p ../topol.top -o 7.tpr -maxwarn 2
ibrun mdrun_mpi -s 7.tpr -deffnm 7 -g 7.log
# Step 8
gmx grompp -f 8.mdp -c 7.gro -t 7.cpt -p ../topol.top -o 8.tpr -maxwarn 2
ibrun mdrun_mpi -s 8.tpr -deffnm 8 -g 8.log
# Step 9
gmx grompp -f 9.mdp -c 8.gro -t 8.cpt -p ../topol.top -o 9.tpr -maxwarn 2
ibrun mdrun_mpi -s 9.tpr -deffnm 9 -g 9.log
# Step 10
gmx grompp -f 10.mdp -c 9.gro -t 9.cpt -p ../topol.top -o 10.tpr -maxwarn 2
ibrun mdrun_mpi -s 10.tpr -deffnm 10 -g 10.log
# Step 11
gmx grompp -f 11.mdp -c 10.gro -t 10.cpt -p ../topol.top -o 11.tpr -maxwarn 2
ibrun mdrun_mpi -s 11.tpr -deffnm 11 -g 11.log
# Step 12
gmx grompp -f 12.mdp -c 11.gro -t 11.cpt -p ../topol.top -o 12.tpr -maxwarn 2
ibrun mdrun_mpi -s 12.tpr -deffnm 12 -g 12.log
# Step 13
gmx grompp -f 13.mdp -c 12.gro -t 12.cpt -p ../topol.top -o 13.tpr -maxwarn 2
ibrun mdrun_mpi -s 13.tpr -deffnm 13 -g 13.log
# Step 14
gmx grompp -f 14.mdp -c 13.gro -t 13.cpt -p ../topol.top -o 14.tpr -maxwarn 2
ibrun mdrun_mpi -s 14.tpr -deffnm 14 -g 14.log
# Step 15
gmx grompp -f 15.mdp -c 14.gro -t 14.cpt -p ../topol.top -o 15.tpr -maxwarn 2
ibrun mdrun_mpi -s 15.tpr -deffnm 15 -g 15.log
# Step 16
gmx grompp -f 16.mdp -c 15.gro -t 15.cpt -p ../topol.top -o 16.tpr -maxwarn 2
ibrun mdrun_mpi -s 16.tpr -deffnm 16 -g 16.log
# Step 17
gmx grompp -f 17.mdp -c 16.gro -t 16.cpt -p ../topol.top -o 17.tpr -maxwarn 2
ibrun mdrun_mpi -s 17.tpr -deffnm 17 -g 17.log
# Step 18
gmx grompp -f 18.mdp -c 17.gro -t 17.cpt -p ../topol.top -o 18.tpr -maxwarn 2
ibrun mdrun_mpi -s 18.tpr -deffnm 18 -g 18.log
# Step 19
gmx grompp -f 19.mdp -c 18.gro -t 18.cpt -p ../topol.top -o 19.tpr -maxwarn 2
ibrun mdrun_mpi -s 19.tpr -deffnm 19 -g 19.log
# Step 20
gmx grompp -f 20.mdp -c 19.gro -t 19.cpt -p ../topol.top -o 20.tpr -maxwarn 2
ibrun mdrun_mpi -s 20.tpr -deffnm 20 -g 20.log
# Step 21
gmx grompp -f 21.mdp -c 20.gro -t 20.cpt -p ../topol.top -o 21.tpr -maxwarn 2
ibrun mdrun_mpi -s 21.tpr -deffnm 21 -g 21.log

# get a structure for futher use
echo 0 | gmx trjconv -f 21.cpt -s 21.tpr -o last_step.gro
echo 0 | gmx trjconv -f 21.cpt -s 21.tpr -o last_step_unwrapped.gro -pbc mol
