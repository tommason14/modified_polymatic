#!/bin/bash
#SBATCH -J CR61-test                
#SBATCH -A your_account
#SBATCH -o slurm.o%j
#SBATCH -e slurm.e%j            
#SBATCH -N 1                   
#SBATCH -n 48
#SBATCH -p skx-dev            
#SBATCH -t 20:00            

module load intel/18.0.2
module load impi/18.0.2
module load lammps/9Jan20
module load vmd
export LAMMPS_EXEC="ibrun lmp_stampede"

setup(){
  packmol < pack.inp
  lmp_gaff.py pack.xyz gaff.ff # will give incorrect charges
  add_additional_params.py -l pack.lmps -f gaff.ff -p polym.in
  add_correct_charges.py pack.lmps
  # check charge is now correct
  sed -n '/Atoms/,/Bonds/p' pack.lmps | grep '^[0-9]' | awk '{charge+=$4} END {printf "Total charge now = %.5f \n", charge}'
  polymatic_types.py pack.lmps > types.txt
  generate_molecule_id_file.py 
}

polymerise(){
  cp -r ../polymatic/scripts .
  python3 ../polymatic/polym_loop.py --controlled
}

time setup > setup.txt
time polymerise > polymerisation.txt
