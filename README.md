# Modified Polymatic

Here the original [Polymatic](https://nanohub.org/resources/17278) code written
by Lauren Abbott has been modified and made compatible with python 3. Additional
functionality has been added to run the algorithm until one single polymer chain
has been formed, using VMD to count the number of fragments in the LAMMPS
datafile at each step.
Additionally, free energy barriers may be used to partially control the
polymerisation. Reactions of high energy are disfavoured, resulting in a more realistic polymer.

Before using this modified program, I strongly recommend reading the manual distributed with the original Polymatic code, to see how the code works, or read the original [paper](https://link.springer.com/article/10.1007/s00214-013-1334-z).

# Setup process

Several steps are required to set up a system for use with Polymatic.
1) Create a simulation box of monomers
2) Create a LAMMPS datafile
3) Add additional parameters of bonds/angles/dihedrals that are going to be
   formed during the polymerisation
4) Create a types.txt file containing the atom types for each atom, bond, angle,
   dihedral and improper dihedral in the LAMMPS datafile.
5) Add correct partial charges (optional)

Methods to achieve these steps are described below. These methods use the
scripts inside the [setup_scripts](setup_scripts) folder, so make these scripts accessible by modifying your `$PATH` variable.
See [example-with-detailed-setup/run.sh](example-with-detailed-setup/run.sh) for an example shell script.

## Creating a simulation box

The simplest way to do this is to use
[Packmol](http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml). This tool uses
xyz or pdb files and generates a simulation box using an input file such as
this:

```
# pack.inp
tolerance 2.0
output pack.xyz
filetype xyz
structure stysulf.xyz
  number 20
  inside cube 0 0 0 40
end structure
structure dvb.xyz
  number 20
  inside cube 0 0 0 40
end structure
```

Save the file above to `pack.inp`, then create the box with `$ packmol < pack.inp`.

## Creating a datafile

There are many ways to create a LAMMPS datafile, including
[moltemplate](https://www.moltemplate.org/). Another method is using the VMD
Topotools plugin. The `lmp_gaff.py` script reads an xyz file (such as the one
generated by Packmol) then uses VMD and
Topotools to generate a datafile with no parameters. The script then searches a
forcefield file for the appropriate parameters. An example forcefield file is
the [gaff.ff](example-with-detailed-setup/gaff.ff) file used in the examples - this was created by converting the
gaff.lt file distributed with moltemplate.

A datafile will be created by running `$ lmp_gaff.py pack.xyz gaff.ff`.

At this point, any forcefield may be chosen. For example, if you wish to use
OPLS, you can use this [lmp_opls_cvff.py](https://github.com/tommason14/scripts/blob/master/chem/lammps/create_opls_jobs/lmp_opls_cvff.py) script and the forcefield file [here](https://github.com/tommason14/scripts/blob/master/chem/lammps/create_opls_jobs/via-topotools), along with the python script in that folder to generate additional parameters. Note that this forcefield uses the oplsaa.ff force field distributed with
[fftool](https://github.com/paduagroup/fftool), but with cvff improper dihedrals.

Any other method of generating a LAMMPS datafile will also work here.

## Adding additional parameters

Polymatic requires that parameters for bond, angles and dihedrals formed during
the polymerisation are present in the datafile before the script is run. To do
this, you can manually create dimers of each possible reaction that may occur,
and then add in the additional parameters to the original datafile.
However, this is extremely time-consuming and prone to user error.

Instead, we can use the `polym.in` file required by Polymatic that outlines the
atom types that may bond (see examples), and the atom types that will be present after bonding.
Using this file along with the LAMMPS datafile, all atom types that may form new
bonds, angles or dihedrals can be found. All combinations of these types are
then considered using the python `itertools` package, and the `gaff.ff` force
field file is used to find the additional
parameters.
To do this, run
`$ add_addtional_parameters.py -l pack.lmps -f gaff.ff -p polym.in`.

## Creating a types.txt file

Polymatic uses a types.txt file, describing the atom types involved in each
bond, angle and dihedral. To create this file from the LAMMPS datafile, run
`$ polymatic_types.py pack.lmps > types.txt`.

## Adding partial charges (optional)

The GAFF forcefield does not contain partial charges. As a result, here I have
included a script that looks inside the pack.inp file for the order of atoms in
the datafile, and extracts charges from GAMESS geodesic log files of the same
name as each xyz file. 
For example, for a monomer called dvb.xyz, include a GAMESS log file 
called dvb.log.
Then run
`$ add_correct_charges.py pack.lmps`. 

# Setup process using the GAFF forcefield

The setup process has been simplified for the GAFF forcefield.

Simply ensure you have the following files in the current directory:

- gaff.ff (forcefield file, feel free to modify this as necessary)
- polym.in (describing how Polymatic should be used)
- pack.inp (packmol file describing a simulation box of monomers)

Make sure that all xyz files are labelled according to the atom types in
gaff.ff, and then run `polymatic_autogenerate.sh`. This will perform steps 1-5
above, with no user input required.

# Running Polymatic

To perform the polymerisation, include the `pack.lmps` file in a directory along
with the `polym.in` and `types.txt` files.

An example `polym.in` file looks like:

```
link LC1,C3 LC2,C3
charge +0.3 -0.3
cutoff 6.0
intra true 5
```

specifying that LC1 and LC2 atoms within 6 angstroms can bond, and crosslinking between atoms at least 5 bonds apart is allowed.

Polymatic relies on several LAMMPS input files throughout the polymerisation.
These are provided inside [this](polymatic/scripts) folder, and need to copied
to the directory where you wish to run the job.

Polymatic is run by using this [polym_loop.py](polymatic/polym_loop.py) file.
Several options can be supplied. 

1. By default, a minimisation of the pack.lmps file is performed, then a
   random polymer is formed.
2. To retain files for each polymerisation step, add a `--keep` flag.
3. To use the original packmol-generated file with no minimisation, rename
   pack.lmps to data.lmps and add a `--no-minimise` flag.
4. If you are forming a co-polymer and wish to control how the polymer is
   formed, create a file of free energy barriers of each propagation step 
   and name it `barriers.in`. 
   For example: 

```
mol1 stysulf    
mol2 dvb        
k11  77.70758788
k12  77.68702204
k21  83.54565331
k22  91.49469955
   ```
   In this example, the reaction between two styrene sulfonate monomers
   has a free energy barrier of 77.71 kJ/mol, while the energy barrier to two
   divinylbenzene units reacting is 91.49 kJ/mol.

   Then add a `--controlled` flag in the call to
   `polym_loop.py`.
   After all of one type of monomer has been consumed, a 20 ps NVT run at 900 K
   is performed to rearrange the molecules, using `polymatic/scripts/shuffle.in`. To
   turn off this behaviour, add the `--controlled --no-shuffle` flags.

Lastly, the LAMMPS executable is defined using a `LAMMPS_EXEC` environmental
variable. In practice, just export `LAMMPS_EXEC` in the shell script before
calling `python3 polym_loop.py`.
For example:

```
export LAMMPS_EXEC='ibrun lmp_stampede'
python3 polymatic/polym_loop.py --no-minimise
```

# After polymerisation

As the polymer has been created from a low density initial system, a series of compression/decompression schemes are commonly used to increase the density to a more accurate value. A 21-step equilibration scheme was presented in the original Polymatic paper, found in [polymatic/21-step](polymatic/21-step).

After the 21-step equilibration, the polymer is ready for use. 
To convert a polymer created with the GAFF force field to gromacs, see this [repo](https://github.com/tommason14/lammps_gaff_to_gromacs). For conversion with other forcefields, [Intermol](https://github.com/shirtsgroup/InterMol) may be worth a look.
