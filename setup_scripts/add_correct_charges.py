#!/usr/bin/env python3

"""
File: add_charges.py
Author: Tom Mason
mail: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Search pack.inp for the required number of atoms, and takes in a GAMESS
geodesic charge log file of the monomer(s) used to create the lammps system.
"""
from autochem import read_file, eof
import os
from glob import glob
import subprocess as sp
import sys
import re

if not 2 <= len(sys.argv) <= 3 or sys.argv[1] == '-h':
    sys.exit(f'Syntax: {os.path.basename(__file__)} lammps_file [output_file]')

lammps_datafile = sys.argv[1]
if len(sys.argv) == 2:
    output = lammps_datafile
else:
    output = sys.argv[2]

def find_structures():
    # map structures to numbers of each
    structs = (
        sp.check_output("grep '^\s*structure.*xyz' pack.inp | awk '{print $NF}'", shell=True)
        .decode("utf-8")
        .strip()
        .split("\n")
    )
    numbers = (
        sp.check_output("grep 'number' pack.inp | awk '{print $NF}'", shell=True)
        .decode("utf-8")
        .strip()
        .split("\n")
    )

    return {s.replace('.xyz', ''): {"number": int(n)} for s, n in zip(structs, numbers)}

def detect_charges(mols):
    charge_logs = glob('*.log')
    for mol in mols.keys():
        for chargefile in charge_logs:
            if mol in chargefile:
                mols[mol]['charge_log'] = chargefile
                print(f'Charges for {mol} taken from {chargefile}')
    return mols

def ask_user_for_charges(mols):
    """
    Reads cwd for files ending in .log, then asks user to assign them to the molecules
    stored in the mols dictionary
    """
    charge_logs = glob('*.log')
    menu = ''
    for idx, log in enumerate(charge_logs, 1):
        menu += f'({idx}) {log}\n'
    print('Log files found:')
    print(menu)
    acceptable_responses = [str(i) for i, _ in enumerate(charge_logs, 1)]
    for mol in mols.keys():
        acceptable = False
        while not acceptable:
            calc = input(f'Enter # of geodesic calc for {mol}: ')
            if calc in acceptable_responses:
                acceptable = True
        mols[mol]['charge_log'] = charge_logs[int(calc) - 1]
    return mols

def add_charges(mols):
    """
    Search each charge_log in mols for partial charges that are printed in the last 
    20% of the file. Doesn't account for the atom type, so the xyz used to generate the
    partial charges must contain atoms in the same order as the xyz used in the packmol script.
    """
    for data in mols.values():
        charges = []
        for line in eof(data['charge_log'], 0.2):
            if re.search('^\s[A-z]+(\s+-?[0-9]+\.[0-9]+){2}$', line):
                charges.append(float(line.split()[1]))
        data['charges'] = charges
    return mols

def generate_charge_sequence(mols):
    """
    Loops over the mols dict to generate the new charge sequence to be added into the 
    Atoms section of the lammps data file.
    Order is the same as the order that packmol generates (all of type1, then all of type 2 etc...)
    """
    charges = []
    for data in mols.values():
        for _ in range(data["number"]):
            charges += data['charges']
    return charges
    
def change_datafile(charges, datafile, output):
    b4_atoms = []
    atoms = []
    after_atoms = []
    b4_atoms_bool = True
    found_atoms_bool = False
    after_atoms_bool = False
    for line in read_file(datafile):
        if re.search('\s*Atoms', line):
            b4_atoms_bool = False
            found_atoms_bool = True
        if re.search('^\s*Bonds', line):
            found_atoms_bool = False
            after_atoms_bool = True
        if b4_atoms_bool:
            b4_atoms.append(line)
        if found_atoms_bool:
            atoms.append(line)
        if after_atoms_bool:
            after_atoms.append(line)

    newatoms = []
    for line, charge in zip(atoms[2:-1], charges):
        line = line.split()
        line[3] = str(charge) # change column
        newatoms.append(' '.join(line) + '\n')
    atoms = atoms[:2] + newatoms + [atoms[-1]]

    with open(output, 'w') as f:
        f.writelines(b4_atoms)
        f.writelines(atoms)
        f.writelines(after_atoms)

def main():
    mols = find_structures()
    mols = detect_charges(mols)
    # mols = ask_user_for_charges(mols) 
    mols = add_charges(mols)
    charge_sequence = generate_charge_sequence(mols)
    change_datafile(charge_sequence, lammps_datafile, output)

main()
