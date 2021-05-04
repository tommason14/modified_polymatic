#!/usr/bin/env python3

import argparse
import sys
import re

if len(sys.argv) < 2 or sys.argv[1] == '-h':
    sys.exit('Syntax: polymatic_types.py <lammps_input>')
inp = sys.argv[1]

found_masses = False
found_bonds = False
found_angles = False
found_dihedrals = False
found_impropers = False

is_data = lambda line: len(line.split()) > 0

atoms = {}
bonds = {}
angles = {}
dihedrals = {}
impropers = {}

with open(inp) as inpfile:
    for line in inpfile:
        if is_data(line):
            if found_masses and re.search('^\s*[A-Z]', line):
                found_masses = False
            if found_bonds and re.search('^\s*[A-Z]', line):
                found_bonds = False
            if found_angles and re.search('^\s*[A-Z]', line):
                found_angles = False
            if found_dihedrals and re.search('^\s*[A-Z]', line):
                found_dihedrals = False
            if found_impropers and re.search('^\s*[A-Z]', line):
                found_impropers = False
            if "Masses" in line:
                found_masses = True
                continue
            if "Bond Coeffs" in line:
                found_bonds = True
                continue
            if "Angle Coeffs" in line:
                found_angles = True
                continue
            if "Dihedral Coeffs" in line:
                found_dihedrals = True
                continue
            if 'Improper Coeffs' in line:
                found_impropers = True
                continue
            if found_masses:
                if len(line.split()) > 2:
                    count, *_, name = line.split()
                    name = name.replace("#", "").strip()
                    atoms[count] = name
            if found_bonds:
                if len(line.split()) > 2:
                    count, *_, name = line.split()
                    name = name.replace("#", "").strip()
                    bonds[count] = name
            if found_angles:
                if len(line.split()) > 2:
                    count, *_, name = line.split()
                    name = name.replace("#", "").strip()
                    angles[count] = name
            if found_dihedrals:
                if len(line.split()) > 2:
                    count, *_, name = line.split()
                    name = name.replace("#", "").strip()
                    dihedrals[count] = name
            if found_impropers:
                if len(line.split()) > 2:
                    count, *_, name = line.split()
                    name = name.replace("#", "").strip()
                    impropers[count] = name

alldata = {
    "atom types": atoms,
    "bond types": bonds,
    "angle types": angles,
    "dihedral types": dihedrals,
    "improper types": impropers,
}

# if no impropers...
alldata = {k: v for k,v in alldata.items() if len(v) > 0}

# print to stdout
for string, data in alldata.items():
    print(string)
    for count, item in data.items():
        print(f"{count:<5} {item.replace('-', ',')}")
    print("#")
