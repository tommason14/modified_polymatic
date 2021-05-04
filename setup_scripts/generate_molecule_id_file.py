#!/usr/bin/env python3
import subprocess as sp
import sys
# map structures to numbers of each
# check all structures have a number assigned to them- shouldn't be an issue here
structs = (sp.check_output(
    "grep '^\s*structure.*xyz' pack.inp | awk '{print $NF}'",
    shell=True).decode("utf-8").strip().split("\n"))
numbers = (sp.check_output("grep '^\s*number' pack.inp | awk '{print $NF}'",
                           shell=True).decode("utf-8").strip().split("\n"))

if len(structs) != len(numbers):
    sys.exit(
        'Error in pack.inp: Make sure all xyz files have a "number" line in the structure block.\n'
        'Cannot assign correct molecule IDs.')

mols = {i: int(n) for i, n in enumerate(numbers, 1)}

with open('moleculeID.txt', 'w') as f:
    f.write('Molecule    Type\n')
    total = 0
    for mol, num in mols.items():
        for _ in range(num):
            total += 1
            f.write(f'{total:<3}    {mol:>6}\n')
