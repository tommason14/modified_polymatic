#!/usr/bin/env python3
"""
File: add_additional_params.py
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Add missing parameters needed by polymatic for bonding.
Works by reading in the polym.in file, then producing parameters for any
combination of those atoms, or atoms that are bonded to atoms in the polym.in
file, as found by looking through the bonds of the monomer data file.

Warning: this file runs very inefficiently! Use it on a small system, maybe a few 
molecules to generate the required additional parameters, then copy them into the
datafile of the larger system.
"""
from autochem import read_file
import itertools
import sys
import re
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-l", "--lammps", help="LAMMPS datafile containing all parameters of individual molecules"
    )
    parser.add_argument(
        "-f", "--forcefield", help="Forcefield 'dictionary' with all parameters that could be added"
    )
    parser.add_argument(
        "-p",
        "--polym",
        help="polym.in file required by Polymatic. For this code to work, this must include a link command",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Filename of the output file. If not given, the original LAMMPS datafile is overwritten",
    )
    args = parser.parse_args()
    if not args.output:
        args.output = args.lammps
    if len(sys.argv) < 3:
        sys.exit(parser.print_help())
    return args


def linker_atom_types(datafile, polym="polym.in", ff="gaff.ff"):
    atoms_to_consider = []
    linkers = {}
    for line in read_file(polym):
        if "link" in line:
            tmp = line.split()[1:]
            for i in tmp:
                atoms_to_consider += i.split(",")
                linker, after_bond = i.split(",")
                linkers[linker] = after_bond
    atoms_to_consider = list(set(atoms_to_consider))

    atom_types = []
    # find atom types of each atom label,
    # then just add all possible parameters
    for line in read_file(ff):
        if "BONDS" in line:
            break
        line = line.split()
        if len(line) > 6:
            if line[0] in atoms_to_consider:
                atom_types.append(line[1])

    # Add atom types of atoms bonded to the linker atoms.
    # To do this, need to create the original monomer system,
    # and read in bond types
    found = False
    for line in read_file(datafile):
        if "Bond Coeffs" in line:
            found = True
            continue
        if "Atom Coeffs" in line:
            break
        line = line.split()
        if len(line) < 5:
            continue
        if any(x in line[-1] for x in atoms_to_consider):
            atom_types += line[-1].split("-")
    return list(set(atom_types)), linkers


def params_from_initial_datafile(datafile):
    params = {}
    params["atoms"] = []
    params["pairs"] = []
    params["bonds"] = []
    params["angles"] = []
    params["dihs"] = []
    params["imps"] = []

    atoms = False
    pairs = False
    bonds = False
    angles = False
    dihs = False
    imps = False
    for line in read_file(datafile):
        if "Masses" in line:
            atoms = True
            continue
        if "Pair Coeffs" in line:
            atoms = False
            pairs = True
            continue
        if "Bond Coeffs" in line:
            pairs = False
            bonds = True
            continue
        if "Angle Coeffs" in line:
            bonds = False
            angles = True
            continue
        if "Dihedral Coeffs" in line:
            angles = False
            dihs = True
            continue
        # improps not always present
        if re.search("^\s*[A-Z]", line):
            dihs = False
        if "Improper Coeffs" in line:
            imps = True
            continue
        if "Atoms" in line:
            break
        if re.search("^\s*$", line):
            continue
        line = line.split()
        if atoms and re.search("^\s*[0-9]", line[0]):
            params["atoms"].append(line[-1])
        if pairs and re.search("^\s*[0-9]", line[0]):
            params["pairs"].append(line[-1])
        if bonds and re.search("^\s*[0-9]", line[0]):
            params["bonds"].append(line[-1])
        if angles and re.search("^\s*[0-9]", line[0]):
            params["angles"].append(line[-1])
        if dihs and re.search("^\s*[0-9]", line[0]):
            params["dihs"].append(line[-1])
        if imps and re.search("^\s*[0-9]", line[0]):
            params["imps"].append(line[-1])

    return params


def atoms_from_ff(atoms, ff):
    ff_atoms = []
    ff_pairs = []
    found = False
    for line in read_file(ff):
        if "ATOMS" in line:
            found = True
            continue
        if "BONDS" in line:
            break
        if found and not re.search("^\s*$", line):
            line = line.split()
            if line[0] in atoms:
                ff_atoms.append([line[0], line[2]])
                ff_pairs.append([line[0]] + line[5:7])
    return ff_atoms, ff_pairs


def bonds_from_ff(atoms, linkers, ff):
    ff_bonds = []
    found = False
    for line in read_file(ff):
        if "BONDS" in line:
            found = True
        if "ANGLES" in line:
            break
        if found and not re.search("^\s*$", line):
            line = line.split()
            if all(x in atoms for x in line[:2]):
                ff_bonds.append(line[:2] + line[3:5])
    # add on linkers by copying lines that include the type after forming
    # the bond
    additional = []
    for linker, newtype in linkers.items():
        for bond in ff_bonds:
            if newtype in bond[:2]:
                tmp = bond.copy()
                for ind, i in enumerate(tmp):
                    if i == newtype:
                        tmp[ind] = linker
                        additional.append(tmp)
                        break
                        # only want to replace one instance of C3, so we have C3, LC1
    ff_bonds += additional
    ret = []
    for i in ff_bonds:
        if i not in ret:
            ret.append(i)
    return ret


def angles_from_ff(atoms, linkers, ff, params_in_original_datafile):
    ff_angles = []
    found = False
    for line in read_file(ff):
        if "ANGLES" in line:
            found = True
        if "DIHEDRALS" in line:
            break
        if found and not re.search("^\s*$", line):
            line = line.split()
            if all(x in atoms for x in line[:3]):
                ff_angles.append(line[:3] + line[4:6])

    orig_atoms = []
    for ang in params_in_original_datafile:
        for atom in ang.split("-"):
            if atom not in orig_atoms:
                orig_atoms.append(atom)
    all_atoms = []
    for atom in atoms + orig_atoms:
        if atom not in all_atoms:
            all_atoms.append(atom)
    all_combos = list(list(i) for i in itertools.product(all_atoms, repeat=3))
    with_types_that_exist = []
    for combo in all_combos:
        with_types_that_exist.append([linkers[i] if i in linkers else i for i in combo])

    # add on linkers by copying lines that include the type after forming
    # the angle
    angles_with_params = []
    for line in ff_angles:
        for combo, types_that_exist in zip(all_combos, with_types_that_exist):
            if line[:3] == types_that_exist:
                params = line[3:]
                angles_with_params.append(combo + params)
                angles_with_params.append(types_that_exist + params)
                possibility = []
                tmpval = combo
                # find each possible different combination,
                # i.e. for combo = [LC1, LC2, CA] and
                # types_that_exist = [C3, C3, CA],
                # also include LC1, C3, CA and
                #              C3, LC2, CA
                for idx, data in enumerate(zip(combo, types_that_exist)):
                    with_linker, without_linker = data
                    if with_linker != without_linker:
                        # change one value and add it
                        to_add = tmpval.copy()
                        to_add[idx] = without_linker
                        possibility.append(to_add)
                for p in possibility:
                    angles_with_params.append(p + params)

    # prevent duplicates...
    checked = []
    ret = []
    for angle in angles_with_params:
        if angle[:3] not in checked:
            checked.append(angle[:3])
            ret.append(angle)
    return ret

    ff_angles = initial + additional
    ret = []
    for i in ff_angles:
        if i not in ret:
            ret.append(i)
    return ret


def compare_db(ff_dihs, all_combos, with_types_that_exist):
    dihs_with_params = []
    for line in ff_dihs:
        for combo, types_that_exist in zip(all_combos, with_types_that_exist):
            # if all atoms defined, check them against the combos
            if all(val != "X" for val in line[:4]):
                if line[:4] == types_that_exist:
                    # variable number of args
                    stop = 0
                    for idx, val in enumerate(line):
                        if "#" in val:
                            stop = idx
                            break
                    params = line[5:stop]
                    # include combo beforee and after linking
                    dihs_with_params.append(combo + params)
                    dihs_with_params.append(types_that_exist + params)
                    possibility = []
                    tmpval = combo
                    # find each possible different combination,
                    # i.e. for combo = [LC1, LC2, CA, CA] and
                    # types_that_exist = [C3, C3, CA, CA],
                    # also include LC1, C3, CA, CA and
                    #              C3, LC2, CA, CA
                    for idx, data in enumerate(zip(combo, types_that_exist)):
                        with_linker, without_linker = data
                        if with_linker != without_linker:
                            # change one value and add it
                            to_add = tmpval.copy()
                            to_add[idx] = without_linker
                            possibility.append(to_add)
                    for p in possibility:
                        dihs_with_params.append(p + params)
            else:
                # check whether the positions of atoms that aren't X
                # match the item
                types_to_check = {idx: val for idx, val in enumerate(line[:4]) if val != "X"}
                # now check the combo dict
                types = {idx: val for idx, val in enumerate(types_that_exist)}
                # make sure all vals from line in ff are the same as in the combo
                # from itertools
                if all(types[idx] == t for idx, t in types_to_check.items()):
                    # variable number of args
                    stop = 0
                    for idx, val in enumerate(line):
                        if "#" in val:
                            stop = idx
                            break
                    params = line[5:stop]
                    dihs_with_params.append(combo + params)
                    dihs_with_params.append(types_that_exist + params)
                    possibility = []
                    tmpval = combo
                    for idx, data in enumerate(zip(combo, types_that_exist)):
                        with_linker, without_linker = data
                        if with_linker != without_linker:
                            # change one value and add it
                            to_add = tmpval.copy()
                            to_add[idx] = without_linker
                            possibility.append(to_add)
                    for p in possibility:
                        dihs_with_params.append(p + params)
    # some duplicates...
    checked = []
    ret = []
    for dih in dihs_with_params:
        if dih[:4] not in checked:
            checked.append(dih[:4])
            ret.append(dih)
    return ret


def find_combos(atoms, linkers, params_in_original_datafile):
    orig_atoms = []
    for dih in params_in_original_datafile:
        for atom in dih.split("-"):
            if atom not in orig_atoms:
                orig_atoms.append(atom)
    # # the above for loop will give us unnecessary dihedrals, but better to include too many
    # # than not enough
    all_atoms = []
    for atom in atoms + orig_atoms:
        if atom not in all_atoms:
            all_atoms.append(atom)
    # all_combos = list(itertools.combinations_with_replacement(all_atoms, 4))
    # # all_combos = [list(i) for i in all_combos]
    # all_combos = [list(j) for i in all_combos for j in itertools.permutations(i, 4)]
    all_combos = list(list(i) for i in itertools.product(all_atoms, repeat=4))
    with_types_that_exist = []
    for combo in all_combos:
        with_types_that_exist.append([linkers[i] if i in linkers else i for i in combo])
    return all_combos, with_types_that_exist


def dihs_from_ff(atoms, linkers, ff, params_in_original_datafile):
    all_combos, with_types_that_exist = find_combos(atoms, linkers, params_in_original_datafile)

    ff_dihs = []
    found = False
    for line in read_file(ff):
        if "DIHEDRALS" in line:
            found = True
            continue
        if "IMPROPERS" in line:
            break
        if found and not re.search("^\s*$", line):
            ff_dihs.append(line.split())

    return compare_db(ff_dihs, all_combos, with_types_that_exist)


def imps_from_ff(atoms, linkers, ff, params_in_original_datafile):
    """
    Same as dihs from ff, just changing the values in the force field 
    to search for
    """
    all_combos, with_types_that_exist = find_combos(atoms, linkers, params_in_original_datafile)

    ff_dihs = []

    found = False
    for line in read_file(ff):
        if "IMPROPERS" in line:
            found = True
            continue
        if found and not re.search("^\s*$", line):
            ff_dihs.append(line.split())

    return compare_db(ff_dihs, all_combos, with_types_that_exist)


def add_atoms(datafile, atoms, original):
    df = datafile.copy()
    atoms_to_add = [i for i in atoms if i[0] not in original]
    found = False
    num = 0
    total = len(original)
    first_newline_skipped = False
    # loop over original file and then add to the copy in the correct place
    for idx, i in enumerate(datafile):
        if "Masses" in i:
            found = True
            continue
        if "Pair Coeffs" in i:
            break
        if found:
            if i != "\n":
                num += 1
            if i == "\n" and not first_newline_skipped:
                first_newline_skipped = True
                continue
            if i == "\n" and first_newline_skipped:
                for new, val in enumerate(atoms_to_add, 1):
                    total += 1
                    df.insert(
                        idx - 1 + new,
                        "{:<4} {:>9.4f}    # {}\n".format(total, float(val[1]), val[0]),
                    )
    return df, total


def add_pairs(datafile, pairs, original):
    df = datafile.copy()
    pairs_to_add = [i for i in pairs if i[0] not in original]
    found = False
    num = 0
    total = len(original)
    first_newline_skipped = False
    # loop over original file and then add to the copy in the correct place
    for idx, i in enumerate(datafile):
        if "Pair Coeffs" in i:
            found = True
            continue
        if "Bond Coeffs" in i:
            break
        if found:
            if i != "\n":
                num += 1
            if i == "\n" and not first_newline_skipped:
                first_newline_skipped = True
                continue
            if i == "\n" and first_newline_skipped:
                for new, val in enumerate(pairs_to_add, 1):
                    total += 1
                    df.insert(
                        idx - 1 + new,
                        "{:<4} {:>10.3f} {:>9.5f}    # {}\n".format(
                            total, float(val[1]), float(val[2]), val[0]
                        ),
                    )
    return df, total


def add_bonds(datafile, bonds, original):
    df = datafile.copy()
    bonds_to_add = [i for i in bonds if f"{i[0]}-{i[1]}" not in original]
    found = False
    num = 0
    total = len(original)
    first_newline_skipped = False
    # loop over original file and then add to the copy in the correct place
    for idx, i in enumerate(datafile):
        if "Bond Coeffs" in i:
            found = True
            continue
        if "Angle Coeffs" in i:
            break
        if found:
            if i != "\n":
                num += 1
            if i == "\n" and not first_newline_skipped:
                first_newline_skipped = True
                continue
            if i == "\n" and first_newline_skipped:
                for new, val in enumerate(bonds_to_add, 1):
                    total += 1
                    df.insert(
                        idx - 1 + new,
                        "{:<4} {:>9.1f} {:>10.3f}    # {}-{}\n".format(
                            total, float(val[2]), float(val[3]), val[0], val[1]
                        ),
                    )
    return df, total


def add_angles(datafile, angles, original):
    df = datafile.copy()
    angles_to_add = [i for i in angles if f"{i[0]}-{i[1]}-{i[2]}" not in original]
    found = False
    num = 0
    total = len(original)
    first_newline_skipped = False
    # loop over original file and then add to the copy in the correct place
    for idx, i in enumerate(datafile):
        if "Angle Coeffs" in i:
            found = True
            continue
        if "Dihedral Coeffs" in i:
            break
        if found:
            if i != "\n":
                num += 1
            if i == "\n" and not first_newline_skipped:
                first_newline_skipped = True
                continue
            if i == "\n" and first_newline_skipped:
                for new, val in enumerate(angles_to_add, 1):
                    total += 1
                    df.insert(
                        idx - 1 + new,
                        "{:<4} {:>9.1f} {:>10.3f}    # {}-{}-{}\n".format(
                            total, float(val[3]), float(val[4]), val[0], val[1], val[2]
                        ),
                    )
    return df, total


def add_dihs(datafile, dihs, original):
    df = datafile.copy()
    dihs_to_add = [i for i in dihs if f"{i[0]}-{i[1]}-{i[2]}-{i[3]}" not in original]
    found = False
    num = 0
    total = len(original)
    first_newline_skipped = False
    # loop over original file and then add to the copy in the correct place
    for idx, i in enumerate(datafile):
        if "Dihedral Coeffs" in i:
            found = True
            continue
        if "Improper Coeffs" in i or "Atoms" in i:
            break
        if found:
            if i != "\n":
                num += 1
            if i == "\n" and not first_newline_skipped:
                first_newline_skipped = True
                continue
            if i == "\n" and first_newline_skipped:
                for new, val in enumerate(dihs_to_add, 1):
                    total += 1
                    params = " ".join(val[4:])
                    df.insert(
                        idx - 1 + new,
                        "{:<4} {}  # {}-{}-{}-{}\n".format(
                            total, params, val[0], val[1], val[2], val[3]
                        ),
                    )
    return df, total


def add_imps(datafile, imps, original):
    df = datafile.copy()
    imps_to_add = [i for i in imps if f"{i[0]}-{i[1]}-{i[2]}-{i[3]}" not in original]
    found = False
    num = 0
    total = len(original)
    first_newline_skipped = False
    # loop over original file and then add to the copy in the correct place
    for idx, i in enumerate(datafile):
        if "Improper Coeffs" in i:
            found = True
            continue
        if "Atoms" in i:
            break
        if found:
            if i != "\n":
                num += 1
            if i == "\n" and not first_newline_skipped:
                first_newline_skipped = True
                continue
            if i == "\n" and first_newline_skipped:
                for new, val in enumerate(imps_to_add, 1):
                    total += 1
                    params = " ".join(val[4:])
                    df.insert(
                        idx - 1 + new,
                        "{:<4} {}  # {}-{}-{}-{}\n".format(
                            total, params, val[0], val[1], val[2], val[3]
                        ),
                    )
    return df, total


def main():
    args = get_args()
    atom_types_to_add, linkers = linker_atom_types(
        args.lammps, ff=args.forcefield, polym=args.polym
    )
    params = params_from_initial_datafile(args.lammps)
    addn_atoms, addn_pairs = atoms_from_ff(atom_types_to_add, args.forcefield)
    addn_bonds = bonds_from_ff(atom_types_to_add, linkers, args.forcefield)
    addn_angles = angles_from_ff(atom_types_to_add, linkers, args.forcefield, params["angles"])
    addn_dihs = dihs_from_ff(atom_types_to_add, linkers, args.forcefield, params["dihs"])
    addn_imps = imps_from_ff(atom_types_to_add, linkers, args.forcefield, params["imps"])

    with open(args.lammps) as f:
        lmps = f.readlines()

    lmps, num_atoms = add_atoms(lmps, addn_atoms, params["atoms"])
    lmps, num_pairs = add_pairs(lmps, addn_pairs, params["pairs"])
    lmps, num_bonds = add_bonds(lmps, addn_bonds, params["bonds"])
    lmps, num_angles = add_angles(lmps, addn_angles, params["angles"])

    # number of dihs and imps added changes each time...
    lmps, num_dihs = add_dihs(lmps, addn_dihs, params["dihs"])
    lmps, num_imps = add_imps(lmps, addn_imps, params["imps"])

    # now change the numbers at top of file
    new = lmps.copy()
    for idx, line in enumerate(lmps):
        if "atom types" in line:
            new[idx] = f" {num_atoms} atom types\n"
        if "bond types" in line:
            new[idx] = f" {num_bonds} bond types\n"
        if "angle types" in line:
            new[idx] = f" {num_angles} angle types\n"
        if "dihedral types" in line:
            new[idx] = f" {num_dihs} dihedral types\n"
        if "improper types" in line:
            new[idx] = f" {num_imps} improper types\n"

    with open(args.output, "w") as n:
        n.writelines(new)


main()
