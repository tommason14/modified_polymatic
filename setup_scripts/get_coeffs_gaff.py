import re, sys, itertools

kj2kcal = 4.1868
factor = 2

def getMass(Atom, path_to_ff):
    with open(path_to_ff, 'r') as f:
        for line in f:
            if re.search("^" + Atom + "\s", line):
                line = line.split()
                return float(line[2])
                
            # IF GET TO BOND SECTION
            elif re.search("^\s*BONDS", line):
                sys.exit("Could not find atom {}".format(Atom))

# GROUPS FOR BONDS, ANGLES, DIHEDRALS, IMPROPERS
def getGroup(Atom, path_to_ff):

    # GET ATOM GROUP FROM ff
    with open(path_to_ff, "r+") as f:
        for line in f:

            # IF LINE STARTS WITH Atom NAME
            if re.search("^" + Atom + "\s", line):
                group = line.split()[1]
                return group

            # IF GET TO BOND SECTION
            elif re.search("^\s*BONDS", line):
                sys.exit("Could not find atom {}".format(Atom))


def getAtomData(Atom, path_to_ff):
    # GET ATOM DATA FROM ff
    with open(path_to_ff, "r+") as f:
        for line in f:

            # IF LINE STARTS WITH Atom NAME
            if re.search("^" + Atom + "\s", line):
                line = line.split()
                eps, sigma = float(line[5]), float(line[6])
                return eps, sigma

            # IF GET TO BOND SECTION
            elif re.search("^\s*BONDS", line):
                sys.exit("Could not find atom {}".format(Atom))


def getAtomPartialCharge(Atom, path_to_ff):

    # GET ATOM DATA FROM ff
    with open(path_to_ff, "r+") as f:
        for line in f:

            # IF LINE STARTS WITH Atom NAME
            if re.search("^" + Atom + "\s", line):
                line = line.split()
                q = float(line[3])
                return q

            # IF GET TO BOND SECTION
            elif re.search("^\s*BONDS", line):
                sys.exit("Could not find atom {}".format(Atom))


def getBond(myAtom1, myAtom2, path_to_ff):

    # GET FF ATOM NAME
    Atom1 = getGroup(myAtom1, path_to_ff)
    Atom2 = getGroup(myAtom2, path_to_ff)

    # GET BOND DATA FROM ff
    found = False
    with open(path_to_ff, "r+") as f:
        for line in f:

            # IF FOUND BOND SECTION
            if re.search("BONDS", line):
                found = True

            # IF LINE STARTS WITH Atom1 Atom2, take params as they are in gaff.ff
            elif found and re.search("^" + Atom1 + "\s*" + Atom2 + "\s", line):
                line = line.split()
                kr, Re = float(line[3]), float(line[4])
                return kr, Re

            # IF LINE STARTS WITH Atom2 Atom1
            elif found and re.search("^" + Atom2 + "\s*" + Atom1 + "\s", line):
                line = line.split()
                kr, Re = float(line[3]), float(line[4])
                return kr, Re

            # IF GET TO ANGLE SECTION
            elif re.search("^\s*ANGLES", line):
                print("Could not find bond {} {}".format(Atom1, Atom2))
                return 0.0, 0.0

def getAngle(myAtom1, myAtom2, myAtom3, path_to_ff):

    # GET FF ATOM NAME
    Atom1 = getGroup(myAtom1, path_to_ff)
    Atom2 = getGroup(myAtom2, path_to_ff)
    Atom3 = getGroup(myAtom3, path_to_ff)

    # GET BOND DATA FROM ff
    found = False
    with open(path_to_ff, "r+") as f:
        for line in f:

            # IF FOUND ANGLES SECTION
            if re.search("^\s*ANGLES", line):
                found = True

            # IF LINE STARTS WITH Atom1 Atom2 Atom3
            elif found and re.search(
                "^" + Atom1 + "\s*" + Atom2 + "\s*" + Atom3 + "\s", line
            ):
                line = line.split()
                ka = float(line[4])
                th = float(line[5])
                return ka, th

            # IF LINE STARTS WITH Atom3 Atom2 Atom1
            elif found and re.search(
                "^" + Atom3 + "\s*" + Atom2 + "\s*" + Atom1 + "\s", line
            ):
                line = line.split()
                ka = float(line[4])
                th = float(line[5])
                return ka, th

            # IF GET TO DIHEDRALS SECTION
            elif re.search("^\s*DIHEDRALS", line):
                print("Could not find angle {} {} {}".format(Atom1, Atom2, Atom3))
                return 0.00, 0.00

def getDihedral(myAtom1, myAtom2, myAtom3, myAtom4, path_to_ff):

    # GET FF ATOM NAME
    Atom1 = getGroup(myAtom1, path_to_ff)
    Atom2 = getGroup(myAtom2, path_to_ff)
    Atom3 = getGroup(myAtom3, path_to_ff)
    Atom4 = getGroup(myAtom4, path_to_ff)

    # GET BOND DATA FROM ff
    # just check middle 2 atoms
    found = False
    with open(path_to_ff, "r+") as f:
        for line in f:

            # IF FOUND DIHEDRALS SECTION
            if re.search("^\s*DIHEDRALS", line):
                found = True

            # IF LINE STARTS WITH Atom1 Atom2 Atom3 Atom4
            elif found and re.search("^.*\s*" + Atom2 + "\s*" + Atom3, line):
                line = line.split()
                # variable number of args
                for idx, val in enumerate(line):
                    if "#" in val:
                        stop = idx
                return line[5:stop]

            # IF LINE STARTS WITH Atom4 Atom3 Atom2 Atom1
            elif found and re.search("^.*\s*" + Atom3 + "\s*" + Atom2, line):
                line = line.split()
                # variable number of args
                for idx, val in enumerate(line):
                    if "#" in val:
                        stop = idx
                return line[5:stop]

            # IF GET TO DIHEDRALS SECTION
            elif re.search("^\s*IMPROPER", line):
                print(
                    "Could not find dihedral {} {} {} {}".format(
                        Atom1, Atom2, Atom3, Atom4
                    )
                )
                return ['1', '0.0', '0', '0.0']

def getImproper(myAtom1, myAtom2, myAtom3, myAtom4, path_to_ff):

    # GET FF ATOM NAME
    Atom1 = getGroup(myAtom1, path_to_ff)
    Atom2 = getGroup(myAtom2, path_to_ff)
    Atom3 = getGroup(myAtom3, path_to_ff)
    Atom4 = getGroup(myAtom4, path_to_ff)

    atoms = [a for a in [Atom1, Atom2, Atom3, Atom4] if a != "X"]

    perms = list(itertools.permutations(atoms, len(atoms)))
    # GET BOND DATA FROM ff
    found = False
    with open(path_to_ff, "r+") as f:
        for line in f:
            # IF FOUND DIHEDRALS SECTION
            if re.search("^\s*IMPROPER", line):
                found = True
            elif found:
                for perm in perms:
                    a1, a2, a3, a4 = perm
                    # now we have CA CA CA HA, but any of these atoms could be an X
                    # want a partial match of any two atoms
                    if (
                        re.search(f"{a1}.*{a2}", line)
                        or re.search(f"{a2}.*{a3}", line)
                        or re.search(f"{a3}.*{a4}", line)
                    ):
                        return line.split()[5:8]

    print("Could not find improper {} {} {} {}".format(Atom1, Atom2, Atom3, Atom4))
    return ['0.0', '0', '0']
