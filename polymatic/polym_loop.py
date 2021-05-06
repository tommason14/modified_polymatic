#!/usr/bin/env python3

################################################################################
#
# polym_loop.py
# This file is part of the Polymatic distribution.
#
# Author: Lauren J. Abbott
# Version: 1.1
# Date: August 16, 2015
#
# Description: Controls the simulated polymerization loop of the Polymatic
# algorithm. Polymerization steps are performed in cycles. After each bond is
# made, an energy minimization is performed in LAMMPS. After each cycle, a
# molecular dynamics step is performed in LAMMPS.
#
# Syntax:
#  ./polym_loop.py
#
# User parameters and file paths necessary for the polymerization should be
# specified at the beginning of the script.
#
################################################################################
#
# Polymatic: a general simulated polymerization algorithm
# Copyright (C) 2013, 2015 Lauren J. Abbott
#
# Polymatic is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Polymatic is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License (COPYING)
# along with Polymatic. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import sys
import os
import subprocess
import shutil
import glob
import argparse

#
# User Setup
#

parser = argparse.ArgumentParser()
parser.add_argument(
    "--no-minimise",
    help=
    "By default, minimsisation of pack.lmps is performed. If you already have a data.lmps file that "
    "you do not want to minimise, pass in this flag",
    action="store_true",
)
parser.add_argument(
    "--no-shuffle",
    help=
    "By default, an quick NVT MD simulation is performed after all of one monomer has been consumed, as early "
    "tests showed this to be necessary for polymerisation to progress. Turn off this behaviour by "
    "using this flag",
    action="store_true",
)
parser.add_argument(
    "--controlled",
    help=
    "Pass in this flag to used a modified version of the Polmatic code that accounts for free "
    "energy barriers of each possible reaction, stored in a barriers.in file. A file containing "
    "molecule IDs is also required, named moleculeID.txt",
    action="store_true",
)
parser.add_argument("-k",
                    "--keep",
                    help="Retain folders of each intermediate step",
                    action="store_true")
args = parser.parse_args()

# Parameters
bonds = 0
bonds_cyc = 5
md_cyc = 3
md_max = 20

# Input Scripts
input_polym = "polym.in"
input_barriers = "barriers.in"
init_min = "init.in"
input_min = "min.in"
input_md0 = "md0.in"
input_md1 = "md1.in"
input_md2 = "md2.in"

# Polymatic Scripts
script_step = "polym.pl"
script_step_no_mods = "polym_no_mods.pl"
script_init = "polym_init.pl"
script_init_no_mods = "polym_init_no_mods.pl"
script_final = "polym_final.pl"

cwd = os.path.dirname(os.path.realpath(__file__))

lmps = os.environ["LAMMPS_EXEC"]

# Main
#


def main():
    check_files_in_dir()
    if not args.no_minimise:
        minimisation_of_packmol()
    print_header()
    polym_loop()
    print_footer()


#
# Functions
#


def polym_loop():

    # Variables
    global bonds

    # Initialization
    polym_init()

    # Step 0
    setup_step()
    em()

    write_vmd_frag_file()

    shuffled = False
    # Steps 1-N
    frags = frags_from_data_lmps()
    while frags > 1:
        if args.controlled:
            # should the atoms be shuffled?
            if not shuffled and shuffle_required() and not args.no_shuffle:
                setup_shuffle()
                run_shuffle()
                shuffled = True

        # Polymerization step
        bonds += 1
        setup_step()
        code = polym_step()
        if code == 1:
            break
        em()

        # # Stop if finished
        # if (bonds == bonds_tot):
        #     break

        # Molecular dynamics
        if bonds % bonds_cyc == 0:
            setup_md(0)
            if (bonds / bonds_cyc) % md_cyc == 0:
                md(2)
            else:
                md(1)
        frags = frags_from_vmd()
        print(f"  Fragments present: {frags}")

    # Finalization
    polym_final()


def polym_step():

    # Variables
    global bonds
    attempts = 1

    # Attempt until successful or max attempts
    while 1:

        # Polymerization step
        if args.controlled:
            cmd = f"perl {cwd}/{script_step} -i init.lmps -t ../types.txt -s ../{input_polym} -b ../barriers.in -n ../moleculeID.txt -m ../bonded.tmp -o data.lmps"
        else:
            cmd = f"perl {cwd}/{script_step_no_mods} -i init.lmps -t ../types.txt -s ../{input_polym} -o data.lmps"
        sys.stdout.flush()
        code = subprocess.call(cmd, shell=True)

        # Bond made
        if code == 0:
            print(f"  Attempts: {attempts}")
            if not os.path.isfile("data.lmps"):
                err_exit(
                    "Polymerization step script did not complete properly.")
            return 0

        # No pair found
        elif code == 3:

            # Stop if maximum attempts reached
            if attempts > md_max:
                print(
                    "  No pair was found within the maximum number of attempts."
                )
                bonds -= 1
                os.chdir("..")
                if not args.keep:
                    for path in glob.glob("step_*"):
                        shutil.rmtree(path)
                return 1

            # Molecular dynamics
            setup_md(attempts)
            md(0)
            attempts += 1

        # Error
        else:
            err_exit("Polymerization step script did not complete properly.")


def minimisation_of_packmol():
    """
    Run initial minimisation from the top dir
    """
    print("Minimising packmol stucture...", end=" ")
    subprocess.call("mkdir init_min".split())
    cmd = f"{lmps} -i scripts/{init_min}"
    subprocess.call(
        cmd,
        shell=True,
        stdout=open("init_min/minimisation_of_initial_structure.out", "w"),
    )
    subprocess.call("mv log.lammps init_min", shell=True)
    print("complete")


def polym_init():
    print("Initialization:")
    if script_init != 0:
        if args.controlled:
            cmd = f"perl {cwd}/{script_init} -i data.lmps -t types.txt -s {input_polym} -o temp.lmps -b {input_barriers}"
        else:
            cmd = f"perl {cwd}/{script_init_no_mods} -i data.lmps -t types.txt -s {input_polym} -o temp.lmps"
        sys.stdout.flush()
        code = subprocess.call(cmd, shell=True)
    else:
        shutil.copy("data.lmps", "temp.lmps")
    if code != 0 or not os.path.isfile("temp.lmps"):
        # err_exit(
        # "Polymerization initialization script did not complete properly.")
        num_attempts = 0
        successful = False
        while not successful:
            num_attempts += 1
            if num_attempts == 5:
                # limit number of restarts- if simulation reaches here, bad starting trajectory
                sys.exit(
                    "Incorrect structure- try changing the initial box dimensions"
                )
            print(
                "Polymerization initialization script did not complete properly."
            )
            print(
                "Minimising trajectory and then attempting initialisation step again"
            )
            os.mkdir("restart")
            shutil.copy("data.lmps", "restart/data.lmps")
            shutil.copy("../scripts/restart_min.in", "restart/min.in")
            os.chdir("restart")
            cmd = f"{lmps} -i min.in"
            subprocess.call(cmd, shell=True, stdout=open("out", "w"))
            # replace 'old' data.lmps with the minimised trajectory
            os.system("rm ../data.lmps")
            shutil.copy("md.lmps", "../data.lmps")
            os.chdir("..")
            shutil.rmtree("restart")
            # try again
            cmd = f"perl {cwd}/{script_init} -i data.lmps -t types.txt -s {input_polym} -o temp.lmps -b {input_barriers}"
            sys.stdout.flush()
            code = subprocess.call(cmd, shell=True)
            successful = os.path.isfile(
                "md.lmps")  # if not, then it will keep minimising


def polym_final():
    print("Finalization:")
    if script_final != 0:
        cmd = f"perl {cwd}/{script_final} -i temp.lmps -t types.txt -s {input_polym} -o final.lmps"
        sys.stdout.flush()
        code = subprocess.call(cmd, shell=True)
    else:
        shutil.copy("temp.lmps", "final.lmps")
    if code != 0 or not os.path.isfile("final.lmps"):
        err_exit(
            "Polymerization finalization script did not complete properly.")
    os.remove("temp.lmps")


def em():

    # Input script
    shutil.copy("../scripts/" + input_min, input_min)

    # LAMMPS EM
    cmd = f"{lmps} -i {input_min}"
    code = subprocess.call(cmd, shell=True, stdout=open("out", "w"))
    if not os.path.isfile("min.lmps"):
        err_exit("LAMMPS energy minimization did not complete properly.")

    # Data file
    shutil.copy("min.lmps", "../temp.lmps")
    os.chdir("..")

    # Keep files?
    if not args.keep:
        for path in glob.glob("step_*"):
            shutil.rmtree(path)


def md(num):

    # Input script
    if num == 0:
        inp = input_md0
        shutil.copy("../../scripts/" + inp, inp)
    elif num == 1:
        inp = input_md1
        shutil.copy("../scripts/" + inp, inp)
    elif num == 2:
        inp = input_md2
        shutil.copy("../scripts/" + inp, inp)

    # LAMMPS MD
    cmd = f"{lmps} -i {inp}"
    code = subprocess.call(cmd, shell=True, stdout=open("out", "w"))
    if not os.path.isfile("md.lmps"):
        # err_exit("LAMMPS molecular dynamics did not complete properly.")
        num_attempts = 0
        successful = False
        while not successful:
            num_attempts += 1
            if num_attempts == 5:
                # limit number of restarts- if simulation reaches here, bad starting trajectory
                sys.exit(
                    "Incorrect structure- try changing the initial box dimensions"
                )
            print("LAMMPS molecular dynamics did not complete properly.")
            print("Minimising trajectory and then attempting md step again")
            os.mkdir("restart")
            shutil.copy("data.lmps", "restart/data.lmps")
            shutil.copy("../scripts/restart_min.in", "restart/min.in")
            os.chdir("restart")
            cmd = f"{lmps} -i min.in"
            subprocess.call(cmd, shell=True, stdout=open("out", "w"))
            # replace 'old' data.lmps with the minimised trajectory
            os.system("rm ../data.lmps")
            shutil.copy("md.lmps", "../data.lmps")
            os.chdir("..")
            shutil.rmtree("restart")
            # try again
            cmd = f"{lmps} -i {inp}"
            code = subprocess.call(cmd, shell=True, stdout=open("out", "w"))
            successful = os.path.isfile(
                "md.lmps")  # if not, then it will keep minimising

    # Data file
    if num == 0:
        shutil.copy("md.lmps", "../init.lmps")
    else:
        shutil.copy("md.lmps", "../temp.lmps")
    os.chdir("..")

    # Keep files?
    if not args.keep:
        if num == 0:
            for path in glob.glob("md_*"):
                shutil.rmtree(path)
        else:
            for path in glob.glob("step_*"):
                shutil.rmtree(path)


def setup_step():

    # Directory
    print(f"Step {bonds}")
    directory = "step_" + "{0:03d}".format(bonds)
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        err_exit(f"Directory '{directory}' already exists.")
    os.chdir(directory)

    # Data file
    if bonds == 0:
        shutil.copy("../temp.lmps", "data.lmps")
    else:
        shutil.copy("../temp.lmps", "init.lmps")


def setup_md(num):

    # Directory
    if num == 0:
        directory = "step_" + "{0:03d}".format(bonds) + "_md"
    else:
        directory = "md_" + "{0:03d}".format(num)
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        err_exit(f"Directory '{directory}' already exists.")
    os.chdir(directory)

    # Data file
    if num == 0:
        shutil.copy("../temp.lmps", "data.lmps")
    else:
        shutil.copy("../init.lmps", "data.lmps")


def print_header():
    print("Polymatic Simulated Polymerization")
    print("Parameters\n----------")
    print(f"Bonds per cycle:           {bonds_cyc}")
    print(f"Frequency of MD type 2:    {md_cyc}")
    print(f"Maximum bond attempts:     {md_max}")
    print("Polymerization Loop\n-------------------")


def print_footer():
    print("\nSummary\n-------")
    print(f"Bonds made:                {bonds}")
    if args.controlled:
        os.system("rm bonded.tmp read_data.tcl")
    else:
        os.system("rm read_data.tcl")


def err_exit(error):
    print(f"Error: {error}")
    sys.exit(1)


def shuffle_required(filename="bonded.tmp"):
    """
    Returns True if all of one monomer has been consumed in the polymerisation.
    Giving the atoms a 'shuffle' by running NVT at high temperature might help
    bring the atoms into closer proximity and form a bond.
    """
    with open(filename) as f:
        for line in f:
            if "total mol1" in line:
                tot_1 = float(line.split()[-1])
            if "total mol2" in line:
                tot_2 = float(line.split()[-1])
            if "mol1 included" in line:
                mol_1 = float(line.split()[-1])
            if "mol2 included" in line:
                mol_2 = float(line.split()[-1])
    return mol_1 >= tot_1 or mol_2 >= tot_2


def setup_shuffle():

    # Directory
    directory = "md_shuffle"
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        err_exit(f"Directory '{directory}' already exists.")
    os.chdir(directory)

    # Data file
    shutil.copy("../temp.lmps", "data.lmps")
    shutil.copy("../scripts/shuffle.in", "shuffle.in")


def run_shuffle():

    print("Shuffling atoms")
    cmd = f"{lmps} -i shuffle.in"
    code = subprocess.call(cmd, shell=True, stdout=open("out", "w"))

    # Data file
    shutil.copy("md.lmps", "../temp.lmps")
    os.chdir("..")

    # Keep files?
    if not args.keep:
        shutil.rmtree("md_shuffle")


def write_vmd_frag_file():
    cmds = [
        "package require topotools\n", "topo readlammpsdata temp.lmps\n",
        "exit"
    ]
    with open("read_data.tcl", "w") as f:
        f.writelines(cmds)


def frags_from_vmd():
    cmd = "vmd -dispdev text -e read_data.tcl 2> /dev/null | grep Fragments | awk '{print $3}'"
    return int(
        subprocess.check_output(cmd,
                                shell=True).decode("utf-8").split("\n")[0])


def frags_from_data_lmps():
    """
    Only happens on first step
    """
    cmds = [
        "package require topotools\n", "topo readlammpsdata data.lmps\n",
        "exit"
    ]
    with open("tmp.tcl", "w") as f:
        f.writelines(cmds)

    cmd = (
        "vmd -dispdev text -e tmp.tcl 2> /dev/null | grep Fragments | awk '{print $3}'"
    )
    ret = int(
        subprocess.check_output(cmd,
                                shell=True).decode("utf-8").split("\n")[0])
    os.system("rm tmp.tcl")
    return ret


def check_files_in_dir():
    if args.no_minimise:
        required = ("data.lmps", "polym.in", "types.txt")
    else:
        required = ("pack.lmps", "polym.in", "types.txt")
    no_file = [f for f in required if not os.path.isfile(f)]
    if len(no_file) > 0:
        print("Required file(s) not present:")
        for f in no_file:
            print(f"- {f}")
        sys.exit()


if __name__ == "__main__":
    main()
