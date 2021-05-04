#!/usr/bin/perl

################################################################################
#
# file_convert.pl
# This file is part of the Polymatic distribution.
#
# Author: Lauren J. Abbott
# Version: 1.1
# Date: August 16, 2015
#
# Description: File conversion using the Polymatic.pm module. The molecular
# system can be read in from multiple file types. Connectivity definition can
# be updated if bonds are given. The molecular system then can be written to
# multiple file types. The Polymatic.pm module is used in this code. It must be
# in the same directory or at a file path recognized by the script (e.g., with
# 'use lib').
#
# Syntax:
#  ./file_convert.pl <flags>
#
# File input flags:
#  --ilmp file, read LAMMPS data file named 'file'
#  --ipdb file, read PDB file named 'file'
#  --igro file, read GRO file named 'file'
#  --ixyz file, read XYZ file named 'file'
#  --ixsd file, read XSD file named 'file'
#  --ipsf file, read PSF file named 'file'
#  --ityp file, read types file named 'file'
#
# Connectivity flags:
#  --mol, define all molecules
#  --ang, define all angles
#  --dih, define all dihedrals
#  --imp, define all impropers
#
# File output flags:
#  --olmp file, write LAMMPS data file named 'file'
#  --otrj file, write LAMMPS dump file named 'file'
#  --opdb file, write PDB file named 'file'
#  --ogro file, write GRO file named 'file'
#  --oxyz file, write XYZ file named 'file'
#  --opsf file, write PSF file named 'file'
#  --ocar file, write CAR file named 'file'
#  --omdf file, write MDF file named 'file'
#  --otyp file, write types file named 'file'
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

use strict;
use File::Basename;
use lib dirname (__FILE__);
use Polymatic();
use Getopt::Long();

# Variables
my %opt;
my %sys;

# Command-line options
Getopt::Long::GetOptions(\%opt,
    'ilmp=s', 'ipdb=s', 'igro=s', 'ixyz=s', 'ipsf=s', 'ityp=s', 'ixsd=s',
    'olmp=s', 'opdb=s', 'ogro=s', 'oxyz=s', 'opsf=s', 'otyp=s', 'otrj=s',
    'ocar=s', 'omdf=s', 'ang', 'dih', 'imp', 'mol'
    );

# Input LAMMPS data file
if ($opt{'ilmp'}) {
    printf "Reading LAMMPS data file '%s'\n", $opt{'ilmp'};
    %sys = Polymatic::readLammps($opt{'ilmp'});
}

# Input PDB file
if ($opt{'ipdb'}) {
    printf "Reading PDB file '%s'\n", $opt{'ipdb'};
    %sys = Polymatic::readPdb($opt{'ipdb'});
}

# Input GRO file
if ($opt{'igro'}) {
    printf "Reading GRO file '%s'\n", $opt{'igro'};
    %sys = Polymatic::readGro($opt{'igro'});
}

# Input XYZ file
if ($opt{'ixyz'}) {
    printf "Reading XYZ file '%s'\n", $opt{'ixyz'};
    %sys = Polymatic::readXyz($opt{'ixyz'});
}

# Input XSD file
if ($opt{'ixsd'}) {
    printf "Reading XSD file '%s'\n", $opt{'ixsd'};
    %sys = Polymatic::readXsd($opt{'ixsd'});
}

# Input PSF file
if ($opt{'ipsf'}) {
    printf "Reading PSF file '%s'\n", $opt{'ipsf'};
    Polymatic::readPsf($opt{'ipsf'}, \%sys);
}

# Input types file
if ($opt{'ityp'}) {
    printf "Reading types file '%s'\n", $opt{'ityp'};
    Polymatic::readTypes($opt{'ityp'}, \%sys);
}

# Define angles
if ($opt{'ang'}) {
    printf "Defining all angles in system\n";
    Polymatic::defineAngles(\%sys);
}

# Define dihedrals
if ($opt{'dih'}) {
    printf "Defining all dihedrals in system\n";
    Polymatic::defineDiheds(\%sys);
}

# Define impropers
if ($opt{'imp'}) {
    printf "Defining all impropers in system\n";
    Polymatic::defineImprops(\%sys);
}

# Define molecules
if ($opt{'mol'}) {
    printf "Defining all molecules in system\n";
    Polymatic::defineMols(\%sys);
}

# Output LAMMPS data file
if ($opt{'olmp'}) {
    printf "Writing LAMMPS data file '%s'\n", $opt{'olmp'};
    Polymatic::writeLammps($opt{'olmp'}, \%sys);
}

# Output PDB file
if ($opt{'opdb'}) {
    printf "Writing PDB file '%s'\n", $opt{'opdb'};
    Polymatic::writePdb($opt{'opdb'}, \%sys);
}

# Output GRO file
if ($opt{'ogro'}) {
    printf "Writing GRO file '%s'\n", $opt{'ogro'};
    Polymatic::writeGro($opt{'ogro'}, \%sys);
}

# Output XYZ file
if ($opt{'oxyz'}) {
    printf "Writing XYZ file '%s'\n", $opt{'oxyz'};
    Polymatic::writeXyz($opt{'oxyz'}, \%sys);
}

# Output PSF file
if ($opt{'opsf'}) {
    printf "Writing PSF file '%s'\n", $opt{'opsf'};
    Polymatic::writePsf($opt{'opsf'}, \%sys);
}

# Output types file
if ($opt{'otyp'}) {
    printf "Writing types file '%s'\n", $opt{'otyp'};
    Polymatic::writeTypes($opt{'otyp'}, \%sys);
}

# Output LAMMPS dump file
if ($opt{'otrj'}) {
    printf "Writing LAMMPS dump file '%s'\n", $opt{'otrj'};
    Polymatic::writeLmpsTrj($opt{'otrj'}, \%sys);
}

# Output CAR file
if ($opt{'ocar'}) {
    printf "Writing CAR file '%s'\n", $opt{'ocar'};
    Polymatic::writeCar($opt{'ocar'}, \%sys);
}

# Output MDF file
if ($opt{'omdf'}) {
    printf "Writing MDF file '%s'\n", $opt{'omdf'};
    Polymatic::writeMdf($opt{'omdf'}, \%sys);
}
