#!/usr/bin/perl

################################################################################
#
# polym_final.pl
# This file is part of the Polymatic distribution.
#
# Author: Lauren J. Abbott
# Version: 1.1
# Date: August 16, 2013
#
# Description: Performs a polymerization finalization for use with the Polymatic
# code. Finds all linker atoms as defined by the 'link' command in the input
# script and subtracts the appropriate artificial charges as given by the
# 'charge' command in the input script. Reads and writes LAMMPS data files. The
# Polymatic.pm module is used in this code. It must be in the same directory or
# at a file path recognized by the script (e.g., with 'use lib').
#
# Syntax:
#  ./polym.pl -i data.lmps
#             -t types.txt
#             -s polym.in
#             -o new.lmps
#
# Arguments:
#  1. data.lmps, LAMMPS data file of initial system (-i)
#  2. types.txt, data types text file (-t)
#  3. polym.in, input script specifying polymerization options (-l)
#  4. new.lmps, updated LAMMPS data file after polymerization step (-o)
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
use Math::Trig();
use POSIX();

# Variables
my ($fileData, $fileTypes, $fileInput, $fileOut);
my (%sys, %inp);

# Main
readArgs();
%sys = Polymatic::readLammps($fileData);
Polymatic::readTypes($fileTypes, \%sys);
%inp = Polymatic::readPolym($fileInput);
removeCharges();
unwrap();
Polymatic::writeLammps($fileOut, \%sys);

################################################################################
# Subroutines

# errExit( $error )
# Exit program and print error
sub errExit
{
    printf "Error: %s\n", $_[0];
    exit 2;
}

# readArgs( )
# Read command line arguments
sub readArgs
{
    # Variables
    my @args = @ARGV;
    my $flag;

    # Read by flag
    while(scalar(@args) > 0)
    {
        $flag = shift(@args);

        # Input file
        if ($flag eq "-i")
        {
            $fileData = shift(@args);
            errExit("LAMMPS data file '$fileData' does not exist.")
                if (! -e $fileData);
        }

        # Types file
        elsif ($flag eq "-t")
        {
            $fileTypes = shift(@args);
            errExit("Data types file '$fileTypes' does not exist.")
                if (! -e $fileTypes);
        }

        # Input script
        elsif ($flag eq "-s")
        {
            $fileInput = shift(@args);
            errExit("Input script '$fileInput' does not exist.")
                if (! -e $fileInput);
        }

        # Output file
        elsif ($flag eq "-o")
        {
            $fileOut = shift(@args);
        }

        # Help/syntax
        elsif ($flag eq "-h")
        {
            printf "Syntax: ./polym_final.pl -i data.lmps -t types.txt ".
                "-s polym.in -o new.lmps\n";
            exit 2;
        }

        # Error
        else
        {
            errExit("Command-line flag '$flag' not recognized.\n".
                "Syntax: ./polym_final.pl -i data.lmps -t types.txt ".
                "-s polym.in -o new.lmps\n");
        }
    }

    # Check values are defined
    errExit("Output file is not defined.") if (!defined($fileOut));
    errExit("Data file is not defined.") if (!defined($fileData));
    errExit("Types file is not defined.") if (!defined($fileTypes));
    errExit("Input script is not defined.") if (!defined($fileInput));
}

# removeCharges( )
# Subtract charges from linker atoms
sub removeCharges
{
    # Variables
    my ($type1, $type2, $q1, $q2, @link1, @link2);

    # Charges
    return if (!defined($inp{'charge'}));
    $q1 = $inp{'charge'}[0];
    $q2 = $inp{'charge'}[1];

    # Numeric atom types of linking atoms
    $type1 = $sys{'atomTypes'}{'num'}{$inp{'link'}[0][0]};
    $type2 = $sys{'atomTypes'}{'num'}{$inp{'link'}[1][0]};
    errExit("At least one linker atom does not have a ".
        "corresponding type in the given types file.")
        if (!defined($type1) || !defined($type2));

    # Find linking atoms in system
    @link1 = Polymatic::group(\@{$sys{'atoms'}{'type'}}, $type1);
    @link2 = Polymatic::group(\@{$sys{'atoms'}{'type'}}, $type2);
    errExit("Linker atoms of both types do not exist in system.")
        if (scalar(@link1) == 0 || scalar(@link2) == 0);

    # Subtract charges to linking atoms
    for (my $i=0; $i < scalar(@link1); $i++) {
        $sys{'atoms'}{'q'}[$link1[$i]] -= $q1;
    }

    for (my $i=0; $i < scalar(@link2); $i++) {
        $sys{'atoms'}{'q'}[$link2[$i]] -= $q2;
    }
}

# unwrap( )
# Unwrap molecules so bonds to span periodic boundaries
sub unwrap
{
    # Variables
    my ($lx, $ly, $lz, $xc, $yc, $zc, $x1, $y1, $z1, $x, $y, $z);
    my (@mols, @atoms, $sep, $init, $closest);

    # Box dimensions
    $lx = @{$sys{'boxDims'}{'x'}}[2];
    $ly = @{$sys{'boxDims'}{'y'}}[2];
    $lz = @{$sys{'boxDims'}{'z'}}[2];

    # Per molecule
    @mols = @{$sys{'mols'}{'atoms'}};
    for (my $i=1; $i < scalar(@mols); $i++)
    {
        # Molecule center
        @atoms = @{$mols[$i]};
        for (my $j=0; $j < scalar(@atoms); $j++)
        {
            $xc += $sys{'atoms'}{'pos'}[$atoms[$j]][0];
            $yc += $sys{'atoms'}{'pos'}[$atoms[$j]][1];
            $zc += $sys{'atoms'}{'pos'}[$atoms[$j]][2];
        }

        $xc = $xc/scalar(@atoms);
        $yc = $yc/scalar(@atoms);
        $zc = $zc/scalar(@atoms);

        # Atom closest to center
        for (my $j=0; $j < scalar(@atoms); $j++)
        {
            ($x1, $y1, $z1) = @{$sys{'atoms'}{'pos'}[$atoms[$j]]};
            $x = $x1 - $xc;
            $y = $y1 - $yc;
            $z = $z1 - $zc;
            $sep = sqrt($x*$x + $y*$y + $z*$z);

            if ($j==0 || $sep < $closest)
            {
                $init = $atoms[$j];
                $closest = $sep;
            }
        }

        Polymatic::unwrapMol(\%sys, $init);
    }
}
