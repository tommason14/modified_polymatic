#!/usr/bin/perl

################################################################################
#
# polym.pl
# This file is part of the Polymatic distribution.
#
# Author: Lauren J. Abbott
# Version: 1.1
# Date: August 16, 2015
#
# Description: Performs a polymerization step for use within the Polymatic code.
# Finds the closest pair of 'linker' atoms satisfying all bonding criteria and
# adds all new bonds, angles, dihedrals, and impropers. Bonding criteria that
# are implemented include a maximum cutoff distance, oriential checks of angles
# between vectors and normal vectors to planes for given atoms, and intra-
# molecular bonding. Options are provided in an input script. Reads and writes
# LAMMPS data files. The Polymatic.pm module is used in this code. It must be in
# the same directory or at a file path recognized by the script (e.g., with 'use
# lib').
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
use List::Util qw( min max );

# Variables
my ($fileData, $fileTypes, $fileInput, $fileOut, $fileMols, $fileBarriers, $molsIncluded);
my (%sys, %inp, @pair, %numberedMols, %barriers);

# Main
readArgs();
%sys = Polymatic::readLammps($fileData);
Polymatic::readTypes($fileTypes, \%sys);
%inp = Polymatic::readPolym($fileInput);
# %numberedMols = Polymatic::readNumberedMols($fileMols);
# %barriers = Polymatic::readBarriers($fileBarriers);
@pair = findPair();
makeUpdates(\@pair);
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
    while (scalar(@args) > 0)
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
            printf "Syntax: ./polym.pl -i data.lmps -t types.txt ".
                "-s polym.in -o new.lmps\n";
            exit 2;
        }

        # Error
        else
        {
            errExit("Command-line flag '$flag' not recognized.\n".
                "Syntax: ./polym.pl -i data.lmps -t types.txt ".
                "-s polym.in -o new.lmps");
        }
    }

    # Check all defined
    errExit("Output file is not defined.") if (!defined($fileOut));
    errExit("Data file is not defined.") if (!defined($fileData));
    errExit("Types file is not defined.") if (!defined($fileTypes));
    errExit("Input script is not defined.") if (!defined($fileInput));
}

# bond_after_considering_propagation_rates(\%sys, $a1, $a2)
sub bond_after_considering_propagation_rates{
  my $sys = $_[0];
  my $a1 = $_[1];
  my $a2 = $_[2];

  # molecule ID from lammps (1..#molecules)
  my $mol_a1 = $sys->{'atoms'}{'mol'}[$a1];
  my $mol_a2 = $sys->{'atoms'}{'mol'}[$a2];
  # molecule number from barriers.in (1 or 2)
  my $mol1 = $numberedMols{$mol_a1};
  my $mol2 = $numberedMols{$mol_a2};

  # read rate from hash created from barriers.in
  my $abbreviation = "k"."$mol1"."$mol2";
  my $rate = $barriers{$abbreviation};

  my $k11 = $barriers{'k11'};
  my $k12 = $barriers{'k12'};
  my $k21 = $barriers{'k21'};
  my $k22 = $barriers{'k22'};

  my @rates = ($k11, $k12, $k21, $k22);

  # reactivity ratios- increase rate of polymerisation for bonds unlikely to react
  my $r1 = $k11 / $k12;
  my $r2 = $k22 / $k21;

  my $bond;

  # check against max rate
  # my $random = rand(1);
  # my $exponent = (exp(-(0.5*($rate - min @rates)) / (1000 * 8.3145 * 298)) - (rand(1) * 0.1 * $r1/$r2));
  # 1 - rate/max rate would be 0 for the max rate, but by doing
  # 1 - (0.99*rate)/max(rate), it is just v.close to 0
  # add some randomness, so do 1-(0.99*rate)/(max rate) * (rand/2)
  my $ratio = (1 - (($rate *0.99) / max @rates)) * rand(1);

  # if all of 1 molecule already in polymer, let any pair bond; the desired
  # composition has already been obtained
  my $curr = Polymatic::molsCurrentlyIncluded($mol1, $molsIncluded);
  $curr = unpack("Q>", substr(("\0"x8) . $curr, -8)); # str -> int
  my $tot = Polymatic::totalMolsIncluded($mol1, $molsIncluded);
  $tot = unpack("Q>", substr(("\0"x8) . $tot, -8));

  if ($curr < $tot){
    $bond = $ratio > 0.1;  # early on
  } else {
    $bond = 1;  # when all of one type have polymerised, automatically bond
  }

  if ($bond){
    # print "Bonding: atom $a1 of mol $mol_a1 (type $mol1) with atom $a2 of mol $mol_a2 (type $mol2)\n";
    return 1;
  } else {
    return 0;
  }
}

# findPair( )
# Find closest pair of linker atoms meeting all bonding criteria
sub findPair
{
    # Variables
    my ($t1, $t2, $a1, $a2, $sep);
    my (@link1, @link2, @pair);
    my $closest = 0;
    my $checks = 0;

    # Numeric atom types of linking atoms
    $t1 = $sys{'atomTypes'}{'num'}{$inp{'link'}[0][0]};
    $t2 = $sys{'atomTypes'}{'num'}{$inp{'link'}[1][0]};
    errExit("At least one linker atom does not have a ".
        "corresponding type in the given types file.")
        if (!defined($t1) || !defined($t2));

    # Find linking atoms in system
    @link1 = Polymatic::group(\@{$sys{'atoms'}{'type'}}, $t1);
    @link2 = Polymatic::group(\@{$sys{'atoms'}{'type'}}, $t2);
    errExit("Linker atoms of both types do not exist in system.")
        if (scalar(@link1) == 0 || scalar(@link2) == 0);

    # Check bonding criteria for all pairs
    for (my $i=0; $i < scalar(@link1); $i++)
    {
        $a1 = $link1[$i];

        for (my $j=0; $j < scalar(@link2); $j++)
        {
            $a2 = $link2[$j];

            # Intra check
            if (defined($inp{'intra'})) {
                next if (intra([$a1, $a2]));
            } else {
                next if ($sys{'atoms'}{'mol'}[$a1] ==
                    $sys{'atoms'}{'mol'}[$a2]);
            }

            # Cutoff check
            $sep = Polymatic::getSep(\%sys, $a1, $a2);
            next if ($sep > $inp{'cutoff'});

            # Alignment check
            if (defined($inp{'align'})) {
                next if (!aligned([$a1, $a2]));
            }

            # # Check propagation rate
            # if (!bond_after_considering_propagation_rates(\%sys, $a1, $a2)){
            #   next;
            # }

            # Save pair if closest
            if ($closest == 0 || $sep < $closest) {

                $closest = $sep;
                @pair = ($a1, $a2);
            }
            
        }
    }

    # Return closest pair
    if ($closest == 0) {
        exit 3;
    } else {
        printf "  Pair: %.2f A (%d,%d)\n", $closest, $pair[0], $pair[1];
        return @pair;
    }
}

# makeUpdates( \@pair )
# Update connectivity of system with new bond between given pair
sub makeUpdates
{
    # Variables
    my @pair = @{$_[0]};
    my ($t1, $t1n, $t2, $t2n, $a1, $a2, $a3, $a4);
    my ($type, $order, $num, $min, $max);
    my (@addBonds, @addAngles, @addDiheds, @addImprops);
    my (@update, @a, @bonded1, @bonded2, @bonded3);
    my $at = 'atomTypes';

    # New bonds
    push(@addBonds, [@pair]);
    push(@addBonds, getExtraBonds(\@pair))
        if (defined($inp{'bond'}));

    # Numeric atom types for linker atoms
    $t1 = $sys{$at}{'num'}{$inp{'link'}[0][0]};
    $t1n = $sys{$at}{'num'}{$inp{'link'}[0][1]};
    $t2 = $sys{$at}{'num'}{$inp{'link'}[1][0]};
    $t2n = $sys{$at}{'num'}{$inp{'link'}[1][1]};
    errExit("At least one linker atom does not have a ".
        "corresponding type in the given types file.")
        if (!defined($t1n) || !defined($t2n));

    # Update atoms in bonds
    for (my $i=0; $i < scalar(@addBonds); $i++)
    {
        # Atoms in bond
        ($a1, $a2) = @{$addBonds[$i]};
        push(@{$sys{'atoms'}{'bonded'}[$a1]}, $a2);
        push(@{$sys{'atoms'}{'bonded'}[$a2]}, $a1);

        # Charges
        if (defined($inp{'charge'}) &&
            $sys{'atoms'}{'type'}[$a1] == $t1 &&
            $sys{'atoms'}{'type'}[$a2] == $t2)
        {
            $sys{'atoms'}{'q'}[$a1] -= $inp{'charge'}[0];
            $sys{'atoms'}{'q'}[$a2] -= $inp{'charge'}[1];
        }

        # Atom types
        $sys{'atoms'}{'type'}[$a1] = $t1n
            if ($sys{'atoms'}{'type'}[$a1] == $t1);
        $sys{'atoms'}{'type'}[$a2] = $t2n
            if ($sys{'atoms'}{'type'}[$a2] == $t2);
    }

    # Update affected bonded terms
    for (my $i=0; $i < scalar(@addBonds); $i++)
    {
        # Atoms in bond
        ($a1, $a2) = @{$addBonds[$i]};

        # Bonds
        @update = Polymatic::uniqueArray((@{$sys{'atoms'}{'bonds'}[$a1]},
            @{$sys{'atoms'}{'bonds'}[$a2]}));
        for (my $j=0; $j < scalar(@update); $j++)
        {
            @a = @{$sys{'bonds'}{'atoms'}[$update[$j]]};
            $type = Polymatic::getBondType(\%sys, \@a);

            if ($type > 0)
            {
                $sys{'bonds'}{'type'}[$update[$j]] = $type;
            }
            elsif ($type < 0)
            {
                $sys{'bonds'}{'type'}[$update[$j]] = -1*$type;
                $sys{'bonds'}{'atoms'}[$update[$j]] = [reverse(@a)];
            }
            else
            {
                errExit("Bond type '".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[0]]].",".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[1]]].
                    "' is not defined.");
            }
        }

        # Angles
        @update = Polymatic::uniqueArray((@{$sys{'atoms'}{'angles'}[$a1]},
            @{$sys{'atoms'}{'angles'}[$a2]}));
        for (my $j=0; $j < scalar(@update); $j++)
        {
            @a = @{$sys{'angles'}{'atoms'}[$update[$j]]};
            $type = Polymatic::getAngleType(\%sys, \@a);

            if ($type > 0)
            {
                $sys{'angles'}{'type'}[$update[$j]] = $type;
            }
            elsif ($type < 0)
            {
                $sys{'angles'}{'type'}[$update[$j]] = -1*$type;
                $sys{'angles'}{'atoms'}[$update[$j]] = [reverse(@a)];
            }
            else
            {
                errExit("Angle type '".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[0]]].",".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[1]]].",".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[2]]].
                    "' is not defined.");
            }
        }

        # Dihedrals
        @update = Polymatic::uniqueArray((@{$sys{'atoms'}{'diheds'}[$a1]},
            @{$sys{'atoms'}{'diheds'}[$a2]}));
        for (my $j=0; $j < scalar(@update); $j++)
        {
            @a = @{$sys{'diheds'}{'atoms'}[$update[$j]]};
            $type = Polymatic::getDihedType(\%sys, \@a);

            if ($type > 0)
            {
                $sys{'diheds'}{'type'}[$update[$j]] = $type;
            }
            elsif ($type < 0)
            {
                $sys{'diheds'}{'type'}[$update[$j]] = -1*$type;
                $sys{'diheds'}{'atoms'}[$update[$j]] = [reverse(@a)];
            }
            else
            {
                errExit("Dihedral type '".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[0]]].",".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[1]]].",".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[2]]].",".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[3]]].
                    "' is not defined.");
            }
        }

        # Impropers
        @update = Polymatic::uniqueArray((@{$sys{'atoms'}{'improps'}[$a1]},
            @{$sys{'atoms'}{'improps'}[$a2]}));
        for (my $j=0; $j < scalar(@update); $j++)
        {
            @a = @{$sys{'improps'}{'atoms'}[$update[$j]]};
            ($type, $order) = Polymatic::getImpropType(\%sys, \@a);

            if ($type > 0)
            {
                $sys{'improps'}{'type'}[$update[$j]] = $type;
                if ($order == 1) {
                    $sys{'improps'}{'atoms'}[$update[$j]] =
                        [$a[0], $a[1], $a[3], $a[2]];
                } elsif ($order == 2) {
                    $sys{'improps'}{'atoms'}[$update[$j]] =
                        [$a[2], $a[1], $a[0], $a[3]];
                } elsif ($order == 3) {
                    $sys{'improps'}{'atoms'}[$update[$j]] =
                        [$a[2], $a[1], $a[3], $a[0]];
                } elsif ($order == 4) {
                    $sys{'improps'}{'atoms'}[$update[$j]] =
                        [$a[3], $a[1], $a[0], $a[2]];
                } elsif ($order == 5) {
                    $sys{'improps'}{'atoms'}[$update[$j]] =
                        [$a[3], $a[1], $a[2], $a[0]];
                }
            }
            else
            {
                errExit("Improper type '".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[0]]].",".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[1]]].",".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[2]]].",".
                    $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[3]]].
                    "' is not defined.");
            }
        }
    }

    # Add new bonded terms
    for (my $i=0; $i < scalar(@addBonds); $i++)
    {
        # Skip if no angles/diheds/improps
        last if ($sys{'angleTypes'}{'count'} == 0 &&
            $sys{'dihedTypes'}{'count'} == 0 &&
            $sys{'impropTypes'}{'count'} == 0);

        # Atoms in bond
        ($a1, $a2) = @{$addBonds[$i]};
        @bonded1 = @{$sys{'atoms'}{'bonded'}[$a1]};
        @bonded2 = @{$sys{'atoms'}{'bonded'}[$a2]};

        for (my $j=0; $j < scalar(@bonded2); $j++)
        {
            # Angle 1,2,x
            $a3 = $bonded2[$j];
            next if ($a3 == $a1);
            push(@addAngles, [$a1, $a2, $a3])
                if ($sys{'angleTypes'}{'count'} > 0);

            # Dihedral 1,2,x,y
            @bonded3 = @{$sys{'atoms'}{'bonded'}[$a3]};
            for (my $k=0; $k < scalar(@bonded3); $k++)
            {
                $a4 = $bonded3[$k];
                next if ($a4 == $a2);
                push(@addDiheds, [$a1, $a2, $a3, $a4])
                    if ($sys{'dihedTypes'}{'count'} > 0);
            }
        }

        for (my $j=0; $j < scalar(@bonded1); $j++)
        {
            # Angle 2,1,x
            $a3 = $bonded1[$j];
            next if ($a3 == $a2);
            push(@addAngles, [$a2, $a1, $a3])
                if ($sys{'angleTypes'}{'count'} > 0);

            # Dihedral 2,1,x,y
            @bonded3 = @{$sys{'atoms'}{'bonded'}[$a3]};
            for (my $k=0; $k < scalar(@bonded3); $k++)
            {
                $a4 = $bonded3[$k];
                next if ($a4 == $a1);
                push(@addDiheds, [$a2, $a1, $a3, $a4])
                    if ($sys{'dihedTypes'}{'count'} > 0);
            }

            # Dihedral x,1,2,y
            for (my $k=0; $k < scalar(@bonded2); $k++)
            {
                $a4 = $bonded2[$k];
                next if ($a4 == $a1);
                push(@addDiheds, [$a3, $a1, $a2, $a4])
                    if ($sys{'dihedTypes'}{'count'} > 0);
            }
        }

        if ($sys{'impropTypes'}{'count'} > 0 && scalar(@bonded1) > 1)
        {
            # Improper 2,1,x,y
            for (my $j=0; $j < scalar(@bonded1)-1; $j++)
            {
                $a3 = $bonded1[$j];
                next if ($a3 == $a2);

                for (my $k=$j+1; $k < scalar(@bonded1); $k++)
                {
                    $a4 = $bonded1[$k];
                    next if ($a4 == $a2);
                    push(@addImprops, [$a2, $a1, $a3, $a4]);
                }
            }
        }

        if ($sys{'impropTypes'}{'count'} > 0 && scalar(@bonded2) > 1)
        {
            # Improper 1,2,x,y
            for (my $j=0; $j < scalar(@bonded1)-1; $j++)
            {
                $a3 = $bonded2[$j];
                next if ($a3 == $a1);

                for (my $k=$j+1; $k < scalar(@bonded2); $k++)
                {
                    $a4 = $bonded2[$k];
                    next if ($a4 == $a1);
                    push(@addImprops, [$a1, $a2, $a3, $a4]);
                }
            }
        }
    }

    # Delete duplicates
    @addBonds = Polymatic::delDupBonds(\@addBonds);
    @addAngles = Polymatic::delDupBonds(\@addAngles);
    @addDiheds = Polymatic::delDupBonds(\@addDiheds);
    @addImprops = Polymatic::delDupImprops(\@addImprops);

    # New bonds
    $num = $sys{'bonds'}{'count'};
    for (my $i=0; $i < scalar(@addBonds); $i++)
    {
        $num++;
        @a = @{$addBonds[$i]};
        $type = Polymatic::getBondType(\%sys, \@a);

        if ($type > 0)
        {
            $sys{'bonds'}{'atoms'}[$num] = [@a];
            $sys{'bonds'}{'type'}[$num] = $type;
        }
        elsif ($type < 0)
        {
            $sys{'bonds'}{'atoms'}[$num] = [reverse(@a)];
            $sys{'bonds'}{'type'}[$num] = -1*$type;
        }
        else
        {
            errExit("Bond type '".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[0]]].",".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[1]]].
                "' is not defined.");
        }
    }
    $sys{'bonds'}{'count'} = $num;

    # New angles
    $num = $sys{'angles'}{'count'};
    for (my $i=0; $i < scalar(@addAngles); $i++)
    {
        $num++;
        @a = @{$addAngles[$i]};
        $type = Polymatic::getAngleType(\%sys, \@a);

        if ($type > 0)
        {
            $sys{'angles'}{'atoms'}[$num] = [@a];
            $sys{'angles'}{'type'}[$num] = $type;
        }
        elsif ($type < 0)
        {
            $sys{'angles'}{'atoms'}[$num] = [reverse(@a)];
            $sys{'angles'}{'type'}[$num] = -1*$type;
        }
        else
        {
            errExit("Angle type '".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[0]]].",".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[1]]].",".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[2]]].
                "' is not defined.");
        }
    }
    $sys{'angles'}{'count'} = $num;

    # New dihedrals
    $num = $sys{'diheds'}{'count'};
    for (my $i=0; $i < scalar(@addDiheds); $i++)
    {
        $num++;
        @a = @{$addDiheds[$i]};
        $type = Polymatic::getDihedType(\%sys, \@a);

        if ($type > 0)
        {
            $sys{'diheds'}{'atoms'}[$num] = [@a];
            $sys{'diheds'}{'type'}[$num] = $type;
        }
        elsif ($type < 0)
        {
            $sys{'diheds'}{'atoms'}[$num] = [reverse(@a)];
            $sys{'diheds'}{'type'}[$num] = -1*$type;
        }
        else
        {
            errExit("Dihedral type '".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[0]]].",".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[1]]].",".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[2]]].",".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[3]]].
                "' is not defined.");
        }
    }
    $sys{'diheds'}{'count'} = $num;

    # New impropers
    $num = $sys{'improps'}{'count'};
    for (my $i=0; $i < scalar(@addImprops); $i++)
    {
        $num++;
        @a = @{$addImprops[$i]};
        ($type, $order) = Polymatic::getImpropType(\%sys, \@a);

        if ($type > 0)
        {
            $sys{'improps'}{'type'}[$num] = $type;
            if ($order == 0) {
                $sys{'improps'}{'atoms'}[$num] = [@a];
            } elsif ($order == 1) {
                $sys{'improps'}{'atoms'}[$num] = [$a[0], $a[1], $a[3], $a[2]];
            } elsif ($order == 2) {
                $sys{'improps'}{'atoms'}[$num] = [$a[2], $a[1], $a[0], $a[3]];
            } elsif ($order == 3) {
                $sys{'improps'}{'atoms'}[$num] = [$a[2], $a[1], $a[3], $a[0]];
            } elsif ($order == 4) {
                $sys{'improps'}{'atoms'}[$num] = [$a[3], $a[1], $a[0], $a[2]];
            } elsif ($order == 5) {
                $sys{'improps'}{'atoms'}[$num] = [$a[3], $a[1], $a[2], $a[0]];
            }
        }
        else
        {
            errExit("Improper type '".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[0]]].",".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[1]]].",".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[2]]].",".
                $sys{$at}{'name'}[$sys{'atoms'}{'type'}[$a[3]]].
                "' is not defined.");
        }
    }
    $sys{'improps'}{'count'} = $num;

    # Molecule numbers
    if ($sys{'atoms'}{'mol'}[$pair[0]] != $sys{'atoms'}{'mol'}[$pair[1]])
    {
        $num = $sys{'mols'}{'count'};

        # Max and min
        if ($sys{'atoms'}{'mol'}[$pair[0]] > $sys{'atoms'}{'mol'}[$pair[1]]) {
            $min = $sys{'atoms'}{'mol'}[$pair[1]];
            $max = $sys{'atoms'}{'mol'}[$pair[0]];
        } else {
            $min = $sys{'atoms'}{'mol'}[$pair[0]];
            $max = $sys{'atoms'}{'mol'}[$pair[1]];
        }

        # Move max to min
        @a = @{$sys{'mols'}{'atoms'}[$max]};
        for (my $i=0; $i < scalar(@a); $i++)
        {
            $sys{'atoms'}{'mol'}[$a[$i]] = $min;
            push(@{$sys{'mols'}{'atoms'}[$min]}, $a[$i]);
        }
        $sys{'mols'}{'atoms'}[$max] = [];

        # Move last to max
        @a = @{$sys{'mols'}{'atoms'}[$num]};
        for (my $i=0; $i < scalar(@a); $i++)
        {
            $sys{'atoms'}{'mol'}[$a[$i]] = $max;
            push(@{$sys{'mols'}{'atoms'}[$max]}, $a[$i]);
        }
        pop(@{$sys{'mols'}{'atoms'}});
        $sys{'mols'}{'count'} = $num-1;
    }
}

# getConnect( \@pair )
# Define connectivity from input script for the given pair
sub getConnect
{
    # Variables
    my @pair = @{$_[0]};
    my ($n, $a1, $a2);
    my (@atoms, @queue, @connect, @bonded);

    # Make sure connect and types are defined
    errExit("The 'connect' and 'types' commands must be defined ".
        "to use the 'vector', 'plane', or 'bond' commands.")
        if (!defined($inp{'connect'}) || !defined($inp{'types'}));

    # Initial atoms
    $atoms[1] = $pair[0];
    $atoms[2] = $pair[1];
    push(@queue, (1,2));

    # Follow connectivity
    while (scalar(@queue) > 0)
    {
        $n = shift(@queue);
        @connect = @{$inp{'connect'}[$n]};
        @bonded = @{$sys{'atoms'}{'bonded'}[$atoms[$n]]};

        for (my $i=0; $i < scalar(@bonded); $i++)
        {
            $a1 = $bonded[$i];
            for (my $j=0; $j < scalar(@connect); $j++)
            {
                $a2 = $connect[$j];
                if ($sys{'atomTypes'}{'name'}[$sys{'atoms'}{'type'}[$a1]] eq
                    $inp{'types'}[$a2])
                {
                    errExit("Atom connectivity in input script is not unique.")
                        if (defined($atoms[$a2]));
                    $atoms[$a2] = $a1;
                    push(@queue, $a2) if (defined($inp{'connect'}[$a2]));
                }
            }
        }
    }

    # Return atom definitions
    return @atoms;
}

# intra( \@pair )
# Check if given pair is within N bonds
sub intra
{
    # Variables
    my ($a1, $a2) = @{$_[0]};
    my ($b1, $b2, $b3, $b4, $b5);
    my (@bonded1, @bonded2, @bonded3, @bonded4, @bonded5);
    my $num = $inp{'intra'};

    # One bond
    @bonded1 = @{$sys{'atoms'}{'bonded'}[$a1]};
    for (my $i1=0; $i1 < scalar(@bonded1); $i1++)
    {
        $b1 = $bonded1[$i1];
        return 1 if ($b1 == $a2);
        next if ($num == 1);

        # Two bonds
        @bonded2 = @{$sys{'atoms'}{'bonded'}[$b1]};
        for (my $i2=0; $i2 < scalar(@bonded2); $i2++)
        {
            $b2 = $bonded2[$i2];
            return 1 if ($b2 == $a2);
            next if ($b2 == $a1);
            next if ($num == 2);

            # Three bonds
            @bonded3 = @{$sys{'atoms'}{'bonded'}[$b2]};
            for (my $i3=0; $i3 < scalar(@bonded3); $i3++)
            {
                $b3 = $bonded3[$i3];
                return 1 if ($b3 == $a2);
                next if ($b3 == $b1);
                next if ($num == 3);

                # Four bonds
                @bonded4 = @{$sys{'atoms'}{'bonded'}[$b3]};
                for (my $i4=0; $i4 < scalar(@bonded4); $i4++)
                {
                    $b4 = $bonded4[$i4];
                    return 1 if ($b4 == $a2);
                    next if ($b4 == $b2);
                    next if ($num == 4);

                    # Five bonds
                    @bonded5 = @{$sys{'atoms'}{'bonded'}[$b4]};
                    for (my $i5=0; $i5 < scalar(@bonded5); $i5++)
                    {
                        $b5 = $bonded5[$i5];
                        return 1 if ($b5 == $a2);
                        next if ($b5 == $b3);
                    }
                }
            }
        }
    }

    # Not bonded within N bonds
    return 0;
}

# aligned( \@pair )
# Perform alignment checks for given pair
sub aligned
{
    # Variables
    my @pair = @{$_[0]};
    my ($ang, @atoms, @cond);
    my (@a1, @a2, @v1, @v2, @pos1, @pos2);

    # Atom definitions
    @atoms = getConnect(\@pair);

    # Loop through checks
    for (my $i=0; $i < scalar(@{$inp{'align'}{'atoms'}}); $i++)
    {
        @a1 = @{$inp{'align'}{'atoms'}[$i][0]};
        @a2 = @{$inp{'align'}{'atoms'}[$i][1]};
        @cond = @{$inp{'align'}{'cond'}[$i]};
        @pos1 = ();
        @pos2 = ();

        # Get atoms based on atom definition
        for (my $j=0; $j < scalar(@a1); $j++)
        {
            errExit("Atom connectivity in input script ".
                "is not properly defined.")
                if (!defined($atoms[$a1[$j]]));
            $a1[$j] = $atoms[$a1[$j]];
            push(@pos1, [@{$sys{'atoms'}{'pos'}[$a1[$j]]}]);
        }

        for (my $j=0; $j < scalar(@a2); $j++)
        {
            errExit("Atom connectivity in input script ".
                "is not properly defined.")
                if (!defined($atoms[$a2[$j]]));
            $a2[$j] = $atoms[$a2[$j]];
            push(@pos2, [@{$sys{'atoms'}{'pos'}[$a2[$j]]}]);
        }

        # Vector or plane
        if (scalar(@a1) == 2) {
            @v1 = Polymatic::vectorSub(\@{$pos1[0]}, \@{$pos1[1]});
        } else {
            @v1 = Polymatic::normalPlane(\@pos1);
        }

        if (scalar(@a2) == 2) {
            @v2 = Polymatic::vectorSub(\@{$pos2[0]}, \@{$pos2[1]});
        } else {
            @v2 = Polymatic::normalPlane(\@pos2);
        }

        # Angle between vectors
        $ang = Polymatic::vectorAng(\@v1, \@v2) * 180/Math::Trig::pi;

        # Return 0 if not aligned
        if (scalar(@cond) == 1) {
            return 0 if (!eval("$ang $cond[0]"));
        } elsif (scalar(@cond) == 3) {
            return 0 if (!eval("$ang $cond[0] $cond[2] $ang $cond[1]"));
        }
    }

    # Return 1 if aligned
    return 1;
}

# getExtraBonds( \@pair )
# Define extra bonds based on the given pair
sub getExtraBonds
{
    # Variables
    my @pair = @{$_[0]};
    my ($a1, $a2, $sep, @atoms, @addBonds);

    # Atom definitions
    @atoms = getConnect(\@pair);

    # Extra bonds
    for (my $i=0; $i < scalar(@{$inp{'bond'}}); $i++)
    {
        $a1 = $atoms[$inp{'bond'}[$i][0]];
        $a2 = $atoms[$inp{'bond'}[$i][1]];
        errExit("Atom connectivity in input script is not properly defined.")
            if (!defined($a1) || !defined($a2));

        $sep = Polymatic::getSep(\%sys, $a1, $a2);
        printf "  Extra bond: %.2f A (%d,%d)\n", $sep, $a1, $a2;
        push(@addBonds, [$a1, $a2]);
    }

    # Return extra bonds
    return @addBonds;
}
