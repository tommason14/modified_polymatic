#!/usr/bin/perl

################################################################################
#
# pack.pl
# This file is part of the Polymatic distribution.
#
# Author: Lauren J. Abbott
# Version: 1.1
# Date: August 16, 2015
#
# Description: Creates a random packing of molecules in a periodic cubic cell.
# Random insertions are made for each molecule to avoid overlap of the atomic
# radii. The molecule structures to be packed are defined in reference LAMMPS
# data file and the final packed box is also given as a LAMMPS data file. The
# Polymatic.pm module is used in this code. It must be in the same directory or
# at a file path recognized by the script (e.g., with 'use lib').
#
# Syntax:
#  ./pack.pl -i num F1.lmps N1 F2.lmps N2 (...)
#            -l boxL
#            -o pack.lmps
#
# Arguments:
#  1. num, number of reference molecules (-i)
#  2. Fi.lmps, LAMMPS data file of the ith reference molecule (-i)
#  3. Ni, number of reference molecule Fi.lmps to pack (-i)
#  4. boxL, box length of cubic box (-l)
#  5. pack.lmps, name of LAMMPS data file to output packed box (-o)
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
use Clone();

# User Variables
my $maxAttempts = 100000;
my $scale = 1.0;
my $radius = 2.0;

# Global Variables
my ($numRefMols, $filePack, $boxL);
my (@files, @numbers, @refs, %sys);

# Main
readArgs();
defineSystems();
packMols();
Polymatic::writeLammps($filePack, \%sys);

################################################################################
# Subroutines

# errExit( $error )
# Exit program and print error
sub errExit
{
    printf "Error: %s\n", $_[0];
    exit 2;
}

# warning( $warning )
# Print warning but do not exit
sub warning
{
    printf "Warning: %s\n", $_[0];
}

# readArgs( )
# Read command-line arguments
sub readArgs
{
    # Variables
    my @args = @ARGV;
    my $flag;

    while (scalar(@args) > 0)
    {
        $flag = shift(@args);

        # Input files
        if ($flag eq "-i")
        {
            $numRefMols = shift(@args);
            errExit("Number of reference molecules must be greater than 0")
                if ($numRefMols <= 0);

            for (my $i=0; $i < $numRefMols; $i++)
            {
                push(@files, shift(@args));
                push(@numbers, shift(@args));
                errExit("Reference molecule file '$files[$i]' does not exist.")
                    if (! -e $files[$i]);
                errExit("Number of molecules to pack must be greater than 0.")
                    if ($numbers[$i] < 0);
                errExit("Number of molecules to pack must be greater than 0")
                    if ($i > 0 && $numbers[$i] <= 0);
            }
        }

        # Output file
        elsif ($flag eq "-o")
        {
            $filePack = shift(@args);
        }

        # Box length
        elsif ($flag eq "-l")
        {
            $boxL = shift(@args);
        }

        # Help/syntax
        elsif ($flag eq "-h")
        {
            printf "Syntax: ./pack.pl -i num F1.lmps N1 F2.lmps N2 (...) ".
                "-l boxL -o pack.lmps\n";
            exit 0;
        }

        # Unrecognized flag
        else
        {
            errExit("Command-line flag '$flag' not recognized.\n".
                "Syntax: ./pack.pl -i num F1.lmps N1 F2.lmps N2 (...) ".
                "-l boxL -o pack.lmps");
        }
    }

    # Check required values are defined
    errExit("Output file name is not defined.")
        if (!defined($filePack));
    errExit("Box length is not defined.")
        if (!defined($boxL));
    errExit("Box length must be greater than 0.")
        if ($boxL <= 0 && $numbers[0] > 0);
    errExit("Box length must be greater than 0.")
        if ($boxL < 0 && $numbers[0] == 0);
    errExit("Reference molecule files are not defined.")
        if (!defined($numRefMols));
}

# defineSystems( )
# Define molecular systems from reference files and setup packing system
sub defineSystems
{
    # Variables
    my (%ref, $num, $diam, $cellN, $cellL);

    # Read reference molecules
    for (my $i=0; $i < $numRefMols; $i++)
    {
        %ref = Polymatic::readLammps($files[$i]);
        push(@refs, \%{Clone::clone(\%ref)});
        $num = $numbers[$i];

        # If first molecule is initial existing system
        if ($i==0 && $num == 0)
        {
            %sys = %{Clone::clone(\%ref)};
            $sys{'header'} = "Random packing of '$files[$i]'";

            if ($boxL == 0)
            {
                $boxL = $sys{'boxDims'}{'x'}[2];
                errExit("Only cubic boxes are currently supported.")
                    if ($sys{'boxDims'}{'x'}[2] != $sys{'boxDims'}{'y'}[2] ||
                        $sys{'boxDims'}{'x'}[2] != $sys{'boxDims'}{'z'}[2]);
            }
            else
            {
                $sys{'boxDims'}{'x'} = [0, $boxL, $boxL];
                $sys{'boxDims'}{'y'} = [0, $boxL, $boxL];
                $sys{'boxDims'}{'z'} = [0, $boxL, $boxL];
            }
        }

        # If first molecule is not initial existing system
        elsif ($i==0 && $num > 0)
        {
            $sys{'header'} = "Random packing of $numbers[$i] x '$files[$i]'";
            $sys{'boxDims'}{'x'} = [0, $boxL, $boxL];
            $sys{'boxDims'}{'y'} = [0, $boxL, $boxL];
            $sys{'boxDims'}{'z'} = [0, $boxL, $boxL];
            $sys{'mols'}{'count'} = 0;

            $sys{'atomTypes'} = \%{Clone::clone(\%{$ref{'atomTypes'}})}
                if (defined($ref{'atomTypes'}));
            $sys{'bondTypes'} = \%{Clone::clone(\%{$ref{'bondTypes'}})}
                if (defined($ref{'bondTypes'}));
            $sys{'angleTypes'} = \%{Clone::clone(\%{$ref{'angleTypes'}})}
                if (defined($ref{'angleTypes'}));
            $sys{'dihedTypes'} = \%{Clone::clone(\%{$ref{'dihedTypes'}})}
                if (defined($ref{'dihedTypes'}));
            $sys{'impropTypes'} = \%{Clone::clone(\%{$ref{'impropTypes'}})}
                if (defined($ref{'impropTypes'}));
            $sys{'flags'} = \%{Clone::clone(\%{$ref{'flags'}})}
                if (defined($ref{'flags'}));
        }

        # After first molecule
        else
        {
            $sys{'header'} = $sys{'header'}.", $numbers[$i] x '$files[$i]'";
            errExit("Number of atom types don't match in reference files.")
                if ($sys{'atomTypes'}{'count'} != $ref{'atomTypes'}{'count'});
            errExit("Number of bond types don't match in reference files.")
                if ($sys{'bondTypes'}{'count'} != $ref{'bondTypes'}{'count'});
            errExit("Number of angle types don't match in reference files.")
                if ($sys{'angleTypes'}{'count'} != $ref{'angleTypes'}{'count'});
            errExit("Number of dihedral types don't match in reference files.")
                if ($sys{'dihedTypes'}{'count'} != $ref{'dihedTypes'}{'count'});
            errExit("Number of improper types don't match in reference files.")
                if ($sys{'impropTypes'}{'count'} !=
                    $ref{'impropTypes'}{'count'});
        }
    }

    # Initialize neighbor list
    $diam = 0.0;
    for (my $i=1; $i <= $sys{'atomTypes'}{'count'}; $i++)
    {
        $diam = $sys{'atomTypes'}{'coeffs'}{'pair'}[$i][1]
            if (defined($sys{'atomTypes'}{'coeffs'}{'pair'}[$i][1]) &&
                $sys{'atomTypes'}{'coeffs'}{'pair'}[$i][1] > $diam);
    }

    if ($diam == 0.0)
    {
        $diam = $radius*2;
        warning("Could not find sigma values, using default radius.");
    }

    $cellN = POSIX::floor($boxL/$diam);
    $cellL = $boxL/$cellN;
    Polymatic::initNeighList(\%sys, $cellN, $cellL);
}

# packMols( )
# Loop to control random packing of reference molecules into system
sub packMols
{
    # Variables
    my (%ref, %test, $num, $count);
    my $attempts = 0;

    # Pack each molecule type
    for (my $i=0; $i < $numRefMols; $i++)
    {
        %ref = %{$refs[$i]};
        $num = $numbers[$i];
        centerMol(\%ref) if ($num > 0);
        printf "Packing %d molecules of reference '%s':\n", $num, $files[$i];

        $count = 0;
        while ($count < $num)
        {
            %test = generateRandomMol(\%ref);
            $attempts++;

            if (checkOverlap(\%test))
            {
                addMol(\%test);
                $count++;
                printf "  %d, %d attempts\n", $count, $attempts;
                $attempts = 0;
            }

            errExit("Reached maximum number of attempts. ".
                "Try adjusting box size or scaling parameter.")
                if ($attempts == $maxAttempts);
        }
    }
}

# centerMol( \%mol )
# Translate molecule to be centered at origin
sub centerMol
{
    # Variables
    my $mol = $_[0];
    my ($num, $xc, $yc, $zc, $xt, $yt, $zt);

    # Calculate geometrical center
    $num = $mol->{'atoms'}{'count'};
    for (my $i=1; $i <= $num; $i++)
    {
        $xt += $mol->{'atoms'}{'pos'}[$i][0];
        $yt += $mol->{'atoms'}{'pos'}[$i][1];
        $zt += $mol->{'atoms'}{'pos'}[$i][2];
    }

    $xc = $xt/$num;
    $yc = $yt/$num;
    $zc = $zt/$num;

    # Translate atom positions
    for (my $i=1; $i <= $num; $i++)
    {
        $mol->{'atoms'}{'pos'}[$i][0] -= $xc;
        $mol->{'atoms'}{'pos'}[$i][1] -= $yc;
        $mol->{'atoms'}{'pos'}[$i][2] -= $zc;
    }
}

# generateRandomMol( \%mol )
# Generate molecule from reference with random translation and rotation
sub generateRandomMol
{
    # Variables
    my $mol = $_[0];
    my (%new, $xt, $yt, $zt, $xi, $yi, $zi, $x, $y, $z, $num);
    my ($a, $b, $c, $m1, $m2, $m3, $n1, $n2, $n3, $p1, $p2, $p3);

    # Initialize new molecule
    %new = %{Clone::clone($mol)};
    $new{'neighList'}{'num'} = $sys{'neighList'}{'num'};
    $new{'neighList'}{'width'} = $sys{'neighList'}{'width'};

    # Translation vector
    $xt = rand($boxL) + $sys{'boxDims'}{'x'}[0];
    $yt = rand($boxL) + $sys{'boxDims'}{'y'}[0];
    $zt = rand($boxL) + $sys{'boxDims'}{'z'}[0];

    # Rotation matrix
    $a = rand(2*Math::Trig::pi);
    $b = rand(2*Math::Trig::pi);
    $c = rand(2*Math::Trig::pi);
    $m1 = cos($b)*cos($c);
    $m2 = cos($b)*sin($c);
    $m3 = -sin($b);
    $n1 = -cos($a)*sin($c) + sin($a)*sin($b)*cos($c);
    $n2 = cos($a)*cos($c) + sin($a)*sin($b)*sin($c);
    $n3 = sin($a)*cos($b);
    $p1 = sin($a)*sin($c) + cos($a)*sin($b)*cos($c);
    $p2 = -sin($a)*cos($c) + cos($a)*sin($b)*sin($c);
    $p3 = cos($a)*cos($b);

    # Translate and rotate each atom
    $num = $new{'atoms'}{'count'};
    for (my $i=1; $i <= $num; $i++)
    {
        ($xi, $yi, $zi) = @{$new{'atoms'}{'pos'}[$i]};
        $x = $m1*$xi + $m2*$yi + $m3*$zi + $xt;
        $y = $n1*$xi + $n2*$yi + $n3*$zi + $yt;
        $z = $p1*$xi + $p2*$yi + $p3*$zi + $zt;
        $new{'atoms'}{'pos'}[$i] = [$x, $y, $z];
    }

    # Return random molecule
    return %new;
}

# addMol( \%mol )
# Add molecule to system
sub addMol
{
    # Variables
    my $mol = $_[0];
    my ($initMol, $init, $num, $n, $a1, $a2, $a3, $a4);

    # Initialize
    $initMol = $sys{'mols'}{'count'};
    $init = $sys{'atoms'}{'count'};

    # Add atoms
    $n = $init;
    $num = $mol->{'atoms'}{'count'};
    for (my $i=1; $i <= $num; $i++)
    {
        $n++;
        $sys{'atoms'}{'mol'}[$n] = $initMol + $mol->{'atoms'}{'mol'}[$i];
        $sys{'atoms'}{'type'}[$n] = $mol->{'atoms'}{'type'}[$i];
        $sys{'atoms'}{'q'}[$n] = $mol->{'atoms'}{'q'}[$i];
        $sys{'atoms'}{'pos'}[$n] = [@{$mol->{'atoms'}{'pos'}[$i]}];
        Polymatic::addNeighList(\%sys, $n);
    }

    $sys{'atoms'}{'count'} = $n;
    $sys{'mols'}{'count'} += $mol->{'mols'}{'count'};

    # Add bonds
    $n = $sys{'bonds'}{'count'};
    $num = $mol->{'bonds'}{'count'};
    for (my $i=1; $i <= $num; $i++)
    {
        $n++;
        $sys{'bonds'}{'type'}[$n] = $mol->{'bonds'}{'type'}[$i];
        ($a1, $a2) = @{$mol->{'bonds'}{'atoms'}[$i]};
        $sys{'bonds'}{'atoms'}[$n] = [$a1+$init, $a2+$init];
    }

    $sys{'bonds'}{'count'} = $n;

    # Add angles
    $n = $sys{'angles'}{'count'};
    $num = $mol->{'angles'}{'count'};
    for (my $i=1; $i <= $num; $i++)
    {
        $n++;
        $sys{'angles'}{'type'}[$n] = $mol->{'angles'}{'type'}[$i];
        ($a1, $a2, $a3) = @{$mol->{'angles'}{'atoms'}[$i]};
        $sys{'angles'}{'atoms'}[$n] = [$a1+$init, $a2+$init, $a3+$init];
    }

    $sys{'angles'}{'count'} = $n;

    # Add dihedrals
    $n = $sys{'diheds'}{'count'};
    $num = $mol->{'diheds'}{'count'};
    for (my $i=1; $i <= $num; $i++)
    {
        $n++;
        $sys{'diheds'}{'type'}[$n] = $mol->{'diheds'}{'type'}[$i];
        ($a1, $a2, $a3, $a4) = @{$mol->{'diheds'}{'atoms'}[$i]};
        $sys{'diheds'}{'atoms'}[$n] =
            [$a1+$init, $a2+$init, $a3+$init, $a4+$init];
    }

    $sys{'diheds'}{'count'} = $n;

    # Add impropers
    $n = $sys{'improps'}{'count'};
    $num = $mol->{'improps'}{'count'};
    for (my $i=1; $i <= $num; $i++)
    {
        $n++;
        $sys{'improps'}{'type'}[$n] = $mol->{'improps'}{'type'}[$i];
        ($a1, $a2, $a3, $a4) = @{$mol->{'improps'}{'atoms'}[$i]};
        $sys{'improps'}{'atoms'}[$n] =
            [$a1+$init, $a2+$init, $a3+$init, $a4+$init];
    }

    $sys{'improps'}{'count'} = $n;
}

# checkOverlap( %\mol )
# Check overlap for new molecule in system, returning 1 if no overlap
sub checkOverlap
{
    # Variables
    my $mol = $_[0];
    my (@neigh, $cellN, $u, $v, $w, $num, $sep);
    my (@p1, @p2, $t1, $t2, $r1, $r2, $x, $y, $z, $xc, $yc, $zc);

    # Check each atom in new molecule for overlap
    $cellN = $mol->{'neighList'}{'num'};
    $num = $mol->{'atoms'}{'count'};
    for (my $n=1; $n <= $num; $n++)
    {
        ($u, $v, $w) = Polymatic::getCell($mol, $n);
        @p1 = @{$mol->{'atoms'}{'pos'}[$n]};
        $t1 = $mol->{'atoms'}{'type'}[$n];
        $r1 = $mol->{'atomTypes'}{'coeffs'}{'pair'}[$t1][1]/2;
        $r1 = $radius if (!defined($r1) || $r1 == 0);

        # Check only surrounding cells in neighbor list
        for (my $i=$u-1; $i <= $u+1; $i++) {
            $xc = $i - $cellN * POSIX::floor($i/$cellN);
        for (my $j=$v-1; $j <= $v+1; $j++) {
            $yc = $j - $cellN * POSIX::floor($j/$cellN);
        for (my $k=$w-1; $k <= $w+1; $k++) {
            $zc = $k - $cellN * POSIX::floor($k/$cellN);

            @neigh = @{$sys{'neighList'}{'cells'}[$xc][$yc][$zc]};
            for (my $n=0; $n < scalar(@neigh); $n++)
            {
                @p2 = @{$sys{'atoms'}{'pos'}[$neigh[$n]]};
                $t2 = $sys{'atoms'}{'type'}[$neigh[$n]];
                $r2 = $sys{'atomTypes'}{'coeffs'}{'pair'}[$t2][1]/2;
                $r2 = $radius if (!defined($r2) || $r2 == 0);

                $x = $p1[0] - $p2[0];
                $x = $x - $boxL * POSIX::floor($x/$boxL + 0.5);
                $y = $p1[1] - $p2[1];
                $y = $y - $boxL * POSIX::floor($y/$boxL + 0.5);
                $z = $p1[2] - $p2[2];
                $z = $z - $boxL * POSIX::floor($z/$boxL + 0.5);

                # Return 0 if overlap
                $sep = sqrt($x*$x + $y*$y + $z*$z);
                return 0 if ($sep < $scale*($r1+$r2));
            }
        }}}
    }

    # Return 1 i no overlaps
    return 1;
}
