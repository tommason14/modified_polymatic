################################################################################
#
# Polymatic.pm
# This file is part of the Polymatic distribution.
#
# Author: Lauren J. Abbott
# Version: 1.1
# Date: August 16, 2015
#
# Description: Perl module to handle tasks related to defining and using
# molecular systems. Included are subroutines to read in files of different
# types to define a molecular system, write out a molecular system to different
# file types, and perform other actions on molecular systems. To use any of
# these subroutines in other codes, load the module with 'use Polymatic'. Note
# that Polymatic.pm must be in the same directory or at a file path recognized
# by the script (e.g., 'use lib').
#
# Reading files:
#  > readLammps, LAMMPS data file
#  > readTypes, types file, molecular system should be defined
#  > readPdb, PDB file
#  > readGro, GRO file
#  > readXyz, XYZ file
#  > readPsf, PSF file, only reads bonds, atoms should be defined
#  > readXsd, Materials studio XSD file
#  > readPolym, Polymatic input script
#
# Writing files:
#  > writeLammps, LAMMPS data file
#  > writeTypes, types file
#  > writePdb, PDB file
#  > writeGro, GRO file
#  > writeXyz, XYZ file
#  > writePsf, PSF file, meant for use in VMD with coordinate file
#  > writeLmpsTrj, LAMMPS trajectory file, using 'atom style', single frame
#  > writeCar, CAR file for use in Materials Studio
#  > writeMdf, MDF file for use in Materials Studio
#
# Other:
#  > swapUnits, swap between nm and angstroms
#  > getBondType, get bond type for given atoms
#  > getAngleType, get angle type for given atoms
#  > getDihedType, get dihedral type for given atoms
#  > getImpropType, get improper type for given atoms
#  > defineAngles, define angles in system with known bonding
#  > defineDiheds, define dihedrals in system with known bonding
#  > defineImprops, define impropers in system with known bonding
#  > defineMols, define molecules in system with known bonding
#  > getCell, get neighbor list cell for given atom
#  > initNeighList, initialize neighbor list
#  > addNeighList, add atom to neighbor list
#  > getSep, calculate separation between given atoms
#  > unwrapMol, unwrap molecule coordinates so no bonds spand PB
#  > delDupBond, delete duplicate bonds from array
#  > delDupImprops, delete duplicate impropers from array
#  > vectorSub, vector subtraction
#  > dot, vector dot product
#  > norm, vector norm
#  > vectorAng, angle between vectors
#  > normalPlane, normal vector to best fit plane
#  > eqArray, element-by-element array equals
#  > uniqueArray, keep only unique elements in array
#  > group, find value in array
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

package Polymatic;

use strict;
use Math::Trig();
use POSIX();
# use XML::Simple();

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

# chompSplit( $str )
# Chomp string, remove leading/trailing whitespace, and split by spaces
sub chompSplit
{
  my $str = $_[0];
  chomp($str);
  $str =~ s/^\s+|\s+$//g;
  return split(' ', $str);
}

################################################################################
# Reading files

# readLammps( $file )
# Read LAMMPS data file and return molecular system
sub readLammps
{
# Variables
  my $file = $_[0];
  my ($line, $num, $n);
  my (@temp, %sys);

# Open file
  open FILE, "< $file" or die "Error opening file '$file': $!";

# Read in header
  $line = <FILE>;
  @temp = chompSplit($line);
  $sys{'header'} = join(' ', @temp);
  $line = <FILE>;

# Read in remaining lines
  while ($line = <FILE>)
  {
    @temp = chompSplit($line);

# Counts
    if ($temp[1] eq "atoms") {
      $sys{'atoms'}{'count'} = $temp[0];
    } elsif ($temp[1] eq "bonds") {
      $sys{'bonds'}{'count'} = $temp[0];
    } elsif ($temp[1] eq "angles") {
      $sys{'angles'}{'count'} = $temp[0];
    } elsif ($temp[1] eq "dihedrals") {
      $sys{'diheds'}{'count'} = $temp[0];
    } elsif ($temp[1] eq "impropers") {
      $sys{'improps'}{'count'} = $temp[0];
    } elsif ($temp[1]." ".$temp[2] eq "atom types") {
      $sys{'atomTypes'}{'count'} = $temp[0];
    } elsif ($temp[1]." ".$temp[2] eq "bond types") {
      $sys{'bondTypes'}{'count'} = $temp[0];
    } elsif ($temp[1]." ".$temp[2] eq "angle types") {
      $sys{'angleTypes'}{'count'} = $temp[0];
    } elsif ($temp[1]." ".$temp[2] eq "dihedral types") {
      $sys{'dihedTypes'}{'count'} = $temp[0];
    } elsif ($temp[1]." ".$temp[2] eq "improper types") {
      $sys{'impropTypes'}{'count'} = $temp[0];
    }

# Box dimensions
    elsif ($temp[2] eq "xlo") {
      $sys{'boxDims'}{'x'} = [$temp[0], $temp[1], $temp[1]-$temp[0]];
    } elsif ($temp[2] eq "ylo") {
      $sys{'boxDims'}{'y'} = [$temp[0], $temp[1], $temp[1]-$temp[0]];
    } elsif ($temp[2] eq "zlo") {
      $sys{'boxDims'}{'z'} = [$temp[0], $temp[1], $temp[1]-$temp[0]];
    }

# Masses
    elsif ($temp[0] eq "Masses")
    {
      if (defined($sys{'atomTypes'}{'count'})) {
        $num = $sys{'atomTypes'}{'count'};
      } else {
        errExit("Cannot read masses, ".
            "number of atom types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'atomTypes'}{'coeffs'}{'mass'}[$n] = [@temp];
      }
    }

# Pair coeffs
    elsif ($temp[0]." ".$temp[1] eq "Pair Coeffs")
    {
      if (defined($sys{'atomTypes'}{'count'})) {
        $num = $sys{'atomTypes'}{'count'};
      } else {
        errExit("Cannot read pair coefficients, ".
            "number of atom types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'atomTypes'}{'coeffs'}{'pair'}[$n] = [@temp];
      }
    }

# Bond coeffs
    elsif ($temp[0]." ".$temp[1] eq "Bond Coeffs")
    {
      if (defined($sys{'bondTypes'}{'count'})) {
        $num = $sys{'bondTypes'}{'count'};
      } else {
        errExit("Cannot read bond coefficients, ".
            "number of bond types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'bondTypes'}{'coeffs'}{'bond'}[$n] = [@temp];
      }
    }

# Angle coeffs
    elsif ($temp[0]." ".$temp[1] eq "Angle Coeffs")
    {
      if (defined($sys{'angleTypes'}{'count'})) {
        $num = $sys{'angleTypes'}{'count'};
      } else {
        errExit("Cannot read angle coefficients, ".
            "number of angle types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'angleTypes'}{'coeffs'}{'angle'}[$n] = [@temp];
      }
    }

# Bond bond coeffs
    elsif ($temp[0]." ".$temp[1] eq "BondBond Coeffs")
    {
      if (defined($sys{'angleTypes'}{'count'})) {
        $num = $sys{'angleTypes'}{'count'};
      } else {
        errExit("Cannot read bond bond coefficients, ".
            "number of angle types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'angleTypes'}{'coeffs'}{'bondBond'}[$n] = [@temp];
      }
    }

# Bond angle coeffs
    elsif ($temp[0]." ".$temp[1] eq "BondAngle Coeffs")
    {
      if (defined($sys{'angleTypes'}{'count'})) {
        $num = $sys{'angleTypes'}{'count'};
      } else {
        errExit("Cannot read bond angle coefficients, ".
            "number of angle types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'angleTypes'}{'coeffs'}{'bondAng'}[$n] = [@temp];
      }
    }

# Dihedral coeffs
    elsif ($temp[0]." ".$temp[1] eq "Dihedral Coeffs")
    {
      if (defined($sys{'dihedTypes'}{'count'})) {
        $num = $sys{'dihedTypes'}{'count'};
      } else {
        errExit("Cannot read dihedral coefficients, ".
            "number of dihedral types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'dihedTypes'}{'coeffs'}{'dihed'}[$n] = [@temp];
      }
    }

# Middle bond torsion coeffs
    elsif ($temp[0]." ".$temp[1] eq "MiddleBondTorsion Coeffs")
    {
      if (defined($sys{'dihedTypes'}{'count'})) {
        $num = $sys{'dihedTypes'}{'count'};
      } else {
        errExit("Cannot read middle bond torsion coefficients, ".
            "number of dihedral types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'dihedTypes'}{'coeffs'}{'midBondTors'}[$n] = [@temp];
      }
    }

# End bond torsion coeffs
    elsif ($temp[0]." ".$temp[1] eq "EndBondTorsion Coeffs")
    {
      if (defined($sys{'dihedTypes'}{'count'})) {
        $num = $sys{'dihedTypes'}{'count'};
      } else {
        errExit("Cannot read end bond torsion coefficients, ".
            "number of dihedral types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'dihedTypes'}{'coeffs'}{'endBondTors'}[$n] = [@temp];
      }
    }

# Angle torsion coeffs
    elsif ($temp[0]." ".$temp[1] eq "AngleTorsion Coeffs")
    {
      if (defined($sys{'dihedTypes'}{'count'})) {
        $num = $sys{'dihedTypes'}{'count'};
      } else {
        errExit("Cannot read angle torsion coefficients, ".
            "number of dihedral types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'dihedTypes'}{'coeffs'}{'angTors'}[$n] = [@temp];
      }
    }

# Angle angle torsion coeffs
    elsif ($temp[0]." ".$temp[1] eq "AngleAngleTorsion Coeffs")
    {
      if (defined($sys{'dihedTypes'}{'count'})) {
        $num = $sys{'dihedTypes'}{'count'};
      } else {
        errExit("Cannot read angle angle torsion coefficients, ".
            "number of dihedral types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'dihedTypes'}{'coeffs'}{'angAngTors'}[$n] = [@temp];
      }
    }

# Bond bond 1-3 coeffs
    elsif ($temp[0]." ".$temp[1] eq "BondBond13 Coeffs")
    {
      if (defined($sys{'dihedTypes'}{'count'})) {
        $num = $sys{'dihedTypes'}{'count'};
      } else {
        errExit("Cannot read bond bond 1-3 coefficients, ".
            "number of dihedral types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'dihedTypes'}{'coeffs'}{'bondBond13'}[$n] = [@temp];
      }
    }

# Improper coeffs
    elsif ($temp[0]." ".$temp[1] eq "Improper Coeffs")
    {
      if (defined($sys{'impropTypes'}{'count'})) {
        $num = $sys{'impropTypes'}{'count'};
      } else {
        errExit("Cannot read improper coefficients, ".
            "number of improper types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'impropTypes'}{'coeffs'}{'improp'}[$n] = [@temp];
      }
    }

# Angle angle coeffs
    elsif ($temp[0]." ".$temp[1] eq "AngleAngle Coeffs")
    {
      if (defined($sys{'impropTypes'}{'count'})) {
        $num = $sys{'impropTypes'}{'count'};
      } else {
        errExit("Cannot read angle angle coefficients, ".
            "number of improper types is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'impropTypes'}{'coeffs'}{'angAng'}[$n] = [@temp];
      }
    }

# Atoms
    elsif ($temp[0] eq "Atoms")
    {
      if (defined($sys{'atoms'}{'count'})) {
        $num = $sys{'atoms'}{'count'};
      } else {
        errExit("Cannot read atoms, "."
            number of atoms is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'atoms'}{'mol'}[$n] = shift(@temp);
        $sys{'atoms'}{'type'}[$n] = shift(@temp);
        $sys{'atoms'}{'q'}[$n] = shift(@temp);
        $sys{'atoms'}{'pos'}[$n] = [$temp[0], $temp[1], $temp[2]];
        push(@{$sys{'mols'}{'atoms'}[$sys{'atoms'}{'mol'}[$n]]}, $n);
      }

    }

# Bonds
    elsif ($temp[0] eq "Bonds")
    {
      if (defined($sys{'bonds'}{'count'})) {
        $num = $sys{'bonds'}{'count'};
      } else {
        errExit("Cannot read bonds, "."
            number of bonds is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'bonds'}{'type'}[$n] = shift(@temp);
        $sys{'bonds'}{'atoms'}[$n] = [@temp];

        push(@{$sys{'atoms'}{'bonds'}[$temp[0]]}, $n);
        push(@{$sys{'atoms'}{'bonds'}[$temp[1]]}, $n);
        push(@{$sys{'atoms'}{'bonded'}[$temp[0]]}, $temp[1]);
        push(@{$sys{'atoms'}{'bonded'}[$temp[1]]}, $temp[0]);
      }
    }

# Angles
    elsif ($temp[0] eq "Angles")
    {
      if (defined($sys{'angles'}{'count'})) {
        $num = $sys{'angles'}{'count'};
      } else {
        errExit("Cannot read angles, "."
            number of angles is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'angles'}{'type'}[$n] = shift(@temp);
        $sys{'angles'}{'atoms'}[$n] = [@temp];

        push(@{$sys{'atoms'}{'angles'}[$temp[0]]}, $n);
        push(@{$sys{'atoms'}{'angles'}[$temp[1]]}, $n);
        push(@{$sys{'atoms'}{'angles'}[$temp[2]]}, $n);
      }
    }

# Dihedrals
    elsif ($temp[0] eq "Dihedrals")
    {
      if (defined($sys{'diheds'}{'count'})) {
        $num = $sys{'diheds'}{'count'};
      } else {
        errExit("Cannot read dihedrals, "."
            number of dihedrals is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'diheds'}{'type'}[$n] = shift(@temp);
        $sys{'diheds'}{'atoms'}[$n] = [@temp];

        push(@{$sys{'atoms'}{'diheds'}[$temp[0]]}, $n);
        push(@{$sys{'atoms'}{'diheds'}[$temp[1]]}, $n);
        push(@{$sys{'atoms'}{'diheds'}[$temp[2]]}, $n);
        push(@{$sys{'atoms'}{'diheds'}[$temp[3]]}, $n);
      }
    }

# Impropers
    elsif ($temp[0] eq "Impropers")
    {
      if (defined($sys{'improps'}{'count'})) {
        $num = $sys{'improps'}{'count'};
      } else {
        errExit("Cannot read impropers, "."
            number of impropers is not defined.");
      }

      $line = <FILE>;
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        $n = shift(@temp);
        $sys{'improps'}{'type'}[$n] = shift(@temp);
        $sys{'improps'}{'atoms'}[$n] = [@temp];

        push(@{$sys{'atoms'}{'improps'}[$temp[0]]}, $n);
        push(@{$sys{'atoms'}{'improps'}[$temp[1]]}, $n);
        push(@{$sys{'atoms'}{'improps'}[$temp[2]]}, $n);
        push(@{$sys{'atoms'}{'improps'}[$temp[3]]}, $n);
      }
    }
  }

# Close file
  close FILE;

# Number of molecules
  $sys{'mols'}{'count'} = scalar(@{$sys{'mols'}{'atoms'}}) - 1;

# System flags
  $sys{'flags'}{'ang'} = 1;

# Return molecular system
  return %sys;
}

# readTypes( $file, \%sys )
# Read types file and store info in given molecular system
sub readTypes
{
# Variables
  my $file = $_[0];
  my $sys = $_[1];
  my (@temp, $num);

# Open file
  open FILE, "< $file" or die "Error opening file '$file': $!";

# Read in each line
  while(my $line = <FILE>)
  {
    @temp = chompSplit($line);
    next if (substr($line,0,1) eq "#");

# Atom types
    if ($temp[0] eq "atom")
    {
      if (defined($sys->{'atomTypes'}{'count'})) {
        $num = $sys->{'atomTypes'}{'count'};
      } else {
        errExit("Number of atom types is not defined, ".
            "cannot read in atom types section.");
      }

      for (my $i=1; $i <= $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        if (substr($line,0,1) eq "#") {
          $i--;
          next;
        }

        $sys->{'atomTypes'}{'name'}[$temp[0]] = $temp[1];
        $sys->{'atomTypes'}{'num'}{$temp[1]} = $temp[0];
      }
    }

# Bond types
    elsif ($temp[0] eq "bond")
    {
      if (defined($sys->{'bondTypes'}{'count'})) {
        $num = $sys->{'bondTypes'}{'count'};
      } else {
        errExit("Number of bond types is not defined, ".
            "cannot read in bond types section.");
      }

      for (my $i=1; $i <= $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        if (substr($line,0,1) eq "#") {
          $i--;
          next;
        }

        $sys->{'bondTypes'}{'name'}[$temp[0]] = $temp[1];
        $sys->{'bondTypes'}{'num'}{$temp[1]} = $temp[0];
      }
    }

# Angle types
    elsif ($temp[0] eq "angle")
    {
      if (defined($sys->{'angleTypes'}{'count'})) {
        $num = $sys->{'angleTypes'}{'count'};
      } else {
        errExit("Number of angle types is not defined, ".
            "cannot read in angle types section.");
      }

      for (my $i=1; $i <= $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        if (substr($line,0,1) eq "#") {
          $i--;
          next;
        }

        $sys->{'angleTypes'}{'name'}[$temp[0]] = $temp[1];
        $sys->{'angleTypes'}{'num'}{$temp[1]} = $temp[0];
      }
    }

# Dihedral types
    elsif ($temp[0] eq "dihedral")
    {
      if (defined($sys->{'dihedTypes'}{'count'})) {
        $num = $sys->{'dihedTypes'}{'count'};
      } else {
        errExit("Number of dihedral types is not defined, ".
            "cannot read in dihedral types section.");
      }

      for (my $i=1; $i <= $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        if (substr($line,0,1) eq "#") {
          $i--;
          next;
        }

        $sys->{'dihedTypes'}{'name'}[$temp[0]] = $temp[1];
        $sys->{'dihedTypes'}{'num'}{$temp[1]} = $temp[0];
      }
    }

# Improper types
    elsif ($temp[0] eq "improper")
    {
      if (defined($sys->{'impropTypes'}{'count'})) {
        $num = $sys->{'impropTypes'}{'count'};
      } else {
        errExit("Number of improper types is not defined, ".
            "cannot read in improper types section.");
      }

      for (my $i=1; $i <= $num; $i++)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        if (substr($line,0,1) eq "#") {
          $i--;
          next;
        }

        $sys->{'impropTypes'}{'name'}[$temp[0]] = $temp[1];
        $sys->{'impropTypes'}{'num'}{$temp[1]} = $temp[0];
      }
    }
  }

# Close file
  close FILE
}

# readPsf( $file, \%sys )
# Read PSF file and apply bonding to given molecular system
sub readPsf
{
# Variables
  my $file = $_[0];
  my $sys = $_[1];
  my ($i, $a1, $a2, $type, $t1, $t2, $str);
  my ($numBonds, $numBondTypes);
  my (@temp);
  my $at = 'atomTypes';

# Open file
  open FILE, "< $file" or die "Error opening file '$file': $!";

# Read lines
  while (my $line = <FILE>)
  {
    @temp = chompSplit($line);

# Bonds
    if ($temp[1] eq '!NBOND:' && $temp[2] eq 'bonds')
    {
      $i = 0;
      $numBonds = $temp[0];
      while ($i < $numBonds)
      {
        $line = <FILE>;
        @temp = chompSplit($line);

        while(scalar(@temp) > 0)
        {
          $a1 = shift(@temp);
          $a2 = shift(@temp);
          $i++;

          $type = getBondType($sys, [$a1, $a2]);
          if ($type == 0)
          {
            $numBondTypes++;
            $t1 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$a1]];
            $t2 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$a2]];
            $str = $t1.','.$t2;
            $sys->{'bondTypes'}{'name'}[$numBondTypes] = $str;
            $sys->{'bondTypes'}{'num'}{$str} = $numBondTypes;
            $sys->{'bonds'}{'type'}[$i] = $numBondTypes;
            $sys->{'bonds'}{'atoms'}[$i] = [$a1, $a2];
          }
          elsif ($type < 0)
          {
            $sys->{'bonds'}{'type'}[$i] = -1*$type;
            $sys->{'bonds'}{'atoms'}[$i] = [$a2, $a1];
          }
          else
          {
            $sys->{'bonds'}{'type'}[$i] = $type;
            $sys->{'bonds'}{'atoms'}[$i] = [$a1, $a2];
          }

          push(@{$sys->{'atoms'}{'bonds'}[$a1]}, $i);
          push(@{$sys->{'atoms'}{'bonds'}[$a2]}, $i);
          push(@{$sys->{'atoms'}{'bonded'}[$a1]}, $a2);
          push(@{$sys->{'atoms'}{'bonded'}[$a2]}, $a1);
        }
      }
    }
  }

# System counts
  $sys->{'bonds'}{'count'} = $numBonds;
  $sys->{'bondTypes'}{'count'} = $numBondTypes;

# Close file
  close FILE;
}

# readXsd( $file )
# Read Materials Studio XSD file and return molecular system
# sub readXsd
# {
# # Variables
#   my $file = $_[0];
#   my ($lenX, $lenY, $lenZ, $id, $type, $str, $x, $y, $z, $a1, $a2);
#   my ($numAtoms, $numAtomTypes, $numBonds, $numBondTypes);
#   my (%atom, %atomImage, %atomKey, %bond, %bondImage, %bondKey);
#   my (%sys, %xml, @temp, @bonds);
#
# # Parse XML file
#   %xml = %{XML::Simple::XMLin(join('', ("<system>", `awk '/Atom3d/' $file`,
#           `awk '/Bond/' $file`, `awk '/SpaceGroup/' $file`, "</system>")))};
#
# # Box dimensions
#   if (defined($xml{'SpaceGroup'}))
#   {
#     @temp = split(',', $xml{'SpaceGroup'}{'AVector'});
#     $lenX = $temp[0];
#     $sys{'boxDims'}{'x'} = [0, $lenX, $lenX];
#
#     @temp = split(',', $xml{'SpaceGroup'}{'BVector'});
#     $lenY = $temp[1];
#     $sys{'boxDims'}{'y'} = [0, $lenY, $lenY];
#
#     @temp = split(',', $xml{'SpaceGroup'}{'CVector'});
#     $lenZ = $temp[2];
#     $sys{'boxDims'}{'z'} = [0, $lenZ, $lenZ];
#   }
#
# # Atoms
#   for (my $i=0; $i < scalar(@{$xml{'Atom3d'}}); $i++)
#   {
#     %atom = %{$xml{'Atom3d'}[$i]};
#
#     if (defined($atom{'ID'}))
#     {
#       $id = $atom{'ID'};
#     }
#     else
#     {
#       errExit("Atoms in XSD file are not properly defined.");
#     }
#
#     if (defined($atom{'ImageOf'}))
#     {
#       $atomImage{$id} = $atom{'ImageOf'};
#     }
#     else
#     {
#       $numAtoms++;
#       $atomKey{$id} = $numAtoms;
#
#       if (defined($atom{'ForcefieldType'}))
#       {
#         $type = $atom{'ForcefieldType'};
#         if (defined($sys{'atomTypes'}{'num'}{$type}))
#         {
#           $sys{'atoms'}{'type'}[$numAtoms] =
#             $sys{'atomTypes'}{'num'}{$type};
#         }
#         else
#         {
#           $numAtomTypes++;
#           $sys{'atomTypes'}{'name'}[$numAtomTypes] = $type;
#           $sys{'atomTypes'}{'num'}{$type} = $numAtomTypes;
#           $sys{'atoms'}{'type'}[$numAtoms] = $numAtomTypes;
#         }
#       }
#       else
#       {
#         errExit("Type for atom $id is not defined.");
#       }
#
#       $sys{'atoms'}{'q'}[$numAtoms] = $atom{'Charge'};
#       $sys{'atoms'}{'bonds'}[$numAtoms] =
#         [split(',', $atom{'Connections'})];
#
#       ($x, $y, $z) = split(',', $atom{'XYZ'});
#       if (defined($xml{'SpaceGroup'}))
#       {
#         $sys{'atoms'}{'pos'}[$numAtoms] =
#           [$x*$lenX, $y*$lenY, $z*$lenZ];
#       }
#       else
#       {
#         $sys{'atoms'}{'pos'}[$numAtoms] = [$x, $y, $z];
#       }
#     }
#   }
#
#   $sys{'atoms'}{'count'} = $numAtoms;
#   $sys{'atomTypes'}{'count'} = $numAtomTypes;
#
# # Bonds
#   for (my $i=0; $i < scalar(@{$xml{'Bond'}}); $i++)
#   {
#     %bond = %{$xml{'Bond'}[$i]};
#
#     if (defined($bond{'ID'}))
#     {
#       $id = $bond{'ID'};
#     }
#     else
#     {
#       errExit("Bonds in XSD file are not properly defined.");
#     }
#
#     if (defined($bond{'ImageOf'}))
#     {
#       $bondImage{$id} = $bond{'ImageOf'};
#     }
#     else
#     {
#       $numBonds++;
#       $bondKey{$id} = $numBonds;
#
#       ($a1, $a2) = split(',', $bond{'Connects'});
#       $a1 = $atomImage{$a1} while (defined($atomImage{$a1}));
#       $a2 = $atomImage{$a2} while (defined($atomImage{$a2}));
#       errExit("Atoms in bond '$id' are not properly defined.")
#         if (!defined($a1) || !defined($a2));
#
#       if (defined($atomKey{$a1}) && defined($atomKey{$a2}))
#       {
#         $a1 = $atomKey{$a1};
#         $a2 = $atomKey{$a2};
#       }
#       else
#       {
#         errExit("Atoms in bond '$id' are not properly defined.");
#       }
#
#       $type = getBondType(\%sys, [$a1, $a2]);
#       if ($type == 0)
#       {
#         $numBondTypes++;
#         $str = $sys{'atomTypes'}{'name'}[$sys{'atoms'}{'type'}[$a1]].
#           ','.$sys{'atomTypes'}{'name'}[$sys{'atoms'}{'type'}[$a2]];
#         $sys{'bondTypes'}{'name'}[$numBondTypes] = $str;
#         $sys{'bondTypes'}{'num'}{$str} = $numBondTypes;
#         $sys{'bonds'}{'type'}[$numBonds] = $numBondTypes;
#         $sys{'bonds'}{'atoms'}[$numBonds] = [$a1, $a2];
#       }
#       elsif ($type < 0)
#       {
#         $sys{'bonds'}{'type'}[$numBonds] = -1*$type;
#         $sys{'bonds'}{'atoms'}[$numBonds] = [$a2, $a1];
#       }
#       else
#       {
#         $sys{'bonds'}{'type'}[$numBonds] = $type;
#         $sys{'bonds'}{'atoms'}[$numBonds] = [$a1, $a2];
#       }
#     }
#   }
#
#   $sys{'bonds'}{'count'} = $numBonds;
#   $sys{'bondTypes'}{'count'} = $numBondTypes;
#
# # Atom bonds
#   for (my $i=1; $i <= $numAtoms; $i++)
#   {
#     @bonds = @{$sys{'atoms'}{'bonds'}[$i]};
#     for (my $j=0; $j < scalar(@bonds); $j++)
#     {
#       $id = $bonds[$j];
#       $id = $bondImage{$id} while (defined($bondImage{$id}));
#       $bonds[$j] = $bondKey{$id};
#       ($a1, $a2) = @{$sys{'bonds'}{'atoms'}[$bondKey{$id}]};
#
#       if ($a1 == $i)
#       {
#         push(@{$sys{'atoms'}{'bonded'}[$i]}, $a2);
#       }
#       elsif ($a2 == $i)
#       {
#         push(@{$sys{'atoms'}{'bonded'}[$i]}, $a1);
#       }
#       else
#       {
#         errExit("Atom bonds are not properly defined.");
#       }
#     }
#
#     $sys{'atoms'}{'bonds'}[$i] = [@bonds];
#   }
#
# # System flags
#   $sys{'flags'}{'ang'} = 1;
#
# # Return system
#   return %sys;
# }

# readPdb( $file )
# Read PDB file and return molecular system
sub readPdb
{
# Variables
  my $file = $_[0];
  my ($cmd, $str, $num, $id, $x, $y, $z, $a1, $a2, $type);
  my ($numAtoms, $numAtomTypes, $numBonds, $numBondTypes);
  my (%sys);

# Open file
  open FILE, "< $file" or die "Error opening file '$file': $!";

# Read by line
  while (my $line = <FILE>)
  {
    chomp($line);
    $cmd = substr($line,0,6);

# Header
    if ($cmd eq "TITLE ")
    {
      $str = substr($line, 10);
      $str =~ s/^\s+|\s+$//g;
      $sys{'header'} = $str;
    }

# Box size
    elsif ($cmd eq "CRYST1")
    {
      $x = substr($line, 6, 9);
      $x =~ s/^\s+|\s+$//g;
      $sys{'boxDims'}{'x'} = [0, $x, $x];

      $y = substr($line, 15, 9);
      $y =~ s/^\s+|\s+$//g;
      $sys{'boxDims'}{'y'} = [0, $y, $y];

      $z = substr($line, 24, 9);
      $z =~ s/^\s+|\s+$//g;
      $sys{'boxDims'}{'z'} = [0, $z, $z];
    }

# Atom
    elsif ($cmd eq "ATOM  ")
    {
      $numAtoms++;
      $id = substr($line, 6, 5);
      $id =~ s/^\s+|\s+$//g;
      $sys{'atoms'}{'mol'}[$id] = 1;

      $str = substr($line, 12, 4);
      $str =~ s/^\s+|\s+$//g;
      $sys{'atoms'}{'name'}[$id] = $str;

      if (defined($sys{'atomTypes'}{'num'}{$str}))
      {
        $sys{'atoms'}{'type'}[$id] = $sys{'atomTypes'}{'num'}{$str};
      }
      else
      {
        $numAtomTypes++;
        $sys{'atomTypes'}{'num'}{$str} = $numAtomTypes;
        $sys{'atomTypes'}{'name'}[$numAtomTypes] = $str;
        $sys{'atoms'}{'type'}[$id] = $numAtomTypes;
      }

      $str = substr($line, 17, 3);
      $str =~ s/^\s+|\s+$//g;
      $sys{'atoms'}{'res'}[$id] = $str;

      $num = substr($line, 22, 4);
      $num =~ s/^\s+|\s+$//g;
      $sys{'atoms'}{'resNum'}[$id] = $num;

      $x = substr($line, 30, 8);
      $x =~ s/^\s+|\s+$//g;
      $y = substr($line, 38, 8);
      $y =~ s/^\s+|\s+$//g;
      $z = substr($line, 46, 8);
      $z =~ s/^\s+|\s+$//g;
      $sys{'atoms'}{'pos'}[$id] = [$x, $y, $z];
    }

# Bonds
    elsif ($cmd eq "CONECT")
    {
      $a1 = substr($line, 6, 5);
      $a1 =~ s/^\s+|\s+$//g;

      for (my $i=11; $i <= 26; $i+=5)
      {
        $a2 = substr($line, $i, 5);
        $a2 =~ s/^\s+|\s+$//g;

        next if ($a2 eq '');
        $numBonds++;

        $type = getBondType(\%sys, [$a1, $a2]);
        if ($type == 0)
        {
          $numBondTypes++;
          $str = $sys{'atomTypes'}{'name'}[$sys{'atoms'}{'type'}[$a1]].
            ','.$sys{'atomTypes'}{'name'}[$sys{'atoms'}{'type'}[$a2]];
          $sys{'bondTypes'}{'name'}[$numBondTypes] = $str;
          $sys{'bondTypes'}{'num'}{$str} = $numBondTypes;
          $sys{'bonds'}{'type'}[$numBonds] = $numBondTypes;
          $sys{'bonds'}{'atoms'}[$numBonds] = [$a1, $a2];
        }
        elsif ($type < 0)
        {
          $sys{'bonds'}{'type'}[$numBonds] = -1*$type;
          $sys{'bonds'}{'atoms'}[$numBonds] = [$a2, $a1];
        }
        else
        {
          $sys{'bonds'}{'type'}[$numBonds] = $type;
          $sys{'bonds'}{'atoms'}[$numBonds] = [$a1, $a2];
        }

        push(@{$sys{'atoms'}{'bonds'}[$a1]}, $numBonds);
        push(@{$sys{'atoms'}{'bonds'}[$a2]}, $numBonds);
        push(@{$sys{'atoms'}{'bonded'}[$a1]}, $a2);
        push(@{$sys{'atoms'}{'bonded'}[$a2]}, $a1);
      }
    }
  }

# System setup
  $sys{'atoms'}{'count'} = $numAtoms;
  $sys{'atomTypes'}{'count'} = $numAtomTypes;
  $sys{'bonds'}{'count'} = $numBonds;
  $sys{'bondTypes'}{'count'} = $numBondTypes;
  $sys{'flags'}{'ang'} = 1;

# Close file
  close FILE;

# Return system
  return %sys;
}

# readGro( $file )
# Read GRO file and return molecular system
sub readGro
{
# Variables
  my $file = $_[0];
  my ($line, $id, $num, $str, $x, $y, $z);
  my ($numAtoms, $numAtomTypes);
  my (%sys, @temp);

# Open file
  open FILE, "< $file" or die "Error opening file '$file': $!";

# Header
  $line = <FILE>;
  chomp($line);
  $line =~ s/^\s+|\s+$//g;
  $sys{'header'} = $line;

# Number of atoms
  $line = <FILE>;
  chomp($line);
  $line =~ s/^\s+|\s+$//g;
  $numAtoms = $line;
  $sys{'atoms'}{'count'} = $numAtoms;

# Atoms
  for (my $i=0; $i < $numAtoms; $i++)
  {
    $line = <FILE>;
    chomp($line);

    $id = substr($line, 15, 5);
    $id =~ s/^\s+|\s+$//g;

    $num = substr($line, 0, 5);
    $num =~ s/^\s+|\s+$//g;
    $sys{'atoms'}{'resNum'}[$id] = $num;

    $str = substr($line, 5, 5);
    $str =~ s/^\s+|\s+$//g;
    $sys{'atoms'}{'resName'}[$id] = $str;

    $str = substr($line, 10, 5);
    $str =~ s/^\s+|\s+$//g;
    $sys{'atoms'}{'name'}[$id] = $str;

    if (defined($sys{'atomTypes'}{'num'}{$str}))
    {
      $sys{'atoms'}{'type'}[$id] = $sys{'atomTypes'}{'num'}{$str};
    }
    else
    {
      $numAtomTypes++;
      $sys{'atomTypes'}{'num'}{$str} = $numAtomTypes;
      $sys{'atomTypes'}{'name'}[$numAtomTypes] = $str;
      $sys{'atoms'}{'type'}[$id] = $numAtomTypes;
    }

    $x = substr($line, 20, 8);
    $x =~ s/^\s+|\s+$//g;
    $y = substr($line, 28, 8);
    $y =~ s/^\s+|\s+$//g;
    $z = substr($line, 36, 8);
    $z =~ s/^\s+|\s+$//g;
    $sys{'atoms'}{'pos'}[$id] = [$x, $y, $z];
  }

# Box size
  $line = <FILE>;
  @temp = chompSplit($line);
  $sys{'boxDims'}{'x'} = [0, $temp[0], $temp[0]];
  $sys{'boxDims'}{'y'} = [0, $temp[1], $temp[1]];
  $sys{'boxDims'}{'z'} = [0, $temp[2], $temp[2]];

# System flags and counts
  $sys{'flags'}{'ang'} = 0;
  $sys{'atomTypes'}{'count'} = $numAtomTypes;

# Close file
  close FILE;

# Return system
  return %sys;
}

# readXyz( $file )
# Read XYZ file and return molecular system
sub readXyz
{
# Variables
  my $file = $_[0];
  my ($line, $numAtoms, $str, $x, $y, $z);
  my (%sys);

# Open file
  open FILE, "< $file" or die "Error opening file '$file': $!";

# Number of atoms
  $line = <FILE>;
  chomp($line);
  $line =~ s/^\s+|\s+$//g;
  $numAtoms = $line;
  $sys{'atoms'}{'count'} = $numAtoms;

# Header
  $line = <FILE>;
  chomp($line);
  $line =~ s/^\s+|\s+$//g;
  $sys{'header'} = $line;

# Atoms
  for (my $i=1; $i <= $numAtoms; $i++)
  {
    $line = <FILE>;
    ($str, $x, $y, $z) = chompSplit($line);
    $sys{'atoms'}{'type'}[$i] = $str;
    $sys{'atoms'}{'pos'}[$i] = [$x, $y, $z];
  }

# System flags
  $sys{'flags'}{'ang'} = 1;

# Close file
  close FILE;

# Return system
  return %sys;
}

# readPolym( $file )
# Read Polymatic input script
sub readPolym
{
# Variables
  my $file = $_[0];
  my (%inp, $command, $num, @params, @temp1, @temp2, @temp3);

# Open file
  open FILE, "< $file" or die "Error opening file '$file': $!";

# Read lines
  while (my $line = <FILE>)
  {
    @params = chompSplit($line);
    $command = shift(@params);

# Link
    if ($command eq "link")
    {
      @temp1 = split(',', $params[0]);
      @temp2 = split(',', $params[1]);
      errExit("The 'link' command is not properly defined.")
        if (scalar(@temp1) != 2 || scalar(@temp2) != 2);
      $inp{'link'} = [[@temp1], [@temp2]];
    }

# Cutoff
    elsif ($command eq "cutoff")
    {
      $inp{'cutoff'} = $params[0];
      errExit("The bonding radius must be greater than 0.")
        if ($inp{'cutoff'} <= 0);
    }

# Charge
    elsif ($command eq "charge")
    {
      $inp{'charge'} = [$params[0], $params[1]];
    }

# Intra
    elsif ($command eq "intra")
    {
      $inp{'intra'} = $params[1] if ($params[0] eq 'true');
      errExit("The intramolecular check must be a value 1 to 5.")
        if ($inp{'intra'} != 1 && $inp{'intra'} != 2 &&
            $inp{'intra'} != 3 && $inp{'intra'} != 4 &&
            $inp{'intra'} != 5);
    }

# Bond
    elsif ($command eq "bond")
    {
      errExit("The 'bond' command is not properly defined.")
        if (scalar(@params) != 2);
      push(@{$inp{'bond'}}, [@params]);
    }

# Align
    elsif ($command eq "align" ||
        $command eq "vector" ||
        $command eq "plane")
    {
      warning("The 'vector' command is deprecated. ".
          "Please use the 'align' command instead.")
        if ($command eq "vector");

      warning("The 'plane' command is deprecated. ".
          "Please use the 'align' command instead.")
        if ($command eq "plane");

      @temp1 = split(',', $params[0]);
      @temp2 = split(',', $params[1]);
      push(@{$inp{'align'}{'atoms'}}, [[@temp1], [@temp2]]);
      errExit("The vector/plane definitions in the 'align' command ".
          "must have more than 2 atoms.")
        if (scalar(@temp1) < 2 || scalar(@temp2) < 2);

      @temp3 = split(',', $params[2]);
      push(@{$inp{'align'}{'cond'}}, [@temp3]);
      errExit("The conditions in the 'align' command are ".
          "not properly defined.")
        if (scalar(@temp3) != 1 && scalar(@temp3) != 3);
    }

# Connect
    elsif ($command eq "connect")
    {
      $num = $params[0];
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @params = chompSplit($line);
        $inp{'connect'}[$params[0]] = [split(',', $params[1])];
      }
    }

# Types
    elsif ($command eq "types")
    {
      $num = $params[0];
      for (my $i=0; $i < $num; $i++)
      {
        $line = <FILE>;
        @params = chompSplit($line);
        $inp{'types'}[$params[0]] = $params[1];
      }
    }

# Error
    else
    {
      errExit("Input script command '$command' not recognized.");
    }
  }

# Close file
  close FILE;

# Check for required
  errExit("Cutoff radius is not defined.")
    if (!defined($inp{'cutoff'}));
  errExit("Linking atoms are not defined.")
    if (!defined($inp{'link'}[0][0]) || !defined($inp{'link'}[0][1]) ||
        !defined($inp{'link'}[1][0]) || !defined($inp{'link'}[1][1]));

# Return input settings
  return %inp;
}

################################################################################
# Writing files

# writeLammps( $file, \%sys )
# Write LAMMPS data file for given molecular system
sub writeLammps
{
# Variables
  my $file = $_[0];
  my $sys = $_[1];
  my ($num, $mol, $type, $q, @temp);
  my $flag = 0;

# Adjust distance units if needed
  if ($sys->{'flags'}{'ang'} == 0) {
    swapUnits($sys);
  }

# Open file
  open FILE, "> $file" or die "Error opening file '$file': $!";

# Header
  printf FILE "%s\n\n", $sys->{'header'};

# Counts
  printf FILE "%d atoms\n", $sys->{'atoms'}{'count'};
  printf FILE "%d bonds\n", $sys->{'bonds'}{'count'};
  printf FILE "%d angles\n", $sys->{'angles'}{'count'};
  printf FILE "%d dihedrals\n", $sys->{'diheds'}{'count'};
  printf FILE "%d impropers\n\n", $sys->{'improps'}{'count'};
  printf FILE "%d atom types\n", $sys->{'atomTypes'}{'count'}
  if ($sys->{'atomTypes'}{'count'} > 0);
  printf FILE "%d bond types\n", $sys->{'bondTypes'}{'count'}
  if ($sys->{'bondTypes'}{'count'} > 0);
  printf FILE "%d angle types\n", $sys->{'angleTypes'}{'count'}
  if ($sys->{'angleTypes'}{'count'} > 0);
  printf FILE "%d dihedral types\n", $sys->{'dihedTypes'}{'count'}
  if ($sys->{'dihedTypes'}{'count'} > 0);
  printf FILE "%d improper types\n", $sys->{'impropTypes'}{'count'}
  if ($sys->{'impropTypes'}{'count'} > 0);
  printf FILE "\n";

# Box dimensions
  printf FILE "%f %f xlo xhi\n", $sys->{'boxDims'}{'x'}[0],
  $sys->{'boxDims'}{'x'}[1];
  printf FILE "%f %f ylo yhi\n", $sys->{'boxDims'}{'y'}[0],
  $sys->{'boxDims'}{'y'}[1];
  printf FILE "%f %f zlo zhi\n\n", $sys->{'boxDims'}{'z'}[0],
  $sys->{'boxDims'}{'z'}[1];

# Masses
  if (defined($sys->{'atomTypes'}{'count'}) &&
      defined($sys->{'atomTypes'}{'coeffs'}{'mass'}))
  {
    $num = $sys->{'atomTypes'}{'count'};
    printf FILE "Masses\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'atomTypes'}{'coeffs'}{'mass'}[$i])) {
        @temp = @{$sys->{'atomTypes'}{'coeffs'}{'mass'}[$i]};
      } else {
        errExit("Mass for atom type $i is not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Pair coeffs
  if (defined($sys->{'atomTypes'}{'count'}) &&
      defined($sys->{'atomTypes'}{'coeffs'}{'pair'}))
  {
    $num = $sys->{'atomTypes'}{'count'};
    printf FILE "Pair Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'atomTypes'}{'coeffs'}{'pair'}[$i])) {
        @temp = @{$sys->{'atomTypes'}{'coeffs'}{'pair'}[$i]};
      } else {
        errExit("Pair coefficients for atom type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Bond coeffs
  if (defined($sys->{'bondTypes'}{'count'}) &&
      defined($sys->{'bondTypes'}{'coeffs'}{'bond'}))
  {
    $num = $sys->{'bondTypes'}{'count'};
    printf FILE "Bond Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'bondTypes'}{'coeffs'}{'bond'}[$i])) {
        @temp = @{$sys->{'bondTypes'}{'coeffs'}{'bond'}[$i]};
      } else {
        errExit("Bond coefficients for bond type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Angle coeffs
  if (defined($sys->{'angleTypes'}{'count'}) &&
      defined($sys->{'angleTypes'}{'coeffs'}{'angle'}))
  {
    $num = $sys->{'angleTypes'}{'count'};
    printf FILE "Angle Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'angleTypes'}{'coeffs'}{'angle'}[$i])) {
        @temp = @{$sys->{'angleTypes'}{'coeffs'}{'angle'}[$i]};
      } else {
        errExit("Angle coefficients for ".
            "angle type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Bond bond coeffs
  if (defined($sys->{'angleTypes'}{'count'}) &&
      defined($sys->{'angleTypes'}{'coeffs'}{'bondBond'}))
  {
    $num = $sys->{'angleTypes'}{'count'};
    printf FILE "BondBond Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'angleTypes'}{'coeffs'}{'bondBond'}[$i])) {
        @temp = @{$sys->{'angleTypes'}{'coeffs'}{'bondBond'}[$i]};
      } else {
        errExit("Bond bond coefficients for ".
            "angle type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Bond angle coeffs
  if (defined($sys->{'angleTypes'}{'count'}) &&
      defined($sys->{'angleTypes'}{'coeffs'}{'bondAng'}))
  {
    $num = $sys->{'angleTypes'}{'count'};
    printf FILE "BondAngle Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'angleTypes'}{'coeffs'}{'bondAng'}[$i])) {
        @temp = @{$sys->{'angleTypes'}{'coeffs'}{'bondAng'}[$i]};
      } else {
        errExit("Bond angle coefficients for ".
            "angle type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Dihedral coeffs
  if (defined($sys->{'dihedTypes'}{'count'}) &&
      defined($sys->{'dihedTypes'}{'coeffs'}{'dihed'}))
  {
    $num = $sys->{'dihedTypes'}{'count'};
    printf FILE "Dihedral Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'dihedTypes'}{'coeffs'}{'dihed'}[$i])) {
        @temp = @{$sys->{'dihedTypes'}{'coeffs'}{'dihed'}[$i]};
      } else {
        errExit("Dihedral coefficients for ".
            "dihedral type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Middle bond torsion coeffs
  if (defined($sys->{'dihedTypes'}{'count'}) &&
      defined($sys->{'dihedTypes'}{'coeffs'}{'midBondTors'}))
  {
    $num = $sys->{'dihedTypes'}{'count'};
    printf FILE "MiddleBondTorsion Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'dihedTypes'}{'coeffs'}{'midBondTors'}[$i])) {
        @temp = @{$sys->{'dihedTypes'}{'coeffs'}{'midBondTors'}[$i]};
      } else {
        errExit("Middle bond torsion coefficients for ".
            "dihedral type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# End bond torsion coeffs
  if (defined($sys->{'dihedTypes'}{'count'}) &&
      defined($sys->{'dihedTypes'}{'coeffs'}{'endBondTors'}))
  {
    $num = $sys->{'dihedTypes'}{'count'};
    printf FILE "EndBondTorsion Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'dihedTypes'}{'coeffs'}{'endBondTors'}[$i])) {
        @temp = @{$sys->{'dihedTypes'}{'coeffs'}{'endBondTors'}[$i]};
      } else {
        errExit("End bond torsion coefficients for ".
            "dihedral type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Angle torsion coeffs
  if (defined($sys->{'dihedTypes'}{'count'}) &&
      defined($sys->{'dihedTypes'}{'coeffs'}{'angTors'}))
  {
    $num = $sys->{'dihedTypes'}{'count'};
    printf FILE "AngleTorsion Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'dihedTypes'}{'coeffs'}{'angTors'}[$i])) {
        @temp = @{$sys->{'dihedTypes'}{'coeffs'}{'angTors'}[$i]};
      } else {
        errExit("Angle torsion coefficients for ".
            "dihedral type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Angle angle torsion coeffs
  if (defined($sys->{'dihedTypes'}{'count'}) &&
      defined($sys->{'dihedTypes'}{'coeffs'}{'angAngTors'}))
  {
    $num = $sys->{'dihedTypes'}{'count'};
    printf FILE "AngleAngleTorsion Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'dihedTypes'}{'coeffs'}{'angAngTors'}[$i])) {
        @temp = @{$sys->{'dihedTypes'}{'coeffs'}{'angAngTors'}[$i]};
      } else {
        errExit("Angle angle torsion coefficients for ".
            "dihedral type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Bond bond 1-3 coeffs
  if (defined($sys->{'dihedTypes'}{'count'}) &&
      defined($sys->{'dihedTypes'}{'coeffs'}{'bondBond13'}))
  {
    $num = $sys->{'dihedTypes'}{'count'};
    printf FILE "BondBond13 Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'dihedTypes'}{'coeffs'}{'bondBond13'}[$i])) {
        @temp = @{$sys->{'dihedTypes'}{'coeffs'}{'bondBond13'}[$i]};
      } else {
        errExit("Bond bond 1-3 coefficients for ".
            "dihedral type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Improper coeffs
  if (defined($sys->{'impropTypes'}{'count'}) &&
      defined($sys->{'impropTypes'}{'coeffs'}{'improp'}))
  {
    $num = $sys->{'impropTypes'}{'count'};
    printf FILE "Improper Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'impropTypes'}{'coeffs'}{'improp'}[$i])) {
        @temp = @{$sys->{'impropTypes'}{'coeffs'}{'improp'}[$i]};
      } else {
        errExit("Improper coefficients for ".
            "improper type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Angle angle coeffs
  if (defined($sys->{'impropTypes'}{'count'}) &&
      defined($sys->{'impropTypes'}{'coeffs'}{'angAng'}))
  {
    $num = $sys->{'impropTypes'}{'count'};
    printf FILE "AngleAngle Coeffs\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'impropTypes'}{'coeffs'}{'angAng'}[$i])) {
        @temp = @{$sys->{'impropTypes'}{'coeffs'}{'angAng'}[$i]};
      } else {
        errExit("Angle angle coefficients for ".
            "improper type $i are not defined.");
      }

      printf FILE " $i @temp\n";
    }
    printf FILE "\n";
  }

# Atoms
  if (defined($sys->{'atoms'}{'count'}) &&
      defined($sys->{'atoms'}{'type'}) &&
      defined($sys->{'atoms'}{'pos'}))
  {
    $num = $sys->{'atoms'}{'count'};
    printf FILE "Atoms\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'atoms'}{'type'}[$i]) &&
          defined($sys->{'atoms'}{'pos'}[$i]))
      {
        $mol = $sys->{'atoms'}{'mol'}[$i];
        $type = $sys->{'atoms'}{'type'}[$i];
        $q = $sys->{'atoms'}{'q'}[$i];
        @temp = @{$sys->{'atoms'}{'pos'}[$i]};
      }
      else
      {
        errExit("Data for atom $i is not defined.");
      }

      if (defined($sys->{'atoms'}{'mol'}[$i]))
      {
        $mol = $sys->{'atoms'}{'mol'}[$i];
      }
      else
      {
        $mol = 1;
        $flag = 1;
      }

      if (defined($sys->{'atoms'}{'q'}[$i]))
      {
        $q = $sys->{'atoms'}{'q'}[$i];
      }
      else
      {
        $q = 0.0;
        $flag = 1;
      }

      printf FILE " $i $mol $type $q @temp\n";
    }
    printf FILE "\n";
  }

# Bonds
  if (defined($sys->{'bonds'}{'count'}) &&
      defined($sys->{'bonds'}{'type'}) &&
      defined($sys->{'bonds'}{'atoms'}))
  {
    $num = $sys->{'bonds'}{'count'};
    printf FILE "Bonds\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'bonds'}{'type'}[$i]) &&
          defined($sys->{'bonds'}{'atoms'}[$i]))
      {
        $type = $sys->{'bonds'}{'type'}[$i];
        @temp = @{$sys->{'bonds'}{'atoms'}[$i]};
      }
      else
      {
        errExit("Data for bond $i is not defined.");
      }

      printf FILE " $i $type @temp\n";
    }
    printf FILE "\n";
  }

# Angles
  if (defined($sys->{'angles'}{'count'}) &&
      defined($sys->{'angles'}{'type'}) &&
      defined($sys->{'angles'}{'atoms'}))
  {
    $num = $sys->{'angles'}{'count'};
    printf FILE "Angles\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'angles'}{'type'}[$i]) &&
          defined($sys->{'angles'}{'atoms'}[$i]))
      {
        $type = $sys->{'angles'}{'type'}[$i];
        @temp = @{$sys->{'angles'}{'atoms'}[$i]};
      }
      else
      {
        errExit("Data for angle $i is not defined.");
      }

      printf FILE " $i $type @temp\n";
    }
    printf FILE "\n";
  }

# Dihedrals
  if (defined($sys->{'diheds'}{'count'}) &&
      defined($sys->{'diheds'}{'type'}) &&
      defined($sys->{'diheds'}{'atoms'}))
  {
    $num = $sys->{'diheds'}{'count'};
    printf FILE "Dihedrals\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'diheds'}{'type'}[$i]) &&
          defined($sys->{'diheds'}{'atoms'}[$i]))
      {
        $type = $sys->{'diheds'}{'type'}[$i];
        @temp = @{$sys->{'diheds'}{'atoms'}[$i]};
      }
      else
      {
        errExit("Data for dihedral $i is not defined.");
      }

      printf FILE " $i $type @temp\n";
    }
    printf FILE "\n";
  }

# Impropers
  if (defined($sys->{'improps'}{'count'}) &&
      defined($sys->{'improps'}{'type'}) &&
      defined($sys->{'improps'}{'atoms'}))
  {
    $num = $sys->{'improps'}{'count'};
    printf FILE "Impropers\n\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'improps'}{'type'}[$i]) &&
          defined($sys->{'improps'}{'atoms'}[$i]))
      {
        $type = $sys->{'improps'}{'type'}[$i];
        @temp = @{$sys->{'improps'}{'atoms'}[$i]};
      }
      else
      {
        errExit("Data for improper $i is not defined.");
      }

      printf FILE " $i $type @temp\n";
    }
    printf FILE "\n";
  }

# Close file
  close FILE;

# Warning for default values
  if ($flag) {
    warning("Not all data needed for LAMMPS data format is defined, ".
        "using some default values.");
  }
}

# writeTypes( $file, \%sys )
# Write types file for given molecular system
sub writeTypes
{
# Variables
  my $file = $_[0];
  my $sys = $_[1];
  my ($num, @types);

# Open file
  open FILE, "> $file" or die "Error opening file '$file': $!";

# Atom
  if (defined($sys->{'atomTypes'}{'count'}) &&
      defined($sys->{'atomTypes'}{'name'}))
  {
    $num = $sys->{'atomTypes'}{'count'};
    @types = @{$sys->{'atomTypes'}{'name'}};

    printf FILE "atom types\n";
    for (my $i=1; $i <= $num; $i++) {
      printf FILE "%-5d %s\n", $i, $types[$i];
    }
  }

# Bond
  if (defined($sys->{'bondTypes'}{'count'}) &&
      defined($sys->{'bondTypes'}{'name'}))
  {
    $num = $sys->{'bondTypes'}{'count'};
    @types = @{$sys->{'bondTypes'}{'name'}};

    printf FILE "bond types\n";
    for (my $i=1; $i <= $num; $i++) {
      printf FILE "%-5d %s\n", $i, $types[$i];
    }
  }

# Angle
  if (defined($sys->{'angleTypes'}{'count'}) &&
      defined($sys->{'angleTypes'}{'name'}))
  {
    $num = $sys->{'angleTypes'}{'count'};
    @types = @{$sys->{'angleTypes'}{'name'}};

    printf FILE "angle types\n";
    for (my $i=1; $i <= $num; $i++) {
      printf FILE "%-5d %s\n", $i, $types[$i];
    }
  }

# Dihedral
  if (defined($sys->{'dihedTypes'}{'count'}) &&
      defined($sys->{'dihedTypes'}{'name'}))
  {
    $num = $sys->{'dihedTypes'}{'count'};
    @types = @{$sys->{'dihedTypes'}{'name'}};

    printf FILE "dihedral types\n";
    for (my $i=1; $i <= $num; $i++) {
      printf FILE "%-5d %s\n", $i, $types[$i];
    }
  }

# Improper
  if (defined($sys->{'impropTypes'}{'count'}) &&
      defined($sys->{'impropTypes'}{'name'}))
  {
    $num = $sys->{'impropTypes'}{'count'};
    @types = @{$sys->{'impropTypes'}{'name'}};

    printf FILE "improper types\n";
    for (my $i=1; $i <= $num; $i++) {
      printf FILE "%-5d %s\n", $i, $types[$i];
    }
  }

# Close file
  close FILE;
}

# writePdb( $file, \%sys )
# Write PDB file for given molecular system
sub writePdb
{
# Variables
  my $file = $_[0];
  my $sys = $_[1];
  my ($num, $type, $name, $res, $resNum, @pos, @bonded);
  my $flag = 0;

# Adjust distance units if needed
  if ($sys->{'flags'}{'ang'} == 0) {
    swapUnits($sys);
  }

# Open file
  open FILE, "> $file" or die "Error opening file '$file': $!";

# Header
  printf FILE "TITLE     %s\n", $sys->{'header'};
  printf FILE "REMARK     Generated using Polymatic\n";
  printf FILE "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
         $sys->{'boxDims'}{'x'}[2], $sys->{'boxDims'}{'y'}[2],
         $sys->{'boxDims'}{'z'}[2], 90, 90, 90, "P 1", 1;

# Atoms
         printf FILE "MODEL     1\n";
         if (defined($sys->{'atoms'}{'count'}) &&
             defined($sys->{'atoms'}{'pos'}))
         {
           $num = $sys->{'atoms'}{'count'};
           for (my $i=1; $i <= $num; $i++)
           {
             if (defined($sys->{'atoms'}{'pos'}[$i]))
             {
               @pos = @{$sys->{'atoms'}{'pos'}[$i]};
             }
             else
             {
               errExit("Position of atom $i is not defined.");
             }

             if (defined($sys->{'atoms'}{'name'}[$i]))
             {
               $name = $sys->{'atoms'}{'name'}[$i];
             }
             elsif (defined($sys->{'atoms'}{'type'}[$i]))
             {
               $type = $sys->{'atoms'}{'type'}[$i];
               if (defined($sys->{'atomTypes'}{'name'}[$type]))
               {
                 $name = $sys->{'atomTypes'}{'name'}[$type];
                 $name = substr($name,0,4);
               }
               else
               {
                 $name = $type;
               }
             }
             else
             {
               $name = 'C';
               $flag = 1;
             }

             if (defined($sys->{'atoms'}{'res'}[$i]))
             {
               $res = substr($sys->{'atoms'}{'res'}[$i], 0, 3);
             }
             else
             {
               $res = "MOL";
               $flag = 1;
             }

             if (defined($sys->{'atoms'}{'resNum'}[$i]))
             {
               $resNum = $sys->{'atoms'}{'resNum'}[$i];
             }
             elsif (defined($sys->{'atoms'}{'mol'}[$i]))
             {
               $resNum = $sys->{'atoms'}{'mol'}[$i];
             }
             else
             {
               $resNum = 1;
               $flag = 1;
             }

             printf FILE "ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f".
               "  1.00  0.00\n", $i, $name, $res, $resNum,
               $pos[0], $pos[1], $pos[2];
           }
         }
         else
         {
           errExit("Atoms are not defined, will not write PDB file.");
         }

# Bonds
         if (defined($sys->{'atoms'}{'bonded'}))
         {
           for (my $i=1; $i <= $num; $i++)
           {
             if (defined($sys->{'atoms'}{'bonded'}[$i]))
             {
               @bonded = @{$sys->{'atoms'}{'bonded'}[$i]};
               printf FILE "CONECT%5d", $i;
               for (my $j=0; $j < scalar(@bonded); $j++) {
                 printf FILE "%5d", $bonded[$j];
               }

               printf FILE "\n";
             }
           }
         }
         printf FILE "TER\n";
         printf FILE "ENDMDL\n";

# Warning for default values
         if ($flag) {
           warning("Not all data needed for PDB format is defined, ".
               "using some default values.");
         }

# Close file
         close FILE
}

# writePsf( $file, \%sys )
# Write PSF file for given molecular system
sub writePsf
{
# Variables
  my $file = $_[0];
  my $sys = $_[1];
  my ($num, $resNum, $res, $name, $type, $q, $mass, $n, $a1, $a2);
  my $flag = 0;

# Open file
  open FILE, "> $file" or die "Error opening file '$file': $!";

# Header
  printf FILE "PSF\n\n";
  printf FILE " %7d !TITLE\n", 3;
  printf FILE " REMARKS %s\n", $sys->{'header'};
  printf FILE " REMARKS Generated with Polymatic\n";
  printf FILE " REMARKS To be used with coordinates file in VMD\n\n";

# Atoms
  if (defined($sys->{'atoms'}{'count'}))
  {
    $num = $sys->{'atoms'}{'count'};
    printf FILE " %7d !NATOM\n", $num;
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'atoms'}{'resNum'}[$i]))
      {
        $resNum = $sys->{'atoms'}{'resNum'}[$i];
      }
      elsif (defined($sys->{'atoms'}{'mol'}[$i]))
      {
        $resNum = $sys->{'atoms'}{'mol'}[$i];
      }
      else
      {
        $resNum = 1;
        $flag = 1;
      }

      if (defined($sys->{'atoms'}{'res'}[$i]))
      {
        $res = $sys->{'atoms'}{'res'}[$i];
      }
      else
      {
        $res = "MOL";
        $flag = 1;
      }

      if (defined($sys->{'atoms'}{'type'}[$i]))
      {
        $type = $sys->{'atoms'}{'type'}[$i];
      }
      else
      {
        $type = 0;
        $flag = 1;
      }

      if (defined($sys->{'atoms'}{'name'}[$i]))
      {
        $name = $sys->{'atoms'}{'name'}[$i];
      }
      elsif (defined($sys->{'atomTypes'}{'name'}[$type]))
      {
        $name = $sys->{'atomTypes'}{'name'}[$type];
        $name = substr($name,0,4);
      }
      else
      {
        $name = 'C';
        $flag = 1;
      }

      if (defined($sys->{'atoms'}{'q'}[$i]))
      {
        $q = $sys->{'atoms'}{'q'}[$i];
      }
      else
      {
        $q = 0.0;
        $flag = 1;
      }

      if (defined($sys->{'atomTypes'}{'coeffs'}{'mass'}[$type][0]))
      {
        $mass = $sys->{'atomTypes'}{'coeffs'}{'mass'}[$type][0];
      }
      else
      {
        $mass = 1.0;
        $flag = 1;
      }

      printf FILE " %7d MOL  %-4d %-4s %-4s %-4s  ".
        "%9.6f       %7.4f           0\n",
        $i, $resNum, $res, $name, $type, $q, $mass;
    }
    printf FILE "\n";
  }
  else
  {
    errExit("Atoms are not defined, will not write PSF file.");
  }

# Bonds
  if (defined($sys->{'bonds'}{'count'}) &&
      defined($sys->{'bonds'}{'atoms'}))
  {
    $num = $sys->{'bonds'}{'count'};
    printf FILE " %7d !NBOND: bonds\n", $num;
    for (my $i=1; $i <= $num; $i+=4)
    {
      for (my $j=0; $j < 4; $j++)
      {
        $n = $i + $j;
        last if ($n > $num);

        if (defined($sys->{'bonds'}{'atoms'}[$n])) {
          ($a1, $a2) = @{$sys->{'bonds'}{'atoms'}[$n]};
        } else {
          errExit("Atoms in bond $n are not defined.");
        }

        printf FILE " %7d %7d", $a1, $a2;
      }
      printf FILE "\n";
    }
  }
  else
  {
    warning("Bonds are not defined, bond section skipped in PSF file.");
  }

# Warning for default values
  if ($flag) {
    warning("Not all data needed for PSF format is defined, ".
        "using some default values.");
  }

# Close file
  close FILE;
}

# writeGro( $file, \%sys )
# Write GRO file for given molecular system
sub writeGro
{
# Variables
  my $file = $_[0];
  my $sys = $_[1];
  my ($num, $resNum, $res, $name, $type, @pos);
  my $flag = 0;

# Adjust distance units if needed
  if ($sys->{'flags'}{'ang'} == 1) {
    swapUnits($sys);
  }

# Open file
  open FILE, "> $file" or die "Error opening file '$file': $!";

# Header
  printf FILE "%s\n", $sys->{'header'};

# Atoms
  if (defined($sys->{'atoms'}{'count'}) &&
      defined($sys->{'atoms'}{'pos'}))
  {
    $num = $sys->{'atoms'}{'count'};
    printf FILE "%5d\n", $num;
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'atoms'}{'pos'}[$i]))
      {
        @pos = @{$sys->{'atoms'}{'pos'}[$i]};
      }
      else
      {
        errExit("Position for atom $i is not defined.");
      }

      if (defined($sys->{'atoms'}{'resNum'}[$i]))
      {
        $resNum = $sys->{'atoms'}{'resNum'}[$i];
      }
      elsif(defined($sys->{'atoms'}{'mol'}[$i]))
      {
        $resNum = $sys->{'atoms'}{'mol'}[$i];
      }
      else
      {
        $resNum = 1;
        $flag = 1;
      }

      if (defined($sys->{'atoms'}{'res'}[$i]))
      {
        $res = $sys->{'atoms'}{'res'}[$i];
      }
      else
      {
        $res = "MOL";
        $flag = 1;
      }

      if (defined($sys->{'atoms'}{'name'}[$i]))
      {
        $name = $sys->{'atoms'}{'name'}[$i];
      }
      elsif (defined($sys->{'atoms'}{'type'}[$i]))
      {
        $type = $sys->{'atoms'}{'type'}[$i];
        if (defined($sys->{'atomTypes'}{'name'}[$type]))
        {
          $name = $sys->{'atomTypes'}{'name'}[$type];
          $name = substr($name,0,5);
        }
        else
        {
          $name = $type;
        }
      }
      else
      {
        $name = 'C';
        $flag = 1;
      }

      printf FILE "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
             $resNum, $res, $name, $i, $pos[0], $pos[1], $pos[2];
    }
  }
  else
  {
    errExit("Atoms are not defined, will not write GRO file.");
  }

# Box dimensions
  printf FILE " %f %f %f\n", $sys->{'boxDims'}{'x'}[2],
  $sys->{'boxDims'}{'y'}[2], $sys->{'boxDims'}{'z'}[2];

# Close file
  close FILE;

# Warning for default values
  if ($flag) {
    warning("Not all data needed for GRO format is defined, ".
        "using some default values.");
  }
}

# writeXyz( $file, \%sys )
# Write XYZ file for given molecular system
sub writeXyz
{
# Variables
  my $file = $_[0];
  my $sys = $_[1];
  my ($num, $type, @pos);
  my $flag = 0;

# Adjust distance units if needed
  if ($sys->{'flags'}{'ang'} == 0) {
    swapUnits($sys);
  }

# Open file
  open FILE, "> $file" or die "Error opening file '$file': $!";

  if (defined($sys->{'atoms'}{'count'}) &&
      defined($sys->{'atoms'}{'pos'}))
  {
# Header
    $num = $sys->{'atoms'}{'count'};
    printf FILE "%d\n", $num;
    printf FILE "%s\n", $sys->{'header'};

# Atoms
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'atoms'}{'pos'}[$i]))
      {
        @pos = @{$sys->{'atoms'}{'pos'}[$i]};
      }
      else
      {
        errExit("Position for atom $i is not defined.");
      }

      if (defined($sys->{'atoms'}{'name'}[$i]))
      {
        $type = $sys->{'atoms'}{'name'}[$i];
      }
      elsif (defined($sys->{'atoms'}{'type'}[$i]))
      {
        $type = $sys->{'atoms'}{'type'}[$i];
        if (defined($sys->{'atomTypes'}{'name'}[$type])) {
          $type = $sys->{'atomTypes'}{'name'}[$type];
        }
      }
      else
      {
        $type = 'C';
        $flag = 1;
      }

      printf FILE "%s %f %f %f\n", $type, $pos[0], $pos[1], $pos[2];
    }
  }
  else
  {
    errExit("Atoms are not defined, will not write XYZ file.");
  }

# Close file
  close FILE;

# Warning for default values
  if ($flag) {
    warning("Not all data needed for XYZ format is defined, ".
        "using some default values.");
  }
}

# writeLmpsTrj( $file, \%sys )
# Write LAMMPS trajectory file for given molecular system
sub writeLmpsTrj
{
# Variables
  my $file = $_[0];
  my $sys = $_[1];
  my ($num, $type, @pos, @x, @y, @z);
  my $flag = 0;

# Adjust distance units if needed
  if ($sys->{'flags'}{'ang'} == 0) {
    swapUnits($sys);
  }

# Open file
  open FILE, "> $file" or die "Error opening file '$file': $!";

# Header
  printf FILE "ITEM: TIMESTEP\n";
  printf FILE " %d\n", 0;
  printf FILE "ITEM: NUMBER OF ATOMS\n";
  printf FILE " %d\n", $sys->{'atoms'}{'count'};

# Box dimensions
  if (defined($sys->{'boxDims'}{'x'}) &&
      defined($sys->{'boxDims'}{'y'}) &&
      defined($sys->{'boxDims'}{'z'}))
  {
    @x = @{$sys->{'boxDims'}{'x'}};
    @y = @{$sys->{'boxDims'}{'y'}};
    @z = @{$sys->{'boxDims'}{'z'}};

    if (scalar(@x) == 0) {
      @x = (0, 100, 100);
      $flag = 1;
    }

    if (scalar(@y) == 0) {
      @y = (0, 100, 100);
      $flag = 1;
    }

    if (scalar(@z) == 0) {
      @z = (0, 100, 100);
      $flag = 1;
    }
  }
  else
  {
    @x = (0, 100, 100);
    @y = (0, 100, 100);
    @z = (0, 100, 100);
    $flag = 1;
  }

  printf FILE "ITEM: BOX BOUNDS pp pp pp\n";
  printf FILE " %f %f\n", $x[0], $x[1];
  printf FILE " %f %f\n", $y[0], $y[1];
  printf FILE " %f %f\n", $z[0], $z[1];

# Atoms
  if (defined($sys->{'atoms'}{'count'}) &&
      defined($sys->{'atoms'}{'pos'}))
  {
    $num = $sys->{'atoms'}{'count'};
    printf FILE "ITEM: ATOMS id type xs ys zs\n";
    for (my $i=1; $i <= $num; $i++)
    {
      if (defined($sys->{'atoms'}{'pos'}[$i])) {
        @pos = @{$sys->{'atoms'}{'pos'}[$i]};
      } else {
        errExit("Position of atom $i is not defined.");
      }

      if (defined($sys->{'atoms'}{'type'}[$i])) {
        $type = $sys->{'atoms'}{'type'}[$i];
      } else {
        $type = 1;
        $flag = 1;
      }

      printf FILE " %d %d %f %f %f\n", $i, $type, ($pos[0]-$x[0])/$x[2],
      ($pos[1]-$y[0])/$y[2], ($pos[2]-$z[0])/$z[2];
    }
  }
  else
  {
    errExit("Atoms are not defined, will not ".
        "write LAMMPS dump file.");
  }

# Close file
  close FILE;

# Warning for default values
  if ($flag) {
    warning("Not all data needed for LAMMPS dump format ".
        "is defined, using some default values.");
  }
}

# writeCar( $file, \%sys )
# Write CAR file for given molecular system
sub writeCar
{
# Variables
  my $file = $_[0];
  my $sys = $_[1];
  my ($num, $lenX, $lenY, $lenZ, @mols, @atoms);
  my ($id, $x, $y, $z, $type, $element, $q);
  my $flag = 0;
  my $time = localtime();

# Adjust distance units if needed
  if ($sys->{'flags'}{'ang'} == 0) {
    swapUnits($sys);
  }

# Open file
  open FILE, "> $file" or die "Error opening file '$file': $!";

# Header
  printf FILE "!BIOSYM archive 3\n";
  printf FILE "PBC=ON\n";
  printf FILE "CAR file generated by Polymatic for Materials Studio\n";
  printf FILE "!DATE %s\n", $time;

# Box dimensions
  if (defined($sys->{'boxDims'}{'x'}) &&
      defined($sys->{'boxDims'}{'y'}) &&
      defined($sys->{'boxDims'}{'z'}))
  {
    $lenX = @{$sys->{'boxDims'}{'x'}}[2];
    $lenY = @{$sys->{'boxDims'}{'y'}}[2];
    $lenZ = @{$sys->{'boxDims'}{'z'}}[2];

    if ($lenX == 0) {
      $lenX = 100;
      $flag = 1;
    }

    if ($lenY == 0) {
      $lenY = 100;
      $flag = 1;
    }

    if ($lenZ == 0) {
      $lenZ = 100;
      $flag = 1;
    }
  }
  else
  {
    $lenX = 100;
    $lenY = 100;
    $lenZ = 100;
    $flag = 1;
  }

  printf FILE "PBC %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f (P1)\n",
         $lenX, $lenY, $lenZ, 90, 90, 90;

# Atoms defined?
  if (!defined($sys->{'atoms'}{'count'}) ||
      !defined($sys->{'atoms'}{'pos'}))
  {
    errExit("Atoms are not defined, will not write CAR file.");
  }

# Molecule definitions
  if (defined($sys->{'mols'}{'count'}))
  {
    $num = $sys->{'mols'}{'count'};
    @mols = @{$sys->{'mols'}{'atoms'}};
  }
  else
  {
    $num = $sys->{'atoms'}{'count'};
    @atoms = (1..$num);
    $mols[1] = [@atoms];
    $num = 1;
  }

# Per molecule
  for (my $i=1; $i <= $num; $i++)
  {
    @atoms = @{$mols[$i]};
    for (my $j=0; $j < scalar(@atoms); $j++)
    {
# Atom
      $id = $atoms[$j];

# Coordinates
      if (defined($sys->{'atoms'}{'pos'}[$id]))
      {
        ($x, $y, $z) = @{$sys->{'atoms'}{'pos'}[$id]};
      }
      else
      {
        errExit("Position of atom $id is not defined.");
      }

# Apply PBC
      $x = $x - $lenX*POSIX::floor($x/$lenX);
      $y = $y - $lenY*POSIX::floor($y/$lenY);
      $z = $z - $lenZ*POSIX::floor($z/$lenZ);

# Atom type
      if (defined($sys->{'atoms'}{'type'}[$id]))
      {
        $type = $sys->{'atoms'}{'type'}[$id];
        if (defined($sys->{'atomTypes'}{'name'}[$type]))
        {
          $type = $sys->{'atomTypes'}{'name'}[$type];
          $element = substr($type,0,2);
        }
        else
        {
          $element = $type;
          $flag = 1;
        }
      }

# Atom charge
      if (defined($sys->{'atoms'}{'q'}[$id]))
      {
        $q = $sys->{'atoms'}{'q'}[$id];
      }
      else
      {
        $q = 0.0;
        $flag = 1;
      }

# Atom line
      printf FILE "%-5s %14.9f %14.9f %14.9f XXXX %-6d %-7s %-2s %6.3f\n",
             $id, $x, $y, $z, $i, $type, $element, $q;
    }

# End of molecule
    printf FILE "end\n";
  }

# End of system
  printf FILE "end\n";

# Close file
  close FILE;

# Warning for default values
  if ($flag) {
    warning("Not all data needed for CAR format is ".
        "defined, using some default values.");
  }
}

# writeMdf( $file, \%sys )
# Write MDF file for given molecular system
sub writeMdf
{
# Variables
  my $file = $_[0];
  my $sys = $_[1];
  my ($num, $id, $element, $type, $q, $x, $y, $z);
  my ($a2, $x2, $y2, $z2, $cx, $cy, $cz, $cell);
  my ($lenX, $lenY, $lenZ, @mols, @atoms, @bonded);
  my $flag = 0;
  my $time = localtime();

# Adjust distance units if needed
  if ($sys->{'flags'}{'ang'} == 0) {
    swapUnits($sys);
  }

# Open file
  open FILE, "> $file" or die "Error opening file '$file': $!";

# Box dimensions
  if (defined($sys->{'boxDims'}{'x'}) &&
      defined($sys->{'boxDims'}{'y'}) &&
      defined($sys->{'boxDims'}{'z'}))
  {
    $lenX = @{$sys->{'boxDims'}{'x'}}[2];
    $lenY = @{$sys->{'boxDims'}{'y'}}[2];
    $lenZ = @{$sys->{'boxDims'}{'z'}}[2];

    if ($lenX == 0) {
      $lenX = 100;
      $flag = 1;
    }

    if ($lenY == 0) {
      $lenY = 100;
      $flag = 1;
    }

    if ($lenZ == 0) {
      $lenZ = 100;
      $flag = 1;
    }
  }
  else
  {
    $lenX = 100;
    $lenY = 100;
    $lenZ = 100;
    $flag = 1;
  }

# Header
  printf FILE "!BIOSYM molecular_data 4\n\n";
  printf FILE "!DATE %s  MDF file generated by ".
    "Polymatic for Materials Studio\n", $time;
  printf FILE "#topology\n\n";
  printf FILE "\@column 1 element\n";
  printf FILE "\@column 2 atom_type\n";
  printf FILE "\@column 3 charge_group\n";
  printf FILE "\@column 4 isotope\n";
  printf FILE "\@column 5 formal_charge\n";
  printf FILE "\@column 6 charge\n";
  printf FILE "\@column 7 switching_atom\n";
  printf FILE "\@column 8 oop_flag\n";
  printf FILE "\@column 9 chirality_flag\n";
  printf FILE "\@column 10 occupancy\n";
  printf FILE "\@column 11 xray_temp_factor\n";
  printf FILE "\@column 12 connections\n\n";

# Atoms defined?
  if (!defined($sys->{'atoms'}{'count'}) ||
      !defined($sys->{'atoms'}{'pos'}) ||
      !defined($sys->{'atoms'}{'bonded'}))
  {
    errExit("Atoms or bonded atoms are not defined, "."
        will not write CAR file.");
  }

# Molecule definitions
  if (defined($sys->{'mols'}{'count'}))
  {
    $num = $sys->{'mols'}{'count'};
    @mols = @{$sys->{'mols'}{'atoms'}};
  }
  else
  {
    $num = $sys->{'atoms'}{'count'};
    @atoms = (1..$num);
    $mols[1] = [@atoms];
    $num = 1;
  }

# Per molecule
  for (my $i=1; $i <= $num; $i++)
  {
    printf FILE "\@molecule Molecule%d\n\n", $i;
    @atoms = @{$mols[$i]};
    for (my $j=0; $j < scalar(@atoms); $j++)
    {
# Atom number
      $id = $atoms[$j];

# Coordinates
      if (defined($sys->{'atoms'}{'pos'}[$id]))
      {
        ($x, $y, $z) = @{$sys->{'atoms'}{'pos'}[$id]};
      }
      else
      {
        errExit("Position of atom $id is not defined.");
      }

# Apply PBC
      $x = $x - $lenX*POSIX::floor($x/$lenX);
      $y = $y - $lenY*POSIX::floor($y/$lenY);
      $z = $z - $lenZ*POSIX::floor($z/$lenZ);

# Atom type
      if (defined($sys->{'atoms'}{'type'}[$id]))
      {
        $type = $sys->{'atoms'}{'type'}[$id];
        if (defined($sys->{'atomTypes'}{'name'}[$type]))
        {
          $type = $sys->{'atomTypes'}{'name'}[$type];
          $element = substr($type,0,2);
        }
        else
        {
          $element = $type;
          $flag = 1;
        }
      }

# Atom charge
      if (defined($sys->{'atoms'}{'q'}[$id]))
      {
        $q = $sys->{'atoms'}{'q'}[$id];
      }
      else
      {
        $q = 0.0;
        $flag = 1;
      }

# Bonded atoms
      if (defined($sys->{'atoms'}{'bonded'}[$id]))
      {
        @bonded = @{$sys->{'atoms'}{'bonded'}[$id]};
      }
      else
      {
        @bonded = ();
      }

# Atom record
      printf FILE "XXXX_%-14s %-2s %-7s ?     0  ".
        "0    %7.4f 0 0 8 1.0000  0.0000",
        $i.":".$id, $element, $type, $q;

# Connectivity
      for (my $k=0; $k < scalar(@bonded); $k++)
      {
# Atom number
        $a2 = $bonded[$k];

# Coordinates
        if (defined($sys->{'atoms'}{'pos'}[$a2])) {
          ($x2, $y2, $z2) = @{$sys->{'atoms'}{'pos'}[$a2]};
        } else {
          errExit("Position of atom $a2 is not defined.");
        }

# Apply PBC
        $x2 = $x2 - $lenX*POSIX::floor($x2/$lenX);
        $y2 = $y2 - $lenY*POSIX::floor($y2/$lenY);
        $z2 = $z2 - $lenZ*POSIX::floor($z2/$lenZ);

# Cell
        if (abs($x-$x2) > 0.5*$lenX) {
          $cx = ($x-$x2)/abs($x-$x2);
        } else {
          $cx = 0;
        }

        if (abs($y-$y2) > 0.5*$lenY) {
          $cy = ($y-$y2)/abs($y-$y2);
        } else {
          $cy = 0;
        }

        if (abs($z-$z2) > 0.5*$lenZ) {
          $cz = ($z-$z2)/abs($z-$z2);
        } else {
          $cz = 0;
        }

        if ($cx != 0 || $cy != 0 || $cz != 0) {
          $cell = "%".$cx.$cy.$cz;
        } else {
          $cell = "";
        }

# Connection
        printf FILE " ".$a2.$cell."/1.0";
      }

# End atom record
      printf FILE "\n";
    }

# End molecule
    printf FILE "\n";
  }

# Footer
  printf FILE "\n!\n";
  printf FILE "#symmetry\n";
  printf FILE "\@periodicity 3 xyz\n";
  printf FILE "\@group (P1)\n\n";
  printf FILE "#end\n";

# Close file
  close FILE;

# Warning for default values
  if ($flag) {
    warning("Not all data needed for MDF format is ".
        "defined, using some default values.");
  }
}

################################################################################
# Other

# swapUnits( \%sys )
# Swap units (nm <=> ang) for given molecular system
sub swapUnits
{
# Variables
  my $sys = $_[0];
  my ($scale, $num, @pos);

  if ($sys->{'flags'}{'ang'})
  {
    $scale = 0.1;
    $sys->{'flags'}{'ang'} = 0;
  }
  else
  {
    $scale = 10;
    $sys->{'flags'}{'ang'} = 1;
  }

  $num = $sys->{'atoms'}{'count'};
  for (my $i=1; $i <= $num; $i++)
  {
    @pos = @{$sys->{'atoms'}{'pos'}[$i]};
    $sys->{'atoms'}{'pos'}[$i] =
      [$pos[0]*$scale, $pos[1]*$scale, $pos[2]*$scale];
  }
}

# getBondType( \%sys, \@atoms )
# Get numerical bond type for given atoms in given system
sub getBondType
{
# Variables
  my $sys = $_[0];
  my @a = @{$_[1]};
  my ($t1, $t2, $str, $type);

# Atom string types
  $t1 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[0]]];
  $t2 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[1]]];
  errExit("Atom types for atoms in bond are not defined.")
    if (!defined($t1) || !defined($t2));

# Bond type, original order
  $str = join(',', ($t1, $t2));
  $type = $sys->{'bondTypes'}{'num'}{$str};
  return $type if (defined($type));

# Bond type, reverse order
  $str = join(',', ($t2, $t1));
  $type = $sys->{'bondTypes'}{'num'}{$str};
  return (-1*$type) if (defined($type));

# No type found
  return 0;
}

# getAngleType( \%sys, \@atoms )
# Get numerical angle type for given atoms in given system
sub getAngleType
{
# Variables
  my $sys = $_[0];
  my @a = @{$_[1]};
  my ($t1, $t2, $t3, $str, $type);

# Atom string types
# print Dumper($sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}]);
  $t1 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[0]]];
  $t2 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[1]]];
  $t3 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[2]]];
  errExit("Atom types for atoms in angle are not defined.")
    if (!defined($t1) || !defined($t2) || !defined($t3));

# Angle type, original order
  $str = join(',', ($t1, $t2, $t3));
  $type = $sys->{'angleTypes'}{'num'}{$str};
  return $type if (defined($type));

# Angle type, reverse order
  $str = join(',', ($t3, $t2, $t1));
  $type = $sys->{'angleTypes'}{'num'}{$str};
  return (-1*$type) if (defined($type));

# No type found
  return 0;
}

# getDihedType( \%sys, \@atoms )
# Get numerical dihedral type for given atoms in given system
sub getDihedType
{
# Variables
  my $sys = $_[0];
  my @a = @{$_[1]};
  my ($t1, $t2, $t3, $t4, $str, $type);

# Atom string types
# print "$a[0],$a[1],$a[2], $a[3]\n";
# print Dumper($sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[0]]]);
# print "@{[$sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[0]]]]}\n";
  $t1 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[0]]];
  $t2 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[1]]];
  $t3 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[2]]];
  $t4 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[3]]];
  errExit("Atom types for atoms in dihedral are not defined.")
    if (!defined($t1) || !defined($t2) || !defined($t3) || !defined($t4));

# Dihedral type, original order
  $str = join(',', ($t1, $t2, $t3, $t4));
  $type = $sys->{'dihedTypes'}{'num'}{$str};
  return $type if (defined($type));

# Dihedral type, reverse order
  $str = join(',', ($t4, $t3, $t2, $t1));
  $type = $sys->{'dihedTypes'}{'num'}{$str};
  return (-1*$type) if (defined($type));

# No type found
  return 0;
}

# getImpropType( \%sys, \@atoms )
# Get numerical improper type for given atoms in given system
sub getImpropType
{
# Variables
  my $sys = $_[0];
  my @a = @{$_[1]};
  my ($t1, $t2, $t3, $t4, $str, $type);

# Atom string types
  $t1 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[0]]];
  $t2 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[1]]];
  $t3 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[2]]];
  $t4 = $sys->{'atomTypes'}{'name'}[$sys->{'atoms'}{'type'}[$a[3]]];
  errExit("Atom types for atoms in improper are not defined.")
    if (!defined($t1) || !defined($t2) || !defined($t3) || !defined($t4));

# Improper type, original order
  $str = join(',', ($t1, $t2, $t3, $t4));
  $type = $sys->{'impropTypes'}{'num'}{$str};
  return ($type, 0) if (defined($type));

# Improper type, order 1
  $str = join(',', ($t1, $t2, $t4, $t3));
  $type = $sys->{'impropTypes'}{'num'}{$str};
  return ($type, 1) if (defined($type));

# Improper type, order 2
  $str = join(',', ($t3, $t2, $t1, $t4));
  $type = $sys->{'impropTypes'}{'num'}{$str};
  return ($type, 2) if (defined($type));

# Improper type, order 3
  $str = join(',', ($t3, $t2, $t4, $t1));
  $type = $sys->{'impropTypes'}{'num'}{$str};
  return ($type, 3) if (defined($type));

# Improper type, order 4
  $str = join(',', ($t4, $t2, $t1, $t3));
  $type = $sys->{'impropTypes'}{'num'}{$str};
  return ($type, 4) if (defined($type));

# Improper type, order 5
  $str = join(',', ($t4, $t2, $t3, $t1));
  $type = $sys->{'impropTypes'}{'num'}{$str};
  return ($type, 5) if (defined($type));

# No type found
  return (0, 0);
}

# defineAngles( \%sys )
# Define all angles in given system
sub defineAngles
{
# Variables
  my $sys = $_[0];
  my ($numAtoms, $numAngles, $numAngleTypes);
  my ($a1, $a2, $t1, $t2, $t3, $str, $type, @bonded);
  my $at = 'atomTypes';

# Check atoms and bonds are defined
  errExit("Atoms and bonds must be defined to define angles.")
    if (!defined($sys->{'atoms'}{'count'}) || !defined($sys->{'atoms'}) ||
        !defined($sys->{'bonds'}{'count'}) || !defined($sys->{'bonds'}));

# Check for angles already defined
  if (defined($sys->{'angles'}))
  {
    warning("Angles are already defined in the system ".
        "and will be overwritten.");
    $sys->{'angles'} = ();
  }

# Angle: bonded1-atom-bonded2
  $numAtoms = $sys->{'atoms'}{'count'};
  for (my $i=1; $i <= $numAtoms; $i++)
  {
    if (defined($sys->{'atoms'}{'bonded'}[$i])) {
      @bonded = @{$sys->{'atoms'}{'bonded'}[$i]};
    } else {
      @bonded = ();
    }

    next if (scalar(@bonded) < 2);

    for (my $j=0; $j < scalar(@bonded)-1; $j++)
    {
      $a1 = $bonded[$j];
      for (my $k=$j+1; $k < scalar(@bonded); $k++)
      {
        $a2 = $bonded[$k];
        $type = getAngleType($sys, [$a1, $i, $a2]);
        $numAngles++;

        if ($type == 0)
        {
          $numAngleTypes++;
          $t1 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$i]];
          $t2 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$a1]];
          $t3 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$a2]];
          $str = $t2.','.$t1.','.$t3;
          $sys->{'angleTypes'}{'name'}[$numAngleTypes] = $str;
          $sys->{'angleTypes'}{'num'}{$str} = $numAngleTypes;
          $sys->{'angles'}{'type'}[$numAngles] = $numAngleTypes;
          $sys->{'angles'}{'atoms'}[$numAngles] = [$a1, $i, $a2];
        }
        elsif ($type < 0)
        {
          $sys->{'angles'}{'type'}[$numAngles] = -1*$type;
          $sys->{'angles'}{'atoms'}[$numAngles] = [$a2, $i, $a1];
        }
        else
        {
          $sys->{'angles'}{'type'}[$numAngles] = $type;
          $sys->{'angles'}{'atoms'}[$numAngles] = [$a1, $i, $a2];
        }
      }
    }
  }

  $sys->{'angles'}{'count'} = $numAngles;
  $sys->{'angleTypes'}{'count'} = $numAngleTypes;
}

# defineDiheds( \%sys )
# Define all dihedrals in given system
sub defineDiheds
{
# Variables
  my $sys = $_[0];
  my ($numBonds, $numDiheds, $numDihedTypes);
  my ($a1, $a2, $a3, $a4, $t1, $t2, $t3, $t4, $str, $type);
  my (@bonded1, @bonded2);
  my $at = 'atomTypes';

# Check atoms and bonds are defined
  errExit("Atoms and bonds must be defined to define dihedrals.")
    if (!defined($sys->{'atoms'}{'count'}) || !defined($sys->{'atoms'}) ||
        !defined($sys->{'bonds'}{'count'}) || !defined($sys->{'bonds'}));

# Check for dihedrals already defined
  if (defined($sys->{'diheds'}))
  {
    warning("Dihedrals are already defined in the system ".
        "and will be overwritten.");
    $sys->{'diheds'} = ();
  }

# Dihedral: bond1-atom1-atom2-bond2
  $numBonds = $sys->{'bonds'}{'count'};
  for (my $i=1; $i <= $numBonds; $i++)
  {
    ($a1, $a2) = @{$sys->{'bonds'}{'atoms'}[$i]};
    if (defined($sys->{'atoms'}{'bonded'}[$a1])) {
      @bonded1 = @{$sys->{'atoms'}{'bonded'}[$a1]};
    } else {
      @bonded1 = ();
    }

    if (defined($sys->{'atoms'}{'bonded'}[$a2])) {
      @bonded2 = @{$sys->{'atoms'}{'bonded'}[$a2]};
    } else {
      @bonded2 = ();
    }

    next if (scalar(@bonded1) < 2 || scalar(@bonded2) < 2);

    for (my $j=0; $j < scalar(@bonded1); $j++)
    {
      $a3 = $bonded1[$j];
      next if ($a3 == $a2);
      for (my $k = 0; $k < scalar(@bonded2); $k++)
      {
        $a4 = $bonded2[$k];
        next if ($a4 == $a1 || $a3 == $a4);
        $type = getDihedType($sys, [$a3, $a1, $a2, $a4]);
        $numDiheds++;

        if ($type == 0)
        {
          $numDihedTypes++;
          $t1 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$a1]];
          $t2 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$a2]];
          $t3 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$a3]];
          $t4 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$a4]];
          $str = $t3.','.$t1.','.$t2.','.$t4;
          $sys->{'dihedTypes'}{'name'}[$numDihedTypes] = $str;
          $sys->{'dihedTypes'}{'num'}{$str} = $numDihedTypes;
          $sys->{'diheds'}{'type'}[$numDiheds] = $numDihedTypes;
          $sys->{'diheds'}{'atoms'}[$numDiheds] =
            [$a3, $a1, $a2, $a4];
        }
        elsif ($type < 0)
        {
          $sys->{'diheds'}{'type'}[$numDiheds] = -1*$type;
          $sys->{'diheds'}{'atoms'}[$numDiheds] =
            [$a4, $a2, $a1, $a3];
        }
        else
        {
          $sys->{'diheds'}{'type'}[$numDiheds] = $type;
          $sys->{'diheds'}{'atoms'}[$numDiheds] =
            [$a3, $a1, $a2, $a4];
        }
      }
    }
  }

  $sys->{'diheds'}{'count'} = $numDiheds;
  $sys->{'dihedTypes'}{'count'} = $numDihedTypes;
}

# defineImprops( \%sys )
# Define all impropers in given system
sub defineImprops
{
# Variables
  my $sys = $_[0];
  my ($numAtoms, $numImprops, $numImpropTypes);
  my ($a2, $a3, $a4, $t1, $t2, $t3, $t4, $str, $type, $order);
  my (@bonded);
  my $at = 'atomTypes';

# Check atoms and bonds are defined
  errExit("Atoms and bonds must be defined to define impropers.")
    if (!defined($sys->{'atoms'}{'count'}) || !defined($sys->{'atoms'}) ||
        !defined($sys->{'bonds'}{'count'}) || !defined($sys->{'bonds'}));

# Check for impropers already defined
  if (defined($sys->{'improps'})) {
    warning("Impropers are already defined in the system ".
        "and will be overwritten.");
    $sys->{'improps'} = ();
  }

# Improper: bond1-atom-bond2-bond3
  $numAtoms = $sys->{'atoms'}{'count'};
  for (my $i=1; $i <= $numAtoms; $i++)
  {
    if (defined($sys->{'atoms'}{'bonded'}[$i])) {
      @bonded = @{$sys->{'atoms'}{'bonded'}[$i]};
    } else {
      @bonded = ();
    }

    next if (scalar(@bonded) < 3);

    for (my $j=0; $j < scalar(@bonded)-2; $j++)
    {
      $a2 = $bonded[$j];
      for (my $k=$j+1; $k < scalar(@bonded)-1; $k++)
      {
        $a3 = $bonded[$k];
        for (my $l=$k+1; $l < scalar(@bonded); $l++)
        {
          $a4 = $bonded[$l];
          ($type, $order) = getImpropType($sys, [$a2, $i, $a3, $a4]);
          $numImprops++;

          if ($type == 0)
          {
            $numImpropTypes++;
            $t1 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$i]];
            $t2 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$a2]];
            $t3 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$a3]];
            $t4 = $sys->{$at}{'name'}[$sys->{'atoms'}{'type'}[$a4]];
            $str = $t2.','.$t1.','.$t3.','.$t4;
            $sys->{'impropTypes'}{'name'}[$numImpropTypes] = $str;
            $sys->{'impropTypes'}{'num'}{$str} = $numImpropTypes;
            $sys->{'improps'}{'type'}[$numImprops] =
              $numImpropTypes;
            $sys->{'improps'}{'atoms'}[$numImprops] =
              [$a2, $i, $a3, $a4];
          }
          else
          {
            $sys->{'improps'}{'type'}[$numImprops] = $type;
            if ($order == 0) {
              $sys->{'improps'}{'atoms'}[$numImprops] =
                [$a2, $i, $a3, $a4];
            } elsif ($order == 1) {
              $sys->{'improps'}{'atoms'}[$numImprops] =
                [$a2, $i, $a4, $a3];
            } elsif ($order == 2) {
              $sys->{'improps'}{'atoms'}[$numImprops] =
                [$a3, $i, $a2, $a4];
            } elsif ($order == 3) {
              $sys->{'improps'}{'atoms'}[$numImprops] =
                [$a3, $i, $a4, $a2];
            } elsif ($order == 4) {
              $sys->{'improps'}{'atoms'}[$numImprops] =
                [$a4, $i, $a2, $a3];
            } else {
              $sys->{'improps'}{'atoms'}[$numImprops] =
                [$a4, $i, $a3, $a2];
            }
          }
        }
      }
    }
  }

  $sys->{'improps'}{'count'} = $numImprops;
  $sys->{'impropTypes'}{'count'} = $numImpropTypes;
}

# defineMols( \%sys )
# Define molecules for given system
sub defineMols
{
# Variables
  my $sys = $_[0];
  my ($numAtoms, $numMols, $a);
  my (@mol, @bonded, @queue);

# Check atoms and bonds are defined
  errExit("Atoms and bonds must be defined to define molecules.")
    if (!defined($sys->{'atoms'}{'count'}) || !defined($sys->{'atoms'}) ||
        !defined($sys->{'bonds'}{'count'}) || !defined($sys->{'bonds'}));

# Check for molecules already defined
  if (defined($sys->{'mols'})) {
    warning("Molecules are already defined in the system ".
        "and will be overwritten.");
    $sys->{'mols'} = ();
  }

# Recursively add bonded molecules
  $numAtoms = $sys->{'atoms'}{'count'};
  for (my $i=1; $i <= $numAtoms; $i++)
  {
    next if (defined($sys->{'atoms'}{'mol'}[$i]));
    $numMols++;
    $sys->{'atoms'}{'mol'}[$i] = $numMols;
    push(@mol, $i);
    push(@queue, $i);

    while (scalar(@queue) > 0)
    {
      $a = shift(@queue);
      if (defined($sys->{'atoms'}{'bonded'}[$a])) {
        @bonded = @{$sys->{'atoms'}{'bonded'}[$a]};
      } else {
        @bonded = ();
      }

      for (my $j=0; $j < scalar(@bonded); $j++)
      {
        if (!defined($sys->{'atoms'}{'mol'}[$bonded[$j]]))
        {
          $sys->{'atoms'}{'mol'}[$bonded[$j]] = $numMols;;
          push(@mol, $bonded[$j]);
          push(@queue, $bonded[$j]);
        }
      }
    }

    $sys->{'mols'}{'atoms'}[$numMols] = [@mol];
    @mol = ();
  }

  $sys->{'mols'}{'count'} = $numMols;
}

# getCell( \%sys, $atom )
# Calculate neighbor list cell for given atom in given system
sub getCell
{
# Variables
  my $sys = $_[0];
  my $atom = $_[1];
  my ($cellN, $cellL, $x, $y, $z, $xc, $yc, $zc);

  ($x, $y, $z) = @{$sys->{'atoms'}{'pos'}[$atom]};
  $cellN = $sys->{'neighList'}{'num'};
  $cellL = $sys->{'neighList'}{'width'};

  $xc = POSIX::floor($x/$cellL);
  $xc = $xc - $cellN * POSIX::floor($xc/$cellN);
  $yc = POSIX::floor($y/$cellL);
  $yc = $yc - $cellN * POSIX::floor($yc/$cellN);
  $zc = POSIX::floor($z/$cellL);
  $zc = $zc - $cellN * POSIX::floor($zc/$cellN);

  return ($xc, $yc, $zc);
}

# initNeighList( \%sys, $cellN, $cellL )
# Initialize neighbor list for given system with given cell number and width
sub initNeighList
{
# Variables
  my $sys = $_[0];
  my $cellN = $_[1];
  my $cellL = $_[2];
  my (@neigh, $u, $v, $w, $num);

  $sys->{'neighList'}{'num'} = $cellN;
  $sys->{'neighList'}{'width'} = $cellL;

  for (my $i=0; $i < $cellN; $i++) {
    for (my $j=0; $j < $cellN; $j++) {
      for (my $k=0; $k < $cellN; $k++) {
        @{$neigh[$i][$j][$k]} = ();
      }}}

  $num = $sys->{'atoms'}{'count'};
  for (my $i=1; $i <= $num; $i++)
  {
    ($u, $v, $w) = getCell($sys, $i);
    push(@{$neigh[$u][$v][$w]}, $i);
  }
  $sys->{'neighList'}{'cells'} = [@neigh];
}

# addNeighList( \%sys, $atom )
# Add given atom in given system to neighbor list
sub addNeighList
{
# Variables
  my $sys = $_[0];
  my $atom = $_[1];
  my ($u, $v, $w);

  ($u, $v, $w) = getCell($sys, $atom);
  push(@{$sys->{'neighList'}{'cells'}[$u][$v][$w]}, $atom)
}

# getSep( \%sys, $atom1, $atom2 )
# Calcualte separation between given atoms with nearest image convention
sub getSep
{
# Variables
  my $sys = $_[0];
  my $atom1 = $_[1];
  my $atom2 = $_[2];
  my (@pos1, @pos2, $lenX, $lenY, $lenZ, $x, $y, $z);

  @pos1 = @{$sys->{'atoms'}{'pos'}[$atom1]};
  @pos2 = @{$sys->{'atoms'}{'pos'}[$atom2]};
  $lenX = $sys->{'boxDims'}{'x'}[2];
  $lenY = $sys->{'boxDims'}{'y'}[2];
  $lenZ = $sys->{'boxDims'}{'z'}[2];

  $x = $pos1[0] - $pos2[0];
  $y = $pos1[1] - $pos2[1];
  $z = $pos1[2] - $pos2[2];

  $x = $x - $lenX * POSIX::floor($x/$lenX + 0.5);
  $y = $y - $lenY * POSIX::floor($y/$lenY + 0.5);
  $z = $z - $lenZ * POSIX::floor($z/$lenZ + 0.5);

  return sqrt($x*$x + $y*$y + $z*$z);
}

# unwrapMol( \%sys, $ref )
# Unwrap molecule coordinates starting from given reference atom
sub unwrapMol
{
# Variables
  my $sys = $_[0];
  my $ref = $_[1];
  my ($a1, $a2, $x1, $y1, $z1, $x2, $y2, $z2);
  my ($lx, $ly, $lz, $flag, $x, $y, $z);
  my (@todo, @bonded);

# Box dimensions
  if (defined(@{$sys->{'boxDims'}{'x'}}[2]) &&
      defined(@{$sys->{'boxDims'}{'y'}}[2]) &&
      defined(@{$sys->{'boxDims'}{'z'}}[2]))
  {
    $lx = @{$sys->{'boxDims'}{'x'}}[2];
    $ly = @{$sys->{'boxDims'}{'y'}}[2];
    $lz = @{$sys->{'boxDims'}{'z'}}[2];
  }
  else
  {
    warning("Box dimensions are undefined, can't unwrap molecules.");
    return 0;
  }

# Initial atom
  push(@todo, $ref);
  $sys->{'atoms'}{'checked'}[$ref] = 1;

# Loop through atoms in molecule
  while (scalar(@todo) > 0)
  {
    $a1 = shift(@todo);
    ($x1, $y1, $z1) = @{$sys->{'atoms'}{'pos'}[$a1]};

    if (defined($sys->{'atoms'}{'bonded'}[$a1])) {
      @bonded = @{$sys->{'atoms'}{'bonded'}[$a1]};
    } else {
      @bonded = ()
    }

    for (my $i=0; $i < scalar(@bonded); $i++)
    {
      $a2 = $bonded[$i];
      next if ($sys->{'atoms'}{'checked'}[$a2] == 1);
      push(@todo, $a2);
      $sys->{'atoms'}{'checked'}[$a2] = 1;
      ($x2, $y2, $z2) = @{$sys->{'atoms'}{'pos'}[$a2]};

# Adjust coordinates
      $flag = 0;
      while ($flag == 0)
      {
        $x = $x1 - $x2;

        if (abs($x) > $lx/2)
        {
          if ($x > 0) {
            $x2 += $lx;
          } else {
            $x2 -= $lx;
          }
        }
        else
        {
          $flag = 1;
        }
      }

      $flag = 0;
      while ($flag == 0)
      {
        $y = $y1 - $y2;

        if (abs($y) > $ly/2)
        {
          if ($y > 0) {
            $y2 += $ly;
          } else {
            $y2 -= $ly;
          }
        }
        else
        {
          $flag = 1;
        }
      }

      $flag = 0;
      while ($flag == 0)
      {
        $z = $z1 - $z2;

        if (abs($z) > $lz/2)
        {
          if ($z > 0) {
            $z2 += $lz;
          } else {
            $z2 -= $lz;
          }
        }
        else
        {
          $flag = 1;
        }
      }

      $sys->{'atoms'}{'pos'}[$a2] = [$x2, $y2, $z2];
    }
  }
}

# delDupBonds( \@array )
# Delete duplicate bonds from given array, considering both orders
sub delDupBonds
{
# Variables
  my @array = @{$_[0]};
  my (@unique, @a1, @a2, $flag);

  return @array if (scalar(@array) < 2);
  push(@unique, $array[0]);

  for (my $i=1; $i < scalar(@array); $i++)
  {
    $flag = 0;
    @a1 = @{$array[$i]};
    for (my $j=0; $j < scalar(@unique); $j++)
    {
      @a2 = @{$unique[$j]};
      if (eqArray(\@a1, \@a2) || eqArray(\@a1, [reverse(@a2)]))
      {
        $flag = 1;
        last;
      }
    }
    push(@unique, [@a1]) if ($flag == 0);
  }

  return @unique;
}

# delDupImprops( \@array )
# Delete duplicate impropers from given array, considering all orders
sub delDupImprops
{
# Variables
  my @array = @{$_[0]};
  my (@unique, @a1, @a2, $flag);

  return @array if (scalar(@array) < 2);
  push(@unique, $array[0]);

  for (my $i=1; $i < scalar(@array); $i++)
  {
    $flag = 0;
    @a1 = @{$array[$i]};
    for (my $j=0; $j < scalar(@unique); $j++)
    {
      @a2 = @{$unique[$j]};
      if (eqArray(\@a1, \@a2) ||
          eqArray(\@a1, [$a2[0], $a2[1], $a2[3], $a2[2]]) ||
          eqArray(\@a1, [$a2[2], $a2[1], $a2[0], $a2[3]]) ||
          eqArray(\@a1, [$a2[2], $a2[1], $a2[3], $a2[0]]) ||
          eqArray(\@a1, [$a2[3], $a2[1], $a2[0], $a2[2]]) ||
          eqArray(\@a1, [$a2[3], $a2[1], $a2[2], $a2[0]]))
      {
        $flag = 1;
        last;
      }
    }
    push(@unique, [@a1]) if ($flag == 0);
  }

  return @unique;
}

# vectorSub( \@vec1, \@vec2 )
# Vector subtraction between given vectors, @v1 - @v2
sub vectorSub
{
  my @v1 = @{$_[0]};
  my @v2 = @{$_[1]};
  my @v;

  errExit("Vectors must be the same length for subtraction.")
    if (scalar(@v1) != scalar(@v2));

  for (my $i=0; $i < scalar(@v1); $i++) {
    $v[$i] = $v1[$i] - $v2[$i];
  }

  return @v;
}

# dot( \@vec1, \@vec2 )
# Dot product for given vectors
sub dot
{
  my @v1 = @{$_[0]};
  my @v2 = @{$_[1]};
  my $sum=0;

  errExit("Vectors must be the same length for dot product.")
    if (scalar(@v1) != scalar(@v2));

  for (my $i=0; $i < scalar(@v1); $i++) {
    $sum += ($v1[$i] * $v2[$i]);
  }

  return $sum;
}

# norm( \@vec )
# Normal of given vector
sub norm
{
  my @v = @{$_[0]};
  my $sum=0;

  for (my $i=0; $i < scalar(@v); $i++) {
    $sum += ($v[$i] * $v[$i]);
  }

  return sqrt($sum);
}

# vectorAng( \@vec1, \@vec2 )
# Angle between given vectors
sub vectorAng
{
  my @v1 = @{$_[0]};
  my @v2 = @{$_[1]};
  return (Math::Trig::acos(dot(\@v1, \@v2) / (norm(\@v1) * norm(\@v2))));
}

# normalPlane( \@coords )
# Normal vector to best fit plane for given set of coordinates
sub normalPlane
{
  my @coords = @{$_[0]};
  my ($xi, $yi, $zi, $x, $y, $z, $a, $b, $temp);
  my ($xx, $yy, $zz, $xy, $yz, $xz);
  my $n = scalar(@coords);

  for (my $i=0; $i < $n; $i++)
  {
    ($xi, $yi, $zi) = @{$coords[$i]};

    $x += $xi;
    $y += $yi;
    $z += $zi;
    $xx += ($xi*$xi);
    $yy += ($yi*$yi);
    $zz += ($zi*$zi);
    $xy += ($xi*$yi);
    $xz += ($xi*$zi);
    $yz += ($yi*$zi);
  }

# www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
  $temp = $n*$xy*$xy - 2*$x*$xy*$y + $xx*$y*$y + $x*$x*$yy - $n*$xx*$yy;
  $a = -(-$xz*$y*$y + $n*$xz*$yy - $n*$xy*$yz + $x*$y*$yz + $xy*$y*$z
      - $x*$yy*$z) / $temp;
  $b = -(-$n*$xy*$xz + $x*$xz*$y - $x*$x*$yz + $n*$xx*$yz + $x*$xy*$z
      - $xx*$y*$z) / $temp;

  return ($a, $b, -1);
}

# eqArray( \@array1, \@array2 )
# Determine if given arrays are equal element by element
sub eqArray
{
  my @a1 = @{$_[0]};
  my @a2 = @{$_[1]};

  return 0 if (scalar(@a1) != scalar(@a2));
  for (my $i=0; $i < scalar(@a1); $i++) {
    return 0 if ($a1[$i] != $a2[$i] && $a1[$i] ne $a2[$i]);
  }

  return 1;
}

# uniqueArray( @array )
# Keep only unique values in array (remove duplicates)
sub uniqueArray
{
  my %hash = map { $_ => 1 } @_;
  my @unique = keys %hash;
  return @unique;
}

# group( \@array, $val)
# Indices of given array that have given value, all i such that array[i] = val
sub group
{
  my @array = @{$_[0]};
  my $val = $_[1];
  my @match;

  for (my $i=0; $i < scalar(@array); $i++) {
    push (@match, $i) if ($array[$i] == $val);
  }

  return @match;
}

#############################################################
#  Consider propagation rates in the form of free energies  #
#############################################################

# readBarriers( $file )
sub readBarriers
{
# Variables
  my $file = $_[0];
  my %barriers;
  my (@temp, $num);

# Open file
  open FILE, "< $file" or die "Error opening file '$file': $!";

# Read in each line
  while(my $line = <FILE>)
  {
    @temp = chompSplit($line);
    next if (substr($line,0,3) eq "mol");
    $barriers{$temp[0]} = $temp[1];
  }
  return %barriers;
}

# atomcounts(@atoms)
sub atomcounts
{
  my @molatoms = @{$_[0]};
  my %molcounts;

  foreach my $atom (@molatoms){
    $molcounts{$atom}++
  }
# remove empty keys
  my %ret;
  foreach my $key (keys %molcounts){
    $ret{$key} = $molcounts{$key} if (!$key eq "");
  }
  return %ret;
};

# read in xyz filename from the barriers hash,
# then return atom counts
# molCountFromBarriersFile($barriers->{molref})
sub molCountFromBarriersFile{
  my $mol = $_[0];
  my %molsys = readXyz($mol);
  my @molatoms = @{$molsys{'atoms'}{'type'}};
# first element is undefined
  @molatoms = grep defined, @molatoms;
  my %molcounts = atomcounts(\@molatoms);
  return %molcounts
}
# Search the system as read in initially, before any polymerisation
# has occurred, for atom types
# search_initial_sys_for_atom_types (\%barriers, \%sys)
sub search_initial_sys_for_molecules
{
  my $barriers = $_[0];
  my $sys = $_[1];
  my %mol1counts = molCountFromBarriersFile($barriers->{'mol1'});
  my %mol2counts = molCountFromBarriersFile($barriers->{'mol2'});
  print "mol1 = $barriers->{'mol1'}\n";
  print "mol2 = $barriers->{'mol2'}\n";
  print "mol1\n"; 
  foreach my $key (keys %mol1counts) {
    print "$key, $mol1counts{$key}\n";
  }
  print "mol2\n"; 
  foreach my $key (keys %mol2counts) {
    print "$key, $mol2counts{$key}\n";
  }
   
  ### need to do this kind of thing to check all atoms of a certain molecule ID
  my ($num,$id,$type,$element);
  my (@atoms,@mols);

  # Molecule definitions
  if (defined($sys->{'mols'}{'count'}))
  {
    $num = $sys->{'mols'}{'count'};
    @mols = @{$sys->{'mols'}{'atoms'}};
  }
  else
  {
    $num = $sys->{'atoms'}{'count'};
    @atoms = (1..$num);
    $mols[1] = [@atoms];
    $num = 1;
  }
  $num = $sys->{'mols'}{'count'};
  # Per molecule
  my %typesOfSystem;
  for (my $i=1; $i <= $num; $i++)
  {
    my @atomTypesPerMol;
    @atoms = @{$mols[$i]};
    for (my $j=0; $j < scalar(@atoms); $j++)
    {
      # Atom type
      $id = $atoms[$j];
      $type = $sys->{'atoms'}{'type'}[$id];
      $type = $sys->{'atomTypes'}{'name'}[$type];
      push(@atomTypesPerMol, $type);
    }
    $typesOfSystem{$i} = \@atomTypesPerMol;
  }
  # for each molecule, store counts of atom types for reference
  # later in the polym loop
  my %molsOfSystem;
  foreach my $molnum (keys %typesOfSystem){
    my @mol = @{$typesOfSystem{$molnum}};
    my %count = atomcounts(\@mol);
    # Now check that the count equals a count of a molecule
    my $mol_is_mol1 = 1;
    foreach my $type (keys %count){
      # compare count to mol1counts and mol2counts, return 1 or 2
      # if any number of types is different, then the molecule in
      # question is not mol1, and must be mol2
      # in future, may expand to more molecules- this check would be different in
      # that case
      if (!exists $mol1counts{$type}){
        $mol_is_mol1 = 0;
        print "$type not in mol1\n";
        last;
      }
      # atom type in mol2 might not exist in mol1
      if (!$count{$type} eq $mol1counts{$type}){
        $mol_is_mol1 = 0;
        print "$type not in mol1 in the same quantity\n";
        # check molecule is 2
        if(!$count{$type} eq $mol2counts{$type}){
          errExit("Molecule $type not found in barriers.in")
        }
        last;
      }
    }
    if ($mol_is_mol1){
      $molsOfSystem{$molnum} = 1;
    } else {
      $molsOfSystem{$molnum} = 2;
    }
  }
  return %molsOfSystem
}

# written during initial check (polym_init.pl)
sub writeNumberedMols
{
  my $file = $_[0];
  my %molIDs = %{$_[1]};

  open FILE, "> $file" or die "Error opening file '$file': $!";
  print FILE "Molecule   Type\n";
  foreach my $id (keys %molIDs){
    printf FILE "%-10d %d\n", $id, $molIDs{$id};
  }
  close FILE;
}


# readNumberedMols( $file ) # return hash
sub readNumberedMols
{
  # Variables
  my $file = $_[0];
  my %mols;
  my @temp;

  # Open file
  open FILE, "< $file" or die "Error opening file '$file': $!";

  # Read in each line
  while(my $line = <FILE>)
  {
    @temp = chompSplit($line);
    next if ($temp[0] eq "Molecule");
    $mols{$temp[0]} = $temp[1];
  }
  close FILE;
  return %mols;

}

1;

# countMols( $file )
sub countMols
{
  my $file = $_[0];
  my %counts;
  my @temp;
  open FILE, "< $file" or die "Error opening file '$file': $!";

  # Read in each line
  while(my $line = <FILE>)
  {
    @temp = chompSplit($line);
    next if ($temp[0] eq "Molecule");
    if (exists($counts{$temp[1]})){
      $counts{$temp[1]}++;
    } else {
      $counts{$temp[1]} = 1;
    }
  }
  close FILE;
  return %counts;
}

#######################################################
#  tracking monomers included in the growing polymer  #
#######################################################

# writeMolInfo ( \%counts , $file)
sub writeMolInfo
{
  my $counts = $_[0];
  my $file = $_[1];
  open FILE, "> $file" or die "Error opening file '$file': $!";
  if (defined($counts->{'1'})) {
    print FILE "total mol1 $counts->{'1'}\n"
  } else {
    print FILE "total mol1 0\n"
  }
  if (defined($counts->{'2'})) {
    print FILE "total mol2 $counts->{'2'}\n"
  } else {
    print FILE "total mol2 0\n"
  }
  print FILE "mol1 included 0\n"; # increment these in polymerisation
  print FILE "mol2 included 0\n";
  close FILE;
}

###############################################
#  the next 3 functions are used in polym.pl  #
###############################################

# incrementMolIncluded ( $mol, $file ) # 1 or 2
sub incrementMolIncluded
{
  my $mol = $_[0];
  my $file = $_[1];
  my ($string, $line, $num, $newnum);
  my @temp;
  $string = "mol$mol included";
  $line = `grep "$string" $file`;
  @temp = chompSplit($line);
  $num = $temp[2];
  $newnum = $num + 1;

  rename($file, $file . '.bak');
  open (my $input, '<' . $file . '.bak') or die ('File err: '. $!);
  open (my $output, '>' . $file) or die $!;
  while (<$input>)
  {
      $_ =~ s/mol$mol included $num/mol$mol included $newnum/;
      print $output $_;
  }
  close($input);
  close($output);
  unlink($file . '.bak') or die "Can't remove $file.bak: $!";
}

# molsCurrentlyIncluded ($mol, $file) # 1 or 2, bonded.tmp
sub molsCurrentlyIncluded
{
  my $mol = $_[0];
  my $file = $_[1];
  my ($sed, $string, $line, $num);
  my @temp;
  $line = `grep "mol$mol included" $file`;
  @temp = chompSplit($line);
  $num = $temp[2];
  return $num;
}

# totalIncluded ($mol, $file) # 1 or 2
sub totalMolsIncluded
{
  my $mol = $_[0];
  my $file = $_[1];
  my ($sed, $string, $line, $num);
  my @temp;
  $line = `grep "total mol$mol" $file`;
  @temp = chompSplit($line);
  $num = $temp[2];
  return $num;
}
