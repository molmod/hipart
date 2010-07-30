# HiPart is a software toolkit to analyse molecular densities with the hirshfeld partitioning scheme.
# Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of HiPart.
#
# HiPart is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HiPart is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --






from molmod import angstrom, Rotation
from molmod.periodic import periodic

import numpy, sys, os


__all__ = [
    "get_atom_grid",
    "load_atom_scalars", "load_atom_vectors", "dump_atom_scalars", "dump_atom_vectors"
]


def get_atom_grid(lebedev_xyz, center, radii):
    num_lebedev = len(lebedev_xyz)
    grid_points = numpy.zeros((num_lebedev*len(radii),3), float)
    counter = 0
    for r in radii:
        rot = Rotation.random()
        grid_points[counter:counter+num_lebedev] = r*numpy.dot(lebedev_xyz,rot.r)+center
        counter += num_lebedev
    return grid_points


def load_atom_scalars(filename):
    f = file(filename)
    # read the number of atoms
    line = f.next()
    N = int(line[line.rfind(" "):])
    # read the charges
    charges = numpy.zeros(N, float)
    f.next()
    f.next()
    for i in xrange(N):
        line = f.next()
        words = line.split()
        charges[i] = float(words[3])
    f.close()
    return charges


def load_atom_vectors(filename):
    f = file(filename)
    # read the number of atoms
    line = f.next()
    N = int(line[line.rfind(" "):])
    # read the dipoles
    dipoles = numpy.zeros((N,3), float)
    f.next()
    f.next()
    for i in xrange(N):
        line = f.next()
        words = line.split()
        dipoles[i, 0] = float(words[3])
        dipoles[i, 1] = float(words[4])
        dipoles[i, 2] = float(words[5])
    f.close()
    return dipoles


def dump_atom_scalars(filename, charges, numbers=None, name="Charge"):
    f = file(filename, "w")
    print >> f, "number of atoms: %i" % len(charges)
    print >> f, "  i        Z  %s" % name.rjust(10)
    print >> f, "-----------------------------"
    for i in xrange(len(charges)):
        if numbers is None:
            number = 0
            symbol = "?"
        else:
            number = numbers[i]
            symbol = periodic[number].symbol
        print >> f, "% 3i  %2s  % 3i   % 10.5f" % (
            i+1, symbol, number, charges[i]
        )
    print >> f, "-----------------------------"


def dump_atom_vectors(filename, dipoles, numbers=None, name="Dipole"):
    name = name.rjust(10)
    f = file(filename, "w")
    print >> f, "number of atoms:", len(dipoles)
    print >> f, "  i        Z %(name)s-X %(name)s-Y %(name)s-Z  %(name)s" % {"name": name}
    print >> f, "------------------------------------------------------------------"
    for i in xrange(len(dipoles)):
        if numbers is None:
            number = 0
            symbol = "?"
        else:
            number = numbers[i]
            symbol = periodic[number].symbol
        print >> f, "% 3i  %2s  % 3i   % 10.5f   % 10.5f   % 10.5f   % 10.5f" % (
            i+1, symbol, number, dipoles[i,0], dipoles[i,1],
            dipoles[i,2], numpy.linalg.norm(dipoles[i]),
        )
    f.close()
