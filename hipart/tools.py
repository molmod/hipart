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
    "load_atom_scalars", "load_atom_vectors",
    "dump_atom_scalars", "dump_atom_vectors", "dump_atom_matrix"
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
    scalars = numpy.zeros(N, float)
    f.next()
    f.next()
    for i in xrange(N):
        line = f.next()
        words = line.split()
        scalars[i] = float(words[3])
    f.close()
    return scalars


def load_atom_vectors(filename):
    f = file(filename)
    # read the number of atoms
    line = f.next()
    N = int(line[line.rfind(" "):])
    # read the dipoles
    vectors = numpy.zeros((N,3), float)
    f.next()
    f.next()
    for i in xrange(N):
        line = f.next()
        words = line.split()
        vectors[i, 0] = float(words[3])
        vectors[i, 1] = float(words[4])
        vectors[i, 2] = float(words[5])
    f.close()
    return vectors


def load_atom_matrix(filename):
    f = file(filename)
    # read the number of atoms
    line = f.next()
    N = int(line[line.rfind(" "):])
    # read the dipoles
    matrix = numpy.zeros((N,N), float)
    f.next()
    f.next()
    for i in xrange(N):
        line = f.next()
        words = line.split()[4:]
        for j in xrange(N):
            matrix[i,j] = float(words[j])
    f.close()
    return matrix


def _iter_symbols_numbers(numbers, N):
    for i in xrange(N):
        if numbers is None:
            number = 0
            symbol = "?"
        else:
            number = numbers[i]
            symbol = periodic[number].symbol
        yield i, symbol, number


def dump_atom_scalars(filename, scalars, numbers=None, name="Charge"):
    f = file(filename, "w")
    print >> f, "number of atoms: %i" % len(scalars)
    print >> f, "  i        Z  %s" % name.rjust(10)
    print >> f, "-----------------------------"
    for i, symbol, number in _iter_symbols_numbers(numbers, len(scalars)):
        print >> f, "% 3i  %2s  % 3i   % 10.5f" % (
            i+1, symbol, number, scalars[i]
        )
    print >> f, "-----------------------------"


def dump_atom_vectors(filename, vectors, numbers=None, name="Dipole"):
    name = name.rjust(10)
    f = file(filename, "w")
    print >> f, "number of atoms:", len(vectors)
    print >> f, "  i        Z %(name)s-X %(name)s-Y %(name)s-Z  %(name)s" % {"name": name}
    print >> f, "------------------------------------------------------------------"
    for i, symbol, number in _iter_symbols_numbers(numbers, len(vectors)):
        print >> f, "% 3i  %2s  % 3i   % 10.5f   % 10.5f   % 10.5f   % 10.5f" % (
            i+1, symbol, number, vectors[i,0], vectors[i,1],
            vectors[i,2], numpy.linalg.norm(vectors[i]),
        )
    f.close()


def dump_atom_matrix(filename, matrix, numbers=None, name="Matrix"):
    name = name.center(15)
    f = file(filename, "w")
    print >> f, "number of atoms:", len(matrix)
    print >> f, "%s | %s" % (name, " ".join(
        "  % 3i %2s  " % (key[0]+1, key[1]) for key
        in _iter_symbols_numbers(numbers, len(matrix))
    ))
    print >> f, "----------------+-"+"-"*(8+10*len(matrix))
    for i, symbol, number in _iter_symbols_numbers(numbers, len(matrix)):
        print >> f, "% 3i  %2s  % 3i    | %s" % (
            i+1, symbol, number, " ".join("% 10.5f" % val for val in matrix[i])
        )
    f.close()
