# HiPart is a software toolkit to analyse molecular densities with the hirshfeld partitioning scheme.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from molmod.units import angstrom
from molmod.transformations import random_rotation

import numpy, sys, os


__all__ = [
    "Error", "ProgressBar",
    "load_cube", "guess_density_type", "write_cube_in", "get_atom_grid",
    "compute_hirshfeld_weights"
]

class Error(Exception):
    pass


class ProgressBar(object):
    def __init__(self, label, n, f=sys.stdout):
        self.i = 0
        self.n = n
        self.f = f
        self.f.write("%s: " % label)
        self.f.flush()

    def __call__(self):
        if self.i % 10 == 0 or self.i == self.n:
            self.f.write(" %i%% " % ((100*self.i)/self.n))
        else:
            self.f.write(".")
        if self.i == self.n:
            self.f.write("\n")
        self.f.flush()
        self.i += 1


def load_cube(filename):
    f = file(filename)
    data = numpy.array([float(line.split()[3]) for line in f])
    if len(data) == 0:
        raise Error("Could not load cube file. File is empty.")
    f.close()
    return data


def guess_density_type(lot):
    if "mp2" in lot.lower():
        return "mp2"
    elif "mp3" in lot.lower():
        return "mp3"
    elif "mp4" in lot.lower():
        return "mp4"
    else:
        return "scf"


def write_cube_in(filename, grid_points):
    f = file(filename, "w")
    for point in grid_points:
        print >> f, "%15.10e %15.10e %15.10e" % tuple(point/angstrom)
    f.close()


def get_atom_grid(lebedev_xyz, center, radii):
    num_lebedev = len(lebedev_xyz)
    grid_points = numpy.zeros((num_lebedev*len(radii),3), float)
    counter = 0
    for r in radii:
        rot = random_rotation()
        grid_points[counter:counter+num_lebedev] = r*numpy.dot(lebedev_xyz,rot.r)+center
        counter += num_lebedev
    return grid_points


def compute_hirshfeld_weights(i, atom_fns, num_lebedev, distances):
    N = len(atom_fns)
    # construct the pro-atom and pro-molecule on this atomic grid
    pro_atom = numpy.array([atom_fns[i].density.y]*num_lebedev).transpose().ravel()
    pro_mol = numpy.zeros(len(pro_atom), float)
    for j in xrange(N):
        if i==j:
            pro_mol += pro_atom
        else:
            pro_mol += atom_fns[j].density(distances[(i,j)])

    # avoid division by zero
    pro_atom[pro_mol < 1e-40] = 1e-40
    pro_mol[pro_mol < 1e-40] = 1e-40
    # multiply the density on the grid by the weight function
    return pro_atom/pro_mol


