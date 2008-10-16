# HiPart is a tool to analyse molecular densities with the hirshfeld partitioning scheme
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

import numpy


__all__ = ["load_cube", "guess_density", "write_atom_grid"]


def load_cube(filename, values_only=False):
    f = file(filename)
    if values_only:
        tmp = numpy.array([float(line.split()[3]) for line in f])
        if len(tmp) == 0:
            raise Error("Could not load cube file. File is empty.")
    else:
        tmp = [
            [float(word) for word in line.split()]
            for line in f
        ]
        tmp = numpy.array([row for row in tmp if len(row) == 4])
        if len(tmp) == 0:
            raise Error("Could not load cube file. File is empty.")
        tmp[:,:3] *= angstrom
    f.close()
    return tmp


def guess_density(lot):
    if "mp2" in lot.lower():
        return "mp2"
    elif "mp3" in lot.lower():
        return "mp3"
    elif "mp4" in lot.lower():
        return "mp4"
    else:
        return "scf"


def write_atom_grid(prefix, lebedev_xyz, center, radii):
    f = file("%s.txt" % prefix, "w")
    grid_points = numpy.zeros((len(lebedev_xyz)*len(radii),3), float)
    counter = 0
    for r in radii:
        rot = random_rotation()
        for xyz in lebedev_xyz:
            point = r*numpy.dot(rot.r,xyz)+center
            grid_points[counter] = point
            print >> f, "%15.10e %15.10e %15.10e" % tuple(point/angstrom)
            counter += 1
    f.close()
    grid_points.tofile("%s.bin" % prefix)

