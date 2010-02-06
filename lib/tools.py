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


from molmod import angstrom, Rotation
from molmod.periodic import periodic

import numpy, sys, os


__all__ = [
    "Error", "load_cube", "guess_density_type", "write_cube_in",
    "cubegen_density", "cubegen_potential", "cubegen_orbital",
    "get_atom_grid", "compute_stockholder_weights",
    "load_charges", "load_dipoles", "dump_charges"
]


class Error(Exception):
    pass


def load_cube(filename, N):
    f = file(filename)
    data = numpy.zeros(N, float)
    counter = 0
    for line in f:
        if counter == N:
            raise Error("Cube file has wrong size. Expecting %i points. Got more." % (N, counter))
        line = line.strip()
        data[counter] = float(line[line.rfind(" "):])
        counter += 1
    if counter != N:
        raise Error("Cube file has wrong size. Expecting %i points. Got %i." % (N, counter))
    if numpy.isnan(data[0]): # ugly workarond for stupid cubegen
        data[0] = data[1]
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

def cubegen_density(grid_fn, den_fn, fchk, density_type, grid_size):
    if fchk.lot.startswith("ROHF"):
        # ugly hack: Workaround for stupid cubegen that does not work
        # property on ROHF calculations. pfff...
        densities = 0.0
        num_elec = fchk.fields.get("Number of electrons")
        for j in xrange((num_elec+1)/2):
            orb_fn = den_fn.replace("dens", "orb%05i" % j)
            orb_fn_bin = "%s.bin" % orb_fn
            os.system(". ~/g03.profile; cubegen 0 MO=%s %s %s -5 < %s" % (
                j+1, fchk.filename, orb_fn, grid_fn
            ))
            orb = load_cube(orb_fn, grid_size)
            orb.tofile(orb_fn_bin)
            if 2*j+1==num_elec:
                densities += orb**2
            else:
                densities += 2*orb**2
    elif fchk.lot.startswith("RO"):
        raise ComputeError("Can not cope with RO calculation, except ROHF. Cubegen can not compute the density properly for RO calculations!")
    else:
        os.system(". ~/g03.profile; cubegen 0 fdensity=%s %s %s -5 < %s" % (
            density_type, fchk.filename, den_fn, grid_fn
        ))
        densities = load_cube(den_fn, grid_size)

    return densities

def cubegen_potential(grid_fn, pot_fn, fchk, density_type, grid_size):
    os.system(". ~/g03.profile; cubegen 0 potential=%s %s %s -5 < %s" % (
        density_type, fchk.filename, pot_fn, grid_fn,
    ))
    return load_cube(pot_fn, grid_size)

def cubegen_orbital(grid_fn, orb_fn, fchk, orb_index, grid_size):
    os.system(". ~/g03.profile; cubegen 0 MO=%i %s %s -5 < %s" % (
        orb_index+1, fchk.filename, orb_fn, grid_fn,
    ))
    return load_cube(orb_fn, grid_size)


def get_atom_grid(lebedev_xyz, center, radii):
    num_lebedev = len(lebedev_xyz)
    grid_points = numpy.zeros((num_lebedev*len(radii),3), float)
    counter = 0
    for r in radii:
        rot = Rotation.random()
        grid_points[counter:counter+num_lebedev] = r*numpy.dot(lebedev_xyz,rot.r)+center
        counter += num_lebedev
    return grid_points


def compute_stockholder_weights(i, atom_fns, num_lebedev, distances, k=None):
    """ Evaulates the weights of atom i in the grid of atom k.

    If k is not given, it defaults to i.
    """
    N = len(atom_fns)
    # construct the pro-atom and pro-molecule on this atomic grid
    if k is None:
        k = i

    if k == i:
        pro_atom = numpy.array([atom_fns[i].density.y]*num_lebedev).transpose().ravel()
    else:
        pro_atom = atom_fns[i].density(distances[(k,i)])

    pro_mol = numpy.zeros(len(pro_atom), float)
    for j in xrange(N):
        if i==j:
            pro_mol += pro_atom
        else:
            if k==j:
                pro_mol += numpy.array([atom_fns[j].density.y]*num_lebedev).transpose().ravel()
            else:
                pro_mol += atom_fns[j].density(distances[(k,j)])

    # avoid division by zero
    pro_atom[pro_mol < 1e-40] = 1e-40
    pro_mol[pro_mol < 1e-40] = 1e-40
    # multiply the density on the grid by the weight function
    return pro_atom/pro_mol


def load_charges(filename):
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


def load_dipoles(filename):
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


def dump_charges(filename, charges, numbers=None):
    f = file(filename, "w")
    print >> f, "number of atoms: %i" % len(charges)
    print >> f, "  i        Z      Charge"
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



