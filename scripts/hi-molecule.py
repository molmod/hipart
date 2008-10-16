#!/usr/bin/python
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


from molmod.io.gaussian03.fchk import FCHKFile
from molmod.data.periodic import periodic
from molmod.units import angstrom

from hipart.core import *
from hipart.lebedev_laikov import get_grid, grid_fns
from hipart.tools import load_cube, guess_density, write_atom_grid
from hipart.fit import compute_mol_esp

from optparse import OptionParser
import os, shutil, copy, sys, numpy


parser = OptionParser("%prog [options] densities.txt gaussian.fchk")
parser.add_option(
    "--density",
    help="The density field to use from the gaussian fchk file (scf, mp2, mp3, "
    "...). If not given, the program will guess it based on the level of "
    "theory."
)
parser.add_option(
    "-l", "--lebedev", default=110, type='int',
    help="The number of grid points for the spherical averaging. "
    "[default=%default]. Select from: " + (", ".join(str(i) for i in sorted(grid_fns)))
)
parser.add_option(
    "--clone", default=None,
    help="Link the grid coordinates (and weights) from the given directory."
)
parser.add_option(
    "-n", "--no-fix-total-charge", dest="fix_total_charge", default="True",
    action="store_false", help="Do not correct the total charge."
)
parser.add_option(
    "-t", "--threshold", default=1e-4, type='float',
    help="When the maximum change in the charges drops below this threshold "
    "value, the iteration stops. [default=%default]"
)
(options, args) = parser.parse_args()
densities_filename, fchk_filename = args

rootdir = os.path.dirname(fchk_filename)
prefix = os.path.basename(fchk_filename).replace(".fchk", "")
workdir = os.path.join(rootdir, "%s-hipart-work" % prefix)
out_fn = os.path.join(rootdir, "%s-hipart.out" % prefix)
esp_out_fn = os.path.join(rootdir, "%s-esp.out" % prefix)

if os.path.isfile(out_fn):
    print "Nothing todo. Bye bye."
    sys.exit()


# A) Load FCHK file and atomic densities
fchk = FCHKFile(fchk_filename, field_labels=["Charge"])
at = AtomTable(densities_filename)

# A.1) Guess options.density
if options.density is None:
    options.density = guess_density(fchk.lot)

# B) Compute densities on atomic grids
if not os.path.isdir(workdir):
    os.mkdir(workdir)

# B.0) make the lebedev grid
num_lebedev = options.lebedev
lebedev_xyz, lebedev_weights = get_grid(num_lebedev)

# B.1) call cubegen a few times
pb = ProgressBar("Density on atomic grids", fchk.molecule.size)
for i, number in enumerate(fchk.molecule.numbers):
    pb()
    den_fn = os.path.join(workdir, "atom%05idens.cube" % i)
    den_fn_bin = "%s.bin" % den_fn
    if not os.path.isfile(den_fn_bin):
        center = fchk.molecule.coordinates[i]
        grid_prefix = os.path.join(workdir, "atom%05igrid" % i)
        write_atom_grid(grid_prefix, lebedev_xyz, center, at.records[number].rs)
        os.system(". ~/g03.profile; cubegen 0 fdensity=%s %s %s -5 < %s" % (
            options.density, fchk_filename, den_fn, "%s.txt" % grid_prefix
        ))
        tmp = load_cube(den_fn, values_only=True)
        tmp.tofile(den_fn_bin)
        os.remove("%s.txt" % grid_prefix)
        os.remove(den_fn)
pb()


# C) Precompute distances and load density data
distances = {}
densities = []
pb = ProgressBar("Precomputing distances", fchk.molecule.size**2)
for i, number_i in enumerate(fchk.molecule.numbers):
    den_fn_bin = os.path.join(workdir, "atom%05idens.cube.bin" % i)
    densities.append(numpy.fromfile(den_fn_bin, float))
    grid_fn_bin = os.path.join(workdir, "atom%05igrid.bin" % i)
    grid_points = numpy.fromfile(grid_fn_bin, float).reshape((-1,3))

    for j, number_j in enumerate(fchk.molecule.numbers):
        pb()
        if i!=j:
            distances[(i,j)] = numpy.sqrt(((grid_points - fchk.molecule.coordinates[j])**2).sum(axis=1))
pb()


# D) compute (iterative) hirshfeld charges
counter = 0
old_charges = numpy.zeros(fchk.molecule.size, float)
while True:
    charges = []

    # construct the pro-atom density functions, using the densities from the
    # previous iteration
    atom_fns = []
    for i, number_i in enumerate(fchk.molecule.numbers):
        atom_fns.append(at.records[number_i].get_atom_fn(old_charges[i]))

    # Run over each atom and ...
    for i, number_i in enumerate(fchk.molecule.numbers):
        # construct the pro-atom and pro-molecule on this atomic grid
        atom_weights = numpy.array([atom_fns[i].density.y]*num_lebedev).transpose().ravel()
        promol_weights = numpy.zeros(len(atom_weights), float)
        for j, number_j in enumerate(fchk.molecule.numbers):
            if i==j:
                promol_weights += atom_weights
            else:
                promol_weights += atom_fns[j].density(distances[(i,j)])

        # avoid division by zero
        atom_weights[promol_weights < 1e-40] = 1e-40
        promol_weights[promol_weights < 1e-40] = 1e-40
        # multiply the density on the grid by the weight function
        fn = densities[i]*atom_weights/promol_weights

        # integrate over the spherical degrees of freedom, using lebedev and
        # numpy tricks
        radial_int = (fn.reshape((-1, num_lebedev))*lebedev_weights).sum(axis=1)*4*numpy.pi
        # integrate over the radial axis, using our simple integration routine
        rs = at.records[number_i].rs
        total_int = integrate(rs, radial_int*rs**2)
        # append the integral to the list of charges
        charges.append(number_i - total_int)
    # ordinary blablabla below ...
    charges = numpy.array(charges)
    if counter == 0:
        charges0 = charges
    max_change = abs(charges-old_charges).max()
    print "Hirshfeld (%03i)    max change = %10.5e    total charge = %10.5e" % (
        counter, max_change, charges.sum()
    )
    if max_change < options.threshold:
        break
    counter += 1
    old_charges = charges


if options.fix_total_charge:
    print "Ugly step: Adding constant to charges so that the total charge is zero."
    charges0 -= (charges0.sum() - fchk.fields.get("Charge"))/fchk.molecule.size
    charges -= (charges.sum() - fchk.fields.get("Charge"))/fchk.molecule.size


# print some nice output
f = file(out_fn, "w")
print >> f, "  i        Z      Hirsh          Hirsh-I"
print >> f, "--------------------------------------------"
for i, number in enumerate(fchk.molecule.numbers):
    print >> f, "% 3i  %2s  % 3i   % 10.5f     % 10.5f" % (
        i+1, periodic[number].symbol, number, charges0[i], charges[i]
    )
print >> f, "--------------------------------------------"
print >> f, "   Q SUM       % 10.5f     % 10.5f" % (charges0.sum(), charges.sum())
print >> f, "   Q RMS       % 10.5f     % 10.5f" % (
    numpy.sqrt((charges0**2).mean()),
    numpy.sqrt((charges**2).mean()),
)
mol_esp_cost = compute_mol_esp(fchk, at, workdir, esp_out_fn, options.lebedev, options.density, options.clone)
print >> f, " ESP RMS         % 10.5e   % 10.5e" % (
    mol_esp_cost.rms,
    mol_esp_cost.rms,
)
print >> f, "ESP RMSD         % 10.5e   % 10.5e" % (
    mol_esp_cost.rmsd(charges0),
    mol_esp_cost.rmsd(charges),
)
print >> f
f.close()

