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
from hipart.fit import *
from hipart.lebedev_laikov import grid_fns
from hipart.tools import guess_density

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
    "-t", "--threshold", default=1e-4, type='float',
    help="When the maximum change in the charges drops below this threshold "
    "value, the iteration stops. [default=%default]"
)
(options, args) = parser.parse_args()
densities_filename, fchk_filename = args


prefix = os.path.basename(fchk_filename).replace(".fchk", "")
rootdir = os.path.dirname(fchk_filename)
workdir = os.path.join(rootdir, "%s-hipart-work" % prefix)
hipart_esp_out_fn = os.path.join(rootdir, "%s-hipart-esp.out" % prefix)
esp_out_fn = os.path.join(rootdir, "%s-esp.out" % prefix)

if os.path.isfile(hipart_esp_out_fn) and os.path.isfile(esp_out_fn):
    print "Nothing todo. Bye bye."
    sys.exit()


# A) Load FCHK file and atomic densities
fchk = FCHKFile(fchk_filename, field_labels=["Charge"])
at = AtomTable(densities_filename)

# A.1) Guess options.density
if options.density is None:
    options.density = guess_density(fchk.lot)


# B) Compute density and the esp
if not os.path.isdir(workdir):
    os.mkdir(workdir)

mol_esp_cost = compute_mol_esp(fchk, at, workdir, esp_out_fn, options.lebedev, options.density, options.clone)
print "Conventional ESP condition number: % 10.5e" % mol_esp_cost.condition_number

if not os.path.isfile(hipart_esp_out_fn):
    distances = mol_esp_cost.distances

    counter = 0
    old_charges = numpy.zeros(fchk.molecule.size, float)

    original_expected_values = mol_esp_cost.potentials.copy()
    for i, number in enumerate(fchk.molecule.numbers):
        original_expected_values -= number/distances[i]

    design_matrix = numpy.zeros((len(mol_esp_cost.grid_points), fchk.molecule.size), float)
    model_density = numpy.zeros(len(mol_esp_cost.grid_points), float)

    def tail(x):
        return numpy.exp(-abs(x))
        #return 1/(1+abs(x))

    while True:
        design_matrix[:] = 0.0
        model_density[:] = 0

        for i, number_i in enumerate(fchk.molecule.numbers):
            atom_fn = at.records[number_i].get_atom_fn(old_charges[i])
            model_density += atom_fn.density(distances[i])
            design_matrix[:,i] = atom_fn.potential(distances[i])/atom_fn.num_elec
            #design_matrix[:,i] = 1/distances[i]

        my_weights = mol_esp_cost.grid_weights.copy()
        my_weights *= tail((model_density - mol_esp_cost.densities)/1e-1)
        my_weights *= tail(mol_esp_cost.densities/1e-2)
        expected_values = original_expected_values
        expected_values = numpy.dot(design_matrix, fchk.molecule.numbers) - expected_values

        cost_fn = ChargeLSCostFunction(design_matrix, expected_values, my_weights, fchk.fields.get("Charge"))
        charges = cost_fn.solution

        if counter == 0:
            charges0 = charges
        max_change = abs(charges-old_charges).max()
        print "Hirshfeld-ESP (%03i)    max change = %10.5e    condition number = %10.5e" % (
            counter, max_change, cost_fn.condition_number
        )
        if max_change < options.threshold:
            break
        counter += 1
        old_charges = charges

    # To visualize the most relevant part of the grid ...
    #mask = (my_weights > my_weights.max()*0.01)
    #f = file(os.path.join(workdir, "molecule_relevant.xyz"), "w")
    #print >> f, sum(mask)
    #print >> f, "No interesting title today"
    #for cor in mol_esp_cost.grid_points[mask]/angstrom:
    #    print >> f, "X %s %s %s" % tuple(cor)
    #f.close()

    # Write output file
    f = file(hipart_esp_out_fn, "w")
    print >> f, "  i        Z      Hirsh-ESP      Hirsh-I-ESP"
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
    print >> f, " ESP RMS         % 10.5e   % 10.5e" % (
        mol_esp_cost.rms,
        mol_esp_cost.rms,
    )
    print >> f, "ESP RMSD         % 10.5e   % 10.5e" % (
        mol_esp_cost.rmsd(charges0),
        mol_esp_cost.rmsd(charges),
    )
    print >> f
    print >> f
    cost_fn.dump_matrices_to_file(f, mol_esp_cost.other_potentials)
    f.close()

