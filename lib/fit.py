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


from hipart.lebedev_laikov import get_grid
from hipart.core import ProgressBar
from hipart.tools import load_cube

from molmod.transformations import random_rotation
from molmod.units import angstrom

import os, numpy


__all__ = ["compute_mol_esp", "ChargeLSCostFunction", "ESPCostFunction"]


def compute_mol_esp(fchk, atom_table, workdir, esp_out_fn, num_lebedev, density_type, clonedir=None, scale_min=2.0, scale_max=15.0, scale_steps=40):
    """Wrapper function that computes grid and constructs ESP cost function"""
    grid_points, grid_weights, densities, potentials = _compute_mol_grid(
        fchk, atom_table, workdir, num_lebedev, density_type, clonedir,
        scale_min, scale_max, scale_steps
    )
    if clonedir is not None:
        other_potentials = load_cube(os.path.join(clone_dir, "molecule_potential.cube"), value_only=True)
        other_densities = load_cube(os.path.join(clone_dir, "molecule_density.cube"), value_only=True)
        mol_esp_cost = ESPCostFunction(fchk.molecule.coordinates, grid_points, grid_weights, other_densities, potentials, fchk.fields.get("Charge"))
        mol_esp_cost.write_matrices_to_file(esp_out_fn, other_potentials)
        mol_esp_cost.densities = densities
        mol_esp_cost.other_densities = other_densities
        mol_esp_cost.other_potentials = other_potentials
    else:
        mol_esp_cost = ESPCostFunction(fchk.molecule.coordinates, grid_points, grid_weights, densities, potentials, fchk.fields.get("Charge"))
        mol_esp_cost.write_matrices_to_file(esp_out_fn)
        mol_esp_cost.other_densities = None
        mol_esp_cost.other_potentials = None
    return mol_esp_cost


def _compute_mol_grid(fchk, atom_table, workdir, num_lebedev, density_type, clonedir=None, scale_min=1.5, scale_max=30.0, scale_steps=60):
    """Central function that computes grid and constructs ESP cost function"""
    lebedev_xyz, lebedev_weights = get_grid(num_lebedev)
    if not os.path.isfile(os.path.join(workdir, "molecule.points")):
        if clonedir is None:
            # we have to generate a new grid. The grid is constructed taking
            # into account the following considerations:
            # 1) Grid points within the cusp region are discarded
            # 2) The rest of the molecular and surrounding volume is sampled
            #    with spherical grids centered on the atoms. Around each atom,
            #    'scale_steps' of shells are placed with lebedev grid points
            #    (num_lebedev).
            # 3) The radii of the shells start from scale_min*(cusp_radius+0.2)
            #    and go up to scale_max*(cusp_radius+0.2).
            # 4) Each shell will be randomly rotated around the atom to avoid
            #    preferential directions in the grid.
            # 5) The default parameters for the grid should be sufficient for
            #    sane ESP fitting.

            scale_factor = (scale_max/scale_min)**(1.0/(scale_steps-1))
            scales = scale_min*scale_factor**numpy.arange(scale_steps)

            atom_table.init_cusp_cutoffs()
            atom_radii = numpy.array([atom_table.records[number].cusp_cutoff+0.2 for number in fchk.molecule.numbers])

            grid_points = []
            grid_weights = []
            pb = ProgressBar("Computing grid points", scale_steps)
            for scale in scales:
                pb()
                radii = scale*atom_radii
                #print radii/angstrom
                for i in xrange(fchk.molecule.size):
                    rot = random_rotation()
                    for j in xrange(len(lebedev_xyz)):
                        my_point = radii[i]*numpy.dot(rot.r, lebedev_xyz[j]) + fchk.molecule.coordinates[i]
                        distances = numpy.sqrt(((fchk.molecule.coordinates - my_point)**2).sum(axis=1))
                        if (distances < scales[0]*atom_radii).any():
                            continue
                        grid_points.append(my_point)
                        grid_weights.append(lebedev_weights[j])
            pb()
            grid_points = numpy.array(grid_points)
            grid_weights = numpy.array(grid_weights)

            # write the grid data to file
            f = file(os.path.join(workdir, "molecule.points"), "w")
            for cor in grid_points/angstrom:
                print >> f, "%s %s %s" % tuple(cor)
            f.close()
            f = file(os.path.join(workdir, "molecule.xyz"), "w")
            print >> f, len(grid_points)
            print >> f, "No interesting title today"
            for cor in grid_points/angstrom:
                print >> f, "X %s %s %s" % tuple(cor)
            f.close()
            grid_weights.tofile(os.path.join(workdir, "molecule.weights"), " ")
            del grid_points
            del grid_weights
        else:
            # Refer to another grid. Remove existing files first.
            print "Linking grid points"
            abs_clonedir = os.path.abspath(os.path.realpath(clonedir))
            for filename in ["molecule.points","molecule.weights","molecule.xyz"]:
                src = os.path.join(abs_clonedir, filename)
                dest = os.path.join(workdir, filename)
                if os.path.exists(dest):
                    os.path.remove(dest)
                os.symlink(src, dest)

    if not os.path.isfile(os.path.join(workdir, "molecule_density.cube")):
        print "Computing density molecule"
        os.system(". ~/g03.profile; cubegen 0 fdensity=%s %s %s -5 < %s" % (
            density_type, fchk.filename,
            os.path.join(workdir, "molecule_density.cube"),
            os.path.join(workdir, "molecule.points"),
        ))
        print "Computing potential molecule"
        os.system(". ~/g03.profile; cubegen 0 potential=%s %s %s -5 < %s" % (
            density_type, fchk.filename,
            os.path.join(workdir, "molecule_potential.cube"),
            os.path.join(workdir, "molecule.points"),
        ))


    print "Loading data"
    grid_weights = numpy.fromfile(os.path.join(workdir, "molecule.weights"), float, sep=" ")
    tmp = load_cube(os.path.join(workdir, "molecule_density.cube"))
    densities = tmp[:,-1]
    grid_points = tmp[:,:-1]
    potentials = load_cube(os.path.join(workdir, "molecule_potential.cube"), values_only=True)

    return grid_points, grid_weights, densities, potentials


class ChargeLSCostFunction(object):
    def __init__(self, design_matrix, expected_values, weights, total_charge):
        self.design_matrix = design_matrix*weights.reshape((-1,1))
        self.expected_values = expected_values*weights
        self.weights = weights
        self.weights_sqsum = (weights**2).sum()

        self.A = numpy.dot(self.design_matrix.transpose(), self.design_matrix)
        self.B = numpy.dot(self.design_matrix.transpose(), self.expected_values)
        self.C = numpy.dot(self.expected_values.transpose(), self.expected_values)
        self.rms = numpy.sqrt(self.C/self.weights_sqsum)

        Ap = numpy.zeros((len(self.A)+1,len(self.A)+1), float)
        Ap[:-1,:-1] = self.A
        Ap[:-1,-1] = 1
        Ap[-1,:-1] = 1
        Bp = numpy.zeros(len(self.B)+1, float)
        Bp[:-1] = self.B
        Bp[-1] = total_charge
        self.solution = numpy.linalg.solve(Ap, Bp)[:-1]

        A_evals = numpy.linalg.eigvalsh(self.A)
        self.condition_number = A_evals[-1]/A_evals[0]

    def rmsd(self, charges):
        return numpy.sqrt((numpy.dot(charges, numpy.dot(self.A, charges) - 2*self.B) + self.C)/self.weights_sqsum)

    def write_matrices_to_file(self, filename, other_expected_values=None):
        f = file(filename, "w")
        self.dump_matrices_to_file(f, other_expected_values)
        f.close()

    def dump_matrices_to_file(self, f, other_expected_values=None):
        print >> f, "A: Quadratic part of the ESP Cost function"
        for row in self.A:
            print >> f, " ".join("% 15.10e" % value for value in row)

        print >> f, "B: Linear part of the ESP Cost function"
        print >> f, " ".join("% 15.10e" % value for value in self.B)

        print >> f, "C: Constant term of the ESP Cost function"
        print >> f, "% 15.10e" % self.C

        print >> f, "W: Sum of the squared weights"
        print >> f, "% 15.10e" % self.weights_sqsum

        print >> f, "Ccross: Cross term from the reference ESP Cost function"
        if other_expected_values is None:
            Ccross = 0
        else:
            Ccross = numpy.dot(self.expected_values.transpose(), other_expected_values*self.weights)
        print >> f, "% 15.10e" % Ccross


class ESPCostFunction(ChargeLSCostFunction):
    def __init__(self, coordinates, grid_points, grid_weights, densities, potentials, total_charge):
        self.grid_points = grid_points
        self.grid_weights = grid_weights
        self.densities = densities
        self.potentials = potentials

        expected_values = potentials.copy()
        design_matrix = numpy.zeros((len(expected_values), len(coordinates)), float)
        self.distances = {}
        for i in xrange(len(coordinates)):
            distances = numpy.sqrt(((grid_points - coordinates[i])**2).sum(axis=1))
            self.distances[i] = distances
            design_matrix[:,i] = 1/distances
        weights = grid_weights * numpy.exp(-densities/1e-5)
        ChargeLSCostFunction.__init__(self, design_matrix, expected_values, weights, total_charge)


