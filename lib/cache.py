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


from hipart.tools import ProgressBar, write_cube_in, load_cube, get_atom_grid, \
    compute_hirshfeld_weights
from hipart.integrate import integrate_log, integrate_lebedev

from molmod.transformations import random_rotation
from molmod.data.periodic import periodic

import os, numpy


__all__ = ["Error", "Cache"]


class ComputeError(Exception):
    pass


class Cache(object):
    def __init__(self, context):
        self.context = context

    def do_atom_grids(self):
        if hasattr(self, "atom_grid_distances"):
            return

        molecule = self.context.fchk.molecule
        workdir = self.context.workdir

        self.atom_grid_distances = {}
        pb = ProgressBar("Atomic grids (and distances)", molecule.size**2)
        for i, number_i in enumerate(molecule.numbers):
            grid_prefix = os.path.join(workdir, "atom%05igrid" % i)
            grid_fn_bin = "%s.hipart.bin" % grid_prefix
            grid_fn_txt = "%s.txt" % grid_prefix
            if os.path.isfile(grid_fn_bin):
                grid_points = numpy.fromfile(grid_fn_bin).reshape((-1,3))
            else:
                center = molecule.coordinates[i]
                grid_points = get_atom_grid(
                    self.context.lebedev_xyz, center,
                    self.context.atom_table.records[number_i].rs
                )
                grid_points.tofile(grid_fn_bin)
            if not os.path.isfile(grid_fn_txt):
                write_cube_in(grid_fn_txt, grid_points)

            for j, number_j in enumerate(molecule.numbers):
                pb()
                if i!=j:
                    distances = numpy.sqrt(((grid_points - molecule.coordinates[j])**2).sum(axis=1))
                    self.atom_grid_distances[(i,j)] = distances
        pb()

    def do_atom_densities(self):
        if hasattr(self, "atom_densities"):
            return
        self.do_atom_grids()
        molecule = self.context.fchk.molecule
        workdir = self.context.workdir

        pb = ProgressBar("Molecular density on atomic grids", molecule.size)
        self.atom_densities = []
        for i in xrange(molecule.size):
            pb()
            den_fn = os.path.join(workdir, "atom%05idens.cube" % i)
            den_fn_bin = "%s.bin" % den_fn
            if os.path.isfile(den_fn_bin):
                densities = numpy.fromfile(den_fn_bin, float)
            else:
                if not os.path.isfile(den_fn):
                    grid_fn = os.path.join(workdir, "atom%05igrid.txt" % i)
                    os.system(". ~/g03.profile; cubegen 0 fdensity=%s %s %s -5 < %s" % (
                        self.context.options.density,
                        self.context.fchk.filename, den_fn, grid_fn
                    ))
                densities = load_cube(den_fn)
                densities.tofile(den_fn_bin)
            self.atom_densities.append(densities)
        pb()

    def do_iterative_hirshfeld(self):
        if hasattr(self, "hirshi_charges") and hasattr(self, "hirshi_weights"):
            return
        print "Iterative hirshfeld"

        molecule = self.context.fchk.molecule
        atom_table = self.context.atom_table

        hirshi_charges_fn_bin = os.path.join(self.context.workdir, "hirshi_charges.out.bin")
        hirshi_weights_tpl_bin = os.path.join(self.context.workdir, "hirshi_weights%05i.hipart.bin")

        # Try to read the data from the workdir
        if os.path.isfile(hirshi_charges_fn_bin):
            self.hirshi_charges = numpy.fromfile(hirshi_charges_fn_bin, float)
            self.hirshi_weights = []

            atom_fns = [] # construct the atom fns, just in case...
            def construct_atom_fns():
                if len(atom_fns) > 0:
                    return
                self.do_atom_grids()
                for i, number_i in enumerate(molecule.numbers):
                    atom_fns.append(atom_table.records[number_i].get_atom_fn(self.hirshi_charges[i]))

            for i in xrange(molecule.size):
                if os.path.isfile(hirshi_weights_tpl_bin % i):
                    hw = numpy.fromfile(hirshi_weights_tpl_bin % i, float)
                else:
                    construct_atom_fns()
                    hw = compute_hirshfeld_weights(
                        i, atom_fns, self.context.num_lebedev,
                        self.atom_grid_distances
                    )
                    hw.tofile(hirshi_weights_tpl_bin % i)
                self.hirshi_weights.append(hw)
            print "Read iterative hirsfeld results from workdir."
        else:
            self.do_atom_densities()
            counter = 0
            old_charges = numpy.zeros(molecule.size, float)
            while True:
                # construct the pro-atom density functions, using the densities
                # from the previous iteration.
                atom_fns = []
                for i, number_i in enumerate(molecule.numbers):
                    atom_fns.append(atom_table.records[number_i].get_atom_fn(old_charges[i]))

                # compute the hirshfeld charges
                charges = []
                hirshfeld_weights = []
                for i, number_i in enumerate(molecule.numbers):
                    hw = compute_hirshfeld_weights(
                        i, atom_fns, self.context.num_lebedev,
                        self.atom_grid_distances
                    )
                    fn = self.atom_densities[i]*hw
                    radfun = integrate_lebedev(self.context.lebedev_weights, fn)
                    rs = atom_table.records[number_i].rs
                    num_electrons = integrate_log(rs, radfun*rs**2)

                    charges.append(number_i - num_electrons)
                    hirshfeld_weights.append(hw)

                # ordinary blablabla below ...
                charges = numpy.array(charges)
                max_change = abs(charges-old_charges).max()
                print "Iteration %03i    max change = %10.5e    total charge = %10.5e" % (
                    counter, max_change, charges.sum()
                )
                if max_change < self.context.options.threshold:
                    break
                counter += 1
                old_charges = charges

            if self.context.options.fix_total_charge:
                print "Ugly step: Adding constant to charges so that the total charge is zero."
                charges -= (charges.sum() - self.context.fchk.fields.get("Charge"))/molecule.size

            self.hirshi_charges = charges
            self.hirshi_weights = hirshfeld_weights

            print "Writing iterative hirshfeld charges to workdir"
            self.hirshi_charges.tofile(hirshi_charges_fn_bin)
            for i, hw in enumerate(self.hirshi_weights):
                hw.tofile(hirshi_weights_tpl_bin % i)

        # now some nice output
        self.do_esp_costfunction()
        f = file(os.path.join(self.context.outdir, "hirshi_charges.txt"), "w")
        print >> f, "  i        Z      Charge"
        print >> f, "-----------------------------"
        for i, number in enumerate(molecule.numbers):
            print >> f, "% 3i  %2s  % 3i   % 10.5f" % (
                i+1, periodic[number].symbol, number, self.hirshi_charges[i]
            )
        print >> f, "-----------------------------"
        print >> f, "   Q SUM       % 10.5f" % self.hirshi_charges.sum()
        print >> f, "   Q RMS       % 10.5f" % numpy.sqrt((self.hirshi_charges**2).mean())
        print >> f, " ESP RMS         % 10.5e" % self.mol_esp_cost.rms
        print >> f, "ESP RMSD         % 10.5e" % self.mol_esp_cost.rmsd(self.hirshi_charges)
        print >> f
        f.close()


    def do_esp_costfunction(self):
        clonedir = self.context.options.clone
        if clonedir is not None:
            clone_dens_fn = os.path.join(clone_dir, "molecule_dens.cube.bin")
            if not os.path.isfile(clone_dens_fn):
                raise ComputeError("Could not find reference data: %s" % clone_dens_fn)
            clone_pot_fn = os.path.join(clone_dir, "molecule_pot.cube.bin")
            if not os.path.isfile(clone_pot_fn):
                raise ComputeError("Could not find reference data: %s" % clone_pot_fn)

        self._do_molecular_grid()

        total_charge = self.context.fchk.fields.get("Charge")
        coordinates = self.context.fchk.molecule.coordinates
        if clonedir is not None:
            from hipart.fit import ESPCrossCostFunction
            other_densities = numpy.fromfile(clone_dens_fn)
            other_potentials = numpy.fromfile(clone_pot_fn)
            self.mol_esp_cost = ESPCrossCostFunction(
                coordinates, self.mol_points, self.mol_weights,
                self.mol_densities, self.mol_potentials,
                other_densities, other_potentials, total_charge,
            )
        else:
            from hipart.fit import ESPCostFunction
            self.mol_esp_cost = ESPCostFunction(
                coordinates, self.mol_points, self.mol_weights,
                self.mol_densities, self.mol_potentials, total_charge,
            )

        outfn = os.path.join(self.context.outdir, "mol_esp_cost.txt")
        self.mol_esp_cost.write_matrices_to_file(outfn)

    def _do_molecular_grid(self):
        lebedev_xyz = self.context.lebedev_xyz
        lebedev_weights = self.context.lebedev_weights
        molecule = self.context.fchk.molecule
        clonedir = self.context.options.clone
        workdir = self.context.workdir
        atom_table = self.context.atom_table

        dens_fn = os.path.join(workdir, "molecule_dens.cube")
        pot_fn = os.path.join(workdir, "molecule_pot.cube")
        dens_fn_bin = "%s.bin" % dens_fn
        pot_fn_bin = "%s.bin" % pot_fn
        points_fn = os.path.join(workdir, "molecule_points.txt")
        points_fn_bin = os.path.join(workdir, "molecule_points.bin")
        weights_fn_bin = os.path.join(workdir, "molecule_weights.bin")

        if not (os.path.isfile(dens_fn_bin) and os.path.isfile(pot_fn_bin) and
                os.path.isfile(points_fn_bin) and os.path.isfile(weights_fn_bin)):
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

                scale_min = 1.5
                scale_max = 30.0
                scale_steps = 30
                scale_factor = (scale_max/scale_min)**(1.0/(scale_steps-1))
                scales = scale_min*scale_factor**numpy.arange(scale_steps)

                atom_table.init_cusp_cutoffs()
                atom_radii = numpy.array([atom_table.records[number].cusp_cutoff+0.2 for number in molecule.numbers])

                grid_points = []
                grid_weights = []
                pb = ProgressBar("Molecular grid", scale_steps)
                for scale in scales:
                    pb()
                    radii = scale*atom_radii
                    for i in xrange(molecule.size):
                        rot = random_rotation()
                        for j in xrange(len(lebedev_xyz)):
                            my_point = radii[i]*numpy.dot(rot.r, lebedev_xyz[j]) + molecule.coordinates[i]
                            distances = numpy.sqrt(((molecule.coordinates - my_point)**2).sum(axis=1))
                            if (distances < scales[0]*atom_radii).any():
                                continue
                            grid_points.append(my_point)
                            grid_weights.append(lebedev_weights[j])
                pb()
                grid_points = numpy.array(grid_points)
                grid_weights = numpy.array(grid_weights)

                # write the grid data to file
                grid_points.tofile(points_fn_bin)
                grid_weights.tofile(weights_fn_bin)
            else:
                # Refer to another grid. Remove existing files first.
                print "Linking grid from clonedir"
                abs_clonedir = os.path.abspath(os.path.realpath(clonedir))
                for filename in ["molecule_points.bin","molecule_weights.bin"]:
                    src = os.path.join(abs_clonedir, filename)
                    dest = os.path.join(workdir, filename)
                    if os.path.isfile(dest):
                        os.path.remove(dest)
                    os.symlink(src, dest)
                grid_points = numpy.fromfile(points_fn_bin, float).reshape((-1,3))
                grid_weights = numpy.fromfile(weights_fn_bin, float)

            write_cube_in(points_fn, grid_points) # prepare for cubegen
            print "Molecular density on moleculer grid"
            os.system(". ~/g03.profile; cubegen 0 fdensity=%s %s %s -5 < %s" % (
                self.context.options.density, self.context.fchk.filename,
                dens_fn, points_fn,
            ))
            densities = load_cube(dens_fn)
            densities.tofile(dens_fn_bin)
            print "Molecular potential on moleculer grid"
            os.system(". ~/g03.profile; cubegen 0 potential=%s %s %s -5 < %s" % (
                self.context.options.density, self.context.fchk.filename,
                pot_fn, points_fn,
            ))
            potentials = load_cube(pot_fn)
            potentials.tofile(pot_fn_bin)
        else:
            print "Loading molecular data from workdir"
            grid_points = numpy.fromfile(points_fn_bin, float).reshape((-1,3))
            grid_weights = numpy.fromfile(weights_fn_bin, float)
            densities = numpy.fromfile(dens_fn_bin, float)
            potentials = numpy.fromfile(pot_fn_bin, float)

        self.mol_points = grid_points
        self.mol_weights = grid_weights
        self.mol_densities = densities
        self.mol_potentials = potentials

    def do_atom_orbitals(self):
        if hasattr(self, "atom_orbitals"):
            return
        self.do_atom_grids()
        molecule = self.context.fchk.molecule
        num_orbitals = self.context.fchk.fields.get("Number of basis functions")
        workdir = self.context.workdir

        pb = ProgressBar("Molecular orbitals on atomic grids", molecule.size*num_orbitals)
        self.atom_orbitals = []
        for i in xrange(molecule.size):
            orbitals = []
            self.atom_orbitals.append(orbitals)
            grid_fn = os.path.join(workdir, "atom%05igrid.txt" % i)
            for j in xrange(num_orbitals):
                pb()
                cube_fn = os.path.join(workdir, "atom%05iorb%05i.cube" % (i,j))
                cube_fn_bin = "%s.bin" % cube_fn
                if os.path.isfile(cube_fn_bin):
                    wavefn = numpy.fromfile(cube_fn_bin, float)
                else:
                    if not os.path.isfile(cube_fn):
                        os.system(". ~/g03.profile; cubegen 0 MO=%i %s %s -5 < %s" % (
                            j+1, self.context.fchk.filename, cube_fn, grid_fn,
                        ))
                    wavefn = load_cube(cube_fn)
                    wavefn.tofile(cube_fn_bin)
                orbitals.append(wavefn)
        pb()

    def do_atom_matrices(self):
        if hasattr(self, "atom_matrices"):
            return
        self.do_atom_orbitals()
        self.do_iterative_hirshfeld()
        num_orbitals = self.context.fchk.fields.get("Number of basis functions")
        molecule = self.context.fchk.molecule

        pb = ProgressBar("Atom centered matrices", molecule.size)
        self.atom_matrices = []
        for i, number_i in enumerate(molecule.numbers):
            pb()
            orbitals = self.atom_orbitals[i]
            hw = self.hirshi_weights[i]
            matrix = numpy.zeros((num_orbitals,num_orbitals), float)
            self.atom_matrices.append(matrix)
            for j1 in xrange(num_orbitals):
                for j2 in xrange(j1+1):
                    integrand = orbitals[j1]*orbitals[j2]*hw
                    radfun = integrate_lebedev(self.context.lebedev_weights, integrand)
                    rs = self.context.atom_table.records[number_i].rs
                    value = integrate_log(rs, radfun*rs**2)
                    matrix[j1,j2] = value
                    matrix[j2,j1] = value
        pb()

        f = file(os.path.join(self.context.outdir, "atom_matrices.txt"), "w")
        print >> f, "number of orbitas:", num_orbitals
        print >> f, "number of atoms: ", molecule.size
        for i, number_i in enumerate(molecule.numbers):
            print >> f, "Atom %i: %s" % (i, periodic[number_i].symbol)
            for row in self.atom_matrices[i]:
                print >> f, " ".join("% 15.10e" % value for value in row)
        f.close()


