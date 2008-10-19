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
from hipart.fit import ESPCostFunction

from molmod.transformations import random_rotation
from molmod.data.periodic import periodic

import os, numpy


__all__ = ["Error", "Cache"]


class ComputeError(Exception):
    pass


class Cache(object):
    def __init__(self, context):
        self.context = context
        if self.context.reference is not None:
            self.reference = self.context.reference.cache
        else:
            self.reference = None

    def work_link(self, filename):
        if self.reference is None:
            raise Error("Can only link if there is a reference state.")
        if not os.path.isfile(os.path.join(self.context.workdir, filename)):
            workdir = os.path.abspath(self.context.workdir)
            ref_workdir = os.path.abspath(self.reference.context.workdir)
            common = os.path.commonprefix([workdir, ref_workdir])
            l = len(common)
            l = common.rfind("/")+1
            workdir = workdir[l:]
            ref_workdir = ref_workdir[l:]
            relpath = os.path.normpath(os.path.join("/".join(".." for i in xrange(workdir.count("/")+1)), ref_workdir))
            os.symlink(
                os.path.join(relpath, filename),
                os.path.join(self.context.workdir, filename),
            )

    def do_atom_grids(self):
        if hasattr(self, "atom_grid_distances") and hasattr(self, "atom_grid_points"):
            return

        molecule = self.context.fchk.molecule
        if self.reference is not None:
            self.reference.do_atom_grids()
            self.atom_grid_distances = self.reference.atom_grid_distances
            self.atom_grid_points = self.reference.atom_grid_points
            for i in xrange(molecule.size):
                self.work_link("atom%05igrid.txt" % i)
        else:
            workdir = self.context.workdir

            self.atom_grid_distances = {}
            self.atom_grid_points = []
            pb = ProgressBar("Atomic grids (and distances)", molecule.size**2)
            for i, number_i in enumerate(molecule.numbers):
                grid_prefix = os.path.join(workdir, "atom%05igrid" % i)
                grid_fn_bin = "%s.bin" % grid_prefix
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

                self.atom_grid_points.append(grid_points)
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

        hirshi_charges_fn_bin = os.path.join(self.context.workdir, "hirshi_charges.bin")
        hirshi_weights_tpl_bin = os.path.join(self.context.workdir, "hirshi_weights%05i.bin")

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

        self.do_esp_costfunction()

        # now some nice output
        def output(filename, charges, esp_cost):
            f = file(os.path.join(self.context.outdir, filename), "w")
            print >> f, "  i        Z      Charge"
            print >> f, "-----------------------------"
            for i, number in enumerate(molecule.numbers):
                print >> f, "% 3i  %2s  % 3i   % 10.5f" % (
                    i+1, periodic[number].symbol, number, charges[i]
                )
            print >> f, "-----------------------------"
            print >> f, "   Q SUM       % 10.5f" % charges.sum()
            print >> f, "   Q RMS       % 10.5f" % numpy.sqrt((charges**2).mean())
            print >> f, " ESP RMS         % 10.5e" % esp_cost.rms
            print >> f, "ESP RMSD         % 10.5e" % esp_cost.rmsd(charges)
            print >> f
            f.close()

        output("hirshi_charges.txt", self.hirshi_charges, self.mol_esp_cost)

        if self.reference is not None:
            relative_charges = self.hirshi_charges - self.reference.hirshi_charges
            output("hirshi_charges_relative.txt", relative_charges, self.relative_mol_esp_cost)

    def do_esp_costfunction(self):
        if hasattr(self, "mol_esp_cost"):
            return
        if self.reference is not None:
            self.reference.do_esp_costfunction()

        self._do_molecular_grid()

        total_charge = self.context.fchk.fields.get("Charge")
        coordinates = self.context.fchk.molecule.coordinates
        self.mol_esp_cost = ESPCostFunction(
            coordinates, self.mol_points, self.mol_weights,
            self.mol_densities, self.mol_potentials, total_charge,
        )
        outfn = os.path.join(self.context.outdir, "mol_esp_cost.txt")
        self.mol_esp_cost.write_to_file(outfn)
        if self.reference is not None:
            relative_total_charge = total_charge - self.reference.context.fchk.fields.get("Charge")
            self.relative_mol_esp_cost = ESPCostFunction(
                coordinates, self.mol_points, self.mol_weights,
                self.reference.mol_densities,
                self.mol_potentials - self.reference.mol_potentials,
                relative_total_charge
            )
            outfn = os.path.join(self.context.outdir, "mol_esp_cost_relative.txt")
            self.mol_esp_cost.write_to_file(outfn)


    def _do_molecular_grid(self):
        if hasattr(self, "mol_points") and hasattr(self, "mol_weights") and \
           hasattr(self, "mol_densities") and hasattr(self, "mol_potentials"):
            return

        lebedev_xyz = self.context.lebedev_xyz
        lebedev_weights = self.context.lebedev_weights
        molecule = self.context.fchk.molecule
        workdir = self.context.workdir
        atom_table = self.context.atom_table

        dens_fn = os.path.join(workdir, "molecule_dens.cube")
        pot_fn = os.path.join(workdir, "molecule_pot.cube")
        dens_fn_bin = "%s.bin" % dens_fn
        pot_fn_bin = "%s.bin" % pot_fn
        points_fn = os.path.join(workdir, "molecule_points.txt")
        points_fn_bin = os.path.join(workdir, "molecule_points.bin")
        weights_fn_bin = os.path.join(workdir, "molecule_weights.bin")


        if not (os.path.isfile(dens_fn_bin) and os.path.isfile(pot_fn_bin)):
            if self.reference is None:
                if not (os.path.isfile(points_fn_bin) and os.path.isfile(weights_fn_bin)):
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
                    grid_points = numpy.fromfile(points_fn_bin, float).reshape((-1,3))
                    grid_weights = numpy.fromfile(weights_fn_bin, float)
            else:
                print "Linking grid from reference state"
                grid_points = self.reference.mol_points
                grid_weights = self.reference.mol_weights
                self.work_link("molecule_points.bin")
                self.work_link("molecule_weights.bin")

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
            densities = numpy.fromfile(dens_fn_bin, float)
            potentials = numpy.fromfile(pot_fn_bin, float)
            # REMARK: We assume that the grid files exist
            if self.reference is None:
                grid_points = numpy.fromfile(points_fn_bin, float).reshape((-1,3))
                grid_weights = numpy.fromfile(weights_fn_bin, float)
            else:
                grid_points = self.reference.mol_points
                grid_weights = self.reference.mol_weights

        self.mol_points = grid_points
        self.mol_weights = grid_weights
        self.mol_densities = densities
        self.mol_potentials = potentials

    def do_atom_orbitals(self):
        if hasattr(self, "atom_orbitals"):
            return
        if self.reference is not None:
            print "Warning: The orbital analysis ignores the reference state."

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
        if self.reference is not None:
            print "Warning: The orbital analysis ignores the reference state."

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

    def do_dipoles(self):
        if hasattr(self, "hirshi_dipoles"):
            return
        if self.reference is not None:
            self.reference.do_dipoles()

        self.do_iterative_hirshfeld()
        hirshi_dipoles_fn_bin = os.path.join(self.context.workdir, "hirshi_dipoles.bin")
        molecule = self.context.fchk.molecule

        if os.path.isfile(hirshi_dipoles_fn_bin):
            print "Loading dipoles based in iterative hirshfeld partitioning."
            self.hirshi_dipoles = numpy.fromfile(hirshi_dipoles_fn_bin, float).reshape((-1,3))
        else:
            self.do_atom_grids()
            self.do_atom_densities()

            print "Computing dipoles based in iterative hirshfeld partitioning."
            pb = ProgressBar("Molecular dipoles with iterative hirshfeld weights", molecule.size)
            self.hirshi_dipoles = numpy.zeros((molecule.size,3), float)
            for i, number_i in enumerate(molecule.numbers):
                pb()
                hw = self.hirshi_weights[i]
                grid_points = self.atom_grid_points[i]
                densities = self.atom_densities[i]
                center = molecule.coordinates[i]

                for j in 0,1,2:
                    integrand = -(grid_points[:,j] - center[j])*densities*hw
                    radfun = integrate_lebedev(self.context.lebedev_weights, integrand)
                    rs = self.context.atom_table.records[number_i].rs
                    self.hirshi_dipoles[i,j] = integrate_log(rs, radfun*rs**2)
            pb()
            self.hirshi_dipoles.tofile(hirshi_dipoles_fn_bin)

        # now some nice output

        def output(filename, charges, dipoles, esp_cost, dipole_qm):
            f = file(os.path.join(self.context.outdir, filename), "w")
            print >> f, "  i        Z     Dipole-X     Dipole-Y     Dipole-Z      Dipole"
            print >> f, "------------------------------------------------------------------"
            for i, number in enumerate(molecule.numbers):
                print >> f, "% 3i  %2s  % 3i   % 10.5f   % 10.5f   % 10.5f   % 10.5f" % (
                    i+1, periodic[number].symbol, number, dipoles[i,0],
                    dipoles[i,1], dipoles[i,2], numpy.linalg.norm(dipoles[i]),
                )
            print >> f, "------------------------------------------------------------------"

            dipole_q = (molecule.coordinates*charges.reshape((-1,1))).sum(axis=0)
            dipole_p = dipoles.sum(axis=0)
            dipole_qp = dipole_q + dipole_p
            print >> f, "Molecular dipoles due to ..."
            print >> f, "charges (q)    % 10.5f   % 10.5f   % 10.5f   % 10.5f" % (
                dipole_q[0], dipole_q[1], dipole_q[2], numpy.linalg.norm(dipole_q),
            )
            print >> f, "dipoles (p)    % 10.5f   % 10.5f   % 10.5f   % 10.5f" % (
                dipole_p[0], dipole_p[1], dipole_p[2], numpy.linalg.norm(dipole_p),
            )
            print >> f, "q and p        % 10.5f   % 10.5f   % 10.5f   % 10.5f" % (
                dipole_qp[0], dipole_qp[1], dipole_qp[2], numpy.linalg.norm(dipole_qp),
            )
            print >> f, "total density  % 10.5f   % 10.5f   % 10.5f   % 10.5f" % (
                dipole_qm[0], dipole_qm[1], dipole_qm[2], numpy.linalg.norm(dipole_qm),
            )
            print >> f, "------------------------------------------------------------------"

            print >> f, "Reproduction of the external molecular ESP ..."
            print >> f, "                     RMSD             RMS"
            print >> f, "charges (q)      % 10.5e    % 10.5e" % (
                esp_cost.rmsd(charges),
                esp_cost.model_rms(charges),
            )
            print >> f, "dipoles (p)      % 10.5e    % 10.5e" % (
                esp_cost.rmsd(None, dipoles),
                esp_cost.model_rms(None, dipoles),
            )
            print >> f, "q and p          % 10.5e    % 10.5e" % (
                esp_cost.rmsd(charges, dipoles),
                esp_cost.model_rms(charges, dipoles),
            )
            print >> f, "total densitty                   % 10.5e" % esp_cost.rms
            f.close()

        dipole_qm = self.context.fchk.fields.get("Dipole Moment")
        output(
            "hirshi_dipoles.txt", self.hirshi_charges, self.hirshi_dipoles,
            self.mol_esp_cost, dipole_qm
        )
        if self.reference is not None:
            relative_dipole_qm = dipole_qm - self.reference.context.fchk.fields.get("Dipole Moment")
            output(
                "hirshi_dipoles_relative.txt",
                self.hirshi_charges - self.reference.hirshi_charges,
                self.hirshi_dipoles - self.reference.hirshi_dipoles,
                self.relative_mol_esp_cost, relative_dipole_qm,
            )

