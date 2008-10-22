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


from hipart.log import log
from hipart.tools import write_cube_in, cubegen_density, cubegen_potential, \
    cubegen_orbital, get_atom_grid, compute_hirshfeld_weights
from hipart.integrate import cumul_integrate_log, integrate_log, integrate_lebedev
from hipart.fit import ESPCostFunction
from hipart.lebedev_laikov import get_grid
from hipart.atoms import AtomTable, AtomFn

from molmod.transformations import random_rotation
from molmod.data.periodic import periodic

import os, numpy


__all__ = ["ComputeError", "BaseCache", "HirshfeldICache", "cache_classes"]


noble_numbers = numpy.array([0,2,10,18,36,54,86,118])
core_sizes = dict((number, noble_numbers[noble_numbers<=number].max()) for number in periodic.yield_numbers())


class ComputeError(Exception):
    pass


class BaseCache(object):
    num_args = None

    @classmethod
    def new_from_args(cls, context, args):
        raise NotImplementedError

    def __init__(self, context, prefix):
        self.prefix = prefix
        self.context = context
        if self.context.reference is not None:
            self.reference = self.clone(self.context.reference)
        else:
            self.reference = None

    def get_rs(self, number):
        raise NotImplementedError

    def clone(self, other_context):
        raise NotImplementedError

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
        log.begin("Atomic grids")
        molecule = self.context.fchk.molecule
        if self.reference is not None:
            log("Using atomic grids from reference context.")
            self.reference.do_atom_grids()
            self.atom_grid_distances = self.reference.atom_grid_distances
            self.atom_grid_points = self.reference.atom_grid_points
            for i in xrange(molecule.size):
                self.work_link("atom%05igrid.txt" % i)
        else:
            workdir = self.context.workdir

            self.atom_grid_distances = {}
            self.atom_grid_points = []
            pb = log.pb("Computing/Loading atomic grids (and distances)", molecule.size**2)
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
                        self.get_rs(number_i),
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
        log.end("Atom grids")

    def do_atom_densities(self):
        if hasattr(self, "atom_densities"):
            return
        log.begin("Molecular densities on atomic grids")
        self.do_atom_grids()
        molecule = self.context.fchk.molecule
        workdir = self.context.workdir

        pb = log.pb("Computing/Loading densities", molecule.size)
        self.atom_densities = []
        for i, number_i in enumerate(molecule.numbers):
            pb()
            grid_size = self.context.num_lebedev*len(self.get_rs(number_i))
            den_fn = os.path.join(workdir, "atom%05idens.cube" % i)
            den_fn_bin = "%s.bin" % den_fn
            if os.path.isfile(den_fn_bin):
                densities = numpy.fromfile(den_fn_bin, float)
            else:
                if os.path.isfile(den_fn):
                    densities = load_cube(den_fn, cube_size)
                else:
                    grid_fn = os.path.join(workdir, "atom%05igrid.txt" % i)
                    densities = cubegen_density(
                        grid_fn, den_fn, self.context.fchk,
                        self.context.options.density, grid_size
                    )
                densities.tofile(den_fn_bin)
            self.atom_densities.append(densities)
        pb()
        log.end("Molecular densities on atomic grids")

    def do_cusp_radii(self):
        if hasattr(self, "cusp_radii"):
            return
        if self.reference is not None:
            self.reference.do_cusp_radii()
            self.cusp_radii = self.reference.cusp_radii
            return

        molecule = self.context.fchk.molecule
        workdir = self.context.workdir

        cusp_radii_fn_bin = os.path.join(workdir, "cusp_radii.bin")
        if os.path.isfile(cusp_radii_fn_bin):
            log("Loading cusp radii")
            self.cusp_radii = numpy.fromfile(cusp_radii_fn_bin, float)
        else:
            self.do_atom_densities()
            log("Computing cusp radii")
            self.cusp_radii = numpy.zeros(molecule.size, float)
            for i, number_i in enumerate(molecule.numbers):
                if number_i < 3:
                    self.cusp_radii[i] = 0.2
                else:
                    densities = self.atom_densities[i]
                    radfun = integrate_lebedev(self.context.lebedev_weights, densities)
                    rs = self.get_rs(number_i)
                    charge_int = cumul_integrate_log(rs, radfun*rs**2)
                    j = charge_int.searchsorted([core_sizes[number_i]])[0]
                    self.cusp_radii[i] = rs[j]
            self.cusp_radii.tofile(cusp_radii_fn_bin)

    def do_esp_costfunction(self):
        if hasattr(self, "mol_esp_cost"):
            return
        log.begin("ESP Cost function")
        if self.reference is not None:
            log("Also construct reference ESP")
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
        log("Written %s" % outfn)
        if self.reference is not None:
            relative_total_charge = total_charge - self.reference.context.fchk.fields.get("Charge")
            self.relative_mol_esp_cost = ESPCostFunction(
                coordinates, self.mol_points, self.mol_weights,
                self.reference.mol_densities,
                self.mol_potentials - self.reference.mol_potentials,
                relative_total_charge
            )
            outfn = os.path.join(self.context.outdir, "mol_esp_cost_relative.txt")
            self.relative_mol_esp_cost.write_to_file(outfn)
            log("Written %s" % outfn)
        log.end("ESP Cost function")

    def _do_molecular_grid(self):
        if hasattr(self, "mol_points") and hasattr(self, "mol_weights") and \
           hasattr(self, "mol_densities") and hasattr(self, "mol_potentials"):
            return

        log.begin("Molecular grid + density and potential on this grid")
        lebedev_xyz, lebedev_weights = get_grid(self.context.options.mol_lebedev)
        molecule = self.context.fchk.molecule
        workdir = self.context.workdir

        den_fn = os.path.join(workdir, "molecule_dens.cube")
        pot_fn = os.path.join(workdir, "molecule_pot.cube")
        den_fn_bin = "%s.bin" % den_fn
        pot_fn_bin = "%s.bin" % pot_fn
        points_fn = os.path.join(workdir, "molecule_points.txt")
        points_fn_bin = os.path.join(workdir, "molecule_points.bin")
        weights_fn_bin = os.path.join(workdir, "molecule_weights.bin")

        self.do_cusp_radii()

        if not (os.path.isfile(den_fn_bin) and os.path.isfile(pot_fn_bin)):
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

                    grid_points = []
                    grid_weights = []
                    pb = log.pb("Constructing molecular grid", scale_steps)
                    for scale in scales:
                        pb()
                        radii = scale*self.cusp_radii
                        for i in xrange(molecule.size):
                            rot = random_rotation()
                            for j in xrange(len(lebedev_xyz)):
                                my_point = radii[i]*numpy.dot(rot.r, lebedev_xyz[j]) + molecule.coordinates[i]
                                distances = numpy.sqrt(((molecule.coordinates - my_point)**2).sum(axis=1))
                                if (distances < scales[0]*self.cusp_radii).any():
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
                log("Linking grid from reference state")
                grid_points = self.reference.mol_points
                grid_weights = self.reference.mol_weights
                self.work_link("molecule_points.bin")
                self.work_link("molecule_weights.bin")

            write_cube_in(points_fn, grid_points) # prepare for cubegen
            grid_size = len(grid_weights)
            log("Molecular density on moleculer grid.")
            densities = cubegen_density(
                points_fn, den_fn, self.context.fchk,
                self.context.options.density, grid_size
            )
            densities.tofile(den_fn_bin)
            log("Molecular potential on moleculer grid.")
            potentials = cubegen_potential(
                points_fn, pot_fn, self.context.fchk,
                self.context.options.density, grid_size
            )
            potentials.tofile(pot_fn_bin)
        else:
            log("Loading molecular data from workdir")
            densities = numpy.fromfile(den_fn_bin, float)
            potentials = numpy.fromfile(pot_fn_bin, float)
            if self.reference is None:
                # REMARK: We assume that the grid files exist
                grid_points = numpy.fromfile(points_fn_bin, float).reshape((-1,3))
                grid_weights = numpy.fromfile(weights_fn_bin, float)
            else:
                grid_points = self.reference.mol_points
                grid_weights = self.reference.mol_weights

        self.mol_points = grid_points
        self.mol_weights = grid_weights
        self.mol_densities = densities
        self.mol_potentials = potentials
        log.end("Molecular grid + density and potential on this grid")

    def do_partitions(self):
        if hasattr(self, "partition_weights") and hasattr(self, "pro_atom_fns"):
            return

        log.begin("Pro-atoms")
        weights_tpl_bin = os.path.join(self.context.workdir, "%s_weights%%05i.bin" % self.prefix)
        pro_atoms_tpl_bin = os.path.join(self.context.workdir, "%s_proatom%%05i.bin" % self.prefix)

        # Try to read the data from the workdir
        self._load_partitions(weights_tpl_bin, pro_atoms_tpl_bin)
        if self.partition_weights is None or self.pro_atom_fns is None:
            log("Could not load partitions from workdir. Computing them...")
            self._compute_partitions()

            log("Writing pro atoms to workdir")
            for i, pafn in enumerate(self.pro_atom_fns):
                pafn.density.y.tofile(pro_atoms_tpl_bin % i)
            log("Writing partition weights to workdir")
            for i, pw in enumerate(self.partition_weights):
                pw.tofile(weights_tpl_bin % i)
        log.end("Pro-atoms")

    def _load_partitions(self, weights_tpl_bin, pro_atoms_tpl_bin):
        log("Try to load partition weights")
        molecule = self.context.fchk.molecule

        self.partition_weights = []
        for i in xrange(molecule.size):
            if os.path.isfile(weights_tpl_bin % i):
                pw = numpy.fromfile(weights_tpl_bin % i, float)
                self.partition_weights.append(pw)
            else:
                self.partition_weights = None
                return
        self.pro_atom_fns = []
        for i, number_i in enumerate(molecule.numbers):
            if os.path.isfile(pro_atoms_tpl_bin % i):
                rs = self.get_rs(number_i)
                rhos = numpy.fromfile(pro_atoms_tpl_bin % i, float)
                self.pro_atom_fns.append(AtomFn(rs, rhos))
            else:
                self.pro_atom_fns = None
                return

    def _compute_partitions(self, weights_tpl_bin, pro_atoms_tpl_bin):
        raise NotImplementedError

    def do_charges(self):
        if hasattr(self, "charges"):
            return
        log.begin("Charges")

        if self.reference is not None:
            self.reference.do_charges()

        self.do_partitions()
        self.do_esp_costfunction()
        charges_fn_bin = os.path.join(self.context.workdir, "%s_charges.bin" % self.prefix)
        molecule = self.context.fchk.molecule

        if os.path.isfile(charges_fn_bin):
            log("Loading charges.")
            self.charges = numpy.fromfile(charges_fn_bin, float)
        else:
            self.do_atom_grids()
            self.do_atom_densities()

            pb = log.pb("Computing charges", molecule.size)
            self.charges = numpy.zeros(molecule.size, float)
            for i, number_i in enumerate(molecule.numbers):
                pb()
                pw = self.partition_weights[i]
                grid_points = self.atom_grid_points[i]
                densities = self.atom_densities[i]
                center = molecule.coordinates[i]

                integrand = densities*pw
                radfun = integrate_lebedev(self.context.lebedev_weights, integrand)
                rs = self.get_rs(number_i)
                self.charges[i] = number_i - integrate_log(rs, radfun*rs**2)
            pb()
            if self.context.options.fix_total_charge:
                log("Ugly step: Adding constant to charges so that the total charge is zero.")
                self.charges -= (self.charges.sum() - self.context.fchk.fields.get("Charge"))/molecule.size
            self.charges.tofile(charges_fn_bin)

        # now some nice output
        def output(filename, charges, esp_cost):
            filename = os.path.join(self.context.outdir, filename)
            f = file(filename, "w")
            print >> f, "number of atoms:", molecule.size
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
            log("Written %s" % filename)

        output("%s_charges.txt" % self.prefix, self.charges, self.mol_esp_cost)

        if self.reference is not None:
            relative_charges = self.charges - self.reference.charges
            output(
                "%s_charges_relative.txt" % self.prefix, relative_charges,
                self.relative_mol_esp_cost
            )
        log.end("Charges")

    def do_dipoles(self):
        if hasattr(self, "dipoles"):
            return
        if self.reference is not None:
            self.reference.do_dipoles()

        log.begin("Dipoles")
        self.do_partitions()
        self.do_esp_costfunction()
        self.do_charges()
        dipoles_fn_bin = os.path.join(self.context.workdir, "%s_dipoles.bin" % self.prefix)
        molecule = self.context.fchk.molecule

        if os.path.isfile(dipoles_fn_bin):
            log("Loading dipoles.")
            self.dipoles = numpy.fromfile(dipoles_fn_bin, float).reshape((-1,3))
        else:
            self.do_atom_grids()
            self.do_atom_densities()

            pb = log.pb("Computing dipoles", molecule.size)
            self.dipoles = numpy.zeros((molecule.size,3), float)
            for i, number_i in enumerate(molecule.numbers):
                pb()
                pw = self.partition_weights[i]
                grid_points = self.atom_grid_points[i]
                densities = self.atom_densities[i]
                center = molecule.coordinates[i]

                for j in 0,1,2:
                    integrand = -(grid_points[:,j] - center[j])*densities*pw
                    radfun = integrate_lebedev(self.context.lebedev_weights, integrand)
                    rs = self.get_rs(number_i)
                    self.dipoles[i,j] = integrate_log(rs, radfun*rs**2)
            pb()
            self.dipoles.tofile(dipoles_fn_bin)

        # now some nice output
        def output(filename, charges, dipoles, esp_cost, dipole_qm):
            filename = os.path.join(self.context.outdir, filename)
            f = file(filename, "w")
            print >> f, "number of atoms:", molecule.size
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
            print >> f, "Molecular dipole due to ..."
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
            print >> f, "total density                    % 10.5e" % esp_cost.rms
            f.close()
            log("Written %s" % filename)

        dipole_qm = self.context.fchk.fields.get("Dipole Moment")
        output(
            "%s_dipoles.txt" % self.prefix, self.charges, self.
            dipoles, self.mol_esp_cost, dipole_qm
        )
        if self.reference is not None:
            relative_dipole_qm = dipole_qm - self.reference.context.fchk.fields.get("Dipole Moment")
            output(
                "%s_dipoles_relative.txt" % self.prefix,
                self.charges - self.reference.charges,
                self.dipoles - self.reference.dipoles,
                self.relative_mol_esp_cost, relative_dipole_qm,
            )
        log.end("Dipoles")

    def do_atom_orbitals(self):
        if hasattr(self, "atom_orbitals"):
            return

        log.begin("Orbitals on atomic grids")
        if self.reference is not None:
            log("Warning: The orbital analysis ignores the reference state.")

        self.do_atom_grids()
        molecule = self.context.fchk.molecule
        num_orbitals = self.context.fchk.fields.get("Number of basis functions")
        workdir = self.context.workdir

        pb = log.pb("Computing/Loading orbitals", molecule.size*num_orbitals)
        self.atom_orbitals = []
        for i, number_i in enumerate(molecule.numbers):
            orbitals = []
            self.atom_orbitals.append(orbitals)
            grid_size = self.context.num_lebedev*len(self.get_rs(number_i))
            grid_fn = os.path.join(workdir, "atom%05igrid.txt" % i)
            for j in xrange(num_orbitals):
                pb()
                cube_fn = os.path.join(workdir, "atom%05iorb%05i.cube" % (i,j))
                cube_fn_bin = "%s.bin" % cube_fn
                if os.path.isfile(cube_fn_bin):
                    wavefn = numpy.fromfile(cube_fn_bin, float)
                else:
                    wavefn = cubegen_orbital(grid_fn, cube_fn, self.context.fchk, j, grid_size)
                    wavefn.tofile(cube_fn_bin)
                orbitals.append(wavefn)
        pb()
        log.end("Orbitals on atomic grids")

    def do_atom_matrices(self):
        if hasattr(self, "atom_matrices"):
            return
        log.begin("Partinioning density matrix")
        if self.reference is not None:
            log("Warning: The orbital analysis ignores the reference state.")

        self.do_atom_orbitals()
        self.do_partitions()
        num_orbitals = self.context.fchk.fields.get("Number of basis functions")
        molecule = self.context.fchk.molecule

        pb = log.pb("Computing matrices", molecule.size)
        self.atom_matrices = []
        for i, number_i in enumerate(molecule.numbers):
            pb()
            orbitals = self.atom_orbitals[i]
            pw = self.partition_weights[i]
            matrix = numpy.zeros((num_orbitals,num_orbitals), float)
            self.atom_matrices.append(matrix)
            for j1 in xrange(num_orbitals):
                for j2 in xrange(j1+1):
                    integrand = orbitals[j1]*orbitals[j2]*pw
                    radfun = integrate_lebedev(self.context.lebedev_weights, integrand)
                    rs = self.get_rs(number_i)
                    value = integrate_log(rs, radfun*rs**2)
                    matrix[j1,j2] = value
                    matrix[j2,j1] = value
        pb()

        filename = os.path.join(self.context.outdir, "%s_weight_elements.txt" % self.prefix)
        f = file(filename, "w")
        print >> f, "number of orbitas:", num_orbitals
        print >> f, "number of atoms: ", molecule.size
        for i, number_i in enumerate(molecule.numbers):
            print >> f, "Atom %i: %s" % (i, periodic[number_i].symbol)
            for row in self.atom_matrices[i]:
                print >> f, " ".join("% 15.10e" % value for value in row)
        f.close()
        log("Written %s" % filename)
        log.end("Partinioning density matrix")


class TableBaseCache(BaseCache):
    num_args = 1

    @classmethod
    def new_from_args(cls, context, args):
        atom_table = AtomTable(args[0])
        return cls(context, atom_table)

    def __init__(self, context, atom_table, prefix):
        BaseCache.__init__(self, context, prefix)
        self.atom_table = atom_table
        # write the rs to the workdir for plotting purposes:
        atom_table.rs.tofile(os.path.join(self.context.workdir, "%s_rs.bin" % self.prefix))

    def get_rs(self, number):
        return self.atom_table.records[number].rs

    def clone(self, other_context):
        return self.__class__(other_context, self.atom_table)


class HirshfeldCache(TableBaseCache):
    def __init__(self, context, atom_table):
        TableBaseCache.__init__(self, context, atom_table, "hirsh")

    def _compute_partitions(self):
        log.begin("Normal Hirshfeld (with neutral pro-atoms)")
        self.do_atom_grids()

        molecule = self.context.fchk.molecule
        self.pro_atom_fns = []
        for i, number_i in enumerate(molecule.numbers):
            self.pro_atom_fns.append(self.atom_table.records[number_i].get_atom_fn(0.0))

        self.partition_weights = []
        for i, number_i in enumerate(molecule.numbers):
            hw = compute_hirshfeld_weights(
                i, self.pro_atom_fns, self.context.num_lebedev,
                self.atom_grid_distances
            )
            self.partition_weights.append(hw)

        log.end("Normal Hirshfeld (with neutral pro-atoms)")


class HirshfeldICache(TableBaseCache):
    def __init__(self, context, atom_table):
        TableBaseCache.__init__(self, context, atom_table, "hirshi")

    def _compute_partitions(self):
        log.begin("Iterative Hirshfeld")
        molecule = self.context.fchk.molecule
        self.do_atom_densities()

        counter = 0
        old_charges = numpy.zeros(molecule.size, float)
        while True:
            # construct the pro-atom density functions, using the densities
            # from the previous iteration.
            atom_fns = []
            for i, number_i in enumerate(molecule.numbers):
                atom_fns.append(self.atom_table.records[number_i].get_atom_fn(old_charges[i]))

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
                rs = self.atom_table.records[number_i].rs
                num_electrons = integrate_log(rs, radfun*rs**2)

                charges.append(number_i - num_electrons)
                hirshfeld_weights.append(hw)

            # ordinary blablabla ...
            charges = numpy.array(charges)
            max_change = abs(charges-old_charges).max()
            log("Iteration %03i    max change = %10.5e    total charge = %10.5e" % (
                counter, max_change, charges.sum()
            ))
            if max_change < self.context.options.threshold:
                break
            counter += 1
            old_charges = charges

        self.partition_weights = hirshfeld_weights
        self.pro_atom_fns = atom_fns
        log.end("Iterative Hirshfeld")


cache_classes = {
    "hirsh": HirshfeldCache,
    "hirshi": HirshfeldICache,
}
