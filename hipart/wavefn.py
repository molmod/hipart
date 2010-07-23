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


from molmod.io import FCHKFile
from molmod import angstrom
import os, numpy

from hipart.gint.basis import GaussianBasis
from hipart.gint.gint1_fn import gint1_fn_basis
from hipart.gint.ctools import reorder_density_matrix


__all__ = ["load_wavefunction", "FchkWaveFunction"]


def load_wavefunction(filename):
    if filename.endswith(".fchk"):
        return MyFCHKWaveFunction(filename)
    else:
        raise ValueError("File extension of %s not recognized" % filename)


class FCHKWaveFunction(object):
    def __init__(self, filename):
        self.filename = filename
        self.prefix = filename[:-5]
        fchk = FCHKFile(filename, field_labels=[
            "Charge", "Number of basis functions", "Dipole Moment",
            "Number of electrons", "Number of alpha electrons",
            "Number of beta electrons", "Alpha Orbital Energies",
            "Beta Orbital Energies",
        ])
        # for usage by the rest of the program
        self.charge = fchk.fields["Charge"]
        self.num_orbitals = fchk.fields["Number of basis functions"]
        self.dipole = fchk.fields["Dipole Moment"]
        self.num_electrons = fchk.fields["Number of electrons"]
        self.num_alpha = fchk.fields["Number of alpha electrons"]
        self.num_beta = fchk.fields["Number of beta electrons"]
        self.restricted = "Beta Orbital Energies" not in fchk.fields
        self.molecule = fchk.molecule
        # for internal usage
        self._hack_cubegen = self.restricted and (self.num_alpha != self.num_beta)
        if "mp2" in fchk.lot.lower():
            self._density_type = "mp2"
        elif "mp3" in fchk.lot.lower():
            self._density_type = "mp3"
        elif "mp4" in fchk.lot.lower():
            self._density_type = "mp4"
        else:
            self._density_type = "scf"

    def _load_cube(self, filename, N):
        f = file(filename)
        data = numpy.zeros(N, float)
        counter = 0
        for line in f:
            if counter == N:
                raise IOError("Cube file has wrong size. Expecting %i points. Got more." % (N, counter))
            line = line.strip()
            data[counter] = float(line[line.rfind(" "):])
            counter += 1
        if counter != N:
            raise IOError("Cube file has wrong size. Expecting %i points. Got %i." % (N, counter))
        if numpy.isnan(data[0]): # ugly workarond for stupid cubegen
            data[0] = data[1]
        f.close()
        return data

    def _write_cube_in(self, filename, grid_points):
        if os.path.isfile(filename):
            return
        f = file(filename, "w")
        for point in grid_points:
            print >> f, "%15.10e %15.10e %15.10e" % tuple(point/angstrom)
        f.close()

    def compute_density(self, grid):
        moldens = grid.load("moldens")
        if moldens is None:
            points_fn = "%s.txt" % grid.prefix
            den_fn = "%s_moldens.txt" % grid.prefix

            self._write_cube_in(points_fn, grid.points)

            if self._hack_cubegen:
                # ugly hack: Workaround for stupid cubegen that does not work
                # properly on ROHF calculations. pfff...
                self.compute_orbitals(grid)
                moldens = 0.0
                for j in xrange(max(self.num_alpha, self.num_beta)):
                    orb = grid.alpha_orbitals[j]
                    occup = (j < self.num_alpha) + (j < self.num_beta)
                    moldens += occup*orb**2
            else:
                if not os.path.isfile(den_fn):
                    os.system("cubegen 0 fdensity=%s %s %s -5 < %s" % (
                        self._density_type, self.filename, den_fn, points_fn
                    ))
                moldens = self._load_cube(den_fn, len(grid.points))
            grid.dump("moldens", moldens)
        grid.moldens = moldens

    def compute_spin_density(self, grid):
        molspindens = grid.load("molspindens")
        if molspindens is None:
            points_fn = "%s.txt" % grid.prefix
            sden_fn = "%s_molspindens.txt" % grid.prefix

            self._write_cube_in(points_fn, grid.points)

            if self._hack_cubegen:
                # ugly hack: Workaround for stupid cubegen that does not work
                # properly on ROHF calculations. pfff...
                self.compute_orbitals(grid)
                molspindens = 0.0
                n_min = min(self.num_alpha, self.num_beta)
                n_max = max(self.num_alpha, self.num_beta)
                for j in xrange(n_min, n_max):
                    orb = grid.alpha_orbitals[j]
                    molspindens += orb**2
            else:
                if not os.path.isfile(sden_fn):
                    os.system("cubegen 0 spin=%s %s %s -5 < %s" % (
                        self._density_type, self.filename, sden_fn, points_fn
                    ))
                molspindens = self._load_cube(sden_fn, len(grid.points))
            grid.dump("molspindens", molspindens)
        grid.molspindens = molspindens

    def compute_potential(self, grid):
        molpot = grid.load("molpot")
        if molpot is None:
            points_fn = "%s.txt" % grid.prefix
            pot_fn = "%s_pot.txt" % grid.prefix

            self._write_cube_in(points_fn, grid.points)

            if not os.path.isfile(pot_fn):
                os.system("cubegen 0 potential=%s %s %s -5 < %s" % (
                    self._density_type, self.filename, pot_fn, points_fn,
                ))
            molpot = self._load_cube(pot_fn, len(grid.points))
            grid.dump("molpot", molpot)
        grid.molpot = molpot

    def compute_orbitals(self, grid):
        alpha_orbitals = []
        beta_orbitals = []
        for i in xrange(self.num_orbitals):
            alpha_suffix = "alpha_orb%05i" % i
            alpha_orb = grid.load(alpha_suffix)
            if alpha_orb is None:
                points_fn = "%s.txt" % grid.prefix
                alpha_orb_fn = "%s_%s.txt" % (grid.prefix, alpha_suffix)

                self._write_cube_in(points_fn, grid.points)

                if not os.path.isfile(alpha_orb_fn):
                    os.system("cubegen 0 AMO=%i %s %s -5 < %s" % (
                        i+1, self.filename, alpha_orb_fn, points_fn,
                    ))
                alpha_orb = self._load_cube(alpha_orb_fn, len(grid.points))
                grid.dump(alpha_suffix, alpha_orb)
            alpha_orbitals.append(alpha_orb)
            if self.restricted:
                beta_orbitals.append(alpha_orb)
            else:
                beta_suffix = "beta_orb%05i" % i
                beta_orb = grid.load(beta_suffix)
                if beta_orb is None:
                    points_fn = "%s.txt" % grid.prefix
                    beta_orb_fn = "%s_%s.txt" % (grid.prefix, beta_suffix)

                    self._write_cube_in(points_fn, grid.points)

                    if not os.path.isfile(beta_orb_fn):
                        os.system("cubegen 0 BMO=%i %s %s -5 < %s" % (
                            i+1, self.filename, beta_orb_fn, points_fn,
                        ))
                    beta_orb = self._load_cube(beta_orb_fn, len(grid.points))
                    grid.dump(beta_suffix, beta_orb)
                beta_orbitals.append(beta_orb)

        grid.alpha_orbitals = alpha_orbitals
        grid.beta_orbitals = beta_orbitals


class MyFCHKWaveFunction(FCHKWaveFunction):
    def __init__(self, filename):
        self.filename = filename
        self.prefix = filename[:-5]
        fchk = FCHKFile(filename)
        # for use by the rest of the program
        self.charge = fchk.fields["Charge"]
        self.num_orbitals = fchk.fields["Number of basis functions"]
        self.dipole = fchk.fields["Dipole Moment"]
        self.num_electrons = fchk.fields["Number of electrons"]
        self.num_alpha = fchk.fields["Number of alpha electrons"]
        self.num_beta = fchk.fields["Number of beta electrons"]
        self.restricted = "Beta Orbital Energies" not in fchk.fields
        self.molecule = fchk.molecule
        # for internal usage
        self.basis = GaussianBasis.from_fchk(fchk)
        self.alpha_orbital_energies = fchk.fields["Alpha Orbital Energies"]
        self.beta_orbital_energies = fchk.fields.get("Beta Orbital Energies", self.alpha_orbital_energies)
        self.alpha_orbitals = fchk.fields["Alpha MO coefficients"].reshape((-1,self.basis.num_dof))
        self.alpha_orbitals = self.alpha_orbitals[:,self.basis.g03_permutation]
        self.beta_orbitals = fchk.fields.get("Beta MO coefficients")
        if self.beta_orbitals is None:
            self.beta_orbitals = self.alpha_orbitals
        else:
            self.beta_orbitals = self.beta_orbitals.reshape((-1,self.basis.num_dof))[:,self.basis.g03_permutation]
        self.density_matrices = {}
        for key in fchk.fields:
            if key.startswith("Total") and key.endswith("Density"):
                dmat = fchk.fields[key]
                reorder_density_matrix(dmat, self.basis.g03_permutation)
                self.density_matrices[key[6:-8].lower()] = dmat
        # for compatibility
        self._hack_cubegen = self.restricted and (self.num_alpha != self.num_beta)
        if "mp2" in fchk.lot.lower():
            self._density_type = "mp2"
        elif "mp3" in fchk.lot.lower():
            self._density_type = "mp3"
        elif "mp4" in fchk.lot.lower():
            self._density_type = "mp4"
        else:
            self._density_type = "scf"

    def compute_orbitals(self, grid):
        alpha_orbitals = []
        beta_orbitals = []
        for i in xrange(self.num_orbitals):
            alpha_suffix = "alpha_orb%05i" % i
            alpha_orb = grid.load(alpha_suffix)
            if alpha_orb is None:
                weights = self.alpha_orbitals[i]
                alpha_orb = self.basis.call_gint1(gint1_fn_basis, weights, grid.points)
                grid.dump(alpha_suffix, alpha_orb)
            alpha_orbitals.append(alpha_orb)
            if self.restricted:
                beta_orbitals.append(alpha_orb)
            else:
                beta_suffix = "beta_orb%05i" % i
                beta_orb = grid.load(beta_suffix)
                if beta_orb is None:
                    weights = self.beta_orbitals[i]
                    beta_orb = self.basis.call_gint1(gint1_fn_basis, weights, grid.points)
                    grid.dump(beta_suffix, beta_orb)
                    grid.dump(beta_suffix, beta_orb)
                beta_orbitals.append(beta_orb)

        grid.alpha_orbitals = alpha_orbitals
        grid.beta_orbitals = beta_orbitals
