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


__all__ = ["load_wavefunction", "FchkWaveFunction"]


def load_wavefunction(filename):
    if filename.endswith(".fchk"):
        return FCHKWaveFunction(filename)
    else:
        raise ValueError("File extension of %s not recognized" % filename)


class FCHKWaveFunction(object):
    def __init__(self, filename):
        self.filename = filename
        self.prefix = filename[:-5]
        fchk = FCHKFile(filename, field_labels=[
            "Charge", "Number of basis functions", "Dipole Moment",
            "Number of electrons", "Number of alpha electrons",
            "Number of beta electrons"
        ])
        # for usage by the rest of the program
        self.charge = fchk.fields["Charge"]
        self.num_orbitals = fchk.fields["Number of basis functions"]
        self.dipole = fchk.fields["Dipole Moment"]
        self.num_electrons = fchk.fields["Number of electrons"]
        self.num_alpha = fchk.fields["Number of alpha electrons"]
        self.num_beta = fchk.fields["Number of beta electrons"]
        self.molecule = fchk.molecule
        # for internal usage
        self._hack_cubegen = fchk.lot.startswith("ROHF")
        if not self._hack_cubegen and fchk.lot.startswith("RO"):
            raise RuntimeError("Can not deal with restricted open calculations except for ROHF.")
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
                    orb = grid.orbitals[i]
                    occup = (j < self.num_alpha) + (j < self.num_beta)
                    moldens += occup*orb**2
            else:
                if not os.path.isfile(den_fn):
                    os.system(". ~/g03.profile; cubegen 0 fdensity=%s %s %s -5 < %s" % (
                        self._density_type, self.filename, den_fn, points_fn
                    ))
                moldens = self._load_cube(den_fn, len(grid.points))
            grid.dump("moldens", moldens)
        grid.moldens = moldens

    def compute_potential(self, grid):
        molpot = grid.load("molpot")
        if molpot is None:
            points_fn = "%s.txt" % grid.prefix
            pot_fn = "%s_pot.txt" % grid.prefix

            self._write_cube_in(points_fn, grid.points)

            if not os.path.isfile(pot_fn):
                os.system(". ~/g03.profile; cubegen 0 potential=%s %s %s -5 < %s" % (
                    self._density_type, self.filename, pot_fn, points_fn,
                ))
            molpot = self._load_cube(pot_fn, len(grid.points))
            grid.dump("molpot", molpot)
        grid.molpot = molpot

    def compute_orbitals(self, grid):
        orbitals = []
        for i in xrange(self.num_orbitals):
            suffix = "orb%05i" % i
            orb = grid.load(suffix)
            if orb is None:
                points_fn = "%s.txt" % grid.prefix
                orb_fn = "%s_%s.txt" % (grid.prefix, suffix)

                self._write_cube_in(points_fn, grid.points)

                if not os.path.isfile(orb_fn):
                    os.system(". ~/g03.profile; cubegen 0 MO=%i %s %s -5 < %s" % (
                        orb_index+1, self.filename, orb_fn, points_fn,
                    ))
                orb = self._load_cube(orb_fn, len(grid.points))
                grid.dump(suffix, orb)
            orbitals.append(orb)
        grid.orbitals = orbitals
