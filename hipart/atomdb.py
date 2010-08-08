#!/usr/bin/python
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


from hipart.lebedev_laikov import get_grid, grid_fns
from hipart.grids import AtomicGrid, RLogIntGrid, ALebedevIntGrid
from hipart.wavefn import FCHKWaveFunction
from hipart.work import Work
from hipart.log import log

from molmod.periodic import periodic
from molmod.units import angstrom, electronvolt

from optparse import OptionParser
from glob import glob
import numpy, os, sys


__all__ = ["main"]


class Gaussian(object):
    def __init__(self, executable, options):
        self.executable = executable
        self.qc = options.qc
        log("Computing atomic database with program Gaussian (%s,qc=%s)" % (self.executable,self.qc))

    def make_atom_input(self, dirname, number, charge, mult, lot):
        command = "sp guess=indo density=current"
        if self.qc:
            command += " scf(qc,tight)"
        else:
            command += " scf(fermi,tight)"
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        f = file("%s/gaussian.com" % dirname, "w")
        print >> f, "%mem=500MB"
        print >> f, "%chk=gaussian.chk"
        print >> f, "%nproc=1"
        print >> f, "#p %s %s" % (lot, command)
        print >> f, ""
        print >> f, "Foo"
        print >> f, ""
        print >> f, "%i %i" % (charge, mult)
        symbol = periodic[number].symbol
        print >> f, "%2s 0.0 0.0 0.0" % symbol
        print >> f, ""
        f.close()

    def run(self, dirname):
        fn_fchk = os.path.join(dirname, "gaussian.fchk")
        if not os.path.isfile(fn_fchk):
            os.system(
                ("(cd %s; %s < gaussian.com > gaussian.log 2> /dev/null &&"
                "formchk gaussian.chk gaussian.fchk > /dev/null 2> /dev/null;"
                "rm -f gaussian.chk)") % (dirname, self.executable)
            )
            if not os.path.isfile(fn_fchk):
                return False
        return True

    def get_energies(self):
        os.system("grep 'Total Energy' 0*/*/*/gaussian.fchk | grep -v gs > energies.txt")
        f = file("energies.txt")
        # dictionary: all_energies[number][charge][mult] = energy
        all_energies = {}
        for line in f:
            words = line.split()
            energy = float(words[-1])
            tags = words[0][:words[0].find(":")].split("/")
            number = int(tags[0][:3])
            charge = label_to_charge(tags[1])
            mult = int(tags[2][-1])
            atom_energies = all_energies.setdefault(number, {})
            charge_energies = atom_energies.setdefault(charge, {})
            charge_energies[mult] = energy
        f.close()
        return all_energies

    def compute_density(self, grid, dirname):
        fchk_fn = os.path.join(dirname, "gaussian.fchk")
        wavefn = FCHKWaveFunction(fchk_fn, [])
        wavefn.compute_density(grid)


def charge_to_label(charge):
    if charge == 0:
        return "neut"
    elif charge < 0:
        return "neg%i" % (-charge)
    else:
        return "pos%i" % charge


def label_to_charge(label):
    if label == "neut":
        return 0
    elif label.startswith("neg"):
        return -int(label[3:])
    else:
        return int(label[3:])


def iter_states(atom_numbers, max_ion):
    noble_numbers = [0] + [atom.number for atom in periodic.atoms_by_number.itervalues() if atom.col == 18]
    for i in atom_numbers:
        atom = periodic[i]
        next_noble = min(n for n in noble_numbers if n >= atom.number)
        last_noble = max(n for n in noble_numbers if n < atom.number)

        if atom.row == 1:
            num_states = 1
        elif atom.row == 2 or atom.row == 3:
            num_states = 2
        elif atom.row == 4 or atom.row == 5:
            num_states = 3
        else:
            num_states = 4
        for charge in xrange(-max_ion, max_ion+1):
            num_elec = atom.number - charge
            if num_elec <= 0: continue
            charge_num_states = max(1, min(
                num_states,
                (num_elec - last_noble)/2+1,
                (next_noble - num_elec)/2+1
            ))
            if num_elec % 2 == 0:
                mults = list(2*i+1 for i in xrange(charge_num_states))
            else:
                mults = list(2*i+2 for i in xrange(charge_num_states))
            for mult in mults:
                yield atom, charge, mult


def make_inputs(program, lot, atom_numbers, max_ion):
    states = list(iter_states(atom_numbers, max_ion))
    pb = log.pb("Creating input files", len(states))
    dirnames = []

    for atom, charge, mult in states:
        pb()
        charge_label = charge_to_label(charge)
        dirname = os.path.join("%03i%s" % (atom.number, atom.symbol), charge_label, "mult%i" % mult)
        program.make_atom_input(dirname, atom.number, charge, mult, lot)
        dirnames.append(dirname)
    pb()
    return dirnames


def run_jobs(program, dirnames):
    pb = log.pb("Atomic computations", len(dirnames))
    failed = []
    for dirname in dirnames:
        pb()
        succes = program.run(dirname)
        if not succes:
            failed.append(dirname)
    pb()
    if len(failed) == len(dirnames):
        log("Could not execute any job. Is %s in the PATH?" % program.executable)
        sys.exit(-1)
    if len(failed) > 0:
        log("Some jobs failed:")
        for dirname in failed:
            log("  %s" % dirname)


def select_ground_states(program, max_ion):
    log("Selecting ground states.")
    all_energies = program.get_energies()

    if 1 in all_energies:
        # Computation on a proton without electrons: solve manually.
        # crunch crunch ...
        all_energies[1][1] = {1: 0.0}

    f_au = file("chieta_au.txt", "w")
    f_ev = file("chieta_ev.txt", "w")

    print >> f_au, "All values below are in atomic units."
    print >> f_au, "             A         I          chi        eta     mult(neg,neut,pos)"
    print >> f_ev, "All values below are in electron volts."
    print >> f_ev, "             A         I          chi        eta     mult(neg,neut,pos)"

    for number, atom_energies in sorted(all_energies.iteritems()):
        energies = {}
        mult = {}
        symbol = periodic[number].symbol
        for charge in xrange(-max_ion, max_ion+1):
            charge_label = charge_to_label(charge)
            if charge not in atom_energies:
                continue
            energies[charge] = min(atom_energies[charge].itervalues())
            mult[charge] = min((energy, mult) for mult, energy in atom_energies[charge].iteritems())[1]
            newlink = os.path.join("%03i%s" % (number, symbol), charge_label, "gs")
            if os.path.isdir(os.path.dirname(newlink)):
                if os.path.exists(newlink): os.remove(newlink)
                os.symlink("mult%i" % mult[charge], newlink)

        if -1 in energies and 0 in energies and 1 in energies:
            chi = (energies[+1] - energies[-1])/2
            eta = (energies[+1] + energies[-1] - 2*energies[0])
            values = [energies[0] - energies[-1], energies[+1] - energies[0], chi, eta]
            print >> f_au, "% 3s" % symbol, "%2i" % number, " ".join("% 10.5f" % v for v in values), "    ", mult[-1],mult[0],mult[1]
            print >> f_ev, "% 3s" % symbol, "%2i" % number, " ".join("% 10.5f" % (v/electronvolt) for v in values), "    ", mult[-1],mult[0],mult[1]

    f_au.close()
    f_ev.close()


def make_density_profiles(program, num_lebedev, r_low, r_high, steps, atom_numbers, max_ion, do_work):
    # generate lebedev grid
    lebedev_xyz, lebedev_weights = get_grid(num_lebedev)

    # define radii
    rgrid = RLogIntGrid(r_low, r_high, steps)
    agrid = ALebedevIntGrid(num_lebedev)

    f_pro = file("densities.txt", "w")
    print >> f_pro, rgrid.get_description()
    charges = []

    # run over all directories, run cubegen, load cube data
    pb = log.pb("Density profiles", len(atom_numbers)*(2*max_ion+1))
    for number in atom_numbers:
        symbol = periodic[number].symbol
        for charge in xrange(-max_ion, max_ion+1):
            charge_label = charge_to_label(charge)
            pb()
            dirname = os.path.join("%03i%s" % (number, symbol), charge_label, "gs")
            # get the grid
            if not os.path.isdir(dirname): continue
            if do_work:
                work = Work(dirname)
            else:
                work = Work()
            grid = AtomicGrid.from_prefix("grid", work)
            if grid is None:
                center = numpy.zeros(3,float)
                grid = AtomicGrid.from_parameters("grid", work, center, rgrid, agrid)
            # compute densities
            program.compute_density(grid, dirname)
            # this is spherical averaging, i.e. integral/(4*pi)
            radrhos = agrid.integrate(grid.moldens)/(4*numpy.pi)
            print >> f_pro, "%3i %+2i" % (number, charge),
            # leave out near zeros to save space and time
            print >> f_pro, " ".join("%22.16e" % rho for rho in radrhos if rho > 1e-100)
            check = rgrid.integrate(4*numpy.pi*rgrid.rs*rgrid.rs*radrhos)
            charges.append((number, symbol, charge, check))


    pb()
    f_pro.close()

    counter = 0
    for number, symbol, charge, real_charge in charges:
        log("Total charge error: %3i %2s %+2i    % 10.5e" % (number, symbol, charge, -real_charge+number-charge))
        counter += 1


def parse_numbers(atom_str):
    numbers = []
    for item in atom_str.split(','):
        if '-' in item:
            words = item.split("-")
            if len(words) != 2:
                raise ValueError("Each item should contain at most one dash.")
            begin = int(words[0])
            end = int(words[1])
            numbers.extend(xrange(begin,end+1))
        else:
            numbers.append(int(item))
    return numbers


usage = """\
%prog [options] executable lot atoms

%prog computes a database of pro-atomic densities.

The following arguments are mandatory:
  * executable  --  the name of the Gaussian binary (g03 or g09)
  * lot  --  The level of theory to be used in Gaussian input notation.
  * atoms  -- The atoms to be computed. One can specify ranges, e.g 1,2-5'
              (avoid whitespace)

It is recommended to run this script in a directory that is initially empty and
that will contain nothing but the generated atom database. This script will
generate quite a few files and subdirectories.

Examples:

%prog g03 MP2/Aug-CC-pVDZ 1-10,17
%prog g09 HF/3-21G 1,6,7,8 -l 110
"""

def parse_args(args):
    parser = OptionParser(usage)
    parser.add_option(
        "-l", "--lebedev", default=350, type='int',
        help="The number of grid points for the spherical averaging. "
        "[default=%default]. Select from: " + (", ".join(str(i) for i in sorted(grid_fns)))
    )
    parser.add_option(
        "--rlow", default=2e-5, type='float',
        help="The smallest radius for the density profile (in angstroms). [default=%default]"
    )
    parser.add_option(
        "--rhigh", default=20.0, type='float',
        help="The largest radius for the density profile (in angstroms). [default=%default]"
    )
    parser.add_option(
        "--num-steps", default=100, type='int',
        help="The number of steps in density profile. [default=%default]"
    )
    parser.add_option(
        "--max-ion", default=2, type='int',
        help="The maximum ionization to consider. [default=%default]"
    )
    parser.add_option(
        "--qc", default=False, action="store_true",
        help="Specify the qc convergence scheme in Gaussian input. [default=%default]"
    )
    parser.add_option(
        "--no-work", default=True, action='store_false', dest='do_work',
        help="Do not save intermediate results in work directory for later reuse."
    )
    parser.add_option(
        "-q", "--quiet", default=True, action='store_false', dest='verbose',
        help="Do not write any screen output."
    )
    (options, args) = parser.parse_args(args)
    if len(args) != 3:
        parser.error("Expecting three arguments: executable, level of theory (+ basis set) and the atom specification.")
    executable, lot, atom_str = args
    atom_numbers = parse_numbers(atom_str)
    return options, executable, lot, atom_numbers


def main(args=None):
    options, executable, lot, atom_numbers = parse_args(args=None)
    log.verbose = options.verbose
    program = Gaussian(executable, options)
    dirnames = make_inputs(program, lot, atom_numbers, options.max_ion)
    run_jobs(program, dirnames)
    select_ground_states(program, options.max_ion)
    make_density_profiles(
        program,
        options.lebedev, options.rlow*angstrom, options.rhigh*angstrom,
        options.num_steps, atom_numbers, options.max_ion, options.do_work
    )
