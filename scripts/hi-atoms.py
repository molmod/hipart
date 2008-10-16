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


from hipart.core import *
from hipart.tools import guess_density, load_cube
from hipart.lebedev_laikov import get_grid, grid_fns

from molmod.data.periodic import periodic
from molmod.io.gaussian03.mkinput import mkinput
from molmod.molecules import Molecule
from molmod.units import angstrom, eV

from optparse import OptionParser
from glob import glob
import numpy, os


def get_charge_label(charge):
    if charge == 0:
        return "neut"
    elif charge < 0:
        return "neg%i" % (-charge)
    else:
        return "pos%i" % charge


def get_charge_from_label(label):
    if label == "neut":
        return 0
    elif label.startswith("neg"):
        return -int(label[3:])
    else:
        return int(label[3:])


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


def make_inputs(lot, atom_numbers, max_ion):
    pb = ProgressBar("Creating input files", len(atom_numbers)*(2*max_ion+1))
    noble_numbers = [0] + [atom.number for atom in periodic.atoms_by_number.itervalues() if atom.col == 18]
    molecule = Molecule(numpy.zeros(1, int), numpy.zeros((1,3), float), "Computer says nooooo...")

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
        molecule.numbers[0] = atom.number
        for charge in xrange(-max_ion, max_ion+1):
            charge_label = get_charge_label(charge)
            pb()
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
                dirname = os.path.join("%03i%s" % (atom.number, atom.symbol), charge_label, "mult%i" % mult)
                command = "sp scf(fermi,tight) guess=core density=current"
                if num_elec - last_noble > 1:
                    command += " polar"
                mkinput(
                    molecule, charge, mult, lot, command, "", 1,
                    "500MB", "10GB", os.path.join(dirname, "gaussian.com")
                )
    pb()


def run_jobs():
    dirnames = glob("0*/*/*/")
    dirnames.sort()
    pb = ProgressBar("Gaussian calculations", len(dirnames))
    for dirname in dirnames:
        pb()
        if not os.path.isfile(os.path.join(dirname, "gaussian.fchk")):
            os.system("(cd %s; . ~/g03.profile; g03 < gaussian.com > gaussian.log 2> /dev/null; formchk gaussian.chk gaussian.fchk > /dev/null 2> /dev/null; rm gaussian.chk)" % dirname)
    pb()


def select_ground_states(energy_field, max_ion):
    os.system("grep '%s' 0*/*/*/gaussian.fchk | grep -v gs > energies.txt" % (energy_field))

    f = file("energies.txt")

    all_energies = {
        1: {1: {1: 0.0}}
    }

    for line in f:
        words = line.split()
        energy = float(words[-1])
        tags = words[0][:words[0].find(":")].split("/")
        number = int(tags[0][:3])
        charge = get_charge_from_label(tags[1])
        mult = int(tags[2][-1])

        atom_energies = all_energies.setdefault(number, {})
        charge_energies = atom_energies.setdefault(charge, {})
        charge_energies[mult] = energy

    f.close()


    f_au = file("chieta_au.txt", "w")
    f_ev = file("chieta_ev.txt", "w")

    print >> f_au, "All values below are in atomic units."
    print >> f_au, "             A         I          chi        eta     mult(neg,neut,pos)"
    print >> f_ev, "All values below are in electron volts."
    print >> f_ev, "             A         I          chi        eta     mult(neg,neut,pos)"

    for number, atom_energies in all_energies.iteritems():
        energies = {}
        mult = {}
        symbol = periodic[number].symbol
        for charge in xrange(-max_ion, max_ion+1):
            charge_label = get_charge_label(charge)
            if charge not in atom_energies:
                continue
            energies[charge] = min(atom_energies[charge].itervalues())
            mult[charge] = min((energy, mult) for mult, energy in atom_energies[charge].iteritems())[1]
            newlink = os.path.join("%03i%s" % (number, symbol), charge_label, "gs")
            if os.path.isdir(os.path.dirname(newlink)):
                if os.path.exists(newlink): os.remove(newlink)
                os.symlink("mult%i" % mult[charge], newlink)

        chi = (energies[+1] - energies[-1])/2
        eta = (energies[+1] + energies[-1] - 2*energies[0])
        values = [energies[0] - energies[-1], energies[+1] - energies[0], chi, eta]
        print >> f_au, "% 3s" % symbol, "%2i" % number, " ".join("% 10.5f" % v for v in values), "    ", mult[-1],mult[0],mult[1]
        print >> f_ev, "% 3s" % symbol, "%2i" % number, " ".join("% 10.5f" % (v/eV) for v in values), "    ", mult[-1],mult[0],mult[1]

    f_au.close()
    f_ev.close()


def make_density_profile(density, num_lebedev, r_low, r_high, steps, atom_numbers, max_ion):
    # generate lebedev grid
    lebedev_xyz, lebedev_weights = get_grid(num_lebedev)

    # define radii
    ratio = (r_high/r_low)**(1.0/(steps-1))
    alpha = numpy.log(ratio)
    rs = r_low*numpy.exp(alpha*numpy.arange(0,steps))

    f_pro = file("densities.txt", "w")
    print >> f_pro, "Radii [bohr]", " ".join("%12.7e" % r for r in rs)
    charges = []

    # run over all directories, run cubegen, load cube data and plot
    pb = ProgressBar("Density profiles", len(atom_numbers)*(2*max_ion+1))
    for number in atom_numbers:
        atom = periodic[number]
        symbol = atom.symbol
        for charge in xrange(-max_ion, max_ion+1):
            charge_label = get_charge_label(charge)
            pb()
            dirname = os.path.join("%03i%s" % (number, symbol), charge_label, "gs")
            if not os.path.isdir(dirname): continue
            points_filename = os.path.join(dirname, "grid_points.txt")
            if not os.path.isfile(points_filename):
                # write points
                f = file(points_filename, "w")
                for r in rs:
                    for lp in lebedev_xyz:
                        x, y, z = r*lp/angstrom
                        print >> f, "% 10.5f % 10.5f % 10.5f" % (x, y, z)
                f.close()
                os.system("cd %s; . ~/g03.profile; cubegen 0 fdensity=%s gaussian.fchk grid_density.cube -5 > /dev/null 2> /dev/null < grid_points.txt" % (dirname, density))
            # load densities
            den_filename = os.path.join(dirname, "grid_density.cube")
            if os.path.isfile(den_filename):
                rhos = load_cube(den_filename, values_only=True)
                rhos = 4*numpy.pi*(rhos.reshape((-1,num_lebedev))*lebedev_weights).sum(axis=1)
                print >> f_pro, "Densities %3i %2s %+2i [a.u.]" % (number, symbol, charge), " ".join("%12.7e" % rho for rho in rhos)
                charges.append((number, symbol, charge, integrate(rs, rs**2*rhos)))
            else:
                print "Skipping %s" % den_filename

    pb()
    f_pro.close()

    counter = 0
    for number, symbol, charge, real_charge in charges:
        print "Charge check %3i %2s %+2i    %10.5e" % (number, symbol, charge, -real_charge+number-charge)
        counter += 1


parser = OptionParser("%prog [options] lot")
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
    "--rlow", default=2e-5, type='float',
    help="The smallest radius for the density profile (in angstroms)."
)
parser.add_option(
    "--rhigh", default=20.0, type='float',
    help="The largest radius for the density profile (in angstroms)."
)
parser.add_option(
    "--num-steps", default=100, type='int',
    help="The number of steps in density profile."
)
parser.add_option(
    "--max-ion", default=2, type='int',
    help="The maximum ionization to consider."
)
(options, args) = parser.parse_args()
lot, atom_str = args

if options.density is None:
    options.density = guess_density(lot)
atom_numbers = parse_numbers(atom_str)

make_inputs(lot, atom_numbers, options.max_ion)
run_jobs()
select_ground_states('%s Energy' % options.density.upper(), options.max_ion)
make_density_profile(
    options.density, options.lebedev, options.rlow*angstrom, options.rhigh*angstrom,
    options.num_steps, atom_numbers, options.max_ion
)

