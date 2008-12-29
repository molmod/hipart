#!/usr/bin/python
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
from hipart.integrate import integrate_log
from hipart.tools import guess_density_type, write_cube_in, cubegen_density, get_atom_grid
from hipart.lebedev_laikov import get_grid, grid_fns

from molmod.data.periodic import periodic
from molmod.io.gaussian03.mkinput import mkinput
from molmod.io.gaussian03.fchk import FCHKFile
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


def make_inputs(lot, atom_numbers, max_ion, use_qc):
    pb = log.pb("Creating input files", len(atom_numbers)*(2*max_ion+1))
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
        molecule = Molecule(numpy.array([atom.number], int), numpy.zeros((1,3), float), "Computer says nooooo...")
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
                command = "sp guess=indo density=current"
                if use_qc:
                    command += " scf(qc,tight)"
                else:
                    command += " scf(fermi,tight)"
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
    pb = log.pb("Gaussian calculations", len(dirnames))
    failed = []
    for dirname in dirnames:
        pb()
        if not os.path.isfile(os.path.join(dirname, "gaussian.fchk")):
            os.system("(cd %s; . ~/g03.profile; g03 < gaussian.com > gaussian.log 2> /dev/null; formchk gaussian.chk gaussian.fchk > /dev/null 2> /dev/null; rm gaussian.chk)" % dirname)
            if not os.path.isfile(os.path.join(dirname, "gaussian.fchk")):
                failed.append(dirname)
    pb()
    if len(failed) > 0:
        print "Some jobs failed:"
        for dirname in failed:
            print "  %s" % dirname


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

    for number, atom_energies in sorted(all_energies.iteritems()):
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


def make_density_profile(density_type, num_lebedev, r_low, r_high, steps, atom_numbers, max_ion):
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
    pb = log.pb("Density profiles", len(atom_numbers)*(2*max_ion+1))
    for number in atom_numbers:
        atom = periodic[number]
        symbol = atom.symbol
        for charge in xrange(-max_ion, max_ion+1):
            charge_label = get_charge_label(charge)
            pb()
            workdir = os.path.join("%03i%s" % (number, symbol), charge_label, "gs")
            if not os.path.isdir(workdir): continue
            den_fn = os.path.join(workdir, "densities.cube")
            den_bin = os.path.join(workdir, "densities.cube.bin")
            if os.path.isfile(den_bin):
                radrhos = numpy.fromfile(den_bin, float)
            else:
                grid_fn = os.path.join(workdir, "grid.txt")
                atom_grid = get_atom_grid(lebedev_xyz, numpy.zeros(3,float), rs)
                write_cube_in(grid_fn, atom_grid)
                fchk_fn = os.path.join(workdir, "gaussian.fchk")
                fchk = FCHKFile(fchk_fn, field_labels=["Number of electrons"])
                rhos = cubegen_density(grid_fn, den_fn, fchk, options.density, num_lebedev*len(rs))
                radrhos = (rhos.reshape((-1,num_lebedev))*lebedev_weights).sum(axis=1) # this is averaging, i.e. integral/(4*pi)
                radrhos.tofile(den_bin)
                #else:
                #    print "Skipping %s" % den_fn
                #    continue
            print >> f_pro, "Densities %3i %2s %+2i [a.u.]" % (number, symbol, charge),
            print >> f_pro, " ".join("%12.7e" % rho for rho in radrhos)
            charges.append((number, symbol, charge, integrate_log(rs, 4*numpy.pi*rs*rs*radrhos)))


    pb()
    f_pro.close()

    counter = 0
    for number, symbol, charge, real_charge in charges:
        log("Total charge error: %3i %2s %+2i    %10.5e" % (number, symbol, charge, -real_charge+number-charge))
        counter += 1


parser = OptionParser("""%prog [options] lot atoms

%prog computs a database of pro-atom densities.

The following arguments are mandatory:
  * lot  --  The level of theory to be used in gaussian input notation.
  * atoms  -- The atoms to be computed. One can specify ranges, e.g 1,2-5'
              (avoid whitespace)

It is wise to run this script in a directory that is initially empty and that
will contain nothing but the generated atom database. This script will generate
quite a few files and subdirectories.

Examples:

%prog MP2/Aug-CC-pVDZ 1-10,17
%prog HF/3-21G 1,6,7,8 -l 110
""")
parser.add_option(
    "--density",
    help="The density field to use from the gaussian fchk file (scf, mp2, mp3, "
    "...). If not given, the program will guess it based on the level of "
    "theory."
)
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
    help="Specify the qc convergence scheme in gaussian input. [default=%default]"
)
(options, args) = parser.parse_args()
if len(args) != 2:
    parser.error("Expecting two arguments: level of theory (+ basis set) and the atom specification.")
lot, atom_str = args

if options.density is None:
    options.density = guess_density_type(lot)
atom_numbers = parse_numbers(atom_str)

make_inputs(lot, atom_numbers, options.max_ion, options.qc)
run_jobs()
select_ground_states('%s Energy' % options.density.upper(), options.max_ion)
make_density_profile(
    options.density, options.lebedev, options.rlow*angstrom, options.rhigh*angstrom,
    options.num_steps, atom_numbers, options.max_ion
)

