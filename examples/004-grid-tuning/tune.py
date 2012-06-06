#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HiPart is a program to analyze the electronic structure of molecules with
# fuzzy-atom partitioning methods.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
#--
#!/usr/bin/env python


def run_atom_dbs():
    import os
    from hipart.atomdb import run_atomdb, Options
    from molmod.periodic import periodic

    # atom numbers of interest
    atom_numbers = [1, 6, 8]

    # auxiliary function to setup an atomic database with a given radial grid
    # size
    def my_run_atomdb(size):
        directory = "atoms%04i" % size
        if os.path.isfile(os.path.join(directory, "densities.txt")):
            return
        lot = 'B3LYP/6-31g(d)'
        options = Options(num_steps=size, do_work=False, do_random=True)
        run_atomdb("g03", lot, atom_numbers, options, directory)


    # array with different sizes for the radial grid
    sizes = [50, 100, 200, 400, 800, 1600]

    # first run with the least radial grid points
    my_run_atomdb(sizes[0])

    # then make symlinks to reuse the atomic computations
    for size in sizes[1:]:
        directory = "atoms%04i" % size
        if not os.path.isdir(directory):
            os.mkdir(directory)
            for number in atom_numbers:
                name = "%03i%s" % (number, periodic[number].symbol)
                source = os.path.join("../atoms0050", name)
                link_name = os.path.join(directory, name)
                os.symlink(source, link_name)

    # run atomdb in other directories
    for size in sizes[1:]:
        my_run_atomdb(size)

    # load the tables with atomic density profiles
    from hipart.atoms import AtomTable
    atom_tables = {}
    for size in sizes:
        densities_fn = os.path.join("atoms%04i" % size, "densities.txt")
        atom_tables[size] = AtomTable(densities_fn)
    return atom_tables


class Config(object):
    def __init__(self, size, lebedev, charges, cost):
        self.size = size
        self.lebedev = lebedev
        self.charges = charges
        self.error = 0.0
        self.cost = cost
        self.optimal = False

    #def prnt(self):
    #    array_str = "numpy.array([%s])" % (",".join(str(c) for c in charges))
    #    print "Config(size=%i, lebedev=%i, charges=%s, error=%.e, cost=%.e)," % \
    #        (self.size, self.lebedev, array_str, self.error, self.cost)

    def __str__(self):
        return "size = %3i    lebedev = %3i    error = %.2e    cost = %.2e" % \
            (self.size, self.lebedev, self.error, self.cost)


def estimate_errors(atom_tables, fns_fchk):
    from hipart.context import Context, Options
    from hipart.schemes import HirshfeldScheme
    from hipart.log import log
    import numpy, time

    configs = []
    ref_config = None
    # The loops below one jobs for each grid configuration, i.e. combination
    # of radial and angular grid.
    for size, atom_table in sorted(atom_tables.iteritems()):
        for lebedev in 26, 50, 110, 230, 434, 770, 1454, 2702:
            all_charges = []
            start = time.clock()
            for fn_fchk in fns_fchk:
                log.begin("Config size=%i lebedev=%i fchk=%s" % (size, lebedev, fn_fchk))
                # Make a new working environment for every sample
                options = Options(lebedev=lebedev, do_work=False,
                                  do_output=False, do_random=True)
                context = Context(fn_fchk, options)
                scheme = HirshfeldScheme(context, atom_table)
                scheme.do_charges()
                all_charges.append(scheme.charges)
                del scheme
                log.end()
            cost = time.clock() - start
            log("cost=%.2e" % cost)
            config = Config(size, lebedev, numpy.concatenate(all_charges), cost)
            configs.append(config)
            if ref_config is None or (ref_config.size <= config.size and ref_config.lebedev <= config.lebedev):
                ref_config = config

    # Compute the errors
    ref_charges = configs[-1].charges
    for config in configs:
        config.error = (ref_charges - config.charges).std()

    # Log all
    log.begin("All configs")
    for config in configs:
        log(str(config))
    log.end()
    return configs


def pareto_front(configs):
    from hipart.log import log
    # take a copy of the list, so that we can modify locally.
    configs = list(configs)
    # eliminate all configs that are not pareto-optimal
    configs.sort(key=(lambda c: c.error))
    i = 0
    while i < len(configs):
        ref_cost = configs[i].cost
        j = len(configs) -1
        while j > i:
            if configs[j].cost >= configs[i].cost:
                del configs[j]
            j -= 1
        i += 1
    # print the pareto front.
    log.begin("Pareto front")
    for config in configs:
        log(str(config))
        config.optimal = True
    log.end()


def write_table(configs):
    # This is an ulgy routine, and it is not so important as an illustration of
    # using HiPart as a library. It is usefull to write documentation though.
    # 1) Arrange everything in a convenient data format for making a table
    config_map = {}
    lebedevs = set([])
    sizes = set([])
    for config in configs:
        config_map[(config.size, config.lebedev)] = config
        lebedevs.add(config.lebedev)
        sizes.add(config.size)
    lebedevs = sorted(lebedevs)
    sizes = sorted(sizes)
    # 2) Write the table header
    print "="*9, ("="*10 + " ")*(len(lebedevs)*2)
    print "lebedev  ", " ".join(("l=%i" % l).center(21) for l in lebedevs)
    print "-"*9, ("-"*21 + " ")*(len(lebedevs))
    for s in sizes:
        words = []
        for l in lebedevs:
            c = config_map[(s,l)]
            if c.optimal:
                word = "**%7.1e** " % c.error
                word += ("**%.1f**" % c.cost).center(10)
            else:
                word = "  %7.1e   " % c.error
                word += ("%.1f" % c.cost).center(10)
            word = word.replace("e+0", "e+").replace("e-0", "e-")
            words.append(word)
        print "%9s" % ("s=%i" % s), " ".join(words)
    print "="*9, ("="*10 + " ")*(len(lebedevs)*2)


if __name__ == '__main__':
    import glob
    # 1) Generate databases of spherically averaged atomic densities at
    # different grid resolutions.
    atom_tables = run_atom_dbs()
    # 2) Use each radial grid to estimate the error on the Hirshfeld-I charges.
    configs = estimate_errors(atom_tables, glob.glob("*.fchk"))
    # 3) Find the pareto front
    pareto_front(configs)
    # 4) Write a table in reStructuredText.
    write_table(configs)
