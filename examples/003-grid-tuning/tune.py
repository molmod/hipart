#!/usr/bin/env python


def run_atom_dbs():
    import os
    from hipart.atomdb import run_atomdb, Options
    from molmod.periodic import periodic

    # auxiliary function to setup an atomic database with a given radial grid
    # size
    def my_run_atomdb(size):
        directory = "atoms%03i" % size
        if os.path.isfile(os.path.join(directory, "densities.txt")):
            return
        lot = 'B3LYP/6-31g(d)'
        atom_numbers = [1, 6, 7, 8]
        options = Options(num_steps=size, do_work=False)
        run_atomdb("g03", lot, atom_numbers, options, "atoms050")


    # array with different sizes for the radial grid
    sizes = [50, 100, 200, 400, 600]

    # first run with the least radial grid points
    my_run_atomdb(sizes[0])

    # then make symlinks to reuse the atomic computations
    for size in sizes[1:]:
        directory = "atoms%03i" % size
        if not os.path.isdir(directory):
            os.mkdir(directory)
            for number in atom_numbers:
                name = "%03i%s" % (number, periodic[number].symbol)
                source = os.path.join("../atoms050", name)
                link_name = os.path.join(directory, name)
                os.symlink(source, link_name)

    # run atomdb in other directories
    for size in sizes[1:]:
        my_run_atomdb(size)

    # load the tables with atomic density profiles
    from hipart.atoms import AtomTable
    atom_tables = {}
    for size in sizes:
        densities_fn = os.path.join("atoms%03i" % size, "densities.txt")
        atom_tables[size] = AtomTable(densities_fn)
    return atom_tables


class Config(object):
    def __init__(self, size, lebedev, error, cost):
        self.size = size
        self.lebedev = lebedev
        self.error = error
        self.cost = cost
        self.optimal = False

    def __str__(self):
        return "size = %3i    lebedev = %3i    error = %.2e    cost = %.2e" % \
            (self.size, self.lebedev, self.error, self.cost)


def estimate_errors(atom_tables, fn_fchk):
    from hipart.context import Context, Options
    from hipart.schemes import HirshfeldScheme
    from hipart.log import log
    import numpy, time

    configs = []
    # The loops below run 10 jobs for each grid configuration, i.e. combination
    # of radial and angular grid. Due to the random rotations of the angular
    # integration grids, the resulting charges will differ in each run. The
    # standard devation on each charge over the 10 runs is computed, and then
    # the error over all atoms is averaged. This error is used as 'the' error on
    # the charges. The cost is the cpu time consumed for computing this error.
    for size, atom_table in sorted(atom_tables.iteritems()):
        for lebedev in 26, 50, 110, 170, 266:
            log.begin("Config size=%i lebedev=%i" % (size, lebedev))
            start = time.clock()
            all_charges = []
            for counter in xrange(10):
                log.begin("Sample counter=%i" % counter)
                # Make a new working environment for every sample
                options = Options(lebedev=lebedev, do_work=False,
                                  do_output=False)
                context = Context(fn_fchk, options)
                scheme = HirshfeldScheme(context, atom_table)
                scheme.do_charges()
                all_charges.append(scheme.charges)
                del scheme
                log.end()
            all_charges = numpy.array(all_charges)
            error = all_charges.std(axis=0).mean()
            cost = time.clock() - start
            log("error=%.2e cost=%.2e" % (error, cost))
            configs.append(Config(size, lebedev, error, cost))
            log.end()
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
    print "="*9, ("="*13 + " ")*(len(lebedevs)*2)
    print "lebedev  ", " ".join(("l=%i" % l).center(27) for l in lebedevs)
    print "-"*9, ("-"*27 + " ")*(len(lebedevs))
    for s in sizes:
        words = []
        for l in lebedevs:
            c = config_map[(s,l)]
            if c.optimal:
                word = "**%8.3e**" % c.error
                word += " **%8.3e**" % c.cost
            else:
                word = "  %8.3e  " % c.error
                word += "   %8.3e  " % c.cost
            word.replace("e+0", "e+").replace("e-0", "e-")
            words.append(word)
        print "%9s" % ("s=%i" % s), " ".join(words)
    print "="*9, ("="*13 + " ")*(len(lebedevs)*2)


if __name__ == '__main__':
    # 1) Generate databases of spherically averaged atomic densities at
    # different grid resolutions.
    atom_tables = run_atom_dbs()
    # 2) Use each radial grid to estimate the error on the Hirshfeld-I charges.
    configs = estimate_errors(atom_tables, "../001-usage/alanine/gaussian.fchk")
    # 3) Find the pareto front
    pareto_front(configs)
    # 4) Write a table in reStructuredText.
    write_table(configs)
