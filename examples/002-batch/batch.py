#!/usr/bin/env python


def main():
    import glob
    from molmod import angstrom
    from hipart.atoms import AtomTable
    from hipart.context import Context, Options
    from hipart.schemes import HirshfeldScheme

    # load the atom table from the first example
    atom_table = AtomTable("../001-usage/atoms/densities.txt")

    # select formatted checkpoint files from the first example
    fns_fchk = glob.glob("../001-usage/*.fchk")

    f = file("all.txt", "w")
    # loop over all formatted checkpoint files
    for fn_fchk in fns_fchk:
        # setup a HiPart computation
        options = Options(lebedev=266, do_work=False, do_output=False)
        context = Context(fn_fchk, options)
        scheme = HirshfeldScheme(context, atom_table)
        # the following line does the computation
        scheme.do_charges()
        # print output to the file all.txt
        m = scheme.molecule
        print >> f, "XXX", m.size
        for i in xrange(m.size):
            print >> f, "%3i % 15.10f % 15.10f % 15.10f % 9.5f" % (
                m.numbers[i], m.coordinates[i,0]/angstrom,
                m.coordinates[i,1]/angstrom, m.coordinates[i,2]/angstrom,
                scheme.charges[i]
            )
    f.close()


if __name__ == "__main__":
    main()
