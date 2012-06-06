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
