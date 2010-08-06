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


from writer import *


class Gint1Fn(GaussianIntegral):
    def __init__(self):
        a = ShellArgGroup("a")
        p = PointArgGroup("p")
        self.a = a.args[0].symbols
        self.a_a = a.args[1].symbol
        self.p =  p.args[0].symbols
        name = "gint1_fn"
        GaussianIntegral.__init__(self, name, [a, p], [])

    def get_key(self, st_row):
        return get_shell_label(st_row[0])

    def add_expressions(self, st_row, routine):
        self.out_counter = 0

        v = symbol_vector("v")
        routine.add(v[0], self.p[0] - self.a[0], "local")
        routine.add(v[1], self.p[1] - self.a[1], "local")
        routine.add(v[2], self.p[2] - self.a[2], "local")

        rsq = Symbol("rsq")
        routine.add(rsq, v[0]*v[0] + v[1]*v[1] + v[2]*v[2], "local")

        for poly, wfn_norm in get_polys(st_row[0], self.a_a, v):
            fn = mypowsimp(simplify(poly/wfn_norm))*C.Function("exp")(-self.a_a*rsq)
            out_symbol = self.get_out_symbol()
            routine.add(out_symbol, fn, "final")


def main():
    Gint1Fn().write(max_shell=3)


if __name__ == "__main__":
    main()
