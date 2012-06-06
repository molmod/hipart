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


from writer import *


class Gint2NAI(GaussianIntegral):
    def __init__(self):
        a = ShellArgGroup("a")
        b = ShellArgGroup("b")
        c = PointArgGroup("c")
        self.a = a.args[0].symbols
        self.a_a = a.args[1].symbol
        self.b = b.args[0].symbols
        self.b_a = b.args[1].symbol
        self.c = c.args[0].symbols
        name = "gint2_nai"
        includes = ["gaux.h"]
        GaussianIntegral.__init__(self, name, [a, b, c], includes)

    def get_SS_integral(self, routine, ab_a, ab_overlap, usq, m):
        result = 2*sqrt(ab_a/pi)*ab_overlap*C.Function("gaux")(usq, m)
        symbol = Symbol("nai_000_000_%i" % m)
        routine.add(symbol, result, "local")
        return symbol

    def get_integral(self, routine, a_n, b_n, u, v, w, ab_a, ab_overlap, usq, m):
        #print "get_integral: %s %s" % (a_n, b_n)
        # This the recurrence relation from the paper of Obara and Saika
        # (see http://dx.doi.org/10.1063/1.450106)
        if (a_n[0] < 0) or (a_n[1] < 0) or (a_n[2] < 0) or (b_n[0] < 0) or (b_n[1] < 0) or (b_n[2] < 0):
            # This is where the recurrence relations leads to zero terms
            return 0
        elif (a_n[0] == 0) and (a_n[1] == 0) and (a_n[2] == 0) and (b_n[0] == 0) and (b_n[1] == 0) and (b_n[2] == 0):
            # This is where the recurrence relations leads to the SS term
            return self.get_SS_integral(routine, ab_a, ab_overlap, usq, m)
        else:
            # apply the recurrence relation
            # Note: we must decide which of the six possible paths we take.
            # Current strategy: get rid of the largest power first
            nonzero_a = [n for n in a_n if n > 0]
            if len(nonzero_a) == 0:
                low_a = 0
            else:
                low_a = min(nonzero_a)
            nonzero_b = [n for n in b_n if n > 0]
            if len(nonzero_b) == 0:
                low_b = 0
            else:
                low_b = min(nonzero_b)
            if (low_a < low_b and low_a > 0) or low_b == 0:
                assert low_a > 0
                index = [i for i, n in enumerate(a_n) if n==low_a][0]
                a_n1 = list(a_n)
                a_n1[index] -= 1
                a_n1 = tuple(a_n1)
                b_n1 = list(b_n)
                b_n1[index] -= 1
                b_n1 = tuple(b_n1)
                a_n2 = list(a_n)
                a_n2[index] -= 2
                a_n2 = tuple(a_n2)
                result = mypowsimp((
                      v[index]*self.get_integral(routine, a_n1, b_n, u, v, w, ab_a, ab_overlap, usq, m)
                    - u[index]*self.get_integral(routine, a_n1, b_n, u, v, w, ab_a, ab_overlap, usq, m+1)
                ) + a_n1[index]/(2*ab_a)*(
                      self.get_integral(routine, a_n2, b_n, u, v, w, ab_a, ab_overlap, usq, m)
                    - self.get_integral(routine, a_n2, b_n, u, v, w, ab_a, ab_overlap, usq, m+1)
                ) + b_n[index]/(2*ab_a)*(
                      self.get_integral(routine, a_n1, b_n1, u, v, w, ab_a, ab_overlap, usq, m)
                    - self.get_integral(routine, a_n1, b_n1, u, v, w, ab_a, ab_overlap, usq, m+1)
                ))
            else:
                assert low_b > 0
                index = [i for i, n in enumerate(b_n) if n==low_b][0]
                b_n1 = list(b_n)
                b_n1[index] -= 1
                b_n1 = tuple(b_n1)
                a_n1 = list(a_n)
                a_n1[index] -= 1
                a_n1 = tuple(a_n1)
                b_n2 = list(b_n)
                b_n2[index] -= 2
                b_n2 = tuple(b_n2)
                result = mypowsimp((
                      w[index]*self.get_integral(routine, a_n, b_n1, u, v, w, ab_a, ab_overlap, usq, m)
                    - u[index]*self.get_integral(routine, a_n, b_n1, u, v, w, ab_a, ab_overlap, usq, m+1)
                ) + b_n1[index]/(2*ab_a)*(
                      self.get_integral(routine, a_n, b_n2, u, v, w, ab_a, ab_overlap, usq, m)
                    - self.get_integral(routine, a_n, b_n2, u, v, w, ab_a, ab_overlap, usq, m+1)
                ) + a_n[index]/(2*ab_a)*(
                      self.get_integral(routine, a_n1, b_n1, u, v, w, ab_a, ab_overlap, usq, m)
                    - self.get_integral(routine, a_n1, b_n1, u, v, w, ab_a, ab_overlap, usq, m+1)
                ))
        symbol = Symbol("nai_%i%i%i_%i%i%i_%i" % (
            a_n[0], a_n[1], a_n[2], b_n[0], b_n[1], b_n[2], m
        ))
        routine.add(symbol, result, "local")
        return symbol

    def get_key(self, st_row):
        return "_".join(get_shell_label(st) for st in st_row[:-1])

    def add_expressions(self, st_row, routine):
        self.out_counter = 0

        ab_a = Symbol("ab_a")
        routine.add(ab_a, self.a_a + self.b_a, "local")
        ab_b = self.a_a*self.b_a/ab_a

        p = symbol_vector("p")
        routine.add(p[0], (self.a_a*self.a[0] + self.b_a*self.b[0])/ab_a, "local")
        routine.add(p[1], (self.a_a*self.a[1] + self.b_a*self.b[1])/ab_a, "local")
        routine.add(p[2], (self.a_a*self.a[2] + self.b_a*self.b[2])/ab_a, "local")

        d = symbol_vector("d")
        routine.add(d[0], self.a[0] - self.b[0], "local")
        routine.add(d[1], self.a[1] - self.b[1], "local")
        routine.add(d[2], self.a[2] - self.b[2], "local")

        dsq = Symbol("dsq")
        routine.add(dsq, d[0]**2 + d[1]**2 + d[2]**2, "local")

        myexp = Symbol("myexp")
        routine.add(myexp, C.Function("exp")(-ab_b*dsq), "local")
        ab_overlap = sqrt(pi/ab_a)**3*myexp

        u = symbol_vector("u")
        routine.add(u[0], p[0] - self.c[0], "local")
        routine.add(u[1], p[1] - self.c[1], "local")
        routine.add(u[2], p[2] - self.c[2], "local")

        usq = Symbol("usq")
        routine.add(usq, ab_a*(u[0]**2 + u[1]**2 + u[2]**2), "local")

        v = symbol_vector("v")
        routine.add(v[0], p[0] - self.a[0], "local")
        routine.add(v[1], p[1] - self.a[1], "local")
        routine.add(v[2], p[2] - self.a[2], "local")

        w = symbol_vector("w")
        routine.add(w[0], p[0] - self.b[0], "local")
        routine.add(w[1], p[1] - self.b[1], "local")
        routine.add(w[2], p[2] - self.b[2], "local")

        if st_row == (0,0,0):
            nai = self.get_SS_integral(routine, ab_a, ab_overlap, usq, 0)
            nai /= get_cartesian_wfn_norm(self.a_a, 0, 0, 0).evalf()
            nai /= get_cartesian_wfn_norm(self.b_a, 0, 0, 0).evalf()
            out_symbol = self.get_out_symbol()
            routine.add(out_symbol, nai, "final")
        elif st_row[0] >= -1 and st_row[1] >= -1:
            # generate cartesian integral expressions
            a_ns = list(iter_cartesian_powers(st_row[0]))
            b_ns = list(iter_cartesian_powers(st_row[1]))
            nais = numpy.zeros((len(a_ns),len(b_ns)), dtype=object)
            for a_i, a_n in enumerate(a_ns):
                for b_i, b_n in enumerate(b_ns):
                    nai = self.get_integral(routine, a_n, b_n, u, v, w, ab_a, ab_overlap, usq, 0)
                    nai /= get_cartesian_wfn_norm(self.a_a, a_n[0], a_n[1], a_n[2]).evalf()
                    nai /= get_cartesian_wfn_norm(self.b_a, b_n[0], b_n[1], b_n[2]).evalf()
                    out_symbol = self.get_out_symbol()
                    routine.add(out_symbol, nai, "final")
                    nais[a_i,b_i] = nai
        else:
            raise NotImplementedError("No integrals with pure functions.")

    def write_routine(self, f_pyf, f_c, f_h, st_row):
        if st_row[0] < -1 or st_row[1] < -1:
            # convert the result of a Cartesian gaussian to pure one
            if st_row[0] < -1:
                # make linear combinations of rows
                other_st_row = (-st_row[0], st_row[1], 0)
                size_a = get_shell_dof(-st_row[0])
                size_b = get_shell_dof(st_row[1])
                nais = numpy.zeros((size_a,size_b), dtype=object)
                self.out_counter = 0
                for ia in xrange(size_a):
                    for ib in xrange(size_b):
                        nais[ia,ib] = self.get_out_symbol()
                lcs = get_poly_conversion(st_row[0])
                nais = numpy.dot(lcs.transpose(), nais)
                self.out_counter = 0
                routine = InplaceRoutine([])
                for nai in nais.flat:
                    routine.append(self.get_out_symbol(), nai, "final")
            else:
                # make linear combintions of columns
                other_st_row = (st_row[0], -st_row[1], 0)
                size_a = get_shell_dof(st_row[0])
                size_b = get_shell_dof(-st_row[1])
                nais = numpy.zeros((size_a,size_b), dtype=object)
                self.out_counter = 0
                for ia in xrange(size_a):
                    for ib in xrange(size_b):
                        nais[ia,ib] = self.get_out_symbol()
                lcs = get_poly_conversion(st_row[1])
                nais = numpy.dot(nais, lcs)
                self.out_counter = 0
                routine = InplaceRoutine([])
                for nai in nais.flat:
                    routine.append(self.get_out_symbol(), nai, "final")
            return self.write_inplace_routine(f_pyf, f_c, f_h, st_row, other_st_row, routine)
        elif st_row[0] > st_row[1]:
            # transpose the result of another routine
            other_st_row = (st_row[1], st_row[0], 0)
            size_a = get_shell_dof(st_row[0])
            size_b = get_shell_dof(st_row[1])
            out_permutation = numpy.zeros(size_a*size_b, int)
            for ia in xrange(size_a):
                for ib in xrange(size_b):
                    new = ib + size_b*ia
                    old = ia + size_a*ib
                    out_permutation[new] = old
            arg_permutation = [1, 0, 2]
            return self.write_permutation_routine(f_pyf, f_c, f_h, st_row, other_st_row, out_permutation, arg_permutation)
        else:
            return self.write_cse_routine(f_pyf, f_c, f_h, st_row)


def main():
    Gint2NAI().write(max_shell=3)


if __name__ == "__main__":
    main()
