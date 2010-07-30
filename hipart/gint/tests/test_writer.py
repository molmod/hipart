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


from sympy import sin

from hipart.gint.writer import *


def test_cartesian_powers():
    assert(list(iter_cartesian_powers(0))==[(0,0,0)])
    assert(list(iter_cartesian_powers(1))==[(1,0,0), (0,1,0), (0,0,1)])
    assert(list(iter_cartesian_powers(2))==[(2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2)])

def test_norms():
    alpha = Symbol("alpha")
    a = simplify(mypowsimp(get_cartesian_wfn_norm(alpha, 0, 0, 0)))
    b = simplify(mypowsimp(get_pure_wfn_norm(alpha, 0)))*sqrt(4*pi)
    assert(a==b)
    a = simplify(mypowsimp(get_cartesian_wfn_norm(alpha, 1, 0, 0)))
    b = simplify(mypowsimp(get_pure_wfn_norm(alpha, 1)))*sqrt(4*pi/3)
    assert(a==b)

def test_commands1():
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    routine = CSERoutine()
    routine.add(Symbol("out"), x*x+y*y+z*z, "final")
    routine.substitute(Symbol("tmp1"), x*x, "locked")
    routine.substitute(Symbol("tmp2"), y*y, "locked")
    assert(len(routine.commands)==3)
    assert(routine.commands[0].symbol == Symbol("tmp1"))
    assert(routine.commands[1].symbol == Symbol("tmp2"))
    assert(routine.commands[2].symbol == Symbol("out"))
    assert(routine.commands[0].expr == x*x)
    assert(routine.commands[1].expr == y*y)
    assert(routine.commands[2].expr == Symbol("tmp1")+Symbol("tmp2")+z*z)
    assert(routine.commands[0].tag == "locked")
    assert(routine.commands[1].tag == "locked")
    assert(routine.commands[2].tag == "final")

def test_commands2():
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    routine = CSERoutine()
    routine.add(Symbol("out"), x*x+y*y+z*z, "final")
    routine.full()
    assert(len(routine.commands)==1)

def test_commands3():
    x = Symbol("x")
    routine = CSERoutine()
    routine.add(Symbol("out"), x*x+sin(x*x), "final")
    routine.full()
    for record in routine.commands:
        print record
    assert(len(routine.commands)==2)

def test_permutation_to_cycles():
    permutation = [0, 3, 6, 9, 1, 4, 7, 10, 2, 5, 8, 11]
    cycles = permutation_to_cycles(permutation)
    assert(len(cycles)==2)
    assert([1, 3, 9, 5, 4] in cycles)
    assert([2, 6, 7, 10, 8] in cycles)
