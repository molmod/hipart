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


from sympy import Symbol, Wild, sqrt, pi, C, S, simplify
from sympy.printing.ccode import CCodePrinter
from sympy.utilities.iterables import numbered_symbols
from sympy.printing.precedence import precedence
from sympy.mpmath import fac2
from sympy.utilities.iterables import preorder_traversal

import numpy

from hipart.gint.basis import get_shell_dof
from hipart.gint.solid_harmonics import get_solid_harmonics


# Sympy stuff


class Command(object):
    def __init__(self, symbol, expr, tag):
        self.symbol = symbol
        self.expr = expr
        self.tag = tag

    def __str__(self):
        return "%s = %s   {%s}" % (self.symbol, self.expr, self.tag)

    def copy(self):
        return Command(self.symbol, self.expr, self.tag)


def iter_subsets(l):
    from itertools import combinations
    for size in xrange(2,len(l)):
        for comb in combinations(l, size):
            yield comb


class MyPreOrder(object):
    def __init__(self, expr):
        self.expr = expr
        self.state = -1
        self.counter = 0
        self.sub_iter = None

    def next(self, deeper=True):
        if self.state == -1:
            if deeper:
                self.state = 0
                return self.expr
            else:
                return None
        elif len(self.expr.args) > 0:
            if self.sub_iter is None:
                self.sub_iter = MyPreOrder(self.expr.args[0])
            while True:
                result = self.sub_iter.next(deeper)
                if result is None:
                    deeper = True
                    self.counter += 1
                    if self.counter >= len(self.expr.args):
                        break
                    self.sub_iter = MyPreOrder(self.expr.args[self.counter])
                else:
                    return result


class BaseRoutine(object):
    def __init__(self, commands):
        self.commands = commands
        self.get_symbol = numbered_symbols("tmp")

    size = property(lambda self: len(self.commands))

    def get_num_out(self):
        return sum(command.tag == "final" for command in self.commands)

    def clean(self):
        # go in reverse order through all commands and remove unused ones.
        counter = self.size-1
        while counter >= 0:
            command = self.commands[counter]
            if command.tag != "final":
                # do not remove end results
                # check for usage
                used = False
                for later_command in self.commands[counter+1:]:
                    if command.symbol in later_command.expr:
                        used = True
                        break
                if not used:
                    print "CLEAN", command
                    del self.commands[counter]
            counter -= 1


class CSERoutine(BaseRoutine):
    def __init__(self):
        self.commands_by_symbol = {}
        BaseRoutine.__init__(self, [])

    def add(self, symbol, expr, tag):
        if symbol in self.commands_by_symbol:
            assert(expr == self.commands_by_symbol[symbol].expr)
        else:
            command = Command(symbol, expr, tag)
            print "RECORD", command
            self.commands.append(command)
            self.commands_by_symbol[command.symbol] = command

    def _iter_subs(self, old_expr, sub, symbol):
        yield old_expr.subs(sub, symbol)
        if isinstance(sub, C.Pow):
            for power in old_expr.atoms(C.Pow):
                if power.base == sub.base:
                    if power.exp - sub.exp == 1 and power.exp > 1:
                        yield old_expr.subs(power, sub.base*symbol)
                    elif power.exp == 2*sub.exp:
                        yield old_expr.subs(power, sub*sub)
        yield old_expr.subs(-sub, -symbol)

    def substitute(self, symbol, expr, tag, force=False):
        # Make substitutions
        first = None
        for p in xrange(self.size):
            old_expr = self.commands[p].expr
            for new_expr in self._iter_subs(old_expr, expr, symbol):
                if old_expr == new_expr:
                    continue
                if weigh(old_expr) > weigh(new_expr) or force:
                    self.commands[p].expr = new_expr
                    if first is None:
                        first = p
                    break
        if first is None:
            raise ValueError("Failed to substitute")
        # insert new expression in commands
        self.commands.insert(first, Command(symbol, expr, tag))

    def autosub(self, level=0):
        def iter_parts(expr, level=0):
            yield expr
            if level > 0 and (isinstance(expr, C.Add) or isinstance(expr, C.Mul)):
                for new_args in iter_subsets(expr.args):
                    good = True
                    for arg in new_args:
                        if isinstance(arg, C.Number):
                            good = False
                            break
                    if good:
                        yield expr.__class__(*new_args)

        all_parts = {}
        highest_weight = None
        iterators = []
        for command in self.commands:
            #print "   ", command.expr
            my_pre_order = MyPreOrder(command.expr)
            iterators.append((my_pre_order, True))

        while len(iterators) > 0:
            my_pre_order, deeper = iterators.pop(0)
            subtree = my_pre_order.next(deeper)
            if subtree is None:
                continue
            #print "---   ", subtree
            deeper = (highest_weight is None)
            for part in iter_parts(subtree, level):
                if part == command.expr:
                    continue
                if isinstance(part, C.Mul) and part.args[0] == -1:
                    continue
                if isinstance(part, C.Pow) and part.exp == -1:
                    continue
                if level==0 and isinstance(part, C.Mul) and isinstance(part.args[0], C.Number):
                    continue
                weight = weigh(part)
                if weight > highest_weight:
                    deeper = True
                if weight > 0:
                    count, weight = all_parts.get(part, (0, weight))
                    count += 1
                    #print "------   ", part, weight, count
                    all_parts[part] = count, weight
                    if count > 1:
                        highest_weight = max(weight, highest_weight)
            #print deeper
            iterators.append((my_pre_order, deeper))
        all_parts = [
            (weight, count, part)
            for (part, (count, weight))
            in all_parts.iteritems()
            if count > 1
        ]
        all_parts.sort(reverse=True)
        success = False
        for weight, count, part in all_parts:
            if weight < highest_weight:
                break
            try:
                symbol = self.get_symbol.next()
                self.substitute(symbol, part, "auto")
                print "AUTOSUB-%i" % level, count, weight, part
                success = True
            except ValueError:
                pass
        return success

    def singles(self):
        counter = self.size-1
        while counter >= 0:
            command = self.commands[counter]
            if command.tag != "final":
                # do not consider end results
                # check for usage
                used = 0
                tmp = None
                for later_command in self.commands[counter+1:]:
                    if command.symbol in later_command.expr:
                        for subtree in preorder_traversal(later_command.expr):
                            if subtree == command.symbol:
                                used += 1
                                tmp = later_command
                            if isinstance(subtree, C.Pow) and subtree.exp == 2 and subtree.base == command.symbol:
                                used += 1
                                tmp = later_command
                if used == 1:
                    tmp.expr = tmp.expr.subs(command.symbol, command.expr)
                    print "SINGLE", command
                    del self.commands[counter]
            counter -= 1

    def full(self):
        # remove unused commands
        self.clean()
        # evalf
        for command in self.commands:
            command.expr = mypowsimp(mycollectsimp(command.expr.evalf()))
            command.expr = command.expr.subs(C.Real(-1.0), -1)
            command.expr = command.expr.subs(C.Real(-2.0), -2)
        # substitute as much as possible
        while self.autosub(level=0):
            pass
        while self.autosub(level=1):
            pass
        # substitute back temporary variables that are used only once
        self.singles()
        # mypowsimp
        for command in self.commands:
            command.expr = mypowsimp(command.expr)
            command.expr = command.expr.subs(C.Real(-1.0), -1)
            command.expr = command.expr.subs(C.Real(-2.0), -2)

    def recycled(self):
        commands = [command.copy() for command in self.commands]
        inuse = []
        for i, command in enumerate(commands):
            if command.tag == "final":
                continue
            # Determine a symbol that is no longer used after this command.
            avail = False
            for j in xrange(len(inuse)):
                pivot = inuse.pop(0)
                avail = True
                for later_command in commands[i+1:]:
                    if pivot in later_command.expr:
                        avail = False
                        break
                if avail:
                    break
                else:
                    inuse.append(pivot)
            # If there are some, use it and substitute in later commands
            if avail:
                for later_command in commands[i+1:]:
                    later_command.expr = later_command.expr.subs(command.symbol, pivot)
                command.symbol = pivot
                command.tag += "+recycle"
                print "RECYCLE", command
            inuse.append(command.symbol)
        return BaseRoutine(commands)


class InplaceRoutine(BaseRoutine):
    def append(self, symbol, expr, tag):
        self.commands.append(Command(symbol, expr, tag))

    def sort(self):
        i = 0
        old_penalty = None
        penalty = 0
        while True:
            command1 = self.commands[i]
            command2 = self.commands[i+1]
            # check if command i can be swapped with command i+1
            if command2.symbol not in command1.expr:
                # check if command i wants to move after command i+1
                move = False
                for j, command in enumerate(self.commands[i+1:]):
                    if command1.symbol in command.expr:
                        move = True
                        penalty += j+1
                if move:
                    # perform the swap
                    #print "SWAPPING", i, i+1
                    self.commands[i+1] = command1
                    self.commands[i] = command2
                    changed = True
            i = (i + 1) % (self.size - 1)
            if i == 0:
                #print "PENALTY", penalty
                if old_penalty is None or penalty != old_penalty:
                    old_penalty = penalty
                    penalty = 0
                else:
                    break

    def non_overlapping(self):
        self.sort()
        # introduce temporary vairables so that every in-place operation is
        # based on the original output array, and not on the already partially
        # modified output.
        commands = [command.copy() for command in self.commands]
        avail_tmp_vars = []
        i = 0
        while i < len(commands):
            command = commands[i]
            if command.tag == "temporary":
                avail_tmp_vars.append(command.expr)
                command.tag = "final"
            else:
                assign_to = command.symbol
                last_used = None
                for j, later_command in enumerate(commands[i+1:]):
                    if assign_to in later_command.expr:
                        last_used = i+j+1
                if last_used is not None:
                    # get a temporary var
                    if len(avail_tmp_vars) > 0:
                        tmp = avail_tmp_vars.pop(0)
                    else:
                        tmp = self.get_symbol.next()
                    command.symbol = tmp
                    command.tag = "stacked"
                    commands.insert(last_used+1, Command(assign_to, tmp, "temporary"))
            i += 1
        # remove assignment to itself
        commands = [command for command in commands if (command.symbol!=command.expr)]
        return BaseRoutine(commands)



def weigh(expr):
    ops = expr.count_ops()
    weight = 0
    if ops == 0:
        return weight
    if isinstance(ops, C.Add):
        ops = ops.args
    else:
        ops = [ops]
    for op in ops:
        if isinstance(op, C.Mul):
            weight += op.args[0]
        else:
            weight += 1
    for pow in expr.atoms(C.Pow):
        if not (pow.exp==-1 or pow.exp==2):
            weight += 1
    return int(weight)


class MyCCodePrinter(CCodePrinter):
    _default_settings = {
        "order": None,
        "full_prec": "auto",
        "lookup": {}
    }

    def _print_Symbol(self, expr):
        result = self._settings["lookup"].get(expr.name)
        if result is None:
            result = expr.name
        return result

    def _print_Pow(self, expr):
        PREC = precedence(expr)
        if expr.exp is S.NegativeOne:
            return '1.0/%s'%(self.parenthesize(expr.base, PREC))
        elif isinstance(expr.base, Symbol) and expr.exp == 2:
            tmp = self.parenthesize(expr.base, PREC)
            return tmp+"*"+tmp
        else:
            return 'pow(%s,%s)'%(self.parenthesize(expr.base, PREC),
                                 self.parenthesize(expr.exp, PREC))


def ccode(expr, **settings):
    return MyCCodePrinter(settings).doprint(expr)


def mypowsimp(expr):
    def find_double_pow(expr):
        for sub in preorder_traversal(expr):
            if isinstance(sub, C.Pow) and isinstance(sub.base, C.Pow):
                return sub
    while True:
        sub = find_double_pow(expr)
        if sub is None:
            break
        expr = expr.subs(sub, C.Pow(sub.base.base, sub.exp*sub.base.exp))
    return expr


def mycollectsimp(expr):
    from sympy import collect
    counts = {}
    for subtree in preorder_traversal(expr):
        if isinstance(subtree, Symbol):
            counts[subtree] = counts.get(subtree, 0)+1
    counts = [(count, symbol) for symbol, count in counts.iteritems()]
    counts.sort()
    for count, symbol in counts:
        expr = collect(expr, symbol)
    return expr


# Representation of Gaussian Integral Routines


class BaseArg(object):
    def __init__(self, name):
        self.name = name

    def get_c_type_name(self):
        return "%s %s" % (self.get_c_type(), self.name)


class ArrayArg(BaseArg):
    def __init__(self, name, symbols):
        BaseArg.__init__(self, name)
        self.symbols = symbols

    def update_lookup(self, lookup):
        for i, symbol in enumerate(self.symbols):
            lookup[symbol.name] = "%s[%i]" % (self.name, i)

    def get_c_type(self):
        return "double*"

    def get_pyf_name(self):
        return "%s(%i)" % (self.name, len(self.symbols))


class ScalarArg(BaseArg):
    def __init__(self, symbol):
        BaseArg.__init__(self, symbol.name)
        self.symbol = symbol

    def update_lookup(self, lookup):
        lookup[self.symbol.name] = self.symbol.name

    def get_c_type(self):
        return "double"

    def get_pyf_name(self):
        return self.name


class BaseArgGroup(object):
    def __init__(self, prefix):
        self.prefix = prefix
        self.args = []

    def get_c_types_names(self):
        return ", ".join(arg.get_c_type_name() for arg in self.args)

    def get_c_names(self):
        return ", ".join(arg.name for arg in self.args)

    def get_pyf_names(self):
        return ", ".join(arg.get_pyf_name() for arg in self.args)

    def iter_shell_types(self, max_shell):
        raise NotImplementedError

    def get_switch_name(self):
        raise NotImplementedError


class PointArgGroup(BaseArgGroup):
    def __init__(self, prefix):
        BaseArgGroup.__init__(self, prefix)
        self.args.append(ArrayArg(self.prefix, [
            Symbol("%s_0" % prefix),
            Symbol("%s_1" % prefix),
            Symbol("%s_2" % prefix),
        ]))

    def iter_shell_types(self, max_shell):
        yield 0

    def get_switch_name(self):
        return None


class ShellArgGroup(PointArgGroup):
    def __init__(self, prefix):
        PointArgGroup.__init__(self, prefix)
        self.args.append(ScalarArg(Symbol("%s_a" % prefix, positive=True)))

    def iter_shell_types(self, max_shell):
        for st in xrange(-max_shell, max_shell+1):
            yield st

    def get_switch_name(self):
        return "%s_s" % self.prefix


class GaussianIntegral():
    def __init__(self, name, arg_groups, interface_fns, includes):
        self.name = name
        self.arg_groups = arg_groups
        self.interface_fns = interface_fns
        self.includes = includes
        self.lookup = {}
        self.out_counter = 0
        for ag in self.arg_groups:
            for arg in ag.args:
                arg.update_lookup(self.lookup)

    def get_out_symbol(self):
        i = self.out_counter
        symbol = Symbol("out_%i" % i)
        self.lookup[symbol.name] = "out[%i]" % i
        self.out_counter += 1
        return symbol

    def iter_shell_types(self, max_shell, arg_groups=None):
        if arg_groups is None:
            arg_groups = self.arg_groups
        if len(arg_groups) == 0:
            yield tuple()
        else:
            for result in self.iter_shell_types(max_shell, arg_groups[1:]):
                for st in arg_groups[0].iter_shell_types(max_shell):
                    yield (st,) + result

    def get_key(self, st_row):
        raise NotImplementedError

    def add_expressions(self, st_row, routine):
        raise NotImplementedError

    def symbol_to_c(self, symbol):
        return self.lookup.get(symbol.name, symbol.name)

    def get_c_types_names(self):
        return ", ".join(ag.get_c_types_names() for ag in self.arg_groups)

    def get_c_names(self, permutation=None):
        if permutation is None:
            return ", ".join(ag.get_c_names() for ag in self.arg_groups)
        else:
            return ", ".join(self.arg_groups[p].get_c_names() for p in permutation)

    def get_pyf_names(self):
        return ", ".join(ag.get_pyf_names() for ag in self.arg_groups)

    def write(self, max_shell=3):
        # C source code
        f_c = open("%s.c" % self.name, "w")
        f_h = open("%s.h" % self.name, "w")
        f_pyf = open("%s.pyf.inc" % self.name, "w")
        f_c_header = open("../../HEADER.c")
        f_c.write(f_c_header.read())
        f_c_header.seek(0)
        f_h.write(f_c_header.read())
        f_c_header.close()
        f_f_header = open("../../HEADER.f")
        f_pyf.write(f_f_header.read())
        f_f_header.close()
        fn_names = []
        print >> f_c
        print >> f_c, "#include <math.h>"
        print >> f_c, "#include <stdlib.h>"
        print >> f_c, "#include \"%s.h\"" % self.name
        for include in self.includes:
            print >> f_c, "#include \"%s\"" % include
        print >> f_c, "#define MAX_SHELL %i" % max_shell
        print >> f_c, "#define NUM_SHELL_TYPES %i" % (2*max_shell+1)
        print >> f_c, "#define MAX_SHELL_DOF %i" % max(get_shell_dof(shell_type) for shell_type in xrange(-max_shell, max_shell+1))
        print >> f_c, "#define CHECK_ALLOC(pointer) if (pointer==NULL) {result = -1; goto EXIT; }"
        print >> f_c, "#define CHECK_SHELL(shell_type) if (abs(shell_type) > MAX_SHELL) { result = -2; goto EXIT; }"
        print >> f_c
        for st_row in self.iter_shell_types(max_shell):
            print st_row
            fn_name = self.write_routine(f_pyf, f_c, f_h, st_row)
            fn_names.append(fn_name)
        # add some sugar
        arg_c_types = []
        for ag in self.arg_groups:
            for arg in ag.args:
                arg_c_types.append(arg.get_c_type())
        print >> f_c, "typedef void (*fntype)(%s, double*);" % (", ".join(arg_c_types))
        print >> f_c, "static const fntype fns[%i] = {%s};" % (len(fn_names), ", ".join(fn_names))
        print >> f_c
        all_c_types_names = []
        c_names = []
        switches = []
        for ag in self.arg_groups:
            switch = ag.get_switch_name()
            if switch is not None:
                all_c_types_names.append("int %s" % switch)
                switches.append(switch)
            for arg in ag.args:
                c_names.append(arg.name)
            all_c_types_names.append(ag.get_c_types_names())
        print >> f_c, "static void %s_dispatch(%s, double* out)" % (self.name, ", ".join(all_c_types_names))
        print >> f_c, "{"
        factor = 1
        offsets = []
        for switch in switches:
            if factor == 1:
                offsets.append("%i+%s" % (max_shell, switch))
            else:
                offsets.append("%i*(%i+%s)" % (factor, max_shell, switch))
            factor *= 2*max_shell + 1
        print >> f_c, "  fns[%s](%s, out);" % ("+".join(offsets), ", ".join(c_names))
        print >> f_c, "}"
        for interface_fn in self.interface_fns:
            f_c.write("\n")
            f_c.write(interface_fn[0])
            f_c.write("\n")
            print >> f_pyf, interface_fn[1]
        f_c.close()
        f_h.close()
        f_pyf.close()

    def write_routine_proto(self, f_pyf, f_c, f_h, fn_name, num_out):
        c_names = self.get_c_names()
        pyf_names = self.get_pyf_names()
        c_types_names = self.get_c_types_names()
        # prototype definitions in pyf
        print >> f_pyf, "  subroutine %s(%s, out)" % (fn_name, c_names)
        print >> f_pyf, "    intent(c) %s" % fn_name
        print >> f_pyf, "    intent(c)"
        print >> f_pyf, "    double precision intent(in) :: %s" % pyf_names
        print >> f_pyf, "    double precision intent(inout) :: out(%i)" % num_out
        print >> f_pyf, "  end subroutine %s" % fn_name
        print >> f_pyf
        # prototype definitions in h
        print >> f_h, "void %s(%s, double* out);" % (fn_name, c_types_names)
        # routine itself
        print >> f_c, "void %s(%s, double* out)" % (fn_name, c_types_names)
        print >> f_c, "{"

    def write_local_variables(self, f_c, routine):
        # variables
        variables = set([
            command.symbol.name for command in routine.commands
            if not command.tag == "final"
        ])
        print "VARIABLES", len(variables)
        if len(variables) > 0:
            print >> f_c, "  // Number of local variables: %i" % len(variables)
            print >> f_c, "  double %s;" % (", ".join(sorted(variables)))

    def write_commands(self, f_c, commands):
        total_weight = 0
        for command in commands:
            print "CCODE", command
            weight = weigh(command.expr)
            total_weight += weight
            print >> f_c, "  %s = %s; // %s, weighs %i" % (
                self.symbol_to_c(command.symbol),
                ccode(command.expr, lookup=self.lookup),
                command.tag, weight
            )
        print >> f_c, "  // total weight = %i" % total_weight

    def write_routine(self, f_pyf, f_c, f_h, st_row):
        return self.write_cse_routine(f_pyf, f_c, f_h, st_row)

    def write_cse_routine(self, f_pyf, f_c, f_h, st_row):
        fn_name = "%s_%s" % (self.name, self.get_key(st_row))
        routine = CSERoutine()
        self.add_expressions(st_row, routine)
        # use routine thing to do cse
        routine.full()
        routine = routine.recycled()
        num_out = routine.get_num_out()
        self.write_routine_proto(f_pyf, f_c, f_h, fn_name, num_out)
        print >> f_c, "  // Generated code based on 'Common SubExpression' analysis."
        self.write_local_variables(f_c, routine)
        self.write_commands(f_c, routine.commands)
        print >> f_c, "}"
        print >> f_c
        return fn_name

    def write_permutation_routine(self, f_pyf, f_c, f_h, st_row, other_st_row, out_permutation, arg_permutation):
        """Permutation of the output of another routine."""
        fn_name = "%s_%s" % (self.name, self.get_key(st_row))
        other_fn_name = "%s_%s" % (self.name, self.get_key(other_st_row))
        num_out = len(out_permutation)
        self.write_routine_proto(f_pyf, f_c, f_h, fn_name, num_out)
        cycles = permutation_to_cycles(out_permutation)
        print >> f_c, "  // In-place permutation of the output of another routine."
        if len(cycles) > 0:
            print >> f_c, "  double tmp;"
        print >> f_c, "  %s(%s, out);" % (other_fn_name, self.get_c_names(arg_permutation))
        if len(cycles) > 0:
            for cycle in cycles:
                print >> f_c, "  tmp = out[%i];" % cycle[0]
                for i in xrange(1, len(cycle)):
                    print >> f_c, "  out[%i] = out[%i];" % (cycle[i-1], cycle[i])
                print >> f_c, "  out[%i] = tmp;" % cycle[-1]
        print >> f_c, "}"
        print >> f_c
        return fn_name

    def write_inplace_routine(self, f_pyf, f_c, f_h, st_row, other_st_row, inplace_routine):
        fn_name = "%s_%s" % (self.name, self.get_key(st_row))
        other_fn_name = "%s_%s" % (self.name, self.get_key(other_st_row))
        num_out = inplace_routine.get_num_out()
        self.write_routine_proto(f_pyf, f_c, f_h, fn_name, num_out)
        print >> f_c, "  // In-place modification of the output of another routine."
        routine = inplace_routine.non_overlapping()
        self.write_local_variables(f_c, routine)
        # call the other function
        print >> f_c, "  %s(%s, out);" % (other_fn_name, self.get_c_names())
        self.write_commands(f_c, routine.commands)
        print >> f_c, "}"
        print >> f_c
        return fn_name


# Auxiliary functions


def permutation_to_cycles(permutation):
    # the permutation is an array with the indexes of the old positions.
    # the cycles are lists of indexes where every old index is preceded by its
    # new index in a cyclic fashion.
    # turn the permutation into a dictionary
    permutation = dict((new, old) for new, old in enumerate(permutation))
    # turn the permutation into cycles
    cycles = []
    current = None
    while len(permutation) > 0:
        if current is None:
            new, old = permutation.popitem()
            current = [new]
        else:
            old = permutation.pop(current[-1])
        if old == current[0]:
            if len(current) > 1:
                cycles.append(current)
            current = None
        else:
            current.append(old)
    return cycles


def symbol_vector(prefix):
    return (Symbol("%s_0" % prefix), Symbol("%s_1" % prefix), Symbol("%s_2" % prefix))


def get_shell_label(shell_type):
    shell_labels = "SPDFGHIJKLMNO"
    if shell_type > 1:
        return "c" + shell_labels[shell_type]
    elif shell_type >= 0:
        return shell_labels[shell_type]
    elif shell_type < -1:
        return "p" + shell_labels[-shell_type]
    else:
        return "SP"


def get_cartesian_wfn_norm(alpha, l, m, n):
    return sqrt(
        int(fac2(2*l-1)*fac2(2*m-1)*fac2(2*n-1))/
        sqrt((2*alpha/pi)**3*2**(4*(l+m+n))*alpha**(2*(l+m+n)))
    )


def get_pure_wfn_norm(alpha, l):
    return sqrt(int(fac2(2*l+1))*sqrt(pi/alpha/2)**3/(4*alpha)**l/(4*pi))


def iter_cartesian_powers(order):
    if order == -1:
        yield (0,0,0)
        yield (1,0,0)
        yield (0,1,0)
        yield (0,0,1)
    else:
        for l in xrange(order,-1,-1):
            for m in xrange(order-l,-1,-1):
                yield l, m, order-l-m


def get_polys(shell_type, alpha, xyz):
    result = []
    if shell_type < -1:
        shell_type = abs(shell_type)
        wfn_norm = get_pure_wfn_norm(alpha, shell_type)
        polys = get_solid_harmonics(shell_type, xyz)
        result = [(poly, wfn_norm) for poly in polys]
    elif shell_type == -1:
        return get_polys(0, alpha, xyz) + get_polys(1, alpha, xyz)
    else:
        x, y, z = xyz
        for l, m, n in iter_cartesian_powers(shell_type):
            result.append((x**l*y**m*z**n, get_cartesian_wfn_norm(alpha, l, m, n)))
    return result


def get_poly_conversion(shell_type):
    # conversion from normalised Cartesian to normalized pure
    assert shell_type < -1
    alpha = Symbol("alpha")
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    xyz = (x,y,z)

    px = Wild("px")
    py = Wild("py")
    pz = Wild("pz")
    c = Wild("c")

    num_dof_in = get_shell_dof(-shell_type)
    num_dof_out = get_shell_dof(shell_type)
    lcs = numpy.zeros((num_dof_in, num_dof_out), dtype=object)

    for i_out, (poly, pure_wfn_norm) in enumerate(get_polys(shell_type, alpha, xyz)):
        poly = poly.expand()
        if isinstance(poly, C.Add):
            terms = poly.args
        else:
            terms = [poly]
        coeffs = {}
        for term in terms:
            d = term.evalf(20).match(c*x**px*y**py*z**pz)
            key = (int(d[px]), int(d[py]), int(d[pz]))
            coeffs[key] = d[c]
        for i_in, key in enumerate(iter_cartesian_powers(abs(shell_type))):
            cart_wfn_norm = get_cartesian_wfn_norm(alpha, key[0], key[1], key[2])
            norm_ratio = mypowsimp((cart_wfn_norm/pure_wfn_norm).evalf(20))
            lc = float(coeffs.get(key, 0)*norm_ratio)
            if abs(lc-int(lc)) < 1e-12:
                lc = int(lc)
            lcs[i_in,i_out] = lc
    return lcs
