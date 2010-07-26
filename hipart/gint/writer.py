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



from sympy import simplify, Symbol, sqrt, cos, sin, exp, Ylm, I, Mul, pi, C, S
from sympy.printing.ccode import CCodePrinter
from sympy.printing.tree import print_tree
from sympy.utilities.iterables import numbered_symbols
from sympy.printing.precedence import precedence
from sympy.core.basic import S
from sympy.mpmath import fac2
from sympy.utilities.iterables import preorder_traversal

from hipart.gint.basis import get_shell_dof

# Sympy stuff

x = Symbol("x", real=True)
y = Symbol("y", real=True)
z = Symbol("z", real=True)
r = Symbol("r", real=True)


class Record(object):
    def __init__(self, symbol, expr, tag):
        self.symbol = symbol
        self.expr = expr
        self.tag = tag

    def __str__(self):
        return "%s = %s   {%s}" % (self.symbol, self.expr, self.tag)

    def copy(self):
        return Record(self.symbol, self.expr, self.tag)


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


class Commands(object):
    def __init__(self):
        self.records = []
        self.records_by_symbol = {}
        self.get_symbol = numbered_symbols("tmp")

    def add(self, record):
        if record.symbol in self.records_by_symbol:
            assert(record.expr == self.records_by_symbol[record.symbol].expr)
        else:
            print "RECORD", record
            self.records.append(record)
            self.records_by_symbol[record.symbol] = record

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

    def substitute(self, record, force=False):
        # Make substitutions
        first = None
        for p in xrange(len(self.records)):
            old_expr = self.records[p].expr
            for new_expr in self._iter_subs(old_expr, record.expr, record.symbol):
                if old_expr == new_expr:
                    continue
                if weigh(old_expr) > weigh(new_expr) or force:
                    self.records[p].expr = new_expr
                    if first is None:
                        first = p
                    break
        if first is None:
            raise ValueError("Failed to substitute: %s" % record)
        # insert new expression in records
        self.records.insert(first, record)

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
        highest_weight = 0
        iterators = []
        for record in self.records:
            #print "   ", record.expr
            my_pre_order = MyPreOrder(record.expr)
            iterators.append((my_pre_order, True))

        while len(iterators) > 0:
            my_pre_order, deeper = iterators.pop(0)
            subtree = my_pre_order.next(deeper)
            if subtree is None:
                continue
            #print "---   ", subtree,
            deeper = True
            for part in iter_parts(subtree, level):
                if part == record.expr:
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
                    #print "------   ", part, weight
                    count, weight = all_parts.get(part, (0, weight))
                    count += 1
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
                self.substitute(Record(symbol, part, "auto"))
                print "AUTOSUB-%i" % level, count, weight, part
                success = True
            except ValueError:
                pass
        return success

    def clean(self):
        # go in reverse order through all records
        counter = len(self.records)-1
        while counter >= 0:
            record = self.records[counter]
            if record.tag != "final":
                # do not remove end results
                # check for usage
                used = False
                for later_record in self.records[counter+1:]:
                    if record.symbol in later_record.expr:
                        used = True
                        break
                if not used:
                    print "CLEAN", record
                    del self.records[counter]
            counter -= 1

    def singles(self):
        counter = len(self.records)-1
        while counter >= 0:
            record = self.records[counter]
            if record.tag != "final":
                # do not consider end results
                # check for usage
                used = 0
                tmp = None
                for later_record in self.records[counter+1:]:
                    if record.symbol in later_record.expr:
                        # TODO: count the number of occurences and treat squares
                        # as double
                        for subtree in preorder_traversal(later_record.expr):
                            if subtree == record.symbol:
                                used += 1
                                tmp = later_record
                            if isinstance(subtree, C.Pow) and subtree.exp == 2 and subtree.base == record.symbol:
                                used += 1
                                tmp = later_record
                if used == 1:
                    tmp.expr = tmp.expr.subs(record.symbol, record.expr)
                    print "SINGLE", record
                    del self.records[counter]
            counter -= 1

    def full(self):
        # remove unused commands
        self.clean()
        # evalf
        for record in self.records:
            record.expr = mypowsimp(mycollectsimp(record.expr.evalf()))
            record.expr = record.expr.subs(C.Real(-1.0), -1)
            record.expr = record.expr.subs(C.Real(-2.0), -2)
        # substitute as much as possible
        while self.autosub(level=0):
            pass
        while self.autosub(level=1):
            pass
        # substitute back temporary variables that are used only once
        self.singles()
        # mypowsimp
        for record in self.records:
            record.expr = mypowsimp(record.expr)
            record.expr = record.expr.subs(C.Real(-1.0), -1)
            record.expr = record.expr.subs(C.Real(-2.0), -2)

    def get_recycle_records(self):
        result = [record.copy() for record in self.records]
        inuse = set([])
        for i, record in enumerate(result):
            if record.tag == "final":
                continue
            # Determine which symbols are no longer used after this record.
            avail = inuse.copy()
            for later_record in result[i+1:]:
                for symbol in list(avail):
                    if symbol in later_record.expr:
                        avail.discard(symbol)
            # If there are some, use it and substitute in later records
            if len(avail) > 0:
                symbol = sorted(avail)[0]
                for later_record in result[i+1:]:
                    later_record.expr = later_record.expr.subs(record.symbol, symbol)
                record.symbol = symbol
                record.tag += "+recycle"
                print "RECYCLE", record
            inuse.add(record.symbol)
        return result


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
    return weight


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

# Representation of Gaussian integral algorithm

class BaseArg(object):
    def __init__(self, name):
        self.name = name

    def get_c_arg(self):
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


class ScalarArg(BaseArg):
    def __init__(self, symbol):
        BaseArg.__init__(self, symbol.name)
        self.symbol = symbol

    def update_lookup(self, lookup):
        lookup[self.symbol.name] = self.symbol.name

    def get_c_type(self):
        return "double"


class BaseArgGroup(object):
    def __init__(self, prefix):
        self.prefix = prefix
        self.args = []

    def get_c_args(self):
        return ", ".join(arg.get_c_arg() for arg in self.args)

    def iter_shell_types(self, max_shell):
        raise NotImplementedError

    def get_switch_name(self):
        raise NotImplementedError


class PointArgGroup(BaseArgGroup):
    def __init__(self, prefix):
        BaseArgGroup.__init__(self, prefix)
        self.args.append(ArrayArg(self.prefix, [
            Symbol("%s_x" % prefix),
            Symbol("%s_y" % prefix),
            Symbol("%s_z" % prefix),
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

    def add_expressions(self, st_row, commands):
        raise NotImplementedError

    def get_c_args(self):
        return ", ".join(ag.get_c_args() for ag in self.arg_groups)


# Auxiliary functions

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


def get_polys(shell_type, alpha):
    result = []
    if shell_type < -1:
        shell_type = abs(shell_type)
        theta = Symbol("theta", real=True)
        phi = Symbol("phi", real=True)
        wfn_norm = get_pure_wfn_norm(alpha, shell_type)
        for m in range(shell_type+1):
            if m > 0:
                part_plus = (Ylm(shell_type,m,theta,phi)*r**shell_type).expand()
                part_plus = part_plus.subs(exp(I*phi),(x+I*y)/sin(theta)/r)
                part_min = (Ylm(shell_type,-m,theta,phi)*r**shell_type).expand()
                part_min = part_min.subs(exp(-I*phi),(x-I*y)/sin(theta)/r)
                sym = simplify(((-1)**m*part_plus+part_min)/sqrt(2))
                sym = sym.subs(r*cos(theta),z)
                sym = sym.subs(r**2,x**2+y**2+z**2)
                sym = sym.subs(cos(theta)**2,1-sin(theta)**2)
                sym = simplify(sym)
                asym = simplify(((-1)**m*part_plus-part_min)/I/sqrt(2))
                asym = asym.subs(r*cos(theta),z)
                asym = asym.subs(r**2,x**2+y**2+z**2)
                asym = asym.subs(cos(theta)**2,1-sin(theta)**2)
                asym = simplify(asym)
                result.append((sym, wfn_norm))
                result.append((asym, wfn_norm))
            else:
                part = (Ylm(shell_type,0,theta,phi)*r**shell_type).expand()
                part = part.subs(r*cos(theta),z)
                part = part.subs(r**2,x**2+y**2+z**2)
                part = simplify(part)
                result.append((part, wfn_norm))
    elif shell_type == -1:
        return get_polys(0, alpha) + get_polys(1, alpha)
    else:
        for l, m, n in iter_cartesian_powers(shell_type):
            result.append((x**l*y**m*z**n, get_cartesian_wfn_norm(alpha, l, m, n)))
    return result


# Simple example

code_gint1_fn_basis_c = """\
int gint1_fn_basis(double* weights, double* fns, double* points,
  double* centers, int* shell_types, int* shell_map,  int* num_primitives,
  double* ccoeffs, double* exponents, int num_weights, int num_points,
  int num_centers, int num_shells, int num_ccoeffs, int num_exponents)
{
  int shell, shell_type, primitive, dof, num_dof, result, i_point;
  double *work, *out, *center, *ccoeff, *exponent, *weight, *shell_weights;

  result = 0;

  work = malloc(MAX_SHELL_DOF*sizeof(double));
  if (work==NULL) {result = -1; goto EXIT;}

  for (i_point=0; i_point<num_points; i_point++) {
    *fns = 0.0;
    shell_weights = weights;
    ccoeff = ccoeffs;
    exponent = exponents;
    for (shell=0; shell<num_shells; shell++) {
      center = centers + (3*shell_map[shell]);
      shell_type = shell_types[shell];
      for (primitive=0; primitive<num_primitives[shell]; primitive++) {
        gint_fn_dispatch(shell_type, center, *exponent, points, work);
        //printf("shell_type=%d  primitive=%d  exponent=%f\\n", shell_type, primitive, *exponent);
        out = work;
        exponent++;
        weight = shell_weights;
        if (shell_type==-1) {
          *fns += (*weight)*(*out)*(*ccoeff);
          //printf("weight=%f  out=%f  ccoeff=%f  contrib=%f  fn=%f\\n", *weight, *out, *ccoeff, (*weight)*(*out)*(*ccoeff), *fns);
          weight++;
          out++;
          ccoeff++;
          num_dof = 3;
        } else if (shell_type > 0) {
          num_dof = ((shell_type+1)*(shell_type+2))/2;
        } else {
          num_dof = -2*shell_type+1;
        }
        for (dof=0; dof<num_dof; dof++) {
          *fns += (*weight)*(*out)*(*ccoeff);
          //printf("weight=%f  out=%f  ccoeff=%f  contrib=%f  fn=%f\\n", *weight, *out, *ccoeff, (*weight)*(*out)*(*ccoeff), *fns);
          weight++;
          out++;
        }
        ccoeff++;
      }
      if (shell_type==-1) {
        num_dof = 4;
      } else if (shell_type > 0) {
        num_dof = ((shell_type+1)*(shell_type+2))/2;
      } else {
        num_dof = -2*shell_type+1;
      }
      shell_weights += num_dof;
    }
    points += 3;
    fns++;
    //printf("\\n");
  }

EXIT:
  free(work);
  return result;
}"""


code_gint1_fn_basis_pyf = """\
  integer function gint1_fn_basis(weights, fns, points, centers, shell_types, shell_map, num_primitives, ccoeffs, exponents, num_weights, num_points, num_centers, num_shells, num_ccoeffs, num_exponents)
    intent(c) gint1_fn_basis
    intent(c)
    double precision intent(in) :: weights(num_weights)
    double precision intent(inout) :: fns(num_points)
    double precision intent(in) :: points(num_points,3)
    double precision intent(in)  :: centers(num_centers,3)
    integer intent(int) :: shell_types(num_shells)
    integer intent(int) :: shell_map(num_shells)
    integer intent(int) :: num_primitives(num_shells)
    double precision intent(in) :: ccoeffs(num_ccoeffs)
    double precision intent(in) :: exponents(num_exponents)
    integer intent(hide), depend(weights) :: num_weights=len(weights)
    integer intent(hide), depend(points) :: num_points=len(points)
    integer intent(hide), depend(centers) :: num_centers=len(centers)
    integer intent(hide), depend(shell_types) :: num_shells=len(shell_types)
    integer intent(hide), depend(ccoeffs) :: num_ccoeffs=len(ccoeffs)
    integer intent(hide), depend(exponents) :: num_exponents=len(exponents)
  end function gint1_fn_basis
"""


code_gint1_fn_dmat_c = """\
int gint1_fn_dmat(double* dmat, double* density, double* points,
  double* centers, int* shell_types, int* shell_map,  int* num_primitives,
  double* ccoeffs, double* exponents, int num_dmat, int num_points,
  int num_centers, int num_shells, int num_ccoeffs, int num_exponents)
{
  int shell, shell_type, primitive, dof1, dof2, num_dof, num_shell_dof, result, i_point;
  double *work, *out, *center, *ccoeff, *exponent, *basis_fns, *shell_fns, *fn, *dmat_element;

  result = 0;

  work = malloc(MAX_SHELL_DOF*sizeof(double));
  if (work==NULL) {result = -1; goto EXIT;}
  num_dof = ((int)(sqrt(1.0+8.0*num_dmat)-1.0))/2;
  basis_fns = malloc(num_dof*sizeof(double));
  if (basis_fns==NULL) {result = -1; goto EXIT;}

  for (i_point=0; i_point<num_points; i_point++) {
    // A) clear the basis functions.
    fn = basis_fns;
    for (dof1=0; dof1<num_dof; dof1++) {
      *fn = 0.0;
      fn++;
    }
    // B) evaluate the basis functions in the current point.
    ccoeff = ccoeffs;
    exponent = exponents;
    shell_fns = basis_fns;
    for (shell=0; shell<num_shells; shell++) {
      center = centers + (3*shell_map[shell]);
      shell_type = shell_types[shell];
      for (primitive=0; primitive<num_primitives[shell]; primitive++) {
        gint1_fn(shell_type, center, *exponent, points, work);
        //printf("shell_type=%d  primitive=%d  exponent=%f\\n", shell_type, primitive, *exponent);
        out = work;
        exponent++;
        fn = shell_fns;
        if (shell_type==-1) {
          *fn += (*out)*(*ccoeff);
          //printf("out=%f  ccoeff=%f  fn=%f\\n", *out, *ccoeff, (*out)*(*ccoeff));
          fn++;
          out++;
          ccoeff++;
          num_shell_dof = 3;
        } else if (shell_type > 0) {
          num_shell_dof = ((shell_type+1)*(shell_type+2))/2;
        } else {
          num_shell_dof = -2*shell_type+1;
        }
        for (dof1=0; dof1<num_shell_dof; dof1++) {
          *fn += (*out)*(*ccoeff);
          //printf("out=%f  ccoeff=%f  fn=%f\\n", *out, *ccoeff, (*out)*(*ccoeff));
          fn++;
          out++;
        }
        ccoeff++;
      }
      if (shell_type==-1) {
        num_shell_dof = 4;
      } else if (shell_type > 0) {
        num_shell_dof = ((shell_type+1)*(shell_type+2))/2;
      } else {
        num_shell_dof = -2*shell_type+1;
      }
      shell_fns += num_shell_dof;
    }
    //printf("\\n");
    // C) Make dot product of basis functions with density matrix.
    *density = 0.0;
    dmat_element = dmat;
    for (dof1=0; dof1<num_dof; dof1++) {
      for (dof2=0; dof2<=dof1; dof2++) {
        if (dof1==dof2) {
          *density += basis_fns[dof1]*basis_fns[dof2]*(*dmat_element);
        } else {
          *density += 2*basis_fns[dof1]*basis_fns[dof2]*(*dmat_element);
        }
        //printf("dof1=%d  dof2=%d  basis1=%f  basis2=%f  dmat_element=%f  density=%f\\n", dof1, dof2, basis_fns[dof1], basis_fns[dof2], *dmat_element, *density);
        dmat_element++;
      }
    }
    // D) Prepare for next iteration
    density++;
    points += 3;
    //printf("\\n");
  }

EXIT:
  free(work);
  free(basis_fns);
  return result;
}"""


code_gint1_fn_dmat_pyf = """\
  integer function gint1_fn_dmat(dmat, density, points, centers, shell_types, shell_map, num_primitives, ccoeffs, exponents, num_dmat, num_points, num_centers, num_shells, num_ccoeffs, num_exponents)
    intent(c) gint1_fn_dmat
    intent(c)
    double precision intent(in) :: dmat(num_dmat)
    double precision intent(inout) :: density(num_points)
    double precision intent(in) :: points(num_points,3)
    double precision intent(in)  :: centers(num_centers,3)
    integer intent(int) :: shell_types(num_shells)
    integer intent(int) :: shell_map(num_shells)
    integer intent(int) :: num_primitives(num_shells)
    double precision intent(in) :: ccoeffs(num_ccoeffs)
    double precision intent(in) :: exponents(num_exponents)
    integer intent(hide), depend(dmat) :: num_dmat=len(dmat)
    integer intent(hide), depend(points) :: num_points=len(points)
    integer intent(hide), depend(centers) :: num_centers=len(centers)
    integer intent(hide), depend(shell_types) :: num_shells=len(shell_types)
    integer intent(hide), depend(ccoeffs) :: num_ccoeffs=len(ccoeffs)
    integer intent(hide), depend(exponents) :: num_exponents=len(exponents)
  end function gint1_fn_dmat
"""


class Gint1Fn(GaussianIntegral):
    def __init__(self):
        a = ShellArgGroup("a")
        p = PointArgGroup("p")
        self.a_x = a.args[0].symbols[0]
        self.a_y = a.args[0].symbols[1]
        self.a_z = a.args[0].symbols[2]
        self.a_a = a.args[1].symbol
        self.p_x = p.args[0].symbols[0]
        self.p_y = p.args[0].symbols[1]
        self.p_z = p.args[0].symbols[2]
        name = "gint1_fn"
        interface_fns = [
            (code_gint1_basis_fn_c, code_gint1_basis_fn_pyf),
            (code_gint1_basis_dmat_c, code_gint1_basis_dmat_pyf),
        ]
        GaussianIntegral.__init__(self, name, [a, p], interface_fns)

    def get_expressions(self, st_row):
        results = []
        st = st_row[0]
        for poly, wfn_norm in get_polys(st, self.a_a):
            rsq = x*x + y*y + z*z
            expr = mypowsimp(simplify(poly/wfn_norm))*C.Function("exp")(-self.a_a*rsq)
            expr = expr.subs(x, self.p_x - self.a_x)
            expr = expr.subs(y, self.p_y - self.a_y)
            expr = expr.subs(z, self.p_z - self.a_z)
            results.append(expr)
        key = get_shell_label(st)
        return key, results, [self.p_x - self.a_x, self.p_y - self.a_y, self.p_z - self.a_z]


def write(gint, max_shell=3):
    # C source code
    f = open("%s.c" % gint.name, "w")
    f_header = open("../../HEADER.c")
    f.write(f_header.read())
    f_header.close()
    fn_names = []
    print >> f
    print >> f, "#include <math.h>"
    print >> f, "#include <stdlib.h>"
    for include in gint.includes:
        print >> f, "#include \"%s\"" % include
    print >> f, "#define MAX_SHELL %i" % max_shell
    print >> f, "#define NUM_SHELL_TYPES %i" % (2*max_shell+1)
    print >> f, "#define MAX_SHELL_DOF %i" % max(get_shell_dof(shell_type) for shell_type in xrange(-max_shell, max_shell+1))
    print >> f
    for st_row in gint.iter_shell_types(max_shell):
        key = gint.get_key(st_row)
        fn_name = "%s_%s" % (gint.name, key)
        fn_names.append(fn_name)
        print "Starting", fn_name
        commands = Commands()
        gint.add_expressions(st_row, commands)
        # use commands thing to do cse
        commands.full()
        records = commands.get_recycle_records()
        print "VARIABLES"
        print >> f, "static void %s(%s, double* out)" % (fn_name, gint.get_c_args())
        print >> f, "{"
        # variables
        variables = set([
            record.symbol.name for record in records
            if not record.symbol.name.startswith("out")
        ])
        if len(variables) > 0:
            print >> f, "  // Number of local variables: %i" % len(variables)
            print >> f, "  double %s;" % (", ".join(sorted(variables)))
        # code lines
        total_weight = 0
        for record in records:
            print "CCODE", record
            weight = weigh(record.expr)
            total_weight += weight
            print >> f, "  %s = %s; // %s, weighs %i" % (
                gint.lookup.get(record.symbol.name, record.symbol.name),
                ccode(record.expr, lookup=gint.lookup),
                record.tag, weight
            )
        print >> f, "  // total weight = %i" % total_weight
        print >> f, "}"
        print >> f
    # add some sugar
    arg_c_types = []
    for ag in gint.arg_groups:
        for arg in ag.args:
            arg_c_types.append(arg.get_c_type())
    print >> f, "typedef void (*fntype)(%s, double*);" % (", ".join(arg_c_types))
    print >> f, "const fntype fns[%i] = {%s};" % (len(fn_names), ", ".join(fn_names))
    print >> f
    all_args = []
    c_arg_names = []
    switches = []
    for ag in gint.arg_groups:
        switch = ag.get_switch_name()
        if switch is not None:
            all_args.append("int %s" % switch)
            switches.append(switch)
        for arg in ag.args:
            c_arg_names.append(arg.name)
        all_args.append(ag.get_c_args())
    print >> f, "void %s_dispatch(%s, double* out)" % (gint.name, ", ".join(all_args))
    print >> f, "{"
    factor = 1
    offsets = []
    for switch in switches:
        if factor == 1:
            offsets.append("%i+%s" % (max_shell, switch))
        else:
            offsets.append("%i*(%i+%s)" % (factor, max_shell, switch))
        factor *= 2*max_shell + 1
    print >> f, "  fns[%s](%s, out);" % ("+".join(offsets), ", ".join(c_arg_names))
    print >> f, "}"
    for interface_fn in gint.interface_fns:
        f.write("\n")
        f.write(interface_fn[0])
        f.write("\n")
    f.close()

    # F2PY Interface
    f = open("%s.pyf" % gint.name, "w")
    print >> f, "python module %s" % gint.name
    print >> f, "interface"
    print >> f
    for interface_fn in gint.interface_fns:
        print >> f, interface_fn[1]
    print >> f, "end interface"
    print >> f, "end python module %s" % gint.name
    f.close()


def main():
    write(Gint1Fn())


if __name__ == "__main__":
    main()
