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



from sympy import simplify, Symbol, sqrt, cos, sin, exp, Ylm, I, Mul, cse, pi, Rational, pprint
from sympy.printing.ccode import CCodePrinter
from sympy.utilities.iterables import numbered_symbols
from sympy.printing.precedence import precedence
from sympy.core.basic import S
from sympy.mpmath import fac2

from hipart.gint.basis import get_shell_dof

# Sympy stuff

x = Symbol("x", real=True)
y = Symbol("y", real=True)
z = Symbol("z", real=True)
r = Symbol("r", real=True)


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
    from sympy.utilities.iterables import postorder_traversal
    from sympy import Pow
    def find_double_pow(expr):
        for sub in postorder_traversal(expr):
            if isinstance(sub, Pow) and isinstance(sub.base, Pow):
                return sub
    while True:
        sub = find_double_pow(expr)
        if sub is None:
            break
        expr = expr.subs(sub, Pow(sub.base.base, sub.exp*sub.base.exp))
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
    def __init__(self, name, arg_groups, interface_fns):
        self.name = name
        self.arg_groups = arg_groups
        self.interface_fns = interface_fns
        self.lookup = {}
        for ag in self.arg_groups:
            for arg in ag.args:
                arg.update_lookup(self.lookup)

    def iter_shell_types(self, max_shell, arg_groups=None):
        if arg_groups is None:
            arg_groups = self.arg_groups
        if len(arg_groups) == 0:
            yield tuple()
        else:
            for result in self.iter_shell_types(max_shell, arg_groups[1:]):
                for st in arg_groups[0].iter_shell_types(max_shell):
                    yield (st,) + result

    def get_expressions(self, st_row):
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


def get_cartesian_powers(order):
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
                print "sym", m
                pprint(sym)
                asym = simplify(((-1)**m*part_plus-part_min)/I/sqrt(2))
                asym = asym.subs(r*cos(theta),z)
                asym = asym.subs(r**2,x**2+y**2+z**2)
                asym = asym.subs(cos(theta)**2,1-sin(theta)**2)
                asym = simplify(asym)
                print "asym", m
                pprint(asym)
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
        for l, m, n in get_cartesian_powers(shell_type):
            result.append((x**l*y**m*z**n, get_cartesian_wfn_norm(alpha, l, m, n)))
    return result


# Simple example

template_shell1_point_sum_c = """\
int %(gint_name)s_basis(double* weights, double* fns, double* points,
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
        %(gint_name)s(shell_type, center, *exponent, points, work);
        //printf("shell_type=%%d  primitive=%%d  exponent=%%f\\n", shell_type, primitive, *exponent);
        out = work;
        exponent++;
        weight = shell_weights;
        if (shell_type==-1) {
          *fns += (*weight)*(*out)*(*ccoeff);
          //printf("weight=%%f  out=%%f  ccoeff=%%f  contrib=%%f  fn=%%f\\n", *weight, *out, *ccoeff, (*weight)*(*out)*(*ccoeff), *fns);
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
          //printf("weight=%%f  out=%%f  ccoeff=%%f  contrib=%%f  fn=%%f\\n", *weight, *out, *ccoeff, (*weight)*(*out)*(*ccoeff), *fns);
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


template_shell1_point_sum_pyf = """\
  integer function %(gint_name)s_basis(weights, fns, points, centers, shell_types, shell_map, num_primitives, ccoeffs, exponents, num_weights, num_points, num_centers, num_shells, num_ccoeffs, num_exponents)
    intent(c) %(gint_name)s_basis
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
  end function %(gint_name)s_basis
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
        interface_fns = [(
            template_shell1_point_sum_c % {"gint_name": name},
            template_shell1_point_sum_pyf % {"gint_name": name},
        )]
        GaussianIntegral.__init__(self, name, [a, p], interface_fns)

    def get_expressions(self, st_row):
        results = []
        st = st_row[0]
        print st
        for poly, wfn_norm in get_polys(st, self.a_a):
            rsq = x*x + y*y + z*z
            #print simplify(poly/wfn_norm)
            expr = (mypowsimp(simplify(poly/wfn_norm))).evalf()*exp(-self.a_a*rsq)
            expr = expr.subs(x, self.p_x - self.a_x)
            expr = expr.subs(y, self.p_y - self.a_y)
            expr = expr.subs(z, self.p_z - self.a_z)
            results.append(expr)
        key = get_shell_label(st)
        return key, results


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
    print >> f, "#define MAX_SHELL %i" % max_shell
    print >> f, "#define NUM_SHELL_TYPES %i" % (2*max_shell+1)
    print >> f, "#define MAX_SHELL_DOF %i" % (((max_shell+1)*(max_shell+2))/2)
    print >> f
    for st_row in gint.iter_shell_types(max_shell):
        key, results = gint.get_expressions(st_row)
        fn_name = "fn_%s" % key
        fn_names.append(fn_name)
        print >> f, "static void %s(%s, double* out)" % (fn_name, gint.get_c_args())
        print >> f, "{"
        # some cse comments
        temporaries, final_results = cse(results, numbered_symbols("tmp"))
        if len(temporaries) > 0:
            print >> f, "  double %s;" % (", ".join(var.name for var, expr in temporaries))
            for var, expr in temporaries:
                print >> f, "  %s = %s;" % (var.name, ccode(expr, lookup=gint.lookup))
        for i, expr in enumerate(final_results):
            print >> f, "  out[%i] = %s;" % (i, ccode(expr, lookup=gint.lookup))
        ## actual code
        #for i, expr in enumerate(results):
        #    print >> f, "  out[%i] = %s;" % (i, ccode(expr, lookup=gint.lookup))
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
    print >> f, "void gint1_fn(%s, double* out)" % (", ".join(all_args))
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
    print >> f
    for interface_fn in gint.interface_fns:
        f.write(interface_fn[0])
        f.write("\n")
    f.close()

    # F2PY Interface
    f = open("%s.pyf" % gint.name, "w")
    print >> f, "python module gint1_fn"
    print >> f, "interface"
    print >> f
    for interface_fn in gint.interface_fns:
        print >> f, interface_fn[1]
    print >> f, "end interface"
    print >> f, "end python module gint1_fn"
    f.close()


def main():
    write(Gint1Fn())


if __name__ == "__main__":
    main()
