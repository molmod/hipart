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
  CHECK_ALLOC(work);

  for (i_point=0; i_point<num_points; i_point++) {
    *fns = 0.0;
    shell_weights = weights;
    ccoeff = ccoeffs;
    exponent = exponents;
    for (shell=0; shell<num_shells; shell++) {
      center = centers + (3*shell_map[shell]);
      shell_type = shell_types[shell];
      CHECK_SHELL(shell_type);
      for (primitive=0; primitive<num_primitives[shell]; primitive++) {
        gint1_fn_dispatch(shell_type, center, *exponent, points, work);
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
  CHECK_ALLOC(work);
  num_dof = ((int)(sqrt(1.0+8.0*num_dmat)-1.0))/2;
  basis_fns = malloc(num_dof*sizeof(double));
  CHECK_ALLOC(basis_fns);

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
      CHECK_SHELL(shell_type);
      for (primitive=0; primitive<num_primitives[shell]; primitive++) {
        gint1_fn_dispatch(shell_type, center, *exponent, points, work);
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


code_gint1_fn_overlap_c="""\
int gint1_fn_overlap(double* overlap, double* points, double* weights,
  double* centers, int* shell_types, int* shell_map, int* num_primitives,
  double* ccoeffs, double* exponents, int num_overlap, int num_points,
  int num_centers, int num_shells, int num_ccoeffs, int num_exponents)
{
  int shell, shell_type, primitive, dof1, dof2, num_dof, num_shell_dof, result, i_point;
  double *work, *out, *center, *ccoeff, *exponent, *basis_fns, *shell_fns, *fn, *overlap_element;

  result = 0;

  work = malloc(MAX_SHELL_DOF*sizeof(double));
  CHECK_ALLOC(work);
  num_dof = (int)(sqrt(num_overlap));
  basis_fns = malloc(num_dof*sizeof(double));
  CHECK_ALLOC(basis_fns);

  // Clear the output
  out = overlap;
  for (dof1=0; dof1<num_overlap; dof1++) {
    *out = 0.0;
    out++;
  }

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
      CHECK_SHELL(shell_type);
      for (primitive=0; primitive<num_primitives[shell]; primitive++) {
        gint1_fn_dispatch(shell_type, center, *exponent, points, work);
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
    // C) Multiply overlap elements with the wieght and add to the atomic
    //    overlap matrix
    overlap_element = overlap;
    for (dof1=0; dof1<num_dof; dof1++) {
      for (dof2=0; dof2<num_dof; dof2++) {
        *overlap_element += (*weights)*basis_fns[dof1]*basis_fns[dof2];
        //printf("dof1=%d  dof2=%d  basis1=%f  basis2=%f  overlap_element=%f  weight=%f\\n", dof1, dof2, basis_fns[dof1], basis_fns[dof2], *overlap_element, *weight);
        overlap_element++;
      }
    }
    // D) Prepare for next iteration
    points += 3;
    weights++;
    //printf("\\n");
  }

EXIT:
  free(work);
  free(basis_fns);
  return result;
}"""


code_gint1_fn_overlap_pyf = """\
  integer function gint1_fn_overlap(overlap, points, weights, centers, shell_types, shell_map, num_primitives, ccoeffs, exponents, num_overlap, num_points, num_centers, num_shells, num_ccoeffs, num_exponents)
    intent(c) gint1_fn_overlap
    intent(c)
    double precision intent(inout) :: overlap(num_overlap)
    double precision intent(in) :: points(num_points,3)
    double precision intent(in) :: weights(num_points)
    double precision intent(in)  :: centers(num_centers,3)
    integer intent(int) :: shell_types(num_shells)
    integer intent(int) :: shell_map(num_shells)
    integer intent(int) :: num_primitives(num_shells)
    double precision intent(in) :: ccoeffs(num_ccoeffs)
    double precision intent(in) :: exponents(num_exponents)
    integer intent(hide), depend(overlap) :: num_overlap=len(overlap)
    integer intent(hide), depend(points) :: num_points=len(points)
    integer intent(hide), depend(centers) :: num_centers=len(centers)
    integer intent(hide), depend(shell_types) :: num_shells=len(shell_types)
    integer intent(hide), depend(ccoeffs) :: num_ccoeffs=len(ccoeffs)
    integer intent(hide), depend(exponents) :: num_exponents=len(exponents)
  end function gint1_fn_overlap
"""


class Gint1Fn(GaussianIntegral):
    def __init__(self):
        a = ShellArgGroup("a")
        p = PointArgGroup("p")
        self.a = a.args[0].symbols
        self.a_a = a.args[1].symbol
        self.p =  p.args[0].symbols
        name = "gint1_fn"
        interface_fns = [
            (code_gint1_fn_basis_c, code_gint1_fn_basis_pyf),
            (code_gint1_fn_dmat_c, code_gint1_fn_dmat_pyf),
            (code_gint1_fn_overlap_c, code_gint1_fn_overlap_pyf),
        ]
        GaussianIntegral.__init__(self, name, [a, p], interface_fns, [])

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
