#!/usr/bin/env python


from writer import *

# TODO: take into account normalization constants
# TODO: reorder in density tests

code_gint2_nai_dmat_c = """\
static int get_shell_dof(int shell_type) {
  if (shell_type==-1) {
    return 4;
  } else if (shell_type > 0) {
    return ((shell_type+1)*(shell_type+2))/2;
  } else {
    return -2*shell_type+1;
  }
}

int gint2_nai_dmat(double* dmat, double* potentials, double* points,
  double* centers, int* shell_types, int* shell_map,  int* num_primitives,
  double* ccoeffs, double* exponents, int num_dmat, int num_points,
  int num_centers, int num_shells, int num_ccoeffs, int num_exponents)
{
  int result, size, i_point, k, k1, k2, primitive1, primitive2, dof1, dof2;
  int shell1, shell2, shell_type1, shell_type2, shell_dof1, shell_dof2, shell_offset1, shell_offset2;
  double *work, *work_sum, *out, *out_sum, *center1, *center2;
  double *shell_ccoeffs1, *shell_ccoeffs2, *ccoeff1, *ccoeff2, c1, c2;
  double *shell_exponents1, *shell_exponents2, *exponent1, *exponent2;

  result = 0;

  size = MAX_SHELL_DOF*MAX_SHELL_DOF;
  work = malloc(size*sizeof(double));
  if (work==NULL) {result = -1; goto EXIT;}
  work_sum = malloc(size*sizeof(double));
  if (work_sum==NULL) {result = -1; goto EXIT;}

  for (i_point=0; i_point<num_points; i_point++) {
    // A) clear the result
    *potentials = 0.0;
    // B) evaluate the potential in the current point.
    // prep inner loop.
    shell_ccoeffs1 = ccoeffs;
    shell_exponents1 = exponents;
    shell_offset1 = 0;
    for (shell1=0; shell1<num_shells; shell1++) {
      center1 = centers + (3*shell_map[shell1]);
      shell_type1 = shell_types[shell1];
      shell_dof1 = get_shell_dof(shell_type1);
      //printf("shell1=%d  type=%d  dof=%d  offset=%d\\n", shell1, shell_type1, shell_dof1, shell_offset1);
      // prep inner loop.
      shell_ccoeffs2 = ccoeffs;
      shell_exponents2 = exponents;
      shell_offset2 = 0;
      for (shell2=0; shell2<num_shells; shell2++) {
        center2 = centers + (3*shell_map[shell2]);
        shell_type2 = shell_types[shell2];
        shell_dof2 = get_shell_dof(shell_type2);
        //printf("  shell2=%d  type=%d  dof=%d  offset=%d\\n", shell2, shell_type2, shell_dof2, shell_offset2);
        // Clear the worksum
        out_sum = work_sum;
        for (k=0; k<size; k++) { *out_sum = 0.0; out_sum++; }
        // Build up the worksum by combining the expectation values of the
        // primitives into the worksum. (density matrix is not yet included).
        // prep inner loop.
        exponent1 = shell_exponents1;
        ccoeff1 = shell_ccoeffs1;
        for (primitive1=0; primitive1<num_primitives[shell1]; primitive1++) {
          // prep inner loop.
          //printf("    primitive1=%d  exponent1=%f\\n", primitive1, *exponent1);
          exponent2 = shell_exponents2;
          ccoeff2 = shell_ccoeffs2;
          for (primitive2=0; primitive2<num_primitives[shell2]; primitive2++) {
            //printf("      primitive2=%d  exponent2=%f\\n", primitive2, *exponent2);
            gint2_nai_dispatch(shell_type1, center1, *exponent1, shell_type2, center2, *exponent2, points, work);
            // add to work sum
            out = work;
            out_sum = work_sum;
            for (dof1=0; dof1<shell_dof1; dof1++) {
              c1 = ccoeff1[(shell_type1==-1)&&(dof1>0)];
              for (dof2=0; dof2<shell_dof2; dof2++) {
                c2 = ccoeff2[(shell_type2==-1)&&(dof2>0)];
                *out_sum += c1*c2*(*out);
                //printf("        dof1=%d  dof2=%d  c1=%f  c2=%f  out=%f  out_sum=%f\\n", dof1, dof2, c1, c2, *out, *out_sum);
                out++;
                out_sum++;
              }
            }
            exponent2++;
            ccoeff2 += 1+(shell_type2==-1);
          }
          exponent1++;
          ccoeff1 += 1+(shell_type1==-1);
        }
        // Compute the trace of the product of the work sum and part of the
        // density matrix that belongs to the (shell1,shell2) part
        out_sum = work_sum;
        for (dof1=0; dof1<shell_dof1; dof1++) {
          for (dof2=0; dof2<shell_dof2; dof2++) {
            k1 = dof1 + shell_offset1;
            k2 = dof2 + shell_offset2;
            if (k1>k2) {
              k = (k1*(k1+1))/2+k2;
            } else {
              k = (k2*(k2+1))/2+k1;
            }
            *potentials += (*out_sum)*dmat[k];
            //printf("++++dof1=%d  dof2=%d  k1=%d  k2=%d  k=%d  out_sum=%f  dmat[k]=%f  potential=%f\\n", dof1, dof2, k1, k2, k, *out_sum, dmat[k], *potentials);
            out_sum++;
          }
        }
        // Move on.
        shell_exponents2 += num_primitives[shell2];
        shell_ccoeffs2 += num_primitives[shell2]*(1+(shell_type2==-1));
        shell_offset2 += shell_dof2;
      }
      // Move on.
      shell_exponents1 += num_primitives[shell1];
      shell_ccoeffs1 += num_primitives[shell1]*(1+(shell_type1==-1));
      shell_offset1 += shell_dof1;
    }
    // C) Move on
    points += 3;
    potentials++;
    //printf("\\n");
  }

EXIT:
  free(work);
  free(work_sum);
  return result;
}"""

code_gint2_nai_dmat_pyf = """\
  integer function gint2_nai_dmat(dmat, potentials, points, centers, shell_types, shell_map, num_primitives, ccoeffs, exponents, num_dmat, num_points, num_centers, num_shells, num_ccoeffs, num_exponents)
    intent(c) gint2_nai_dmat
    intent(c)
    double precision intent(in) :: dmat(num_dmat)
    double precision intent(inout) :: potentials(num_points)
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
  end function gint2_nai_dmat
"""


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
        interface_fns = [(code_gint2_nai_dmat_c, code_gint2_nai_dmat_pyf)]
        includes = ["gaux.h"]
        GaussianIntegral.__init__(self, name, [a, b, c], interface_fns, includes)

    def get_SS_integral(self, ab_a, ab_overlap, u2, m):
        return 2*sqrt(ab_a/pi)*ab_overlap*C.Function("gaux")(u2, m)

    def get_integral(self, a_n, b_n, u, v, w, ab_a, ab_overlap, u2, m):
        # This the recurrence relation from the paper of Obara and Saika
        # (see http://dx.doi.org/10.1063/1.450106)
        if (a_n[0] < 0) or (a_n[1] < 0) or (a_n[2] < 0) or (b_n[0] < 0) or (b_n[1] < 0) or (b_n[2] < 0):
            # This is where the recurrence relations leads to zero terms
            return 0
        elif (a_n[0] == 0) and (a_n[1] == 0) and (a_n[2] == 0) and (b_n[0] == 0) and (b_n[1] == 0) and (b_n[2] == 0):
            # This is where the recurrence relations leads to the SS term
            return self.get_SS_integral(ab_a, ab_overlap, u2, m)
        else:
            # apply the recurrence relation
            # Note: we must decide which of the six possible paths we take.
            # Current strategy: get rid of the largest power first
            if max(a_n) > max(b_n):
                index = [i for i, n in enumerate(a_n) if n==max(a_n)][0]
                a_n1 = list(a_n)
                a_n1[index] -= 1
                a_n1 = tuple(a_n1)
                b_n1 = list(b_n)
                b_n1[index] -= 1
                b_n1 = tuple(b_n1)
                a_n2 = list(a_n)
                a_n2[index] -= 2
                a_n2 = tuple(a_n2)
                return (
                      v[index]*self.get_integral(a_n1, b_n, u, v, w, ab_a, ab_overlap, u2, m)
                    - u[index]*self.get_integral(a_n1, b_n, u, v, w, ab_a, ab_overlap, u2, m+1)
                ) + a_n1[index]/(2*ab_a)*(
                      self.get_integral(a_n2, b_n, u, v, w, ab_a, ab_overlap, u2, m)
                    - self.get_integral(a_n2, b_n, u, v, w, ab_a, ab_overlap, u2, m+1)
                ) + b_n[index]/(2*ab_a)*(
                      self.get_integral(a_n1, b_n1, u, v, w, ab_a, ab_overlap, u2, m)
                    - self.get_integral(a_n1, b_n1, u, v, w, ab_a, ab_overlap, u2, m+1)
                )
            else:
                index = [i for i, n in enumerate(b_n) if n==max(b_n)][0]
                b_n1 = list(b_n)
                b_n1[index] -= 1
                b_n1 = tuple(b_n1)
                a_n1 = list(a_n)
                a_n1[index] -= 1
                a_n1 = tuple(a_n1)
                b_n2 = list(b_n)
                b_n2[index] -= 2
                b_n2 = tuple(b_n2)
                return (
                      w[index]*self.get_integral(a_n, b_n1, u, v, w, ab_a, ab_overlap, u2, m)
                    - u[index]*self.get_integral(a_n, b_n1, u, v, w, ab_a, ab_overlap, u2, m+1)
                ) + b_n1[index]/(2*ab_a)*(
                      self.get_integral(a_n, b_n2, u, v, w, ab_a, ab_overlap, u2, m)
                    - self.get_integral(a_n, b_n2, u, v, w, ab_a, ab_overlap, u2, m+1)
                ) + a_n[index]/(2*ab_a)*(
                      self.get_integral(a_n1, b_n1, u, v, w, ab_a, ab_overlap, u2, m)
                    - self.get_integral(a_n1, b_n1, u, v, w, ab_a, ab_overlap, u2, m+1)
                )

    def get_expressions(self, st_row):
        ab_a = self.a_a + self.b_a
        ab_b = self.a_a*self.b_a/ab_a
        p = ((self.a_a*self.a[0] + self.b_a*self.b[0])/ab_a,
             (self.a_a*self.a[1] + self.b_a*self.b[1])/ab_a,
             (self.a_a*self.a[2] + self.b_a*self.b[2])/ab_a,)
        d = (self.a[0] - self.b[0],
             self.a[1] - self.b[1],
             self.a[2] - self.b[2],)
        d2 = (d[0])**2 + (d[1])**2 + (d[2])**2
        u = (p[0] - self.c[0],
             p[1] - self.c[1],
             p[2] - self.c[2],)
        v = (p[0] - self.a[0],
             p[1] - self.a[1],
             p[2] - self.a[2],)
        w = (p[0] - self.b[0],
             p[1] - self.b[1],
             p[2] - self.b[2],)
        u2 = ab_a*((u[0])**2 + (u[1])**2 + (u[2])**2)
        ab_overlap = sqrt(pi/ab_a)**3*C.Function("exp")(-ab_b*d2)
        obliges = [d[0], d[1], d[2], u[0], u[1], u[2], v[0], v[1], v[2], w[0], w[1], w[2]]
        if st_row == (0,0,0):
            nai = self.get_SS_integral(ab_a, ab_overlap, u2, 0)
            nai /= get_cartesian_wfn_norm(self.a_a, 0, 0, 0).evalf()
            nai /= get_cartesian_wfn_norm(self.b_a, 0, 0, 0).evalf()
            results = [nai]
        elif st_row[0] >= -1 and st_row[1] >= -1:
            results = []
            for a_n in iter_cartesian_powers(st_row[0]):
                for b_n in iter_cartesian_powers(st_row[1]):
                    nai = self.get_integral(a_n, b_n, u, v, w, ab_a, ab_overlap, u2, 0)
                    nai /= get_cartesian_wfn_norm(self.a_a, a_n[0], a_n[1], a_n[2]).evalf()
                    nai /= get_cartesian_wfn_norm(self.b_a, b_n[0], b_n[1], b_n[2]).evalf()
                    results.append(nai)
        else:
            results = []
            obliges = []
        key = "_".join(get_shell_label(st) for st in st_row[:-1])
        return key, results, obliges


def main():
    write(Gint2NAI(), max_shell=1)


if __name__ == "__main__":
    main()
