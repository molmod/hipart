// HiPart is a software toolkit to analyse molecular densities with the hirshfeld partitioning scheme.
// Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>
//
// This file is part of HiPart.
//
// HiPart is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// HiPart is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
// --


#include "gint2_nai.h"
#include "gint2_nai.inc.c"


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
  CHECK_ALLOC(work);
  work_sum = malloc(size*sizeof(double));
  CHECK_ALLOC(work_sum);

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
      CHECK_SHELL(shell_type1);
      shell_dof1 = GET_SHELL_DOF(shell_type1);
      //printf("shell1=%d  type=%d  dof=%d  offset=%d\n", shell1, shell_type1, shell_dof1, shell_offset1);
      // prep inner loop.
      shell_ccoeffs2 = ccoeffs;
      shell_exponents2 = exponents;
      shell_offset2 = 0;
      for (shell2=0; shell2<num_shells; shell2++) {
        center2 = centers + (3*shell_map[shell2]);
        shell_type2 = shell_types[shell2];
        CHECK_SHELL(shell_type2);
        shell_dof2 = GET_SHELL_DOF(shell_type2);
        //printf("  shell2=%d  type=%d  dof=%d  offset=%d\n", shell2, shell_type2, shell_dof2, shell_offset2);
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
          //printf("    primitive1=%d  exponent1=%f\n", primitive1, *exponent1);
          exponent2 = shell_exponents2;
          ccoeff2 = shell_ccoeffs2;
          for (primitive2=0; primitive2<num_primitives[shell2]; primitive2++) {
            //printf("      primitive2=%d  exponent2=%f\n", primitive2, *exponent2);
            gint2_nai_dispatch(shell_type1, center1, *exponent1, shell_type2, center2, *exponent2, points, work);
            // add to work sum
            out = work;
            out_sum = work_sum;
            for (dof1=0; dof1<shell_dof1; dof1++) {
              c1 = ccoeff1[(shell_type1==-1)&&(dof1>0)];
              for (dof2=0; dof2<shell_dof2; dof2++) {
                c2 = ccoeff2[(shell_type2==-1)&&(dof2>0)];
                *out_sum += c1*c2*(*out);
                //printf("        dof1=%d  dof2=%d  c1=%f  c2=%f  out=%f  out_sum=%f\n", dof1, dof2, c1, c2, *out, *out_sum);
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
            //printf("++++dof1=%d  dof2=%d  k1=%d  k2=%d  k=%d  out_sum=%f  dmat[k]=%f  potential=%f\n", dof1, dof2, k1, k2, k, *out_sum, dmat[k], *potentials);
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
    //printf("\n");
  }

EXIT:
  free(work);
  free(work_sum);
  return result;
}
