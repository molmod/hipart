// HiPart is a program to analyze the electronic structure of molecules with
// fuzzy-atom partitioning methods.
// Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
//--


#include <math.h>
#include <stdlib.h>
#include "gint1_fn.h"
#include "gint1_fn_auto.h" // auto-generated code


int gint1_fns_basis(double* weights, double* orbs, double* points,
  double* centers, int* shell_types, int* shell_map,  int* num_primitives,
  double* ccoeffs, double* exponents, int num_orbs, int num_dof, int num_points,
  int num_centers, int num_shells, int num_ccoeffs, int num_exponents)
{
  int shell, shell_type, primitive, dof, num_shell_dof, result, i_point, i_orb;
  double *work, *out, *center, *ccoeff, *exponent, *weight, *shell_fns, *basis_fns, *fn;
  double tmp;

  result = 0;

  work = malloc(MAX_SHELL_DOF*sizeof(double));
  CHECK_ALLOC(work);
  basis_fns = malloc(num_dof*sizeof(double));
  CHECK_ALLOC(basis_fns);

  for (i_point=0; i_point<num_points; i_point++) {
    // A) clear the basis functions.
    fn = basis_fns;
    for (dof=0; dof<num_dof; dof++) {
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
      num_shell_dof = GET_SHELL_DOF(shell_type);
      for (primitive=0; primitive<num_primitives[shell]; primitive++) {
        gint1_fn_dispatch(shell_type, center, *exponent, points, work);
        //printf("shell_type=%d  primitive=%d  exponent=%f\n", shell_type, primitive, *exponent);
        out = work;
        exponent++;
        fn = shell_fns;
        // Take care of exceptional contraction rules for SP shells
        if (shell_type==-1) {
          *fn += (*out)*(*ccoeff);
          //printf("out=%f  ccoeff=%f  fn=%f\n", *out, *ccoeff, (*out)*(*ccoeff));
          fn++;
          out++;
          ccoeff++;
        }
        for (dof=num_shell_dof-(shell_type==-1); dof>0; dof--) {
          *fn += (*out)*(*ccoeff);
          //printf("out=%f  ccoeff=%f  fn=%f\n", *out, *ccoeff, (*out)*(*ccoeff));
          fn++;
          out++;
        }
        ccoeff++;
      }
      shell_fns += num_shell_dof;
    }
    //printf("\n");
    // C) Make proper linear combinations with weights and assign to orb
    for (i_orb=0; i_orb<num_orbs; i_orb++) {
      tmp = 0.0;
      fn = basis_fns;
      weight = weights+i_orb*num_dof;
      for (dof=0; dof<num_dof; dof++) {
        tmp += (*fn)*(*weight);
        fn++;
        weight++;
      }
      *(orbs+i_orb*num_points) = tmp;
    }
    points += 3;
    orbs++;
    //printf("\n");
  }

EXIT:
  free(work);
  free(basis_fns);
  return result;
}


int gint1_fn_dmat(double* dmat, double* density, double* points,
  double* centers, int* shell_types, int* shell_map,  int* num_primitives,
  double* ccoeffs, double* exponents, int num_dmat, int num_points,
  int num_centers, int num_shells, int num_ccoeffs, int num_exponents)
{
  int shell, shell_type, primitive, dof1, dof2, num_dof, num_shell_dof, result, i_point;
  double *work, *out, *center, *ccoeff, *exponent, *basis_fns, *shell_fns, *fn1, *fn2, *dmat_element;
  double tmp;

  result = 0;

  work = malloc(MAX_SHELL_DOF*sizeof(double));
  CHECK_ALLOC(work);
  num_dof = ((int)(sqrt(1.0+8.0*num_dmat)-1.0))/2;
  basis_fns = malloc(num_dof*sizeof(double));
  CHECK_ALLOC(basis_fns);

  for (i_point=0; i_point<num_points; i_point++) {
    // A) clear the basis functions.
    fn1 = basis_fns;
    for (dof1=0; dof1<num_dof; dof1++) {
      *fn1 = 0.0;
      fn1++;
    }
    // B) evaluate the basis functions in the current point.
    ccoeff = ccoeffs;
    exponent = exponents;
    shell_fns = basis_fns;
    for (shell=0; shell<num_shells; shell++) {
      center = centers + (3*shell_map[shell]);
      shell_type = shell_types[shell];
      CHECK_SHELL(shell_type);
      num_shell_dof = GET_SHELL_DOF(shell_type);
      for (primitive=0; primitive<num_primitives[shell]; primitive++) {
        gint1_fn_dispatch(shell_type, center, *exponent, points, work);
        //printf("shell_type=%d  primitive=%d  exponent=%f\n", shell_type, primitive, *exponent);
        out = work;
        exponent++;
        fn1 = shell_fns;
        // Take care of exceptional contraction rules for SP shells
        if (shell_type==-1) {
          *fn1 += (*out)*(*ccoeff);
          //printf("out=%f  ccoeff=%f  fn=%f\n", *out, *ccoeff, (*out)*(*ccoeff));
          fn1++;
          out++;
          ccoeff++;
        }
        for (dof1=num_shell_dof-(shell_type==-1); dof1>0; dof1--) {
          *fn1 += (*out)*(*ccoeff);
          //printf("out=%f  ccoeff=%f  fn=%f\n", *out, *ccoeff, (*out)*(*ccoeff));
          fn1++;
          out++;
        }
        ccoeff++;
      }
      shell_fns += num_shell_dof;
    }
    //printf("\n");
    // C) Make dot product of basis functions with density matrix.
    *density = 0.0;
    dmat_element = dmat;
    fn1 = basis_fns;
    for (dof1=0; dof1<num_dof; dof1++) {
      fn2 = basis_fns;
      tmp = 0.0;
      // off-diagonal times 2
      for (dof2=0; dof2<dof1; dof2++) {
        tmp += (*fn2)*(*dmat_element);
        //printf("dof1=%d  dof2=%d  basis1=%f  basis2=%f  dmat_element=%f  density=%f\n", dof1, dof2, basis_fns[dof1], basis_fns[dof2], *dmat_element, *density);
        dmat_element++;
        fn2++;
      }
      tmp *= 2;
      // diagonal
      tmp += (*fn1)*(*dmat_element);
      //printf("dof1=%d  dof2=%d  basis1=%f  basis2=%f  dmat_element=%f  density=%f\n", dof1, dof2, basis_fns[dof1], basis_fns[dof2], *dmat_element, *density);
      dmat_element++;
      *density += tmp*(*fn1);
      fn1++;
    }
    // D) Prepare for next iteration
    density++;
    points += 3;
    //printf("\n");
  }

EXIT:
  free(work);
  free(basis_fns);
  return result;
}


int gint1_fn_overlap(double* overlap, double* points, double* weights,
  double* centers, int* shell_types, int* shell_map, int* num_primitives,
  double* ccoeffs, double* exponents, int num_overlap, int num_points,
  int num_centers, int num_shells, int num_ccoeffs, int num_exponents)
{
  int shell, shell_type, primitive, dof1, dof2, num_dof, num_shell_dof, result, i_point;
  double *work, *out, *center, *ccoeff, *exponent, *basis_fns, *shell_fns, *fn1, *fn2, *overlap_element;
  double tmp;

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
    fn1 = basis_fns;
    for (dof1=0; dof1<num_dof; dof1++) {
      *fn1 = 0.0;
      fn1++;
    }
    // B) evaluate the basis functions in the current point.
    ccoeff = ccoeffs;
    exponent = exponents;
    shell_fns = basis_fns;
    for (shell=0; shell<num_shells; shell++) {
      center = centers + (3*shell_map[shell]);
      shell_type = shell_types[shell];
      CHECK_SHELL(shell_type);
      num_shell_dof = GET_SHELL_DOF(shell_type);
      for (primitive=0; primitive<num_primitives[shell]; primitive++) {
        gint1_fn_dispatch(shell_type, center, *exponent, points, work);
        //printf("shell_type=%d  primitive=%d  exponent=%f\n", shell_type, primitive, *exponent);
        out = work;
        exponent++;
        fn1 = shell_fns;
        // Take care of exceptional contraction rules for SP shells
        if (shell_type==-1) {
          *fn1 += (*out)*(*ccoeff);
          //printf("out=%f  ccoeff=%f  fn=%f\n", *out, *ccoeff, (*out)*(*ccoeff));
          fn1++;
          out++;
          ccoeff++;
        }
        for (dof1=num_shell_dof-(shell_type==-1); dof1>0; dof1--) {
          *fn1 += (*out)*(*ccoeff);
          //printf("out=%f  ccoeff=%f  fn=%f\n", *out, *ccoeff, (*out)*(*ccoeff));
          fn1++;
          out++;
        }
        ccoeff++;
      }
      shell_fns += num_shell_dof;
    }
    //printf("\n");
    // C) Multiply overlap elements with the wieght and add to the atomic
    //    overlap matrix
    overlap_element = overlap;
    fn1 = basis_fns;
    for (dof1=0; dof1<num_dof; dof1++) {
      fn2 = basis_fns;
      tmp = (*weights)*(*fn1);
      for (dof2=0; dof2<num_dof; dof2++) {
        *overlap_element += tmp*(*fn2);
        //printf("dof1=%d  dof2=%d  basis1=%f  basis2=%f  overlap_element=%f  weight=%f\n", dof1, dof2, basis_fns[dof1], basis_fns[dof2], *overlap_element, *weight);
        overlap_element++;
        fn2++;
      }
      fn1++;
    }
    // D) Prepare for next iteration
    points += 3;
    weights++;
    //printf("\n");
  }

EXIT:
  free(work);
  free(basis_fns);
  return result;
}
