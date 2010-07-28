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



#include <math.h>
#include <stdlib.h>
#include "gint1_fn.h"
#define MAX_SHELL 3
#define NUM_SHELL_TYPES 7
#define MAX_SHELL_DOF 10

void gint1_fn_pF(double* a, double a_a, double* p, double* out)
{
  // Number of local variables: 9
  double tmp11, tmp12, tmp13, tmp2, tmp3, tmp4, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp13 = v_1*v_1; // auto, weighs 1
  tmp2 = tmp13; // auto, weighs 0
  tmp3 = v_0*v_0; // auto, weighs 1
  tmp4 = v_2*v_2; // auto, weighs 1
  tmp12 = pow(a_a,2.25)*exp(-a_a*(tmp2 + tmp3 + tmp4)); // auto, weighs 7
  tmp3 = tmp12*tmp3; // auto+recycle, weighs 1
  tmp2 = tmp12*tmp2; // auto+recycle, weighs 1
  out[0] = v_2*(-2.20823713394864*tmp2 - 2.20823713394864*tmp3) + 1.47215808929909*tmp12*pow(v_2,3); // final, weighs 8
  tmp4 = 3.60603613949341*tmp12*tmp4; // auto+recycle, weighs 2
  tmp11 = tmp12*pow(v_0,3); // auto, weighs 2
  out[1] = -0.901509034873353*tmp11 + v_0*(tmp4 - 0.901509034873353*tmp2); // final, weighs 5
  tmp12 = tmp12*v_1; // auto+recycle, weighs 1
  tmp13 = tmp12*tmp13; // auto+recycle, weighs 1
  out[2] = -0.901509034873353*tmp13 + v_1*(tmp4 - 0.901509034873353*tmp3); // final, weighs 5
  out[3] = v_2*(2.85082188141996*tmp3 - 2.85082188141996*tmp2); // final, weighs 4
  out[4] = 5.70164376283992*tmp12*v_0*v_2; // final, weighs 3
  out[5] = 1.16384315950667*tmp11 - 3.49152947852002*tmp2*v_0; // final, weighs 4
  out[6] = -1.16384315950667*tmp13 + 3.49152947852002*tmp3*v_1; // final, weighs 4
  // total weight = 57
}

void gint1_fn_pD(double* a, double a_a, double* p, double* out)
{
  // Number of local variables: 7
  double tmp10, tmp2, tmp3, tmp4, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp2 = v_1*v_1; // auto, weighs 1
  tmp3 = v_0*v_0; // auto, weighs 1
  tmp4 = v_2*v_2; // auto, weighs 1
  tmp10 = pow(a_a,1.75)*exp(-a_a*(tmp2 + tmp3 + tmp4)); // auto, weighs 7
  tmp3 = tmp10*tmp3; // auto+recycle, weighs 1
  tmp2 = tmp10*tmp2; // auto+recycle, weighs 1
  out[0] = -0.822961390324745*tmp2 - 0.822961390324745*tmp3 + 1.64592278064949*tmp10*tmp4; // final, weighs 6
  tmp4 = tmp10*v_2; // auto+recycle, weighs 1
  out[1] = 2.85082188141996*tmp4*v_0; // final, weighs 2
  out[2] = 2.85082188141996*tmp4*v_1; // final, weighs 2
  out[3] = 1.42541094070998*tmp3 - 1.42541094070998*tmp2; // final, weighs 3
  out[4] = 2.85082188141996*tmp10*v_0*v_1; // final, weighs 3
  // total weight = 35
}

void gint1_fn_SP(double* a, double a_a, double* p, double* out)
{
  // Number of local variables: 4
  double tmp0, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp0 = exp(-a_a*(v_0*v_0 + v_1*v_1 + v_2*v_2)); // auto, weighs 8
  out[0] = 0.71270547035499*tmp0*pow(a_a,0.75); // final, weighs 3
  tmp0 = tmp0*pow(a_a,1.25); // auto+recycle, weighs 2
  out[1] = 1.42541094070998*tmp0*v_0; // final, weighs 2
  out[2] = 1.42541094070998*tmp0*v_1; // final, weighs 2
  out[3] = 1.42541094070998*tmp0*v_2; // final, weighs 2
  // total weight = 25
}

void gint1_fn_S(double* a, double a_a, double* p, double* out)
{
  // Number of local variables: 3
  double v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  out[0] = 0.71270547035499*pow(a_a,0.75)*exp(-a_a*(v_0*v_0 + v_1*v_1 + v_2*v_2)); // final, weighs 11
  // total weight = 17
}

void gint1_fn_P(double* a, double a_a, double* p, double* out)
{
  // Number of local variables: 4
  double tmp2, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp2 = pow(a_a,1.25)*exp(-a_a*(v_0*v_0 + v_1*v_1 + v_2*v_2)); // auto, weighs 10
  out[0] = 1.42541094070998*tmp2*v_0; // final, weighs 2
  out[1] = 1.42541094070998*tmp2*v_1; // final, weighs 2
  out[2] = 1.42541094070998*tmp2*v_2; // final, weighs 2
  // total weight = 22
}

void gint1_fn_cD(double* a, double a_a, double* p, double* out)
{
  // Number of local variables: 7
  double tmp2, tmp3, tmp4, tmp8, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp2 = v_2*v_2; // auto, weighs 1
  tmp3 = v_1*v_1; // auto, weighs 1
  tmp4 = v_0*v_0; // auto, weighs 1
  tmp8 = pow(a_a,1.75)*exp(-a_a*(tmp2 + tmp3 + tmp4)); // auto, weighs 7
  out[0] = 1.64592278064949*tmp4*tmp8; // final, weighs 2
  out[1] = 2.85082188141996*tmp8*v_0*v_1; // final, weighs 3
  tmp4 = tmp8*v_2; // auto+recycle, weighs 1
  out[2] = 2.85082188141996*tmp4*v_0; // final, weighs 2
  out[3] = 1.64592278064949*tmp3*tmp8; // final, weighs 2
  out[4] = 2.85082188141996*tmp4*v_1; // final, weighs 2
  out[5] = 1.64592278064949*tmp2*tmp8; // final, weighs 2
  // total weight = 30
}

void gint1_fn_cF(double* a, double a_a, double* p, double* out)
{
  // Number of local variables: 11
  double tmp11, tmp12, tmp13, tmp14, tmp2, tmp3, tmp4, tmp7, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp12 = v_2*v_2; // auto, weighs 1
  tmp2 = tmp12; // auto, weighs 0
  tmp13 = v_1*v_1; // auto, weighs 1
  tmp3 = tmp13; // auto, weighs 0
  tmp14 = v_0*v_0; // auto, weighs 1
  tmp4 = tmp14; // auto, weighs 0
  tmp11 = pow(a_a,2.25)*exp(-a_a*(tmp2 + tmp3 + tmp4)); // auto, weighs 7
  tmp7 = tmp11*v_0; // auto, weighs 1
  out[0] = 1.47215808929909*tmp14*tmp7; // final, weighs 2
  tmp14 = tmp11*v_1; // auto+recycle, weighs 1
  out[1] = 3.29184556129898*tmp14*tmp4; // final, weighs 2
  tmp11 = tmp11*v_2; // auto+recycle, weighs 1
  out[2] = 3.29184556129898*tmp11*tmp4; // final, weighs 2
  out[3] = 3.29184556129898*tmp3*tmp7; // final, weighs 2
  out[4] = 5.70164376283992*tmp11*v_0*v_1; // final, weighs 3
  out[5] = 3.29184556129898*tmp2*tmp7; // final, weighs 2
  out[6] = 1.47215808929909*tmp13*tmp14; // final, weighs 2
  out[7] = 3.29184556129898*tmp11*tmp3; // final, weighs 2
  out[8] = 3.29184556129898*tmp14*tmp2; // final, weighs 2
  out[9] = 1.47215808929909*tmp11*tmp12; // final, weighs 2
  // total weight = 40
}

typedef void (*fntype)(double*, double, double*, double*);
static const fntype fns[7] = {gint1_fn_pF, gint1_fn_pD, gint1_fn_SP, gint1_fn_S, gint1_fn_P, gint1_fn_cD, gint1_fn_cF};

static void gint1_fn_dispatch(int a_s, double* a, double a_a, double* p, double* out)
{
  fns[3+a_s](a, a_a, p, out);
}

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
        gint1_fn_dispatch(shell_type, center, *exponent, points, work);
        //printf("shell_type=%d  primitive=%d  exponent=%f\n", shell_type, primitive, *exponent);
        out = work;
        exponent++;
        weight = shell_weights;
        if (shell_type==-1) {
          *fns += (*weight)*(*out)*(*ccoeff);
          //printf("weight=%f  out=%f  ccoeff=%f  contrib=%f  fn=%f\n", *weight, *out, *ccoeff, (*weight)*(*out)*(*ccoeff), *fns);
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
          //printf("weight=%f  out=%f  ccoeff=%f  contrib=%f  fn=%f\n", *weight, *out, *ccoeff, (*weight)*(*out)*(*ccoeff), *fns);
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
    //printf("\n");
  }

EXIT:
  free(work);
  return result;
}

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
        gint1_fn_dispatch(shell_type, center, *exponent, points, work);
        //printf("shell_type=%d  primitive=%d  exponent=%f\n", shell_type, primitive, *exponent);
        out = work;
        exponent++;
        fn = shell_fns;
        if (shell_type==-1) {
          *fn += (*out)*(*ccoeff);
          //printf("out=%f  ccoeff=%f  fn=%f\n", *out, *ccoeff, (*out)*(*ccoeff));
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
          //printf("out=%f  ccoeff=%f  fn=%f\n", *out, *ccoeff, (*out)*(*ccoeff));
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
    //printf("\n");
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
        //printf("dof1=%d  dof2=%d  basis1=%f  basis2=%f  dmat_element=%f  density=%f\n", dof1, dof2, basis_fns[dof1], basis_fns[dof2], *dmat_element, *density);
        dmat_element++;
      }
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
