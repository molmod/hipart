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
#define MAX_SHELL 3
#define NUM_SHELL_TYPES 7
#define MAX_SHELL_DOF 10

static void fn_pF(double* a, double a_a, double* p, double* out)
{
  double fix0, fix1, fix2, tmp0, tmp1, tmp10, tmp11, tmp12, tmp13, tmp5, tmp6, tmp9;
  fix0 = p[0] - a[0]; // oblige
  fix1 = p[1] - a[1]; // oblige
  fix2 = p[2] - a[2]; // oblige
  tmp11 = fix2*fix2; // auto
  tmp12 = fix1*fix1; // auto
  tmp13 = fix0*fix0; // auto
  tmp0 = exp(-a_a*(tmp11 + tmp12 + tmp13)); // auto
  tmp1 = pow(a_a,2.25); // auto
  tmp12 = tmp1*tmp12; // auto+recycle
  tmp13 = tmp1*tmp13; // auto+recycle
  tmp11 = tmp1*tmp11; // auto+recycle
  tmp5 = fix2*tmp13; // auto
  tmp6 = fix2*tmp12; // auto
  out[0] = tmp0*(-2.20823713394864*tmp5 - 2.20823713394864*tmp6 + 1.47215808929909*fix2*tmp11); // final
  tmp9 = fix0*tmp13; // auto
  tmp10 = fix0*tmp12; // auto
  out[1] = tmp0*(-0.901509034873353*tmp10 - 0.901509034873353*tmp9 + 3.60603613949341*fix0*tmp11); // final
  tmp13 = fix1*tmp13; // auto+recycle
  tmp12 = fix1*tmp12; // auto+recycle
  out[2] = tmp0*(-0.901509034873353*tmp12 - 0.901509034873353*tmp13 + 3.60603613949341*fix1*tmp11); // final
  out[3] = tmp0*(2.85082188141996*tmp5 - 2.85082188141996*tmp6); // final
  out[4] = 5.70164376283992*fix0*fix1*fix2*tmp0*tmp1; // final
  out[5] = tmp0*(1.16384315950667*tmp9 - 3.49152947852002*tmp10); // final
  out[6] = tmp0*(3.49152947852002*tmp13 - 1.16384315950667*tmp12); // final
}

static void fn_pD(double* a, double a_a, double* p, double* out)
{
  double fix0, fix1, fix2, tmp0, tmp1, tmp6, tmp7, tmp8;
  fix0 = p[0] - a[0]; // oblige
  fix1 = p[1] - a[1]; // oblige
  fix2 = p[2] - a[2]; // oblige
  tmp6 = fix2*fix2; // auto
  tmp7 = fix1*fix1; // auto
  tmp8 = fix0*fix0; // auto
  tmp0 = exp(-a_a*(tmp6 + tmp7 + tmp8)); // auto
  tmp1 = pow(a_a,1.75); // auto
  tmp7 = tmp1*tmp7; // auto+recycle
  tmp8 = tmp1*tmp8; // auto+recycle
  out[0] = tmp0*(-0.822961390324745*tmp7 - 0.822961390324745*tmp8 + 1.64592278064949*tmp1*tmp6); // final
  tmp1 = tmp0*tmp1; // auto+recycle
  tmp6 = fix2*tmp1; // auto+recycle
  out[1] = 2.85082188141996*fix0*tmp6; // final
  out[2] = 2.85082188141996*fix1*tmp6; // final
  out[3] = tmp0*(1.42541094070998*tmp8 - 1.42541094070998*tmp7); // final
  out[4] = 2.85082188141996*fix0*fix1*tmp1; // final
}

static void fn_SP(double* a, double a_a, double* p, double* out)
{
  double fix0, fix1, fix2, tmp0;
  fix0 = p[0] - a[0]; // oblige
  fix1 = p[1] - a[1]; // oblige
  fix2 = p[2] - a[2]; // oblige
  tmp0 = exp(-a_a*(fix0*fix0 + fix1*fix1 + fix2*fix2)); // auto
  out[0] = 0.71270547035499*tmp0*pow(a_a,0.75); // final
  tmp0 = tmp0*pow(a_a,1.25); // auto+recycle
  out[1] = 1.42541094070998*fix0*tmp0; // final
  out[2] = 1.42541094070998*fix1*tmp0; // final
  out[3] = 1.42541094070998*fix2*tmp0; // final
}

static void fn_S(double* a, double a_a, double* p, double* out)
{
  double fix0, fix1, fix2;
  fix0 = p[0] - a[0]; // oblige
  fix1 = p[1] - a[1]; // oblige
  fix2 = p[2] - a[2]; // oblige
  out[0] = 0.71270547035499*pow(a_a,0.75)*exp(-a_a*(fix0*fix0 + fix1*fix1 + fix2*fix2)); // final
}

static void fn_P(double* a, double a_a, double* p, double* out)
{
  double fix0, fix1, fix2, tmp0;
  fix0 = p[0] - a[0]; // oblige
  fix1 = p[1] - a[1]; // oblige
  fix2 = p[2] - a[2]; // oblige
  tmp0 = pow(a_a,1.25)*exp(-a_a*(fix0*fix0 + fix1*fix1 + fix2*fix2)); // auto
  out[0] = 1.42541094070998*fix0*tmp0; // final
  out[1] = 1.42541094070998*fix1*tmp0; // final
  out[2] = 1.42541094070998*fix2*tmp0; // final
}

static void fn_cD(double* a, double a_a, double* p, double* out)
{
  double fix0, fix1, fix2, tmp0, tmp2, tmp3;
  fix0 = p[0] - a[0]; // oblige
  fix1 = p[1] - a[1]; // oblige
  fix2 = p[2] - a[2]; // oblige
  tmp2 = fix1*fix1; // auto
  tmp3 = fix0*fix0; // auto
  tmp0 = pow(a_a,1.75)*exp(-a_a*(tmp2 + tmp3 + fix2*fix2)); // auto
  out[0] = 1.64592278064949*tmp0*tmp3; // final
  out[1] = 2.85082188141996*fix0*fix1*tmp0; // final
  tmp3 = fix2*tmp0; // auto+recycle
  out[2] = 2.85082188141996*fix0*tmp3; // final
  out[3] = 1.64592278064949*tmp0*tmp2; // final
  out[4] = 2.85082188141996*fix1*tmp3; // final
  out[5] = 1.64592278064949*fix2*tmp3; // final
}

static void fn_cF(double* a, double a_a, double* p, double* out)
{
  double fix0, fix1, fix2, tmp0, tmp4, tmp5, tmp6;
  fix0 = p[0] - a[0]; // oblige
  fix1 = p[1] - a[1]; // oblige
  fix2 = p[2] - a[2]; // oblige
  tmp4 = fix2*fix2; // auto
  tmp5 = fix1*fix1; // auto
  tmp6 = fix0*fix0; // auto
  tmp0 = pow(a_a,2.25)*exp(-a_a*(tmp4 + tmp5 + tmp6)); // auto
  tmp6 = tmp0*tmp6; // auto+recycle
  out[0] = 1.47215808929909*fix0*tmp6; // final
  out[1] = 3.29184556129898*fix1*tmp6; // final
  out[2] = 3.29184556129898*fix2*tmp6; // final
  tmp5 = tmp0*tmp5; // auto+recycle
  out[3] = 3.29184556129898*fix0*tmp5; // final
  out[4] = 5.70164376283992*fix0*fix1*fix2*tmp0; // final
  tmp4 = tmp0*tmp4; // auto+recycle
  out[5] = 3.29184556129898*fix0*tmp4; // final
  out[6] = 1.47215808929909*fix1*tmp5; // final
  out[7] = 3.29184556129898*fix2*tmp5; // final
  out[8] = 3.29184556129898*fix1*tmp4; // final
  out[9] = 1.47215808929909*fix2*tmp4; // final
}

typedef void (*fntype)(double*, double, double*, double*);
const fntype fns[7] = {fn_pF, fn_pD, fn_SP, fn_S, fn_P, fn_cD, fn_cF};

void gint1_fn(int a_s, double* a, double a_a, double* p, double* out)
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
        gint1_fn(shell_type, center, *exponent, points, work);
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
        gint1_fn(shell_type, center, *exponent, points, work);
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
