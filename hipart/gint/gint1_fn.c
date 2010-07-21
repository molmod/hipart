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
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11;
  tmp0 = p[2] - a[2];
  tmp1 = pow(a_a,2.25);
  tmp2 = p[1] - a[1];
  tmp3 = tmp2*tmp2;
  tmp4 = p[0] - a[0];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp0*tmp0;
  tmp7 = tmp3 + tmp5 + tmp6;
  tmp8 = -a_a*tmp7;
  tmp9 = exp(tmp8);
  tmp10 = pow(tmp5,(3.0/2.0));
  tmp11 = pow(tmp3,(3.0/2.0));
  out[0] = tmp9*(1.47215808929909*tmp0*tmp1*tmp6 - 2.20823713394864*tmp0*tmp1*tmp3 - 2.20823713394864*tmp0*tmp1*tmp5);
  out[1] = tmp9*(3.60603613949341*tmp1*tmp4*tmp6 - 0.901509034873353*tmp1*tmp3*tmp4 - 0.901509034873353*tmp1*tmp4*tmp5);
  out[2] = tmp9*(3.60603613949341*tmp1*tmp2*tmp6 - 0.901509034873353*tmp1*tmp2*tmp3 - 0.901509034873353*tmp1*tmp2*tmp5);
  out[3] = tmp9*(2.85082188141996*tmp0*tmp1*tmp5 - 2.85082188141996*tmp0*tmp1*tmp3);
  out[4] = 5.70164376283992*tmp0*tmp1*tmp2*tmp4*tmp9;
  out[5] = tmp9*(1.16384315950667*tmp1*tmp4*tmp5 - 3.49152947852002*tmp1*tmp3*tmp4);
  out[6] = tmp9*(3.49152947852002*tmp1*tmp2*tmp5 - 1.16384315950667*tmp1*tmp2*tmp3);
}

static void fn_pD(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  tmp0 = p[2] - a[2];
  tmp1 = tmp0*tmp0;
  tmp2 = pow(a_a,1.75);
  tmp3 = p[1] - a[1];
  tmp4 = tmp3*tmp3;
  tmp5 = p[0] - a[0];
  tmp6 = tmp5*tmp5;
  tmp7 = tmp1 + tmp4 + tmp6;
  tmp8 = -a_a*tmp7;
  tmp9 = exp(tmp8);
  out[0] = tmp9*(1.64592278064949*tmp1*tmp2 - 0.822961390324745*tmp2*tmp4 - 0.822961390324745*tmp2*tmp6);
  out[1] = 2.85082188141996*tmp0*tmp2*tmp5*tmp9;
  out[2] = 2.85082188141996*tmp0*tmp2*tmp3*tmp9;
  out[3] = tmp9*(1.42541094070998*tmp2*tmp6 - 1.42541094070998*tmp2*tmp4);
  out[4] = 2.85082188141996*tmp2*tmp3*tmp5*tmp9;
}

static void fn_SP(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  tmp0 = p[0] - a[0];
  tmp1 = p[1] - a[1];
  tmp2 = tmp1*tmp1;
  tmp3 = tmp0*tmp0;
  tmp4 = p[2] - a[2];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp2 + tmp3 + tmp5;
  tmp7 = -a_a*tmp6;
  tmp8 = exp(tmp7);
  tmp9 = pow(a_a,1.25);
  out[0] = 0.71270547035499*tmp8*pow(tmp9,0.6);
  out[1] = 1.42541094070998*tmp0*tmp8*tmp9;
  out[2] = 1.42541094070998*tmp1*tmp8*tmp9;
  out[3] = 1.42541094070998*tmp4*tmp8*tmp9;
}

static void fn_S(double* a, double a_a, double* p, double* out)
{
  out[0] = 0.71270547035499*pow(a_a,0.75)*exp(-a_a*(pow((p[0] - a[0]),2) + pow((p[1] - a[1]),2) + pow((p[2] - a[2]),2)));
}

static void fn_P(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  tmp0 = p[0] - a[0];
  tmp1 = pow(a_a,1.25);
  tmp2 = p[1] - a[1];
  tmp3 = tmp2*tmp2;
  tmp4 = tmp0*tmp0;
  tmp5 = p[2] - a[2];
  tmp6 = tmp5*tmp5;
  tmp7 = tmp3 + tmp4 + tmp6;
  tmp8 = -a_a*tmp7;
  tmp9 = exp(tmp8);
  out[0] = 1.42541094070998*tmp0*tmp1*tmp9;
  out[1] = 1.42541094070998*tmp1*tmp2*tmp9;
  out[2] = 1.42541094070998*tmp1*tmp5*tmp9;
}

static void fn_cD(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  tmp0 = p[0] - a[0];
  tmp1 = tmp0*tmp0;
  tmp2 = pow(a_a,1.75);
  tmp3 = p[1] - a[1];
  tmp4 = tmp3*tmp3;
  tmp5 = p[2] - a[2];
  tmp6 = tmp5*tmp5;
  tmp7 = tmp1 + tmp4 + tmp6;
  tmp8 = -a_a*tmp7;
  tmp9 = exp(tmp8);
  out[0] = 1.64592278064949*tmp1*tmp2*tmp9;
  out[1] = 2.85082188141996*tmp0*tmp2*tmp3*tmp9;
  out[2] = 2.85082188141996*tmp0*tmp2*tmp5*tmp9;
  out[3] = 1.64592278064949*tmp2*tmp4*tmp9;
  out[4] = 2.85082188141996*tmp2*tmp3*tmp5*tmp9;
  out[5] = 1.64592278064949*tmp2*tmp6*tmp9;
}

static void fn_cF(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  tmp0 = p[0] - a[0];
  tmp1 = pow(a_a,2.25);
  tmp2 = p[1] - a[1];
  tmp3 = tmp0*tmp0;
  tmp4 = tmp2*tmp2;
  tmp5 = p[2] - a[2];
  tmp6 = tmp5*tmp5;
  tmp7 = tmp3 + tmp4 + tmp6;
  tmp8 = -a_a*tmp7;
  tmp9 = exp(tmp8);
  out[0] = 1.47215808929909*tmp1*tmp9*pow(tmp3,(3.0/2.0));
  out[1] = 3.29184556129898*tmp1*tmp2*tmp3*tmp9;
  out[2] = 3.29184556129898*tmp1*tmp3*tmp5*tmp9;
  out[3] = 3.29184556129898*tmp0*tmp1*tmp4*tmp9;
  out[4] = 5.70164376283992*tmp0*tmp1*tmp2*tmp5*tmp9;
  out[5] = 3.29184556129898*tmp0*tmp1*tmp6*tmp9;
  out[6] = 1.47215808929909*tmp1*tmp9*pow(tmp4,(3.0/2.0));
  out[7] = 3.29184556129898*tmp1*tmp4*tmp5*tmp9;
  out[8] = 3.29184556129898*tmp1*tmp2*tmp6*tmp9;
  out[9] = 1.47215808929909*tmp1*tmp5*tmp6*tmp9;
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
