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
#include "gint1_fn_auto.h"

void gint1_fn_pF(double* a, double a_a, double* p, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 10
  double tmp10, tmp15, tmp4, tmp5, tmp6, tmp8, tmp9, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp4 = v_1*v_1; // auto, weighs 1
  tmp5 = v_0*v_0; // auto, weighs 1
  tmp6 = v_2*v_2; // auto, weighs 1
  tmp15 = pow(a_a,2.25)*exp(-a_a*(tmp4 + tmp5 + tmp6)); // auto, weighs 8
  tmp8 = tmp15*tmp5; // auto, weighs 1
  tmp9 = tmp15*tmp4; // auto, weighs 1
  tmp10 = tmp15*v_2; // auto, weighs 1
  out[0] = v_2*(-2.20823713394864*tmp8 - 2.20823713394864*tmp9) + 1.47215808929909*tmp10*tmp6; // final, weighs 7
  tmp6 = 3.60603613949341*tmp15*tmp6; // auto+recycle, weighs 2
  tmp5 = tmp15*tmp5*v_0; // auto+recycle, weighs 2
  out[1] = -0.901509034873353*tmp5 + v_0*(tmp6 - 0.901509034873353*tmp9); // final, weighs 5
  tmp15 = tmp15*tmp4*v_1; // auto+recycle, weighs 2
  out[2] = -0.901509034873353*tmp15 + v_1*(tmp6 - 0.901509034873353*tmp8); // final, weighs 5
  out[3] = v_2*(2.85082188141996*tmp8 - 2.85082188141996*tmp9); // final, weighs 4
  out[4] = 5.70164376283992*tmp10*v_0*v_1; // final, weighs 3
  out[5] = 1.16384315950667*tmp5 - 3.49152947852002*tmp9*v_0; // final, weighs 4
  out[6] = -1.16384315950667*tmp15 + 3.49152947852002*tmp8*v_1; // final, weighs 4
  // total weight = 58
}

void gint1_fn_pD(double* a, double a_a, double* p, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 7
  double tmp10, tmp2, tmp3, tmp4, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp2 = v_1*v_1; // auto, weighs 1
  tmp3 = v_0*v_0; // auto, weighs 1
  tmp4 = v_2*v_2; // auto, weighs 1
  tmp10 = pow(a_a,1.75)*exp(-a_a*(tmp2 + tmp3 + tmp4)); // auto, weighs 8
  tmp3 = tmp10*tmp3; // auto+recycle, weighs 1
  tmp2 = tmp10*tmp2; // auto+recycle, weighs 1
  out[0] = -0.822961390324745*tmp2 - 0.822961390324745*tmp3 + 1.64592278064949*tmp10*tmp4; // final, weighs 6
  tmp4 = tmp10*v_2; // auto+recycle, weighs 1
  out[1] = 2.85082188141996*tmp4*v_0; // final, weighs 2
  out[2] = 2.85082188141996*tmp4*v_1; // final, weighs 2
  out[3] = 1.42541094070998*tmp3 - 1.42541094070998*tmp2; // final, weighs 3
  out[4] = 2.85082188141996*tmp10*v_0*v_1; // final, weighs 3
  // total weight = 36
}

void gint1_fn_SP(double* a, double a_a, double* p, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 4
  double tmp0, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp0 = exp(-a_a*(v_0*v_0 + v_1*v_1 + v_2*v_2)); // auto, weighs 8
  out[0] = 0.71270547035499*tmp0*pow(a_a,0.75); // final, weighs 4
  tmp0 = tmp0*pow(a_a,1.25); // auto+recycle, weighs 3
  out[1] = 1.42541094070998*tmp0*v_0; // final, weighs 2
  out[2] = 1.42541094070998*tmp0*v_1; // final, weighs 2
  out[3] = 1.42541094070998*tmp0*v_2; // final, weighs 2
  // total weight = 27
}

void gint1_fn_S(double* a, double a_a, double* p, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 3
  double v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  out[0] = 0.71270547035499*pow(a_a,0.75)*exp(-a_a*(v_0*v_0 + v_1*v_1 + v_2*v_2)); // final, weighs 12
  // total weight = 18
}

void gint1_fn_P(double* a, double a_a, double* p, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 4
  double tmp2, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp2 = pow(a_a,1.25)*exp(-a_a*(v_0*v_0 + v_1*v_1 + v_2*v_2)); // auto, weighs 11
  out[0] = 1.42541094070998*tmp2*v_0; // final, weighs 2
  out[1] = 1.42541094070998*tmp2*v_1; // final, weighs 2
  out[2] = 1.42541094070998*tmp2*v_2; // final, weighs 2
  // total weight = 23
}

void gint1_fn_cD(double* a, double a_a, double* p, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 7
  double tmp2, tmp3, tmp4, tmp8, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp2 = v_2*v_2; // auto, weighs 1
  tmp3 = v_1*v_1; // auto, weighs 1
  tmp4 = v_0*v_0; // auto, weighs 1
  tmp8 = pow(a_a,1.75)*exp(-a_a*(tmp2 + tmp3 + tmp4)); // auto, weighs 8
  out[0] = 1.64592278064949*tmp4*tmp8; // final, weighs 2
  out[1] = 2.85082188141996*tmp8*v_0*v_1; // final, weighs 3
  v_2 = tmp8*v_2; // auto+recycle, weighs 1
  out[2] = 2.85082188141996*v_0*v_2; // final, weighs 2
  out[3] = 1.64592278064949*tmp3*tmp8; // final, weighs 2
  out[4] = 2.85082188141996*v_1*v_2; // final, weighs 2
  out[5] = 1.64592278064949*tmp2*tmp8; // final, weighs 2
  // total weight = 31
}

void gint1_fn_cF(double* a, double a_a, double* p, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 9
  double tmp11, tmp2, tmp3, tmp4, tmp6, tmp7, v_0, v_1, v_2;
  v_0 = p[0] - a[0]; // local, weighs 2
  v_1 = p[1] - a[1]; // local, weighs 2
  v_2 = p[2] - a[2]; // local, weighs 2
  tmp2 = v_2*v_2; // auto, weighs 1
  tmp3 = v_1*v_1; // auto, weighs 1
  tmp4 = v_0*v_0; // auto, weighs 1
  tmp11 = pow(a_a,2.25)*exp(-a_a*(tmp2 + tmp3 + tmp4)); // auto, weighs 8
  tmp7 = tmp11*v_0; // auto, weighs 1
  out[0] = 1.47215808929909*tmp4*tmp7; // final, weighs 2
  tmp6 = tmp11*v_1; // auto, weighs 1
  out[1] = 3.29184556129898*tmp4*tmp6; // final, weighs 2
  v_2 = tmp11*v_2; // auto+recycle, weighs 1
  out[2] = 3.29184556129898*tmp4*v_2; // final, weighs 2
  out[3] = 3.29184556129898*tmp3*tmp7; // final, weighs 2
  out[4] = 5.70164376283992*v_0*v_1*v_2; // final, weighs 3
  out[5] = 3.29184556129898*tmp2*tmp7; // final, weighs 2
  out[6] = 1.47215808929909*tmp3*tmp6; // final, weighs 2
  out[7] = 3.29184556129898*tmp3*v_2; // final, weighs 2
  out[8] = 3.29184556129898*tmp2*tmp6; // final, weighs 2
  out[9] = 1.47215808929909*tmp2*v_2; // final, weighs 2
  // total weight = 41
}

typedef void (*fntype)(double*, double, double*, double*);
static const fntype fns[7] = {gint1_fn_pF, gint1_fn_pD, gint1_fn_SP, gint1_fn_S, gint1_fn_P, gint1_fn_cD, gint1_fn_cF};

void gint1_fn_dispatch(int a_s, double* a, double a_a, double* p, double* out)
{
  fns[3+a_s](a, a_a, p, out);
}
