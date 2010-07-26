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
#include "gaux.h"
#define MAX_SHELL 3
#define NUM_SHELL_TYPES 7
#define MAX_SHELL_DOF 10

static void gint2_nai_pF_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_pD_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_SP_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_S_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_P_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_cD_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_cF_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_pF_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_pD_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_SP_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_S_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_P_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_cD_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_cF_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_pF_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_pD_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_SP_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 20
  double ab_a, d_0, d_1, d_2, nai_000_000_1, p_0, p_1, p_2, tmp1, tmp10, tmp11, tmp4, tmp9, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = pow(b_a,0.75); // auto+recycle, weighs 1
  tmp4 = pow(a_a,0.75); // auto, weighs 1
  out[0] = 0.507949087473928*d_1*d_2*tmp4; // final, weighs 3
  nai_000_000_1 = 6.28318530717959*d_0*gaux(usq, 1); // local, weighs 3
  tmp11 = nai_000_000_1*u_0; // auto, weighs 1
  tmp1 = pow(b_a,1.25); // auto, weighs 1
  tmp4 = tmp1*tmp4; // auto+recycle, weighs 1
  out[1] = 1.01589817494786*tmp4*(-tmp11 + d_1*p_0); // final, weighs 5
  tmp10 = nai_000_000_1*u_1; // auto, weighs 1
  out[2] = 1.01589817494786*tmp4*(-tmp10 + d_1*p_1); // final, weighs 5
  tmp9 = nai_000_000_1*u_2; // auto, weighs 1
  out[3] = 1.01589817494786*tmp4*(-tmp9 + d_1*p_2); // final, weighs 5
  tmp11 = -tmp11 + d_1*v_0; // local+recycle, weighs 3
  tmp4 = pow(a_a,1.25); // auto+recycle, weighs 1
  d_2 = d_2*tmp4; // auto+recycle, weighs 1
  out[4] = 1.01589817494786*d_2*tmp11; // final, weighs 2
  d_0 = 6.28318530717959*d_0*gaux(usq, 2); // local+recycle, weighs 3
  usq = nai_000_000_1*v_0 - d_0*u_0; // local+recycle, weighs 4
  ab_a = (0.5*d_1 - 0.5*nai_000_000_1)/ab_a; // auto+recycle, weighs 5
  tmp1 = tmp1*tmp4; // auto+recycle, weighs 1
  out[5] = 2.03179634989571*tmp1*(ab_a + p_0*tmp11 - u_0*usq); // final, weighs 7
  out[6] = 2.03179634989571*tmp1*(p_1*tmp11 - u_1*usq); // final, weighs 6
  out[7] = 2.03179634989571*tmp1*(p_2*tmp11 - u_2*usq); // final, weighs 6
  tmp10 = -tmp10 + d_1*v_1; // local+recycle, weighs 3
  out[8] = 1.01589817494786*d_2*tmp10; // final, weighs 2
  tmp11 = nai_000_000_1*v_1 - d_0*u_1; // local+recycle, weighs 4
  out[9] = 2.03179634989571*tmp1*(p_0*tmp10 - tmp11*u_0); // final, weighs 6
  out[10] = 2.03179634989571*tmp1*(ab_a + p_1*tmp10 - tmp11*u_1); // final, weighs 7
  out[11] = 2.03179634989571*tmp1*(p_2*tmp10 - tmp11*u_2); // final, weighs 6
  d_1 = -tmp9 + d_1*v_2; // local+recycle, weighs 3
  out[12] = 1.01589817494786*d_1*d_2; // final, weighs 2
  d_0 = nai_000_000_1*v_2 - d_0*u_2; // local+recycle, weighs 4
  out[13] = 2.03179634989571*tmp1*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[14] = 2.03179634989571*tmp1*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[15] = 2.03179634989571*tmp1*(ab_a + d_1*p_2 - d_0*u_2); // final, weighs 7
  // total weight = 185
}

static void gint2_nai_S_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 11
  double ab_a, d_0, d_1, d_2, p_0, p_1, p_2, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  ab_a = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_0 = 6.28318530717959*ab_a*gaux(usq, 0); // local+recycle, weighs 3
  d_1 = pow(a_a,0.75); // auto+recycle, weighs 1
  out[0] = 0.507949087473928*d_0*d_1*pow(b_a,0.75); // final, weighs 4
  ab_a = 6.28318530717959*ab_a*gaux(usq, 1); // local+recycle, weighs 3
  d_1 = d_1*pow(b_a,1.25); // auto+recycle, weighs 2
  out[1] = 1.01589817494786*d_1*(d_0*(p_0 - b[0]) - ab_a*u_0); // final, weighs 8
  out[2] = 1.01589817494786*d_1*(d_0*(p_1 - b[1]) - ab_a*u_1); // final, weighs 8
  out[3] = 1.01589817494786*d_1*(d_0*(p_2 - b[2]) - ab_a*u_2); // final, weighs 8
  // total weight = 84
}

static void gint2_nai_P_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 17
  double ab_a, d_0, d_1, d_2, nai_100_000_0, p_0, p_1, p_2, tmp1, tmp6, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  nai_100_000_0 = d_1*v_0 - d_2*u_0; // local, weighs 4
  tmp1 = pow(a_a,1.25); // auto, weighs 1
  tmp6 = tmp1*pow(b_a,0.75); // auto, weighs 2
  out[0] = 1.01589817494786*nai_100_000_0*tmp6; // final, weighs 2
  d_0 = 6.28318530717959*d_0*gaux(usq, 2); // local+recycle, weighs 3
  usq = d_2*v_0 - d_0*u_0; // local+recycle, weighs 4
  ab_a = (0.5*d_1 - 0.5*d_2)/ab_a; // auto+recycle, weighs 5
  tmp1 = tmp1*pow(b_a,1.25); // auto+recycle, weighs 2
  out[1] = 2.03179634989571*tmp1*(ab_a + nai_100_000_0*p_0 - u_0*usq); // final, weighs 7
  out[2] = 2.03179634989571*tmp1*(nai_100_000_0*p_1 - u_1*usq); // final, weighs 6
  out[3] = 2.03179634989571*tmp1*(nai_100_000_0*p_2 - u_2*usq); // final, weighs 6
  nai_100_000_0 = d_1*v_1 - d_2*u_1; // local+recycle, weighs 4
  out[4] = 1.01589817494786*nai_100_000_0*tmp6; // final, weighs 2
  usq = d_2*v_1 - d_0*u_1; // local+recycle, weighs 4
  out[5] = 2.03179634989571*tmp1*(nai_100_000_0*p_0 - u_0*usq); // final, weighs 6
  out[6] = 2.03179634989571*tmp1*(ab_a + nai_100_000_0*p_1 - u_1*usq); // final, weighs 7
  out[7] = 2.03179634989571*tmp1*(nai_100_000_0*p_2 - u_2*usq); // final, weighs 6
  d_1 = d_1*v_2 - d_2*u_2; // local+recycle, weighs 4
  out[8] = 1.01589817494786*d_1*tmp6; // final, weighs 2
  d_0 = d_2*v_2 - d_0*u_2; // local+recycle, weighs 4
  out[9] = 2.03179634989571*tmp1*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[10] = 2.03179634989571*tmp1*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[11] = 2.03179634989571*tmp1*(ab_a + d_1*p_2 - d_0*u_2); // final, weighs 7
  // total weight = 165
}

static void gint2_nai_cD_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 26
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_010_000_2, nai_100_000_0, nai_100_000_1, nai_110_000_0, nai_110_000_1, nai_200_000_0, nai_200_000_1, p_0, p_1, p_2, tmp1, tmp10, tmp4, tmp5, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  nai_100_000_0 = d_1*v_0 - d_2*u_0; // local, weighs 4
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  nai_100_000_1 = d_2*v_0 - nai_000_000_2*u_0; // local, weighs 4
  tmp1 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  nai_200_000_0 = tmp1 + nai_100_000_0*v_0 - nai_100_000_1*u_0; // local, weighs 5
  tmp5 = pow(a_a,1.75); // auto, weighs 1
  tmp10 = tmp5*pow(b_a,0.75); // auto, weighs 2
  out[0] = 1.17305816955079*nai_200_000_0*tmp10; // final, weighs 2
  d_0 = 6.28318530717959*d_0*gaux(usq, 3); // local+recycle, weighs 3
  usq = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto+recycle, weighs 5
  nai_200_000_1 = usq + nai_100_000_1*v_0 - u_0*(nai_000_000_2*v_0 - d_0*u_0); // local, weighs 9
  tmp5 = tmp5*pow(b_a,1.25); // auto+recycle, weighs 2
  out[1] = 2.34611633910158*tmp5*(nai_200_000_0*p_0 + (nai_100_000_0 - nai_100_000_1)/ab_a - nai_200_000_1*u_0); // final, weighs 11
  out[2] = 2.34611633910158*tmp5*(nai_200_000_0*p_1 - nai_200_000_1*u_1); // final, weighs 6
  out[3] = 2.34611633910158*tmp5*(nai_200_000_0*p_2 - nai_200_000_1*u_2); // final, weighs 6
  nai_200_000_0 = d_1*v_1 - d_2*u_1; // local+recycle, weighs 4
  nai_200_000_1 = d_2*v_1 - nai_000_000_2*u_1; // local+recycle, weighs 4
  nai_110_000_0 = nai_200_000_0*v_0 - nai_200_000_1*u_0; // local, weighs 4
  out[4] = 2.03179634989571*nai_110_000_0*tmp10; // final, weighs 2
  nai_010_000_2 = nai_000_000_2*v_1 - d_0*u_1; // local, weighs 4
  nai_110_000_1 = nai_200_000_1*v_0 - nai_010_000_2*u_0; // local, weighs 4
  tmp4 = (0.5*nai_200_000_0 - 0.5*nai_200_000_1)/ab_a; // auto, weighs 5
  out[5] = 4.06359269979142*tmp5*(tmp4 + nai_110_000_0*p_0 - nai_110_000_1*u_0); // final, weighs 7
  nai_100_000_0 = (0.5*nai_100_000_0 - 0.5*nai_100_000_1)/ab_a; // auto+recycle, weighs 5
  out[6] = 4.06359269979142*tmp5*(nai_100_000_0 + nai_110_000_0*p_1 - nai_110_000_1*u_1); // final, weighs 7
  out[7] = 4.06359269979142*tmp5*(nai_110_000_0*p_2 - nai_110_000_1*u_2); // final, weighs 6
  d_1 = d_1*v_2 - d_2*u_2; // local+recycle, weighs 4
  d_2 = d_2*v_2 - nai_000_000_2*u_2; // local+recycle, weighs 4
  nai_100_000_1 = d_1*v_0 - d_2*u_0; // local+recycle, weighs 4
  out[8] = 2.03179634989571*nai_100_000_1*tmp10; // final, weighs 2
  d_0 = nai_000_000_2*v_2 - d_0*u_2; // local+recycle, weighs 4
  nai_000_000_2 = d_2*v_0 - d_0*u_0; // local+recycle, weighs 4
  nai_110_000_0 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto+recycle, weighs 5
  out[9] = 4.06359269979142*tmp5*(nai_110_000_0 + nai_100_000_1*p_0 - nai_000_000_2*u_0); // final, weighs 7
  out[10] = 4.06359269979142*tmp5*(nai_100_000_1*p_1 - nai_000_000_2*u_1); // final, weighs 6
  out[11] = 4.06359269979142*tmp5*(nai_100_000_0 + nai_100_000_1*p_2 - nai_000_000_2*u_2); // final, weighs 7
  nai_000_000_2 = tmp1 + nai_200_000_0*v_1 - nai_200_000_1*u_1; // local+recycle, weighs 5
  out[12] = 1.17305816955079*nai_000_000_2*tmp10; // final, weighs 2
  nai_010_000_2 = usq + nai_200_000_1*v_1 - nai_010_000_2*u_1; // local+recycle, weighs 5
  out[13] = 2.34611633910158*tmp5*(nai_000_000_2*p_0 - nai_010_000_2*u_0); // final, weighs 6
  out[14] = 2.34611633910158*tmp5*(nai_000_000_2*p_1 + (nai_200_000_0 - nai_200_000_1)/ab_a - nai_010_000_2*u_1); // final, weighs 11
  out[15] = 2.34611633910158*tmp5*(nai_000_000_2*p_2 - nai_010_000_2*u_2); // final, weighs 6
  nai_000_000_2 = d_1*v_1 - d_2*u_1; // local+recycle, weighs 4
  out[16] = 2.03179634989571*nai_000_000_2*tmp10; // final, weighs 2
  nai_010_000_2 = d_2*v_1 - d_0*u_1; // local+recycle, weighs 4
  out[17] = 4.06359269979142*tmp5*(nai_000_000_2*p_0 - nai_010_000_2*u_0); // final, weighs 6
  out[18] = 4.06359269979142*tmp5*(nai_110_000_0 + nai_000_000_2*p_1 - nai_010_000_2*u_1); // final, weighs 7
  out[19] = 4.06359269979142*tmp5*(tmp4 + nai_000_000_2*p_2 - nai_010_000_2*u_2); // final, weighs 7
  nai_000_000_2 = tmp1 + d_1*v_2 - d_2*u_2; // local+recycle, weighs 5
  out[20] = 1.17305816955079*nai_000_000_2*tmp10; // final, weighs 2
  d_0 = usq + d_2*v_2 - d_0*u_2; // local+recycle, weighs 5
  out[21] = 2.34611633910158*tmp5*(nai_000_000_2*p_0 - d_0*u_0); // final, weighs 6
  out[22] = 2.34611633910158*tmp5*(nai_000_000_2*p_1 - d_0*u_1); // final, weighs 6
  out[23] = 2.34611633910158*tmp5*(nai_000_000_2*p_2 + (d_1 - d_2)/ab_a - d_0*u_2); // final, weighs 11
  // total weight = 332
}

static void gint2_nai_cF_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 35
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_010_000_1, nai_010_000_2, nai_020_000_2, nai_100_000_0, nai_100_000_1, nai_100_000_2, nai_110_000_0, nai_110_000_1, nai_120_000_1, nai_200_000_0, nai_200_000_1, nai_200_000_2, p_0, p_1, p_2, tmp1, tmp14, tmp2, tmp3, tmp4, tmp7, tmp9, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  nai_100_000_0 = d_1*v_0 - d_2*u_0; // local, weighs 4
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  nai_100_000_1 = d_2*v_0 - nai_000_000_2*u_0; // local, weighs 4
  tmp2 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  nai_200_000_0 = tmp2 + nai_100_000_0*v_0 - nai_100_000_1*u_0; // local, weighs 5
  nai_000_000_3 = 6.28318530717959*d_0*gaux(usq, 3); // local, weighs 3
  nai_100_000_2 = nai_000_000_2*v_0 - nai_000_000_3*u_0; // local, weighs 4
  tmp1 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto, weighs 5
  nai_200_000_1 = tmp1 + nai_100_000_1*v_0 - nai_100_000_2*u_0; // local, weighs 5
  nai_100_000_0 = nai_200_000_0*v_0 + (nai_100_000_0 - nai_100_000_1)/ab_a - nai_200_000_1*u_0; // local+recycle, weighs 9
  tmp9 = pow(a_a,2.25); // auto, weighs 1
  tmp14 = tmp9*pow(b_a,0.75); // auto, weighs 2
  out[0] = 1.04921512347081*nai_100_000_0*tmp14; // final, weighs 2
  d_0 = 6.28318530717959*d_0*gaux(usq, 4); // local+recycle, weighs 3
  usq = (0.5*nai_000_000_2 - 0.5*nai_000_000_3)/ab_a; // auto+recycle, weighs 5
  nai_200_000_2 = usq + nai_100_000_2*v_0 - u_0*(nai_000_000_3*v_0 - d_0*u_0); // local, weighs 9
  nai_100_000_1 = nai_200_000_1*v_0 + (nai_100_000_1 - nai_100_000_2)/ab_a - nai_200_000_2*u_0; // local+recycle, weighs 9
  nai_100_000_2 = tmp9*pow(b_a,1.25); // auto+recycle, weighs 2
  out[1] = 2.09843024694163*nai_100_000_2*(nai_100_000_0*p_0 + (1.5*nai_200_000_0 - 1.5*nai_200_000_1)/ab_a - nai_100_000_1*u_0); // final, weighs 12
  out[2] = 2.09843024694163*nai_100_000_2*(nai_100_000_0*p_1 - nai_100_000_1*u_1); // final, weighs 6
  out[3] = 2.09843024694163*nai_100_000_2*(nai_100_000_0*p_2 - nai_100_000_1*u_2); // final, weighs 6
  nai_100_000_0 = nai_200_000_0*v_1 - nai_200_000_1*u_1; // local+recycle, weighs 4
  out[4] = 2.34611633910158*nai_100_000_0*tmp14; // final, weighs 2
  nai_100_000_1 = nai_200_000_1*v_1 - nai_200_000_2*u_1; // local+recycle, weighs 4
  tmp9 = d_1*v_1 - d_2*u_1; // local+recycle, weighs 4
  nai_010_000_1 = d_2*v_1 - nai_000_000_2*u_1; // local, weighs 4
  nai_110_000_0 = tmp9*v_0 - nai_010_000_1*u_0; // local, weighs 4
  nai_010_000_2 = nai_000_000_2*v_1 - nai_000_000_3*u_1; // local, weighs 4
  nai_110_000_1 = nai_010_000_1*v_0 - nai_010_000_2*u_0; // local, weighs 4
  tmp7 = (nai_110_000_0 - nai_110_000_1)/ab_a; // auto, weighs 4
  out[5] = 4.69223267820315*nai_100_000_2*(tmp7 + nai_100_000_0*p_0 - nai_100_000_1*u_0); // final, weighs 7
  tmp3 = (0.5*nai_200_000_0 - 0.5*nai_200_000_1)/ab_a; // auto, weighs 5
  out[6] = 4.69223267820315*nai_100_000_2*(tmp3 + nai_100_000_0*p_1 - nai_100_000_1*u_1); // final, weighs 7
  out[7] = 4.69223267820315*nai_100_000_2*(nai_100_000_0*p_2 - nai_100_000_1*u_2); // final, weighs 6
  nai_100_000_0 = nai_200_000_0*v_2 - nai_200_000_1*u_2; // local+recycle, weighs 4
  out[8] = 2.34611633910158*nai_100_000_0*tmp14; // final, weighs 2
  nai_100_000_1 = nai_200_000_1*v_2 - nai_200_000_2*u_2; // local+recycle, weighs 4
  d_1 = d_1*v_2 - d_2*u_2; // local+recycle, weighs 4
  d_2 = d_2*v_2 - nai_000_000_2*u_2; // local+recycle, weighs 4
  nai_200_000_0 = d_1*v_0 - d_2*u_0; // local+recycle, weighs 4
  nai_000_000_2 = nai_000_000_2*v_2 - nai_000_000_3*u_2; // local+recycle, weighs 4
  nai_200_000_1 = d_2*v_0 - nai_000_000_2*u_0; // local+recycle, weighs 4
  nai_200_000_2 = (nai_200_000_0 - nai_200_000_1)/ab_a; // auto+recycle, weighs 4
  out[9] = 4.69223267820315*nai_100_000_2*(nai_200_000_2 + nai_100_000_0*p_0 - nai_100_000_1*u_0); // final, weighs 7
  out[10] = 4.69223267820315*nai_100_000_2*(nai_100_000_0*p_1 - nai_100_000_1*u_1); // final, weighs 6
  out[11] = 4.69223267820315*nai_100_000_2*(tmp3 + nai_100_000_0*p_2 - nai_100_000_1*u_2); // final, weighs 7
  nai_100_000_0 = tmp2 + tmp9*v_1 - nai_010_000_1*u_1; // local+recycle, weighs 5
  nai_100_000_1 = tmp1 + nai_010_000_1*v_1 - nai_010_000_2*u_1; // local+recycle, weighs 5
  tmp3 = nai_100_000_0*v_0 - nai_100_000_1*u_0; // local+recycle, weighs 4
  out[12] = 2.34611633910158*tmp14*tmp3; // final, weighs 2
  nai_020_000_2 = usq + nai_010_000_2*v_1 - u_1*(nai_000_000_3*v_1 - d_0*u_1); // local, weighs 9
  nai_120_000_1 = nai_100_000_1*v_0 - nai_020_000_2*u_0; // local, weighs 4
  tmp4 = (0.5*nai_100_000_0 - 0.5*nai_100_000_1)/ab_a; // auto, weighs 5
  out[13] = 4.69223267820315*nai_100_000_2*(tmp4 + p_0*tmp3 - nai_120_000_1*u_0); // final, weighs 7
  out[14] = 4.69223267820315*nai_100_000_2*(tmp7 + p_1*tmp3 - nai_120_000_1*u_1); // final, weighs 7
  out[15] = 4.69223267820315*nai_100_000_2*(p_2*tmp3 - nai_120_000_1*u_2); // final, weighs 6
  nai_120_000_1 = d_1*v_1 - d_2*u_1; // local+recycle, weighs 4
  tmp3 = d_2*v_1 - nai_000_000_2*u_1; // local+recycle, weighs 4
  tmp7 = nai_120_000_1*v_0 - tmp3*u_0; // local+recycle, weighs 4
  out[16] = 4.06359269979142*tmp14*tmp7; // final, weighs 2
  d_0 = nai_000_000_3*v_2 - d_0*u_2; // local+recycle, weighs 4
  nai_000_000_3 = tmp3*v_0 - u_0*(nai_000_000_2*v_1 - d_0*u_1); // local+recycle, weighs 8
  out[17] = 8.12718539958285*nai_100_000_2*(p_0*tmp7 + (0.5*nai_120_000_1 - 0.5*tmp3)/ab_a - nai_000_000_3*u_0); // final, weighs 12
  out[18] = 8.12718539958285*nai_100_000_2*(p_1*tmp7 + (0.5*nai_200_000_0 - 0.5*nai_200_000_1)/ab_a - nai_000_000_3*u_1); // final, weighs 12
  out[19] = 8.12718539958285*nai_100_000_2*(p_2*tmp7 + (0.5*nai_110_000_0 - 0.5*nai_110_000_1)/ab_a - nai_000_000_3*u_2); // final, weighs 12
  nai_000_000_3 = tmp2 + d_1*v_2 - d_2*u_2; // local+recycle, weighs 5
  nai_110_000_0 = tmp1 + d_2*v_2 - nai_000_000_2*u_2; // local+recycle, weighs 5
  nai_110_000_1 = nai_000_000_3*v_0 - nai_110_000_0*u_0; // local+recycle, weighs 4
  out[20] = 2.34611633910158*nai_110_000_1*tmp14; // final, weighs 2
  d_0 = usq + nai_000_000_2*v_2 - d_0*u_2; // local+recycle, weighs 5
  nai_200_000_0 = nai_110_000_0*v_0 - d_0*u_0; // local+recycle, weighs 4
  nai_200_000_1 = (0.5*nai_000_000_3 - 0.5*nai_110_000_0)/ab_a; // auto+recycle, weighs 5
  out[21] = 4.69223267820315*nai_100_000_2*(nai_200_000_1 + nai_110_000_1*p_0 - nai_200_000_0*u_0); // final, weighs 7
  out[22] = 4.69223267820315*nai_100_000_2*(nai_110_000_1*p_1 - nai_200_000_0*u_1); // final, weighs 6
  out[23] = 4.69223267820315*nai_100_000_2*(nai_200_000_2 + nai_110_000_1*p_2 - nai_200_000_0*u_2); // final, weighs 7
  nai_110_000_1 = nai_100_000_0*v_1 + (tmp9 - nai_010_000_1)/ab_a - nai_100_000_1*u_1; // local+recycle, weighs 9
  out[24] = 1.04921512347081*nai_110_000_1*tmp14; // final, weighs 2
  nai_010_000_1 = nai_100_000_1*v_1 + (nai_010_000_1 - nai_010_000_2)/ab_a - nai_020_000_2*u_1; // local+recycle, weighs 9
  out[25] = 2.09843024694163*nai_100_000_2*(nai_110_000_1*p_0 - nai_010_000_1*u_0); // final, weighs 6
  out[26] = 2.09843024694163*nai_100_000_2*(nai_110_000_1*p_1 + (1.5*nai_100_000_0 - 1.5*nai_100_000_1)/ab_a - nai_010_000_1*u_1); // final, weighs 12
  out[27] = 2.09843024694163*nai_100_000_2*(nai_110_000_1*p_2 - nai_010_000_1*u_2); // final, weighs 6
  nai_010_000_1 = nai_100_000_0*v_2 - nai_100_000_1*u_2; // local+recycle, weighs 4
  out[28] = 2.34611633910158*nai_010_000_1*tmp14; // final, weighs 2
  nai_010_000_2 = nai_100_000_1*v_2 - nai_020_000_2*u_2; // local+recycle, weighs 4
  out[29] = 4.69223267820315*nai_100_000_2*(nai_010_000_1*p_0 - nai_010_000_2*u_0); // final, weighs 6
  nai_020_000_2 = (nai_120_000_1 - tmp3)/ab_a; // auto+recycle, weighs 4
  out[30] = 4.69223267820315*nai_100_000_2*(nai_020_000_2 + nai_010_000_1*p_1 - nai_010_000_2*u_1); // final, weighs 7
  out[31] = 4.69223267820315*nai_100_000_2*(tmp4 + nai_010_000_1*p_2 - nai_010_000_2*u_2); // final, weighs 7
  nai_010_000_1 = nai_000_000_3*v_1 - nai_110_000_0*u_1; // local+recycle, weighs 4
  out[32] = 2.34611633910158*nai_010_000_1*tmp14; // final, weighs 2
  nai_010_000_2 = nai_110_000_0*v_1 - d_0*u_1; // local+recycle, weighs 4
  out[33] = 4.69223267820315*nai_100_000_2*(nai_010_000_1*p_0 - nai_010_000_2*u_0); // final, weighs 6
  out[34] = 4.69223267820315*nai_100_000_2*(nai_200_000_1 + nai_010_000_1*p_1 - nai_010_000_2*u_1); // final, weighs 7
  out[35] = 4.69223267820315*nai_100_000_2*(nai_020_000_2 + nai_010_000_1*p_2 - nai_010_000_2*u_2); // final, weighs 7
  d_1 = nai_000_000_3*v_2 + (d_1 - d_2)/ab_a - nai_110_000_0*u_2; // local+recycle, weighs 9
  out[36] = 1.04921512347081*d_1*tmp14; // final, weighs 2
  d_0 = nai_110_000_0*v_2 + (d_2 - nai_000_000_2)/ab_a - d_0*u_2; // local+recycle, weighs 9
  out[37] = 2.09843024694163*nai_100_000_2*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[38] = 2.09843024694163*nai_100_000_2*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[39] = 2.09843024694163*nai_100_000_2*(d_1*p_2 + (1.5*nai_000_000_3 - 1.5*nai_110_000_0)/ab_a - d_0*u_2); // final, weighs 12
  // total weight = 600
}

static void gint2_nai_pF_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_pD_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_SP_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 11
  double ab_a, d_0, d_1, d_2, p_0, p_1, p_2, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  ab_a = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_0 = 6.28318530717959*ab_a*gaux(usq, 0); // local+recycle, weighs 3
  d_1 = pow(b_a,0.75); // auto+recycle, weighs 1
  out[0] = 0.507949087473928*d_0*d_1*pow(a_a,0.75); // final, weighs 4
  ab_a = 6.28318530717959*ab_a*gaux(usq, 1); // local+recycle, weighs 3
  d_1 = d_1*pow(a_a,1.25); // auto+recycle, weighs 2
  out[1] = 1.01589817494786*d_1*(d_0*(p_0 - a[0]) - ab_a*u_0); // final, weighs 8
  out[2] = 1.01589817494786*d_1*(d_0*(p_1 - a[1]) - ab_a*u_1); // final, weighs 8
  out[3] = 1.01589817494786*d_1*(d_0*(p_2 - a[2]) - ab_a*u_2); // final, weighs 8
  // total weight = 84
}

static void gint2_nai_S_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 7
  double ab_a, d_0, d_1, d_2, u_0, u_1, u_2;
  ab_a = a_a + b_a; // local, weighs 1
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = -c[0] + (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 7
  u_1 = -c[1] + (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 7
  u_2 = -c[2] + (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 7
  out[0] = 3.19153824321146*pow(a_a,0.75)*pow(b_a,0.75)*exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)*gaux(ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2), 0)/ab_a; // final, weighs 25
  // total weight = 53
}

static void gint2_nai_P_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 11
  double ab_a, d_0, d_1, d_2, p_0, p_1, p_2, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  ab_a = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_0 = 6.28318530717959*ab_a*gaux(usq, 0); // local+recycle, weighs 3
  ab_a = 6.28318530717959*ab_a*gaux(usq, 1); // local+recycle, weighs 3
  d_1 = pow(a_a,1.25)*pow(b_a,0.75); // auto+recycle, weighs 3
  out[0] = 1.01589817494786*d_1*(d_0*(p_0 - a[0]) - ab_a*u_0); // final, weighs 8
  out[1] = 1.01589817494786*d_1*(d_0*(p_1 - a[1]) - ab_a*u_1); // final, weighs 8
  out[2] = 1.01589817494786*d_1*(d_0*(p_2 - a[2]) - ab_a*u_2); // final, weighs 8
  // total weight = 80
}

static void gint2_nai_cD_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 13
  double ab_a, d_0, d_1, d_2, nai_010_000_0, nai_010_000_1, p_0, p_1, p_2, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  p_0 = p_0 - a[0]; // local+recycle, weighs 2
  p_1 = p_1 - a[1]; // local+recycle, weighs 2
  p_2 = p_2 - a[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  d_0 = 6.28318530717959*d_0*gaux(usq, 2); // local+recycle, weighs 3
  ab_a = (0.5*d_1 - 0.5*d_2)/ab_a; // auto+recycle, weighs 5
  usq = pow(a_a,1.75)*pow(b_a,0.75); // auto+recycle, weighs 3
  out[0] = 1.17305816955079*usq*(ab_a + p_0*(d_1*p_0 - d_2*u_0) - u_0*(d_2*p_0 - d_0*u_0)); // final, weighs 15
  nai_010_000_0 = d_1*p_1 - d_2*u_1; // local, weighs 4
  nai_010_000_1 = d_2*p_1 - d_0*u_1; // local, weighs 4
  out[1] = 2.03179634989571*usq*(nai_010_000_0*p_0 - nai_010_000_1*u_0); // final, weighs 6
  d_1 = d_1*p_2 - d_2*u_2; // local+recycle, weighs 4
  d_0 = d_2*p_2 - d_0*u_2; // local+recycle, weighs 4
  out[2] = 2.03179634989571*usq*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[3] = 1.17305816955079*usq*(ab_a + nai_010_000_0*p_1 - nai_010_000_1*u_1); // final, weighs 7
  out[4] = 2.03179634989571*usq*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[5] = 1.17305816955079*usq*(ab_a + d_1*p_2 - d_0*u_2); // final, weighs 7
  // total weight = 133
}

static void gint2_nai_cF_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 18
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_100_000_0, nai_100_000_1, nai_200_000_0, nai_200_000_1, p_0, p_1, p_2, tmp1, tmp5, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  p_0 = p_0 - a[0]; // local+recycle, weighs 2
  p_1 = p_1 - a[1]; // local+recycle, weighs 2
  p_2 = p_2 - a[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  nai_100_000_0 = d_1*p_0 - d_2*u_0; // local, weighs 4
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  nai_100_000_1 = d_2*p_0 - nai_000_000_2*u_0; // local, weighs 4
  tmp1 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  nai_200_000_0 = tmp1 + nai_100_000_0*p_0 - nai_100_000_1*u_0; // local, weighs 5
  d_0 = 6.28318530717959*d_0*gaux(usq, 3); // local+recycle, weighs 3
  usq = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto+recycle, weighs 5
  nai_200_000_1 = usq + nai_100_000_1*p_0 - u_0*(nai_000_000_2*p_0 - d_0*u_0); // local, weighs 9
  tmp5 = pow(a_a,2.25)*pow(b_a,0.75); // auto, weighs 3
  out[0] = 1.04921512347081*tmp5*(nai_200_000_0*p_0 + (nai_100_000_0 - nai_100_000_1)/ab_a - nai_200_000_1*u_0); // final, weighs 11
  out[1] = 2.34611633910158*tmp5*(nai_200_000_0*p_1 - nai_200_000_1*u_1); // final, weighs 6
  out[2] = 2.34611633910158*tmp5*(nai_200_000_0*p_2 - nai_200_000_1*u_2); // final, weighs 6
  nai_100_000_0 = d_1*p_1 - d_2*u_1; // local+recycle, weighs 4
  nai_100_000_1 = d_2*p_1 - nai_000_000_2*u_1; // local+recycle, weighs 4
  nai_200_000_0 = tmp1 + nai_100_000_0*p_1 - nai_100_000_1*u_1; // local+recycle, weighs 5
  nai_200_000_1 = usq + nai_100_000_1*p_1 - u_1*(nai_000_000_2*p_1 - d_0*u_1); // local+recycle, weighs 9
  out[3] = 2.34611633910158*tmp5*(nai_200_000_0*p_0 - nai_200_000_1*u_0); // final, weighs 6
  d_1 = d_1*p_2 - d_2*u_2; // local+recycle, weighs 4
  d_2 = d_2*p_2 - nai_000_000_2*u_2; // local+recycle, weighs 4
  d_0 = nai_000_000_2*p_2 - d_0*u_2; // local+recycle, weighs 4
  out[4] = 4.06359269979142*tmp5*(p_0*(d_1*p_1 - d_2*u_1) - u_0*(d_2*p_1 - d_0*u_1)); // final, weighs 14
  nai_000_000_2 = tmp1 + d_1*p_2 - d_2*u_2; // local+recycle, weighs 5
  d_0 = usq + d_2*p_2 - d_0*u_2; // local+recycle, weighs 5
  out[5] = 2.34611633910158*tmp5*(nai_000_000_2*p_0 - d_0*u_0); // final, weighs 6
  out[6] = 1.04921512347081*tmp5*(nai_200_000_0*p_1 + (nai_100_000_0 - nai_100_000_1)/ab_a - nai_200_000_1*u_1); // final, weighs 11
  out[7] = 2.34611633910158*tmp5*(nai_200_000_0*p_2 - nai_200_000_1*u_2); // final, weighs 6
  out[8] = 2.34611633910158*tmp5*(nai_000_000_2*p_1 - d_0*u_1); // final, weighs 6
  out[9] = 1.04921512347081*tmp5*(nai_000_000_2*p_2 + (d_1 - d_2)/ab_a - d_0*u_2); // final, weighs 11
  // total weight = 227
}

static void gint2_nai_pF_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_pD_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_SP_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 19
  double ab_a, d_0, d_1, d_2, p_0, p_1, p_2, tmp1, tmp6, tmp7, tmp8, tmp9, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  tmp9 = d_2*u_0; // auto, weighs 1
  tmp1 = pow(b_a,1.25); // auto, weighs 1
  tmp6 = tmp1*pow(a_a,0.75); // auto, weighs 2
  out[0] = 1.01589817494786*tmp6*(-tmp9 + d_1*p_0); // final, weighs 5
  tmp8 = d_2*u_1; // auto, weighs 1
  out[1] = 1.01589817494786*tmp6*(-tmp8 + d_1*p_1); // final, weighs 5
  tmp7 = d_2*u_2; // auto, weighs 1
  out[2] = 1.01589817494786*tmp6*(-tmp7 + d_1*p_2); // final, weighs 5
  tmp6 = -tmp9 + d_1*v_0; // local+recycle, weighs 3
  d_0 = 6.28318530717959*d_0*gaux(usq, 2); // local+recycle, weighs 3
  tmp9 = d_2*v_0 - d_0*u_0; // local+recycle, weighs 4
  ab_a = (0.5*d_1 - 0.5*d_2)/ab_a; // auto+recycle, weighs 5
  tmp1 = tmp1*pow(a_a,1.25); // auto+recycle, weighs 2
  out[3] = 2.03179634989571*tmp1*(ab_a + p_0*tmp6 - tmp9*u_0); // final, weighs 7
  out[4] = 2.03179634989571*tmp1*(p_1*tmp6 - tmp9*u_1); // final, weighs 6
  out[5] = 2.03179634989571*tmp1*(p_2*tmp6 - tmp9*u_2); // final, weighs 6
  tmp6 = -tmp8 + d_1*v_1; // local+recycle, weighs 3
  tmp8 = d_2*v_1 - d_0*u_1; // local+recycle, weighs 4
  out[6] = 2.03179634989571*tmp1*(p_0*tmp6 - tmp8*u_0); // final, weighs 6
  out[7] = 2.03179634989571*tmp1*(ab_a + p_1*tmp6 - tmp8*u_1); // final, weighs 7
  out[8] = 2.03179634989571*tmp1*(p_2*tmp6 - tmp8*u_2); // final, weighs 6
  d_1 = -tmp7 + d_1*v_2; // local+recycle, weighs 3
  d_0 = d_2*v_2 - d_0*u_2; // local+recycle, weighs 4
  out[9] = 2.03179634989571*tmp1*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[10] = 2.03179634989571*tmp1*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[11] = 2.03179634989571*tmp1*(ab_a + d_1*p_2 - d_0*u_2); // final, weighs 7
  // total weight = 174
}

static void gint2_nai_S_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 11
  double ab_a, d_0, d_1, d_2, p_0, p_1, p_2, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  ab_a = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_0 = 6.28318530717959*ab_a*gaux(usq, 0); // local+recycle, weighs 3
  ab_a = 6.28318530717959*ab_a*gaux(usq, 1); // local+recycle, weighs 3
  d_1 = pow(a_a,0.75)*pow(b_a,1.25); // auto+recycle, weighs 3
  out[0] = 1.01589817494786*d_1*(d_0*(p_0 - b[0]) - ab_a*u_0); // final, weighs 8
  out[1] = 1.01589817494786*d_1*(d_0*(p_1 - b[1]) - ab_a*u_1); // final, weighs 8
  out[2] = 1.01589817494786*d_1*(d_0*(p_2 - b[2]) - ab_a*u_2); // final, weighs 8
  // total weight = 80
}

static void gint2_nai_P_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 15
  double ab_a, d_0, d_1, d_2, nai_100_000_0, p_0, p_1, p_2, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  nai_100_000_0 = d_1*v_0 - d_2*u_0; // local, weighs 4
  d_0 = 6.28318530717959*d_0*gaux(usq, 2); // local+recycle, weighs 3
  usq = d_2*v_0 - d_0*u_0; // local+recycle, weighs 4
  ab_a = (0.5*d_1 - 0.5*d_2)/ab_a; // auto+recycle, weighs 5
  v_0 = pow(a_a,1.25)*pow(b_a,1.25); // auto+recycle, weighs 3
  out[0] = 2.03179634989571*v_0*(ab_a + nai_100_000_0*p_0 - u_0*usq); // final, weighs 7
  out[1] = 2.03179634989571*v_0*(nai_100_000_0*p_1 - u_1*usq); // final, weighs 6
  out[2] = 2.03179634989571*v_0*(nai_100_000_0*p_2 - u_2*usq); // final, weighs 6
  nai_100_000_0 = d_1*v_1 - d_2*u_1; // local+recycle, weighs 4
  usq = d_2*v_1 - d_0*u_1; // local+recycle, weighs 4
  out[3] = 2.03179634989571*v_0*(nai_100_000_0*p_0 - u_0*usq); // final, weighs 6
  out[4] = 2.03179634989571*v_0*(ab_a + nai_100_000_0*p_1 - u_1*usq); // final, weighs 7
  out[5] = 2.03179634989571*v_0*(nai_100_000_0*p_2 - u_2*usq); // final, weighs 6
  d_1 = d_1*v_2 - d_2*u_2; // local+recycle, weighs 4
  d_0 = d_2*v_2 - d_0*u_2; // local+recycle, weighs 4
  out[6] = 2.03179634989571*v_0*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[7] = 2.03179634989571*v_0*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[8] = 2.03179634989571*v_0*(ab_a + d_1*p_2 - d_0*u_2); // final, weighs 7
  // total weight = 157
}

static void gint2_nai_cD_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 25
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_010_000_2, nai_100_000_0, nai_100_000_1, nai_110_000_0, nai_110_000_1, nai_200_000_0, nai_200_000_1, p_0, p_1, p_2, tmp1, tmp4, tmp8, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  nai_100_000_0 = d_1*v_0 - d_2*u_0; // local, weighs 4
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  nai_100_000_1 = d_2*v_0 - nai_000_000_2*u_0; // local, weighs 4
  tmp1 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  nai_200_000_0 = tmp1 + nai_100_000_0*v_0 - nai_100_000_1*u_0; // local, weighs 5
  d_0 = 6.28318530717959*d_0*gaux(usq, 3); // local+recycle, weighs 3
  usq = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto+recycle, weighs 5
  nai_200_000_1 = usq + nai_100_000_1*v_0 - u_0*(nai_000_000_2*v_0 - d_0*u_0); // local, weighs 9
  tmp8 = pow(a_a,1.75)*pow(b_a,1.25); // auto, weighs 3
  out[0] = 2.34611633910158*tmp8*(nai_200_000_0*p_0 + (nai_100_000_0 - nai_100_000_1)/ab_a - nai_200_000_1*u_0); // final, weighs 11
  out[1] = 2.34611633910158*tmp8*(nai_200_000_0*p_1 - nai_200_000_1*u_1); // final, weighs 6
  out[2] = 2.34611633910158*tmp8*(nai_200_000_0*p_2 - nai_200_000_1*u_2); // final, weighs 6
  nai_200_000_0 = d_1*v_1 - d_2*u_1; // local+recycle, weighs 4
  nai_200_000_1 = d_2*v_1 - nai_000_000_2*u_1; // local+recycle, weighs 4
  nai_110_000_0 = nai_200_000_0*v_0 - nai_200_000_1*u_0; // local, weighs 4
  nai_010_000_2 = nai_000_000_2*v_1 - d_0*u_1; // local, weighs 4
  nai_110_000_1 = nai_200_000_1*v_0 - nai_010_000_2*u_0; // local, weighs 4
  tmp4 = (0.5*nai_200_000_0 - 0.5*nai_200_000_1)/ab_a; // auto, weighs 5
  out[3] = 4.06359269979142*tmp8*(tmp4 + nai_110_000_0*p_0 - nai_110_000_1*u_0); // final, weighs 7
  nai_100_000_0 = (0.5*nai_100_000_0 - 0.5*nai_100_000_1)/ab_a; // auto+recycle, weighs 5
  out[4] = 4.06359269979142*tmp8*(nai_100_000_0 + nai_110_000_0*p_1 - nai_110_000_1*u_1); // final, weighs 7
  out[5] = 4.06359269979142*tmp8*(nai_110_000_0*p_2 - nai_110_000_1*u_2); // final, weighs 6
  d_1 = d_1*v_2 - d_2*u_2; // local+recycle, weighs 4
  d_2 = d_2*v_2 - nai_000_000_2*u_2; // local+recycle, weighs 4
  nai_100_000_1 = d_1*v_0 - d_2*u_0; // local+recycle, weighs 4
  d_0 = nai_000_000_2*v_2 - d_0*u_2; // local+recycle, weighs 4
  nai_000_000_2 = d_2*v_0 - d_0*u_0; // local+recycle, weighs 4
  nai_110_000_0 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto+recycle, weighs 5
  out[6] = 4.06359269979142*tmp8*(nai_110_000_0 + nai_100_000_1*p_0 - nai_000_000_2*u_0); // final, weighs 7
  out[7] = 4.06359269979142*tmp8*(nai_100_000_1*p_1 - nai_000_000_2*u_1); // final, weighs 6
  out[8] = 4.06359269979142*tmp8*(nai_100_000_0 + nai_100_000_1*p_2 - nai_000_000_2*u_2); // final, weighs 7
  nai_000_000_2 = tmp1 + nai_200_000_0*v_1 - nai_200_000_1*u_1; // local+recycle, weighs 5
  nai_010_000_2 = usq + nai_200_000_1*v_1 - nai_010_000_2*u_1; // local+recycle, weighs 5
  out[9] = 2.34611633910158*tmp8*(nai_000_000_2*p_0 - nai_010_000_2*u_0); // final, weighs 6
  out[10] = 2.34611633910158*tmp8*(nai_000_000_2*p_1 + (nai_200_000_0 - nai_200_000_1)/ab_a - nai_010_000_2*u_1); // final, weighs 11
  out[11] = 2.34611633910158*tmp8*(nai_000_000_2*p_2 - nai_010_000_2*u_2); // final, weighs 6
  nai_000_000_2 = d_1*v_1 - d_2*u_1; // local+recycle, weighs 4
  nai_010_000_2 = d_2*v_1 - d_0*u_1; // local+recycle, weighs 4
  out[12] = 4.06359269979142*tmp8*(nai_000_000_2*p_0 - nai_010_000_2*u_0); // final, weighs 6
  out[13] = 4.06359269979142*tmp8*(nai_110_000_0 + nai_000_000_2*p_1 - nai_010_000_2*u_1); // final, weighs 7
  out[14] = 4.06359269979142*tmp8*(tmp4 + nai_000_000_2*p_2 - nai_010_000_2*u_2); // final, weighs 7
  nai_000_000_2 = tmp1 + d_1*v_2 - d_2*u_2; // local+recycle, weighs 5
  d_0 = usq + d_2*v_2 - d_0*u_2; // local+recycle, weighs 5
  out[15] = 2.34611633910158*tmp8*(nai_000_000_2*p_0 - d_0*u_0); // final, weighs 6
  out[16] = 2.34611633910158*tmp8*(nai_000_000_2*p_1 - d_0*u_1); // final, weighs 6
  out[17] = 2.34611633910158*tmp8*(nai_000_000_2*p_2 + (d_1 - d_2)/ab_a - d_0*u_2); // final, weighs 11
  // total weight = 318
}

static void gint2_nai_cF_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 34
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_010_000_0, nai_010_000_1, nai_010_000_2, nai_020_000_2, nai_100_000_0, nai_100_000_1, nai_100_000_2, nai_110_000_0, nai_110_000_1, nai_120_000_1, nai_200_000_0, nai_200_000_1, nai_200_000_2, p_0, p_1, p_2, tmp1, tmp2, tmp3, tmp4, tmp7, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  nai_100_000_0 = d_1*v_0 - d_2*u_0; // local, weighs 4
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  nai_100_000_1 = d_2*v_0 - nai_000_000_2*u_0; // local, weighs 4
  tmp2 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  nai_200_000_0 = tmp2 + nai_100_000_0*v_0 - nai_100_000_1*u_0; // local, weighs 5
  nai_000_000_3 = 6.28318530717959*d_0*gaux(usq, 3); // local, weighs 3
  nai_100_000_2 = nai_000_000_2*v_0 - nai_000_000_3*u_0; // local, weighs 4
  tmp1 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto, weighs 5
  nai_200_000_1 = tmp1 + nai_100_000_1*v_0 - nai_100_000_2*u_0; // local, weighs 5
  nai_100_000_0 = nai_200_000_0*v_0 + (nai_100_000_0 - nai_100_000_1)/ab_a - nai_200_000_1*u_0; // local+recycle, weighs 9
  d_0 = 6.28318530717959*d_0*gaux(usq, 4); // local+recycle, weighs 3
  usq = (0.5*nai_000_000_2 - 0.5*nai_000_000_3)/ab_a; // auto+recycle, weighs 5
  nai_200_000_2 = usq + nai_100_000_2*v_0 - u_0*(nai_000_000_3*v_0 - d_0*u_0); // local, weighs 9
  nai_100_000_1 = nai_200_000_1*v_0 + (nai_100_000_1 - nai_100_000_2)/ab_a - nai_200_000_2*u_0; // local+recycle, weighs 9
  nai_100_000_2 = pow(a_a,2.25)*pow(b_a,1.25); // auto+recycle, weighs 3
  out[0] = 2.09843024694163*nai_100_000_2*(nai_100_000_0*p_0 + (1.5*nai_200_000_0 - 1.5*nai_200_000_1)/ab_a - nai_100_000_1*u_0); // final, weighs 12
  out[1] = 2.09843024694163*nai_100_000_2*(nai_100_000_0*p_1 - nai_100_000_1*u_1); // final, weighs 6
  out[2] = 2.09843024694163*nai_100_000_2*(nai_100_000_0*p_2 - nai_100_000_1*u_2); // final, weighs 6
  nai_100_000_0 = nai_200_000_0*v_1 - nai_200_000_1*u_1; // local+recycle, weighs 4
  nai_100_000_1 = nai_200_000_1*v_1 - nai_200_000_2*u_1; // local+recycle, weighs 4
  nai_010_000_0 = d_1*v_1 - d_2*u_1; // local, weighs 4
  nai_010_000_1 = d_2*v_1 - nai_000_000_2*u_1; // local, weighs 4
  nai_110_000_0 = nai_010_000_0*v_0 - nai_010_000_1*u_0; // local, weighs 4
  nai_010_000_2 = nai_000_000_2*v_1 - nai_000_000_3*u_1; // local, weighs 4
  nai_110_000_1 = nai_010_000_1*v_0 - nai_010_000_2*u_0; // local, weighs 4
  tmp7 = (nai_110_000_0 - nai_110_000_1)/ab_a; // auto, weighs 4
  out[3] = 4.69223267820315*nai_100_000_2*(tmp7 + nai_100_000_0*p_0 - nai_100_000_1*u_0); // final, weighs 7
  tmp3 = (0.5*nai_200_000_0 - 0.5*nai_200_000_1)/ab_a; // auto, weighs 5
  out[4] = 4.69223267820315*nai_100_000_2*(tmp3 + nai_100_000_0*p_1 - nai_100_000_1*u_1); // final, weighs 7
  out[5] = 4.69223267820315*nai_100_000_2*(nai_100_000_0*p_2 - nai_100_000_1*u_2); // final, weighs 6
  nai_100_000_0 = nai_200_000_0*v_2 - nai_200_000_1*u_2; // local+recycle, weighs 4
  nai_100_000_1 = nai_200_000_1*v_2 - nai_200_000_2*u_2; // local+recycle, weighs 4
  d_1 = d_1*v_2 - d_2*u_2; // local+recycle, weighs 4
  d_2 = d_2*v_2 - nai_000_000_2*u_2; // local+recycle, weighs 4
  nai_200_000_0 = d_1*v_0 - d_2*u_0; // local+recycle, weighs 4
  nai_000_000_2 = nai_000_000_2*v_2 - nai_000_000_3*u_2; // local+recycle, weighs 4
  nai_200_000_1 = d_2*v_0 - nai_000_000_2*u_0; // local+recycle, weighs 4
  nai_200_000_2 = (nai_200_000_0 - nai_200_000_1)/ab_a; // auto+recycle, weighs 4
  out[6] = 4.69223267820315*nai_100_000_2*(nai_200_000_2 + nai_100_000_0*p_0 - nai_100_000_1*u_0); // final, weighs 7
  out[7] = 4.69223267820315*nai_100_000_2*(nai_100_000_0*p_1 - nai_100_000_1*u_1); // final, weighs 6
  out[8] = 4.69223267820315*nai_100_000_2*(tmp3 + nai_100_000_0*p_2 - nai_100_000_1*u_2); // final, weighs 7
  nai_100_000_0 = tmp2 + nai_010_000_0*v_1 - nai_010_000_1*u_1; // local+recycle, weighs 5
  nai_100_000_1 = tmp1 + nai_010_000_1*v_1 - nai_010_000_2*u_1; // local+recycle, weighs 5
  tmp3 = nai_100_000_0*v_0 - nai_100_000_1*u_0; // local+recycle, weighs 4
  nai_020_000_2 = usq + nai_010_000_2*v_1 - u_1*(nai_000_000_3*v_1 - d_0*u_1); // local, weighs 9
  nai_120_000_1 = nai_100_000_1*v_0 - nai_020_000_2*u_0; // local, weighs 4
  tmp4 = (0.5*nai_100_000_0 - 0.5*nai_100_000_1)/ab_a; // auto, weighs 5
  out[9] = 4.69223267820315*nai_100_000_2*(tmp4 + p_0*tmp3 - nai_120_000_1*u_0); // final, weighs 7
  out[10] = 4.69223267820315*nai_100_000_2*(tmp7 + p_1*tmp3 - nai_120_000_1*u_1); // final, weighs 7
  out[11] = 4.69223267820315*nai_100_000_2*(p_2*tmp3 - nai_120_000_1*u_2); // final, weighs 6
  nai_120_000_1 = d_1*v_1 - d_2*u_1; // local+recycle, weighs 4
  tmp3 = d_2*v_1 - nai_000_000_2*u_1; // local+recycle, weighs 4
  tmp7 = nai_120_000_1*v_0 - tmp3*u_0; // local+recycle, weighs 4
  d_0 = nai_000_000_3*v_2 - d_0*u_2; // local+recycle, weighs 4
  nai_000_000_3 = tmp3*v_0 - u_0*(nai_000_000_2*v_1 - d_0*u_1); // local+recycle, weighs 8
  out[12] = 8.12718539958285*nai_100_000_2*(p_0*tmp7 + (0.5*nai_120_000_1 - 0.5*tmp3)/ab_a - nai_000_000_3*u_0); // final, weighs 12
  out[13] = 8.12718539958285*nai_100_000_2*(p_1*tmp7 + (0.5*nai_200_000_0 - 0.5*nai_200_000_1)/ab_a - nai_000_000_3*u_1); // final, weighs 12
  out[14] = 8.12718539958285*nai_100_000_2*(p_2*tmp7 + (0.5*nai_110_000_0 - 0.5*nai_110_000_1)/ab_a - nai_000_000_3*u_2); // final, weighs 12
  nai_000_000_3 = tmp2 + d_1*v_2 - d_2*u_2; // local+recycle, weighs 5
  nai_110_000_0 = tmp1 + d_2*v_2 - nai_000_000_2*u_2; // local+recycle, weighs 5
  nai_110_000_1 = nai_000_000_3*v_0 - nai_110_000_0*u_0; // local+recycle, weighs 4
  d_0 = usq + nai_000_000_2*v_2 - d_0*u_2; // local+recycle, weighs 5
  nai_200_000_0 = nai_110_000_0*v_0 - d_0*u_0; // local+recycle, weighs 4
  nai_200_000_1 = (0.5*nai_000_000_3 - 0.5*nai_110_000_0)/ab_a; // auto+recycle, weighs 5
  out[15] = 4.69223267820315*nai_100_000_2*(nai_200_000_1 + nai_110_000_1*p_0 - nai_200_000_0*u_0); // final, weighs 7
  out[16] = 4.69223267820315*nai_100_000_2*(nai_110_000_1*p_1 - nai_200_000_0*u_1); // final, weighs 6
  out[17] = 4.69223267820315*nai_100_000_2*(nai_200_000_2 + nai_110_000_1*p_2 - nai_200_000_0*u_2); // final, weighs 7
  nai_010_000_0 = nai_100_000_0*v_1 + (nai_010_000_0 - nai_010_000_1)/ab_a - nai_100_000_1*u_1; // local+recycle, weighs 9
  nai_010_000_1 = nai_100_000_1*v_1 + (nai_010_000_1 - nai_010_000_2)/ab_a - nai_020_000_2*u_1; // local+recycle, weighs 9
  out[18] = 2.09843024694163*nai_100_000_2*(nai_010_000_0*p_0 - nai_010_000_1*u_0); // final, weighs 6
  out[19] = 2.09843024694163*nai_100_000_2*(nai_010_000_0*p_1 + (1.5*nai_100_000_0 - 1.5*nai_100_000_1)/ab_a - nai_010_000_1*u_1); // final, weighs 12
  out[20] = 2.09843024694163*nai_100_000_2*(nai_010_000_0*p_2 - nai_010_000_1*u_2); // final, weighs 6
  nai_010_000_0 = nai_100_000_0*v_2 - nai_100_000_1*u_2; // local+recycle, weighs 4
  nai_010_000_1 = nai_100_000_1*v_2 - nai_020_000_2*u_2; // local+recycle, weighs 4
  out[21] = 4.69223267820315*nai_100_000_2*(nai_010_000_0*p_0 - nai_010_000_1*u_0); // final, weighs 6
  nai_010_000_2 = (nai_120_000_1 - tmp3)/ab_a; // auto+recycle, weighs 4
  out[22] = 4.69223267820315*nai_100_000_2*(nai_010_000_2 + nai_010_000_0*p_1 - nai_010_000_1*u_1); // final, weighs 7
  out[23] = 4.69223267820315*nai_100_000_2*(tmp4 + nai_010_000_0*p_2 - nai_010_000_1*u_2); // final, weighs 7
  nai_010_000_0 = nai_000_000_3*v_1 - nai_110_000_0*u_1; // local+recycle, weighs 4
  nai_010_000_1 = nai_110_000_0*v_1 - d_0*u_1; // local+recycle, weighs 4
  out[24] = 4.69223267820315*nai_100_000_2*(nai_010_000_0*p_0 - nai_010_000_1*u_0); // final, weighs 6
  out[25] = 4.69223267820315*nai_100_000_2*(nai_200_000_1 + nai_010_000_0*p_1 - nai_010_000_1*u_1); // final, weighs 7
  out[26] = 4.69223267820315*nai_100_000_2*(nai_010_000_2 + nai_010_000_0*p_2 - nai_010_000_1*u_2); // final, weighs 7
  d_1 = nai_000_000_3*v_2 + (d_1 - d_2)/ab_a - nai_110_000_0*u_2; // local+recycle, weighs 9
  d_0 = nai_110_000_0*v_2 + (d_2 - nai_000_000_2)/ab_a - d_0*u_2; // local+recycle, weighs 9
  out[27] = 2.09843024694163*nai_100_000_2*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[28] = 2.09843024694163*nai_100_000_2*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[29] = 2.09843024694163*nai_100_000_2*(d_1*p_2 + (1.5*nai_000_000_3 - 1.5*nai_110_000_0)/ab_a - d_0*u_2); // final, weighs 12
  // total weight = 578
}

static void gint2_nai_pF_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_pD_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_SP_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 34
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_001_0, nai_000_001_1, nai_000_002_0, nai_000_010_0, nai_000_010_1, nai_000_020_0, nai_000_100_0, nai_000_100_1, nai_000_200_0, nai_000_200_1, p_0, p_1, p_2, tmp1, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp3, tmp8, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  tmp17 = d_2*u_0; // auto, weighs 1
  nai_000_100_0 = -tmp17 + d_1*p_0; // local, weighs 3
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  tmp14 = nai_000_000_2*u_0; // auto, weighs 1
  nai_000_100_1 = -tmp14 + d_2*p_0; // local, weighs 3
  tmp1 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  nai_000_200_0 = tmp1 + nai_000_100_0*p_0 - nai_000_100_1*u_0; // local, weighs 5
  tmp3 = pow(b_a,1.75); // auto, weighs 1
  tmp8 = tmp3*pow(a_a,0.75); // auto, weighs 2
  out[0] = 1.17305816955079*nai_000_200_0*tmp8; // final, weighs 2
  tmp16 = d_2*u_1; // auto, weighs 1
  nai_000_010_0 = -tmp16 + d_1*p_1; // local, weighs 3
  tmp13 = nai_000_000_2*u_1; // auto, weighs 1
  nai_000_010_1 = -tmp13 + d_2*p_1; // local, weighs 3
  out[1] = 2.03179634989571*tmp8*(nai_000_010_0*p_0 - nai_000_010_1*u_0); // final, weighs 6
  tmp15 = d_2*u_2; // auto, weighs 1
  nai_000_001_0 = -tmp15 + d_1*p_2; // local, weighs 3
  tmp12 = nai_000_000_2*u_2; // auto, weighs 1
  nai_000_001_1 = -tmp12 + d_2*p_2; // local, weighs 3
  out[2] = 2.03179634989571*tmp8*(nai_000_001_0*p_0 - nai_000_001_1*u_0); // final, weighs 6
  nai_000_020_0 = tmp1 + nai_000_010_0*p_1 - nai_000_010_1*u_1; // local, weighs 5
  out[3] = 1.17305816955079*nai_000_020_0*tmp8; // final, weighs 2
  out[4] = 2.03179634989571*tmp8*(nai_000_001_0*p_1 - nai_000_001_1*u_1); // final, weighs 6
  nai_000_002_0 = tmp1 + nai_000_001_0*p_2 - nai_000_001_1*u_2; // local, weighs 5
  out[5] = 1.17305816955079*nai_000_002_0*tmp8; // final, weighs 2
  d_0 = 6.28318530717959*d_0*gaux(usq, 3); // local+recycle, weighs 3
  tmp8 = d_0*u_0; // auto+recycle, weighs 1
  usq = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto+recycle, weighs 5
  nai_000_200_1 = usq + nai_000_100_1*p_0 - u_0*(-tmp8 + nai_000_000_2*p_0); // local, weighs 8
  tmp3 = tmp3*pow(a_a,1.25); // auto+recycle, weighs 2
  out[6] = 2.34611633910158*tmp3*(nai_000_200_0*v_0 + (nai_000_100_0 - nai_000_100_1)/ab_a - nai_000_200_1*u_0); // final, weighs 11
  nai_000_100_0 = -tmp17 + d_1*v_0; // local+recycle, weighs 3
  nai_000_100_1 = -tmp14 + d_2*v_0; // local+recycle, weighs 3
  tmp14 = -tmp8 + nai_000_000_2*v_0; // local+recycle, weighs 3
  out[7] = 4.06359269979142*tmp3*(p_0*(nai_000_100_0*p_1 - nai_000_100_1*u_1) + (0.5*nai_000_010_0 - 0.5*nai_000_010_1)/ab_a - u_0*(nai_000_100_1*p_1 - tmp14*u_1)); // final, weighs 20
  nai_000_100_0 = nai_000_100_0*p_2 - nai_000_100_1*u_2; // local+recycle, weighs 4
  nai_000_100_1 = nai_000_100_1*p_2 - tmp14*u_2; // local+recycle, weighs 4
  tmp14 = (0.5*nai_000_001_0 - 0.5*nai_000_001_1)/ab_a; // auto+recycle, weighs 5
  out[8] = 4.06359269979142*tmp3*(tmp14 + nai_000_100_0*p_0 - nai_000_100_1*u_0); // final, weighs 7
  tmp17 = d_0*u_1; // auto+recycle, weighs 1
  tmp8 = usq + nai_000_010_1*p_1 - u_1*(-tmp17 + nai_000_000_2*p_1); // local+recycle, weighs 8
  out[9] = 2.34611633910158*tmp3*(nai_000_020_0*v_0 - tmp8*u_0); // final, weighs 6
  out[10] = 4.06359269979142*tmp3*(nai_000_100_0*p_1 - nai_000_100_1*u_1); // final, weighs 6
  d_0 = d_0*u_2; // auto+recycle, weighs 1
  nai_000_100_0 = usq + nai_000_001_1*p_2 - u_2*(-d_0 + nai_000_000_2*p_2); // local+recycle, weighs 8
  out[11] = 2.34611633910158*tmp3*(nai_000_002_0*v_0 - nai_000_100_0*u_0); // final, weighs 6
  out[12] = 2.34611633910158*tmp3*(nai_000_200_0*v_1 - nai_000_200_1*u_1); // final, weighs 6
  nai_000_100_1 = -tmp16 + d_1*v_1; // local+recycle, weighs 3
  tmp13 = -tmp13 + d_2*v_1; // local+recycle, weighs 3
  tmp16 = -tmp17 + nai_000_000_2*v_1; // local+recycle, weighs 3
  out[13] = 4.06359269979142*tmp3*(p_0*(tmp1 + nai_000_100_1*p_1 - tmp13*u_1) - u_0*(usq + p_1*tmp13 - tmp16*u_1)); // final, weighs 16
  nai_000_100_1 = nai_000_100_1*p_2 - tmp13*u_2; // local+recycle, weighs 4
  tmp13 = p_2*tmp13 - tmp16*u_2; // local+recycle, weighs 4
  out[14] = 4.06359269979142*tmp3*(nai_000_100_1*p_0 - tmp13*u_0); // final, weighs 6
  out[15] = 2.34611633910158*tmp3*(nai_000_020_0*v_1 + (nai_000_010_0 - nai_000_010_1)/ab_a - tmp8*u_1); // final, weighs 11
  out[16] = 4.06359269979142*tmp3*(tmp14 + nai_000_100_1*p_1 - tmp13*u_1); // final, weighs 7
  out[17] = 2.34611633910158*tmp3*(nai_000_002_0*v_1 - nai_000_100_0*u_1); // final, weighs 6
  out[18] = 2.34611633910158*tmp3*(nai_000_200_0*v_2 - nai_000_200_1*u_2); // final, weighs 6
  d_1 = -tmp15 + d_1*v_2; // local+recycle, weighs 3
  d_2 = -tmp12 + d_2*v_2; // local+recycle, weighs 3
  d_0 = -d_0 + nai_000_000_2*v_2; // local+recycle, weighs 3
  out[19] = 4.06359269979142*tmp3*(p_0*(d_1*p_1 - d_2*u_1) - u_0*(d_2*p_1 - d_0*u_1)); // final, weighs 14
  d_1 = tmp1 + d_1*p_2 - d_2*u_2; // local+recycle, weighs 5
  d_0 = usq + d_2*p_2 - d_0*u_2; // local+recycle, weighs 5
  out[20] = 4.06359269979142*tmp3*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[21] = 2.34611633910158*tmp3*(nai_000_020_0*v_2 - tmp8*u_2); // final, weighs 6
  out[22] = 4.06359269979142*tmp3*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[23] = 2.34611633910158*tmp3*(nai_000_002_0*v_2 + (nai_000_001_0 - nai_000_001_1)/ab_a - nai_000_100_0*u_2); // final, weighs 11
  // total weight = 391
}

static void gint2_nai_S_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 13
  double ab_a, d_0, d_1, d_2, nai_000_010_0, nai_000_010_1, p_0, p_1, p_2, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  d_0 = 6.28318530717959*d_0*gaux(usq, 2); // local+recycle, weighs 3
  ab_a = (0.5*d_1 - 0.5*d_2)/ab_a; // auto+recycle, weighs 5
  usq = pow(a_a,0.75)*pow(b_a,1.75); // auto+recycle, weighs 3
  out[0] = 1.17305816955079*usq*(ab_a + p_0*(d_1*p_0 - d_2*u_0) - u_0*(d_2*p_0 - d_0*u_0)); // final, weighs 15
  nai_000_010_0 = d_1*p_1 - d_2*u_1; // local, weighs 4
  nai_000_010_1 = d_2*p_1 - d_0*u_1; // local, weighs 4
  out[1] = 2.03179634989571*usq*(nai_000_010_0*p_0 - nai_000_010_1*u_0); // final, weighs 6
  d_1 = d_1*p_2 - d_2*u_2; // local+recycle, weighs 4
  d_0 = d_2*p_2 - d_0*u_2; // local+recycle, weighs 4
  out[2] = 2.03179634989571*usq*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[3] = 1.17305816955079*usq*(ab_a + nai_000_010_0*p_1 - nai_000_010_1*u_1); // final, weighs 7
  out[4] = 2.03179634989571*usq*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[5] = 1.17305816955079*usq*(ab_a + d_1*p_2 - d_0*u_2); // final, weighs 7
  // total weight = 133
}

static void gint2_nai_P_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 33
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_001_0, nai_000_001_1, nai_000_010_1, nai_000_020_0, nai_000_020_1, nai_000_100_0, nai_000_100_1, nai_000_200_0, nai_000_200_1, p_0, p_1, p_2, tmp0, tmp1, tmp10, tmp11, tmp12, tmp15, tmp2, tmp6, tmp8, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  tmp15 = d_2*u_0; // auto, weighs 1
  nai_000_100_0 = -tmp15 + d_1*p_0; // local, weighs 3
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  tmp12 = nai_000_000_2*u_0; // auto, weighs 1
  nai_000_100_1 = -tmp12 + d_2*p_0; // local, weighs 3
  tmp1 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  nai_000_200_0 = tmp1 + nai_000_100_0*p_0 - nai_000_100_1*u_0; // local, weighs 5
  d_0 = 6.28318530717959*d_0*gaux(usq, 3); // local+recycle, weighs 3
  usq = d_0*u_0; // auto+recycle, weighs 1
  tmp0 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto, weighs 5
  nai_000_200_1 = tmp0 + nai_000_100_1*p_0 - u_0*(-usq + nai_000_000_2*p_0); // local, weighs 8
  tmp6 = pow(a_a,1.25)*pow(b_a,1.75); // auto, weighs 3
  out[0] = 2.34611633910158*tmp6*(nai_000_200_0*v_0 + (nai_000_100_0 - nai_000_100_1)/ab_a - nai_000_200_1*u_0); // final, weighs 11
  nai_000_100_0 = -tmp15 + d_1*v_0; // local+recycle, weighs 3
  nai_000_100_1 = -tmp12 + d_2*v_0; // local+recycle, weighs 3
  tmp12 = -usq + nai_000_000_2*v_0; // local+recycle, weighs 3
  tmp15 = d_2*u_1; // auto+recycle, weighs 1
  usq = -tmp15 + d_1*p_1; // local+recycle, weighs 3
  tmp11 = nai_000_000_2*u_1; // auto, weighs 1
  nai_000_010_1 = -tmp11 + d_2*p_1; // local, weighs 3
  out[1] = 4.06359269979142*tmp6*(p_0*(nai_000_100_0*p_1 - nai_000_100_1*u_1) + (0.5*usq - 0.5*nai_000_010_1)/ab_a - u_0*(nai_000_100_1*p_1 - tmp12*u_1)); // final, weighs 20
  nai_000_100_0 = nai_000_100_0*p_2 - nai_000_100_1*u_2; // local+recycle, weighs 4
  nai_000_100_1 = nai_000_100_1*p_2 - tmp12*u_2; // local+recycle, weighs 4
  tmp12 = d_2*u_2; // auto+recycle, weighs 1
  nai_000_001_0 = -tmp12 + d_1*p_2; // local, weighs 3
  tmp10 = nai_000_000_2*u_2; // auto, weighs 1
  nai_000_001_1 = -tmp10 + d_2*p_2; // local, weighs 3
  tmp2 = (0.5*nai_000_001_0 - 0.5*nai_000_001_1)/ab_a; // auto, weighs 5
  out[2] = 4.06359269979142*tmp6*(tmp2 + nai_000_100_0*p_0 - nai_000_100_1*u_0); // final, weighs 7
  nai_000_020_0 = tmp1 + p_1*usq - nai_000_010_1*u_1; // local, weighs 5
  tmp8 = d_0*u_1; // auto, weighs 1
  nai_000_020_1 = tmp0 + nai_000_010_1*p_1 - u_1*(-tmp8 + nai_000_000_2*p_1); // local, weighs 8
  out[3] = 2.34611633910158*tmp6*(nai_000_020_0*v_0 - nai_000_020_1*u_0); // final, weighs 6
  out[4] = 4.06359269979142*tmp6*(nai_000_100_0*p_1 - nai_000_100_1*u_1); // final, weighs 6
  nai_000_100_0 = tmp1 + nai_000_001_0*p_2 - nai_000_001_1*u_2; // local+recycle, weighs 5
  d_0 = d_0*u_2; // auto+recycle, weighs 1
  nai_000_100_1 = tmp0 + nai_000_001_1*p_2 - u_2*(-d_0 + nai_000_000_2*p_2); // local+recycle, weighs 8
  out[5] = 2.34611633910158*tmp6*(nai_000_100_0*v_0 - nai_000_100_1*u_0); // final, weighs 6
  out[6] = 2.34611633910158*tmp6*(nai_000_200_0*v_1 - nai_000_200_1*u_1); // final, weighs 6
  tmp15 = -tmp15 + d_1*v_1; // local+recycle, weighs 3
  tmp11 = -tmp11 + d_2*v_1; // local+recycle, weighs 3
  tmp8 = -tmp8 + nai_000_000_2*v_1; // local+recycle, weighs 3
  out[7] = 4.06359269979142*tmp6*(p_0*(tmp1 + p_1*tmp15 - tmp11*u_1) - u_0*(tmp0 + p_1*tmp11 - tmp8*u_1)); // final, weighs 16
  tmp15 = p_2*tmp15 - tmp11*u_2; // local+recycle, weighs 4
  tmp11 = p_2*tmp11 - tmp8*u_2; // local+recycle, weighs 4
  out[8] = 4.06359269979142*tmp6*(p_0*tmp15 - tmp11*u_0); // final, weighs 6
  out[9] = 2.34611633910158*tmp6*(nai_000_020_0*v_1 + (usq - nai_000_010_1)/ab_a - nai_000_020_1*u_1); // final, weighs 11
  out[10] = 4.06359269979142*tmp6*(tmp2 + p_1*tmp15 - tmp11*u_1); // final, weighs 7
  out[11] = 2.34611633910158*tmp6*(nai_000_100_0*v_1 - nai_000_100_1*u_1); // final, weighs 6
  out[12] = 2.34611633910158*tmp6*(nai_000_200_0*v_2 - nai_000_200_1*u_2); // final, weighs 6
  d_1 = -tmp12 + d_1*v_2; // local+recycle, weighs 3
  d_2 = -tmp10 + d_2*v_2; // local+recycle, weighs 3
  d_0 = -d_0 + nai_000_000_2*v_2; // local+recycle, weighs 3
  out[13] = 4.06359269979142*tmp6*(p_0*(d_1*p_1 - d_2*u_1) - u_0*(d_2*p_1 - d_0*u_1)); // final, weighs 14
  d_1 = tmp1 + d_1*p_2 - d_2*u_2; // local+recycle, weighs 5
  d_0 = tmp0 + d_2*p_2 - d_0*u_2; // local+recycle, weighs 5
  out[14] = 4.06359269979142*tmp6*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[15] = 2.34611633910158*tmp6*(nai_000_020_0*v_2 - nai_000_020_1*u_2); // final, weighs 6
  out[16] = 4.06359269979142*tmp6*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[17] = 2.34611633910158*tmp6*(nai_000_100_0*v_2 + (nai_000_001_0 - nai_000_001_1)/ab_a - nai_000_100_1*u_2); // final, weighs 11
  // total weight = 365
}

static void gint2_nai_cD_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 48
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_000_002_1, nai_000_002_2, nai_000_020_0, nai_000_020_1, nai_001_010_1, nai_010_000_2, nai_010_000_3, nai_010_001_1, nai_010_010_0, nai_010_010_1, nai_100_000_0, nai_100_000_1, nai_100_000_2, nai_100_001_1, nai_110_000_2, nai_200_000_0, nai_200_000_1, nai_200_000_2, p_0, p_1, p_2, tmp0, tmp1, tmp11, tmp14, tmp17, tmp19, tmp2, tmp21, tmp22, tmp24, tmp29, tmp33, tmp36, tmp39, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  tmp39 = d_2*u_0; // auto, weighs 1
  nai_100_000_0 = -tmp39 + d_1*v_0; // local, weighs 3
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  tmp36 = nai_000_000_2*u_0; // auto, weighs 1
  nai_100_000_1 = -tmp36 + d_2*v_0; // local, weighs 3
  tmp1 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  tmp17 = tmp1 - nai_100_000_1*u_0; // auto, weighs 3
  nai_200_000_0 = tmp17 + nai_100_000_0*v_0; // local, weighs 2
  nai_000_000_3 = 6.28318530717959*d_0*gaux(usq, 3); // local, weighs 3
  tmp33 = nai_000_000_3*u_0; // auto, weighs 1
  nai_100_000_2 = -tmp33 + nai_000_000_2*v_0; // local, weighs 3
  tmp0 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto, weighs 5
  tmp11 = tmp0 - nai_100_000_2*u_0; // auto, weighs 3
  nai_200_000_1 = tmp11 + nai_100_000_1*v_0; // local, weighs 2
  d_0 = 6.28318530717959*d_0*gaux(usq, 4); // local+recycle, weighs 3
  usq = d_0*u_0; // auto+recycle, weighs 1
  tmp2 = (0.5*nai_000_000_2 - 0.5*nai_000_000_3)/ab_a; // auto, weighs 5
  nai_200_000_2 = tmp2 + nai_100_000_2*v_0 - u_0*(-usq + nai_000_000_3*v_0); // local, weighs 8
  tmp14 = 0.5*nai_200_000_0 - 0.5*nai_200_000_1; // auto, weighs 3
  tmp21 = pow(a_a,1.75)*pow(b_a,1.75); // auto, weighs 3
  out[0] = 2.70906179986095*tmp21*(p_0*(nai_200_000_0*p_0 + (nai_100_000_0 - nai_100_000_1)/ab_a - nai_200_000_1*u_0) + (tmp14 + tmp17 - tmp11 + nai_100_000_0*p_0 - nai_100_000_1*p_0)/ab_a - u_0*(nai_200_000_1*p_0 + (nai_100_000_1 - nai_100_000_2)/ab_a - nai_200_000_2*u_0)); // final, weighs 35
  tmp11 = nai_200_000_0*p_1 - nai_200_000_1*u_1; // local+recycle, weighs 4
  tmp17 = nai_200_000_1*p_1 - nai_200_000_2*u_1; // local+recycle, weighs 4
  out[1] = 4.69223267820315*tmp21*(p_0*tmp11 + (nai_100_000_0*p_1 + nai_100_000_2*u_1 - nai_100_000_1*p_1 - nai_100_000_1*u_1)/ab_a - tmp17*u_0); // final, weighs 18
  nai_200_000_0 = nai_200_000_0*p_2 - nai_200_000_1*u_2; // local+recycle, weighs 4
  nai_200_000_1 = nai_200_000_1*p_2 - nai_200_000_2*u_2; // local+recycle, weighs 4
  nai_200_000_2 = nai_100_000_0*p_2 - nai_100_000_1*u_2; // local+recycle, weighs 4
  nai_100_001_1 = nai_100_000_1*p_2 - nai_100_000_2*u_2; // local, weighs 4
  out[2] = 4.69223267820315*tmp21*(nai_200_000_0*p_0 + (nai_200_000_2 - nai_100_001_1)/ab_a - nai_200_000_1*u_0); // final, weighs 11
  tmp14 = tmp14/ab_a; // auto+recycle, weighs 2
  out[3] = 2.70906179986095*tmp21*(tmp14 + p_1*tmp11 - tmp17*u_1); // final, weighs 7
  out[4] = 4.69223267820315*tmp21*(nai_200_000_0*p_1 - nai_200_000_1*u_1); // final, weighs 6
  out[5] = 2.70906179986095*tmp21*(tmp14 + nai_200_000_0*p_2 - nai_200_000_1*u_2); // final, weighs 7
  nai_200_000_0 = -tmp36 + d_2*p_0; // local+recycle, weighs 3
  nai_200_000_1 = tmp1 + p_0*(-tmp39 + d_1*p_0) - nai_200_000_0*u_0; // local+recycle, weighs 8
  tmp11 = -tmp33 + nai_000_000_2*p_0; // local+recycle, weighs 3
  nai_200_000_0 = tmp0 + nai_200_000_0*p_0 - tmp11*u_0; // local+recycle, weighs 5
  tmp11 = tmp2 + p_0*tmp11 - u_0*(-usq + nai_000_000_3*p_0); // local+recycle, weighs 8
  tmp14 = d_2*u_1; // auto+recycle, weighs 1
  tmp17 = -tmp14 + d_1*v_1; // local+recycle, weighs 3
  tmp33 = nai_000_000_2*u_1; // auto+recycle, weighs 1
  tmp36 = -tmp33 + d_2*v_1; // local+recycle, weighs 3
  tmp39 = tmp36*u_0; // auto+recycle, weighs 1
  usq = nai_000_000_3*u_1; // auto+recycle, weighs 1
  nai_010_000_2 = -usq + nai_000_000_2*v_1; // local, weighs 3
  tmp22 = nai_010_000_2*u_0; // auto, weighs 1
  out[6] = 4.69223267820315*tmp21*(v_0*(nai_200_000_1*v_1 - nai_200_000_0*u_1) + (tmp22 - tmp39 + p_0*tmp17 - p_0*tmp36)/ab_a - u_0*(nai_200_000_0*v_1 - tmp11*u_1)); // final, weighs 24
  tmp39 = -tmp39 + tmp17*v_0; // local+recycle, weighs 3
  tmp22 = -tmp22 + tmp36*v_0; // local+recycle, weighs 3
  nai_100_000_0 = (0.5*nai_100_000_0 - 0.5*nai_100_000_1)/ab_a; // auto+recycle, weighs 5
  tmp29 = d_0*u_1; // auto, weighs 1
  nai_010_000_3 = -tmp29 + nai_000_000_3*v_1; // local, weighs 3
  nai_110_000_2 = nai_010_000_2*v_0 - nai_010_000_3*u_0; // local, weighs 4
  nai_100_000_1 = (0.5*nai_100_000_1 - 0.5*nai_100_000_2)/ab_a; // auto+recycle, weighs 5
  nai_100_000_2 = tmp1 - tmp36*u_1; // auto+recycle, weighs 3
  nai_010_010_0 = nai_100_000_2 + p_1*tmp17; // local, weighs 2
  tmp19 = tmp0 - nai_010_000_2*u_1; // auto, weighs 3
  nai_010_010_1 = tmp19 + p_1*tmp36; // local, weighs 2
  out[7] = 8.12718539958285*tmp21*(p_0*(nai_100_000_0 + p_1*tmp39 - tmp22*u_1) + (0.5*nai_010_010_0 - 0.5*nai_010_010_1)/ab_a - u_0*(nai_100_000_1 + p_1*tmp22 - nai_110_000_2*u_1)); // final, weighs 22
  tmp39 = p_2*tmp39 - tmp22*u_2; // local+recycle, weighs 4
  nai_110_000_2 = p_2*tmp22 - nai_110_000_2*u_2; // local+recycle, weighs 4
  tmp22 = p_2*tmp17 - tmp36*u_2; // local+recycle, weighs 4
  nai_010_001_1 = p_2*tmp36 - nai_010_000_2*u_2; // local, weighs 4
  out[8] = 8.12718539958285*tmp21*(p_0*tmp39 + (0.5*tmp22 - 0.5*nai_010_001_1)/ab_a - nai_110_000_2*u_0); // final, weighs 12
  tmp14 = -tmp14 + d_1*p_1; // local+recycle, weighs 3
  tmp33 = -tmp33 + d_2*p_1; // local+recycle, weighs 3
  nai_000_020_0 = tmp1 + p_1*tmp14 - tmp33*u_1; // local, weighs 5
  usq = -usq + nai_000_000_2*p_1; // local+recycle, weighs 3
  nai_000_020_1 = tmp0 + p_1*tmp33 - u_1*usq; // local, weighs 5
  tmp29 = tmp2 + p_1*usq - u_1*(-tmp29 + nai_000_000_3*p_1); // local+recycle, weighs 8
  out[9] = 4.69223267820315*tmp21*(v_0*(nai_000_020_0*v_1 + (tmp14 - tmp33)/ab_a - nai_000_020_1*u_1) - u_0*(nai_000_020_1*v_1 + (tmp33 - usq)/ab_a - tmp29*u_1)); // final, weighs 24
  out[10] = 8.12718539958285*tmp21*(p_1*tmp39 + (0.5*nai_200_000_2 - 0.5*nai_100_001_1)/ab_a - nai_110_000_2*u_1); // final, weighs 12
  nai_100_001_1 = d_2*u_2; // auto+recycle, weighs 1
  nai_110_000_2 = -nai_100_001_1 + d_1*p_2; // local+recycle, weighs 3
  nai_200_000_2 = nai_000_000_2*u_2; // auto+recycle, weighs 1
  tmp14 = -nai_200_000_2 + d_2*p_2; // local+recycle, weighs 3
  tmp33 = tmp1 + nai_110_000_2*p_2 - tmp14*u_2; // local+recycle, weighs 5
  tmp39 = nai_000_000_3*u_2; // auto+recycle, weighs 1
  usq = -tmp39 + nai_000_000_2*p_2; // local+recycle, weighs 3
  nai_000_002_1 = tmp0 + p_2*tmp14 - u_2*usq; // local, weighs 5
  d_0 = d_0*u_2; // auto+recycle, weighs 1
  nai_000_002_2 = tmp2 + p_2*usq - u_2*(-d_0 + nai_000_000_3*p_2); // local, weighs 8
  out[11] = 4.69223267820315*tmp21*(v_0*(tmp33*v_1 - nai_000_002_1*u_1) - u_0*(nai_000_002_1*v_1 - nai_000_002_2*u_1)); // final, weighs 14
  nai_200_000_1 = nai_200_000_1*v_2 - nai_200_000_0*u_2; // local+recycle, weighs 4
  nai_200_000_0 = nai_200_000_0*v_2 - tmp11*u_2; // local+recycle, weighs 4
  d_1 = -nai_100_001_1 + d_1*v_2; // local+recycle, weighs 3
  d_2 = -nai_200_000_2 + d_2*v_2; // local+recycle, weighs 3
  nai_100_001_1 = d_2*u_0; // auto+recycle, weighs 1
  nai_000_000_2 = -tmp39 + nai_000_000_2*v_2; // local+recycle, weighs 3
  nai_200_000_2 = nai_000_000_2*u_0; // auto+recycle, weighs 1
  out[12] = 4.69223267820315*tmp21*(nai_200_000_1*v_0 + (nai_200_000_2 - nai_100_001_1 + d_1*p_0 - d_2*p_0)/ab_a - nai_200_000_0*u_0); // final, weighs 16
  nai_100_001_1 = -nai_100_001_1 + d_1*v_0; // local+recycle, weighs 3
  nai_200_000_2 = -nai_200_000_2 + d_2*v_0; // local+recycle, weighs 3
  d_0 = -d_0 + nai_000_000_3*v_2; // local+recycle, weighs 3
  nai_000_000_3 = nai_000_000_2*v_0 - d_0*u_0; // local+recycle, weighs 4
  tmp11 = d_2*u_1; // auto+recycle, weighs 1
  tmp39 = -tmp11 + d_1*p_1; // local+recycle, weighs 3
  tmp24 = nai_000_000_2*u_1; // auto, weighs 1
  nai_001_010_1 = -tmp24 + d_2*p_1; // local, weighs 3
  out[13] = 8.12718539958285*tmp21*(p_0*(nai_100_001_1*p_1 - nai_200_000_2*u_1) + (0.5*tmp39 - 0.5*nai_001_010_1)/ab_a - u_0*(nai_200_000_2*p_1 - nai_000_000_3*u_1)); // final, weighs 20
  nai_100_000_0 = nai_100_000_0 + nai_100_001_1*p_2 - nai_200_000_2*u_2; // local+recycle, weighs 5
  nai_000_000_3 = nai_100_000_1 + nai_200_000_2*p_2 - nai_000_000_3*u_2; // local+recycle, weighs 5
  nai_100_000_1 = tmp1 - d_2*u_2; // auto+recycle, weighs 3
  nai_100_001_1 = nai_100_000_1 + d_1*p_2; // local+recycle, weighs 2
  nai_200_000_2 = tmp0 - nai_000_000_2*u_2; // auto+recycle, weighs 3
  tmp0 = nai_200_000_2 + d_2*p_2; // local+recycle, weighs 2
  tmp1 = (0.5*nai_100_001_1 - 0.5*tmp0)/ab_a; // auto+recycle, weighs 5
  out[14] = 8.12718539958285*tmp21*(tmp1 + nai_100_000_0*p_0 - nai_000_000_3*u_0); // final, weighs 7
  nai_000_020_0 = nai_000_020_0*v_2 - nai_000_020_1*u_2; // local+recycle, weighs 4
  nai_000_020_1 = nai_000_020_1*v_2 - tmp29*u_2; // local+recycle, weighs 4
  out[15] = 4.69223267820315*tmp21*(nai_000_020_0*v_0 - nai_000_020_1*u_0); // final, weighs 6
  out[16] = 8.12718539958285*tmp21*(nai_100_000_0*p_1 - nai_000_000_3*u_1); // final, weighs 6
  nai_000_000_3 = tmp33*v_2 + (nai_110_000_2 - tmp14)/ab_a - nai_000_002_1*u_2; // local+recycle, weighs 9
  nai_000_002_1 = nai_000_002_1*v_2 + (tmp14 - usq)/ab_a - nai_000_002_2*u_2; // local+recycle, weighs 9
  out[17] = 4.69223267820315*tmp21*(nai_000_000_3*v_0 - nai_000_002_1*u_0); // final, weighs 6
  nai_000_002_2 = nai_100_000_2 + tmp17*v_1; // local+recycle, weighs 2
  nai_100_000_0 = tmp19 + tmp36*v_1; // local+recycle, weighs 2
  nai_010_000_3 = tmp2 + nai_010_000_2*v_1 - nai_010_000_3*u_1; // local+recycle, weighs 5
  nai_100_000_2 = 0.5*nai_000_002_2 - 0.5*nai_100_000_0; // auto+recycle, weighs 3
  nai_110_000_2 = nai_100_000_2/ab_a; // auto+recycle, weighs 2
  out[18] = 2.70906179986095*tmp21*(nai_110_000_2 + p_0*(nai_000_002_2*p_0 - nai_100_000_0*u_0) - u_0*(nai_100_000_0*p_0 - nai_010_000_3*u_0)); // final, weighs 15
  tmp14 = nai_000_002_2*p_1 + (tmp17 - tmp36)/ab_a - nai_100_000_0*u_1; // local+recycle, weighs 9
  tmp19 = nai_100_000_0*p_1 + (tmp36 - nai_010_000_2)/ab_a - nai_010_000_3*u_1; // local+recycle, weighs 9
  out[19] = 4.69223267820315*tmp21*(p_0*tmp14 - tmp19*u_0); // final, weighs 6
  nai_000_002_2 = nai_000_002_2*p_2 - nai_100_000_0*u_2; // local+recycle, weighs 4
  nai_010_000_3 = nai_100_000_0*p_2 - nai_010_000_3*u_2; // local+recycle, weighs 4
  out[20] = 4.69223267820315*tmp21*(nai_000_002_2*p_0 - nai_010_000_3*u_0); // final, weighs 6
  out[21] = 2.70906179986095*tmp21*(p_1*tmp14 + (nai_010_010_0 + nai_100_000_2 - nai_010_010_1)/ab_a - tmp19*u_1); // final, weighs 12
  out[22] = 4.69223267820315*tmp21*(nai_000_002_2*p_1 + (tmp22 - nai_010_001_1)/ab_a - nai_010_000_3*u_1); // final, weighs 11
  out[23] = 2.70906179986095*tmp21*(nai_110_000_2 + nai_000_002_2*p_2 - nai_010_000_3*u_2); // final, weighs 7
  out[24] = 4.69223267820315*tmp21*(nai_200_000_1*v_1 - nai_200_000_0*u_1); // final, weighs 6
  nai_000_002_2 = -tmp11 + d_1*v_1; // local+recycle, weighs 3
  nai_010_000_3 = -tmp24 + d_2*v_1; // local+recycle, weighs 3
  nai_010_001_1 = nai_000_000_2*v_1 - d_0*u_1; // local+recycle, weighs 4
  out[25] = 8.12718539958285*tmp21*(p_0*(nai_000_002_2*p_1 + (0.5*d_1 - 0.5*d_2)/ab_a - nai_010_000_3*u_1) - u_0*(nai_010_000_3*p_1 + (0.5*d_2 - 0.5*nai_000_000_2)/ab_a - nai_010_001_1*u_1)); // final, weighs 26
  nai_000_002_2 = nai_000_002_2*p_2 + (0.5*tmp17 - 0.5*tmp36)/ab_a - nai_010_000_3*u_2; // local+recycle, weighs 10
  nai_010_000_2 = nai_010_000_3*p_2 + (0.5*tmp36 - 0.5*nai_010_000_2)/ab_a - nai_010_001_1*u_2; // local+recycle, weighs 10
  out[26] = 8.12718539958285*tmp21*(nai_000_002_2*p_0 - nai_010_000_2*u_0); // final, weighs 6
  out[27] = 4.69223267820315*tmp21*(nai_000_020_0*v_1 + (tmp39 - nai_001_010_1)/ab_a - nai_000_020_1*u_1); // final, weighs 11
  out[28] = 8.12718539958285*tmp21*(tmp1 + nai_000_002_2*p_1 - nai_010_000_2*u_1); // final, weighs 7
  out[29] = 4.69223267820315*tmp21*(nai_000_000_3*v_1 - nai_000_002_1*u_1); // final, weighs 6
  nai_000_000_3 = nai_100_000_1 + d_1*v_2; // local+recycle, weighs 2
  nai_000_002_1 = nai_200_000_2 + d_2*v_2; // local+recycle, weighs 2
  d_0 = tmp2 + nai_000_000_2*v_2 - d_0*u_2; // local+recycle, weighs 5
  nai_000_002_2 = 0.5*nai_000_000_3 - 0.5*nai_000_002_1; // auto+recycle, weighs 3
  nai_000_020_0 = nai_000_002_2/ab_a; // auto+recycle, weighs 2
  out[30] = 2.70906179986095*tmp21*(nai_000_020_0 + p_0*(nai_000_000_3*p_0 - nai_000_002_1*u_0) - u_0*(nai_000_002_1*p_0 - d_0*u_0)); // final, weighs 15
  nai_000_020_1 = nai_000_000_3*p_1 - nai_000_002_1*u_1; // local+recycle, weighs 4
  nai_001_010_1 = nai_000_002_1*p_1 - d_0*u_1; // local+recycle, weighs 4
  out[31] = 4.69223267820315*tmp21*(nai_000_020_1*p_0 - nai_001_010_1*u_0); // final, weighs 6
  d_1 = nai_000_000_3*p_2 + (d_1 - d_2)/ab_a - nai_000_002_1*u_2; // local+recycle, weighs 9
  d_0 = nai_000_002_1*p_2 + (d_2 - nai_000_000_2)/ab_a - d_0*u_2; // local+recycle, weighs 9
  out[32] = 4.69223267820315*tmp21*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[33] = 2.70906179986095*tmp21*(nai_000_020_0 + nai_000_020_1*p_1 - nai_001_010_1*u_1); // final, weighs 7
  out[34] = 4.69223267820315*tmp21*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[35] = 2.70906179986095*tmp21*(d_1*p_2 + (nai_000_002_2 + nai_100_001_1 - tmp0)/ab_a - d_0*u_2); // final, weighs 12
  // total weight = 926
}

static void gint2_nai_cF_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 81
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_000_000_4, nai_001_000_3, nai_010_000_2, nai_010_000_3, nai_011_000_1, nai_011_000_2, nai_020_000_3, nai_020_001_0, nai_020_001_1, nai_020_010_1, nai_020_200_0, nai_020_200_1, nai_100_000_0, nai_100_000_1, nai_100_000_2, nai_100_000_3, nai_101_000_1, nai_101_000_2, nai_110_000_1, nai_110_000_2, nai_200_000_0, nai_200_000_1, nai_200_000_2, nai_200_000_3, nai_200_001_0, nai_200_001_1, nai_200_010_0, nai_200_010_1, nai_300_000_0, nai_300_000_1, nai_300_000_2, nai_300_010_0, nai_300_010_1, p_0, p_1, p_2, tmp0, tmp1, tmp2, tmp3, tmp30, tmp31, tmp35, tmp40, tmp47, tmp49, tmp51, tmp52, tmp58, tmp6, tmp61, tmp62, tmp64, tmp66, tmp67, tmp7, tmp71, tmp72, tmp79, tmp80, tmp83, tmp85, tmp86, tmp87, tmp89, tmp90, tmp93, tmp96, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  tmp96 = d_2*u_0; // auto, weighs 1
  nai_100_000_0 = -tmp96 + d_1*v_0; // local, weighs 3
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  tmp93 = nai_000_000_2*u_0; // auto, weighs 1
  nai_100_000_1 = -tmp93 + d_2*v_0; // local, weighs 3
  tmp2 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  tmp49 = tmp2 - nai_100_000_1*u_0; // auto, weighs 3
  nai_200_000_0 = tmp49 + nai_100_000_0*v_0; // local, weighs 2
  nai_000_000_3 = 6.28318530717959*d_0*gaux(usq, 3); // local, weighs 3
  tmp90 = nai_000_000_3*u_0; // auto, weighs 1
  nai_100_000_2 = -tmp90 + nai_000_000_2*v_0; // local, weighs 3
  tmp1 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto, weighs 5
  tmp31 = tmp1 - nai_100_000_2*u_0; // auto, weighs 3
  nai_200_000_1 = tmp31 + nai_100_000_1*v_0; // local, weighs 2
  tmp51 = (nai_100_000_0 - nai_100_000_1)/ab_a - nai_200_000_1*u_0; // auto, weighs 7
  nai_300_000_0 = tmp51 + nai_200_000_0*v_0; // local, weighs 2
  nai_000_000_4 = 6.28318530717959*d_0*gaux(usq, 4); // local, weighs 3
  tmp87 = nai_000_000_4*u_0; // auto, weighs 1
  nai_100_000_3 = -tmp87 + nai_000_000_3*v_0; // local, weighs 3
  tmp0 = (0.5*nai_000_000_2 - 0.5*nai_000_000_3)/ab_a; // auto, weighs 5
  tmp30 = tmp0 - nai_100_000_3*u_0; // auto, weighs 3
  nai_200_000_2 = tmp30 + nai_100_000_2*v_0; // local, weighs 2
  tmp52 = (nai_100_000_1 - nai_100_000_2)/ab_a - nai_200_000_2*u_0; // auto, weighs 7
  nai_300_000_1 = tmp52 + nai_200_000_1*v_0; // local, weighs 2
  d_0 = 6.28318530717959*d_0*gaux(usq, 5); // local+recycle, weighs 3
  usq = d_0*u_0; // auto+recycle, weighs 1
  tmp3 = (0.5*nai_000_000_3 - 0.5*nai_000_000_4)/ab_a; // auto, weighs 5
  nai_200_000_3 = tmp3 + nai_100_000_3*v_0 - u_0*(-usq + nai_000_000_4*v_0); // local, weighs 8
  nai_100_000_3 = (nai_100_000_2 - nai_100_000_3)/ab_a - nai_200_000_3*u_0; // auto+recycle, weighs 7
  nai_300_000_2 = nai_100_000_3 + nai_200_000_2*v_0; // local, weighs 2
  tmp51 = tmp51 + nai_200_000_0*p_0; // local+recycle, weighs 2
  tmp52 = tmp52 + nai_200_000_1*p_0; // local+recycle, weighs 2
  tmp35 = 0.5*nai_300_000_0 - 0.5*nai_300_000_1; // auto, weighs 3
  tmp58 = pow(a_a,2.25)*pow(b_a,1.75); // auto, weighs 3
  out[0] = 2.4230585358948*tmp58*(p_0*(nai_300_000_0*p_0 + (1.5*nai_200_000_0 - 1.5*nai_200_000_1)/ab_a - nai_300_000_1*u_0) + (tmp35 + 1.5*tmp51 - 1.5*tmp52)/ab_a - u_0*(nai_300_000_1*p_0 + (1.5*nai_200_000_1 - 1.5*nai_200_000_2)/ab_a - nai_300_000_2*u_0)); // final, weighs 33
  nai_300_010_0 = nai_300_000_0*p_1 - nai_300_000_1*u_1; // local, weighs 4
  nai_300_010_1 = nai_300_000_1*p_1 - nai_300_000_2*u_1; // local, weighs 4
  tmp64 = nai_200_000_1*u_1; // auto, weighs 1
  nai_200_010_0 = -tmp64 + nai_200_000_0*p_1; // local, weighs 3
  tmp62 = nai_200_000_2*u_1; // auto, weighs 1
  nai_200_010_1 = -tmp62 + nai_200_000_1*p_1; // local, weighs 3
  out[1] = 4.19686049388326*tmp58*(nai_300_010_0*p_0 + (1.5*nai_200_010_0 - 1.5*nai_200_010_1)/ab_a - nai_300_010_1*u_0); // final, weighs 12
  nai_300_000_0 = nai_300_000_0*p_2 - nai_300_000_1*u_2; // local+recycle, weighs 4
  nai_300_000_1 = nai_300_000_1*p_2 - nai_300_000_2*u_2; // local+recycle, weighs 4
  nai_300_000_2 = nai_200_000_1*u_2; // auto+recycle, weighs 1
  nai_200_001_0 = -nai_300_000_2 + nai_200_000_0*p_2; // local, weighs 3
  tmp61 = nai_200_000_2*u_2; // auto, weighs 1
  nai_200_001_1 = -tmp61 + nai_200_000_1*p_2; // local, weighs 3
  out[2] = 4.19686049388326*tmp58*(nai_300_000_0*p_0 + (1.5*nai_200_001_0 - 1.5*nai_200_001_1)/ab_a - nai_300_000_1*u_0); // final, weighs 12
  tmp35 = tmp35/ab_a; // auto+recycle, weighs 2
  out[3] = 2.4230585358948*tmp58*(tmp35 + nai_300_010_0*p_1 - nai_300_010_1*u_1); // final, weighs 7
  out[4] = 4.19686049388326*tmp58*(nai_300_000_0*p_1 - nai_300_000_1*u_1); // final, weighs 6
  out[5] = 2.4230585358948*tmp58*(tmp35 + nai_300_000_0*p_2 - nai_300_000_1*u_2); // final, weighs 7
  nai_300_000_0 = tmp31 + nai_100_000_1*p_0; // local+recycle, weighs 2
  nai_300_000_1 = 0.5*nai_200_000_0 - 0.5*nai_200_000_1; // auto+recycle, weighs 3
  nai_300_010_0 = p_0*tmp51 + (nai_300_000_1 + tmp49 - nai_300_000_0 + nai_100_000_0*p_0)/ab_a - tmp52*u_0; // local+recycle, weighs 12
  nai_300_010_1 = 0.5*nai_200_000_1 - 0.5*nai_200_000_2; // auto+recycle, weighs 3
  nai_100_000_3 = p_0*tmp52 + (nai_300_000_0 + nai_300_010_1 - tmp30 - nai_100_000_2*p_0)/ab_a - u_0*(nai_100_000_3 + nai_200_000_2*p_0); // local+recycle, weighs 15
  out[6] = 5.4181235997219*tmp58*(nai_300_010_0*v_1 - nai_100_000_3*u_1); // final, weighs 6
  nai_300_000_0 = -tmp64 + nai_200_000_0*v_1; // local+recycle, weighs 3
  tmp30 = -tmp62 + nai_200_000_1*v_1; // local+recycle, weighs 3
  nai_300_000_1 = nai_300_000_1/ab_a; // auto+recycle, weighs 2
  tmp31 = nai_200_000_3*u_1; // auto+recycle, weighs 1
  tmp35 = -tmp31 + nai_200_000_2*v_1; // local+recycle, weighs 3
  nai_300_010_1 = nai_300_010_1/ab_a; // auto+recycle, weighs 2
  tmp49 = d_2*u_1; // auto+recycle, weighs 1
  tmp51 = -tmp49 + d_1*v_1; // local+recycle, weighs 3
  tmp52 = nai_000_000_2*u_1; // auto+recycle, weighs 1
  tmp62 = -tmp52 + d_2*v_1; // local+recycle, weighs 3
  tmp64 = tmp51*v_0 - tmp62*u_0; // local+recycle, weighs 4
  tmp89 = nai_000_000_3*u_1; // auto, weighs 1
  nai_010_000_2 = -tmp89 + nai_000_000_2*v_1; // local, weighs 3
  nai_110_000_1 = tmp62*v_0 - nai_010_000_2*u_0; // local, weighs 4
  nai_100_000_0 = (0.5*nai_100_000_0 - 0.5*nai_100_000_1)/ab_a; // auto+recycle, weighs 5
  tmp86 = nai_000_000_4*u_1; // auto, weighs 1
  nai_010_000_3 = -tmp86 + nai_000_000_3*v_1; // local, weighs 3
  nai_110_000_2 = nai_010_000_2*v_0 - nai_010_000_3*u_0; // local, weighs 4
  nai_100_000_1 = (0.5*nai_100_000_1 - 0.5*nai_100_000_2)/ab_a; // auto+recycle, weighs 5
  out[7] = 9.3844653564063*tmp58*(p_0*(nai_300_000_1 + nai_300_000_0*p_1 - tmp30*u_1) + (nai_100_000_0 - nai_100_000_1 + nai_110_000_2*u_1 + p_1*tmp64 - nai_110_000_1*p_1 - nai_110_000_1*u_1)/ab_a - u_0*(nai_300_010_1 + p_1*tmp30 - tmp35*u_1)); // final, weighs 31
  nai_100_000_2 = nai_300_000_0*p_2 - tmp30*u_2; // local+recycle, weighs 4
  nai_300_000_0 = p_2*tmp30 - tmp35*u_2; // local+recycle, weighs 4
  tmp30 = (nai_110_000_2*u_2 + p_2*tmp64 - nai_110_000_1*p_2 - nai_110_000_1*u_2)/ab_a; // auto+recycle, weighs 11
  out[8] = 9.3844653564063*tmp58*(tmp30 + nai_100_000_2*p_0 - nai_300_000_0*u_0); // final, weighs 7
  tmp35 = nai_300_000_1 + nai_200_010_0*p_1 - nai_200_010_1*u_1; // local+recycle, weighs 5
  tmp31 = nai_300_010_1 + nai_200_010_1*p_1 - u_1*(-tmp31 + nai_200_000_2*p_1); // local+recycle, weighs 8
  out[9] = 5.4181235997219*tmp58*(tmp35*v_1 + (nai_200_010_0 - nai_200_010_1)/ab_a - tmp31*u_1); // final, weighs 11
  out[10] = 9.3844653564063*tmp58*(nai_100_000_2*p_1 + (0.5*nai_200_001_0 - 0.5*nai_200_001_1)/ab_a - nai_300_000_0*u_1); // final, weighs 12
  nai_100_000_2 = nai_300_000_1 + nai_200_001_0*p_2 - nai_200_001_1*u_2; // local+recycle, weighs 5
  nai_200_000_3 = nai_200_000_3*u_2; // auto+recycle, weighs 1
  nai_200_010_0 = nai_300_010_1 + nai_200_001_1*p_2 - u_2*(-nai_200_000_3 + nai_200_000_2*p_2); // local+recycle, weighs 8
  out[11] = 5.4181235997219*tmp58*(nai_100_000_2*v_1 - nai_200_010_0*u_1); // final, weighs 6
  out[12] = 5.4181235997219*tmp58*(nai_300_010_0*v_2 - nai_100_000_3*u_2); // final, weighs 6
  nai_100_000_3 = -nai_300_000_2 + nai_200_000_0*v_2; // local+recycle, weighs 3
  nai_200_000_0 = -tmp61 + nai_200_000_1*v_2; // local+recycle, weighs 3
  nai_200_000_1 = -nai_200_000_3 + nai_200_000_2*v_2; // local+recycle, weighs 3
  nai_200_000_2 = d_2*u_2; // auto+recycle, weighs 1
  nai_200_000_3 = -nai_200_000_2 + d_1*v_2; // local+recycle, weighs 3
  nai_200_010_1 = nai_000_000_2*u_2; // auto+recycle, weighs 1
  nai_300_000_0 = -nai_200_010_1 + d_2*v_2; // local+recycle, weighs 3
  nai_300_000_2 = nai_200_000_3*v_0 - nai_300_000_0*u_0; // local+recycle, weighs 4
  nai_300_010_0 = nai_000_000_3*u_2; // auto+recycle, weighs 1
  tmp61 = -nai_300_010_0 + nai_000_000_2*v_2; // local+recycle, weighs 3
  nai_101_000_1 = nai_300_000_0*v_0 - tmp61*u_0; // local, weighs 4
  tmp85 = nai_000_000_4*u_2; // auto, weighs 1
  nai_001_000_3 = -tmp85 + nai_000_000_3*v_2; // local, weighs 3
  nai_101_000_2 = tmp61*v_0 - nai_001_000_3*u_0; // local, weighs 4
  out[13] = 9.3844653564063*tmp58*(p_0*(nai_100_000_3*p_1 - nai_200_000_0*u_1) + (nai_101_000_2*u_1 + nai_300_000_2*p_1 - nai_101_000_1*p_1 - nai_101_000_1*u_1)/ab_a - u_0*(nai_200_000_0*p_1 - nai_200_000_1*u_1)); // final, weighs 26
  nai_100_000_3 = nai_300_000_1 + nai_100_000_3*p_2 - nai_200_000_0*u_2; // local+recycle, weighs 5
  nai_200_000_0 = nai_300_010_1 + nai_200_000_0*p_2 - nai_200_000_1*u_2; // local+recycle, weighs 5
  nai_100_000_0 = nai_100_000_0 + nai_300_000_2*p_2 - nai_101_000_1*u_2; // local+recycle, weighs 5
  nai_100_000_1 = nai_100_000_1 + nai_101_000_1*p_2 - nai_101_000_2*u_2; // local+recycle, weighs 5
  out[14] = 9.3844653564063*tmp58*(nai_100_000_3*p_0 + (nai_100_000_0 - nai_100_000_1)/ab_a - nai_200_000_0*u_0); // final, weighs 11
  out[15] = 5.4181235997219*tmp58*(tmp35*v_2 - tmp31*u_2); // final, weighs 6
  out[16] = 9.3844653564063*tmp58*(nai_100_000_3*p_1 - nai_200_000_0*u_1); // final, weighs 6
  out[17] = 5.4181235997219*tmp58*(nai_100_000_2*v_2 + (nai_200_001_0 - nai_200_001_1)/ab_a - nai_200_010_0*u_2); // final, weighs 11
  nai_100_000_2 = tmp2 - tmp62*u_1; // auto+recycle, weighs 3
  nai_100_000_3 = nai_100_000_2 + tmp51*v_1; // local+recycle, weighs 2
  nai_200_000_0 = tmp1 - nai_010_000_2*u_1; // auto+recycle, weighs 3
  nai_200_000_1 = nai_200_000_0 + tmp62*v_1; // local+recycle, weighs 2
  nai_200_001_0 = nai_200_000_1*u_0; // auto+recycle, weighs 1
  nai_200_001_1 = -nai_200_001_0 + nai_100_000_3*p_0; // local+recycle, weighs 3
  nai_200_010_0 = tmp0 - nai_010_000_3*u_1; // auto+recycle, weighs 3
  nai_300_000_1 = nai_200_010_0 + nai_010_000_2*v_1; // local+recycle, weighs 2
  nai_300_010_1 = nai_300_000_1*u_0; // auto+recycle, weighs 1
  tmp31 = -nai_300_010_1 + nai_200_000_1*p_0; // local+recycle, weighs 3
  tmp35 = 0.5*nai_100_000_3 - 0.5*nai_200_000_1; // auto+recycle, weighs 3
  tmp7 = tmp35/ab_a; // auto, weighs 2
  nai_020_200_0 = tmp7 + nai_200_001_1*p_0 - tmp31*u_0; // local, weighs 5
  tmp83 = d_0*u_1; // auto, weighs 1
  nai_020_000_3 = tmp3 + nai_010_000_3*v_1 - u_1*(-tmp83 + nai_000_000_4*v_1); // local, weighs 8
  tmp66 = nai_020_000_3*u_0; // auto, weighs 1
  tmp40 = 0.5*nai_200_000_1 - 0.5*nai_300_000_1; // auto, weighs 3
  tmp6 = tmp40/ab_a; // auto, weighs 2
  nai_020_200_1 = tmp6 + p_0*tmp31 - u_0*(-tmp66 + nai_300_000_1*p_0); // local, weighs 8
  out[18] = 5.4181235997219*tmp58*(nai_020_200_0*v_0 + (nai_200_001_1 - tmp31)/ab_a - nai_020_200_1*u_0); // final, weighs 11
  nai_200_001_0 = -nai_200_001_0 + nai_100_000_3*v_0; // local+recycle, weighs 3
  nai_200_001_1 = -nai_300_010_1 + nai_200_000_1*v_0; // local+recycle, weighs 3
  nai_300_010_1 = -tmp66 + nai_300_000_1*v_0; // local+recycle, weighs 3
  tmp31 = (tmp51 - tmp62)/ab_a - nai_200_000_1*u_1; // auto+recycle, weighs 7
  tmp66 = tmp31 + nai_100_000_3*p_1; // local+recycle, weighs 2
  tmp47 = (tmp62 - nai_010_000_2)/ab_a - nai_300_000_1*u_1; // auto, weighs 7
  nai_020_010_1 = tmp47 + nai_200_000_1*p_1; // local, weighs 2
  out[19] = 9.3844653564063*tmp58*(p_0*(nai_200_001_0*p_1 + (tmp64 - nai_110_000_1)/ab_a - nai_200_001_1*u_1) + (0.5*tmp66 - 0.5*nai_020_010_1)/ab_a - u_0*(nai_200_001_1*p_1 + (nai_110_000_1 - nai_110_000_2)/ab_a - nai_300_010_1*u_1)); // final, weighs 30
  nai_200_001_0 = nai_200_001_0*p_2 - nai_200_001_1*u_2; // local+recycle, weighs 4
  nai_200_001_1 = nai_200_001_1*p_2 - nai_300_010_1*u_2; // local+recycle, weighs 4
  nai_300_010_1 = nai_200_000_1*u_2; // auto+recycle, weighs 1
  nai_020_001_0 = -nai_300_010_1 + nai_100_000_3*p_2; // local, weighs 3
  tmp67 = nai_300_000_1*u_2; // auto, weighs 1
  nai_020_001_1 = -tmp67 + nai_200_000_1*p_2; // local, weighs 3
  out[20] = 9.3844653564063*tmp58*(nai_200_001_0*p_0 + (0.5*nai_020_001_0 - 0.5*nai_020_001_1)/ab_a - nai_200_001_1*u_0); // final, weighs 12
  nai_200_000_0 = nai_200_000_0 + p_1*tmp62; // local+recycle, weighs 2
  nai_100_000_2 = p_1*tmp66 + (nai_100_000_2 + tmp35 - nai_200_000_0 + p_1*tmp51)/ab_a - nai_020_010_1*u_1; // local+recycle, weighs 12
  nai_010_000_3 = (nai_010_000_2 - nai_010_000_3)/ab_a - nai_020_000_3*u_1; // auto+recycle, weighs 7
  nai_200_000_0 = nai_020_010_1*p_1 + (nai_200_000_0 + tmp40 - nai_200_010_0 - nai_010_000_2*p_1)/ab_a - u_1*(nai_010_000_3 + nai_300_000_1*p_1); // local+recycle, weighs 15
  out[21] = 5.4181235997219*tmp58*(nai_100_000_2*v_0 - nai_200_000_0*u_0); // final, weighs 6
  out[22] = 9.3844653564063*tmp58*(tmp30 + nai_200_001_0*p_1 - nai_200_001_1*u_1); // final, weighs 7
  nai_200_001_0 = tmp7 + nai_020_001_0*p_2 - nai_020_001_1*u_2; // local+recycle, weighs 5
  nai_020_000_3 = nai_020_000_3*u_2; // auto+recycle, weighs 1
  nai_200_001_1 = tmp6 + nai_020_001_1*p_2 - u_2*(-nai_020_000_3 + nai_300_000_1*p_2); // local+recycle, weighs 8
  out[23] = 5.4181235997219*tmp58*(nai_200_001_0*v_0 - nai_200_001_1*u_0); // final, weighs 6
  nai_200_010_0 = -tmp93 + d_2*p_0; // local+recycle, weighs 3
  tmp30 = -tmp90 + nai_000_000_2*p_0; // local+recycle, weighs 3
  tmp35 = tmp1 + nai_200_010_0*p_0 - tmp30*u_0; // local+recycle, weighs 5
  tmp40 = -tmp87 + nai_000_000_3*p_0; // local+recycle, weighs 3
  tmp30 = tmp0 + p_0*tmp30 - tmp40*u_0; // local+recycle, weighs 5
  tmp87 = tmp35*v_2 - tmp30*u_2; // local+recycle, weighs 4
  tmp90 = nai_300_000_0*u_1; // auto+recycle, weighs 1
  tmp93 = -tmp90 + nai_200_000_3*v_1; // local+recycle, weighs 3
  tmp80 = tmp61*u_1; // auto, weighs 1
  nai_011_000_1 = -tmp80 + nai_300_000_0*v_1; // local, weighs 3
  tmp72 = nai_011_000_1*u_0; // auto, weighs 1
  tmp79 = nai_001_000_3*u_1; // auto, weighs 1
  nai_011_000_2 = -tmp79 + tmp61*v_1; // local, weighs 3
  tmp71 = nai_011_000_2*u_0; // auto, weighs 1
  out[24] = 9.3844653564063*tmp58*(v_0*(v_1*(v_2*(tmp2 + p_0*(-tmp96 + d_1*p_0) - nai_200_010_0*u_0) - tmp35*u_2) - tmp87*u_1) + (tmp71 - tmp72 + p_0*tmp93 - nai_011_000_1*p_0)/ab_a - u_0*(tmp87*v_1 - u_1*(tmp30*v_2 - u_2*(tmp3 + p_0*tmp40 - u_0*(-usq + nai_000_000_4*p_0))))); // final, weighs 48
  nai_200_010_0 = -tmp72 + tmp93*v_0; // local+recycle, weighs 3
  tmp30 = -tmp71 + nai_011_000_1*v_0; // local+recycle, weighs 3
  d_0 = d_0*u_2; // auto+recycle, weighs 1
  tmp35 = -d_0 + nai_000_000_4*v_2; // local+recycle, weighs 3
  tmp40 = nai_011_000_2*v_0 - u_0*(nai_001_000_3*v_1 - tmp35*u_1); // local+recycle, weighs 8
  out[25] = 16.2543707991657*tmp58*(p_0*(nai_200_010_0*p_1 + (0.5*nai_300_000_2 - 0.5*nai_101_000_1)/ab_a - tmp30*u_1) + (0.5*nai_011_000_2*u_1 + 0.5*p_1*tmp93 + 0.5*(0.5*nai_200_000_3 - 0.5*nai_300_000_0)/ab_a - 0.5*nai_011_000_1*p_1 - 0.5*nai_011_000_1*u_1 - 0.5*(0.5*nai_300_000_0 - 0.5*tmp61)/ab_a)/ab_a - u_0*(p_1*tmp30 + (0.5*nai_101_000_1 - 0.5*nai_101_000_2)/ab_a - tmp40*u_1)); // final, weighs 54
  nai_200_010_0 = nai_200_010_0*p_2 + (0.5*tmp64 - 0.5*nai_110_000_1)/ab_a - tmp30*u_2; // local+recycle, weighs 10
  nai_110_000_1 = p_2*tmp30 + (0.5*nai_110_000_1 - 0.5*nai_110_000_2)/ab_a - tmp40*u_2; // local+recycle, weighs 10
  nai_110_000_2 = p_2*tmp93 + (0.5*tmp51 - 0.5*tmp62)/ab_a - nai_011_000_1*u_2; // local+recycle, weighs 10
  nai_010_000_2 = nai_011_000_1*p_2 + (0.5*tmp62 - 0.5*nai_010_000_2)/ab_a - nai_011_000_2*u_2; // local+recycle, weighs 10
  out[26] = 16.2543707991657*tmp58*(nai_200_010_0*p_0 + (0.5*nai_110_000_2 - 0.5*nai_010_000_2)/ab_a - nai_110_000_1*u_0); // final, weighs 12
  tmp30 = -tmp52 + d_2*p_1; // local+recycle, weighs 3
  tmp40 = -tmp89 + nai_000_000_2*p_1; // local+recycle, weighs 3
  tmp51 = tmp1 + p_1*tmp30 - tmp40*u_1; // local+recycle, weighs 5
  tmp52 = -tmp86 + nai_000_000_3*p_1; // local+recycle, weighs 3
  tmp40 = tmp0 + p_1*tmp40 - tmp52*u_1; // local+recycle, weighs 5
  tmp62 = tmp51*v_2 - tmp40*u_2; // local+recycle, weighs 4
  tmp64 = -tmp80 + nai_300_000_0*p_1; // local+recycle, weighs 3
  out[27] = 9.3844653564063*tmp58*(v_0*(v_1*(v_2*(tmp2 + p_1*(-tmp49 + d_1*p_1) - tmp30*u_1) - tmp51*u_2) + (-tmp64 - tmp90 + nai_200_000_3*p_1)/ab_a - tmp62*u_1) - u_0*(tmp62*v_1 + (tmp64 + tmp79 - p_1*tmp61)/ab_a - u_1*(tmp40*v_2 - u_2*(tmp3 + p_1*tmp52 - u_1*(-tmp83 + nai_000_000_4*p_1))))); // final, weighs 53
  out[28] = 16.2543707991657*tmp58*(nai_200_010_0*p_1 + (0.5*nai_100_000_0 - 0.5*nai_100_000_1)/ab_a - nai_110_000_1*u_1); // final, weighs 12
  d_1 = -nai_200_000_2 + d_1*p_2; // local+recycle, weighs 3
  d_2 = -nai_200_010_1 + d_2*p_2; // local+recycle, weighs 3
  nai_000_000_2 = -nai_300_010_0 + nai_000_000_2*p_2; // local+recycle, weighs 3
  nai_100_000_0 = tmp1 + d_2*p_2 - nai_000_000_2*u_2; // local+recycle, weighs 5
  nai_000_000_3 = -tmp85 + nai_000_000_3*p_2; // local+recycle, weighs 3
  nai_100_000_1 = tmp0 + nai_000_000_2*p_2 - nai_000_000_3*u_2; // local+recycle, weighs 5
  nai_110_000_1 = nai_100_000_0*v_2 + (d_2 - nai_000_000_2)/ab_a - nai_100_000_1*u_2; // local+recycle, weighs 9
  out[29] = 9.3844653564063*tmp58*(v_0*(v_1*(v_2*(tmp2 + d_1*p_2 - d_2*u_2) + (d_1 - d_2)/ab_a - nai_100_000_0*u_2) - nai_110_000_1*u_1) - u_0*(nai_110_000_1*v_1 - u_1*(nai_100_000_1*v_2 + (nai_000_000_2 - nai_000_000_3)/ab_a - u_2*(tmp3 + nai_000_000_3*p_2 - u_2*(-d_0 + nai_000_000_4*p_2))))); // final, weighs 45
  d_0 = tmp2 - nai_300_000_0*u_2; // auto+recycle, weighs 3
  d_1 = d_0 + nai_200_000_3*v_2; // local+recycle, weighs 2
  d_2 = tmp1 - tmp61*u_2; // auto+recycle, weighs 3
  nai_000_000_2 = d_2 + nai_300_000_0*v_2; // local+recycle, weighs 2
  nai_000_000_3 = nai_000_000_2*u_0; // auto+recycle, weighs 1
  nai_000_000_4 = -nai_000_000_3 + d_1*p_0; // local+recycle, weighs 3
  nai_100_000_0 = tmp0 - nai_001_000_3*u_2; // auto+recycle, weighs 3
  nai_100_000_1 = nai_100_000_0 + tmp61*v_2; // local+recycle, weighs 2
  nai_110_000_1 = nai_100_000_1*u_0; // auto+recycle, weighs 1
  nai_200_000_2 = -nai_110_000_1 + nai_000_000_2*p_0; // local+recycle, weighs 3
  nai_200_010_0 = 0.5*d_1 - 0.5*nai_000_000_2; // auto+recycle, weighs 3
  nai_200_010_1 = nai_200_010_0/ab_a; // auto+recycle, weighs 2
  nai_300_010_0 = nai_200_010_1 + nai_000_000_4*p_0 - nai_200_000_2*u_0; // local+recycle, weighs 5
  tmp0 = tmp3 + nai_001_000_3*v_2 - tmp35*u_2; // local+recycle, weighs 5
  tmp1 = tmp0*u_0; // auto+recycle, weighs 1
  tmp2 = 0.5*nai_000_000_2 - 0.5*nai_100_000_1; // auto+recycle, weighs 3
  tmp3 = tmp2/ab_a; // auto+recycle, weighs 2
  tmp30 = tmp3 + nai_200_000_2*p_0 - u_0*(-tmp1 + nai_100_000_1*p_0); // local+recycle, weighs 8
  out[30] = 5.4181235997219*tmp58*(nai_300_010_0*v_0 + (nai_000_000_4 - nai_200_000_2)/ab_a - tmp30*u_0); // final, weighs 11
  nai_000_000_3 = -nai_000_000_3 + d_1*v_0; // local+recycle, weighs 3
  nai_000_000_4 = -nai_110_000_1 + nai_000_000_2*v_0; // local+recycle, weighs 3
  nai_110_000_1 = -tmp1 + nai_100_000_1*v_0; // local+recycle, weighs 3
  nai_200_000_2 = nai_000_000_2*u_1; // auto+recycle, weighs 1
  tmp1 = -nai_200_000_2 + d_1*p_1; // local+recycle, weighs 3
  tmp35 = nai_100_000_1*u_1; // auto+recycle, weighs 1
  tmp40 = -tmp35 + nai_000_000_2*p_1; // local+recycle, weighs 3
  out[31] = 9.3844653564063*tmp58*(p_0*(nai_000_000_3*p_1 - nai_000_000_4*u_1) + (0.5*tmp1 - 0.5*tmp40)/ab_a - u_0*(nai_000_000_4*p_1 - nai_110_000_1*u_1)); // final, weighs 20
  nai_000_000_3 = nai_000_000_3*p_2 + (nai_300_000_2 - nai_101_000_1)/ab_a - nai_000_000_4*u_2; // local+recycle, weighs 9
  nai_000_000_4 = nai_000_000_4*p_2 + (nai_101_000_1 - nai_101_000_2)/ab_a - nai_110_000_1*u_2; // local+recycle, weighs 9
  nai_101_000_1 = (nai_200_000_3 - nai_300_000_0)/ab_a - nai_000_000_2*u_2; // auto+recycle, weighs 7
  nai_101_000_2 = nai_101_000_1 + d_1*p_2; // local+recycle, weighs 2
  nai_110_000_1 = (nai_300_000_0 - tmp61)/ab_a - nai_100_000_1*u_2; // auto+recycle, weighs 7
  nai_300_000_2 = nai_110_000_1 + nai_000_000_2*p_2; // local+recycle, weighs 2
  tmp49 = (0.5*nai_101_000_2 - 0.5*nai_300_000_2)/ab_a; // auto+recycle, weighs 5
  out[32] = 9.3844653564063*tmp58*(tmp49 + nai_000_000_3*p_0 - nai_000_000_4*u_0); // final, weighs 7
  tmp51 = nai_200_010_1 + p_1*tmp1 - tmp40*u_1; // local+recycle, weighs 5
  tmp52 = tmp0*u_1; // auto+recycle, weighs 1
  tmp62 = tmp3 + p_1*tmp40 - u_1*(-tmp52 + nai_100_000_1*p_1); // local+recycle, weighs 8
  out[33] = 5.4181235997219*tmp58*(tmp51*v_0 - tmp62*u_0); // final, weighs 6
  out[34] = 9.3844653564063*tmp58*(nai_000_000_3*p_1 - nai_000_000_4*u_1); // final, weighs 6
  d_2 = d_2 + nai_300_000_0*p_2; // local+recycle, weighs 2
  d_0 = nai_101_000_2*p_2 + (d_0 + nai_200_010_0 - d_2 + nai_200_000_3*p_2)/ab_a - nai_300_000_2*u_2; // local+recycle, weighs 12
  nai_000_000_3 = (tmp61 - nai_001_000_3)/ab_a - tmp0*u_2; // auto+recycle, weighs 7
  d_2 = nai_300_000_2*p_2 + (d_2 + tmp2 - nai_100_000_0 - p_2*tmp61)/ab_a - u_2*(nai_000_000_3 + nai_100_000_1*p_2); // local+recycle, weighs 15
  out[35] = 5.4181235997219*tmp58*(d_0*v_0 - d_2*u_0); // final, weighs 6
  nai_000_000_4 = tmp31 + nai_100_000_3*v_1; // local+recycle, weighs 2
  nai_001_000_3 = tmp47 + nai_200_000_1*v_1; // local+recycle, weighs 2
  nai_010_000_3 = nai_010_000_3 + nai_300_000_1*v_1; // local+recycle, weighs 2
  nai_100_000_0 = 0.5*nai_000_000_4 - 0.5*nai_001_000_3; // auto+recycle, weighs 3
  nai_200_000_3 = nai_100_000_0/ab_a; // auto+recycle, weighs 2
  out[36] = 2.4230585358948*tmp58*(nai_200_000_3 + p_0*(nai_000_000_4*p_0 - nai_001_000_3*u_0) - u_0*(nai_001_000_3*p_0 - nai_010_000_3*u_0)); // final, weighs 15
  nai_200_010_0 = nai_000_000_4*p_1 + (1.5*nai_100_000_3 - 1.5*nai_200_000_1)/ab_a - nai_001_000_3*u_1; // local+recycle, weighs 10
  nai_300_000_0 = nai_001_000_3*p_1 + (1.5*nai_200_000_1 - 1.5*nai_300_000_1)/ab_a - nai_010_000_3*u_1; // local+recycle, weighs 10
  out[37] = 4.19686049388326*tmp58*(nai_200_010_0*p_0 - nai_300_000_0*u_0); // final, weighs 6
  nai_000_000_4 = nai_000_000_4*p_2 - nai_001_000_3*u_2; // local+recycle, weighs 4
  nai_001_000_3 = nai_001_000_3*p_2 - nai_010_000_3*u_2; // local+recycle, weighs 4
  out[38] = 4.19686049388326*tmp58*(nai_000_000_4*p_0 - nai_001_000_3*u_0); // final, weighs 6
  out[39] = 2.4230585358948*tmp58*(nai_200_010_0*p_1 + (nai_100_000_0 + 1.5*tmp66 - 1.5*nai_020_010_1)/ab_a - nai_300_000_0*u_1); // final, weighs 13
  out[40] = 4.19686049388326*tmp58*(nai_000_000_4*p_1 + (1.5*nai_020_001_0 - 1.5*nai_020_001_1)/ab_a - nai_001_000_3*u_1); // final, weighs 12
  out[41] = 2.4230585358948*tmp58*(nai_200_000_3 + nai_000_000_4*p_2 - nai_001_000_3*u_2); // final, weighs 7
  out[42] = 5.4181235997219*tmp58*(nai_020_200_0*v_2 - nai_020_200_1*u_2); // final, weighs 6
  nai_000_000_4 = -nai_300_010_1 + nai_100_000_3*v_2; // local+recycle, weighs 3
  nai_001_000_3 = -tmp67 + nai_200_000_1*v_2; // local+recycle, weighs 3
  nai_010_000_3 = (tmp93 - nai_011_000_1)/ab_a; // auto+recycle, weighs 4
  nai_020_000_3 = -nai_020_000_3 + nai_300_000_1*v_2; // local+recycle, weighs 3
  nai_011_000_1 = (nai_011_000_1 - nai_011_000_2)/ab_a; // auto+recycle, weighs 4
  out[43] = 9.3844653564063*tmp58*(p_0*(nai_010_000_3 + nai_000_000_4*p_1 - nai_001_000_3*u_1) - u_0*(nai_011_000_1 + nai_001_000_3*p_1 - nai_020_000_3*u_1)); // final, weighs 16
  nai_000_000_4 = tmp7 + nai_000_000_4*p_2 - nai_001_000_3*u_2; // local+recycle, weighs 5
  nai_001_000_3 = tmp6 + nai_001_000_3*p_2 - nai_020_000_3*u_2; // local+recycle, weighs 5
  out[44] = 9.3844653564063*tmp58*(nai_000_000_4*p_0 - nai_001_000_3*u_0); // final, weighs 6
  out[45] = 5.4181235997219*tmp58*(nai_100_000_2*v_2 - nai_200_000_0*u_2); // final, weighs 6
  out[46] = 9.3844653564063*tmp58*(nai_000_000_4*p_1 + (nai_110_000_2 - nai_010_000_2)/ab_a - nai_001_000_3*u_1); // final, weighs 11
  out[47] = 5.4181235997219*tmp58*(nai_200_001_0*v_2 + (nai_020_001_0 - nai_020_001_1)/ab_a - nai_200_001_1*u_2); // final, weighs 11
  out[48] = 5.4181235997219*tmp58*(nai_300_010_0*v_1 - tmp30*u_1); // final, weighs 6
  nai_000_000_4 = -nai_200_000_2 + d_1*v_1; // local+recycle, weighs 3
  nai_001_000_3 = -tmp35 + nai_000_000_2*v_1; // local+recycle, weighs 3
  nai_010_000_2 = -tmp52 + nai_100_000_1*v_1; // local+recycle, weighs 3
  out[49] = 9.3844653564063*tmp58*(p_0*(nai_200_010_1 + nai_000_000_4*p_1 - nai_001_000_3*u_1) - u_0*(tmp3 + nai_001_000_3*p_1 - nai_010_000_2*u_1)); // final, weighs 16
  nai_000_000_4 = nai_010_000_3 + nai_000_000_4*p_2 - nai_001_000_3*u_2; // local+recycle, weighs 5
  nai_001_000_3 = nai_011_000_1 + nai_001_000_3*p_2 - nai_010_000_2*u_2; // local+recycle, weighs 5
  out[50] = 9.3844653564063*tmp58*(nai_000_000_4*p_0 - nai_001_000_3*u_0); // final, weighs 6
  out[51] = 5.4181235997219*tmp58*(tmp51*v_1 + (tmp1 - tmp40)/ab_a - tmp62*u_1); // final, weighs 11
  out[52] = 9.3844653564063*tmp58*(tmp49 + nai_000_000_4*p_1 - nai_001_000_3*u_1); // final, weighs 7
  out[53] = 5.4181235997219*tmp58*(d_0*v_1 - d_2*u_1); // final, weighs 6
  d_0 = nai_101_000_1 + d_1*v_2; // local+recycle, weighs 2
  d_2 = nai_110_000_1 + nai_000_000_2*v_2; // local+recycle, weighs 2
  nai_000_000_3 = nai_000_000_3 + nai_100_000_1*v_2; // local+recycle, weighs 2
  nai_000_000_4 = 0.5*d_0 - 0.5*d_2; // auto+recycle, weighs 3
  nai_001_000_3 = nai_000_000_4/ab_a; // auto+recycle, weighs 2
  out[54] = 2.4230585358948*tmp58*(nai_001_000_3 + p_0*(d_0*p_0 - d_2*u_0) - u_0*(d_2*p_0 - nai_000_000_3*u_0)); // final, weighs 15
  nai_010_000_2 = d_0*p_1 - d_2*u_1; // local+recycle, weighs 4
  nai_010_000_3 = d_2*p_1 - nai_000_000_3*u_1; // local+recycle, weighs 4
  out[55] = 4.19686049388326*tmp58*(nai_010_000_2*p_0 - nai_010_000_3*u_0); // final, weighs 6
  d_0 = d_0*p_2 + (1.5*d_1 - 1.5*nai_000_000_2)/ab_a - d_2*u_2; // local+recycle, weighs 10
  d_1 = d_2*p_2 + (1.5*nai_000_000_2 - 1.5*nai_100_000_1)/ab_a - nai_000_000_3*u_2; // local+recycle, weighs 10
  out[56] = 4.19686049388326*tmp58*(d_0*p_0 - d_1*u_0); // final, weighs 6
  out[57] = 2.4230585358948*tmp58*(nai_001_000_3 + nai_010_000_2*p_1 - nai_010_000_3*u_1); // final, weighs 7
  out[58] = 4.19686049388326*tmp58*(d_0*p_1 - d_1*u_1); // final, weighs 6
  out[59] = 2.4230585358948*tmp58*(d_0*p_2 + (nai_000_000_4 + 1.5*nai_101_000_2 - 1.5*nai_300_000_2)/ab_a - d_1*u_2); // final, weighs 13
  // total weight = 1771
}

static void gint2_nai_pF_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_pD_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // total weight = 0
}

static void gint2_nai_SP_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 60
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_000_001_0, nai_000_001_1, nai_000_001_2, nai_000_002_0, nai_000_002_1, nai_000_003_0, nai_000_010_0, nai_000_010_1, nai_000_010_2, nai_000_011_0, nai_000_011_1, nai_000_020_0, nai_000_020_1, nai_000_030_0, nai_000_100_0, nai_000_100_1, nai_000_100_2, nai_000_200_0, nai_000_200_1, nai_000_200_2, nai_000_300_0, nai_100_001_1, nai_100_020_1, p_0, p_1, p_2, tmp1, tmp12, tmp19, tmp2, tmp23, tmp24, tmp25, tmp26, tmp27, tmp28, tmp29, tmp33, tmp34, tmp35, tmp36, tmp37, tmp38, tmp39, tmp40, tmp41, tmp5, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  tmp41 = d_2*u_0; // auto, weighs 1
  nai_000_100_0 = -tmp41 + d_1*p_0; // local, weighs 3
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  tmp38 = nai_000_000_2*u_0; // auto, weighs 1
  nai_000_100_1 = -tmp38 + d_2*p_0; // local, weighs 3
  tmp2 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  nai_000_200_0 = tmp2 + nai_000_100_0*p_0 - nai_000_100_1*u_0; // local, weighs 5
  nai_000_000_3 = 6.28318530717959*d_0*gaux(usq, 3); // local, weighs 3
  tmp35 = nai_000_000_3*u_0; // auto, weighs 1
  nai_000_100_2 = -tmp35 + nai_000_000_2*p_0; // local, weighs 3
  tmp1 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto, weighs 5
  nai_000_200_1 = tmp1 + nai_000_100_1*p_0 - nai_000_100_2*u_0; // local, weighs 5
  nai_000_100_0 = (nai_000_100_0 - nai_000_100_1)/ab_a - nai_000_200_1*u_0; // auto+recycle, weighs 7
  nai_000_300_0 = nai_000_100_0 + nai_000_200_0*p_0; // local, weighs 2
  tmp12 = pow(b_a,2.25); // auto, weighs 1
  tmp23 = tmp12*pow(a_a,0.75); // auto, weighs 2
  out[0] = 1.04921512347081*nai_000_300_0*tmp23; // final, weighs 2
  tmp25 = nai_000_200_1*u_1; // auto, weighs 1
  out[1] = 2.34611633910158*tmp23*(-tmp25 + nai_000_200_0*p_1); // final, weighs 5
  tmp24 = nai_000_200_1*u_2; // auto, weighs 1
  out[2] = 2.34611633910158*tmp23*(-tmp24 + nai_000_200_0*p_2); // final, weighs 5
  tmp40 = d_2*u_1; // auto, weighs 1
  nai_000_010_0 = -tmp40 + d_1*p_1; // local, weighs 3
  tmp37 = nai_000_000_2*u_1; // auto, weighs 1
  nai_000_010_1 = -tmp37 + d_2*p_1; // local, weighs 3
  nai_000_020_0 = tmp2 + nai_000_010_0*p_1 - nai_000_010_1*u_1; // local, weighs 5
  tmp34 = nai_000_000_3*u_1; // auto, weighs 1
  nai_000_010_2 = -tmp34 + nai_000_000_2*p_1; // local, weighs 3
  nai_000_020_1 = tmp1 + nai_000_010_1*p_1 - nai_000_010_2*u_1; // local, weighs 5
  tmp27 = nai_000_020_1*u_0; // auto, weighs 1
  out[3] = 2.34611633910158*tmp23*(-tmp27 + nai_000_020_0*p_0); // final, weighs 5
  tmp39 = d_2*u_2; // auto, weighs 1
  nai_000_001_0 = -tmp39 + d_1*p_2; // local, weighs 3
  tmp36 = nai_000_000_2*u_2; // auto, weighs 1
  nai_000_001_1 = -tmp36 + d_2*p_2; // local, weighs 3
  nai_000_011_0 = nai_000_001_0*p_1 - nai_000_001_1*u_1; // local, weighs 4
  tmp33 = nai_000_000_3*u_2; // auto, weighs 1
  nai_000_001_2 = -tmp33 + nai_000_000_2*p_2; // local, weighs 3
  nai_000_011_1 = nai_000_001_1*p_1 - nai_000_001_2*u_1; // local, weighs 4
  out[4] = 4.06359269979142*tmp23*(nai_000_011_0*p_0 - nai_000_011_1*u_0); // final, weighs 6
  nai_000_002_0 = tmp2 + nai_000_001_0*p_2 - nai_000_001_1*u_2; // local, weighs 5
  nai_000_002_1 = tmp1 + nai_000_001_1*p_2 - nai_000_001_2*u_2; // local, weighs 5
  tmp29 = nai_000_002_1*u_0; // auto, weighs 1
  out[5] = 2.34611633910158*tmp23*(-tmp29 + nai_000_002_0*p_0); // final, weighs 5
  nai_000_010_0 = (nai_000_010_0 - nai_000_010_1)/ab_a - nai_000_020_1*u_1; // auto+recycle, weighs 7
  nai_000_030_0 = nai_000_010_0 + nai_000_020_0*p_1; // local, weighs 2
  out[6] = 1.04921512347081*nai_000_030_0*tmp23; // final, weighs 2
  tmp26 = nai_000_020_1*u_2; // auto, weighs 1
  out[7] = 2.34611633910158*tmp23*(-tmp26 + nai_000_020_0*p_2); // final, weighs 5
  tmp28 = nai_000_002_1*u_1; // auto, weighs 1
  out[8] = 2.34611633910158*tmp23*(-tmp28 + nai_000_002_0*p_1); // final, weighs 5
  tmp19 = (nai_000_001_0 - nai_000_001_1)/ab_a - nai_000_002_1*u_2; // auto, weighs 7
  nai_000_003_0 = tmp19 + nai_000_002_0*p_2; // local, weighs 2
  out[9] = 1.04921512347081*nai_000_003_0*tmp23; // final, weighs 2
  d_0 = 6.28318530717959*d_0*gaux(usq, 4); // local+recycle, weighs 3
  tmp23 = d_0*u_0; // auto+recycle, weighs 1
  usq = (0.5*nai_000_000_2 - 0.5*nai_000_000_3)/ab_a; // auto+recycle, weighs 5
  nai_000_200_2 = usq + nai_000_100_2*p_0 - u_0*(-tmp23 + nai_000_000_3*p_0); // local, weighs 8
  nai_000_100_1 = (nai_000_100_1 - nai_000_100_2)/ab_a - nai_000_200_2*u_0; // auto+recycle, weighs 7
  nai_000_100_2 = nai_000_100_1 + nai_000_200_1*p_0; // local+recycle, weighs 2
  tmp12 = tmp12*pow(a_a,1.25); // auto+recycle, weighs 2
  out[10] = 2.09843024694163*tmp12*(nai_000_300_0*v_0 + (1.5*nai_000_200_0 - 1.5*nai_000_200_1)/ab_a - nai_000_100_2*u_0); // final, weighs 12
  nai_000_100_0 = nai_000_100_0 + nai_000_200_0*v_0; // local+recycle, weighs 2
  nai_000_100_1 = nai_000_100_1 + nai_000_200_1*v_0; // local+recycle, weighs 2
  out[11] = 4.69223267820315*tmp12*(nai_000_100_0*p_1 - nai_000_100_1*u_1); // final, weighs 6
  out[12] = 4.69223267820315*tmp12*(nai_000_100_0*p_2 - nai_000_100_1*u_2); // final, weighs 6
  nai_000_100_0 = -tmp27 + nai_000_020_0*v_0; // local+recycle, weighs 3
  nai_000_100_1 = d_0*u_1; // auto+recycle, weighs 1
  tmp27 = usq + nai_000_010_2*p_1 - u_1*(-nai_000_100_1 + nai_000_000_3*p_1); // local+recycle, weighs 8
  nai_100_020_1 = nai_000_020_1*v_0 - tmp27*u_0; // local, weighs 4
  tmp5 = (0.5*nai_000_020_0 - 0.5*nai_000_020_1)/ab_a; // auto, weighs 5
  out[13] = 4.69223267820315*tmp12*(tmp5 + nai_000_100_0*p_0 - nai_100_020_1*u_0); // final, weighs 7
  tmp38 = -tmp38 + d_2*v_0; // local+recycle, weighs 3
  tmp35 = -tmp35 + nai_000_000_2*v_0; // local+recycle, weighs 3
  nai_100_001_1 = p_2*tmp38 - tmp35*u_2; // local, weighs 4
  out[14] = 8.12718539958285*tmp12*(p_0*(p_1*(p_2*(-tmp41 + d_1*v_0) - tmp38*u_2) - nai_100_001_1*u_1) + (0.5*nai_000_011_0 - 0.5*nai_000_011_1)/ab_a - u_0*(nai_100_001_1*p_1 - u_1*(p_2*tmp35 - u_2*(-tmp23 + nai_000_000_3*v_0)))); // final, weighs 34
  nai_000_011_0 = -tmp29 + nai_000_002_0*v_0; // local+recycle, weighs 3
  d_0 = d_0*u_2; // auto+recycle, weighs 1
  nai_000_011_1 = usq + nai_000_001_2*p_2 - u_2*(-d_0 + nai_000_000_3*p_2); // local+recycle, weighs 8
  nai_100_001_1 = nai_000_002_1*v_0 - nai_000_011_1*u_0; // local+recycle, weighs 4
  tmp23 = (0.5*nai_000_002_0 - 0.5*nai_000_002_1)/ab_a; // auto+recycle, weighs 5
  out[15] = 4.69223267820315*tmp12*(tmp23 + nai_000_011_0*p_0 - nai_100_001_1*u_0); // final, weighs 7
  nai_000_010_1 = (nai_000_010_1 - nai_000_010_2)/ab_a - tmp27*u_1; // auto+recycle, weighs 7
  nai_000_010_2 = nai_000_010_1 + nai_000_020_1*p_1; // local+recycle, weighs 2
  out[16] = 2.09843024694163*tmp12*(nai_000_030_0*v_0 - nai_000_010_2*u_0); // final, weighs 6
  out[17] = 4.69223267820315*tmp12*(nai_000_100_0*p_2 - nai_100_020_1*u_2); // final, weighs 6
  out[18] = 4.69223267820315*tmp12*(nai_000_011_0*p_1 - nai_100_001_1*u_1); // final, weighs 6
  nai_000_011_0 = (nai_000_001_1 - nai_000_001_2)/ab_a - nai_000_011_1*u_2; // auto+recycle, weighs 7
  nai_000_100_0 = nai_000_011_0 + nai_000_002_1*p_2; // local+recycle, weighs 2
  out[19] = 2.09843024694163*tmp12*(nai_000_003_0*v_0 - nai_000_100_0*u_0); // final, weighs 6
  out[20] = 2.09843024694163*tmp12*(nai_000_300_0*v_1 - nai_000_100_2*u_1); // final, weighs 6
  nai_100_001_1 = -tmp25 + nai_000_200_0*v_1; // local+recycle, weighs 3
  nai_100_020_1 = nai_000_200_1*v_1 - nai_000_200_2*u_1; // local+recycle, weighs 4
  tmp25 = (0.5*nai_000_200_0 - 0.5*nai_000_200_1)/ab_a; // auto+recycle, weighs 5
  out[21] = 4.69223267820315*tmp12*(tmp25 + nai_100_001_1*p_1 - nai_100_020_1*u_1); // final, weighs 7
  out[22] = 4.69223267820315*tmp12*(nai_100_001_1*p_2 - nai_100_020_1*u_2); // final, weighs 6
  nai_000_010_0 = nai_000_010_0 + nai_000_020_0*v_1; // local+recycle, weighs 2
  nai_000_010_1 = nai_000_010_1 + nai_000_020_1*v_1; // local+recycle, weighs 2
  out[23] = 4.69223267820315*tmp12*(nai_000_010_0*p_0 - nai_000_010_1*u_0); // final, weighs 6
  nai_100_001_1 = -tmp37 + d_2*v_1; // local+recycle, weighs 3
  nai_100_020_1 = -tmp34 + nai_000_000_2*v_1; // local+recycle, weighs 3
  tmp29 = nai_100_001_1*p_2 - nai_100_020_1*u_2; // local+recycle, weighs 4
  out[24] = 8.12718539958285*tmp12*(p_0*(p_1*(p_2*(-tmp40 + d_1*v_1) - nai_100_001_1*u_2) + (0.5*nai_000_001_0 - 0.5*nai_000_001_1)/ab_a - tmp29*u_1) - u_0*(p_1*tmp29 + (0.5*nai_000_001_1 - 0.5*nai_000_001_2)/ab_a - u_1*(nai_100_020_1*p_2 - u_2*(-nai_000_100_1 + nai_000_000_3*v_1)))); // final, weighs 40
  nai_000_001_0 = -tmp28 + nai_000_002_0*v_1; // local+recycle, weighs 3
  nai_000_001_1 = nai_000_002_1*v_1 - nai_000_011_1*u_1; // local+recycle, weighs 4
  out[25] = 4.69223267820315*tmp12*(nai_000_001_0*p_0 - nai_000_001_1*u_0); // final, weighs 6
  out[26] = 2.09843024694163*tmp12*(nai_000_030_0*v_1 + (1.5*nai_000_020_0 - 1.5*nai_000_020_1)/ab_a - nai_000_010_2*u_1); // final, weighs 12
  out[27] = 4.69223267820315*tmp12*(nai_000_010_0*p_2 - nai_000_010_1*u_2); // final, weighs 6
  out[28] = 4.69223267820315*tmp12*(tmp23 + nai_000_001_0*p_1 - nai_000_001_1*u_1); // final, weighs 7
  out[29] = 2.09843024694163*tmp12*(nai_000_003_0*v_1 - nai_000_100_0*u_1); // final, weighs 6
  out[30] = 2.09843024694163*tmp12*(nai_000_300_0*v_2 - nai_000_100_2*u_2); // final, weighs 6
  nai_000_001_0 = -tmp24 + nai_000_200_0*v_2; // local+recycle, weighs 3
  nai_000_001_1 = nai_000_200_1*v_2 - nai_000_200_2*u_2; // local+recycle, weighs 4
  out[31] = 4.69223267820315*tmp12*(nai_000_001_0*p_1 - nai_000_001_1*u_1); // final, weighs 6
  out[32] = 4.69223267820315*tmp12*(tmp25 + nai_000_001_0*p_2 - nai_000_001_1*u_2); // final, weighs 7
  nai_000_001_0 = -tmp26 + nai_000_020_0*v_2; // local+recycle, weighs 3
  nai_000_001_1 = nai_000_020_1*v_2 - tmp27*u_2; // local+recycle, weighs 4
  out[33] = 4.69223267820315*tmp12*(nai_000_001_0*p_0 - nai_000_001_1*u_0); // final, weighs 6
  d_2 = -tmp36 + d_2*v_2; // local+recycle, weighs 3
  nai_000_000_2 = -tmp33 + nai_000_000_2*v_2; // local+recycle, weighs 3
  nai_000_001_2 = tmp1 + d_2*p_2 - nai_000_000_2*u_2; // local+recycle, weighs 5
  out[34] = 8.12718539958285*tmp12*(p_0*(p_1*(tmp2 + p_2*(-tmp39 + d_1*v_2) - d_2*u_2) - nai_000_001_2*u_1) - u_0*(nai_000_001_2*p_1 - u_1*(usq + nai_000_000_2*p_2 - u_2*(-d_0 + nai_000_000_3*v_2)))); // final, weighs 30
  d_0 = tmp19 + nai_000_002_0*v_2; // local+recycle, weighs 2
  d_1 = nai_000_011_0 + nai_000_002_1*v_2; // local+recycle, weighs 2
  out[35] = 4.69223267820315*tmp12*(d_0*p_0 - d_1*u_0); // final, weighs 6
  out[36] = 2.09843024694163*tmp12*(nai_000_030_0*v_2 - nai_000_010_2*u_2); // final, weighs 6
  out[37] = 4.69223267820315*tmp12*(tmp5 + nai_000_001_0*p_2 - nai_000_001_1*u_2); // final, weighs 7
  out[38] = 4.69223267820315*tmp12*(d_0*p_1 - d_1*u_1); // final, weighs 6
  out[39] = 2.09843024694163*tmp12*(nai_000_003_0*v_2 + (1.5*nai_000_002_0 - 1.5*nai_000_002_1)/ab_a - nai_000_100_0*u_2); // final, weighs 12
  // total weight = 687
}

static void gint2_nai_S_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 18
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_100_0, nai_000_100_1, nai_000_200_0, nai_000_200_1, p_0, p_1, p_2, tmp1, tmp5, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  nai_000_100_0 = d_1*p_0 - d_2*u_0; // local, weighs 4
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  nai_000_100_1 = d_2*p_0 - nai_000_000_2*u_0; // local, weighs 4
  tmp1 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  nai_000_200_0 = tmp1 + nai_000_100_0*p_0 - nai_000_100_1*u_0; // local, weighs 5
  d_0 = 6.28318530717959*d_0*gaux(usq, 3); // local+recycle, weighs 3
  usq = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto+recycle, weighs 5
  nai_000_200_1 = usq + nai_000_100_1*p_0 - u_0*(nai_000_000_2*p_0 - d_0*u_0); // local, weighs 9
  tmp5 = pow(a_a,0.75)*pow(b_a,2.25); // auto, weighs 3
  out[0] = 1.04921512347081*tmp5*(nai_000_200_0*p_0 + (nai_000_100_0 - nai_000_100_1)/ab_a - nai_000_200_1*u_0); // final, weighs 11
  out[1] = 2.34611633910158*tmp5*(nai_000_200_0*p_1 - nai_000_200_1*u_1); // final, weighs 6
  out[2] = 2.34611633910158*tmp5*(nai_000_200_0*p_2 - nai_000_200_1*u_2); // final, weighs 6
  nai_000_100_0 = d_1*p_1 - d_2*u_1; // local+recycle, weighs 4
  nai_000_100_1 = d_2*p_1 - nai_000_000_2*u_1; // local+recycle, weighs 4
  nai_000_200_0 = tmp1 + nai_000_100_0*p_1 - nai_000_100_1*u_1; // local+recycle, weighs 5
  nai_000_200_1 = usq + nai_000_100_1*p_1 - u_1*(nai_000_000_2*p_1 - d_0*u_1); // local+recycle, weighs 9
  out[3] = 2.34611633910158*tmp5*(nai_000_200_0*p_0 - nai_000_200_1*u_0); // final, weighs 6
  d_1 = d_1*p_2 - d_2*u_2; // local+recycle, weighs 4
  d_2 = d_2*p_2 - nai_000_000_2*u_2; // local+recycle, weighs 4
  d_0 = nai_000_000_2*p_2 - d_0*u_2; // local+recycle, weighs 4
  out[4] = 4.06359269979142*tmp5*(p_0*(d_1*p_1 - d_2*u_1) - u_0*(d_2*p_1 - d_0*u_1)); // final, weighs 14
  nai_000_000_2 = tmp1 + d_1*p_2 - d_2*u_2; // local+recycle, weighs 5
  d_0 = usq + d_2*p_2 - d_0*u_2; // local+recycle, weighs 5
  out[5] = 2.34611633910158*tmp5*(nai_000_000_2*p_0 - d_0*u_0); // final, weighs 6
  out[6] = 1.04921512347081*tmp5*(nai_000_200_0*p_1 + (nai_000_100_0 - nai_000_100_1)/ab_a - nai_000_200_1*u_1); // final, weighs 11
  out[7] = 2.34611633910158*tmp5*(nai_000_200_0*p_2 - nai_000_200_1*u_2); // final, weighs 6
  out[8] = 2.34611633910158*tmp5*(nai_000_000_2*p_1 - d_0*u_1); // final, weighs 6
  out[9] = 1.04921512347081*tmp5*(nai_000_000_2*p_2 + (d_1 - d_2)/ab_a - d_0*u_2); // final, weighs 11
  // total weight = 227
}

static void gint2_nai_P_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 52
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_000_001_0, nai_000_001_1, nai_000_001_2, nai_000_010_1, nai_000_010_2, nai_000_020_0, nai_000_020_1, nai_000_020_2, nai_000_030_0, nai_000_100_0, nai_000_100_1, nai_000_100_2, nai_000_200_0, nai_000_200_1, nai_000_200_2, nai_000_300_0, nai_010_200_1, nai_100_001_1, nai_100_020_0, nai_100_020_1, p_0, p_1, p_2, tmp0, tmp1, tmp2, tmp21, tmp23, tmp25, tmp26, tmp27, tmp28, tmp29, tmp3, tmp30, tmp31, tmp33, tmp4, tmp5, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  tmp33 = d_2*u_0; // auto, weighs 1
  nai_000_100_0 = -tmp33 + d_1*p_0; // local, weighs 3
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  tmp30 = nai_000_000_2*u_0; // auto, weighs 1
  nai_000_100_1 = -tmp30 + d_2*p_0; // local, weighs 3
  tmp2 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  nai_000_200_0 = tmp2 + nai_000_100_0*p_0 - nai_000_100_1*u_0; // local, weighs 5
  nai_000_000_3 = 6.28318530717959*d_0*gaux(usq, 3); // local, weighs 3
  tmp27 = nai_000_000_3*u_0; // auto, weighs 1
  nai_000_100_2 = -tmp27 + nai_000_000_2*p_0; // local, weighs 3
  tmp1 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto, weighs 5
  nai_000_200_1 = tmp1 + nai_000_100_1*p_0 - nai_000_100_2*u_0; // local, weighs 5
  nai_000_100_0 = (nai_000_100_0 - nai_000_100_1)/ab_a - nai_000_200_1*u_0; // auto+recycle, weighs 7
  nai_000_300_0 = nai_000_100_0 + nai_000_200_0*p_0; // local, weighs 2
  d_0 = 6.28318530717959*d_0*gaux(usq, 4); // local+recycle, weighs 3
  usq = d_0*u_0; // auto+recycle, weighs 1
  tmp0 = (0.5*nai_000_000_2 - 0.5*nai_000_000_3)/ab_a; // auto, weighs 5
  nai_000_200_2 = tmp0 + nai_000_100_2*p_0 - u_0*(-usq + nai_000_000_3*p_0); // local, weighs 8
  nai_000_100_1 = (nai_000_100_1 - nai_000_100_2)/ab_a - nai_000_200_2*u_0; // auto+recycle, weighs 7
  nai_000_100_2 = nai_000_100_1 + nai_000_200_1*p_0; // local+recycle, weighs 2
  tmp21 = pow(a_a,1.25)*pow(b_a,2.25); // auto, weighs 3
  out[0] = 2.09843024694163*tmp21*(nai_000_300_0*v_0 + (1.5*nai_000_200_0 - 1.5*nai_000_200_1)/ab_a - nai_000_100_2*u_0); // final, weighs 12
  nai_000_100_0 = nai_000_100_0 + nai_000_200_0*v_0; // local+recycle, weighs 2
  nai_000_100_1 = nai_000_100_1 + nai_000_200_1*v_0; // local+recycle, weighs 2
  out[1] = 4.69223267820315*tmp21*(nai_000_100_0*p_1 - nai_000_100_1*u_1); // final, weighs 6
  out[2] = 4.69223267820315*tmp21*(nai_000_100_0*p_2 - nai_000_100_1*u_2); // final, weighs 6
  nai_000_100_0 = d_2*u_1; // auto+recycle, weighs 1
  nai_000_100_1 = -nai_000_100_0 + d_1*p_1; // local+recycle, weighs 3
  tmp29 = nai_000_000_2*u_1; // auto, weighs 1
  nai_000_010_1 = -tmp29 + d_2*p_1; // local, weighs 3
  nai_000_020_0 = tmp2 + nai_000_100_1*p_1 - nai_000_010_1*u_1; // local, weighs 5
  tmp26 = nai_000_000_3*u_1; // auto, weighs 1
  nai_000_010_2 = -tmp26 + nai_000_000_2*p_1; // local, weighs 3
  nai_000_020_1 = tmp1 + nai_000_010_1*p_1 - nai_000_010_2*u_1; // local, weighs 5
  nai_100_020_0 = nai_000_020_0*v_0 - nai_000_020_1*u_0; // local, weighs 4
  tmp23 = d_0*u_1; // auto, weighs 1
  nai_000_020_2 = tmp0 + nai_000_010_2*p_1 - u_1*(-tmp23 + nai_000_000_3*p_1); // local, weighs 8
  nai_100_020_1 = nai_000_020_1*v_0 - nai_000_020_2*u_0; // local, weighs 4
  tmp5 = (0.5*nai_000_020_0 - 0.5*nai_000_020_1)/ab_a; // auto, weighs 5
  out[3] = 4.69223267820315*tmp21*(tmp5 + nai_100_020_0*p_0 - nai_100_020_1*u_0); // final, weighs 7
  tmp30 = -tmp30 + d_2*v_0; // local+recycle, weighs 3
  tmp27 = -tmp27 + nai_000_000_2*v_0; // local+recycle, weighs 3
  nai_100_001_1 = p_2*tmp30 - tmp27*u_2; // local, weighs 4
  tmp31 = d_2*u_2; // auto, weighs 1
  nai_000_001_0 = -tmp31 + d_1*p_2; // local, weighs 3
  tmp28 = nai_000_000_2*u_2; // auto, weighs 1
  nai_000_001_1 = -tmp28 + d_2*p_2; // local, weighs 3
  tmp25 = nai_000_000_3*u_2; // auto, weighs 1
  nai_000_001_2 = -tmp25 + nai_000_000_2*p_2; // local, weighs 3
  out[4] = 8.12718539958285*tmp21*(p_0*(p_1*(p_2*(-tmp33 + d_1*v_0) - tmp30*u_2) - nai_100_001_1*u_1) + (0.5*nai_000_001_0*p_1 + 0.5*nai_000_001_2*u_1 - 0.5*nai_000_001_1*p_1 - 0.5*nai_000_001_1*u_1)/ab_a - u_0*(nai_100_001_1*p_1 - u_1*(p_2*tmp27 - u_2*(-usq + nai_000_000_3*v_0)))); // final, weighs 42
  nai_100_001_1 = tmp2 + nai_000_001_0*p_2 - nai_000_001_1*u_2; // local+recycle, weighs 5
  tmp27 = tmp1 + nai_000_001_1*p_2 - nai_000_001_2*u_2; // local+recycle, weighs 5
  tmp30 = nai_100_001_1*v_0 - tmp27*u_0; // local+recycle, weighs 4
  d_0 = d_0*u_2; // auto+recycle, weighs 1
  tmp33 = tmp0 + nai_000_001_2*p_2 - u_2*(-d_0 + nai_000_000_3*p_2); // local+recycle, weighs 8
  usq = tmp27*v_0 - tmp33*u_0; // local+recycle, weighs 4
  tmp3 = (0.5*nai_100_001_1 - 0.5*tmp27)/ab_a; // auto, weighs 5
  out[5] = 4.69223267820315*tmp21*(tmp3 + p_0*tmp30 - u_0*usq); // final, weighs 7
  nai_000_100_1 = (nai_000_100_1 - nai_000_010_1)/ab_a - nai_000_020_1*u_1; // auto+recycle, weighs 7
  nai_000_030_0 = nai_000_100_1 + nai_000_020_0*p_1; // local, weighs 2
  nai_000_010_1 = (nai_000_010_1 - nai_000_010_2)/ab_a - nai_000_020_2*u_1; // auto+recycle, weighs 7
  nai_000_010_2 = nai_000_010_1 + nai_000_020_1*p_1; // local+recycle, weighs 2
  out[6] = 2.09843024694163*tmp21*(nai_000_030_0*v_0 - nai_000_010_2*u_0); // final, weighs 6
  out[7] = 4.69223267820315*tmp21*(nai_100_020_0*p_2 - nai_100_020_1*u_2); // final, weighs 6
  out[8] = 4.69223267820315*tmp21*(p_1*tmp30 - u_1*usq); // final, weighs 6
  nai_100_020_0 = (nai_000_001_0 - nai_000_001_1)/ab_a - tmp27*u_2; // auto+recycle, weighs 7
  nai_100_020_1 = nai_100_020_0 + nai_100_001_1*p_2; // local+recycle, weighs 2
  tmp30 = (nai_000_001_1 - nai_000_001_2)/ab_a - tmp33*u_2; // auto+recycle, weighs 7
  usq = tmp30 + p_2*tmp27; // local+recycle, weighs 2
  out[9] = 2.09843024694163*tmp21*(nai_100_020_1*v_0 - u_0*usq); // final, weighs 6
  out[10] = 2.09843024694163*tmp21*(nai_000_300_0*v_1 - nai_000_100_2*u_1); // final, weighs 6
  v_0 = nai_000_200_0*v_1 - nai_000_200_1*u_1; // local+recycle, weighs 4
  nai_010_200_1 = nai_000_200_1*v_1 - nai_000_200_2*u_1; // local, weighs 4
  tmp4 = (0.5*nai_000_200_0 - 0.5*nai_000_200_1)/ab_a; // auto, weighs 5
  out[11] = 4.69223267820315*tmp21*(tmp4 + p_1*v_0 - nai_010_200_1*u_1); // final, weighs 7
  out[12] = 4.69223267820315*tmp21*(p_2*v_0 - nai_010_200_1*u_2); // final, weighs 6
  nai_000_100_1 = nai_000_100_1 + nai_000_020_0*v_1; // local+recycle, weighs 2
  nai_000_010_1 = nai_000_010_1 + nai_000_020_1*v_1; // local+recycle, weighs 2
  out[13] = 4.69223267820315*tmp21*(nai_000_100_1*p_0 - nai_000_010_1*u_0); // final, weighs 6
  nai_010_200_1 = -tmp29 + d_2*v_1; // local+recycle, weighs 3
  tmp26 = -tmp26 + nai_000_000_2*v_1; // local+recycle, weighs 3
  tmp29 = nai_010_200_1*p_2 - tmp26*u_2; // local+recycle, weighs 4
  out[14] = 8.12718539958285*tmp21*(p_0*(p_1*(p_2*(-nai_000_100_0 + d_1*v_1) - nai_010_200_1*u_2) + (0.5*nai_000_001_0 - 0.5*nai_000_001_1)/ab_a - tmp29*u_1) - u_0*(p_1*tmp29 + (0.5*nai_000_001_1 - 0.5*nai_000_001_2)/ab_a - u_1*(p_2*tmp26 - u_2*(-tmp23 + nai_000_000_3*v_1)))); // final, weighs 40
  nai_000_001_0 = nai_100_001_1*v_1 - tmp27*u_1; // local+recycle, weighs 4
  nai_000_001_1 = tmp27*v_1 - tmp33*u_1; // local+recycle, weighs 4
  out[15] = 4.69223267820315*tmp21*(nai_000_001_0*p_0 - nai_000_001_1*u_0); // final, weighs 6
  out[16] = 2.09843024694163*tmp21*(nai_000_030_0*v_1 + (1.5*nai_000_020_0 - 1.5*nai_000_020_1)/ab_a - nai_000_010_2*u_1); // final, weighs 12
  out[17] = 4.69223267820315*tmp21*(nai_000_100_1*p_2 - nai_000_010_1*u_2); // final, weighs 6
  out[18] = 4.69223267820315*tmp21*(tmp3 + nai_000_001_0*p_1 - nai_000_001_1*u_1); // final, weighs 7
  out[19] = 2.09843024694163*tmp21*(nai_100_020_1*v_1 - u_1*usq); // final, weighs 6
  out[20] = 2.09843024694163*tmp21*(nai_000_300_0*v_2 - nai_000_100_2*u_2); // final, weighs 6
  nai_000_001_0 = nai_000_200_0*v_2 - nai_000_200_1*u_2; // local+recycle, weighs 4
  nai_000_001_1 = nai_000_200_1*v_2 - nai_000_200_2*u_2; // local+recycle, weighs 4
  out[21] = 4.69223267820315*tmp21*(nai_000_001_0*p_1 - nai_000_001_1*u_1); // final, weighs 6
  out[22] = 4.69223267820315*tmp21*(tmp4 + nai_000_001_0*p_2 - nai_000_001_1*u_2); // final, weighs 7
  nai_000_001_0 = nai_000_020_0*v_2 - nai_000_020_1*u_2; // local+recycle, weighs 4
  nai_000_001_1 = nai_000_020_1*v_2 - nai_000_020_2*u_2; // local+recycle, weighs 4
  out[23] = 4.69223267820315*tmp21*(nai_000_001_0*p_0 - nai_000_001_1*u_0); // final, weighs 6
  d_2 = -tmp28 + d_2*v_2; // local+recycle, weighs 3
  nai_000_000_2 = -tmp25 + nai_000_000_2*v_2; // local+recycle, weighs 3
  nai_000_001_2 = tmp1 + d_2*p_2 - nai_000_000_2*u_2; // local+recycle, weighs 5
  out[24] = 8.12718539958285*tmp21*(p_0*(p_1*(tmp2 + p_2*(-tmp31 + d_1*v_2) - d_2*u_2) - nai_000_001_2*u_1) - u_0*(nai_000_001_2*p_1 - u_1*(tmp0 + nai_000_000_2*p_2 - u_2*(-d_0 + nai_000_000_3*v_2)))); // final, weighs 30
  d_0 = nai_100_020_0 + nai_100_001_1*v_2; // local+recycle, weighs 2
  d_1 = tmp30 + tmp27*v_2; // local+recycle, weighs 2
  out[25] = 4.69223267820315*tmp21*(d_0*p_0 - d_1*u_0); // final, weighs 6
  out[26] = 2.09843024694163*tmp21*(nai_000_030_0*v_2 - nai_000_010_2*u_2); // final, weighs 6
  out[27] = 4.69223267820315*tmp21*(tmp5 + nai_000_001_0*p_2 - nai_000_001_1*u_2); // final, weighs 7
  out[28] = 4.69223267820315*tmp21*(d_0*p_1 - d_1*u_1); // final, weighs 6
  out[29] = 2.09843024694163*tmp21*(nai_100_020_1*v_2 + (1.5*nai_100_001_1 - 1.5*tmp27)/ab_a - u_2*usq); // final, weighs 12
  // total weight = 643
}

static void gint2_nai_cD_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 90
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_000_000_4, nai_000_001_1, nai_000_001_2, nai_000_001_3, nai_000_002_0, nai_000_002_1, nai_000_002_2, nai_000_002_3, nai_000_003_2, nai_000_010_1, nai_000_010_2, nai_000_010_3, nai_000_020_0, nai_000_020_1, nai_000_020_2, nai_000_020_3, nai_000_030_0, nai_000_030_1, nai_000_100_0, nai_000_100_1, nai_000_100_2, nai_000_100_3, nai_000_200_0, nai_000_200_1, nai_000_200_2, nai_000_200_3, nai_000_300_0, nai_000_300_1, nai_010_001_0, nai_010_001_1, nai_010_001_2, nai_010_100_1, nai_010_200_0, nai_010_200_1, nai_010_300_0, nai_010_300_1, nai_100_001_1, nai_100_001_2, nai_100_002_0, nai_100_002_1, nai_100_020_0, nai_100_020_1, nai_110_200_0, nai_200_000_0, nai_200_000_1, nai_200_000_2, nai_200_100_1, p_0, p_1, p_2, tmp0, tmp1, tmp15, tmp2, tmp26, tmp27, tmp3, tmp30, tmp31, tmp32, tmp38, tmp43, tmp5, tmp53, tmp54, tmp55, tmp56, tmp64, tmp66, tmp67, tmp68, tmp69, tmp70, tmp71, tmp73, tmp74, tmp77, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  tmp77 = d_2*u_0; // auto, weighs 1
  nai_000_100_0 = -tmp77 + d_1*p_0; // local, weighs 3
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  tmp74 = nai_000_000_2*u_0; // auto, weighs 1
  nai_000_100_1 = -tmp74 + d_2*p_0; // local, weighs 3
  tmp2 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  nai_000_200_0 = tmp2 + nai_000_100_0*p_0 - nai_000_100_1*u_0; // local, weighs 5
  nai_000_000_3 = 6.28318530717959*d_0*gaux(usq, 3); // local, weighs 3
  tmp71 = nai_000_000_3*u_0; // auto, weighs 1
  nai_000_100_2 = -tmp71 + nai_000_000_2*p_0; // local, weighs 3
  tmp1 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto, weighs 5
  nai_000_200_1 = tmp1 + nai_000_100_1*p_0 - nai_000_100_2*u_0; // local, weighs 5
  nai_000_100_0 = (nai_000_100_0 - nai_000_100_1)/ab_a - nai_000_200_1*u_0; // auto+recycle, weighs 7
  nai_000_300_0 = nai_000_100_0 + nai_000_200_0*p_0; // local, weighs 2
  nai_000_000_4 = 6.28318530717959*d_0*gaux(usq, 4); // local, weighs 3
  tmp68 = nai_000_000_4*u_0; // auto, weighs 1
  nai_000_100_3 = -tmp68 + nai_000_000_3*p_0; // local, weighs 3
  tmp0 = (0.5*nai_000_000_2 - 0.5*nai_000_000_3)/ab_a; // auto, weighs 5
  nai_000_200_2 = tmp0 + nai_000_100_2*p_0 - nai_000_100_3*u_0; // local, weighs 5
  nai_000_100_1 = (nai_000_100_1 - nai_000_100_2)/ab_a - nai_000_200_2*u_0; // auto+recycle, weighs 7
  nai_000_300_1 = nai_000_100_1 + nai_000_200_1*p_0; // local, weighs 2
  d_0 = 6.28318530717959*d_0*gaux(usq, 5); // local+recycle, weighs 3
  usq = d_0*u_0; // auto+recycle, weighs 1
  tmp3 = (0.5*nai_000_000_3 - 0.5*nai_000_000_4)/ab_a; // auto, weighs 5
  nai_000_200_3 = tmp3 + nai_000_100_3*p_0 - u_0*(-usq + nai_000_000_4*p_0); // local, weighs 8
  nai_000_100_2 = nai_000_200_2*p_0 + (nai_000_100_2 - nai_000_100_3)/ab_a - nai_000_200_3*u_0; // local+recycle, weighs 9
  nai_000_100_0 = nai_000_100_0 + nai_000_200_0*v_0; // local+recycle, weighs 2
  nai_000_100_1 = nai_000_100_1 + nai_000_200_1*v_0; // local+recycle, weighs 2
  nai_000_100_3 = 0.5*nai_000_300_0 - 0.5*nai_000_300_1; // auto+recycle, weighs 3
  tmp53 = pow(a_a,1.75)*pow(b_a,2.25); // auto, weighs 3
  out[0] = 2.4230585358948*tmp53*(v_0*(nai_000_300_0*v_0 + (1.5*nai_000_200_0 - 1.5*nai_000_200_1)/ab_a - nai_000_300_1*u_0) + (nai_000_100_3 + 1.5*nai_000_100_0 - 1.5*nai_000_100_1)/ab_a - u_0*(nai_000_300_1*v_0 + (1.5*nai_000_200_1 - 1.5*nai_000_200_2)/ab_a - nai_000_100_2*u_0)); // final, weighs 33
  tmp77 = -tmp77 + d_1*v_0; // local+recycle, weighs 3
  tmp74 = -tmp74 + d_2*v_0; // local+recycle, weighs 3
  tmp43 = tmp2 - tmp74*u_0; // auto, weighs 3
  nai_200_000_0 = tmp43 + tmp77*v_0; // local, weighs 2
  tmp71 = -tmp71 + nai_000_000_2*v_0; // local+recycle, weighs 3
  tmp27 = tmp1 - tmp71*u_0; // auto, weighs 3
  nai_200_000_1 = tmp27 + tmp74*v_0; // local, weighs 2
  tmp68 = -tmp68 + nai_000_000_3*v_0; // local+recycle, weighs 3
  tmp26 = tmp0 - tmp68*u_0; // auto, weighs 3
  nai_200_000_2 = tmp26 + tmp71*v_0; // local, weighs 2
  nai_200_100_1 = nai_200_000_1*p_0 + (tmp74 - tmp71)/ab_a - nai_200_000_2*u_0; // local, weighs 9
  tmp27 = tmp27 + p_0*tmp74; // local+recycle, weighs 2
  tmp32 = 0.5*nai_200_000_0 - 0.5*nai_200_000_1; // auto, weighs 3
  tmp43 = p_0*(nai_200_000_0*p_0 + (tmp77 - tmp74)/ab_a - nai_200_000_1*u_0) + (tmp32 + tmp43 - tmp27 + p_0*tmp77)/ab_a - nai_200_100_1*u_0; // local+recycle, weighs 21
  usq = tmp3 + tmp68*v_0 - u_0*(-usq + nai_000_000_4*v_0); // local+recycle, weighs 8
  tmp31 = 0.5*nai_200_000_1 - 0.5*nai_200_000_2; // auto, weighs 3
  nai_200_100_1 = nai_200_100_1*p_0 + (tmp27 + tmp31 - tmp26 - p_0*tmp71)/ab_a - u_0*(nai_200_000_2*p_0 + (tmp71 - tmp68)/ab_a - u_0*usq); // local+recycle, weighs 22
  out[1] = 5.4181235997219*tmp53*(p_1*tmp43 - nai_200_100_1*u_1); // final, weighs 6
  out[2] = 5.4181235997219*tmp53*(p_2*tmp43 - nai_200_100_1*u_2); // final, weighs 6
  nai_200_100_1 = nai_200_000_1*p_1 - nai_200_000_2*u_1; // local+recycle, weighs 4
  tmp26 = tmp32/ab_a; // auto+recycle, weighs 2
  tmp27 = tmp26 + p_1*(nai_200_000_0*p_1 - nai_200_000_1*u_1) - nai_200_100_1*u_1; // local+recycle, weighs 9
  tmp31 = tmp31/ab_a; // auto+recycle, weighs 2
  nai_200_100_1 = tmp31 + nai_200_100_1*p_1 - u_1*(nai_200_000_2*p_1 - u_1*usq); // local+recycle, weighs 9
  tmp32 = d_2*u_1; // auto+recycle, weighs 1
  tmp43 = -tmp32 + d_1*p_1; // local+recycle, weighs 3
  tmp73 = nai_000_000_2*u_1; // auto, weighs 1
  nai_000_010_1 = -tmp73 + d_2*p_1; // local, weighs 3
  nai_000_020_0 = tmp2 + p_1*tmp43 - nai_000_010_1*u_1; // local, weighs 5
  tmp70 = nai_000_000_3*u_1; // auto, weighs 1
  nai_000_010_2 = -tmp70 + nai_000_000_2*p_1; // local, weighs 3
  nai_000_020_1 = tmp1 + nai_000_010_1*p_1 - nai_000_010_2*u_1; // local, weighs 5
  nai_100_020_0 = nai_000_020_0*v_0 - nai_000_020_1*u_0; // local, weighs 4
  tmp67 = nai_000_000_4*u_1; // auto, weighs 1
  nai_000_010_3 = -tmp67 + nai_000_000_3*p_1; // local, weighs 3
  nai_000_020_2 = tmp0 + nai_000_010_2*p_1 - nai_000_010_3*u_1; // local, weighs 5
  nai_100_020_1 = nai_000_020_1*v_0 - nai_000_020_2*u_0; // local, weighs 4
  out[3] = 5.4181235997219*tmp53*(p_0*tmp27 + (nai_100_020_0 - nai_100_020_1)/ab_a - nai_200_100_1*u_0); // final, weighs 11
  nai_200_000_0 = nai_200_000_0*p_2 - nai_200_000_1*u_2; // local+recycle, weighs 4
  nai_200_000_1 = nai_200_000_1*p_2 - nai_200_000_2*u_2; // local+recycle, weighs 4
  nai_200_000_2 = nai_200_000_2*p_2 - u_2*usq; // local+recycle, weighs 4
  usq = p_2*tmp77 - tmp74*u_2; // local+recycle, weighs 4
  nai_100_001_1 = p_2*tmp74 - tmp71*u_2; // local, weighs 4
  nai_100_001_2 = p_2*tmp71 - tmp68*u_2; // local, weighs 4
  out[4] = 9.3844653564063*tmp53*(p_0*(nai_200_000_0*p_1 - nai_200_000_1*u_1) + (nai_100_001_2*u_1 + p_1*usq - nai_100_001_1*p_1 - nai_100_001_1*u_1)/ab_a - u_0*(nai_200_000_1*p_1 - nai_200_000_2*u_1)); // final, weighs 26
  nai_200_000_0 = tmp26 + nai_200_000_0*p_2 - nai_200_000_1*u_2; // local+recycle, weighs 5
  nai_200_000_1 = tmp31 + nai_200_000_1*p_2 - nai_200_000_2*u_2; // local+recycle, weighs 5
  nai_200_000_2 = d_2*u_2; // auto+recycle, weighs 1
  tmp26 = -nai_200_000_2 + d_1*p_2; // local+recycle, weighs 3
  tmp31 = nai_000_000_2*u_2; // auto+recycle, weighs 1
  nai_000_001_1 = -tmp31 + d_2*p_2; // local, weighs 3
  nai_000_002_0 = tmp2 + p_2*tmp26 - nai_000_001_1*u_2; // local, weighs 5
  tmp69 = nai_000_000_3*u_2; // auto, weighs 1
  nai_000_001_2 = -tmp69 + nai_000_000_2*p_2; // local, weighs 3
  nai_000_002_1 = tmp1 + nai_000_001_1*p_2 - nai_000_001_2*u_2; // local, weighs 5
  nai_100_002_0 = nai_000_002_0*v_0 - nai_000_002_1*u_0; // local, weighs 4
  tmp66 = nai_000_000_4*u_2; // auto, weighs 1
  nai_000_001_3 = -tmp66 + nai_000_000_3*p_2; // local, weighs 3
  nai_000_002_2 = tmp0 + nai_000_001_2*p_2 - nai_000_001_3*u_2; // local, weighs 5
  nai_100_002_1 = nai_000_002_1*v_0 - nai_000_002_2*u_0; // local, weighs 4
  out[5] = 5.4181235997219*tmp53*(nai_200_000_0*p_0 + (nai_100_002_0 - nai_100_002_1)/ab_a - nai_200_000_1*u_0); // final, weighs 11
  tmp43 = (tmp43 - nai_000_010_1)/ab_a - nai_000_020_1*u_1; // auto+recycle, weighs 7
  nai_000_030_0 = tmp43 + nai_000_020_0*p_1; // local, weighs 2
  nai_000_010_1 = (nai_000_010_1 - nai_000_010_2)/ab_a - nai_000_020_2*u_1; // auto+recycle, weighs 7
  nai_000_030_1 = nai_000_010_1 + nai_000_020_1*p_1; // local, weighs 2
  tmp64 = d_0*u_1; // auto, weighs 1
  nai_000_020_3 = tmp3 + nai_000_010_3*p_1 - u_1*(-tmp64 + nai_000_000_4*p_1); // local, weighs 8
  nai_000_010_2 = (nai_000_010_2 - nai_000_010_3)/ab_a - nai_000_020_3*u_1; // auto+recycle, weighs 7
  nai_000_010_3 = nai_000_010_2 + nai_000_020_2*p_1; // local+recycle, weighs 2
  tmp30 = 0.5*nai_000_030_0 - 0.5*nai_000_030_1; // auto, weighs 3
  tmp5 = tmp30/ab_a; // auto, weighs 2
  out[6] = 2.4230585358948*tmp53*(tmp5 + v_0*(nai_000_030_0*v_0 - nai_000_030_1*u_0) - u_0*(nai_000_030_1*v_0 - nai_000_010_3*u_0)); // final, weighs 15
  out[7] = 5.4181235997219*tmp53*(p_2*tmp27 - nai_200_100_1*u_2); // final, weighs 6
  out[8] = 5.4181235997219*tmp53*(nai_200_000_0*p_1 - nai_200_000_1*u_1); // final, weighs 6
  nai_200_000_0 = (tmp26 - nai_000_001_1)/ab_a - nai_000_002_1*u_2; // auto+recycle, weighs 7
  nai_200_000_1 = nai_200_000_0 + nai_000_002_0*p_2; // local+recycle, weighs 2
  nai_200_100_1 = (nai_000_001_1 - nai_000_001_2)/ab_a - nai_000_002_2*u_2; // auto+recycle, weighs 7
  tmp27 = nai_200_100_1 + nai_000_002_1*p_2; // local+recycle, weighs 2
  d_0 = d_0*u_2; // auto+recycle, weighs 1
  nai_000_002_3 = tmp3 + nai_000_001_3*p_2 - u_2*(-d_0 + nai_000_000_4*p_2); // local, weighs 8
  nai_000_001_3 = (nai_000_001_2 - nai_000_001_3)/ab_a - nai_000_002_3*u_2; // auto+recycle, weighs 7
  nai_000_003_2 = nai_000_001_3 + nai_000_002_2*p_2; // local, weighs 2
  tmp38 = 0.5*nai_200_000_1 - 0.5*tmp27; // auto, weighs 3
  tmp15 = tmp38/ab_a; // auto, weighs 2
  out[9] = 2.4230585358948*tmp53*(tmp15 + v_0*(nai_200_000_1*v_0 - tmp27*u_0) - u_0*(tmp27*v_0 - nai_000_003_2*u_0)); // final, weighs 15
  nai_010_300_0 = nai_000_300_0*v_1 - nai_000_300_1*u_1; // local, weighs 4
  nai_010_300_1 = nai_000_300_1*v_1 - nai_000_100_2*u_1; // local, weighs 4
  nai_010_200_0 = nai_000_200_0*v_1 - nai_000_200_1*u_1; // local, weighs 4
  nai_010_200_1 = nai_000_200_1*v_1 - nai_000_200_2*u_1; // local, weighs 4
  out[10] = 4.19686049388326*tmp53*(nai_010_300_0*v_0 + (1.5*nai_010_200_0 - 1.5*nai_010_200_1)/ab_a - nai_010_300_1*u_0); // final, weighs 12
  tmp32 = -tmp32 + d_1*v_1; // local+recycle, weighs 3
  tmp73 = -tmp73 + d_2*v_1; // local+recycle, weighs 3
  tmp56 = tmp73*u_0; // auto, weighs 1
  tmp70 = -tmp70 + nai_000_000_2*v_1; // local+recycle, weighs 3
  tmp55 = tmp70*u_0; // auto, weighs 1
  nai_010_100_1 = -tmp55 + p_0*tmp73; // local, weighs 3
  nai_110_200_0 = nai_010_200_0*v_0 + (-nai_010_100_1 - tmp56 + p_0*tmp32)/ab_a - nai_010_200_1*u_0; // local, weighs 12
  tmp67 = -tmp67 + nai_000_000_3*v_1; // local+recycle, weighs 3
  tmp54 = tmp67*u_0; // auto, weighs 1
  nai_010_100_1 = nai_010_200_1*v_0 + (nai_010_100_1 + tmp54 - p_0*tmp70)/ab_a - u_0*(nai_000_200_2*v_1 - nai_000_200_3*u_1); // local+recycle, weighs 15
  nai_000_100_0 = (0.5*nai_000_100_0 - 0.5*nai_000_100_1)/ab_a; // auto+recycle, weighs 5
  out[11] = 9.3844653564063*tmp53*(nai_000_100_0 + nai_110_200_0*p_1 - nai_010_100_1*u_1); // final, weighs 7
  out[12] = 9.3844653564063*tmp53*(nai_110_200_0*p_2 - nai_010_100_1*u_2); // final, weighs 6
  nai_000_100_1 = tmp43 + nai_000_020_0*v_1; // local+recycle, weighs 2
  nai_000_010_1 = nai_000_010_1 + nai_000_020_1*v_1; // local+recycle, weighs 2
  nai_010_100_1 = nai_000_100_1*v_0 - nai_000_010_1*u_0; // local+recycle, weighs 4
  nai_000_010_2 = nai_000_010_1*v_0 - u_0*(nai_000_010_2 + nai_000_020_2*v_1); // local+recycle, weighs 6
  nai_110_200_0 = (0.5*nai_000_100_1 - 0.5*nai_000_010_1)/ab_a; // auto+recycle, weighs 5
  out[13] = 9.3844653564063*tmp53*(nai_110_200_0 + nai_010_100_1*p_0 - nai_000_010_2*u_0); // final, weighs 7
  tmp43 = -tmp55 + tmp73*v_0; // local+recycle, weighs 3
  tmp54 = -tmp54 + tmp70*v_0; // local+recycle, weighs 3
  tmp55 = p_2*tmp43 - tmp54*u_2; // local+recycle, weighs 4
  tmp64 = -tmp64 + nai_000_000_4*v_1; // local+recycle, weighs 3
  nai_010_001_0 = p_2*tmp32 - tmp73*u_2; // local, weighs 4
  nai_010_001_1 = p_2*tmp73 - tmp70*u_2; // local, weighs 4
  nai_010_001_2 = p_2*tmp70 - tmp67*u_2; // local, weighs 4
  out[14] = 16.2543707991657*tmp53*(p_0*(p_1*(p_2*(-tmp56 + tmp32*v_0) - tmp43*u_2) + (0.5*usq - 0.5*nai_100_001_1)/ab_a - tmp55*u_1) + (0.5*nai_010_001_0*p_1 + 0.5*nai_010_001_2*u_1 + 0.5*(0.5*tmp26 - 0.5*nai_000_001_1)/ab_a - 0.5*nai_010_001_1*p_1 - 0.5*nai_010_001_1*u_1 - 0.5*(0.5*nai_000_001_1 - 0.5*nai_000_001_2)/ab_a)/ab_a - u_0*(p_1*tmp55 + (0.5*nai_100_001_1 - 0.5*nai_100_001_2)/ab_a - u_1*(p_2*tmp54 - u_2*(tmp67*v_0 - tmp64*u_0)))); // final, weighs 69
  nai_000_001_1 = nai_000_002_0*v_1 - nai_000_002_1*u_1; // local+recycle, weighs 4
  nai_000_001_2 = nai_000_002_1*v_1 - nai_000_002_2*u_1; // local+recycle, weighs 4
  nai_100_001_1 = nai_000_001_1*v_0 - nai_000_001_2*u_0; // local+recycle, weighs 4
  nai_000_002_3 = nai_000_001_2*v_0 - u_0*(nai_000_002_2*v_1 - nai_000_002_3*u_1); // local+recycle, weighs 8
  out[15] = 9.3844653564063*tmp53*(nai_100_001_1*p_0 + (0.5*nai_000_001_1 - 0.5*nai_000_001_2)/ab_a - nai_000_002_3*u_0); // final, weighs 12
  nai_100_001_2 = nai_000_030_0*v_1 + (1.5*nai_000_020_0 - 1.5*nai_000_020_1)/ab_a - nai_000_030_1*u_1; // local+recycle, weighs 10
  tmp26 = nai_000_030_1*v_1 + (1.5*nai_000_020_1 - 1.5*nai_000_020_2)/ab_a - nai_000_010_3*u_1; // local+recycle, weighs 10
  out[16] = 4.19686049388326*tmp53*(nai_100_001_2*v_0 - tmp26*u_0); // final, weighs 6
  out[17] = 9.3844653564063*tmp53*(nai_010_100_1*p_2 - nai_000_010_2*u_2); // final, weighs 6
  out[18] = 9.3844653564063*tmp53*(nai_100_001_1*p_1 + (0.5*nai_100_002_0 - 0.5*nai_100_002_1)/ab_a - nai_000_002_3*u_1); // final, weighs 12
  nai_000_002_3 = nai_200_000_1*v_1 - tmp27*u_1; // local+recycle, weighs 4
  nai_000_010_2 = tmp27*v_1 - nai_000_003_2*u_1; // local+recycle, weighs 4
  out[19] = 4.19686049388326*tmp53*(nai_000_002_3*v_0 - nai_000_010_2*u_0); // final, weighs 6
  nai_000_300_0 = nai_000_300_0*v_2 - nai_000_300_1*u_2; // local+recycle, weighs 4
  nai_000_100_2 = nai_000_300_1*v_2 - nai_000_100_2*u_2; // local+recycle, weighs 4
  nai_000_200_0 = nai_000_200_0*v_2 - nai_000_200_1*u_2; // local+recycle, weighs 4
  nai_000_200_1 = nai_000_200_1*v_2 - nai_000_200_2*u_2; // local+recycle, weighs 4
  out[20] = 4.19686049388326*tmp53*(nai_000_300_0*v_0 + (1.5*nai_000_200_0 - 1.5*nai_000_200_1)/ab_a - nai_000_100_2*u_0); // final, weighs 12
  d_1 = -nai_200_000_2 + d_1*v_2; // local+recycle, weighs 3
  d_2 = -tmp31 + d_2*v_2; // local+recycle, weighs 3
  nai_000_300_1 = d_2*u_0; // auto+recycle, weighs 1
  nai_000_000_2 = -tmp69 + nai_000_000_2*v_2; // local+recycle, weighs 3
  nai_010_100_1 = nai_000_000_2*u_0; // auto+recycle, weighs 1
  nai_100_001_1 = -nai_010_100_1 + d_2*p_0; // local+recycle, weighs 3
  nai_100_002_0 = nai_000_200_0*v_0 + (-nai_000_300_1 - nai_100_001_1 + d_1*p_0)/ab_a - nai_000_200_1*u_0; // local+recycle, weighs 12
  nai_000_200_2 = nai_000_200_2*v_2 - nai_000_200_3*u_2; // local+recycle, weighs 4
  nai_000_000_3 = -tmp66 + nai_000_000_3*v_2; // local+recycle, weighs 3
  nai_000_200_3 = nai_000_000_3*u_0; // auto+recycle, weighs 1
  nai_100_001_1 = nai_000_200_1*v_0 + (nai_000_200_3 + nai_100_001_1 - nai_000_000_2*p_0)/ab_a - nai_000_200_2*u_0; // local+recycle, weighs 11
  out[21] = 9.3844653564063*tmp53*(nai_100_002_0*p_1 - nai_100_001_1*u_1); // final, weighs 6
  out[22] = 9.3844653564063*tmp53*(nai_000_100_0 + nai_100_002_0*p_2 - nai_100_001_1*u_2); // final, weighs 7
  nai_000_020_0 = nai_000_020_0*v_2 - nai_000_020_1*u_2; // local+recycle, weighs 4
  nai_000_020_1 = nai_000_020_1*v_2 - nai_000_020_2*u_2; // local+recycle, weighs 4
  nai_000_100_0 = nai_000_020_0*v_0 - nai_000_020_1*u_0; // local+recycle, weighs 4
  nai_000_020_2 = nai_000_020_2*v_2 - nai_000_020_3*u_2; // local+recycle, weighs 4
  nai_000_020_3 = nai_000_020_1*v_0 - nai_000_020_2*u_0; // local+recycle, weighs 4
  out[23] = 9.3844653564063*tmp53*(nai_000_100_0*p_0 + (0.5*nai_000_020_0 - 0.5*nai_000_020_1)/ab_a - nai_000_020_3*u_0); // final, weighs 12
  nai_010_100_1 = -nai_010_100_1 + d_2*v_0; // local+recycle, weighs 3
  nai_000_200_3 = -nai_000_200_3 + nai_000_000_2*v_0; // local+recycle, weighs 3
  nai_100_001_1 = nai_010_100_1*p_2 + (0.5*tmp74 - 0.5*tmp71)/ab_a - nai_000_200_3*u_2; // local+recycle, weighs 10
  d_0 = -d_0 + nai_000_000_4*v_2; // local+recycle, weighs 3
  nai_000_000_4 = tmp2 - d_2*u_2; // auto+recycle, weighs 3
  nai_100_002_0 = nai_000_000_4 + d_1*p_2; // local+recycle, weighs 2
  nai_100_002_1 = tmp1 - nai_000_000_2*u_2; // auto+recycle, weighs 3
  nai_200_000_2 = nai_100_002_1 + d_2*p_2; // local+recycle, weighs 2
  tmp31 = tmp0 - nai_000_000_3*u_2; // auto+recycle, weighs 3
  tmp43 = tmp31 + nai_000_000_2*p_2; // local+recycle, weighs 2
  out[24] = 16.2543707991657*tmp53*(p_0*(p_1*(p_2*(-nai_000_300_1 + d_1*v_0) + (0.5*tmp77 - 0.5*tmp74)/ab_a - nai_010_100_1*u_2) - nai_100_001_1*u_1) + (0.5*nai_100_002_0*p_1 + 0.5*tmp43*u_1 - 0.5*nai_200_000_2*p_1 - 0.5*nai_200_000_2*u_1)/ab_a - u_0*(nai_100_001_1*p_1 - u_1*(nai_000_200_3*p_2 + (0.5*tmp71 - 0.5*tmp68)/ab_a - u_2*(nai_000_000_3*v_0 - d_0*u_0)))); // final, weighs 55
  nai_000_200_3 = nai_200_000_0 + nai_000_002_0*v_2; // local+recycle, weighs 2
  nai_000_300_1 = nai_200_100_1 + nai_000_002_1*v_2; // local+recycle, weighs 2
  nai_010_100_1 = nai_000_200_3*v_0 - nai_000_300_1*u_0; // local+recycle, weighs 4
  nai_000_001_3 = nai_000_001_3 + nai_000_002_2*v_2; // local+recycle, weighs 2
  nai_100_001_1 = nai_000_300_1*v_0 - nai_000_001_3*u_0; // local+recycle, weighs 4
  nai_200_000_0 = (0.5*nai_000_200_3 - 0.5*nai_000_300_1)/ab_a; // auto+recycle, weighs 5
  out[25] = 9.3844653564063*tmp53*(nai_200_000_0 + nai_010_100_1*p_0 - nai_100_001_1*u_0); // final, weighs 7
  nai_000_030_0 = nai_000_030_0*v_2 - nai_000_030_1*u_2; // local+recycle, weighs 4
  nai_000_010_3 = nai_000_030_1*v_2 - nai_000_010_3*u_2; // local+recycle, weighs 4
  out[26] = 4.19686049388326*tmp53*(nai_000_030_0*v_0 - nai_000_010_3*u_0); // final, weighs 6
  out[27] = 9.3844653564063*tmp53*(nai_000_100_0*p_2 + (0.5*nai_100_020_0 - 0.5*nai_100_020_1)/ab_a - nai_000_020_3*u_2); // final, weighs 12
  out[28] = 9.3844653564063*tmp53*(nai_010_100_1*p_1 - nai_100_001_1*u_1); // final, weighs 6
  nai_000_002_0 = nai_200_000_1*v_2 + (1.5*nai_000_002_0 - 1.5*nai_000_002_1)/ab_a - tmp27*u_2; // local+recycle, weighs 10
  nai_000_002_1 = tmp27*v_2 + (1.5*nai_000_002_1 - 1.5*nai_000_002_2)/ab_a - nai_000_003_2*u_2; // local+recycle, weighs 10
  out[29] = 4.19686049388326*tmp53*(nai_000_002_0*v_0 - nai_000_002_1*u_0); // final, weighs 6
  nai_000_002_2 = nai_000_100_3/ab_a; // auto+recycle, weighs 2
  out[30] = 2.4230585358948*tmp53*(nai_000_002_2 + nai_010_300_0*v_1 - nai_010_300_1*u_1); // final, weighs 7
  nai_000_003_2 = tmp2 - tmp73*u_1; // auto+recycle, weighs 3
  nai_000_020_3 = nai_000_003_2 + tmp32*v_1; // local+recycle, weighs 2
  nai_000_030_1 = tmp1 - tmp70*u_1; // auto+recycle, weighs 3
  nai_000_100_0 = nai_000_030_1 + tmp73*v_1; // local+recycle, weighs 2
  nai_000_100_3 = tmp0 - tmp67*u_1; // auto+recycle, weighs 3
  nai_010_100_1 = nai_000_100_3 + tmp70*v_1; // local+recycle, weighs 2
  nai_010_300_0 = nai_000_100_0*p_0 - nai_010_100_1*u_0; // local+recycle, weighs 4
  nai_010_300_1 = 0.5*nai_000_020_3 - 0.5*nai_000_100_0; // auto+recycle, weighs 3
  nai_100_001_1 = nai_010_300_1/ab_a; // auto+recycle, weighs 2
  nai_100_020_0 = nai_100_001_1 + p_0*(nai_000_020_3*p_0 - nai_000_100_0*u_0) - nai_010_300_0*u_0; // local+recycle, weighs 9
  nai_100_020_1 = tmp3 + tmp67*v_1 - tmp64*u_1; // local+recycle, weighs 5
  nai_200_000_1 = 0.5*nai_000_100_0 - 0.5*nai_010_100_1; // auto+recycle, weighs 3
  nai_200_100_1 = nai_200_000_1/ab_a; // auto+recycle, weighs 2
  nai_010_300_0 = nai_200_100_1 + nai_010_300_0*p_0 - u_0*(nai_010_100_1*p_0 - nai_100_020_1*u_0); // local+recycle, weighs 9
  out[31] = 5.4181235997219*tmp53*(nai_100_020_0*p_1 + (nai_010_200_0 - nai_010_200_1)/ab_a - nai_010_300_0*u_1); // final, weighs 11
  out[32] = 5.4181235997219*tmp53*(nai_100_020_0*p_2 - nai_010_300_0*u_2); // final, weighs 6
  nai_010_300_0 = nai_000_100_0*p_1 + (tmp73 - tmp70)/ab_a - nai_010_100_1*u_1; // local+recycle, weighs 9
  nai_000_030_1 = nai_000_030_1 + p_1*tmp73; // local+recycle, weighs 2
  nai_000_003_2 = p_1*(nai_000_020_3*p_1 + (tmp32 - tmp73)/ab_a - nai_000_100_0*u_1) + (nai_000_003_2 + nai_010_300_1 - nai_000_030_1 + p_1*tmp32)/ab_a - nai_010_300_0*u_1; // local+recycle, weighs 21
  nai_000_030_1 = nai_010_300_0*p_1 + (nai_000_030_1 + nai_200_000_1 - nai_000_100_3 - p_1*tmp70)/ab_a - u_1*(nai_010_100_1*p_1 + (tmp70 - tmp67)/ab_a - nai_100_020_1*u_1); // local+recycle, weighs 22
  out[33] = 5.4181235997219*tmp53*(nai_000_003_2*p_0 - nai_000_030_1*u_0); // final, weighs 6
  nai_000_020_3 = nai_000_020_3*p_2 - nai_000_100_0*u_2; // local+recycle, weighs 4
  nai_000_100_0 = nai_000_100_0*p_2 - nai_010_100_1*u_2; // local+recycle, weighs 4
  nai_000_100_3 = nai_010_100_1*p_2 - nai_100_020_1*u_2; // local+recycle, weighs 4
  out[34] = 9.3844653564063*tmp53*(p_0*(nai_000_020_3*p_1 + (nai_010_001_0 - nai_010_001_1)/ab_a - nai_000_100_0*u_1) - u_0*(nai_000_100_0*p_1 + (nai_010_001_1 - nai_010_001_2)/ab_a - nai_000_100_3*u_1)); // final, weighs 24
  nai_000_020_3 = nai_100_001_1 + nai_000_020_3*p_2 - nai_000_100_0*u_2; // local+recycle, weighs 5
  nai_000_100_0 = nai_200_100_1 + nai_000_100_0*p_2 - nai_000_100_3*u_2; // local+recycle, weighs 5
  out[35] = 5.4181235997219*tmp53*(nai_000_020_3*p_0 - nai_000_100_0*u_0); // final, weighs 6
  out[36] = 2.4230585358948*tmp53*(nai_100_001_2*v_1 + (tmp30 + 1.5*nai_000_100_1 - 1.5*nai_000_010_1)/ab_a - tmp26*u_1); // final, weighs 13
  out[37] = 5.4181235997219*tmp53*(nai_000_003_2*p_2 - nai_000_030_1*u_2); // final, weighs 6
  out[38] = 5.4181235997219*tmp53*(nai_000_020_3*p_1 + (nai_000_001_1 - nai_000_001_2)/ab_a - nai_000_100_0*u_1); // final, weighs 11
  out[39] = 2.4230585358948*tmp53*(tmp15 + nai_000_002_3*v_1 - nai_000_010_2*u_1); // final, weighs 7
  out[40] = 4.19686049388326*tmp53*(nai_000_300_0*v_1 - nai_000_100_2*u_1); // final, weighs 6
  nai_000_001_1 = nai_000_200_0*v_1 - nai_000_200_1*u_1; // local+recycle, weighs 4
  nai_000_001_2 = nai_000_200_1*v_1 - nai_000_200_2*u_1; // local+recycle, weighs 4
  out[41] = 9.3844653564063*tmp53*(nai_000_001_1*p_1 + (0.5*nai_000_200_0 - 0.5*nai_000_200_1)/ab_a - nai_000_001_2*u_1); // final, weighs 12
  out[42] = 9.3844653564063*tmp53*(nai_000_001_1*p_2 + (0.5*nai_010_200_0 - 0.5*nai_010_200_1)/ab_a - nai_000_001_2*u_2); // final, weighs 12
  nai_000_001_1 = d_2*u_1; // auto+recycle, weighs 1
  nai_000_001_2 = nai_000_000_2*u_1; // auto+recycle, weighs 1
  nai_000_002_3 = -nai_000_001_2 + d_2*p_1; // local+recycle, weighs 3
  nai_000_003_2 = nai_000_020_0*v_1 + (-nai_000_001_1 - nai_000_002_3 + d_1*p_1)/ab_a - nai_000_020_1*u_1; // local+recycle, weighs 12
  nai_000_010_1 = nai_000_000_3*u_1; // auto+recycle, weighs 1
  nai_000_002_3 = nai_000_020_1*v_1 + (nai_000_002_3 + nai_000_010_1 - nai_000_000_2*p_1)/ab_a - nai_000_020_2*u_1; // local+recycle, weighs 11
  out[43] = 9.3844653564063*tmp53*(nai_000_003_2*p_0 - nai_000_002_3*u_0); // final, weighs 6
  nai_000_001_2 = -nai_000_001_2 + d_2*v_1; // local+recycle, weighs 3
  nai_000_010_1 = -nai_000_010_1 + nai_000_000_2*v_1; // local+recycle, weighs 3
  nai_000_010_2 = nai_000_001_2*p_2 + (0.5*tmp73 - 0.5*tmp70)/ab_a - nai_000_010_1*u_2; // local+recycle, weighs 10
  out[44] = 16.2543707991657*tmp53*(p_0*(p_1*(p_2*(-nai_000_001_1 + d_1*v_1) + (0.5*tmp32 - 0.5*tmp73)/ab_a - nai_000_001_2*u_2) + (0.5*nai_100_002_0 - 0.5*nai_200_000_2)/ab_a - nai_000_010_2*u_1) - u_0*(nai_000_010_2*p_1 + (0.5*nai_200_000_2 - 0.5*tmp43)/ab_a - u_1*(nai_000_010_1*p_2 + (0.5*tmp70 - 0.5*tmp67)/ab_a - u_2*(nai_000_000_3*v_1 - d_0*u_1)))); // final, weighs 53
  nai_000_001_1 = nai_000_200_3*v_1 - nai_000_300_1*u_1; // local+recycle, weighs 4
  nai_000_001_2 = nai_000_300_1*v_1 - nai_000_001_3*u_1; // local+recycle, weighs 4
  out[45] = 9.3844653564063*tmp53*(nai_000_001_1*p_0 - nai_000_001_2*u_0); // final, weighs 6
  out[46] = 4.19686049388326*tmp53*(nai_000_030_0*v_1 + (1.5*nai_000_020_0 - 1.5*nai_000_020_1)/ab_a - nai_000_010_3*u_1); // final, weighs 12
  out[47] = 9.3844653564063*tmp53*(nai_110_200_0 + nai_000_003_2*p_2 - nai_000_002_3*u_2); // final, weighs 7
  out[48] = 9.3844653564063*tmp53*(nai_200_000_0 + nai_000_001_1*p_1 - nai_000_001_2*u_1); // final, weighs 7
  out[49] = 4.19686049388326*tmp53*(nai_000_002_0*v_1 - nai_000_002_1*u_1); // final, weighs 6
  out[50] = 2.4230585358948*tmp53*(nai_000_002_2 + nai_000_300_0*v_2 - nai_000_100_2*u_2); // final, weighs 7
  nai_000_000_4 = nai_000_000_4 + d_1*v_2; // local+recycle, weighs 2
  nai_000_001_1 = nai_100_002_1 + d_2*v_2; // local+recycle, weighs 2
  nai_000_001_2 = tmp31 + nai_000_000_2*v_2; // local+recycle, weighs 2
  nai_000_001_3 = nai_000_001_1*p_0 - nai_000_001_2*u_0; // local+recycle, weighs 4
  nai_000_002_2 = 0.5*nai_000_000_4 - 0.5*nai_000_001_1; // auto+recycle, weighs 3
  nai_000_002_3 = nai_000_002_2/ab_a; // auto+recycle, weighs 2
  nai_000_003_2 = nai_000_002_3 + p_0*(nai_000_000_4*p_0 - nai_000_001_1*u_0) - nai_000_001_3*u_0; // local+recycle, weighs 9
  d_0 = tmp3 + nai_000_000_3*v_2 - d_0*u_2; // local+recycle, weighs 5
  nai_000_010_1 = 0.5*nai_000_001_1 - 0.5*nai_000_001_2; // auto+recycle, weighs 3
  nai_000_010_2 = nai_000_010_1/ab_a; // auto+recycle, weighs 2
  nai_000_001_3 = nai_000_010_2 + nai_000_001_3*p_0 - u_0*(nai_000_001_2*p_0 - d_0*u_0); // local+recycle, weighs 9
  out[51] = 5.4181235997219*tmp53*(nai_000_003_2*p_1 - nai_000_001_3*u_1); // final, weighs 6
  out[52] = 5.4181235997219*tmp53*(nai_000_003_2*p_2 + (nai_000_200_0 - nai_000_200_1)/ab_a - nai_000_001_3*u_2); // final, weighs 11
  nai_000_001_3 = nai_000_001_1*p_1 - nai_000_001_2*u_1; // local+recycle, weighs 4
  nai_000_002_3 = nai_000_002_3 + p_1*(nai_000_000_4*p_1 - nai_000_001_1*u_1) - nai_000_001_3*u_1; // local+recycle, weighs 9
  nai_000_001_3 = nai_000_010_2 + nai_000_001_3*p_1 - u_1*(nai_000_001_2*p_1 - d_0*u_1); // local+recycle, weighs 9
  out[53] = 5.4181235997219*tmp53*(nai_000_002_3*p_0 - nai_000_001_3*u_0); // final, weighs 6
  d_1 = nai_000_000_4*p_2 + (d_1 - d_2)/ab_a - nai_000_001_1*u_2; // local+recycle, weighs 9
  d_2 = nai_000_001_1*p_2 + (d_2 - nai_000_000_2)/ab_a - nai_000_001_2*u_2; // local+recycle, weighs 9
  d_0 = nai_000_001_2*p_2 + (nai_000_000_2 - nai_000_000_3)/ab_a - d_0*u_2; // local+recycle, weighs 9
  out[54] = 9.3844653564063*tmp53*(p_0*(d_1*p_1 - d_2*u_1) - u_0*(d_2*p_1 - d_0*u_1)); // final, weighs 14
  d_1 = d_1*p_2 + (nai_000_002_2 + nai_100_002_0 - nai_200_000_2)/ab_a - d_2*u_2; // local+recycle, weighs 10
  d_0 = d_2*p_2 + (nai_000_010_1 + nai_200_000_2 - tmp43)/ab_a - d_0*u_2; // local+recycle, weighs 10
  out[55] = 5.4181235997219*tmp53*(d_1*p_0 - d_0*u_0); // final, weighs 6
  out[56] = 2.4230585358948*tmp53*(tmp5 + nai_000_030_0*v_2 - nai_000_010_3*u_2); // final, weighs 7
  out[57] = 5.4181235997219*tmp53*(nai_000_002_3*p_2 + (nai_000_020_0 - nai_000_020_1)/ab_a - nai_000_001_3*u_2); // final, weighs 11
  out[58] = 5.4181235997219*tmp53*(d_1*p_1 - d_0*u_1); // final, weighs 6
  out[59] = 2.4230585358948*tmp53*(nai_000_002_0*v_2 + (tmp38 + 1.5*nai_000_200_3 - 1.5*nai_000_300_1)/ab_a - nai_000_002_1*u_2); // final, weighs 13
  // total weight = 1830
}

static void gint2_nai_cF_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Number of local variables: 137
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_000_000_4, nai_000_000_5, nai_000_001_1, nai_000_001_2, nai_000_001_3, nai_000_001_4, nai_000_002_0, nai_000_002_1, nai_000_002_2, nai_000_002_3, nai_000_010_1, nai_000_010_2, nai_000_010_3, nai_000_010_4, nai_000_020_0, nai_000_020_1, nai_000_020_2, nai_000_020_3, nai_000_020_4, nai_010_001_1, nai_010_002_1, nai_010_020_0, nai_010_020_1, nai_010_020_2, nai_020_001_0, nai_020_001_1, nai_020_001_2, nai_020_010_2, nai_020_200_2, nai_100_000_0, nai_100_000_1, nai_100_000_2, nai_100_000_3, nai_100_000_4, nai_100_001_1, nai_100_030_1, nai_110_000_0, nai_110_000_1, nai_110_000_2, nai_110_000_3, nai_110_001_0, nai_110_001_1, nai_110_001_2, nai_110_020_0, nai_110_020_1, nai_120_001_1, nai_200_000_0, nai_200_000_1, nai_200_000_2, nai_200_000_3, nai_200_000_4, nai_200_001_0, nai_200_001_1, nai_200_001_2, nai_200_002_1, nai_200_003_0, nai_200_010_0, nai_200_010_1, nai_200_010_2, nai_200_020_0, nai_200_020_1, nai_200_030_0, nai_210_001_1, nai_300_000_0, nai_300_000_1, nai_300_000_2, nai_300_000_3, nai_300_100_0, nai_300_100_1, nai_300_200_0, nai_300_200_1, p_0, p_1, p_2, tmp0, tmp1, tmp106, tmp107, tmp108, tmp11, tmp110, tmp115, tmp12, tmp121, tmp122, tmp123, tmp124, tmp126, tmp129, tmp131, tmp133, tmp138, tmp139, tmp140, tmp157, tmp159, tmp160, tmp161, tmp162, tmp163, tmp164, tmp165, tmp166, tmp167, tmp169, tmp170, tmp173, tmp2, tmp20, tmp21, tmp3, tmp33, tmp34, tmp4, tmp62, tmp63, tmp69, tmp70, tmp74, tmp77, tmp8, tmp80, tmp87, tmp88, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a_a*a[0] + b_a*b[0])/ab_a; // local, weighs 5
  p_1 = (a_a*a[1] + b_a*b[1])/ab_a; // local, weighs 5
  p_2 = (a_a*a[2] + b_a*b[2])/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  v_0 = p_0 - a[0]; // local, weighs 2
  v_1 = p_1 - a[1]; // local, weighs 2
  v_2 = p_2 - a[2]; // local, weighs 2
  p_0 = p_0 - b[0]; // local+recycle, weighs 2
  p_1 = p_1 - b[1]; // local+recycle, weighs 2
  p_2 = p_2 - b[2]; // local+recycle, weighs 2
  d_0 = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_1 = 6.28318530717959*d_0*gaux(usq, 0); // local+recycle, weighs 3
  d_2 = 6.28318530717959*d_0*gaux(usq, 1); // local+recycle, weighs 3
  tmp173 = d_2*u_0; // auto, weighs 1
  nai_100_000_0 = -tmp173 + d_1*v_0; // local, weighs 3
  nai_000_000_2 = 6.28318530717959*d_0*gaux(usq, 2); // local, weighs 3
  tmp170 = nai_000_000_2*u_0; // auto, weighs 1
  nai_100_000_1 = -tmp170 + d_2*v_0; // local, weighs 3
  tmp3 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto, weighs 5
  tmp110 = tmp3 - nai_100_000_1*u_0; // auto, weighs 3
  nai_200_000_0 = tmp110 + nai_100_000_0*v_0; // local, weighs 2
  nai_000_000_3 = 6.28318530717959*d_0*gaux(usq, 3); // local, weighs 3
  tmp167 = nai_000_000_3*u_0; // auto, weighs 1
  nai_100_000_2 = -tmp167 + nai_000_000_2*v_0; // local, weighs 3
  tmp2 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto, weighs 5
  tmp63 = tmp2 - nai_100_000_2*u_0; // auto, weighs 3
  nai_200_000_1 = tmp63 + nai_100_000_1*v_0; // local, weighs 2
  tmp106 = (nai_100_000_0 - nai_100_000_1)/ab_a - nai_200_000_1*u_0; // auto, weighs 7
  nai_300_000_0 = tmp106 + nai_200_000_0*v_0; // local, weighs 2
  nai_000_000_4 = 6.28318530717959*d_0*gaux(usq, 4); // local, weighs 3
  tmp164 = nai_000_000_4*u_0; // auto, weighs 1
  nai_100_000_3 = -tmp164 + nai_000_000_3*v_0; // local, weighs 3
  tmp1 = (0.5*nai_000_000_2 - 0.5*nai_000_000_3)/ab_a; // auto, weighs 5
  tmp62 = tmp1 - nai_100_000_3*u_0; // auto, weighs 3
  nai_200_000_2 = tmp62 + nai_100_000_2*v_0; // local, weighs 2
  tmp107 = (nai_100_000_1 - nai_100_000_2)/ab_a - nai_200_000_2*u_0; // auto, weighs 7
  nai_300_000_1 = tmp107 + nai_200_000_1*v_0; // local, weighs 2
  nai_300_100_0 = nai_300_000_0*p_0 + (1.5*nai_200_000_0 - 1.5*nai_200_000_1)/ab_a - nai_300_000_1*u_0; // local, weighs 10
  nai_000_000_5 = 6.28318530717959*d_0*gaux(usq, 5); // local, weighs 3
  tmp161 = nai_000_000_5*u_0; // auto, weighs 1
  nai_100_000_4 = -tmp161 + nai_000_000_4*v_0; // local, weighs 3
  tmp0 = (0.5*nai_000_000_3 - 0.5*nai_000_000_4)/ab_a; // auto, weighs 5
  tmp115 = tmp0 - nai_100_000_4*u_0; // auto, weighs 3
  nai_200_000_3 = tmp115 + nai_100_000_3*v_0; // local, weighs 2
  tmp108 = (nai_100_000_2 - nai_100_000_3)/ab_a - nai_200_000_3*u_0; // auto, weighs 7
  nai_300_000_2 = tmp108 + nai_200_000_2*v_0; // local, weighs 2
  nai_300_100_1 = nai_300_000_1*p_0 + (1.5*nai_200_000_1 - 1.5*nai_200_000_2)/ab_a - nai_300_000_2*u_0; // local, weighs 10
  tmp106 = tmp106 + nai_200_000_0*p_0; // local+recycle, weighs 2
  tmp107 = tmp107 + nai_200_000_1*p_0; // local+recycle, weighs 2
  tmp70 = 0.5*nai_300_000_0 - 0.5*nai_300_000_1; // auto, weighs 3
  nai_300_200_0 = nai_300_100_0*p_0 + (tmp70 + 1.5*tmp106 - 1.5*tmp107)/ab_a - nai_300_100_1*u_0; // local, weighs 11
  d_0 = 6.28318530717959*d_0*gaux(usq, 6); // local+recycle, weighs 3
  usq = d_0*u_0; // auto+recycle, weighs 1
  tmp4 = (0.5*nai_000_000_4 - 0.5*nai_000_000_5)/ab_a; // auto, weighs 5
  nai_200_000_4 = tmp4 + nai_100_000_4*v_0 - u_0*(-usq + nai_000_000_5*v_0); // local, weighs 8
  nai_100_000_4 = (nai_100_000_3 - nai_100_000_4)/ab_a - nai_200_000_4*u_0; // auto+recycle, weighs 7
  nai_300_000_3 = nai_100_000_4 + nai_200_000_3*v_0; // local, weighs 2
  tmp108 = tmp108 + nai_200_000_2*p_0; // local+recycle, weighs 2
  tmp69 = 0.5*nai_300_000_1 - 0.5*nai_300_000_2; // auto, weighs 3
  nai_300_200_1 = nai_300_100_1*p_0 + (tmp69 + 1.5*tmp107 - 1.5*tmp108)/ab_a - u_0*(nai_300_000_2*p_0 + (1.5*nai_200_000_2 - 1.5*nai_200_000_3)/ab_a - nai_300_000_3*u_0); // local, weighs 21
  tmp63 = tmp63 + nai_100_000_1*p_0; // local+recycle, weighs 2
  tmp77 = 0.5*nai_200_000_0 - 0.5*nai_200_000_1; // auto, weighs 3
  tmp106 = p_0*tmp106 + (tmp110 + tmp77 - tmp63 + nai_100_000_0*p_0)/ab_a - tmp107*u_0; // local+recycle, weighs 12
  tmp110 = tmp62 + nai_100_000_2*p_0; // local+recycle, weighs 2
  tmp62 = 0.5*nai_200_000_1 - 0.5*nai_200_000_2; // auto+recycle, weighs 3
  tmp107 = p_0*tmp107 + (tmp62 + tmp63 - tmp110)/ab_a - tmp108*u_0; // local+recycle, weighs 10
  tmp63 = pow(a_a,2.25)*pow(b_a,2.25); // auto+recycle, weighs 3
  out[0] = 2.16724943988876*tmp63*(nai_300_200_0*p_0 + (nai_300_100_0 - nai_300_100_1 + 1.5*tmp106 - 1.5*tmp107)/ab_a - nai_300_200_1*u_0); // final, weighs 15
  out[1] = 4.84611707178961*tmp63*(nai_300_200_0*p_1 - nai_300_200_1*u_1); // final, weighs 6
  out[2] = 4.84611707178961*tmp63*(nai_300_200_0*p_2 - nai_300_200_1*u_2); // final, weighs 6
  nai_300_100_0 = nai_300_000_0*p_1 - nai_300_000_1*u_1; // local+recycle, weighs 4
  nai_300_100_1 = nai_300_000_1*p_1 - nai_300_000_2*u_1; // local+recycle, weighs 4
  nai_300_200_0 = tmp70/ab_a; // auto+recycle, weighs 2
  nai_300_200_1 = nai_300_200_0 + nai_300_100_0*p_1 - nai_300_100_1*u_1; // local+recycle, weighs 5
  tmp69 = tmp69/ab_a; // auto+recycle, weighs 2
  tmp70 = tmp69 + nai_300_100_1*p_1 - u_1*(nai_300_000_2*p_1 - nai_300_000_3*u_1); // local+recycle, weighs 9
  tmp126 = nai_200_000_1*u_1; // auto, weighs 1
  nai_200_010_0 = -tmp126 + nai_200_000_0*p_1; // local, weighs 3
  tmp124 = nai_200_000_2*u_1; // auto, weighs 1
  nai_200_010_1 = -tmp124 + nai_200_000_1*p_1; // local, weighs 3
  tmp77 = tmp77/ab_a; // auto+recycle, weighs 2
  nai_200_020_0 = tmp77 + nai_200_010_0*p_1 - nai_200_010_1*u_1; // local, weighs 5
  tmp122 = nai_200_000_3*u_1; // auto, weighs 1
  nai_200_010_2 = -tmp122 + nai_200_000_2*p_1; // local, weighs 3
  tmp62 = tmp62/ab_a; // auto+recycle, weighs 2
  nai_200_020_1 = tmp62 + nai_200_010_1*p_1 - nai_200_010_2*u_1; // local, weighs 5
  tmp11 = (1.5*nai_200_020_0 - 1.5*nai_200_020_1)/ab_a; // auto, weighs 5
  out[3] = 4.84611707178961*tmp63*(tmp11 + nai_300_200_1*p_0 - tmp70*u_0); // final, weighs 7
  nai_300_000_0 = nai_300_000_0*p_2 - nai_300_000_1*u_2; // local+recycle, weighs 4
  nai_300_000_1 = nai_300_000_1*p_2 - nai_300_000_2*u_2; // local+recycle, weighs 4
  nai_300_000_2 = nai_300_000_2*p_2 - nai_300_000_3*u_2; // local+recycle, weighs 4
  nai_300_000_3 = nai_200_000_1*u_2; // auto+recycle, weighs 1
  nai_200_001_0 = -nai_300_000_3 + nai_200_000_0*p_2; // local, weighs 3
  tmp123 = nai_200_000_2*u_2; // auto, weighs 1
  nai_200_001_1 = -tmp123 + nai_200_000_1*p_2; // local, weighs 3
  tmp121 = nai_200_000_3*u_2; // auto, weighs 1
  nai_200_001_2 = -tmp121 + nai_200_000_2*p_2; // local, weighs 3
  out[4] = 8.39372098776651*tmp63*(p_0*(nai_300_000_0*p_1 - nai_300_000_1*u_1) + (1.5*nai_200_001_0*p_1 + 1.5*nai_200_001_2*u_1 - 1.5*nai_200_001_1*p_1 - 1.5*nai_200_001_1*u_1)/ab_a - u_0*(nai_300_000_1*p_1 - nai_300_000_2*u_1)); // final, weighs 28
  nai_300_200_0 = nai_300_200_0 + nai_300_000_0*p_2 - nai_300_000_1*u_2; // local+recycle, weighs 5
  nai_300_000_2 = tmp69 + nai_300_000_1*p_2 - nai_300_000_2*u_2; // local+recycle, weighs 5
  tmp69 = tmp77 + nai_200_001_0*p_2 - nai_200_001_1*u_2; // local+recycle, weighs 5
  nai_200_002_1 = tmp62 + nai_200_001_1*p_2 - nai_200_001_2*u_2; // local, weighs 5
  tmp12 = (1.5*tmp69 - 1.5*nai_200_002_1)/ab_a; // auto, weighs 5
  out[5] = 4.84611707178961*tmp63*(tmp12 + nai_300_200_0*p_0 - nai_300_000_2*u_0); // final, weighs 7
  out[6] = 2.16724943988876*tmp63*(nai_300_200_1*p_1 + (nai_300_100_0 - nai_300_100_1)/ab_a - tmp70*u_1); // final, weighs 11
  out[7] = 4.84611707178961*tmp63*(nai_300_200_1*p_2 - tmp70*u_2); // final, weighs 6
  out[8] = 4.84611707178961*tmp63*(nai_300_200_0*p_1 - nai_300_000_2*u_1); // final, weighs 6
  out[9] = 2.16724943988876*tmp63*(nai_300_200_0*p_2 + (nai_300_000_0 - nai_300_000_1)/ab_a - nai_300_000_2*u_2); // final, weighs 11
  nai_300_000_0 = -tmp173 + d_1*p_0; // local+recycle, weighs 3
  nai_300_000_1 = -tmp170 + d_2*p_0; // local+recycle, weighs 3
  nai_300_000_2 = tmp3 + nai_300_000_0*p_0 - nai_300_000_1*u_0; // local+recycle, weighs 5
  nai_300_100_0 = -tmp167 + nai_000_000_2*p_0; // local+recycle, weighs 3
  nai_300_100_1 = tmp2 + nai_300_000_1*p_0 - nai_300_100_0*u_0; // local+recycle, weighs 5
  nai_300_000_0 = (nai_300_000_0 - nai_300_000_1)/ab_a - nai_300_100_1*u_0; // auto+recycle, weighs 7
  nai_300_200_0 = nai_300_000_0 + nai_300_000_2*p_0; // local+recycle, weighs 2
  nai_300_200_1 = -tmp164 + nai_000_000_3*p_0; // local+recycle, weighs 3
  tmp164 = tmp1 + nai_300_100_0*p_0 - nai_300_200_1*u_0; // local+recycle, weighs 5
  nai_300_000_1 = (nai_300_000_1 - nai_300_100_0)/ab_a - tmp164*u_0; // auto+recycle, weighs 7
  tmp167 = nai_300_000_1 + nai_300_100_1*p_0; // local+recycle, weighs 2
  tmp161 = -tmp161 + nai_000_000_4*p_0; // local+recycle, weighs 3
  tmp170 = tmp0 + nai_300_200_1*p_0 - tmp161*u_0; // local+recycle, weighs 5
  nai_300_100_0 = (nai_300_100_0 - nai_300_200_1)/ab_a - tmp170*u_0; // auto+recycle, weighs 7
  tmp173 = nai_300_100_0 + p_0*tmp164; // local+recycle, weighs 2
  tmp70 = tmp167*v_0 + (1.5*nai_300_100_1 - 1.5*tmp164)/ab_a - tmp173*u_0; // local+recycle, weighs 10
  nai_300_000_1 = nai_300_000_1 + nai_300_100_1*v_0; // local+recycle, weighs 2
  tmp87 = 0.5*nai_300_200_0 - 0.5*tmp167; // auto, weighs 3
  nai_300_000_0 = v_0*(nai_300_200_0*v_0 + (1.5*nai_300_000_2 - 1.5*nai_300_100_1)/ab_a - tmp167*u_0) + (tmp87 + 1.5*nai_300_000_0 - 1.5*nai_300_000_1 + 1.5*nai_300_000_2*v_0)/ab_a - tmp70*u_0; // local+recycle, weighs 24
  usq = tmp4 + p_0*tmp161 - u_0*(-usq + nai_000_000_5*p_0); // local+recycle, weighs 8
  nai_300_200_1 = p_0*tmp170 + (nai_300_200_1 - tmp161)/ab_a - u_0*usq; // local+recycle, weighs 9
  tmp161 = 0.5*tmp167 - 0.5*tmp173; // auto+recycle, weighs 3
  nai_300_000_1 = tmp70*v_0 + (tmp161 + 1.5*nai_300_000_1 - 1.5*nai_300_100_0 - 1.5*tmp164*v_0)/ab_a - u_0*(tmp173*v_0 + (1.5*tmp164 - 1.5*tmp170)/ab_a - nai_300_200_1*u_0); // local+recycle, weighs 24
  out[10] = 4.84611707178961*tmp63*(nai_300_000_0*v_1 - nai_300_000_1*u_1); // final, weighs 6
  nai_300_100_0 = tmp106*v_1 - tmp107*u_1; // local+recycle, weighs 4
  tmp70 = 0.5*nai_200_000_2 - 0.5*nai_200_000_3; // auto+recycle, weighs 3
  nai_100_000_4 = p_0*tmp108 + (tmp110 + tmp70 - tmp115 - nai_100_000_3*p_0)/ab_a - u_0*(nai_100_000_4 + nai_200_000_3*p_0); // local+recycle, weighs 15
  tmp108 = tmp107*v_1 - nai_100_000_4*u_1; // local+recycle, weighs 4
  tmp110 = (0.5*tmp106 - 0.5*tmp107)/ab_a; // auto+recycle, weighs 5
  out[11] = 10.8362471994438*tmp63*(tmp110 + nai_300_100_0*p_1 - tmp108*u_1); // final, weighs 7
  out[12] = 10.8362471994438*tmp63*(nai_300_100_0*p_2 - tmp108*u_2); // final, weighs 6
  nai_200_010_0 = nai_200_020_0*v_1 + (nai_200_010_0 - nai_200_010_1)/ab_a - nai_200_020_1*u_1; // local+recycle, weighs 9
  nai_300_100_0 = nai_200_000_4*u_1; // auto+recycle, weighs 1
  tmp108 = tmp70/ab_a; // auto+recycle, weighs 2
  tmp115 = tmp108 + nai_200_010_2*p_1 - u_1*(-nai_300_100_0 + nai_200_000_3*p_1); // local+recycle, weighs 8
  nai_200_010_1 = nai_200_020_1*v_1 + (nai_200_010_1 - nai_200_010_2)/ab_a - tmp115*u_1; // local+recycle, weighs 9
  nai_200_010_2 = d_2*u_1; // auto+recycle, weighs 1
  tmp70 = -nai_200_010_2 + d_1*p_1; // local+recycle, weighs 3
  tmp169 = nai_000_000_2*u_1; // auto, weighs 1
  nai_000_010_1 = -tmp169 + d_2*p_1; // local, weighs 3
  nai_000_020_0 = tmp3 + p_1*tmp70 - nai_000_010_1*u_1; // local, weighs 5
  tmp166 = nai_000_000_3*u_1; // auto, weighs 1
  nai_000_010_2 = -tmp166 + nai_000_000_2*p_1; // local, weighs 3
  nai_000_020_1 = tmp2 + nai_000_010_1*p_1 - nai_000_010_2*u_1; // local, weighs 5
  tmp70 = (tmp70 - nai_000_010_1)/ab_a - nai_000_020_1*u_1; // auto+recycle, weighs 7
  nai_010_020_0 = tmp70 + nai_000_020_0*v_1; // local, weighs 2
  tmp163 = nai_000_000_4*u_1; // auto, weighs 1
  nai_000_010_3 = -tmp163 + nai_000_000_3*p_1; // local, weighs 3
  nai_000_020_2 = tmp1 + nai_000_010_2*p_1 - nai_000_010_3*u_1; // local, weighs 5
  nai_000_010_1 = (nai_000_010_1 - nai_000_010_2)/ab_a - nai_000_020_2*u_1; // auto+recycle, weighs 7
  nai_010_020_1 = nai_000_010_1 + nai_000_020_1*v_1; // local, weighs 2
  nai_110_020_0 = nai_010_020_0*v_0 - nai_010_020_1*u_0; // local, weighs 4
  tmp160 = nai_000_000_5*u_1; // auto, weighs 1
  nai_000_010_4 = -tmp160 + nai_000_000_4*p_1; // local, weighs 3
  nai_000_020_3 = tmp0 + nai_000_010_3*p_1 - nai_000_010_4*u_1; // local, weighs 5
  nai_000_010_2 = (nai_000_010_2 - nai_000_010_3)/ab_a - nai_000_020_3*u_1; // auto+recycle, weighs 7
  nai_010_020_2 = nai_000_010_2 + nai_000_020_2*v_1; // local, weighs 2
  nai_110_020_1 = nai_010_020_1*v_0 - nai_010_020_2*u_0; // local, weighs 4
  out[13] = 10.8362471994438*tmp63*(nai_200_010_0*p_0 + (nai_110_020_0 - nai_110_020_1)/ab_a - nai_200_010_1*u_0); // final, weighs 11
  tmp124 = -tmp124 + nai_200_000_1*v_1; // local+recycle, weighs 3
  tmp122 = -tmp122 + nai_200_000_2*v_1; // local+recycle, weighs 3
  nai_210_001_1 = p_2*tmp124 - tmp122*u_2; // local, weighs 4
  nai_200_010_2 = -nai_200_010_2 + d_1*v_1; // local+recycle, weighs 3
  tmp169 = -tmp169 + d_2*v_1; // local+recycle, weighs 3
  tmp140 = tmp169*u_0; // auto, weighs 1
  nai_110_000_0 = -tmp140 + nai_200_010_2*v_0; // local, weighs 3
  tmp166 = -tmp166 + nai_000_000_2*v_1; // local+recycle, weighs 3
  tmp139 = tmp166*u_0; // auto, weighs 1
  nai_110_000_1 = -tmp139 + tmp169*v_0; // local, weighs 3
  nai_110_001_0 = nai_110_000_0*p_2 - nai_110_000_1*u_2; // local, weighs 4
  tmp163 = -tmp163 + nai_000_000_3*v_1; // local+recycle, weighs 3
  tmp138 = tmp163*u_0; // auto, weighs 1
  nai_110_000_2 = -tmp138 + tmp166*v_0; // local, weighs 3
  nai_110_001_1 = nai_110_000_1*p_2 - nai_110_000_2*u_2; // local, weighs 4
  nai_100_001_1 = nai_100_000_1*p_2 - nai_100_000_2*u_2; // local, weighs 4
  tmp160 = -tmp160 + nai_000_000_4*v_1; // local+recycle, weighs 3
  nai_110_000_3 = tmp163*v_0 - tmp160*u_0; // local, weighs 4
  nai_110_001_2 = nai_110_000_2*p_2 - nai_110_000_3*u_2; // local, weighs 4
  out[14] = 18.7689307128126*tmp63*(p_0*(p_1*(p_2*(-tmp126 + nai_200_000_0*v_1) - tmp124*u_2) + (0.5*nai_200_001_0 - 0.5*nai_200_001_1)/ab_a - nai_210_001_1*u_1) + (nai_110_001_0*p_1 + nai_110_001_2*u_1 + (-0.5*nai_100_001_1 + 0.5*nai_100_000_0*p_2 - 0.5*nai_100_000_1*u_2)/ab_a - nai_110_001_1*p_1 - nai_110_001_1*u_1 - (0.5*nai_100_001_1 + 0.5*nai_100_000_3*u_2 - 0.5*nai_100_000_2*p_2)/ab_a)/ab_a - u_0*(nai_210_001_1*p_1 + (0.5*nai_200_001_1 - 0.5*nai_200_001_2)/ab_a - u_1*(p_2*tmp122 - u_2*(-nai_300_100_0 + nai_200_000_3*v_1)))); // final, weighs 73
  nai_100_001_1 = tmp69*v_1 - nai_200_002_1*u_1; // local+recycle, weighs 4
  nai_200_000_4 = nai_200_000_4*u_2; // auto+recycle, weighs 1
  nai_210_001_1 = tmp108 + nai_200_001_2*p_2 - u_2*(-nai_200_000_4 + nai_200_000_3*p_2); // local+recycle, weighs 8
  nai_300_100_0 = nai_200_002_1*v_1 - nai_210_001_1*u_1; // local+recycle, weighs 4
  tmp122 = d_2*u_2; // auto+recycle, weighs 1
  tmp124 = -tmp122 + d_1*p_2; // local+recycle, weighs 3
  tmp126 = nai_000_000_2*u_2; // auto+recycle, weighs 1
  nai_000_001_1 = -tmp126 + d_2*p_2; // local, weighs 3
  nai_000_002_0 = tmp3 + p_2*tmp124 - nai_000_001_1*u_2; // local, weighs 5
  tmp165 = nai_000_000_3*u_2; // auto, weighs 1
  nai_000_001_2 = -tmp165 + nai_000_000_2*p_2; // local, weighs 3
  nai_000_002_1 = tmp2 + nai_000_001_1*p_2 - nai_000_001_2*u_2; // local, weighs 5
  tmp162 = nai_000_000_4*u_2; // auto, weighs 1
  nai_000_001_3 = -tmp162 + nai_000_000_3*p_2; // local, weighs 3
  nai_000_002_2 = tmp1 + nai_000_001_2*p_2 - nai_000_001_3*u_2; // local, weighs 5
  nai_010_002_1 = nai_000_002_1*v_1 - nai_000_002_2*u_1; // local, weighs 4
  tmp159 = nai_000_000_5*u_2; // auto, weighs 1
  nai_000_001_4 = -tmp159 + nai_000_000_4*p_2; // local, weighs 3
  nai_000_002_3 = tmp0 + nai_000_001_3*p_2 - nai_000_001_4*u_2; // local, weighs 5
  nai_010_002_1 = (u_0*(nai_000_002_2*v_1 - nai_000_002_3*u_1) + v_0*(nai_000_002_0*v_1 - nai_000_002_1*u_1) - nai_010_002_1*u_0 - nai_010_002_1*v_0)/ab_a; // auto+recycle, weighs 19
  out[15] = 10.8362471994438*tmp63*(nai_010_002_1 + nai_100_001_1*p_0 - nai_300_100_0*u_0); // final, weighs 7
  tmp70 = tmp70 + nai_000_020_0*p_1; // local+recycle, weighs 2
  nai_000_010_1 = nai_000_010_1 + nai_000_020_1*p_1; // local+recycle, weighs 2
  nai_000_010_2 = nai_000_010_2 + nai_000_020_2*p_1; // local+recycle, weighs 2
  nai_100_030_1 = nai_000_010_1*v_0 - nai_000_010_2*u_0; // local, weighs 4
  tmp74 = 0.5*tmp70 - 0.5*nai_000_010_1; // auto, weighs 3
  tmp21 = tmp74/ab_a; // auto, weighs 2
  nai_200_030_0 = tmp21 + v_0*(tmp70*v_0 - nai_000_010_1*u_0) - nai_100_030_1*u_0; // local, weighs 9
  tmp157 = d_0*u_1; // auto, weighs 1
  nai_000_020_4 = tmp4 + nai_000_010_4*p_1 - u_1*(-tmp157 + nai_000_000_5*p_1); // local, weighs 8
  nai_000_010_3 = nai_000_020_3*p_1 + (nai_000_010_3 - nai_000_010_4)/ab_a - nai_000_020_4*u_1; // local+recycle, weighs 9
  nai_000_010_4 = 0.5*nai_000_010_1 - 0.5*nai_000_010_2; // auto+recycle, weighs 3
  tmp20 = nai_000_010_4/ab_a; // auto, weighs 2
  nai_100_030_1 = tmp20 + nai_100_030_1*v_0 - u_0*(nai_000_010_2*v_0 - nai_000_010_3*u_0); // local+recycle, weighs 9
  out[16] = 4.84611707178961*tmp63*(tmp11 + nai_200_030_0*v_1 - nai_100_030_1*u_1); // final, weighs 7
  out[17] = 10.8362471994438*tmp63*(nai_200_010_0*p_2 - nai_200_010_1*u_2); // final, weighs 6
  out[18] = 10.8362471994438*tmp63*(nai_100_001_1*p_1 + (0.5*tmp69 - 0.5*nai_200_002_1)/ab_a - nai_300_100_0*u_1); // final, weighs 12
  nai_100_001_1 = (tmp124 - nai_000_001_1)/ab_a - nai_000_002_1*u_2; // auto+recycle, weighs 7
  nai_200_010_0 = nai_100_001_1 + nai_000_002_0*p_2; // local+recycle, weighs 2
  nai_000_001_1 = (nai_000_001_1 - nai_000_001_2)/ab_a - nai_000_002_2*u_2; // auto+recycle, weighs 7
  nai_200_010_1 = nai_000_001_1 + nai_000_002_1*p_2; // local+recycle, weighs 2
  nai_000_001_2 = (nai_000_001_2 - nai_000_001_3)/ab_a - nai_000_002_3*u_2; // auto+recycle, weighs 7
  nai_300_100_0 = nai_000_001_2 + nai_000_002_2*p_2; // local+recycle, weighs 2
  tmp11 = nai_200_010_1*v_0 - nai_300_100_0*u_0; // local+recycle, weighs 4
  tmp124 = 0.5*nai_200_010_0 - 0.5*nai_200_010_1; // auto+recycle, weighs 3
  tmp34 = tmp124/ab_a; // auto, weighs 2
  nai_200_003_0 = tmp34 + v_0*(nai_200_010_0*v_0 - nai_200_010_1*u_0) - tmp11*u_0; // local, weighs 9
  d_0 = d_0*u_2; // auto+recycle, weighs 1
  nai_000_001_3 = (nai_000_001_3 - nai_000_001_4)/ab_a - u_2*(tmp4 + nai_000_001_4*p_2 - u_2*(-d_0 + nai_000_000_5*p_2)); // auto+recycle, weighs 15
  nai_000_001_4 = nai_000_001_3 + nai_000_002_3*p_2; // local+recycle, weighs 2
  tmp88 = 0.5*nai_200_010_1 - 0.5*nai_300_100_0; // auto, weighs 3
  tmp33 = tmp88/ab_a; // auto, weighs 2
  tmp11 = tmp33 + tmp11*v_0 - u_0*(nai_300_100_0*v_0 - nai_000_001_4*u_0); // local+recycle, weighs 9
  out[19] = 4.84611707178961*tmp63*(nai_200_003_0*v_1 - tmp11*u_1); // final, weighs 6
  out[20] = 4.84611707178961*tmp63*(nai_300_000_0*v_2 - nai_300_000_1*u_2); // final, weighs 6
  nai_300_000_0 = tmp106*v_2 - tmp107*u_2; // local+recycle, weighs 4
  nai_100_000_4 = tmp107*v_2 - nai_100_000_4*u_2; // local+recycle, weighs 4
  out[21] = 10.8362471994438*tmp63*(nai_300_000_0*p_1 - nai_100_000_4*u_1); // final, weighs 6
  out[22] = 10.8362471994438*tmp63*(tmp110 + nai_300_000_0*p_2 - nai_100_000_4*u_2); // final, weighs 7
  nai_100_000_4 = nai_200_020_0*v_2 - nai_200_020_1*u_2; // local+recycle, weighs 4
  nai_300_000_0 = nai_200_020_1*v_2 - tmp115*u_2; // local+recycle, weighs 4
  nai_300_000_1 = nai_000_020_0*v_2 - nai_000_020_1*u_2; // local+recycle, weighs 4
  tmp106 = nai_000_020_1*v_2 - nai_000_020_2*u_2; // local+recycle, weighs 4
  tmp107 = nai_000_020_2*v_2 - nai_000_020_3*u_2; // local+recycle, weighs 4
  tmp110 = (nai_300_000_1*v_0 + tmp107*u_0 - tmp106*u_0 - tmp106*v_0)/ab_a; // auto+recycle, weighs 11
  out[23] = 10.8362471994438*tmp63*(tmp110 + nai_100_000_4*p_0 - nai_300_000_0*u_0); // final, weighs 7
  nai_200_000_1 = -tmp123 + nai_200_000_1*v_2; // local+recycle, weighs 3
  nai_200_000_2 = -tmp121 + nai_200_000_2*v_2; // local+recycle, weighs 3
  tmp115 = tmp62 + nai_200_000_1*p_2 - nai_200_000_2*u_2; // local+recycle, weighs 5
  d_1 = -tmp122 + d_1*v_2; // local+recycle, weighs 3
  d_2 = -tmp126 + d_2*v_2; // local+recycle, weighs 3
  tmp121 = d_2*u_0; // auto+recycle, weighs 1
  tmp122 = -tmp121 + d_1*v_0; // local+recycle, weighs 3
  nai_000_000_2 = -tmp165 + nai_000_000_2*v_2; // local+recycle, weighs 3
  tmp123 = nai_000_000_2*u_0; // auto+recycle, weighs 1
  tmp126 = -tmp123 + d_2*v_0; // local+recycle, weighs 3
  nai_100_000_0 = p_2*tmp122 + (0.5*nai_100_000_0 - 0.5*nai_100_000_1)/ab_a - tmp126*u_2; // local+recycle, weighs 10
  nai_000_000_3 = -tmp162 + nai_000_000_3*v_2; // local+recycle, weighs 3
  tmp162 = nai_000_000_3*u_0; // auto+recycle, weighs 1
  tmp165 = -tmp162 + nai_000_000_2*v_0; // local+recycle, weighs 3
  nai_100_000_1 = p_2*tmp126 + (0.5*nai_100_000_1 - 0.5*nai_100_000_2)/ab_a - tmp165*u_2; // local+recycle, weighs 10
  nai_000_000_4 = -tmp159 + nai_000_000_4*v_2; // local+recycle, weighs 3
  tmp159 = nai_000_000_3*v_0 - nai_000_000_4*u_0; // local+recycle, weighs 4
  nai_100_000_2 = p_2*tmp165 + (0.5*nai_100_000_2 - 0.5*nai_100_000_3)/ab_a - tmp159*u_2; // local+recycle, weighs 10
  out[24] = 18.7689307128126*tmp63*(p_0*(p_1*(tmp77 + p_2*(-nai_300_000_3 + nai_200_000_0*v_2) - nai_200_000_1*u_2) - tmp115*u_1) + (nai_100_000_0*p_1 + nai_100_000_2*u_1 - nai_100_000_1*p_1 - nai_100_000_1*u_1)/ab_a - u_0*(p_1*tmp115 - u_1*(tmp108 + nai_200_000_2*p_2 - u_2*(-nai_200_000_4 + nai_200_000_3*v_2)))); // final, weighs 42
  nai_100_000_3 = tmp69*v_2 + (nai_200_001_0 - nai_200_001_1)/ab_a - nai_200_002_1*u_2; // local+recycle, weighs 9
  nai_200_000_0 = nai_200_002_1*v_2 + (nai_200_001_1 - nai_200_001_2)/ab_a - nai_210_001_1*u_2; // local+recycle, weighs 9
  nai_100_001_1 = nai_100_001_1 + nai_000_002_0*v_2; // local+recycle, weighs 2
  nai_000_001_1 = nai_000_001_1 + nai_000_002_1*v_2; // local+recycle, weighs 2
  nai_200_000_1 = nai_100_001_1*v_0 - nai_000_001_1*u_0; // local+recycle, weighs 4
  nai_000_001_2 = nai_000_001_2 + nai_000_002_2*v_2; // local+recycle, weighs 2
  nai_200_000_2 = nai_000_001_1*v_0 - nai_000_001_2*u_0; // local+recycle, weighs 4
  out[25] = 10.8362471994438*tmp63*(nai_100_000_3*p_0 + (nai_200_000_1 - nai_200_000_2)/ab_a - nai_200_000_0*u_0); // final, weighs 11
  out[26] = 4.84611707178961*tmp63*(nai_200_030_0*v_2 - nai_100_030_1*u_2); // final, weighs 6
  out[27] = 10.8362471994438*tmp63*(nai_100_000_4*p_2 + (0.5*nai_200_020_0 - 0.5*nai_200_020_1)/ab_a - nai_300_000_0*u_2); // final, weighs 12
  out[28] = 10.8362471994438*tmp63*(nai_100_000_3*p_1 - nai_200_000_0*u_1); // final, weighs 6
  out[29] = 4.84611707178961*tmp63*(tmp12 + nai_200_003_0*v_2 - tmp11*u_2); // final, weighs 7
  nai_100_000_3 = tmp167*v_1 - tmp173*u_1; // local+recycle, weighs 4
  nai_100_000_4 = tmp87/ab_a; // auto+recycle, weighs 2
  nai_100_030_1 = nai_100_000_4 + v_1*(nai_300_200_0*v_1 - tmp167*u_1) - nai_100_000_3*u_1; // local+recycle, weighs 9
  nai_200_000_0 = tmp161/ab_a; // auto+recycle, weighs 2
  nai_100_000_3 = nai_200_000_0 + nai_100_000_3*v_1 - u_1*(tmp173*v_1 - nai_300_200_1*u_1); // local+recycle, weighs 9
  nai_200_000_3 = tmp3 - tmp169*u_1; // auto+recycle, weighs 3
  nai_200_000_4 = nai_200_000_3 + nai_200_010_2*v_1; // local+recycle, weighs 2
  nai_200_001_0 = tmp2 - tmp166*u_1; // auto+recycle, weighs 3
  nai_200_001_1 = nai_200_001_0 + tmp169*v_1; // local+recycle, weighs 2
  nai_200_001_2 = nai_200_001_1*u_0; // auto+recycle, weighs 1
  nai_200_002_1 = -nai_200_001_2 + nai_200_000_4*p_0; // local+recycle, weighs 3
  nai_200_003_0 = tmp1 - tmp163*u_1; // auto+recycle, weighs 3
  nai_200_020_0 = nai_200_003_0 + tmp166*v_1; // local+recycle, weighs 2
  nai_200_020_1 = nai_200_020_0*u_0; // auto+recycle, weighs 1
  nai_200_030_0 = -nai_200_020_1 + nai_200_001_1*p_0; // local+recycle, weighs 3
  nai_210_001_1 = 0.5*nai_200_000_4 - 0.5*nai_200_001_1; // auto+recycle, weighs 3
  nai_300_000_0 = nai_210_001_1/ab_a; // auto+recycle, weighs 2
  nai_300_000_3 = nai_300_000_0 + nai_200_002_1*p_0 - nai_200_030_0*u_0; // local+recycle, weighs 5
  tmp108 = tmp0 - tmp160*u_1; // auto+recycle, weighs 3
  tmp11 = tmp108 + tmp163*v_1; // local+recycle, weighs 2
  tmp115 = tmp11*u_0; // auto+recycle, weighs 1
  tmp12 = -tmp115 + nai_200_020_0*p_0; // local+recycle, weighs 3
  tmp161 = 0.5*nai_200_001_1 - 0.5*nai_200_020_0; // auto+recycle, weighs 3
  tmp62 = tmp161/ab_a; // auto+recycle, weighs 2
  tmp69 = tmp62 + nai_200_030_0*p_0 - tmp12*u_0; // local+recycle, weighs 5
  tmp77 = (1.5*nai_300_000_3 - 1.5*tmp69)/ab_a; // auto+recycle, weighs 5
  out[30] = 4.84611707178961*tmp63*(tmp77 + nai_100_030_1*v_0 - nai_100_000_3*u_0); // final, weighs 7
  nai_200_002_1 = nai_300_000_3*v_0 + (nai_200_002_1 - nai_200_030_0)/ab_a - tmp69*u_0; // local+recycle, weighs 9
  tmp157 = tmp4 + tmp160*v_1 - u_1*(-tmp157 + nai_000_000_5*v_1); // local+recycle, weighs 8
  tmp87 = tmp157*u_0; // auto+recycle, weighs 1
  tmp80 = 0.5*nai_200_020_0 - 0.5*tmp11; // auto, weighs 3
  tmp8 = tmp80/ab_a; // auto, weighs 2
  nai_020_200_2 = tmp8 + p_0*tmp12 - u_0*(-tmp87 + p_0*tmp11); // local, weighs 8
  nai_200_030_0 = tmp69*v_0 + (nai_200_030_0 - tmp12)/ab_a - nai_020_200_2*u_0; // local+recycle, weighs 9
  tmp12 = nai_300_100_1*v_1 - tmp164*u_1; // local+recycle, weighs 4
  tmp139 = -tmp139 + p_0*tmp169; // local+recycle, weighs 3
  tmp140 = v_0*(nai_300_000_2*v_1 - nai_300_100_1*u_1) + (-tmp139 - tmp140 + nai_200_010_2*p_0)/ab_a - tmp12*u_0; // local+recycle, weighs 16
  tmp12 = tmp12*v_0 + (tmp138 + tmp139 - p_0*tmp166)/ab_a - u_0*(tmp164*v_1 - tmp170*u_1); // local+recycle, weighs 15
  out[31] = 10.8362471994438*tmp63*(nai_200_002_1*p_1 + (tmp140 - tmp12)/ab_a - nai_200_030_0*u_1); // final, weighs 11
  out[32] = 10.8362471994438*tmp63*(nai_200_002_1*p_2 - nai_200_030_0*u_2); // final, weighs 6
  nai_200_002_1 = (nai_200_010_2 - tmp169)/ab_a - nai_200_001_1*u_1; // auto+recycle, weighs 7
  nai_200_030_0 = nai_200_002_1 + nai_200_000_4*p_1; // local+recycle, weighs 2
  tmp138 = (tmp169 - tmp166)/ab_a - nai_200_020_0*u_1; // auto+recycle, weighs 7
  tmp139 = tmp138 + nai_200_001_1*p_1; // local+recycle, weighs 2
  nai_200_001_0 = nai_200_001_0 + p_1*tmp169; // local+recycle, weighs 2
  nai_200_000_3 = nai_200_030_0*p_1 + (nai_200_000_3 + nai_210_001_1 - nai_200_001_0 + nai_200_010_2*p_1)/ab_a - tmp139*u_1; // local+recycle, weighs 12
  nai_210_001_1 = (tmp166 - tmp163)/ab_a - tmp11*u_1; // auto+recycle, weighs 7
  nai_020_010_2 = nai_210_001_1 + nai_200_020_0*p_1; // local, weighs 2
  nai_200_003_0 = nai_200_003_0 + p_1*tmp166; // local+recycle, weighs 2
  nai_200_001_0 = p_1*tmp139 + (nai_200_001_0 + tmp161 - nai_200_003_0)/ab_a - nai_020_010_2*u_1; // local+recycle, weighs 10
  tmp161 = nai_200_000_3*v_0 - nai_200_001_0*u_0; // local+recycle, weighs 4
  tmp160 = (tmp163 - tmp160)/ab_a - tmp157*u_1; // auto+recycle, weighs 7
  nai_200_003_0 = nai_020_010_2*p_1 + (nai_200_003_0 + tmp80 - tmp108 - p_1*tmp163)/ab_a - u_1*(tmp160 + p_1*tmp11); // local+recycle, weighs 15
  tmp108 = nai_200_001_0*v_0 - nai_200_003_0*u_0; // local+recycle, weighs 4
  tmp80 = (0.5*nai_200_000_3 - 0.5*nai_200_001_0)/ab_a; // auto+recycle, weighs 5
  out[33] = 10.8362471994438*tmp63*(tmp80 + p_0*tmp161 - tmp108*u_0); // final, weighs 7
  nai_200_020_1 = -nai_200_020_1 + nai_200_001_1*v_0; // local+recycle, weighs 3
  tmp115 = -tmp115 + nai_200_020_0*v_0; // local+recycle, weighs 3
  nai_120_001_1 = nai_200_020_1*p_2 - tmp115*u_2; // local, weighs 4
  tmp133 = nai_200_001_1*u_2; // auto, weighs 1
  nai_020_001_0 = -tmp133 + nai_200_000_4*p_2; // local, weighs 3
  tmp131 = nai_200_020_0*u_2; // auto, weighs 1
  nai_020_001_1 = -tmp131 + nai_200_001_1*p_2; // local, weighs 3
  nai_010_001_1 = p_2*tmp169 - tmp166*u_2; // local, weighs 4
  tmp129 = tmp11*u_2; // auto, weighs 1
  nai_020_001_2 = -tmp129 + nai_200_020_0*p_2; // local, weighs 3
  out[34] = 18.7689307128126*tmp63*(p_0*(p_1*(p_2*(-nai_200_001_2 + nai_200_000_4*v_0) - nai_200_020_1*u_2) + (nai_110_001_0 - nai_110_001_1)/ab_a - nai_120_001_1*u_1) + (0.5*nai_020_001_0*p_1 + 0.5*nai_020_001_2*u_1 + 0.5*(-nai_010_001_1 + nai_200_010_2*p_2 - tmp169*u_2)/ab_a - 0.5*nai_020_001_1*p_1 - 0.5*nai_020_001_1*u_1 - 0.5*(nai_010_001_1 + tmp163*u_2 - p_2*tmp166)/ab_a)/ab_a - u_0*(nai_120_001_1*p_1 + (nai_110_001_1 - nai_110_001_2)/ab_a - u_1*(p_2*tmp115 - u_2*(-tmp87 + tmp11*v_0)))); // final, weighs 71
  nai_010_001_1 = nai_300_000_0 + nai_020_001_0*p_2 - nai_020_001_1*u_2; // local+recycle, weighs 5
  nai_110_001_0 = tmp62 + nai_020_001_1*p_2 - nai_020_001_2*u_2; // local+recycle, weighs 5
  nai_110_001_1 = nai_010_001_1*v_0 - nai_110_001_0*u_0; // local+recycle, weighs 4
  nai_110_001_2 = tmp157*u_2; // auto+recycle, weighs 1
  nai_120_001_1 = tmp8 + nai_020_001_2*p_2 - u_2*(-nai_110_001_2 + p_2*tmp11); // local+recycle, weighs 8
  nai_200_001_2 = nai_110_001_0*v_0 - nai_120_001_1*u_0; // local+recycle, weighs 4
  out[35] = 10.8362471994438*tmp63*(nai_110_001_1*p_0 + (0.5*nai_010_001_1 - 0.5*nai_110_001_0)/ab_a - nai_200_001_2*u_0); // final, weighs 12
  nai_200_020_1 = nai_000_010_1*v_1 + (1.5*nai_000_020_1 - 1.5*nai_000_020_2)/ab_a - nai_000_010_2*u_1; // local+recycle, weighs 10
  nai_000_020_0 = v_1*(tmp70*v_1 + (1.5*nai_000_020_0 - 1.5*nai_000_020_1)/ab_a - nai_000_010_1*u_1) + (tmp74 + 1.5*nai_010_020_0 - 1.5*nai_010_020_1)/ab_a - nai_200_020_1*u_1; // local+recycle, weighs 21
  nai_000_010_4 = nai_200_020_1*v_1 + (nai_000_010_4 + 1.5*nai_010_020_1 - 1.5*nai_010_020_2)/ab_a - u_1*(nai_000_010_2*v_1 + (1.5*nai_000_020_2 - 1.5*nai_000_020_3)/ab_a - nai_000_010_3*u_1); // local+recycle, weighs 21
  out[36] = 4.84611707178961*tmp63*(nai_000_020_0*v_0 - nai_000_010_4*u_0); // final, weighs 6
  out[37] = 10.8362471994438*tmp63*(p_2*tmp161 - tmp108*u_2); // final, weighs 6
  out[38] = 10.8362471994438*tmp63*(nai_010_002_1 + nai_110_001_1*p_1 - nai_200_001_2*u_1); // final, weighs 7
  nai_000_020_1 = nai_200_010_1*v_1 - nai_300_100_0*u_1; // local+recycle, weighs 4
  nai_000_020_2 = tmp34 + v_1*(nai_200_010_0*v_1 - nai_200_010_1*u_1) - nai_000_020_1*u_1; // local+recycle, weighs 9
  nai_000_020_1 = tmp33 + nai_000_020_1*v_1 - u_1*(nai_300_100_0*v_1 - nai_000_001_4*u_1); // local+recycle, weighs 9
  out[39] = 4.84611707178961*tmp63*(nai_000_020_2*v_0 - nai_000_020_1*u_0); // final, weighs 6
  nai_010_002_1 = nai_300_200_0*v_2 - tmp167*u_2; // local+recycle, weighs 4
  nai_010_020_0 = tmp167*v_2 - tmp173*u_2; // local+recycle, weighs 4
  nai_010_020_1 = tmp173*v_2 - nai_300_200_1*u_2; // local+recycle, weighs 4
  nai_010_020_2 = nai_300_000_2*v_2 - nai_300_100_1*u_2; // local+recycle, weighs 4
  nai_110_001_1 = nai_300_100_1*v_2 - tmp164*u_2; // local+recycle, weighs 4
  nai_200_001_2 = nai_010_020_2*v_1 - nai_110_001_1*u_1; // local+recycle, weighs 4
  nai_200_020_1 = tmp164*v_2 - tmp170*u_2; // local+recycle, weighs 4
  nai_300_000_2 = nai_110_001_1*v_1 - nai_200_020_1*u_1; // local+recycle, weighs 4
  out[40] = 8.39372098776651*tmp63*(v_0*(nai_010_002_1*v_1 - nai_010_020_0*u_1) + (1.5*nai_200_001_2 - 1.5*nai_300_000_2)/ab_a - u_0*(nai_010_020_0*v_1 - nai_010_020_1*u_1)); // final, weighs 20
  nai_300_100_1 = d_2*u_1; // auto+recycle, weighs 1
  nai_300_200_0 = -nai_300_100_1 + d_1*v_1; // local+recycle, weighs 3
  nai_300_200_1 = nai_000_000_2*u_1; // auto+recycle, weighs 1
  tmp108 = -nai_300_200_1 + d_2*v_1; // local+recycle, weighs 3
  tmp115 = tmp108*u_0; // auto+recycle, weighs 1
  tmp157 = nai_000_000_3*u_1; // auto+recycle, weighs 1
  tmp161 = -tmp157 + nai_000_000_2*v_1; // local+recycle, weighs 3
  tmp164 = tmp161*u_0; // auto+recycle, weighs 1
  tmp167 = -tmp164 + p_0*tmp108; // local+recycle, weighs 3
  tmp173 = nai_200_001_2*v_0 + (-tmp115 - tmp167 + nai_300_200_0*p_0)/ab_a - nai_300_000_2*u_0; // local+recycle, weighs 12
  tmp33 = nai_000_000_4*u_1; // auto+recycle, weighs 1
  tmp34 = -tmp33 + nai_000_000_3*v_1; // local+recycle, weighs 3
  tmp74 = tmp34*u_0; // auto+recycle, weighs 1
  tmp167 = nai_300_000_2*v_0 + (tmp167 + tmp74 - p_0*tmp161)/ab_a - u_0*(nai_200_020_1*v_1 - u_1*(tmp170*v_2 - u_2*usq)); // local+recycle, weighs 19
  tmp123 = -tmp123 + d_2*p_0; // local+recycle, weighs 3
  nai_010_020_2 = nai_010_020_2*v_0 + (-tmp121 - tmp123 + d_1*p_0)/ab_a - nai_110_001_1*u_0; // local+recycle, weighs 12
  nai_110_001_1 = nai_110_001_1*v_0 + (tmp123 + tmp162 - nai_000_000_2*p_0)/ab_a - nai_200_020_1*u_0; // local+recycle, weighs 11
  out[41] = 18.7689307128126*tmp63*(p_1*tmp173 + (0.5*nai_010_020_2 - 0.5*nai_110_001_1)/ab_a - tmp167*u_1); // final, weighs 12
  out[42] = 18.7689307128126*tmp63*(p_2*tmp173 + (0.5*tmp140 - 0.5*tmp12)/ab_a - tmp167*u_2); // final, weighs 12
  nai_200_020_1 = -nai_300_200_1 + d_2*p_1; // local+recycle, weighs 3
  nai_300_100_1 = nai_300_000_1*v_1 + (-nai_200_020_1 - nai_300_100_1 + d_1*p_1)/ab_a - tmp106*u_1; // local+recycle, weighs 12
  nai_300_200_1 = -tmp157 + nai_000_000_2*p_1; // local+recycle, weighs 3
  nai_200_020_1 = tmp106*v_1 + (nai_200_020_1 - nai_300_200_1)/ab_a - tmp107*u_1; // local+recycle, weighs 9
  tmp12 = nai_300_100_1*v_0 - nai_200_020_1*u_0; // local+recycle, weighs 4
  nai_000_020_3 = nai_200_020_1*v_0 - u_0*(tmp107*v_1 + (nai_300_200_1 + tmp33 - nai_000_000_3*p_1)/ab_a - u_1*(nai_000_020_3*v_2 - nai_000_020_4*u_2)); // local+recycle, weighs 19
  out[43] = 18.7689307128126*tmp63*(p_0*tmp12 + (0.5*nai_300_100_1 - 0.5*nai_200_020_1)/ab_a - nai_000_020_3*u_0); // final, weighs 12
  nai_000_020_4 = -tmp164 + tmp108*v_0; // local+recycle, weighs 3
  nai_300_200_1 = -tmp74 + tmp161*v_0; // local+recycle, weighs 3
  tmp121 = nai_000_020_4*p_2 + (0.5*nai_110_000_1 - 0.5*nai_110_000_2)/ab_a - nai_300_200_1*u_2; // local+recycle, weighs 10
  d_0 = -d_0 + nai_000_000_5*v_2; // local+recycle, weighs 3
  nai_000_000_5 = nai_300_200_0*p_2 + (0.5*nai_200_010_2 - 0.5*tmp169)/ab_a - tmp108*u_2; // local+recycle, weighs 10
  nai_200_010_2 = p_2*tmp108 + (0.5*tmp169 - 0.5*tmp166)/ab_a - tmp161*u_2; // local+recycle, weighs 10
  tmp123 = tmp3 - d_2*u_2; // auto+recycle, weighs 3
  tmp140 = tmp123 + d_1*p_2; // local+recycle, weighs 2
  tmp157 = tmp2 - nai_000_000_2*u_2; // auto+recycle, weighs 3
  tmp162 = tmp157 + d_2*p_2; // local+recycle, weighs 2
  tmp163 = p_2*tmp161 + (0.5*tmp166 - 0.5*tmp163)/ab_a - tmp34*u_2; // local+recycle, weighs 10
  tmp1 = tmp1 - nai_000_000_3*u_2; // auto+recycle, weighs 3
  tmp164 = tmp1 + nai_000_000_2*p_2; // local+recycle, weighs 2
  out[44] = 32.5087415983314*tmp63*(p_0*(p_1*(p_2*(-tmp115 + nai_300_200_0*v_0) + (0.5*nai_110_000_0 - 0.5*nai_110_000_1)/ab_a - nai_000_020_4*u_2) + (0.5*nai_100_000_0 - 0.5*nai_100_000_1)/ab_a - tmp121*u_1) + (0.5*nai_000_000_5*p_1 + 0.5*tmp163*u_1 + 0.5*(0.5*tmp140 - 0.5*tmp162)/ab_a - 0.5*nai_200_010_2*p_1 - 0.5*nai_200_010_2*u_1 - 0.5*(0.5*tmp162 - 0.5*tmp164)/ab_a)/ab_a - u_0*(p_1*tmp121 + (0.5*nai_100_000_1 - 0.5*nai_100_000_2)/ab_a - u_1*(nai_300_200_1*p_2 + (0.5*nai_110_000_2 - 0.5*nai_110_000_3)/ab_a - u_2*(tmp34*v_0 - u_0*(nai_000_000_4*v_1 - d_0*u_1))))); // final, weighs 85
  nai_000_020_4 = nai_100_001_1*v_1 - nai_000_001_1*u_1; // local+recycle, weighs 4
  nai_100_000_0 = nai_000_001_1*v_1 - nai_000_001_2*u_1; // local+recycle, weighs 4
  nai_100_000_1 = nai_000_020_4*v_0 - nai_100_000_0*u_0; // local+recycle, weighs 4
  nai_000_001_3 = nai_100_000_0*v_0 - u_0*(nai_000_001_2*v_1 - u_1*(nai_000_001_3 + nai_000_002_3*v_2)); // local+recycle, weighs 10
  out[45] = 18.7689307128126*tmp63*(nai_100_000_1*p_0 + (0.5*nai_000_020_4 - 0.5*nai_100_000_0)/ab_a - nai_000_001_3*u_0); // final, weighs 12
  nai_100_000_2 = tmp70*v_2 - nai_000_010_1*u_2; // local+recycle, weighs 4
  nai_000_010_1 = nai_000_010_1*v_2 - nai_000_010_2*u_2; // local+recycle, weighs 4
  nai_000_010_2 = nai_000_010_2*v_2 - nai_000_010_3*u_2; // local+recycle, weighs 4
  out[46] = 8.39372098776651*tmp63*(v_0*(nai_100_000_2*v_1 + (1.5*nai_300_000_1 - 1.5*tmp106)/ab_a - nai_000_010_1*u_1) - u_0*(nai_000_010_1*v_1 + (1.5*tmp106 - 1.5*tmp107)/ab_a - nai_000_010_2*u_1)); // final, weighs 26
  out[47] = 18.7689307128126*tmp63*(p_2*tmp12 + (0.5*nai_110_020_0 - 0.5*nai_110_020_1)/ab_a - nai_000_020_3*u_2); // final, weighs 12
  out[48] = 18.7689307128126*tmp63*(nai_100_000_1*p_1 + (0.5*nai_200_000_1 - 0.5*nai_200_000_2)/ab_a - nai_000_001_3*u_1); // final, weighs 12
  nai_000_001_3 = nai_200_010_0*v_2 + (1.5*nai_000_002_0 - 1.5*nai_000_002_1)/ab_a - nai_200_010_1*u_2; // local+recycle, weighs 10
  nai_000_002_0 = nai_200_010_1*v_2 + (1.5*nai_000_002_1 - 1.5*nai_000_002_2)/ab_a - nai_300_100_0*u_2; // local+recycle, weighs 10
  nai_000_001_4 = nai_300_100_0*v_2 + (1.5*nai_000_002_2 - 1.5*nai_000_002_3)/ab_a - nai_000_001_4*u_2; // local+recycle, weighs 10
  out[49] = 8.39372098776651*tmp63*(v_0*(nai_000_001_3*v_1 - nai_000_002_0*u_1) - u_0*(nai_000_002_0*v_1 - nai_000_001_4*u_1)); // final, weighs 14
  nai_000_002_1 = nai_100_000_4 + nai_010_002_1*v_2 - nai_010_020_0*u_2; // local+recycle, weighs 5
  nai_000_002_2 = nai_200_000_0 + nai_010_020_0*v_2 - nai_010_020_1*u_2; // local+recycle, weighs 5
  nai_000_002_3 = tmp123 + d_1*v_2; // local+recycle, weighs 2
  nai_000_010_3 = tmp157 + d_2*v_2; // local+recycle, weighs 2
  nai_000_020_3 = nai_000_010_3*u_0; // auto+recycle, weighs 1
  nai_010_002_1 = -nai_000_020_3 + nai_000_002_3*p_0; // local+recycle, weighs 3
  nai_010_020_0 = tmp1 + nai_000_000_2*v_2; // local+recycle, weighs 2
  nai_010_020_1 = nai_010_020_0*u_0; // auto+recycle, weighs 1
  nai_100_000_1 = -nai_010_020_1 + nai_000_010_3*p_0; // local+recycle, weighs 3
  nai_100_000_4 = 0.5*nai_000_002_3 - 0.5*nai_000_010_3; // auto+recycle, weighs 3
  nai_110_000_0 = nai_100_000_4/ab_a; // auto+recycle, weighs 2
  nai_110_000_1 = nai_110_000_0 + nai_010_002_1*p_0 - nai_100_000_1*u_0; // local+recycle, weighs 5
  nai_110_000_2 = tmp0 - nai_000_000_4*u_2; // auto+recycle, weighs 3
  nai_110_000_3 = nai_110_000_2 + nai_000_000_3*v_2; // local+recycle, weighs 2
  nai_110_020_0 = nai_110_000_3*u_0; // auto+recycle, weighs 1
  nai_110_020_1 = -nai_110_020_0 + nai_010_020_0*p_0; // local+recycle, weighs 3
  nai_200_000_0 = 0.5*nai_000_010_3 - 0.5*nai_010_020_0; // auto+recycle, weighs 3
  nai_200_000_1 = nai_200_000_0/ab_a; // auto+recycle, weighs 2
  nai_200_000_2 = nai_200_000_1 + nai_100_000_1*p_0 - nai_110_020_1*u_0; // local+recycle, weighs 5
  nai_200_010_0 = (1.5*nai_110_000_1 - 1.5*nai_200_000_2)/ab_a; // auto+recycle, weighs 5
  out[50] = 4.84611707178961*tmp63*(nai_200_010_0 + nai_000_002_1*v_0 - nai_000_002_2*u_0); // final, weighs 7
  nai_010_002_1 = nai_110_000_1*v_0 + (nai_010_002_1 - nai_100_000_1)/ab_a - nai_200_000_2*u_0; // local+recycle, weighs 9
  d_0 = tmp4 + nai_000_000_4*v_2 - d_0*u_2; // local+recycle, weighs 5
  nai_200_010_1 = d_0*u_0; // auto+recycle, weighs 1
  nai_300_000_1 = 0.5*nai_010_020_0 - 0.5*nai_110_000_3; // auto+recycle, weighs 3
  nai_300_100_0 = nai_300_000_1/ab_a; // auto+recycle, weighs 2
  nai_300_200_1 = nai_300_100_0 + nai_110_020_1*p_0 - u_0*(-nai_200_010_1 + nai_110_000_3*p_0); // local+recycle, weighs 8
  nai_100_000_1 = nai_200_000_2*v_0 + (nai_100_000_1 - nai_110_020_1)/ab_a - nai_300_200_1*u_0; // local+recycle, weighs 9
  out[51] = 10.8362471994438*tmp63*(nai_010_002_1*p_1 - nai_100_000_1*u_1); // final, weighs 6
  out[52] = 10.8362471994438*tmp63*(nai_010_002_1*p_2 + (nai_010_020_2 - nai_110_001_1)/ab_a - nai_100_000_1*u_2); // final, weighs 11
  nai_010_002_1 = nai_000_010_3*u_1; // auto+recycle, weighs 1
  nai_010_020_2 = -nai_010_002_1 + nai_000_002_3*p_1; // local+recycle, weighs 3
  nai_100_000_1 = nai_010_020_0*u_1; // auto+recycle, weighs 1
  nai_110_001_1 = -nai_100_000_1 + nai_000_010_3*p_1; // local+recycle, weighs 3
  nai_110_000_0 = nai_110_000_0 + nai_010_020_2*p_1 - nai_110_001_1*u_1; // local+recycle, weighs 5
  nai_110_020_1 = nai_110_000_3*u_1; // auto+recycle, weighs 1
  tmp0 = -nai_110_020_1 + nai_010_020_0*p_1; // local+recycle, weighs 3
  nai_200_000_1 = nai_200_000_1 + nai_110_001_1*p_1 - tmp0*u_1; // local+recycle, weighs 5
  tmp1 = nai_110_000_0*v_0 - nai_200_000_1*u_0; // local+recycle, weighs 4
  tmp106 = d_0*u_1; // auto+recycle, weighs 1
  nai_300_100_0 = nai_300_100_0 + p_1*tmp0 - u_1*(-tmp106 + nai_110_000_3*p_1); // local+recycle, weighs 8
  tmp107 = nai_200_000_1*v_0 - nai_300_100_0*u_0; // local+recycle, weighs 4
  out[53] = 10.8362471994438*tmp63*(p_0*tmp1 + (0.5*nai_110_000_0 - 0.5*nai_200_000_1)/ab_a - tmp107*u_0); // final, weighs 12
  nai_010_020_1 = -nai_010_020_1 + nai_000_010_3*v_0; // local+recycle, weighs 3
  nai_110_020_0 = -nai_110_020_0 + nai_010_020_0*v_0; // local+recycle, weighs 3
  tmp115 = nai_010_020_1*p_2 + (tmp126 - tmp165)/ab_a - nai_110_020_0*u_2; // local+recycle, weighs 9
  d_1 = (d_1 - d_2)/ab_a - nai_000_010_3*u_2; // auto+recycle, weighs 7
  tmp12 = d_1 + nai_000_002_3*p_2; // local+recycle, weighs 2
  d_2 = (d_2 - nai_000_000_2)/ab_a - nai_010_020_0*u_2; // auto+recycle, weighs 7
  tmp121 = d_2 + nai_000_010_3*p_2; // local+recycle, weighs 2
  nai_000_000_2 = (nai_000_000_2 - nai_000_000_3)/ab_a - nai_110_000_3*u_2; // auto+recycle, weighs 7
  tmp123 = nai_000_000_2 + nai_010_020_0*p_2; // local+recycle, weighs 2
  out[54] = 18.7689307128126*tmp63*(p_0*(p_1*(p_2*(-nai_000_020_3 + nai_000_002_3*v_0) + (tmp122 - tmp126)/ab_a - nai_010_020_1*u_2) - tmp115*u_1) + (0.5*p_1*tmp12 + 0.5*tmp123*u_1 - 0.5*p_1*tmp121 - 0.5*tmp121*u_1)/ab_a - u_0*(p_1*tmp115 - u_1*(nai_110_020_0*p_2 + (tmp165 - tmp159)/ab_a - u_2*(-nai_200_010_1 + nai_110_000_3*v_0)))); // final, weighs 52
  nai_000_020_3 = p_2*tmp12 + (nai_100_000_4 + tmp140 - tmp162)/ab_a - tmp121*u_2; // local+recycle, weighs 10
  nai_010_020_1 = p_2*tmp121 + (nai_200_000_0 + tmp162 - tmp164)/ab_a - tmp123*u_2; // local+recycle, weighs 10
  nai_100_000_4 = nai_000_020_3*v_0 - nai_010_020_1*u_0; // local+recycle, weighs 4
  d_0 = (nai_000_000_3 - nai_000_000_4)/ab_a - d_0*u_2; // auto+recycle, weighs 7
  nai_000_000_3 = p_2*tmp123 + (nai_300_000_1 + tmp164 - nai_110_000_2 - nai_000_000_3*p_2)/ab_a - u_2*(d_0 + nai_110_000_3*p_2); // local+recycle, weighs 15
  nai_000_000_4 = nai_010_020_1*v_0 - nai_000_000_3*u_0; // local+recycle, weighs 4
  nai_110_000_2 = (0.5*nai_000_020_3 - 0.5*nai_010_020_1)/ab_a; // auto+recycle, weighs 5
  out[55] = 10.8362471994438*tmp63*(nai_110_000_2 + nai_100_000_4*p_0 - nai_000_000_4*u_0); // final, weighs 7
  nai_100_000_2 = tmp21 + nai_100_000_2*v_2 - nai_000_010_1*u_2; // local+recycle, weighs 5
  nai_000_010_1 = tmp20 + nai_000_010_1*v_2 - nai_000_010_2*u_2; // local+recycle, weighs 5
  out[56] = 4.84611707178961*tmp63*(nai_100_000_2*v_0 - nai_000_010_1*u_0); // final, weighs 6
  out[57] = 10.8362471994438*tmp63*(tmp110 + p_2*tmp1 - tmp107*u_2); // final, weighs 7
  out[58] = 10.8362471994438*tmp63*(nai_100_000_4*p_1 - nai_000_000_4*u_1); // final, weighs 6
  nai_000_000_4 = nai_000_001_3*v_2 + (tmp124 + 1.5*nai_100_001_1 - 1.5*nai_000_001_1)/ab_a - nai_000_002_0*u_2; // local+recycle, weighs 11
  nai_000_001_1 = nai_000_002_0*v_2 + (tmp88 + 1.5*nai_000_001_1 - 1.5*nai_000_001_2)/ab_a - nai_000_001_4*u_2; // local+recycle, weighs 11
  out[59] = 4.84611707178961*tmp63*(nai_000_000_4*v_0 - nai_000_001_1*u_0); // final, weighs 6
  nai_000_001_2 = nai_200_002_1 + nai_200_000_4*v_1; // local+recycle, weighs 2
  nai_000_001_3 = tmp138 + nai_200_001_1*v_1; // local+recycle, weighs 2
  nai_000_001_4 = nai_000_001_2*p_0 - nai_000_001_3*u_0; // local+recycle, weighs 4
  nai_000_002_0 = nai_210_001_1 + nai_200_020_0*v_1; // local+recycle, weighs 2
  nai_000_010_2 = nai_000_001_3*p_0 - nai_000_002_0*u_0; // local+recycle, weighs 4
  nai_100_000_4 = 0.5*nai_000_001_2 - 0.5*nai_000_001_3; // auto+recycle, weighs 3
  nai_100_001_1 = nai_100_000_4/ab_a; // auto+recycle, weighs 2
  nai_110_020_0 = nai_100_001_1 + nai_000_001_4*p_0 - nai_000_010_2*u_0; // local+recycle, weighs 5
  nai_200_000_0 = tmp160 + tmp11*v_1; // local+recycle, weighs 2
  nai_200_002_1 = 0.5*nai_000_001_3 - 0.5*nai_000_002_0; // auto+recycle, weighs 3
  nai_200_010_1 = nai_200_002_1/ab_a; // auto+recycle, weighs 2
  nai_210_001_1 = nai_200_010_1 + nai_000_010_2*p_0 - u_0*(nai_000_002_0*p_0 - nai_200_000_0*u_0); // local+recycle, weighs 9
  out[60] = 2.16724943988876*tmp63*(nai_110_020_0*p_0 + (nai_000_001_4 - nai_000_010_2)/ab_a - nai_210_001_1*u_0); // final, weighs 11
  out[61] = 4.84611707178961*tmp63*(tmp77 + nai_110_020_0*p_1 - nai_210_001_1*u_1); // final, weighs 7
  out[62] = 4.84611707178961*tmp63*(nai_110_020_0*p_2 - nai_210_001_1*u_2); // final, weighs 6
  nai_000_001_4 = nai_000_001_2*p_1 + (1.5*nai_200_000_4 - 1.5*nai_200_001_1)/ab_a - nai_000_001_3*u_1; // local+recycle, weighs 10
  nai_000_010_2 = nai_000_001_3*p_1 + (1.5*nai_200_001_1 - 1.5*nai_200_020_0)/ab_a - nai_000_002_0*u_1; // local+recycle, weighs 10
  nai_100_000_4 = nai_000_001_4*p_1 + (nai_100_000_4 + 1.5*nai_200_030_0 - 1.5*tmp139)/ab_a - nai_000_010_2*u_1; // local+recycle, weighs 11
  nai_020_010_2 = nai_000_010_2*p_1 + (nai_200_002_1 + 1.5*tmp139 - 1.5*nai_020_010_2)/ab_a - u_1*(nai_000_002_0*p_1 + (1.5*nai_200_020_0 - 1.5*tmp11)/ab_a - nai_200_000_0*u_1); // local+recycle, weighs 21
  out[63] = 4.84611707178961*tmp63*(nai_100_000_4*p_0 - nai_020_010_2*u_0); // final, weighs 6
  nai_000_001_2 = nai_000_001_2*p_2 - nai_000_001_3*u_2; // local+recycle, weighs 4
  nai_000_001_3 = nai_000_001_3*p_2 - nai_000_002_0*u_2; // local+recycle, weighs 4
  nai_000_002_0 = nai_000_002_0*p_2 - nai_200_000_0*u_2; // local+recycle, weighs 4
  out[64] = 8.39372098776651*tmp63*(p_0*(nai_000_001_2*p_1 + (1.5*nai_020_001_0 - 1.5*nai_020_001_1)/ab_a - nai_000_001_3*u_1) - u_0*(nai_000_001_3*p_1 + (1.5*nai_020_001_1 - 1.5*nai_020_001_2)/ab_a - nai_000_002_0*u_1)); // final, weighs 26
  nai_100_001_1 = nai_100_001_1 + nai_000_001_2*p_2 - nai_000_001_3*u_2; // local+recycle, weighs 5
  nai_000_002_0 = nai_200_010_1 + nai_000_001_3*p_2 - nai_000_002_0*u_2; // local+recycle, weighs 5
  out[65] = 4.84611707178961*tmp63*(nai_100_001_1*p_0 - nai_000_002_0*u_0); // final, weighs 6
  out[66] = 2.16724943988876*tmp63*(nai_100_000_4*p_1 + (nai_000_001_4 - nai_000_010_2 + 1.5*nai_200_000_3 - 1.5*nai_200_001_0)/ab_a - nai_020_010_2*u_1); // final, weighs 15
  out[67] = 4.84611707178961*tmp63*(nai_100_000_4*p_2 - nai_020_010_2*u_2); // final, weighs 6
  nai_000_001_4 = (1.5*nai_010_001_1 - 1.5*nai_110_001_0)/ab_a; // auto+recycle, weighs 5
  out[68] = 4.84611707178961*tmp63*(nai_000_001_4 + nai_100_001_1*p_1 - nai_000_002_0*u_1); // final, weighs 7
  out[69] = 2.16724943988876*tmp63*(nai_100_001_1*p_2 + (nai_000_001_2 - nai_000_001_3)/ab_a - nai_000_002_0*u_2); // final, weighs 11
  out[70] = 4.84611707178961*tmp63*(nai_100_030_1*v_2 - nai_100_000_3*u_2); // final, weighs 6
  nai_000_001_2 = nai_300_000_3*v_2 - tmp69*u_2; // local+recycle, weighs 4
  nai_000_001_3 = tmp69*v_2 - nai_020_200_2*u_2; // local+recycle, weighs 4
  nai_000_002_0 = (nai_200_001_2 - nai_300_000_2)/ab_a; // auto+recycle, weighs 4
  out[71] = 10.8362471994438*tmp63*(nai_000_002_0 + nai_000_001_2*p_1 - nai_000_001_3*u_1); // final, weighs 7
  out[72] = 10.8362471994438*tmp63*(nai_000_001_2*p_2 + (0.5*nai_300_000_3 - 0.5*tmp69)/ab_a - nai_000_001_3*u_2); // final, weighs 12
  nai_000_001_2 = nai_200_000_3*v_2 - nai_200_001_0*u_2; // local+recycle, weighs 4
  nai_000_001_3 = nai_200_001_0*v_2 - nai_200_003_0*u_2; // local+recycle, weighs 4
  out[73] = 10.8362471994438*tmp63*(nai_000_001_2*p_0 - nai_000_001_3*u_0); // final, weighs 6
  nai_000_010_2 = -tmp131 + nai_200_001_1*v_2; // local+recycle, weighs 3
  nai_020_010_2 = -tmp129 + nai_200_020_0*v_2; // local+recycle, weighs 3
  nai_020_200_2 = tmp62 + nai_000_010_2*p_2 - nai_020_010_2*u_2; // local+recycle, weighs 5
  out[74] = 18.7689307128126*tmp63*(p_0*(p_1*(nai_300_000_0 + p_2*(-tmp133 + nai_200_000_4*v_2) - nai_000_010_2*u_2) + (nai_000_000_5 - nai_200_010_2)/ab_a - nai_020_200_2*u_1) - u_0*(nai_020_200_2*p_1 + (nai_200_010_2 - tmp163)/ab_a - u_1*(tmp8 + nai_020_010_2*p_2 - u_2*(-nai_110_001_2 + tmp11*v_2)))); // final, weighs 40
  nai_000_000_5 = nai_010_001_1*v_2 + (nai_020_001_0 - nai_020_001_1)/ab_a - nai_110_001_0*u_2; // local+recycle, weighs 9
  nai_000_010_2 = nai_110_001_0*v_2 + (nai_020_001_1 - nai_020_001_2)/ab_a - nai_120_001_1*u_2; // local+recycle, weighs 9
  out[75] = 10.8362471994438*tmp63*(nai_000_000_5*p_0 - nai_000_010_2*u_0); // final, weighs 6
  out[76] = 4.84611707178961*tmp63*(nai_000_020_0*v_2 - nai_000_010_4*u_2); // final, weighs 6
  out[77] = 10.8362471994438*tmp63*(tmp80 + nai_000_001_2*p_2 - nai_000_001_3*u_2); // final, weighs 7
  out[78] = 10.8362471994438*tmp63*(nai_000_000_5*p_1 + (nai_000_020_4 - nai_100_000_0)/ab_a - nai_000_010_2*u_1); // final, weighs 11
  out[79] = 4.84611707178961*tmp63*(nai_000_001_4 + nai_000_020_2*v_2 - nai_000_020_1*u_2); // final, weighs 7
  out[80] = 4.84611707178961*tmp63*(nai_000_002_1*v_1 - nai_000_002_2*u_1); // final, weighs 6
  nai_000_000_5 = nai_110_000_1*v_1 - nai_200_000_2*u_1; // local+recycle, weighs 4
  nai_000_001_2 = nai_200_000_2*v_1 - nai_300_200_1*u_1; // local+recycle, weighs 4
  out[81] = 10.8362471994438*tmp63*(nai_000_000_5*p_1 + (0.5*nai_110_000_1 - 0.5*nai_200_000_2)/ab_a - nai_000_001_2*u_1); // final, weighs 12
  out[82] = 10.8362471994438*tmp63*(nai_000_002_0 + nai_000_000_5*p_2 - nai_000_001_2*u_2); // final, weighs 7
  nai_000_000_5 = nai_110_000_0*v_1 + (nai_010_020_2 - nai_110_001_1)/ab_a - nai_200_000_1*u_1; // local+recycle, weighs 9
  nai_000_001_2 = nai_200_000_1*v_1 + (nai_110_001_1 - tmp0)/ab_a - nai_300_100_0*u_1; // local+recycle, weighs 9
  out[83] = 10.8362471994438*tmp63*(nai_000_000_5*p_0 - nai_000_001_2*u_0); // final, weighs 6
  nai_000_001_3 = -nai_100_000_1 + nai_000_010_3*v_1; // local+recycle, weighs 3
  nai_000_001_4 = -nai_110_020_1 + nai_010_020_0*v_1; // local+recycle, weighs 3
  nai_000_002_0 = nai_000_001_3*p_2 + (tmp108 - tmp161)/ab_a - nai_000_001_4*u_2; // local+recycle, weighs 9
  out[84] = 18.7689307128126*tmp63*(p_0*(p_1*(p_2*(-nai_010_002_1 + nai_000_002_3*v_1) + (nai_300_200_0 - tmp108)/ab_a - nai_000_001_3*u_2) + (0.5*tmp12 - 0.5*tmp121)/ab_a - nai_000_002_0*u_1) - u_0*(nai_000_002_0*p_1 + (0.5*tmp121 - 0.5*tmp123)/ab_a - u_1*(nai_000_001_4*p_2 + (tmp161 - tmp34)/ab_a - u_2*(-tmp106 + nai_110_000_3*v_1)))); // final, weighs 50
  nai_000_001_3 = nai_000_020_3*v_1 - nai_010_020_1*u_1; // local+recycle, weighs 4
  nai_000_000_3 = nai_010_020_1*v_1 - nai_000_000_3*u_1; // local+recycle, weighs 4
  out[85] = 10.8362471994438*tmp63*(nai_000_001_3*p_0 - nai_000_000_3*u_0); // final, weighs 6
  nai_000_001_4 = (1.5*nai_110_000_0 - 1.5*nai_200_000_1)/ab_a; // auto+recycle, weighs 5
  out[86] = 4.84611707178961*tmp63*(nai_000_001_4 + nai_100_000_2*v_1 - nai_000_010_1*u_1); // final, weighs 7
  out[87] = 10.8362471994438*tmp63*(nai_000_000_5*p_2 + (nai_300_100_1 - nai_200_020_1)/ab_a - nai_000_001_2*u_2); // final, weighs 11
  out[88] = 10.8362471994438*tmp63*(nai_110_000_2 + nai_000_001_3*p_1 - nai_000_000_3*u_1); // final, weighs 7
  out[89] = 4.84611707178961*tmp63*(nai_000_000_4*v_1 - nai_000_001_1*u_1); // final, weighs 6
  d_1 = d_1 + nai_000_002_3*v_2; // local+recycle, weighs 2
  d_2 = d_2 + nai_000_010_3*v_2; // local+recycle, weighs 2
  nai_000_000_3 = d_1*p_0 - d_2*u_0; // local+recycle, weighs 4
  nai_000_000_2 = nai_000_000_2 + nai_010_020_0*v_2; // local+recycle, weighs 2
  nai_000_000_4 = d_2*p_0 - nai_000_000_2*u_0; // local+recycle, weighs 4
  nai_000_000_5 = 0.5*d_1 - 0.5*d_2; // auto+recycle, weighs 3
  nai_000_001_1 = nai_000_000_5/ab_a; // auto+recycle, weighs 2
  nai_000_001_2 = nai_000_001_1 + nai_000_000_3*p_0 - nai_000_000_4*u_0; // local+recycle, weighs 5
  d_0 = d_0 + nai_110_000_3*v_2; // local+recycle, weighs 2
  nai_000_001_3 = 0.5*d_2 - 0.5*nai_000_000_2; // auto+recycle, weighs 3
  nai_000_002_0 = nai_000_001_3/ab_a; // auto+recycle, weighs 2
  nai_000_002_1 = nai_000_002_0 + nai_000_000_4*p_0 - u_0*(nai_000_000_2*p_0 - d_0*u_0); // local+recycle, weighs 9
  out[90] = 2.16724943988876*tmp63*(nai_000_001_2*p_0 + (nai_000_000_3 - nai_000_000_4)/ab_a - nai_000_002_1*u_0); // final, weighs 11
  out[91] = 4.84611707178961*tmp63*(nai_000_001_2*p_1 - nai_000_002_1*u_1); // final, weighs 6
  out[92] = 4.84611707178961*tmp63*(nai_200_010_0 + nai_000_001_2*p_2 - nai_000_002_1*u_2); // final, weighs 7
  nai_000_000_3 = d_1*p_1 - d_2*u_1; // local+recycle, weighs 4
  nai_000_000_4 = d_2*p_1 - nai_000_000_2*u_1; // local+recycle, weighs 4
  nai_000_001_1 = nai_000_001_1 + nai_000_000_3*p_1 - nai_000_000_4*u_1; // local+recycle, weighs 5
  nai_000_001_2 = nai_000_002_0 + nai_000_000_4*p_1 - u_1*(nai_000_000_2*p_1 - d_0*u_1); // local+recycle, weighs 9
  out[93] = 4.84611707178961*tmp63*(nai_000_001_1*p_0 - nai_000_001_2*u_0); // final, weighs 6
  d_1 = d_1*p_2 + (1.5*nai_000_002_3 - 1.5*nai_000_010_3)/ab_a - d_2*u_2; // local+recycle, weighs 10
  d_2 = d_2*p_2 + (1.5*nai_000_010_3 - 1.5*nai_010_020_0)/ab_a - nai_000_000_2*u_2; // local+recycle, weighs 10
  d_0 = nai_000_000_2*p_2 + (1.5*nai_010_020_0 - 1.5*nai_110_000_3)/ab_a - d_0*u_2; // local+recycle, weighs 10
  out[94] = 8.39372098776651*tmp63*(p_0*(d_1*p_1 - d_2*u_1) - u_0*(d_2*p_1 - d_0*u_1)); // final, weighs 14
  nai_000_000_2 = d_1*p_2 + (nai_000_000_5 + 1.5*tmp12 - 1.5*tmp121)/ab_a - d_2*u_2; // local+recycle, weighs 11
  d_0 = d_2*p_2 + (nai_000_001_3 + 1.5*tmp121 - 1.5*tmp123)/ab_a - d_0*u_2; // local+recycle, weighs 11
  out[95] = 4.84611707178961*tmp63*(nai_000_000_2*p_0 - d_0*u_0); // final, weighs 6
  out[96] = 2.16724943988876*tmp63*(nai_000_001_1*p_1 + (nai_000_000_3 - nai_000_000_4)/ab_a - nai_000_001_2*u_1); // final, weighs 11
  out[97] = 4.84611707178961*tmp63*(nai_000_001_4 + nai_000_001_1*p_2 - nai_000_001_2*u_2); // final, weighs 7
  out[98] = 4.84611707178961*tmp63*(nai_000_000_2*p_1 - d_0*u_1); // final, weighs 6
  out[99] = 2.16724943988876*tmp63*(nai_000_000_2*p_2 + (d_1 - d_2 + 1.5*nai_000_020_3 - 1.5*nai_010_020_1)/ab_a - d_0*u_2); // final, weighs 15
  // total weight = 3598
}

typedef void (*fntype)(double*, double, double*, double, double*, double*);
const fntype fns[49] = {gint2_nai_pF_pF, gint2_nai_pD_pF, gint2_nai_SP_pF, gint2_nai_S_pF, gint2_nai_P_pF, gint2_nai_cD_pF, gint2_nai_cF_pF, gint2_nai_pF_pD, gint2_nai_pD_pD, gint2_nai_SP_pD, gint2_nai_S_pD, gint2_nai_P_pD, gint2_nai_cD_pD, gint2_nai_cF_pD, gint2_nai_pF_SP, gint2_nai_pD_SP, gint2_nai_SP_SP, gint2_nai_S_SP, gint2_nai_P_SP, gint2_nai_cD_SP, gint2_nai_cF_SP, gint2_nai_pF_S, gint2_nai_pD_S, gint2_nai_SP_S, gint2_nai_S_S, gint2_nai_P_S, gint2_nai_cD_S, gint2_nai_cF_S, gint2_nai_pF_P, gint2_nai_pD_P, gint2_nai_SP_P, gint2_nai_S_P, gint2_nai_P_P, gint2_nai_cD_P, gint2_nai_cF_P, gint2_nai_pF_cD, gint2_nai_pD_cD, gint2_nai_SP_cD, gint2_nai_S_cD, gint2_nai_P_cD, gint2_nai_cD_cD, gint2_nai_cF_cD, gint2_nai_pF_cF, gint2_nai_pD_cF, gint2_nai_SP_cF, gint2_nai_S_cF, gint2_nai_P_cF, gint2_nai_cD_cF, gint2_nai_cF_cF};

void gint2_nai_dispatch(int a_s, double* a, double a_a, int b_s, double* b, double b_a, double* c, double* out)
{
  fns[3+a_s+7*(3+b_s)](a, a_a, b, b_a, c, out);
}

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
      //printf("shell1=%d  type=%d  dof=%d  offset=%d\n", shell1, shell_type1, shell_dof1, shell_offset1);
      // prep inner loop.
      shell_ccoeffs2 = ccoeffs;
      shell_exponents2 = exponents;
      shell_offset2 = 0;
      for (shell2=0; shell2<num_shells; shell2++) {
        center2 = centers + (3*shell_map[shell2]);
        shell_type2 = shell_types[shell2];
        shell_dof2 = get_shell_dof(shell_type2);
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
