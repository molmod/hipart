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
#define MAX_SHELL 1
#define NUM_SHELL_TYPES 3
#define MAX_SHELL_DOF 4

static void gint2_nai_SP_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  double fix0, fix1, fix2, fix3, fix4, fix5, tmp12, tmp13, tmp14, tmp20, tmp25, tmp26, tmp27, tmp28, tmp29, tmp3, tmp30, tmp31, tmp32, tmp7;
  tmp31 = -b[0]; // auto
  fix0 = a[0] + tmp31; // oblige
  tmp30 = -b[1]; // auto
  fix1 = a[1] + tmp30; // oblige
  tmp29 = -b[2]; // auto
  fix2 = a[2] + tmp29; // oblige
  tmp25 = a_a + b_a; // auto
  tmp7 = 1.0/tmp25; // auto
  tmp14 = tmp7*(a_a*a[0] + b_a*b[0]); // auto
  fix3 = tmp14 - c[0]; // oblige
  tmp13 = tmp7*(a_a*a[1] + b_a*b[1]); // auto
  fix4 = tmp13 - c[1]; // oblige
  tmp12 = tmp7*(a_a*a[2] + b_a*b[2]); // auto
  fix5 = tmp12 - c[2]; // oblige
  fix1 = M_PI*tmp7*exp(-a_a*b_a*tmp7*(fix0*fix0 + fix1*fix1 + fix2*fix2)); // auto+recycle
  tmp25 = tmp25*(fix3*fix3 + fix4*fix4 + fix5*fix5); // auto+recycle
  fix0 = fix1*gaux(tmp25, 0); // auto+recycle
  tmp7 = pow(b_a,0.75); // auto+recycle
  fix2 = pow(a_a,0.75); // auto+recycle
  out[0] = 1.01589817494786*fix0*fix2*tmp7; // final
  tmp31 = tmp14 + tmp31; // oblige+recycle
  fix0 = 2*fix0; // auto+recycle
  tmp20 = fix1*gaux(tmp25, 1); // auto
  tmp3 = -2*tmp20; // auto
  tmp32 = pow(b_a,1.25); // auto
  fix2 = fix2*tmp32; // auto+recycle
  tmp28 = fix3*tmp3; // auto
  out[1] = 1.01589817494786*fix2*(tmp28 + fix0*tmp31); // final
  tmp30 = tmp13 + tmp30; // oblige+recycle
  tmp27 = fix4*tmp3; // auto
  out[2] = 1.01589817494786*fix2*(tmp27 + fix0*tmp30); // final
  tmp29 = tmp12 + tmp29; // oblige+recycle
  tmp26 = fix5*tmp3; // auto
  out[3] = 1.01589817494786*fix2*(tmp26 + fix0*tmp29); // final
  tmp14 = tmp14 - a[0]; // oblige+recycle
  tmp28 = tmp28 + fix0*tmp14; // auto+recycle
  fix2 = pow(a_a,1.25); // auto+recycle
  tmp7 = fix2*tmp7; // auto+recycle
  out[4] = 1.01589817494786*tmp28*tmp7; // final
  tmp20 = 2*tmp20; // auto+recycle
  tmp25 = -2*fix1*gaux(tmp25, 2); // auto+recycle
  tmp32 = fix2*tmp32; // auto+recycle
  fix1 = (fix0 + tmp3)/(2*a_a + 2*b_a); // auto+recycle
  tmp14 = fix3*tmp25 + tmp14*tmp20; // auto+recycle
  tmp14 = -tmp14; // auto+recycle
  out[5] = 2.03179634989571*tmp32*(fix1 + fix3*tmp14 + tmp28*tmp31); // final
  out[6] = 2.03179634989571*tmp32*(fix4*tmp14 + tmp28*tmp30); // final
  out[7] = 2.03179634989571*tmp32*(fix5*tmp14 + tmp28*tmp29); // final
  tmp14 = tmp13 - a[1]; // oblige+recycle
  tmp28 = tmp27 + fix0*tmp14; // auto+recycle
  out[8] = 1.01589817494786*tmp28*tmp7; // final
  tmp14 = fix4*tmp25 + tmp14*tmp20; // auto+recycle
  tmp14 = -tmp14; // auto+recycle
  out[9] = 2.03179634989571*tmp32*(fix3*tmp14 + tmp28*tmp31); // final
  out[10] = 2.03179634989571*tmp32*(fix1 + fix4*tmp14 + tmp28*tmp30); // final
  out[11] = 2.03179634989571*tmp32*(fix5*tmp14 + tmp28*tmp29); // final
  tmp14 = tmp12 - a[2]; // oblige+recycle
  tmp28 = tmp26 + fix0*tmp14; // auto+recycle
  out[12] = 1.01589817494786*tmp28*tmp7; // final
  tmp25 = fix5*tmp25 + tmp14*tmp20; // auto+recycle
  tmp25 = -tmp25; // auto+recycle
  out[13] = 2.03179634989571*tmp32*(fix3*tmp25 + tmp28*tmp31); // final
  out[14] = 2.03179634989571*tmp32*(fix4*tmp25 + tmp28*tmp30); // final
  out[15] = 2.03179634989571*tmp32*(fix1 + fix5*tmp25 + tmp28*tmp29); // final
}

static void gint2_nai_S_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  double fix0, fix1, fix2, fix3, fix4, fix5, tmp10, tmp11, tmp12, tmp13, tmp2, tmp6, tmp7, tmp8;
  tmp13 = -b[0]; // auto
  fix0 = a[0] + tmp13; // oblige
  tmp12 = -b[1]; // auto
  fix1 = a[1] + tmp12; // oblige
  tmp11 = -b[2]; // auto
  fix2 = a[2] + tmp11; // oblige
  tmp10 = a_a + b_a; // auto
  tmp2 = 1.0/tmp10; // auto
  tmp8 = tmp2*(a_a*a[0] + b_a*b[0]); // auto
  fix3 = tmp8 - c[0]; // oblige
  tmp7 = tmp2*(a_a*a[1] + b_a*b[1]); // auto
  fix4 = tmp7 - c[1]; // oblige
  tmp6 = tmp2*(a_a*a[2] + b_a*b[2]); // auto
  fix5 = tmp6 - c[2]; // oblige
  tmp2 = M_PI*tmp2*exp(-a_a*b_a*tmp2*(fix0*fix0 + fix1*fix1 + fix2*fix2)); // auto+recycle
  fix1 = tmp10*(fix3*fix3 + fix4*fix4 + fix5*fix5); // auto+recycle
  tmp10 = tmp2*gaux(fix1, 0); // auto+recycle
  fix0 = pow(a_a,0.75); // auto+recycle
  out[0] = 1.01589817494786*fix0*tmp10*pow(b_a,0.75); // final
  tmp8 = tmp13 + tmp8; // oblige+recycle
  tmp10 = 2*tmp10; // auto+recycle
  tmp2 = -2*tmp2*gaux(fix1, 1); // auto+recycle
  fix1 = fix0*pow(b_a,1.25); // auto+recycle
  out[1] = 1.01589817494786*fix1*(fix3*tmp2 + tmp10*tmp8); // final
  tmp12 = tmp12 + tmp7; // oblige+recycle
  out[2] = 1.01589817494786*fix1*(fix4*tmp2 + tmp10*tmp12); // final
  tmp12 = tmp11 + tmp6; // oblige+recycle
  out[3] = 1.01589817494786*fix1*(fix5*tmp2 + tmp10*tmp12); // final
}

static void gint2_nai_P_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  double fix0, fix1, fix2, fix3, fix4, fix5, fix6, tmp10, tmp14, tmp15, tmp16, tmp17, tmp23, tmp24, tmp25, tmp26, tmp27, tmp6;
  tmp26 = -b[0]; // auto
  fix0 = a[0] + tmp26; // oblige
  tmp25 = -b[1]; // auto
  fix1 = a[1] + tmp25; // oblige
  tmp24 = -b[2]; // auto
  fix2 = a[2] + tmp24; // oblige
  tmp23 = a_a + b_a; // auto
  tmp6 = 1.0/tmp23; // auto
  tmp16 = tmp6*(a_a*a[0] + b_a*b[0]); // auto
  fix3 = tmp16 - c[0]; // oblige
  tmp15 = tmp6*(a_a*a[1] + b_a*b[1]); // auto
  fix4 = tmp15 - c[1]; // oblige
  tmp14 = tmp6*(a_a*a[2] + b_a*b[2]); // auto
  fix5 = tmp14 - c[2]; // oblige
  fix6 = tmp16 - a[0]; // oblige
  fix1 = M_PI*tmp6*exp(-a_a*b_a*tmp6*(fix0*fix0 + fix1*fix1 + fix2*fix2)); // auto+recycle
  tmp6 = tmp23*(fix3*fix3 + fix4*fix4 + fix5*fix5); // auto+recycle
  fix0 = fix1*gaux(tmp6, 1); // auto+recycle
  tmp23 = 2*fix1*gaux(tmp6, 0); // auto+recycle
  fix2 = -2*fix0; // auto+recycle
  tmp10 = fix2*fix3 + fix6*tmp23; // auto
  tmp27 = pow(a_a,1.25); // auto
  tmp17 = tmp27*pow(b_a,0.75); // auto
  out[0] = 1.01589817494786*tmp10*tmp17; // final
  tmp26 = tmp16 + tmp26; // oblige+recycle
  tmp6 = -2*fix1*gaux(tmp6, 2); // auto+recycle
  fix1 = tmp27*pow(b_a,1.25); // auto+recycle
  tmp27 = (fix2 + tmp23)/(2*a_a + 2*b_a); // auto+recycle
  fix0 = 2*fix0; // auto+recycle
  fix6 = fix0*fix6 + fix3*tmp6; // auto+recycle
  fix6 = -fix6; // auto+recycle
  out[1] = 2.03179634989571*fix1*(tmp27 + fix3*fix6 + tmp10*tmp26); // final
  tmp25 = tmp15 + tmp25; // oblige+recycle
  out[2] = 2.03179634989571*fix1*(fix4*fix6 + tmp10*tmp25); // final
  tmp24 = tmp14 + tmp24; // oblige+recycle
  out[3] = 2.03179634989571*fix1*(fix5*fix6 + tmp10*tmp24); // final
  tmp15 = tmp15 - a[1]; // oblige+recycle
  fix6 = fix2*fix4 + tmp15*tmp23; // auto+recycle
  out[4] = 1.01589817494786*fix6*tmp17; // final
  tmp15 = fix0*tmp15 + fix4*tmp6; // auto+recycle
  tmp15 = -tmp15; // auto+recycle
  out[5] = 2.03179634989571*fix1*(fix3*tmp15 + fix6*tmp26); // final
  out[6] = 2.03179634989571*fix1*(tmp27 + fix4*tmp15 + fix6*tmp25); // final
  out[7] = 2.03179634989571*fix1*(fix5*tmp15 + fix6*tmp24); // final
  tmp15 = tmp14 - a[2]; // oblige+recycle
  tmp14 = fix2*fix5 + tmp15*tmp23; // auto+recycle
  out[8] = 1.01589817494786*tmp14*tmp17; // final
  tmp15 = fix0*tmp15 + fix5*tmp6; // auto+recycle
  tmp15 = -tmp15; // auto+recycle
  out[9] = 2.03179634989571*fix1*(fix3*tmp15 + tmp14*tmp26); // final
  out[10] = 2.03179634989571*fix1*(fix4*tmp15 + tmp14*tmp25); // final
  out[11] = 2.03179634989571*fix1*(tmp27 + fix5*tmp15 + tmp14*tmp24); // final
}

static void gint2_nai_SP_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  double fix0, fix1, fix2, fix3, fix4, fix5, tmp10, tmp2, tmp6, tmp7, tmp8;
  fix0 = a[0] - b[0]; // oblige
  fix1 = a[1] - b[1]; // oblige
  fix2 = a[2] - b[2]; // oblige
  tmp10 = a_a + b_a; // auto
  tmp2 = 1.0/tmp10; // auto
  tmp8 = tmp2*(a_a*a[0] + b_a*b[0]); // auto
  fix3 = tmp8 - c[0]; // oblige
  tmp7 = tmp2*(a_a*a[1] + b_a*b[1]); // auto
  fix4 = tmp7 - c[1]; // oblige
  tmp6 = tmp2*(a_a*a[2] + b_a*b[2]); // auto
  fix5 = tmp6 - c[2]; // oblige
  tmp2 = M_PI*tmp2*exp(-a_a*b_a*tmp2*(fix0*fix0 + fix1*fix1 + fix2*fix2)); // auto+recycle
  fix1 = tmp10*(fix3*fix3 + fix4*fix4 + fix5*fix5); // auto+recycle
  tmp10 = tmp2*gaux(fix1, 0); // auto+recycle
  fix0 = pow(b_a,0.75); // auto+recycle
  out[0] = 1.01589817494786*fix0*tmp10*pow(a_a,0.75); // final
  tmp8 = tmp8 - a[0]; // oblige+recycle
  tmp10 = 2*tmp10; // auto+recycle
  tmp2 = -2*tmp2*gaux(fix1, 1); // auto+recycle
  fix1 = fix0*pow(a_a,1.25); // auto+recycle
  out[1] = 1.01589817494786*fix1*(fix3*tmp2 + tmp10*tmp8); // final
  fix0 = tmp7 - a[1]; // oblige+recycle
  out[2] = 1.01589817494786*fix1*(fix0*tmp10 + fix4*tmp2); // final
  fix4 = tmp6 - a[2]; // oblige+recycle
  out[3] = 1.01589817494786*fix1*(fix4*tmp10 + fix5*tmp2); // final
}

static void gint2_nai_S_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  double fix0, fix1, fix2, fix3, fix4, fix5, tmp0, tmp1;
  fix0 = a[0] - b[0]; // oblige
  fix1 = a[1] - b[1]; // oblige
  fix2 = a[2] - b[2]; // oblige
  tmp1 = a_a + b_a; // auto
  tmp0 = 1.0/tmp1; // auto
  fix3 = -c[0] + tmp0*(a_a*a[0] + b_a*b[0]); // oblige
  fix4 = -c[1] + tmp0*(a_a*a[1] + b_a*b[1]); // oblige
  fix5 = -c[2] + tmp0*(a_a*a[2] + b_a*b[2]); // oblige
  out[0] = 1.01589817494786*M_PI*tmp0*pow(a_a,0.75)*pow(b_a,0.75)*exp(-a_a*b_a*tmp0*(fix0*fix0 + fix1*fix1 + fix2*fix2))*gaux(tmp1*(fix3*fix3 + fix4*fix4 + fix5*fix5), 0); // final
}

static void gint2_nai_P_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  double fix0, fix1, fix2, fix3, fix4, fix5, tmp2, tmp6, tmp7, tmp8, tmp9;
  fix0 = a[0] - b[0]; // oblige
  fix1 = a[1] - b[1]; // oblige
  fix2 = a[2] - b[2]; // oblige
  tmp9 = a_a + b_a; // auto
  tmp2 = 1.0/tmp9; // auto
  tmp8 = tmp2*(a_a*a[0] + b_a*b[0]); // auto
  fix3 = tmp8 - c[0]; // oblige
  tmp7 = tmp2*(a_a*a[1] + b_a*b[1]); // auto
  fix4 = tmp7 - c[1]; // oblige
  tmp6 = tmp2*(a_a*a[2] + b_a*b[2]); // auto
  fix5 = tmp6 - c[2]; // oblige
  tmp8 = tmp8 - a[0]; // oblige+recycle
  tmp2 = M_PI*tmp2*exp(-a_a*b_a*tmp2*(fix0*fix0 + fix1*fix1 + fix2*fix2)); // auto+recycle
  fix1 = tmp9*(fix3*fix3 + fix4*fix4 + fix5*fix5); // auto+recycle
  tmp9 = 2*tmp2*gaux(fix1, 0); // auto+recycle
  tmp2 = -2*tmp2*gaux(fix1, 1); // auto+recycle
  fix1 = pow(a_a,1.25)*pow(b_a,0.75); // auto+recycle
  out[0] = 1.01589817494786*fix1*(fix3*tmp2 + tmp8*tmp9); // final
  fix0 = tmp7 - a[1]; // oblige+recycle
  out[1] = 1.01589817494786*fix1*(fix0*tmp9 + fix4*tmp2); // final
  fix4 = tmp6 - a[2]; // oblige+recycle
  out[2] = 1.01589817494786*fix1*(fix4*tmp9 + fix5*tmp2); // final
}

static void gint2_nai_SP_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  double fix0, fix1, fix2, fix3, fix4, fix5, tmp12, tmp13, tmp14, tmp18, tmp23, tmp24, tmp25, tmp26, tmp27, tmp28, tmp29, tmp30, tmp6;
  tmp29 = -b[0]; // auto
  fix0 = a[0] + tmp29; // oblige
  tmp28 = -b[1]; // auto
  fix1 = a[1] + tmp28; // oblige
  tmp27 = -b[2]; // auto
  fix2 = a[2] + tmp27; // oblige
  tmp23 = a_a + b_a; // auto
  tmp6 = 1.0/tmp23; // auto
  tmp14 = tmp6*(a_a*a[0] + b_a*b[0]); // auto
  fix3 = tmp14 - c[0]; // oblige
  tmp13 = tmp6*(a_a*a[1] + b_a*b[1]); // auto
  fix4 = tmp13 - c[1]; // oblige
  tmp12 = tmp6*(a_a*a[2] + b_a*b[2]); // auto
  fix5 = tmp12 - c[2]; // oblige
  tmp29 = tmp14 + tmp29; // oblige+recycle
  fix1 = M_PI*tmp6*exp(-a_a*b_a*tmp6*(fix0*fix0 + fix1*fix1 + fix2*fix2)); // auto+recycle
  fix0 = tmp23*(fix3*fix3 + fix4*fix4 + fix5*fix5); // auto+recycle
  tmp23 = fix1*gaux(fix0, 1); // auto+recycle
  fix2 = 2*fix1*gaux(fix0, 0); // auto+recycle
  tmp6 = -2*tmp23; // auto+recycle
  tmp30 = pow(b_a,1.25); // auto
  tmp18 = tmp30*pow(a_a,0.75); // auto
  tmp26 = fix3*tmp6; // auto
  out[0] = 1.01589817494786*tmp18*(tmp26 + fix2*tmp29); // final
  tmp28 = tmp13 + tmp28; // oblige+recycle
  tmp25 = fix4*tmp6; // auto
  out[1] = 1.01589817494786*tmp18*(tmp25 + fix2*tmp28); // final
  tmp27 = tmp12 + tmp27; // oblige+recycle
  tmp24 = fix5*tmp6; // auto
  out[2] = 1.01589817494786*tmp18*(tmp24 + fix2*tmp27); // final
  tmp14 = tmp14 - a[0]; // oblige+recycle
  fix1 = -2*fix1*gaux(fix0, 2); // auto+recycle
  fix0 = tmp30*pow(a_a,1.25); // auto+recycle
  tmp6 = (fix2 + tmp6)/(2*a_a + 2*b_a); // auto+recycle
  tmp23 = 2*tmp23; // auto+recycle
  tmp18 = fix1*fix3 + tmp14*tmp23; // auto+recycle
  tmp14 = tmp26 + fix2*tmp14; // auto+recycle
  tmp26 = -tmp18; // auto+recycle
  out[3] = 2.03179634989571*fix0*(tmp6 + fix3*tmp26 + tmp14*tmp29); // final
  out[4] = 2.03179634989571*fix0*(fix4*tmp26 + tmp14*tmp28); // final
  out[5] = 2.03179634989571*fix0*(fix5*tmp26 + tmp14*tmp27); // final
  tmp14 = tmp13 - a[1]; // oblige+recycle
  tmp13 = fix1*fix4 + tmp14*tmp23; // auto+recycle
  tmp25 = tmp25 + fix2*tmp14; // auto+recycle
  tmp14 = -tmp13; // auto+recycle
  out[6] = 2.03179634989571*fix0*(fix3*tmp14 + tmp25*tmp29); // final
  out[7] = 2.03179634989571*fix0*(tmp6 + fix4*tmp14 + tmp25*tmp28); // final
  out[8] = 2.03179634989571*fix0*(fix5*tmp14 + tmp25*tmp27); // final
  tmp25 = tmp12 - a[2]; // oblige+recycle
  tmp14 = fix1*fix5 + tmp23*tmp25; // auto+recycle
  tmp25 = tmp24 + fix2*tmp25; // auto+recycle
  tmp14 = -tmp14; // auto+recycle
  out[9] = 2.03179634989571*fix0*(fix3*tmp14 + tmp25*tmp29); // final
  out[10] = 2.03179634989571*fix0*(fix4*tmp14 + tmp25*tmp28); // final
  out[11] = 2.03179634989571*fix0*(tmp6 + fix5*tmp14 + tmp25*tmp27); // final
}

static void gint2_nai_S_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  double fix0, fix1, fix2, fix3, fix4, fix5, tmp10, tmp11, tmp12, tmp2, tmp6, tmp7, tmp8, tmp9;
  tmp12 = -b[0]; // auto
  fix0 = a[0] + tmp12; // oblige
  tmp11 = -b[1]; // auto
  fix1 = a[1] + tmp11; // oblige
  tmp10 = -b[2]; // auto
  fix2 = a[2] + tmp10; // oblige
  tmp9 = a_a + b_a; // auto
  tmp2 = 1.0/tmp9; // auto
  tmp8 = tmp2*(a_a*a[0] + b_a*b[0]); // auto
  fix3 = tmp8 - c[0]; // oblige
  tmp7 = tmp2*(a_a*a[1] + b_a*b[1]); // auto
  fix4 = tmp7 - c[1]; // oblige
  tmp6 = tmp2*(a_a*a[2] + b_a*b[2]); // auto
  fix5 = tmp6 - c[2]; // oblige
  tmp12 = tmp12 + tmp8; // oblige+recycle
  tmp2 = M_PI*tmp2*exp(-a_a*b_a*tmp2*(fix0*fix0 + fix1*fix1 + fix2*fix2)); // auto+recycle
  fix1 = tmp9*(fix3*fix3 + fix4*fix4 + fix5*fix5); // auto+recycle
  tmp9 = 2*tmp2*gaux(fix1, 0); // auto+recycle
  tmp2 = -2*tmp2*gaux(fix1, 1); // auto+recycle
  fix1 = pow(a_a,0.75)*pow(b_a,1.25); // auto+recycle
  out[0] = 1.01589817494786*fix1*(fix3*tmp2 + tmp12*tmp9); // final
  tmp12 = tmp11 + tmp7; // oblige+recycle
  out[1] = 1.01589817494786*fix1*(fix4*tmp2 + tmp12*tmp9); // final
  tmp12 = tmp10 + tmp6; // oblige+recycle
  out[2] = 1.01589817494786*fix1*(fix5*tmp2 + tmp12*tmp9); // final
}

static void gint2_nai_P_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  double fix0, fix1, fix2, fix3, fix4, fix5, fix6, tmp11, tmp14, tmp15, tmp16, tmp22, tmp23, tmp24, tmp25, tmp6;
  tmp25 = -b[0]; // auto
  fix0 = a[0] + tmp25; // oblige
  tmp24 = -b[1]; // auto
  fix1 = a[1] + tmp24; // oblige
  tmp23 = -b[2]; // auto
  fix2 = a[2] + tmp23; // oblige
  tmp22 = a_a + b_a; // auto
  tmp6 = 1.0/tmp22; // auto
  tmp16 = tmp6*(a_a*a[0] + b_a*b[0]); // auto
  fix3 = tmp16 - c[0]; // oblige
  tmp15 = tmp6*(a_a*a[1] + b_a*b[1]); // auto
  fix4 = tmp15 - c[1]; // oblige
  tmp14 = tmp6*(a_a*a[2] + b_a*b[2]); // auto
  fix5 = tmp14 - c[2]; // oblige
  fix6 = tmp16 - a[0]; // oblige
  tmp25 = tmp16 + tmp25; // oblige+recycle
  fix1 = M_PI*tmp6*exp(-a_a*b_a*tmp6*(fix0*fix0 + fix1*fix1 + fix2*fix2)); // auto+recycle
  tmp6 = tmp22*(fix3*fix3 + fix4*fix4 + fix5*fix5); // auto+recycle
  fix0 = fix1*gaux(tmp6, 1); // auto+recycle
  tmp22 = 2*fix1*gaux(tmp6, 0); // auto+recycle
  fix1 = -2*fix1*gaux(tmp6, 2); // auto+recycle
  tmp6 = pow(a_a,1.25)*pow(b_a,1.25); // auto+recycle
  tmp16 = -2*fix0; // auto+recycle
  fix2 = (tmp16 + tmp22)/(2*a_a + 2*b_a); // auto+recycle
  tmp11 = fix3*tmp16 + fix6*tmp22; // auto
  fix0 = 2*fix0; // auto+recycle
  fix6 = fix0*fix6 + fix1*fix3; // auto+recycle
  fix6 = -fix6; // auto+recycle
  out[0] = 2.03179634989571*tmp6*(fix2 + fix3*fix6 + tmp11*tmp25); // final
  tmp24 = tmp15 + tmp24; // oblige+recycle
  out[1] = 2.03179634989571*tmp6*(fix4*fix6 + tmp11*tmp24); // final
  tmp23 = tmp14 + tmp23; // oblige+recycle
  out[2] = 2.03179634989571*tmp6*(fix5*fix6 + tmp11*tmp23); // final
  tmp15 = tmp15 - a[1]; // oblige+recycle
  fix6 = fix0*tmp15 + fix1*fix4; // auto+recycle
  tmp15 = fix4*tmp16 + tmp15*tmp22; // auto+recycle
  fix6 = -fix6; // auto+recycle
  out[3] = 2.03179634989571*tmp6*(fix3*fix6 + tmp15*tmp25); // final
  out[4] = 2.03179634989571*tmp6*(fix2 + fix4*fix6 + tmp15*tmp24); // final
  out[5] = 2.03179634989571*tmp6*(fix5*fix6 + tmp15*tmp23); // final
  tmp15 = tmp14 - a[2]; // oblige+recycle
  tmp14 = fix5*tmp16 + tmp15*tmp22; // auto+recycle
  tmp15 = fix0*tmp15 + fix1*fix5; // auto+recycle
  tmp15 = -tmp15; // auto+recycle
  out[6] = 2.03179634989571*tmp6*(fix3*tmp15 + tmp14*tmp25); // final
  out[7] = 2.03179634989571*tmp6*(fix4*tmp15 + tmp14*tmp24); // final
  out[8] = 2.03179634989571*tmp6*(fix2 + fix5*tmp15 + tmp14*tmp23); // final
}

typedef void (*fntype)(double*, double, double*, double, double*, double*);
const fntype fns[9] = {gint2_nai_SP_SP, gint2_nai_S_SP, gint2_nai_P_SP, gint2_nai_SP_S, gint2_nai_S_S, gint2_nai_P_S, gint2_nai_SP_P, gint2_nai_S_P, gint2_nai_P_P};

void gint2_nai_dispatch(int a_s, double* a, double a_a, int b_s, double* b, double b_a, double* c, double* out)
{
  fns[1+a_s+3*(1+b_s)](a, a_a, b, b_a, c, out);
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
