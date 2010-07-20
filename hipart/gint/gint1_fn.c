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
#define MAX_SHELL 5
#define NUM_SHELL_TYPES 11
#define MAX_SHELL_DOF 21

static void fn_cH(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17;
  tmp0 = p[0] - a[0];
  tmp1 = p[1] - a[1];
  tmp2 = tmp1*tmp1;
  tmp3 = tmp0*tmp0;
  tmp4 = p[2] - a[2];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp2 + tmp3 + tmp5;
  tmp7 = -a_a*tmp6;
  tmp8 = exp(tmp7);
  tmp9 = pow(tmp3,(5.0/2.0));
  tmp10 = pow(tmp9,(4.0/5.0));
  tmp11 = pow(tmp10,(3.0/4.0));
  tmp12 = pow(tmp2,(3.0/2.0));
  tmp13 = pow(tmp5,(3.0/2.0));
  tmp14 = pow(tmp12,(4.0/3.0));
  tmp15 = pow(tmp13,(4.0/3.0));
  tmp16 = pow(tmp14,(5.0/4.0));
  tmp17 = pow(tmp15,(5.0/4.0));
  out[0] = tmp8*tmp11*tmp11;
  out[1] = tmp1*tmp8*tmp9;
  out[2] = tmp4*tmp8*tmp9;
  out[3] = tmp10*tmp2*tmp8;
  out[4] = tmp1*tmp10*tmp4*tmp8;
  out[5] = tmp10*tmp5*tmp8;
  out[6] = tmp11*tmp12*tmp8;
  out[7] = tmp11*tmp2*tmp4*tmp8;
  out[8] = tmp1*tmp11*tmp5*tmp8;
  out[9] = tmp11*tmp4*tmp5*tmp8;
  out[10] = tmp14*tmp3*tmp8;
  out[11] = tmp12*tmp3*tmp4*tmp8;
  out[12] = tmp2*tmp3*tmp5*tmp8;
  out[13] = tmp1*tmp3*tmp4*tmp5*tmp8;
  out[14] = tmp15*tmp3*tmp8;
  out[15] = tmp0*tmp16*tmp8;
  out[16] = tmp0*tmp14*tmp4*tmp8;
  out[17] = tmp0*tmp12*tmp5*tmp8;
  out[18] = tmp0*tmp2*tmp4*tmp5*tmp8;
  out[19] = tmp0*tmp1*tmp15*tmp8;
  out[20] = tmp0*tmp15*tmp4*tmp8;
  out[21] = tmp8*pow(tmp16,(6.0/5.0));
  out[22] = tmp16*tmp4*tmp8;
  out[23] = tmp14*tmp5*tmp8;
  out[24] = tmp12*tmp4*tmp5*tmp8;
  out[25] = tmp15*tmp2*tmp8;
  out[26] = tmp1*tmp15*tmp4*tmp8;
  out[27] = tmp8*pow(tmp17,(6.0/5.0));
}

static void fn_cG(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14;
  tmp0 = p[0] - a[0];
  tmp1 = p[1] - a[1];
  tmp2 = tmp1*tmp1;
  tmp3 = tmp0*tmp0;
  tmp4 = p[2] - a[2];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp2 + tmp3 + tmp5;
  tmp7 = -a_a*tmp6;
  tmp8 = exp(tmp7);
  tmp9 = tmp3*tmp3;
  tmp10 = pow(tmp9,(3.0/4.0));
  tmp11 = pow(tmp2,(3.0/2.0));
  tmp12 = pow(tmp5,(3.0/2.0));
  tmp13 = pow(tmp11,(4.0/3.0));
  tmp14 = pow(tmp12,(4.0/3.0));
  out[0] = tmp8*pow(tmp10,(5.0/3.0));
  out[1] = tmp1*tmp8*tmp9;
  out[2] = tmp4*tmp8*tmp9;
  out[3] = tmp10*tmp2*tmp8;
  out[4] = tmp1*tmp10*tmp4*tmp8;
  out[5] = tmp10*tmp5*tmp8;
  out[6] = tmp11*tmp3*tmp8;
  out[7] = tmp2*tmp3*tmp4*tmp8;
  out[8] = tmp1*tmp3*tmp5*tmp8;
  out[9] = tmp3*tmp4*tmp5*tmp8;
  out[10] = tmp0*tmp13*tmp8;
  out[11] = tmp0*tmp11*tmp4*tmp8;
  out[12] = tmp0*tmp2*tmp5*tmp8;
  out[13] = tmp0*tmp1*tmp4*tmp5*tmp8;
  out[14] = tmp0*tmp14*tmp8;
  out[15] = tmp8*pow(tmp13,(5.0/4.0));
  out[16] = tmp13*tmp4*tmp8;
  out[17] = tmp11*tmp5*tmp8;
  out[18] = tmp2*tmp4*tmp5*tmp8;
  out[19] = tmp1*tmp14*tmp8;
  out[20] = tmp14*tmp4*tmp8;
}

static void fn_cF(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11;
  tmp0 = p[0] - a[0];
  tmp1 = p[1] - a[1];
  tmp2 = tmp1*tmp1;
  tmp3 = tmp0*tmp0;
  tmp4 = p[2] - a[2];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp2 + tmp3 + tmp5;
  tmp7 = -a_a*tmp6;
  tmp8 = exp(tmp7);
  tmp9 = pow(tmp3,(3.0/2.0));
  tmp10 = pow(tmp2,(3.0/2.0));
  tmp11 = pow(tmp5,(3.0/2.0));
  out[0] = tmp8*pow(tmp9,(4.0/3.0));
  out[1] = tmp1*tmp8*tmp9;
  out[2] = tmp4*tmp8*tmp9;
  out[3] = tmp2*tmp3*tmp8;
  out[4] = tmp1*tmp3*tmp4*tmp8;
  out[5] = tmp3*tmp5*tmp8;
  out[6] = tmp0*tmp10*tmp8;
  out[7] = tmp0*tmp2*tmp4*tmp8;
  out[8] = tmp0*tmp1*tmp5*tmp8;
  out[9] = tmp0*tmp4*tmp5*tmp8;
  out[10] = tmp8*pow(tmp10,(4.0/3.0));
  out[11] = tmp10*tmp4*tmp8;
  out[12] = tmp2*tmp5*tmp8;
  out[13] = tmp1*tmp4*tmp5*tmp8;
  out[14] = tmp8*pow(tmp11,(4.0/3.0));
}

static void fn_cD(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  tmp0 = p[0] - a[0];
  tmp1 = p[1] - a[1];
  tmp2 = tmp0*tmp0;
  tmp3 = tmp1*tmp1;
  tmp4 = p[2] - a[2];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp2 + tmp3 + tmp5;
  tmp7 = -a_a*tmp6;
  tmp8 = exp(tmp7);
  out[0] = tmp8*pow(tmp2,(3.0/2.0));
  out[1] = tmp1*tmp2*tmp8;
  out[2] = tmp2*tmp4*tmp8;
  out[3] = tmp0*tmp3*tmp8;
  out[4] = tmp0*tmp1*tmp4*tmp8;
  out[5] = tmp0*tmp5*tmp8;
  out[6] = tmp8*pow(tmp3,(3.0/2.0));
  out[7] = tmp3*tmp4*tmp8;
  out[8] = tmp1*tmp5*tmp8;
  out[9] = tmp4*tmp5*tmp8;
}

static void fn_SP(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  tmp0 = p[0] - a[0];
  tmp1 = p[1] - a[1];
  tmp2 = tmp1*tmp1;
  tmp3 = tmp0*tmp0;
  tmp4 = p[2] - a[2];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp2 + tmp3 + tmp5;
  tmp7 = -a_a*tmp6;
  tmp8 = exp(tmp7);
  out[0] = 0.282094791773878*tmp8;
  out[1] = 0.48860251190292*tmp0*tmp8;
  out[2] = 0.48860251190292*tmp1*tmp8;
  out[3] = 0.48860251190292*tmp4*tmp8;
}

static void fn_S(double* a, double a_a, double* p, double* out)
{
  out[0] = 0.282094791773878*exp(-a_a*(pow((p[0] - a[0]),2) + pow((p[1] - a[1]),2) + pow((p[2] - a[2]),2)));
}

static void fn_P(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  tmp0 = p[0] - a[0];
  tmp1 = p[1] - a[1];
  tmp2 = tmp1*tmp1;
  tmp3 = tmp0*tmp0;
  tmp4 = p[2] - a[2];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp2 + tmp3 + tmp5;
  tmp7 = -a_a*tmp6;
  tmp8 = exp(tmp7);
  out[0] = 0.48860251190292*tmp0*tmp8;
  out[1] = 0.48860251190292*tmp1*tmp8;
  out[2] = 0.48860251190292*tmp4*tmp8;
}

static void fn_pD(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11;
  tmp0 = p[0] - a[0];
  tmp1 = tmp0*tmp0;
  tmp2 = p[1] - a[1];
  tmp3 = tmp2*tmp2;
  tmp4 = p[2] - a[2];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp1 + tmp3 + tmp5;
  tmp7 = -a_a*tmp6;
  tmp8 = exp(tmp7);
  tmp9 = -1.09254843059208*a[2];
  tmp10 = 1.09254843059208*p[2];
  tmp11 = tmp10 + tmp9;
  out[0] = tmp8*(0.54627421529604*tmp1 - 0.54627421529604*tmp3);
  out[1] = tmp0*tmp8*(1.09254843059208*p[1] - 1.09254843059208*a[1]);
  out[2] = tmp0*tmp11*tmp8;
  out[3] = tmp11*tmp2*tmp8;
  out[4] = -tmp8*(0.31539156525252 - 0.94617469575756*tmp5);
}

static void fn_pF(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  tmp0 = p[0] - a[0];
  tmp1 = p[1] - a[1];
  tmp2 = tmp1*tmp1;
  tmp3 = tmp0*tmp0;
  tmp4 = p[2] - a[2];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp2 + tmp3 + tmp5;
  tmp7 = -a_a*tmp6;
  tmp8 = exp(tmp7);
  out[0] = tmp8*(0.590043589926644*tmp0*tmp3 - 1.77013076977993*tmp0*tmp2);
  out[1] = tmp8*(1.77013076977993*tmp1*tmp3 - 0.590043589926644*tmp1*tmp2);
  out[2] = tmp8*(-tmp2*(1.44530572132028*p[2] - 1.44530572132028*a[2]) + 1.44530572132028*tmp3*tmp4);
  out[3] = tmp0*tmp1*tmp8*(2.89061144264055*p[2] - 2.89061144264055*a[2]);
  out[4] = tmp8*(0.457045799464466*a[0] - 0.457045799464466*p[0] + tmp5*(2.28522899732233*p[0] - 2.28522899732233*a[0]));
  out[5] = tmp8*(0.457045799464466*a[1] - 0.457045799464466*p[1] + tmp5*(2.28522899732233*p[1] - 2.28522899732233*a[1]));
  out[6] = tmp8*(1.11952899777035*a[2] - 1.11952899777035*p[2] + 1.86588166295058*tmp4*tmp5);
}

static void fn_pG(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17;
  tmp0 = p[0] - a[0];
  tmp1 = p[1] - a[1];
  tmp2 = tmp0*tmp0;
  tmp3 = tmp1*tmp1;
  tmp4 = p[2] - a[2];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp2 + tmp3 + tmp5;
  tmp7 = -a_a*tmp6;
  tmp8 = exp(tmp7);
  tmp9 = pow(tmp2,(3.0/2.0));
  tmp10 = -1.77013076977993*a[2];
  tmp11 = 1.77013076977993*p[2];
  tmp12 = tmp10 + tmp11;
  tmp13 = pow(tmp3,(3.0/2.0));
  tmp14 = pow(tmp5,(3.0/2.0));
  tmp15 = 2.00713963067187*p[2];
  tmp16 = -2.00713963067187*a[2];
  tmp17 = tmp15 + tmp16;
  out[0] = tmp8*(-3.75501441269506*tmp2*tmp3 + 0.625835735449176*pow(tmp13,(4.0/3.0)) + 0.625835735449176*pow(tmp9,(4.0/3.0)));
  out[1] = tmp8*(tmp0*tmp2*(2.5033429417967*p[1] - 2.5033429417967*a[1]) - 2.5033429417967*tmp0*tmp1*tmp3);
  out[2] = tmp8*(tmp0*tmp12*tmp2 - 5.31039230933979*tmp0*tmp3*tmp4);
  out[3] = tmp8*(-tmp1*tmp12*tmp3 + 5.31039230933979*tmp1*tmp2*tmp4);
  out[4] = tmp8*(0.47308734787878*tmp3 - 0.47308734787878*tmp2 + 3.31161143515146*tmp2*tmp5 - 3.31161143515146*tmp3*tmp5);
  out[5] = tmp8*(-tmp0*(0.94617469575756*p[1] - 0.94617469575756*a[1]) + tmp0*tmp5*(6.62322287030292*p[1] - 6.62322287030292*a[1]));
  out[6] = tmp8*(-tmp0*tmp17 + 4.68332580490102*tmp0*tmp4*tmp5);
  out[7] = tmp8*(-tmp1*tmp17 + 4.68332580490102*tmp1*tmp4*tmp5);
  out[8] = tmp8*(0.317356640745613 - 3.17356640745613*tmp5 + 3.70249414203215*pow(tmp14,(4.0/3.0)));
}

static void fn_pH(double* a, double a_a, double* p, double* out)
{
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp18, tmp19, tmp20;
  tmp0 = p[1] - a[1];
  tmp1 = p[0] - a[0];
  tmp2 = tmp0*tmp0;
  tmp3 = tmp1*tmp1;
  tmp4 = p[2] - a[2];
  tmp5 = tmp4*tmp4;
  tmp6 = tmp2 + tmp3 + tmp5;
  tmp7 = -a_a*tmp6;
  tmp8 = exp(tmp7);
  tmp9 = tmp2*tmp2;
  tmp10 = 2.07566231488104*p[2];
  tmp11 = -2.07566231488104*a[2];
  tmp12 = tmp10 + tmp11;
  tmp13 = tmp3*tmp3;
  tmp14 = pow(tmp13,(3.0/4.0));
  tmp15 = 8.30264925952416*p[2];
  tmp16 = -8.30264925952416*a[2];
  tmp17 = tmp15 + tmp16;
  tmp18 = pow(tmp9,(3.0/4.0));
  tmp19 = pow(tmp5,(3.0/2.0));
  tmp20 = pow(tmp19,(4.0/3.0));
  out[0] = tmp8*(tmp9*(3.28191028420085*p[0] - 3.28191028420085*a[0]) - 6.5638205684017*tmp1*tmp2*tmp3 + 0.65638205684017*pow(tmp14,(5.0/3.0)));
  out[1] = tmp8*(tmp13*(3.28191028420085*p[1] - 3.28191028420085*a[1]) - 6.5638205684017*tmp0*tmp2*tmp3 + 0.65638205684017*pow(tmp18,(5.0/3.0)));
  out[2] = tmp8*(tmp12*tmp13 + tmp12*tmp9 - tmp2*tmp3*(12.4539738892862*p[2] - 12.4539738892862*a[2]));
  out[3] = tmp8*(tmp0*tmp1*tmp17*tmp3 - tmp0*tmp1*tmp17*tmp2);
  out[4] = tmp8*(tmp2*(1.46771489830575*p[0] - 1.46771489830575*a[0]) + 4.40314469491725*tmp14*tmp5 - 0.48923829943525*tmp1*tmp3 - 13.2094340847518*tmp1*tmp2*tmp5);
  out[5] = tmp8*(-tmp3*(1.46771489830575*p[1] - 1.46771489830575*a[1]) + 0.48923829943525*tmp0*tmp2 - 4.40314469491725*tmp18*tmp5 + tmp3*tmp5*(13.2094340847518*p[1] - 13.2094340847518*a[1]));
  out[6] = tmp8*(-tmp3*(2.39676839248666*p[2] - 2.39676839248666*a[2]) + 2.39676839248666*tmp2*tmp4 + 7.19030517745999*tmp3*tmp4*tmp5 - 7.19030517745999*tmp2*tmp4*tmp5);
  out[7] = tmp8*(-tmp1*tmp4*(4.79353678497332*p[1] - 4.79353678497332*a[1]) + tmp1*tmp4*tmp5*(14.38061035492*p[1] - 14.38061035492*a[1]));
  out[8] = tmp8*(0.452946651195697*p[0] - 0.452946651195697*a[0] + tmp20*(9.51187967510963*p[0] - 9.51187967510963*a[0]) - tmp5*(6.34125311673976*p[0] - 6.34125311673976*a[0]));
  out[9] = tmp8*(0.452946651195697*p[1] - 0.452946651195697*a[1] + tmp20*(9.51187967510963*p[1] - 9.51187967510963*a[1]) - tmp5*(6.34125311673976*p[1] - 6.34125311673976*a[1]));
  out[10] = tmp8*(1.75425483680135*p[2] - 1.75425483680135*a[2] - 8.18652257173965*tmp4*tmp5 + 7.36787031456569*pow(tmp20,(5.0/4.0)));
}

typedef void (*fntype)(double*, double, double*, double*);
const fntype fns[11] = {fn_cH, fn_cG, fn_cF, fn_cD, fn_SP, fn_S, fn_P, fn_pD, fn_pF, fn_pG, fn_pH};

void gint1_fn(int a_s, double* a, double a_a, double* p, double* out)
{
  fns[5+a_s](a, a_a, p, out);
}

int gint1_fn_basis(double* weights, double* fn, double* point, double* centers,
  int* shell_types, int* shell_map,  int* num_primitives, double* ccoeffs,
  double* exponents, int num_shells, int num_centers, int num_weights,
  int num_ccoeffs, int num_exponents)
{
  int shell, shell_type, primitive, dof, num_dof, result;
  double *work, *out, *center, *ccoeff, *exponent, *weight;

  result = 0;
  work = malloc(MAX_SHELL_DOF*sizeof(double));
  if (work==NULL) {result = -1; goto EXIT;}

  *fn = 0.0;
  weight = weights;
  ccoeff = ccoeffs;
  exponent = exponents;

  for (shell=0; shell<num_shells; shell++) {
    center = centers + (3*shell_map[shell]);
    shell_type = shell_types[shell];
    for (primitive=0; primitive<num_primitives[shell]; primitive++) {
      gint1_fn(shell_type, center, *exponent, point, work);
      out = work;
      exponent++;
      if (shell_type==-1) {
        *fn += (*weight)*(*out)*(*ccoeff);
        weight++;
        out++;
        ccoeff++;
        num_dof = 3;
      } else if (shell_type >= 0) {
        num_dof = 2*shell_type+1;
      } else {
        num_dof = ((shell_type+1)*(shell_type+2))/2;
      }
      for (dof=0; dof<num_dof; dof++) {
        *fn += (*weight)*(*out)*(*ccoeff);
        weight++;
        out++;
      }
      ccoeff++;
    }
  }

EXIT:
  free(work);
  return result;
}
