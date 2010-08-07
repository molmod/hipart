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
#include "gint2_nai_auto.h"
#include "gaux.h"

void gint2_nai_pF_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 12
  double tmp0, tmp1, tmp10, tmp11, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  gint2_nai_cF_pF(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[35] - 0.273861278752583*out[21] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = out[63] - 0.670820393249937*out[14] - 0.670820393249937*out[49]; // stacked, weighs 4
  tmp2 = 0.866025403784439*out[14] - 0.866025403784439*out[49]; // stacked, weighs 3
  out[14] = 1.09544511501033*out[56] - 0.273861278752583*out[7] - 0.612372435695794*out[42]; // final, weighs 5
  tmp3 = 1.09544511501033*out[37] - 0.273861278752583*out[23] - 0.612372435695794*out[2]; // stacked, weighs 5
  tmp4 = 0.866025403784439*out[16] - 0.866025403784439*out[51]; // stacked, weighs 3
  tmp5 = out[65] - 0.670820393249937*out[16] - 0.670820393249937*out[51]; // stacked, weighs 4
  out[16] = 1.09544511501033*out[58] - 0.273861278752583*out[9] - 0.612372435695794*out[44]; // final, weighs 5
  tmp6 = 1.09544511501033*out[39] - 0.273861278752583*out[25] - 0.612372435695794*out[4]; // stacked, weighs 5
  tmp7 = out[67] - 0.670820393249937*out[18] - 0.670820393249937*out[53]; // stacked, weighs 4
  tmp8 = 0.866025403784439*out[18] - 0.866025403784439*out[53]; // stacked, weighs 3
  out[18] = 1.09544511501033*out[60] - 0.273861278752583*out[11] - 0.612372435695794*out[46]; // final, weighs 5
  tmp9 = 1.09544511501033*out[41] - 0.273861278752583*out[27] - 0.612372435695794*out[6]; // stacked, weighs 5
  tmp10 = 0.866025403784439*out[20] - 0.866025403784439*out[55]; // stacked, weighs 3
  tmp11 = out[69] - 0.670820393249937*out[20] - 0.670820393249937*out[55]; // stacked, weighs 4
  out[20] = 1.09544511501033*out[62] - 0.273861278752583*out[13] - 0.612372435695794*out[48]; // final, weighs 5
  out[35] = 0.790569415042095*out[0] - 1.06066017177982*out[21]; // final, weighs 3
  out[21] = tmp2; // final, weighs 0
  out[0] = tmp1; // final, weighs 0
  tmp2 = 1.09544511501033*out[36] - 0.273861278752583*out[22] - 0.612372435695794*out[1]; // stacked, weighs 5
  out[36] = 0.790569415042095*out[1] - 1.06066017177982*out[22]; // final, weighs 3
  out[1] = out[64] - 0.670820393249937*out[15] - 0.670820393249937*out[50]; // final, weighs 4
  out[22] = 0.866025403784439*out[15] - 0.866025403784439*out[50]; // final, weighs 3
  out[15] = 1.09544511501033*out[57] - 0.273861278752583*out[8] - 0.612372435695794*out[43]; // final, weighs 5
  out[37] = 0.790569415042095*out[2] - 1.06066017177982*out[23]; // final, weighs 3
  out[2] = tmp5; // final, weighs 0
  out[23] = tmp4; // final, weighs 0
  tmp1 = 1.09544511501033*out[38] - 0.273861278752583*out[24] - 0.612372435695794*out[3]; // stacked, weighs 5
  out[38] = 0.790569415042095*out[3] - 1.06066017177982*out[24]; // final, weighs 3
  out[3] = out[66] - 0.670820393249937*out[17] - 0.670820393249937*out[52]; // final, weighs 4
  out[24] = 0.866025403784439*out[17] - 0.866025403784439*out[52]; // final, weighs 3
  out[17] = 1.09544511501033*out[59] - 0.273861278752583*out[10] - 0.612372435695794*out[45]; // final, weighs 5
  out[39] = 0.790569415042095*out[4] - 1.06066017177982*out[25]; // final, weighs 3
  out[25] = tmp8; // final, weighs 0
  out[4] = tmp7; // final, weighs 0
  tmp5 = 1.09544511501033*out[40] - 0.273861278752583*out[26] - 0.612372435695794*out[5]; // stacked, weighs 5
  out[40] = 0.790569415042095*out[5] - 1.06066017177982*out[26]; // final, weighs 3
  out[5] = out[68] - 0.670820393249937*out[19] - 0.670820393249937*out[54]; // final, weighs 4
  out[26] = 0.866025403784439*out[19] - 0.866025403784439*out[54]; // final, weighs 3
  out[19] = 1.09544511501033*out[61] - 0.273861278752583*out[12] - 0.612372435695794*out[47]; // final, weighs 5
  out[41] = 0.790569415042095*out[6] - 1.06066017177982*out[27]; // final, weighs 3
  out[6] = tmp11; // final, weighs 0
  out[27] = tmp10; // final, weighs 0
  out[42] = 1.06066017177982*out[7] - 0.790569415042095*out[42]; // final, weighs 3
  out[7] = tmp0; // final, weighs 0
  out[43] = 1.06066017177982*out[8] - 0.790569415042095*out[43]; // final, weighs 3
  out[8] = tmp2; // final, weighs 0
  out[44] = 1.06066017177982*out[9] - 0.790569415042095*out[44]; // final, weighs 3
  out[9] = tmp3; // final, weighs 0
  out[45] = 1.06066017177982*out[10] - 0.790569415042095*out[45]; // final, weighs 3
  out[10] = tmp1; // final, weighs 0
  out[46] = 1.06066017177982*out[11] - 0.790569415042095*out[46]; // final, weighs 3
  out[11] = tmp6; // final, weighs 0
  out[47] = 1.06066017177982*out[12] - 0.790569415042095*out[47]; // final, weighs 3
  out[12] = tmp5; // final, weighs 0
  out[48] = 1.06066017177982*out[13] - 0.790569415042095*out[48]; // final, weighs 3
  out[13] = tmp9; // final, weighs 0
  // total weight = 161
}

void gint2_nai_pD_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 8
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  gint2_nai_cD_pF(a, a_a, b, b_a, c, out);
  tmp0 = out[14]; // stacked, weighs 0
  out[14] = out[28]; // final, weighs 0
  tmp1 = out[15]; // stacked, weighs 0
  out[15] = out[29]; // final, weighs 0
  tmp2 = out[16]; // stacked, weighs 0
  out[16] = out[30]; // final, weighs 0
  tmp3 = out[17]; // stacked, weighs 0
  out[17] = out[31]; // final, weighs 0
  tmp4 = out[18]; // stacked, weighs 0
  out[18] = out[32]; // final, weighs 0
  tmp5 = out[19]; // stacked, weighs 0
  out[19] = out[33]; // final, weighs 0
  tmp6 = out[20]; // stacked, weighs 0
  out[20] = out[34]; // final, weighs 0
  tmp7 = out[35] - 0.5*out[0] - 0.5*out[21]; // stacked, weighs 4
  out[21] = 0.866025403784439*out[0] - 0.866025403784439*out[21]; // final, weighs 3
  out[0] = tmp7; // final, weighs 0
  tmp7 = out[36] - 0.5*out[1] - 0.5*out[22]; // stacked, weighs 4
  out[22] = 0.866025403784439*out[1] - 0.866025403784439*out[22]; // final, weighs 3
  out[1] = tmp7; // final, weighs 0
  tmp7 = out[37] - 0.5*out[2] - 0.5*out[23]; // stacked, weighs 4
  out[23] = 0.866025403784439*out[2] - 0.866025403784439*out[23]; // final, weighs 3
  out[2] = tmp7; // final, weighs 0
  tmp7 = out[38] - 0.5*out[24] - 0.5*out[3]; // stacked, weighs 4
  out[24] = 0.866025403784439*out[3] - 0.866025403784439*out[24]; // final, weighs 3
  out[3] = tmp7; // final, weighs 0
  tmp7 = out[39] - 0.5*out[25] - 0.5*out[4]; // stacked, weighs 4
  out[25] = 0.866025403784439*out[4] - 0.866025403784439*out[25]; // final, weighs 3
  out[4] = tmp7; // final, weighs 0
  tmp7 = out[40] - 0.5*out[26] - 0.5*out[5]; // stacked, weighs 4
  out[26] = 0.866025403784439*out[5] - 0.866025403784439*out[26]; // final, weighs 3
  out[5] = tmp7; // final, weighs 0
  tmp7 = out[41] - 0.5*out[27] - 0.5*out[6]; // stacked, weighs 4
  out[27] = 0.866025403784439*out[6] - 0.866025403784439*out[27]; // final, weighs 3
  out[6] = tmp7; // final, weighs 0
  out[28] = out[7]; // final, weighs 0
  out[7] = tmp0; // final, weighs 0
  out[29] = out[8]; // final, weighs 0
  out[8] = tmp1; // final, weighs 0
  out[30] = out[9]; // final, weighs 0
  out[9] = tmp2; // final, weighs 0
  out[31] = out[10]; // final, weighs 0
  out[10] = tmp3; // final, weighs 0
  out[32] = out[11]; // final, weighs 0
  out[11] = tmp4; // final, weighs 0
  out[33] = out[12]; // final, weighs 0
  out[12] = tmp5; // final, weighs 0
  out[34] = out[13]; // final, weighs 0
  out[13] = tmp6; // final, weighs 0
  // total weight = 49
}

void gint2_nai_SP_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 3
  double tmp0, tmp1, tmp2;
  gint2_nai_SP_cF(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[5] - 0.273861278752583*out[3] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = 0.866025403784439*out[2] - 0.866025403784439*out[7]; // stacked, weighs 3
  tmp2 = out[9] - 0.670820393249937*out[2] - 0.670820393249937*out[7]; // stacked, weighs 4
  out[2] = 1.09544511501033*out[8] - 0.273861278752583*out[1] - 0.612372435695794*out[6]; // final, weighs 5
  out[5] = 0.790569415042095*out[0] - 1.06066017177982*out[3]; // final, weighs 3
  out[0] = tmp2; // final, weighs 0
  out[3] = tmp1; // final, weighs 0
  out[6] = 1.06066017177982*out[1] - 0.790569415042095*out[6]; // final, weighs 3
  out[1] = tmp0; // final, weighs 0
  out[7] = out[19] - 0.670820393249937*out[12] - 0.670820393249937*out[17]; // final, weighs 4
  out[8] = 1.09544511501033*out[15] - 0.273861278752583*out[13] - 0.612372435695794*out[10]; // final, weighs 5
  out[9] = 1.09544511501033*out[18] - 0.273861278752583*out[11] - 0.612372435695794*out[16]; // final, weighs 5
  tmp2 = 0.866025403784439*out[12] - 0.866025403784439*out[17]; // stacked, weighs 3
  out[12] = 0.790569415042095*out[10] - 1.06066017177982*out[13]; // final, weighs 3
  out[10] = tmp2; // final, weighs 0
  out[13] = 1.06066017177982*out[11] - 0.790569415042095*out[16]; // final, weighs 3
  out[11] = out[14]; // final, weighs 0
  out[14] = out[29] - 0.670820393249937*out[22] - 0.670820393249937*out[27]; // final, weighs 4
  out[15] = 1.09544511501033*out[25] - 0.273861278752583*out[23] - 0.612372435695794*out[20]; // final, weighs 5
  out[16] = 1.09544511501033*out[28] - 0.273861278752583*out[21] - 0.612372435695794*out[26]; // final, weighs 5
  out[17] = 0.866025403784439*out[22] - 0.866025403784439*out[27]; // final, weighs 3
  out[18] = out[24]; // final, weighs 0
  out[19] = 0.790569415042095*out[20] - 1.06066017177982*out[23]; // final, weighs 3
  out[20] = 1.06066017177982*out[21] - 0.790569415042095*out[26]; // final, weighs 3
  out[21] = out[39] - 0.670820393249937*out[32] - 0.670820393249937*out[37]; // final, weighs 4
  out[22] = 1.09544511501033*out[35] - 0.273861278752583*out[33] - 0.612372435695794*out[30]; // final, weighs 5
  out[23] = 1.09544511501033*out[38] - 0.273861278752583*out[31] - 0.612372435695794*out[36]; // final, weighs 5
  out[24] = 0.866025403784439*out[32] - 0.866025403784439*out[37]; // final, weighs 3
  out[25] = out[34]; // final, weighs 0
  out[26] = 0.790569415042095*out[30] - 1.06066017177982*out[33]; // final, weighs 3
  out[27] = 1.06066017177982*out[31] - 0.790569415042095*out[36]; // final, weighs 3
  // total weight = 92
}

void gint2_nai_S_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 3
  double tmp0, tmp1, tmp2;
  gint2_nai_S_cF(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[5] - 0.273861278752583*out[3] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = out[9] - 0.670820393249937*out[2] - 0.670820393249937*out[7]; // stacked, weighs 4
  tmp2 = 0.866025403784439*out[2] - 0.866025403784439*out[7]; // stacked, weighs 3
  out[2] = 1.09544511501033*out[8] - 0.273861278752583*out[1] - 0.612372435695794*out[6]; // final, weighs 5
  out[5] = 0.790569415042095*out[0] - 1.06066017177982*out[3]; // final, weighs 3
  out[3] = tmp2; // final, weighs 0
  out[0] = tmp1; // final, weighs 0
  out[6] = 1.06066017177982*out[1] - 0.790569415042095*out[6]; // final, weighs 3
  out[1] = tmp0; // final, weighs 0
  // total weight = 23
}

void gint2_nai_P_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 3
  double tmp0, tmp1, tmp2;
  gint2_nai_P_cF(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[5] - 0.273861278752583*out[3] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = 0.866025403784439*out[2] - 0.866025403784439*out[7]; // stacked, weighs 3
  tmp2 = out[9] - 0.670820393249937*out[2] - 0.670820393249937*out[7]; // stacked, weighs 4
  out[2] = 1.09544511501033*out[8] - 0.273861278752583*out[1] - 0.612372435695794*out[6]; // final, weighs 5
  out[5] = 0.790569415042095*out[0] - 1.06066017177982*out[3]; // final, weighs 3
  out[0] = tmp2; // final, weighs 0
  out[3] = tmp1; // final, weighs 0
  out[6] = 1.06066017177982*out[1] - 0.790569415042095*out[6]; // final, weighs 3
  out[1] = tmp0; // final, weighs 0
  out[7] = out[19] - 0.670820393249937*out[12] - 0.670820393249937*out[17]; // final, weighs 4
  out[8] = 1.09544511501033*out[15] - 0.273861278752583*out[13] - 0.612372435695794*out[10]; // final, weighs 5
  out[9] = 1.09544511501033*out[18] - 0.273861278752583*out[11] - 0.612372435695794*out[16]; // final, weighs 5
  tmp2 = 0.866025403784439*out[12] - 0.866025403784439*out[17]; // stacked, weighs 3
  out[12] = 0.790569415042095*out[10] - 1.06066017177982*out[13]; // final, weighs 3
  out[10] = tmp2; // final, weighs 0
  out[13] = 1.06066017177982*out[11] - 0.790569415042095*out[16]; // final, weighs 3
  out[11] = out[14]; // final, weighs 0
  out[14] = out[29] - 0.670820393249937*out[22] - 0.670820393249937*out[27]; // final, weighs 4
  out[15] = 1.09544511501033*out[25] - 0.273861278752583*out[23] - 0.612372435695794*out[20]; // final, weighs 5
  out[16] = 1.09544511501033*out[28] - 0.273861278752583*out[21] - 0.612372435695794*out[26]; // final, weighs 5
  out[17] = 0.866025403784439*out[22] - 0.866025403784439*out[27]; // final, weighs 3
  out[18] = out[24]; // final, weighs 0
  out[19] = 0.790569415042095*out[20] - 1.06066017177982*out[23]; // final, weighs 3
  out[20] = 1.06066017177982*out[21] - 0.790569415042095*out[26]; // final, weighs 3
  // total weight = 69
}

void gint2_nai_cD_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 3
  double tmp0, tmp1, tmp2;
  gint2_nai_cD_cF(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[5] - 0.273861278752583*out[3] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = 0.866025403784439*out[2] - 0.866025403784439*out[7]; // stacked, weighs 3
  tmp2 = out[9] - 0.670820393249937*out[2] - 0.670820393249937*out[7]; // stacked, weighs 4
  out[2] = 1.09544511501033*out[8] - 0.273861278752583*out[1] - 0.612372435695794*out[6]; // final, weighs 5
  out[5] = 0.790569415042095*out[0] - 1.06066017177982*out[3]; // final, weighs 3
  out[0] = tmp2; // final, weighs 0
  out[3] = tmp1; // final, weighs 0
  out[6] = 1.06066017177982*out[1] - 0.790569415042095*out[6]; // final, weighs 3
  out[1] = tmp0; // final, weighs 0
  out[7] = out[19] - 0.670820393249937*out[12] - 0.670820393249937*out[17]; // final, weighs 4
  out[8] = 1.09544511501033*out[15] - 0.273861278752583*out[13] - 0.612372435695794*out[10]; // final, weighs 5
  out[9] = 1.09544511501033*out[18] - 0.273861278752583*out[11] - 0.612372435695794*out[16]; // final, weighs 5
  tmp2 = 0.866025403784439*out[12] - 0.866025403784439*out[17]; // stacked, weighs 3
  out[12] = 0.790569415042095*out[10] - 1.06066017177982*out[13]; // final, weighs 3
  out[10] = tmp2; // final, weighs 0
  out[13] = 1.06066017177982*out[11] - 0.790569415042095*out[16]; // final, weighs 3
  out[11] = out[14]; // final, weighs 0
  out[14] = out[29] - 0.670820393249937*out[22] - 0.670820393249937*out[27]; // final, weighs 4
  out[15] = 1.09544511501033*out[25] - 0.273861278752583*out[23] - 0.612372435695794*out[20]; // final, weighs 5
  out[16] = 1.09544511501033*out[28] - 0.273861278752583*out[21] - 0.612372435695794*out[26]; // final, weighs 5
  out[17] = 0.866025403784439*out[22] - 0.866025403784439*out[27]; // final, weighs 3
  out[18] = out[24]; // final, weighs 0
  out[19] = 0.790569415042095*out[20] - 1.06066017177982*out[23]; // final, weighs 3
  out[20] = 1.06066017177982*out[21] - 0.790569415042095*out[26]; // final, weighs 3
  out[21] = out[39] - 0.670820393249937*out[32] - 0.670820393249937*out[37]; // final, weighs 4
  out[22] = 1.09544511501033*out[35] - 0.273861278752583*out[33] - 0.612372435695794*out[30]; // final, weighs 5
  out[23] = 1.09544511501033*out[38] - 0.273861278752583*out[31] - 0.612372435695794*out[36]; // final, weighs 5
  out[24] = 0.866025403784439*out[32] - 0.866025403784439*out[37]; // final, weighs 3
  out[25] = out[34]; // final, weighs 0
  out[26] = 0.790569415042095*out[30] - 1.06066017177982*out[33]; // final, weighs 3
  out[27] = 1.06066017177982*out[31] - 0.790569415042095*out[36]; // final, weighs 3
  out[28] = out[49] - 0.670820393249937*out[42] - 0.670820393249937*out[47]; // final, weighs 4
  out[29] = 1.09544511501033*out[45] - 0.273861278752583*out[43] - 0.612372435695794*out[40]; // final, weighs 5
  out[30] = 1.09544511501033*out[48] - 0.273861278752583*out[41] - 0.612372435695794*out[46]; // final, weighs 5
  out[31] = 0.866025403784439*out[42] - 0.866025403784439*out[47]; // final, weighs 3
  out[32] = out[44]; // final, weighs 0
  out[33] = 0.790569415042095*out[40] - 1.06066017177982*out[43]; // final, weighs 3
  out[34] = 1.06066017177982*out[41] - 0.790569415042095*out[46]; // final, weighs 3
  out[35] = out[59] - 0.670820393249937*out[52] - 0.670820393249937*out[57]; // final, weighs 4
  out[36] = 1.09544511501033*out[55] - 0.273861278752583*out[53] - 0.612372435695794*out[50]; // final, weighs 5
  out[37] = 1.09544511501033*out[58] - 0.273861278752583*out[51] - 0.612372435695794*out[56]; // final, weighs 5
  out[38] = 0.866025403784439*out[52] - 0.866025403784439*out[57]; // final, weighs 3
  out[39] = out[54]; // final, weighs 0
  out[40] = 0.790569415042095*out[50] - 1.06066017177982*out[53]; // final, weighs 3
  out[41] = 1.06066017177982*out[51] - 0.790569415042095*out[56]; // final, weighs 3
  // total weight = 138
}

void gint2_nai_cF_pF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 3
  double tmp0, tmp1, tmp2;
  gint2_nai_cF_cF(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[5] - 0.273861278752583*out[3] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = 0.866025403784439*out[2] - 0.866025403784439*out[7]; // stacked, weighs 3
  tmp2 = out[9] - 0.670820393249937*out[2] - 0.670820393249937*out[7]; // stacked, weighs 4
  out[2] = 1.09544511501033*out[8] - 0.273861278752583*out[1] - 0.612372435695794*out[6]; // final, weighs 5
  out[5] = 0.790569415042095*out[0] - 1.06066017177982*out[3]; // final, weighs 3
  out[0] = tmp2; // final, weighs 0
  out[3] = tmp1; // final, weighs 0
  out[6] = 1.06066017177982*out[1] - 0.790569415042095*out[6]; // final, weighs 3
  out[1] = tmp0; // final, weighs 0
  out[7] = out[19] - 0.670820393249937*out[12] - 0.670820393249937*out[17]; // final, weighs 4
  out[8] = 1.09544511501033*out[15] - 0.273861278752583*out[13] - 0.612372435695794*out[10]; // final, weighs 5
  out[9] = 1.09544511501033*out[18] - 0.273861278752583*out[11] - 0.612372435695794*out[16]; // final, weighs 5
  tmp2 = 0.866025403784439*out[12] - 0.866025403784439*out[17]; // stacked, weighs 3
  out[12] = 0.790569415042095*out[10] - 1.06066017177982*out[13]; // final, weighs 3
  out[10] = tmp2; // final, weighs 0
  out[13] = 1.06066017177982*out[11] - 0.790569415042095*out[16]; // final, weighs 3
  out[11] = out[14]; // final, weighs 0
  out[14] = out[29] - 0.670820393249937*out[22] - 0.670820393249937*out[27]; // final, weighs 4
  out[15] = 1.09544511501033*out[25] - 0.273861278752583*out[23] - 0.612372435695794*out[20]; // final, weighs 5
  out[16] = 1.09544511501033*out[28] - 0.273861278752583*out[21] - 0.612372435695794*out[26]; // final, weighs 5
  out[17] = 0.866025403784439*out[22] - 0.866025403784439*out[27]; // final, weighs 3
  out[18] = out[24]; // final, weighs 0
  out[19] = 0.790569415042095*out[20] - 1.06066017177982*out[23]; // final, weighs 3
  out[20] = 1.06066017177982*out[21] - 0.790569415042095*out[26]; // final, weighs 3
  out[21] = out[39] - 0.670820393249937*out[32] - 0.670820393249937*out[37]; // final, weighs 4
  out[22] = 1.09544511501033*out[35] - 0.273861278752583*out[33] - 0.612372435695794*out[30]; // final, weighs 5
  out[23] = 1.09544511501033*out[38] - 0.273861278752583*out[31] - 0.612372435695794*out[36]; // final, weighs 5
  out[24] = 0.866025403784439*out[32] - 0.866025403784439*out[37]; // final, weighs 3
  out[25] = out[34]; // final, weighs 0
  out[26] = 0.790569415042095*out[30] - 1.06066017177982*out[33]; // final, weighs 3
  out[27] = 1.06066017177982*out[31] - 0.790569415042095*out[36]; // final, weighs 3
  out[28] = out[49] - 0.670820393249937*out[42] - 0.670820393249937*out[47]; // final, weighs 4
  out[29] = 1.09544511501033*out[45] - 0.273861278752583*out[43] - 0.612372435695794*out[40]; // final, weighs 5
  out[30] = 1.09544511501033*out[48] - 0.273861278752583*out[41] - 0.612372435695794*out[46]; // final, weighs 5
  out[31] = 0.866025403784439*out[42] - 0.866025403784439*out[47]; // final, weighs 3
  out[32] = out[44]; // final, weighs 0
  out[33] = 0.790569415042095*out[40] - 1.06066017177982*out[43]; // final, weighs 3
  out[34] = 1.06066017177982*out[41] - 0.790569415042095*out[46]; // final, weighs 3
  out[35] = out[59] - 0.670820393249937*out[52] - 0.670820393249937*out[57]; // final, weighs 4
  out[36] = 1.09544511501033*out[55] - 0.273861278752583*out[53] - 0.612372435695794*out[50]; // final, weighs 5
  out[37] = 1.09544511501033*out[58] - 0.273861278752583*out[51] - 0.612372435695794*out[56]; // final, weighs 5
  out[38] = 0.866025403784439*out[52] - 0.866025403784439*out[57]; // final, weighs 3
  out[39] = out[54]; // final, weighs 0
  out[40] = 0.790569415042095*out[50] - 1.06066017177982*out[53]; // final, weighs 3
  out[41] = 1.06066017177982*out[51] - 0.790569415042095*out[56]; // final, weighs 3
  out[42] = out[69] - 0.670820393249937*out[62] - 0.670820393249937*out[67]; // final, weighs 4
  out[43] = 1.09544511501033*out[65] - 0.273861278752583*out[63] - 0.612372435695794*out[60]; // final, weighs 5
  out[44] = 1.09544511501033*out[68] - 0.273861278752583*out[61] - 0.612372435695794*out[66]; // final, weighs 5
  out[45] = 0.866025403784439*out[62] - 0.866025403784439*out[67]; // final, weighs 3
  out[46] = out[64]; // final, weighs 0
  out[47] = 0.790569415042095*out[60] - 1.06066017177982*out[63]; // final, weighs 3
  out[48] = 1.06066017177982*out[61] - 0.790569415042095*out[66]; // final, weighs 3
  out[49] = out[79] - 0.670820393249937*out[72] - 0.670820393249937*out[77]; // final, weighs 4
  out[50] = 1.09544511501033*out[75] - 0.273861278752583*out[73] - 0.612372435695794*out[70]; // final, weighs 5
  out[51] = 1.09544511501033*out[78] - 0.273861278752583*out[71] - 0.612372435695794*out[76]; // final, weighs 5
  out[52] = 0.866025403784439*out[72] - 0.866025403784439*out[77]; // final, weighs 3
  out[53] = out[74]; // final, weighs 0
  out[54] = 0.790569415042095*out[70] - 1.06066017177982*out[73]; // final, weighs 3
  out[55] = 1.06066017177982*out[71] - 0.790569415042095*out[76]; // final, weighs 3
  out[56] = out[89] - 0.670820393249937*out[82] - 0.670820393249937*out[87]; // final, weighs 4
  out[57] = 1.09544511501033*out[85] - 0.273861278752583*out[83] - 0.612372435695794*out[80]; // final, weighs 5
  out[58] = 1.09544511501033*out[88] - 0.273861278752583*out[81] - 0.612372435695794*out[86]; // final, weighs 5
  out[59] = 0.866025403784439*out[82] - 0.866025403784439*out[87]; // final, weighs 3
  out[60] = out[84]; // final, weighs 0
  out[61] = 0.790569415042095*out[80] - 1.06066017177982*out[83]; // final, weighs 3
  out[62] = 1.06066017177982*out[81] - 0.790569415042095*out[86]; // final, weighs 3
  out[63] = out[99] - 0.670820393249937*out[92] - 0.670820393249937*out[97]; // final, weighs 4
  out[64] = 1.09544511501033*out[95] - 0.273861278752583*out[93] - 0.612372435695794*out[90]; // final, weighs 5
  out[65] = 1.09544511501033*out[98] - 0.273861278752583*out[91] - 0.612372435695794*out[96]; // final, weighs 5
  out[66] = 0.866025403784439*out[92] - 0.866025403784439*out[97]; // final, weighs 3
  out[67] = out[94]; // final, weighs 0
  out[68] = 0.790569415042095*out[90] - 1.06066017177982*out[93]; // final, weighs 3
  out[69] = 1.06066017177982*out[91] - 0.790569415042095*out[96]; // final, weighs 3
  // total weight = 230
}

void gint2_nai_pF_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 9
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  gint2_nai_cF_pD(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[25] - 0.273861278752583*out[15] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = out[45] - 0.670820393249937*out[10] - 0.670820393249937*out[35]; // stacked, weighs 4
  tmp2 = 0.866025403784439*out[10] - 0.866025403784439*out[35]; // stacked, weighs 3
  out[10] = 1.09544511501033*out[40] - 0.273861278752583*out[5] - 0.612372435695794*out[30]; // final, weighs 5
  tmp3 = 1.09544511501033*out[27] - 0.273861278752583*out[17] - 0.612372435695794*out[2]; // stacked, weighs 5
  tmp4 = 0.866025403784439*out[12] - 0.866025403784439*out[37]; // stacked, weighs 3
  tmp5 = out[47] - 0.670820393249937*out[12] - 0.670820393249937*out[37]; // stacked, weighs 4
  out[12] = 1.09544511501033*out[42] - 0.273861278752583*out[7] - 0.612372435695794*out[32]; // final, weighs 5
  tmp6 = 1.09544511501033*out[29] - 0.273861278752583*out[19] - 0.612372435695794*out[4]; // stacked, weighs 5
  tmp7 = out[49] - 0.670820393249937*out[14] - 0.670820393249937*out[39]; // stacked, weighs 4
  tmp8 = 0.866025403784439*out[14] - 0.866025403784439*out[39]; // stacked, weighs 3
  out[14] = 1.09544511501033*out[44] - 0.273861278752583*out[9] - 0.612372435695794*out[34]; // final, weighs 5
  out[25] = 0.790569415042095*out[0] - 1.06066017177982*out[15]; // final, weighs 3
  out[15] = tmp2; // final, weighs 0
  out[0] = tmp1; // final, weighs 0
  tmp2 = 1.09544511501033*out[26] - 0.273861278752583*out[16] - 0.612372435695794*out[1]; // stacked, weighs 5
  out[26] = 0.790569415042095*out[1] - 1.06066017177982*out[16]; // final, weighs 3
  out[1] = out[46] - 0.670820393249937*out[11] - 0.670820393249937*out[36]; // final, weighs 4
  out[16] = 0.866025403784439*out[11] - 0.866025403784439*out[36]; // final, weighs 3
  out[11] = 1.09544511501033*out[41] - 0.273861278752583*out[6] - 0.612372435695794*out[31]; // final, weighs 5
  out[27] = 0.790569415042095*out[2] - 1.06066017177982*out[17]; // final, weighs 3
  out[2] = tmp5; // final, weighs 0
  out[17] = tmp4; // final, weighs 0
  tmp1 = 1.09544511501033*out[28] - 0.273861278752583*out[18] - 0.612372435695794*out[3]; // stacked, weighs 5
  out[28] = 0.790569415042095*out[3] - 1.06066017177982*out[18]; // final, weighs 3
  out[3] = out[48] - 0.670820393249937*out[13] - 0.670820393249937*out[38]; // final, weighs 4
  out[18] = 0.866025403784439*out[13] - 0.866025403784439*out[38]; // final, weighs 3
  out[13] = 1.09544511501033*out[43] - 0.273861278752583*out[8] - 0.612372435695794*out[33]; // final, weighs 5
  out[29] = 0.790569415042095*out[4] - 1.06066017177982*out[19]; // final, weighs 3
  out[19] = tmp8; // final, weighs 0
  out[4] = tmp7; // final, weighs 0
  out[30] = 1.06066017177982*out[5] - 0.790569415042095*out[30]; // final, weighs 3
  out[5] = tmp0; // final, weighs 0
  out[31] = 1.06066017177982*out[6] - 0.790569415042095*out[31]; // final, weighs 3
  out[6] = tmp2; // final, weighs 0
  out[32] = 1.06066017177982*out[7] - 0.790569415042095*out[32]; // final, weighs 3
  out[7] = tmp3; // final, weighs 0
  out[33] = 1.06066017177982*out[8] - 0.790569415042095*out[33]; // final, weighs 3
  out[8] = tmp1; // final, weighs 0
  out[34] = 1.06066017177982*out[9] - 0.790569415042095*out[34]; // final, weighs 3
  out[9] = tmp6; // final, weighs 0
  // total weight = 115
}

void gint2_nai_pD_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 6
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5;
  gint2_nai_cD_pD(a, a_a, b, b_a, c, out);
  tmp0 = out[10]; // stacked, weighs 0
  out[10] = out[20]; // final, weighs 0
  tmp1 = out[11]; // stacked, weighs 0
  out[11] = out[21]; // final, weighs 0
  tmp2 = out[12]; // stacked, weighs 0
  out[12] = out[22]; // final, weighs 0
  tmp3 = out[13]; // stacked, weighs 0
  out[13] = out[23]; // final, weighs 0
  tmp4 = out[14]; // stacked, weighs 0
  out[14] = out[24]; // final, weighs 0
  tmp5 = out[25] - 0.5*out[0] - 0.5*out[15]; // stacked, weighs 4
  out[15] = 0.866025403784439*out[0] - 0.866025403784439*out[15]; // final, weighs 3
  out[0] = tmp5; // final, weighs 0
  tmp5 = out[26] - 0.5*out[1] - 0.5*out[16]; // stacked, weighs 4
  out[16] = 0.866025403784439*out[1] - 0.866025403784439*out[16]; // final, weighs 3
  out[1] = tmp5; // final, weighs 0
  tmp5 = out[27] - 0.5*out[17] - 0.5*out[2]; // stacked, weighs 4
  out[17] = 0.866025403784439*out[2] - 0.866025403784439*out[17]; // final, weighs 3
  out[2] = tmp5; // final, weighs 0
  tmp5 = out[28] - 0.5*out[18] - 0.5*out[3]; // stacked, weighs 4
  out[18] = 0.866025403784439*out[3] - 0.866025403784439*out[18]; // final, weighs 3
  out[3] = tmp5; // final, weighs 0
  tmp5 = out[29] - 0.5*out[19] - 0.5*out[4]; // stacked, weighs 4
  out[19] = 0.866025403784439*out[4] - 0.866025403784439*out[19]; // final, weighs 3
  out[4] = tmp5; // final, weighs 0
  out[20] = out[5]; // final, weighs 0
  out[5] = tmp0; // final, weighs 0
  out[21] = out[6]; // final, weighs 0
  out[6] = tmp1; // final, weighs 0
  out[22] = out[7]; // final, weighs 0
  out[7] = tmp2; // final, weighs 0
  out[23] = out[8]; // final, weighs 0
  out[8] = tmp3; // final, weighs 0
  out[24] = out[9]; // final, weighs 0
  out[9] = tmp4; // final, weighs 0
  // total weight = 35
}

void gint2_nai_SP_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 2
  double tmp0, tmp1;
  gint2_nai_SP_cD(a, a_a, b, b_a, c, out);
  tmp0 = out[2]; // stacked, weighs 0
  out[2] = out[4]; // final, weighs 0
  tmp1 = out[5] - 0.5*out[0] - 0.5*out[3]; // stacked, weighs 4
  out[3] = 0.866025403784439*out[0] - 0.866025403784439*out[3]; // final, weighs 3
  out[0] = tmp1; // final, weighs 0
  out[4] = out[1]; // final, weighs 0
  out[1] = tmp0; // final, weighs 0
  out[5] = out[11] - 0.5*out[6] - 0.5*out[9]; // final, weighs 4
  tmp1 = out[8]; // stacked, weighs 0
  out[8] = 0.866025403784439*out[6] - 0.866025403784439*out[9]; // final, weighs 3
  out[6] = tmp1; // final, weighs 0
  out[9] = out[7]; // final, weighs 0
  out[7] = out[10]; // final, weighs 0
  out[10] = out[17] - 0.5*out[12] - 0.5*out[15]; // final, weighs 4
  out[11] = out[14]; // final, weighs 0
  tmp0 = 0.866025403784439*out[12] - 0.866025403784439*out[15]; // stacked, weighs 3
  out[12] = out[16]; // final, weighs 0
  out[14] = out[13]; // final, weighs 0
  out[13] = tmp0; // final, weighs 0
  out[15] = out[23] - 0.5*out[18] - 0.5*out[21]; // final, weighs 4
  out[16] = out[20]; // final, weighs 0
  out[17] = out[22]; // final, weighs 0
  out[18] = 0.866025403784439*out[18] - 0.866025403784439*out[21]; // final, weighs 3
  // total weight = 28
}

void gint2_nai_S_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 2
  double tmp0, tmp1;
  gint2_nai_S_cD(a, a_a, b, b_a, c, out);
  tmp0 = out[2]; // stacked, weighs 0
  out[2] = out[4]; // final, weighs 0
  tmp1 = out[5] - 0.5*out[0] - 0.5*out[3]; // stacked, weighs 4
  out[3] = 0.866025403784439*out[0] - 0.866025403784439*out[3]; // final, weighs 3
  out[0] = tmp1; // final, weighs 0
  out[4] = out[1]; // final, weighs 0
  out[1] = tmp0; // final, weighs 0
  // total weight = 7
}

void gint2_nai_P_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 2
  double tmp0, tmp1;
  gint2_nai_P_cD(a, a_a, b, b_a, c, out);
  tmp0 = out[2]; // stacked, weighs 0
  out[2] = out[4]; // final, weighs 0
  tmp1 = out[5] - 0.5*out[0] - 0.5*out[3]; // stacked, weighs 4
  out[3] = 0.866025403784439*out[0] - 0.866025403784439*out[3]; // final, weighs 3
  out[0] = tmp1; // final, weighs 0
  out[4] = out[1]; // final, weighs 0
  out[1] = tmp0; // final, weighs 0
  out[5] = out[11] - 0.5*out[6] - 0.5*out[9]; // final, weighs 4
  tmp1 = out[8]; // stacked, weighs 0
  out[8] = 0.866025403784439*out[6] - 0.866025403784439*out[9]; // final, weighs 3
  out[6] = tmp1; // final, weighs 0
  out[9] = out[7]; // final, weighs 0
  out[7] = out[10]; // final, weighs 0
  out[10] = out[17] - 0.5*out[12] - 0.5*out[15]; // final, weighs 4
  out[11] = out[14]; // final, weighs 0
  tmp0 = 0.866025403784439*out[12] - 0.866025403784439*out[15]; // stacked, weighs 3
  out[12] = out[16]; // final, weighs 0
  out[14] = out[13]; // final, weighs 0
  out[13] = tmp0; // final, weighs 0
  // total weight = 21
}

void gint2_nai_cD_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 2
  double tmp0, tmp1;
  gint2_nai_cD_cD(a, a_a, b, b_a, c, out);
  tmp0 = out[2]; // stacked, weighs 0
  out[2] = out[4]; // final, weighs 0
  tmp1 = out[5] - 0.5*out[0] - 0.5*out[3]; // stacked, weighs 4
  out[3] = 0.866025403784439*out[0] - 0.866025403784439*out[3]; // final, weighs 3
  out[0] = tmp1; // final, weighs 0
  out[4] = out[1]; // final, weighs 0
  out[1] = tmp0; // final, weighs 0
  out[5] = out[11] - 0.5*out[6] - 0.5*out[9]; // final, weighs 4
  tmp1 = out[8]; // stacked, weighs 0
  out[8] = 0.866025403784439*out[6] - 0.866025403784439*out[9]; // final, weighs 3
  out[6] = tmp1; // final, weighs 0
  out[9] = out[7]; // final, weighs 0
  out[7] = out[10]; // final, weighs 0
  out[10] = out[17] - 0.5*out[12] - 0.5*out[15]; // final, weighs 4
  out[11] = out[14]; // final, weighs 0
  tmp0 = 0.866025403784439*out[12] - 0.866025403784439*out[15]; // stacked, weighs 3
  out[12] = out[16]; // final, weighs 0
  out[14] = out[13]; // final, weighs 0
  out[13] = tmp0; // final, weighs 0
  out[15] = out[23] - 0.5*out[18] - 0.5*out[21]; // final, weighs 4
  out[16] = out[20]; // final, weighs 0
  out[17] = out[22]; // final, weighs 0
  out[18] = 0.866025403784439*out[18] - 0.866025403784439*out[21]; // final, weighs 3
  out[20] = out[29] - 0.5*out[24] - 0.5*out[27]; // final, weighs 4
  out[21] = out[26]; // final, weighs 0
  out[22] = out[28]; // final, weighs 0
  out[23] = 0.866025403784439*out[24] - 0.866025403784439*out[27]; // final, weighs 3
  out[24] = out[25]; // final, weighs 0
  out[25] = out[35] - 0.5*out[30] - 0.5*out[33]; // final, weighs 4
  out[26] = out[32]; // final, weighs 0
  out[27] = out[34]; // final, weighs 0
  out[28] = 0.866025403784439*out[30] - 0.866025403784439*out[33]; // final, weighs 3
  out[29] = out[31]; // final, weighs 0
  // total weight = 42
}

void gint2_nai_cF_pD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 2
  double tmp0, tmp1;
  gint2_nai_cF_cD(a, a_a, b, b_a, c, out);
  tmp0 = out[2]; // stacked, weighs 0
  out[2] = out[4]; // final, weighs 0
  tmp1 = out[5] - 0.5*out[0] - 0.5*out[3]; // stacked, weighs 4
  out[3] = 0.866025403784439*out[0] - 0.866025403784439*out[3]; // final, weighs 3
  out[0] = tmp1; // final, weighs 0
  out[4] = out[1]; // final, weighs 0
  out[1] = tmp0; // final, weighs 0
  out[5] = out[11] - 0.5*out[6] - 0.5*out[9]; // final, weighs 4
  tmp1 = out[8]; // stacked, weighs 0
  out[8] = 0.866025403784439*out[6] - 0.866025403784439*out[9]; // final, weighs 3
  out[6] = tmp1; // final, weighs 0
  out[9] = out[7]; // final, weighs 0
  out[7] = out[10]; // final, weighs 0
  out[10] = out[17] - 0.5*out[12] - 0.5*out[15]; // final, weighs 4
  out[11] = out[14]; // final, weighs 0
  tmp0 = 0.866025403784439*out[12] - 0.866025403784439*out[15]; // stacked, weighs 3
  out[12] = out[16]; // final, weighs 0
  out[14] = out[13]; // final, weighs 0
  out[13] = tmp0; // final, weighs 0
  out[15] = out[23] - 0.5*out[18] - 0.5*out[21]; // final, weighs 4
  out[16] = out[20]; // final, weighs 0
  out[17] = out[22]; // final, weighs 0
  out[18] = 0.866025403784439*out[18] - 0.866025403784439*out[21]; // final, weighs 3
  out[20] = out[29] - 0.5*out[24] - 0.5*out[27]; // final, weighs 4
  out[21] = out[26]; // final, weighs 0
  out[22] = out[28]; // final, weighs 0
  out[23] = 0.866025403784439*out[24] - 0.866025403784439*out[27]; // final, weighs 3
  out[24] = out[25]; // final, weighs 0
  out[25] = out[35] - 0.5*out[30] - 0.5*out[33]; // final, weighs 4
  out[26] = out[32]; // final, weighs 0
  out[27] = out[34]; // final, weighs 0
  out[28] = 0.866025403784439*out[30] - 0.866025403784439*out[33]; // final, weighs 3
  out[29] = out[31]; // final, weighs 0
  out[30] = out[41] - 0.5*out[36] - 0.5*out[39]; // final, weighs 4
  out[31] = out[38]; // final, weighs 0
  out[32] = out[40]; // final, weighs 0
  out[33] = 0.866025403784439*out[36] - 0.866025403784439*out[39]; // final, weighs 3
  out[34] = out[37]; // final, weighs 0
  out[35] = out[47] - 0.5*out[42] - 0.5*out[45]; // final, weighs 4
  out[36] = out[44]; // final, weighs 0
  out[37] = out[46]; // final, weighs 0
  out[38] = 0.866025403784439*out[42] - 0.866025403784439*out[45]; // final, weighs 3
  out[39] = out[43]; // final, weighs 0
  out[40] = out[53] - 0.5*out[48] - 0.5*out[51]; // final, weighs 4
  out[41] = out[50]; // final, weighs 0
  out[42] = out[52]; // final, weighs 0
  out[43] = 0.866025403784439*out[48] - 0.866025403784439*out[51]; // final, weighs 3
  out[44] = out[49]; // final, weighs 0
  out[45] = out[59] - 0.5*out[54] - 0.5*out[57]; // final, weighs 4
  out[46] = out[56]; // final, weighs 0
  out[47] = out[58]; // final, weighs 0
  out[48] = 0.866025403784439*out[54] - 0.866025403784439*out[57]; // final, weighs 3
  out[49] = out[55]; // final, weighs 0
  // total weight = 70
}

void gint2_nai_pF_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 6
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5;
  gint2_nai_cF_SP(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[20] - 0.273861278752583*out[12] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = out[36] - 0.670820393249937*out[28] - 0.670820393249937*out[8]; // stacked, weighs 4
  tmp2 = 0.866025403784439*out[8] - 0.866025403784439*out[28]; // stacked, weighs 3
  out[8] = 1.09544511501033*out[32] - 0.273861278752583*out[4] - 0.612372435695794*out[24]; // final, weighs 5
  tmp3 = 1.09544511501033*out[22] - 0.273861278752583*out[14] - 0.612372435695794*out[2]; // stacked, weighs 5
  tmp4 = 0.866025403784439*out[10] - 0.866025403784439*out[30]; // stacked, weighs 3
  tmp5 = out[38] - 0.670820393249937*out[10] - 0.670820393249937*out[30]; // stacked, weighs 4
  out[10] = 1.09544511501033*out[34] - 0.273861278752583*out[6] - 0.612372435695794*out[26]; // final, weighs 5
  out[20] = 0.790569415042095*out[0] - 1.06066017177982*out[12]; // final, weighs 3
  out[12] = tmp2; // final, weighs 0
  out[0] = tmp1; // final, weighs 0
  tmp2 = 1.09544511501033*out[21] - 0.273861278752583*out[13] - 0.612372435695794*out[1]; // stacked, weighs 5
  out[21] = 0.790569415042095*out[1] - 1.06066017177982*out[13]; // final, weighs 3
  out[1] = out[37] - 0.670820393249937*out[29] - 0.670820393249937*out[9]; // final, weighs 4
  out[13] = 0.866025403784439*out[9] - 0.866025403784439*out[29]; // final, weighs 3
  out[9] = 1.09544511501033*out[33] - 0.273861278752583*out[5] - 0.612372435695794*out[25]; // final, weighs 5
  out[22] = 0.790569415042095*out[2] - 1.06066017177982*out[14]; // final, weighs 3
  out[2] = tmp5; // final, weighs 0
  out[14] = tmp4; // final, weighs 0
  tmp1 = 1.09544511501033*out[23] - 0.273861278752583*out[15] - 0.612372435695794*out[3]; // stacked, weighs 5
  out[23] = 0.790569415042095*out[3] - 1.06066017177982*out[15]; // final, weighs 3
  out[3] = out[39] - 0.670820393249937*out[11] - 0.670820393249937*out[31]; // final, weighs 4
  out[15] = 0.866025403784439*out[11] - 0.866025403784439*out[31]; // final, weighs 3
  out[11] = 1.09544511501033*out[35] - 0.273861278752583*out[7] - 0.612372435695794*out[27]; // final, weighs 5
  out[24] = 1.06066017177982*out[4] - 0.790569415042095*out[24]; // final, weighs 3
  out[4] = tmp0; // final, weighs 0
  out[25] = 1.06066017177982*out[5] - 0.790569415042095*out[25]; // final, weighs 3
  out[5] = tmp2; // final, weighs 0
  out[26] = 1.06066017177982*out[6] - 0.790569415042095*out[26]; // final, weighs 3
  out[6] = tmp3; // final, weighs 0
  out[27] = 1.06066017177982*out[7] - 0.790569415042095*out[27]; // final, weighs 3
  out[7] = tmp1; // final, weighs 0
  // total weight = 92
}

void gint2_nai_pD_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 5
  double tmp0, tmp1, tmp2, tmp3, tmp4;
  gint2_nai_cD_SP(a, a_a, b, b_a, c, out);
  tmp0 = out[8]; // stacked, weighs 0
  out[8] = out[16]; // final, weighs 0
  tmp1 = out[9]; // stacked, weighs 0
  out[9] = out[17]; // final, weighs 0
  tmp2 = out[10]; // stacked, weighs 0
  out[10] = out[18]; // final, weighs 0
  tmp3 = out[11]; // stacked, weighs 0
  out[11] = out[19]; // final, weighs 0
  tmp4 = out[20] - 0.5*out[0] - 0.5*out[12]; // stacked, weighs 4
  out[12] = 0.866025403784439*out[0] - 0.866025403784439*out[12]; // final, weighs 3
  out[0] = tmp4; // final, weighs 0
  tmp4 = out[21] - 0.5*out[1] - 0.5*out[13]; // stacked, weighs 4
  out[13] = 0.866025403784439*out[1] - 0.866025403784439*out[13]; // final, weighs 3
  out[1] = tmp4; // final, weighs 0
  tmp4 = out[22] - 0.5*out[14] - 0.5*out[2]; // stacked, weighs 4
  out[14] = 0.866025403784439*out[2] - 0.866025403784439*out[14]; // final, weighs 3
  out[2] = tmp4; // final, weighs 0
  tmp4 = out[23] - 0.5*out[15] - 0.5*out[3]; // stacked, weighs 4
  out[15] = 0.866025403784439*out[3] - 0.866025403784439*out[15]; // final, weighs 3
  out[3] = tmp4; // final, weighs 0
  out[16] = out[4]; // final, weighs 0
  out[4] = tmp0; // final, weighs 0
  out[17] = out[5]; // final, weighs 0
  out[5] = tmp1; // final, weighs 0
  out[18] = out[6]; // final, weighs 0
  out[6] = tmp2; // final, weighs 0
  out[19] = out[7]; // final, weighs 0
  out[7] = tmp3; // final, weighs 0
  // total weight = 28
}

void gint2_nai_SP_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 20
  double ab_a, d_0, d_1, d_2, nai_000_000_1, p_0, p_1, p_2, tmp1, tmp10, tmp11, tmp4, tmp9, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  d_2 = pow(b_a,0.75); // auto+recycle, weighs 2
  tmp4 = pow(a_a,0.75); // auto, weighs 2
  out[0] = 0.507949087473928*d_1*d_2*tmp4; // final, weighs 3
  nai_000_000_1 = 6.28318530717959*d_0*gaux(usq, 1); // local, weighs 3
  tmp11 = nai_000_000_1*u_0; // auto, weighs 1
  tmp1 = pow(b_a,1.25); // auto, weighs 2
  tmp4 = tmp1*tmp4; // auto+recycle, weighs 1
  out[1] = 1.01589817494786*tmp4*(-tmp11 + d_1*p_0); // final, weighs 5
  tmp10 = nai_000_000_1*u_1; // auto, weighs 1
  out[2] = 1.01589817494786*tmp4*(-tmp10 + d_1*p_1); // final, weighs 5
  tmp9 = nai_000_000_1*u_2; // auto, weighs 1
  out[3] = 1.01589817494786*tmp4*(-tmp9 + d_1*p_2); // final, weighs 5
  tmp11 = -tmp11 + d_1*v_0; // local+recycle, weighs 3
  tmp4 = pow(a_a,1.25); // auto+recycle, weighs 2
  d_2 = d_2*tmp4; // auto+recycle, weighs 1
  out[4] = 1.01589817494786*d_2*tmp11; // final, weighs 2
  usq = 6.28318530717959*d_0*gaux(usq, 2); // local+recycle, weighs 3
  v_0 = nai_000_000_1*v_0 - u_0*usq; // local+recycle, weighs 4
  ab_a = (0.5*d_1 - 0.5*nai_000_000_1)/ab_a; // auto+recycle, weighs 5
  d_0 = tmp1*tmp4; // auto+recycle, weighs 1
  out[5] = 2.03179634989571*d_0*(ab_a + p_0*tmp11 - u_0*v_0); // final, weighs 7
  out[6] = 2.03179634989571*d_0*(p_1*tmp11 - u_1*v_0); // final, weighs 6
  out[7] = 2.03179634989571*d_0*(p_2*tmp11 - u_2*v_0); // final, weighs 6
  tmp4 = -tmp10 + d_1*v_1; // local+recycle, weighs 3
  out[8] = 1.01589817494786*d_2*tmp4; // final, weighs 2
  tmp10 = nai_000_000_1*v_1 - u_1*usq; // local+recycle, weighs 4
  out[9] = 2.03179634989571*d_0*(p_0*tmp4 - tmp10*u_0); // final, weighs 6
  out[10] = 2.03179634989571*d_0*(ab_a + p_1*tmp4 - tmp10*u_1); // final, weighs 7
  out[11] = 2.03179634989571*d_0*(p_2*tmp4 - tmp10*u_2); // final, weighs 6
  tmp9 = -tmp9 + d_1*v_2; // local+recycle, weighs 3
  out[12] = 1.01589817494786*d_2*tmp9; // final, weighs 2
  nai_000_000_1 = nai_000_000_1*v_2 - u_2*usq; // local+recycle, weighs 4
  out[13] = 2.03179634989571*d_0*(p_0*tmp9 - nai_000_000_1*u_0); // final, weighs 6
  out[14] = 2.03179634989571*d_0*(p_1*tmp9 - nai_000_000_1*u_1); // final, weighs 6
  out[15] = 2.03179634989571*d_0*(ab_a + p_2*tmp9 - nai_000_000_1*u_2); // final, weighs 7
  // total weight = 189
}

void gint2_nai_S_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place permutation of the output of another routine.
  gint2_nai_SP_S(b, b_a, a, a_a, c, out);
}

void gint2_nai_P_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place permutation of the output of another routine.
  double tmp;
  gint2_nai_SP_P(b, b_a, a, a_a, c, out);
  tmp = out[1];
  out[1] = out[3];
  out[3] = out[9];
  out[9] = out[5];
  out[5] = out[4];
  out[4] = tmp;
  tmp = out[2];
  out[2] = out[6];
  out[6] = out[7];
  out[7] = out[10];
  out[10] = out[8];
  out[8] = tmp;
}

void gint2_nai_cD_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place permutation of the output of another routine.
  double tmp;
  gint2_nai_SP_cD(b, b_a, a, a_a, c, out);
  tmp = out[1];
  out[1] = out[6];
  out[6] = out[13];
  out[13] = out[9];
  out[9] = out[8];
  out[8] = out[2];
  out[2] = out[12];
  out[12] = out[3];
  out[3] = out[18];
  out[18] = out[16];
  out[16] = out[4];
  out[4] = tmp;
  tmp = out[5];
  out[5] = out[7];
  out[7] = out[19];
  out[19] = out[22];
  out[22] = out[17];
  out[17] = out[10];
  out[10] = out[14];
  out[14] = out[15];
  out[15] = out[21];
  out[21] = out[11];
  out[11] = out[20];
  out[20] = tmp;
}

void gint2_nai_cF_SP(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place permutation of the output of another routine.
  double tmp;
  gint2_nai_SP_cF(b, b_a, a, a_a, c, out);
  tmp = out[1];
  out[1] = out[10];
  out[10] = out[22];
  out[22] = out[25];
  out[25] = out[16];
  out[16] = out[4];
  out[4] = tmp;
  tmp = out[2];
  out[2] = out[20];
  out[20] = out[5];
  out[5] = out[11];
  out[11] = out[32];
  out[32] = out[8];
  out[8] = tmp;
  tmp = out[3];
  out[3] = out[30];
  out[30] = out[27];
  out[27] = out[36];
  out[36] = out[9];
  out[9] = out[12];
  out[12] = tmp;
  tmp = out[6];
  out[6] = out[21];
  out[21] = out[15];
  out[15] = out[33];
  out[33] = out[18];
  out[18] = out[24];
  out[24] = tmp;
  tmp = out[7];
  out[7] = out[31];
  out[31] = out[37];
  out[37] = out[19];
  out[19] = out[34];
  out[34] = out[28];
  out[28] = tmp;
  tmp = out[14];
  out[14] = out[23];
  out[23] = out[35];
  out[35] = out[38];
  out[38] = out[29];
  out[29] = out[17];
  out[17] = tmp;
}

void gint2_nai_pF_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 3
  double tmp0, tmp1, tmp2;
  gint2_nai_cF_S(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[5] - 0.273861278752583*out[3] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = out[9] - 0.670820393249937*out[2] - 0.670820393249937*out[7]; // stacked, weighs 4
  tmp2 = 0.866025403784439*out[2] - 0.866025403784439*out[7]; // stacked, weighs 3
  out[2] = 1.09544511501033*out[8] - 0.273861278752583*out[1] - 0.612372435695794*out[6]; // final, weighs 5
  out[5] = 0.790569415042095*out[0] - 1.06066017177982*out[3]; // final, weighs 3
  out[3] = tmp2; // final, weighs 0
  out[0] = tmp1; // final, weighs 0
  out[6] = 1.06066017177982*out[1] - 0.790569415042095*out[6]; // final, weighs 3
  out[1] = tmp0; // final, weighs 0
  // total weight = 23
}

void gint2_nai_pD_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 2
  double tmp0, tmp1;
  gint2_nai_cD_S(a, a_a, b, b_a, c, out);
  tmp0 = out[2]; // stacked, weighs 0
  out[2] = out[4]; // final, weighs 0
  tmp1 = out[5] - 0.5*out[0] - 0.5*out[3]; // stacked, weighs 4
  out[3] = 0.866025403784439*out[0] - 0.866025403784439*out[3]; // final, weighs 3
  out[0] = tmp1; // final, weighs 0
  out[4] = out[1]; // final, weighs 0
  out[1] = tmp0; // final, weighs 0
  // total weight = 7
}

void gint2_nai_SP_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 11
  double ab_a, d_0, d_1, d_2, p_0, p_1, p_2, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  ab_a = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_0 = 6.28318530717959*ab_a*gaux(usq, 0); // local+recycle, weighs 3
  d_1 = pow(b_a,0.75); // auto+recycle, weighs 2
  out[0] = 0.507949087473928*d_0*d_1*pow(a_a,0.75); // final, weighs 5
  d_2 = 6.28318530717959*ab_a*gaux(usq, 1); // local+recycle, weighs 3
  usq = d_1*pow(a_a,1.25); // auto+recycle, weighs 3
  out[1] = 1.01589817494786*usq*(d_0*(p_0 - a[0]) - d_2*u_0); // final, weighs 8
  out[2] = 1.01589817494786*usq*(d_0*(p_1 - a[1]) - d_2*u_1); // final, weighs 8
  out[3] = 1.01589817494786*usq*(d_0*(p_2 - a[2]) - d_2*u_2); // final, weighs 8
  // total weight = 87
}

void gint2_nai_S_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 7
  double ab_a, d_0, d_1, d_2, u_0, u_1, u_2;
  ab_a = a_a + b_a; // local, weighs 1
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = -c[0] + (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 7
  u_1 = -c[1] + (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 7
  u_2 = -c[2] + (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 7
  out[0] = 3.19153824321146*pow(a_a,0.75)*pow(b_a,0.75)*exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)*gaux(ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2), 0)/ab_a; // final, weighs 27
  // total weight = 55
}

void gint2_nai_P_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place permutation of the output of another routine.
  gint2_nai_S_P(b, b_a, a, a_a, c, out);
}

void gint2_nai_cD_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place permutation of the output of another routine.
  gint2_nai_S_cD(b, b_a, a, a_a, c, out);
}

void gint2_nai_cF_S(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place permutation of the output of another routine.
  gint2_nai_S_cF(b, b_a, a, a_a, c, out);
}

void gint2_nai_pF_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 6
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5;
  gint2_nai_cF_P(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[15] - 0.273861278752583*out[9] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = out[27] - 0.670820393249937*out[21] - 0.670820393249937*out[6]; // stacked, weighs 4
  tmp2 = 0.866025403784439*out[6] - 0.866025403784439*out[21]; // stacked, weighs 3
  out[6] = 1.09544511501033*out[24] - 0.273861278752583*out[3] - 0.612372435695794*out[18]; // final, weighs 5
  tmp3 = 1.09544511501033*out[17] - 0.273861278752583*out[11] - 0.612372435695794*out[2]; // stacked, weighs 5
  tmp4 = 0.866025403784439*out[8] - 0.866025403784439*out[23]; // stacked, weighs 3
  tmp5 = out[29] - 0.670820393249937*out[23] - 0.670820393249937*out[8]; // stacked, weighs 4
  out[8] = 1.09544511501033*out[26] - 0.273861278752583*out[5] - 0.612372435695794*out[20]; // final, weighs 5
  out[15] = 0.790569415042095*out[0] - 1.06066017177982*out[9]; // final, weighs 3
  out[9] = tmp2; // final, weighs 0
  out[0] = tmp1; // final, weighs 0
  tmp2 = 1.09544511501033*out[16] - 0.273861278752583*out[10] - 0.612372435695794*out[1]; // stacked, weighs 5
  out[16] = 0.790569415042095*out[1] - 1.06066017177982*out[10]; // final, weighs 3
  out[1] = out[28] - 0.670820393249937*out[22] - 0.670820393249937*out[7]; // final, weighs 4
  out[10] = 0.866025403784439*out[7] - 0.866025403784439*out[22]; // final, weighs 3
  out[7] = 1.09544511501033*out[25] - 0.273861278752583*out[4] - 0.612372435695794*out[19]; // final, weighs 5
  out[17] = 0.790569415042095*out[2] - 1.06066017177982*out[11]; // final, weighs 3
  out[2] = tmp5; // final, weighs 0
  out[11] = tmp4; // final, weighs 0
  out[18] = 1.06066017177982*out[3] - 0.790569415042095*out[18]; // final, weighs 3
  out[3] = tmp0; // final, weighs 0
  out[19] = 1.06066017177982*out[4] - 0.790569415042095*out[19]; // final, weighs 3
  out[4] = tmp2; // final, weighs 0
  out[20] = 1.06066017177982*out[5] - 0.790569415042095*out[20]; // final, weighs 3
  out[5] = tmp3; // final, weighs 0
  // total weight = 69
}

void gint2_nai_pD_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 4
  double tmp0, tmp1, tmp2, tmp3;
  gint2_nai_cD_P(a, a_a, b, b_a, c, out);
  tmp0 = out[6]; // stacked, weighs 0
  out[6] = out[12]; // final, weighs 0
  tmp1 = out[7]; // stacked, weighs 0
  out[7] = out[13]; // final, weighs 0
  tmp2 = out[8]; // stacked, weighs 0
  out[8] = out[14]; // final, weighs 0
  tmp3 = out[15] - 0.5*out[0] - 0.5*out[9]; // stacked, weighs 4
  out[9] = 0.866025403784439*out[0] - 0.866025403784439*out[9]; // final, weighs 3
  out[0] = tmp3; // final, weighs 0
  tmp3 = out[16] - 0.5*out[1] - 0.5*out[10]; // stacked, weighs 4
  out[10] = 0.866025403784439*out[1] - 0.866025403784439*out[10]; // final, weighs 3
  out[1] = tmp3; // final, weighs 0
  tmp3 = out[17] - 0.5*out[11] - 0.5*out[2]; // stacked, weighs 4
  out[11] = 0.866025403784439*out[2] - 0.866025403784439*out[11]; // final, weighs 3
  out[2] = tmp3; // final, weighs 0
  out[12] = out[3]; // final, weighs 0
  out[3] = tmp0; // final, weighs 0
  out[13] = out[4]; // final, weighs 0
  out[4] = tmp1; // final, weighs 0
  out[14] = out[5]; // final, weighs 0
  out[5] = tmp2; // final, weighs 0
  // total weight = 21
}

void gint2_nai_SP_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 19
  double ab_a, d_0, d_1, d_2, p_0, p_1, p_2, tmp1, tmp6, tmp7, tmp8, tmp9, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  tmp1 = pow(b_a,1.25); // auto, weighs 2
  tmp6 = tmp1*pow(a_a,0.75); // auto, weighs 3
  out[0] = 1.01589817494786*tmp6*(-tmp9 + d_1*p_0); // final, weighs 5
  tmp8 = d_2*u_1; // auto, weighs 1
  out[1] = 1.01589817494786*tmp6*(-tmp8 + d_1*p_1); // final, weighs 5
  tmp7 = d_2*u_2; // auto, weighs 1
  out[2] = 1.01589817494786*tmp6*(-tmp7 + d_1*p_2); // final, weighs 5
  tmp9 = -tmp9 + d_1*v_0; // local+recycle, weighs 3
  tmp6 = 6.28318530717959*d_0*gaux(usq, 2); // local+recycle, weighs 3
  usq = d_2*v_0 - tmp6*u_0; // local+recycle, weighs 4
  v_0 = (0.5*d_1 - 0.5*d_2)/ab_a; // auto+recycle, weighs 5
  ab_a = tmp1*pow(a_a,1.25); // auto+recycle, weighs 3
  out[3] = 2.03179634989571*ab_a*(v_0 + p_0*tmp9 - u_0*usq); // final, weighs 7
  out[4] = 2.03179634989571*ab_a*(p_1*tmp9 - u_1*usq); // final, weighs 6
  out[5] = 2.03179634989571*ab_a*(p_2*tmp9 - u_2*usq); // final, weighs 6
  d_0 = -tmp8 + d_1*v_1; // local+recycle, weighs 3
  tmp9 = d_2*v_1 - tmp6*u_1; // local+recycle, weighs 4
  out[6] = 2.03179634989571*ab_a*(d_0*p_0 - tmp9*u_0); // final, weighs 6
  out[7] = 2.03179634989571*ab_a*(v_0 + d_0*p_1 - tmp9*u_1); // final, weighs 7
  out[8] = 2.03179634989571*ab_a*(d_0*p_2 - tmp9*u_2); // final, weighs 6
  tmp1 = -tmp7 + d_1*v_2; // local+recycle, weighs 3
  tmp6 = d_2*v_2 - tmp6*u_2; // local+recycle, weighs 4
  out[9] = 2.03179634989571*ab_a*(p_0*tmp1 - tmp6*u_0); // final, weighs 6
  out[10] = 2.03179634989571*ab_a*(p_1*tmp1 - tmp6*u_1); // final, weighs 6
  out[11] = 2.03179634989571*ab_a*(v_0 + p_2*tmp1 - tmp6*u_2); // final, weighs 7
  // total weight = 177
}

void gint2_nai_S_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 11
  double ab_a, d_0, d_1, d_2, p_0, p_1, p_2, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
  d_0 = a[0] - b[0]; // local, weighs 2
  d_1 = a[1] - b[1]; // local, weighs 2
  d_2 = a[2] - b[2]; // local, weighs 2
  u_0 = p_0 - c[0]; // local, weighs 2
  u_1 = p_1 - c[1]; // local, weighs 2
  u_2 = p_2 - c[2]; // local, weighs 2
  usq = ab_a*(u_0*u_0 + u_1*u_1 + u_2*u_2); // local, weighs 6
  ab_a = exp(-a_a*b_a*(d_0*d_0 + d_1*d_1 + d_2*d_2)/ab_a)/ab_a; // auto+recycle, weighs 13
  d_0 = 6.28318530717959*ab_a*gaux(usq, 0); // local+recycle, weighs 3
  d_1 = 6.28318530717959*ab_a*gaux(usq, 1); // local+recycle, weighs 3
  d_2 = pow(a_a,0.75)*pow(b_a,1.25); // auto+recycle, weighs 5
  out[0] = 1.01589817494786*d_2*(d_0*(p_0 - b[0]) - d_1*u_0); // final, weighs 8
  out[1] = 1.01589817494786*d_2*(d_0*(p_1 - b[1]) - d_1*u_1); // final, weighs 8
  out[2] = 1.01589817494786*d_2*(d_0*(p_2 - b[2]) - d_1*u_2); // final, weighs 8
  // total weight = 82
}

void gint2_nai_P_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 15
  double ab_a, d_0, d_1, d_2, nai_100_000_0, p_0, p_1, p_2, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  usq = 6.28318530717959*d_0*gaux(usq, 2); // local+recycle, weighs 3
  v_0 = d_2*v_0 - u_0*usq; // local+recycle, weighs 4
  ab_a = (0.5*d_1 - 0.5*d_2)/ab_a; // auto+recycle, weighs 5
  d_0 = pow(a_a,1.25)*pow(b_a,1.25); // auto+recycle, weighs 5
  out[0] = 2.03179634989571*d_0*(ab_a + nai_100_000_0*p_0 - u_0*v_0); // final, weighs 7
  out[1] = 2.03179634989571*d_0*(nai_100_000_0*p_1 - u_1*v_0); // final, weighs 6
  out[2] = 2.03179634989571*d_0*(nai_100_000_0*p_2 - u_2*v_0); // final, weighs 6
  nai_100_000_0 = d_1*v_1 - d_2*u_1; // local+recycle, weighs 4
  v_0 = d_2*v_1 - u_1*usq; // local+recycle, weighs 4
  out[3] = 2.03179634989571*d_0*(nai_100_000_0*p_0 - u_0*v_0); // final, weighs 6
  out[4] = 2.03179634989571*d_0*(ab_a + nai_100_000_0*p_1 - u_1*v_0); // final, weighs 7
  out[5] = 2.03179634989571*d_0*(nai_100_000_0*p_2 - u_2*v_0); // final, weighs 6
  v_1 = d_1*v_2 - d_2*u_2; // local+recycle, weighs 4
  v_2 = d_2*v_2 - u_2*usq; // local+recycle, weighs 4
  out[6] = 2.03179634989571*d_0*(p_0*v_1 - u_0*v_2); // final, weighs 6
  out[7] = 2.03179634989571*d_0*(p_1*v_1 - u_1*v_2); // final, weighs 6
  out[8] = 2.03179634989571*d_0*(ab_a + p_2*v_1 - u_2*v_2); // final, weighs 7
  // total weight = 159
}

void gint2_nai_cD_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place permutation of the output of another routine.
  double tmp;
  gint2_nai_P_cD(b, b_a, a, a_a, c, out);
  tmp = out[1];
  out[1] = out[6];
  out[6] = out[2];
  out[2] = out[12];
  out[12] = out[4];
  out[4] = out[7];
  out[7] = out[8];
  out[8] = out[14];
  out[14] = out[16];
  out[16] = out[11];
  out[11] = out[15];
  out[15] = out[5];
  out[5] = out[13];
  out[13] = out[10];
  out[10] = out[9];
  out[9] = out[3];
  out[3] = tmp;
}

void gint2_nai_cF_P(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place permutation of the output of another routine.
  double tmp;
  gint2_nai_P_cF(b, b_a, a, a_a, c, out);
  tmp = out[1];
  out[1] = out[10];
  out[10] = out[13];
  out[13] = out[14];
  out[14] = out[24];
  out[24] = out[8];
  out[8] = out[22];
  out[22] = out[17];
  out[17] = out[25];
  out[25] = out[18];
  out[18] = out[6];
  out[6] = out[2];
  out[2] = out[20];
  out[20] = out[26];
  out[26] = out[28];
  out[28] = out[19];
  out[19] = out[16];
  out[16] = out[15];
  out[15] = out[5];
  out[5] = out[21];
  out[21] = out[7];
  out[7] = out[12];
  out[12] = out[4];
  out[4] = out[11];
  out[11] = out[23];
  out[23] = out[27];
  out[27] = out[9];
  out[9] = out[3];
  out[3] = tmp;
}

void gint2_nai_pF_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 9
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  gint2_nai_cF_cD(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[30] - 0.273861278752583*out[18] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = out[54] - 0.670820393249937*out[12] - 0.670820393249937*out[42]; // stacked, weighs 4
  tmp2 = 0.866025403784439*out[12] - 0.866025403784439*out[42]; // stacked, weighs 3
  out[12] = 1.09544511501033*out[48] - 0.273861278752583*out[6] - 0.612372435695794*out[36]; // final, weighs 5
  tmp3 = 1.09544511501033*out[32] - 0.273861278752583*out[20] - 0.612372435695794*out[2]; // stacked, weighs 5
  tmp4 = 0.866025403784439*out[14] - 0.866025403784439*out[44]; // stacked, weighs 3
  tmp5 = out[56] - 0.670820393249937*out[14] - 0.670820393249937*out[44]; // stacked, weighs 4
  out[14] = 1.09544511501033*out[50] - 0.273861278752583*out[8] - 0.612372435695794*out[38]; // final, weighs 5
  tmp6 = 1.09544511501033*out[34] - 0.273861278752583*out[22] - 0.612372435695794*out[4]; // stacked, weighs 5
  tmp7 = out[58] - 0.670820393249937*out[16] - 0.670820393249937*out[46]; // stacked, weighs 4
  tmp8 = 0.866025403784439*out[16] - 0.866025403784439*out[46]; // stacked, weighs 3
  out[16] = 1.09544511501033*out[52] - 0.273861278752583*out[10] - 0.612372435695794*out[40]; // final, weighs 5
  out[30] = 0.790569415042095*out[0] - 1.06066017177982*out[18]; // final, weighs 3
  out[18] = tmp2; // final, weighs 0
  out[0] = tmp1; // final, weighs 0
  tmp2 = 1.09544511501033*out[31] - 0.273861278752583*out[19] - 0.612372435695794*out[1]; // stacked, weighs 5
  out[31] = 0.790569415042095*out[1] - 1.06066017177982*out[19]; // final, weighs 3
  out[1] = out[55] - 0.670820393249937*out[13] - 0.670820393249937*out[43]; // final, weighs 4
  out[19] = 0.866025403784439*out[13] - 0.866025403784439*out[43]; // final, weighs 3
  out[13] = 1.09544511501033*out[49] - 0.273861278752583*out[7] - 0.612372435695794*out[37]; // final, weighs 5
  out[32] = 0.790569415042095*out[2] - 1.06066017177982*out[20]; // final, weighs 3
  out[2] = tmp5; // final, weighs 0
  out[20] = tmp4; // final, weighs 0
  tmp1 = 1.09544511501033*out[33] - 0.273861278752583*out[21] - 0.612372435695794*out[3]; // stacked, weighs 5
  out[33] = 0.790569415042095*out[3] - 1.06066017177982*out[21]; // final, weighs 3
  out[3] = out[57] - 0.670820393249937*out[15] - 0.670820393249937*out[45]; // final, weighs 4
  out[21] = 0.866025403784439*out[15] - 0.866025403784439*out[45]; // final, weighs 3
  out[15] = 1.09544511501033*out[51] - 0.273861278752583*out[9] - 0.612372435695794*out[39]; // final, weighs 5
  out[34] = 0.790569415042095*out[4] - 1.06066017177982*out[22]; // final, weighs 3
  out[22] = tmp8; // final, weighs 0
  out[4] = tmp7; // final, weighs 0
  tmp5 = 1.09544511501033*out[35] - 0.273861278752583*out[23] - 0.612372435695794*out[5]; // stacked, weighs 5
  out[35] = 0.790569415042095*out[5] - 1.06066017177982*out[23]; // final, weighs 3
  out[5] = out[59] - 0.670820393249937*out[17] - 0.670820393249937*out[47]; // final, weighs 4
  out[23] = 0.866025403784439*out[17] - 0.866025403784439*out[47]; // final, weighs 3
  out[17] = 1.09544511501033*out[53] - 0.273861278752583*out[11] - 0.612372435695794*out[41]; // final, weighs 5
  out[36] = 1.06066017177982*out[6] - 0.790569415042095*out[36]; // final, weighs 3
  out[6] = tmp0; // final, weighs 0
  out[37] = 1.06066017177982*out[7] - 0.790569415042095*out[37]; // final, weighs 3
  out[7] = tmp2; // final, weighs 0
  out[38] = 1.06066017177982*out[8] - 0.790569415042095*out[38]; // final, weighs 3
  out[8] = tmp3; // final, weighs 0
  out[39] = 1.06066017177982*out[9] - 0.790569415042095*out[39]; // final, weighs 3
  out[9] = tmp1; // final, weighs 0
  out[40] = 1.06066017177982*out[10] - 0.790569415042095*out[40]; // final, weighs 3
  out[10] = tmp6; // final, weighs 0
  out[41] = 1.06066017177982*out[11] - 0.790569415042095*out[41]; // final, weighs 3
  out[11] = tmp5; // final, weighs 0
  // total weight = 138
}

void gint2_nai_pD_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 7
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  gint2_nai_cD_cD(a, a_a, b, b_a, c, out);
  tmp0 = out[12]; // stacked, weighs 0
  out[12] = out[24]; // final, weighs 0
  tmp1 = out[13]; // stacked, weighs 0
  out[13] = out[25]; // final, weighs 0
  tmp2 = out[14]; // stacked, weighs 0
  out[14] = out[26]; // final, weighs 0
  tmp3 = out[15]; // stacked, weighs 0
  out[15] = out[27]; // final, weighs 0
  tmp4 = out[16]; // stacked, weighs 0
  out[16] = out[28]; // final, weighs 0
  tmp5 = out[17]; // stacked, weighs 0
  out[17] = out[29]; // final, weighs 0
  tmp6 = out[30] - 0.5*out[0] - 0.5*out[18]; // stacked, weighs 4
  out[18] = 0.866025403784439*out[0] - 0.866025403784439*out[18]; // final, weighs 3
  out[0] = tmp6; // final, weighs 0
  tmp6 = out[31] - 0.5*out[1] - 0.5*out[19]; // stacked, weighs 4
  out[19] = 0.866025403784439*out[1] - 0.866025403784439*out[19]; // final, weighs 3
  out[1] = tmp6; // final, weighs 0
  tmp6 = out[32] - 0.5*out[2] - 0.5*out[20]; // stacked, weighs 4
  out[20] = 0.866025403784439*out[2] - 0.866025403784439*out[20]; // final, weighs 3
  out[2] = tmp6; // final, weighs 0
  tmp6 = out[33] - 0.5*out[21] - 0.5*out[3]; // stacked, weighs 4
  out[21] = 0.866025403784439*out[3] - 0.866025403784439*out[21]; // final, weighs 3
  out[3] = tmp6; // final, weighs 0
  tmp6 = out[34] - 0.5*out[22] - 0.5*out[4]; // stacked, weighs 4
  out[22] = 0.866025403784439*out[4] - 0.866025403784439*out[22]; // final, weighs 3
  out[4] = tmp6; // final, weighs 0
  tmp6 = out[35] - 0.5*out[23] - 0.5*out[5]; // stacked, weighs 4
  out[23] = 0.866025403784439*out[5] - 0.866025403784439*out[23]; // final, weighs 3
  out[5] = tmp6; // final, weighs 0
  out[24] = out[6]; // final, weighs 0
  out[6] = tmp0; // final, weighs 0
  out[25] = out[7]; // final, weighs 0
  out[7] = tmp1; // final, weighs 0
  out[26] = out[8]; // final, weighs 0
  out[8] = tmp2; // final, weighs 0
  out[27] = out[9]; // final, weighs 0
  out[9] = tmp3; // final, weighs 0
  out[28] = out[10]; // final, weighs 0
  out[10] = tmp4; // final, weighs 0
  out[29] = out[11]; // final, weighs 0
  out[11] = tmp5; // final, weighs 0
  // total weight = 42
}

void gint2_nai_SP_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 34
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_001_0, nai_000_001_1, nai_000_002_0, nai_000_010_0, nai_000_010_1, nai_000_020_0, nai_000_100_0, nai_000_100_1, nai_000_200_0, nai_000_200_1, p_0, p_1, p_2, tmp1, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp3, tmp8, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  tmp3 = pow(b_a,1.75); // auto, weighs 2
  tmp8 = tmp3*pow(a_a,0.75); // auto, weighs 3
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
  usq = 6.28318530717959*d_0*gaux(usq, 3); // local+recycle, weighs 3
  d_0 = u_0*usq; // auto+recycle, weighs 1
  tmp8 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto+recycle, weighs 5
  nai_000_200_1 = tmp8 + nai_000_100_1*p_0 - u_0*(-d_0 + nai_000_000_2*p_0); // local, weighs 8
  tmp3 = tmp3*pow(a_a,1.25); // auto+recycle, weighs 3
  out[6] = 2.34611633910158*tmp3*(nai_000_200_0*v_0 + (nai_000_100_0 - nai_000_100_1)/ab_a - nai_000_200_1*u_0); // final, weighs 11
  tmp17 = -tmp17 + d_1*v_0; // local+recycle, weighs 3
  nai_000_100_0 = -tmp14 + d_2*v_0; // local+recycle, weighs 3
  tmp14 = -d_0 + nai_000_000_2*v_0; // local+recycle, weighs 3
  out[7] = 4.06359269979142*tmp3*(p_0*(p_1*tmp17 - nai_000_100_0*u_1) + (0.5*nai_000_010_0 - 0.5*nai_000_010_1)/ab_a - u_0*(nai_000_100_0*p_1 - tmp14*u_1)); // final, weighs 20
  nai_000_100_1 = p_2*tmp17 - nai_000_100_0*u_2; // local+recycle, weighs 4
  d_0 = nai_000_100_0*p_2 - tmp14*u_2; // local+recycle, weighs 4
  tmp17 = (0.5*nai_000_001_0 - 0.5*nai_000_001_1)/ab_a; // auto+recycle, weighs 5
  out[8] = 4.06359269979142*tmp3*(tmp17 + nai_000_100_1*p_0 - d_0*u_0); // final, weighs 7
  nai_000_100_0 = u_1*usq; // auto+recycle, weighs 1
  tmp14 = tmp8 + nai_000_010_1*p_1 - u_1*(-nai_000_100_0 + nai_000_000_2*p_1); // local+recycle, weighs 8
  out[9] = 2.34611633910158*tmp3*(nai_000_020_0*v_0 - tmp14*u_0); // final, weighs 6
  out[10] = 4.06359269979142*tmp3*(nai_000_100_1*p_1 - d_0*u_1); // final, weighs 6
  nai_000_100_1 = u_2*usq; // auto+recycle, weighs 1
  usq = tmp8 + nai_000_001_1*p_2 - u_2*(-nai_000_100_1 + nai_000_000_2*p_2); // local+recycle, weighs 8
  out[11] = 2.34611633910158*tmp3*(nai_000_002_0*v_0 - u_0*usq); // final, weighs 6
  out[12] = 2.34611633910158*tmp3*(nai_000_200_0*v_1 - nai_000_200_1*u_1); // final, weighs 6
  v_0 = -tmp16 + d_1*v_1; // local+recycle, weighs 3
  d_0 = -tmp13 + d_2*v_1; // local+recycle, weighs 3
  nai_000_100_0 = -nai_000_100_0 + nai_000_000_2*v_1; // local+recycle, weighs 3
  out[13] = 4.06359269979142*tmp3*(p_0*(tmp1 + p_1*v_0 - d_0*u_1) - u_0*(tmp8 + d_0*p_1 - nai_000_100_0*u_1)); // final, weighs 16
  tmp16 = p_2*v_0 - d_0*u_2; // local+recycle, weighs 4
  tmp13 = d_0*p_2 - nai_000_100_0*u_2; // local+recycle, weighs 4
  out[14] = 4.06359269979142*tmp3*(p_0*tmp16 - tmp13*u_0); // final, weighs 6
  out[15] = 2.34611633910158*tmp3*(nai_000_020_0*v_1 + (nai_000_010_0 - nai_000_010_1)/ab_a - tmp14*u_1); // final, weighs 11
  out[16] = 4.06359269979142*tmp3*(tmp17 + p_1*tmp16 - tmp13*u_1); // final, weighs 7
  out[17] = 2.34611633910158*tmp3*(nai_000_002_0*v_1 - u_1*usq); // final, weighs 6
  out[18] = 2.34611633910158*tmp3*(nai_000_200_0*v_2 - nai_000_200_1*u_2); // final, weighs 6
  nai_000_010_1 = -tmp15 + d_1*v_2; // local+recycle, weighs 3
  tmp15 = -tmp12 + d_2*v_2; // local+recycle, weighs 3
  tmp12 = -nai_000_100_1 + nai_000_000_2*v_2; // local+recycle, weighs 3
  out[19] = 4.06359269979142*tmp3*(p_0*(nai_000_010_1*p_1 - tmp15*u_1) - u_0*(p_1*tmp15 - tmp12*u_1)); // final, weighs 14
  v_0 = tmp1 + nai_000_010_1*p_2 - tmp15*u_2; // local+recycle, weighs 5
  v_1 = tmp8 + p_2*tmp15 - tmp12*u_2; // local+recycle, weighs 5
  out[20] = 4.06359269979142*tmp3*(p_0*v_0 - u_0*v_1); // final, weighs 6
  out[21] = 2.34611633910158*tmp3*(nai_000_020_0*v_2 - tmp14*u_2); // final, weighs 6
  out[22] = 4.06359269979142*tmp3*(p_1*v_0 - u_1*v_1); // final, weighs 6
  out[23] = 2.34611633910158*tmp3*(nai_000_002_0*v_2 + (nai_000_001_0 - nai_000_001_1)/ab_a - u_2*usq); // final, weighs 11
  // total weight = 394
}

void gint2_nai_S_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 13
  double ab_a, d_0, d_1, d_2, nai_000_010_0, nai_000_010_1, p_0, p_1, p_2, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  usq = 6.28318530717959*d_0*gaux(usq, 2); // local+recycle, weighs 3
  ab_a = (0.5*d_1 - 0.5*d_2)/ab_a; // auto+recycle, weighs 5
  d_0 = pow(a_a,0.75)*pow(b_a,1.75); // auto+recycle, weighs 5
  out[0] = 1.17305816955079*d_0*(ab_a + p_0*(d_1*p_0 - d_2*u_0) - u_0*(d_2*p_0 - u_0*usq)); // final, weighs 15
  nai_000_010_0 = d_1*p_1 - d_2*u_1; // local, weighs 4
  nai_000_010_1 = d_2*p_1 - u_1*usq; // local, weighs 4
  out[1] = 2.03179634989571*d_0*(nai_000_010_0*p_0 - nai_000_010_1*u_0); // final, weighs 6
  d_1 = d_1*p_2 - d_2*u_2; // local+recycle, weighs 4
  d_2 = d_2*p_2 - u_2*usq; // local+recycle, weighs 4
  out[2] = 2.03179634989571*d_0*(d_1*p_0 - d_2*u_0); // final, weighs 6
  out[3] = 1.17305816955079*d_0*(ab_a + nai_000_010_0*p_1 - nai_000_010_1*u_1); // final, weighs 7
  out[4] = 2.03179634989571*d_0*(d_1*p_1 - d_2*u_1); // final, weighs 6
  out[5] = 1.17305816955079*d_0*(ab_a + d_1*p_2 - d_2*u_2); // final, weighs 7
  // total weight = 135
}

void gint2_nai_P_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 33
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_001_0, nai_000_001_1, nai_000_010_1, nai_000_020_0, nai_000_020_1, nai_000_100_0, nai_000_100_1, nai_000_200_0, nai_000_200_1, p_0, p_1, p_2, tmp0, tmp1, tmp10, tmp11, tmp12, tmp15, tmp2, tmp6, tmp8, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  usq = 6.28318530717959*d_0*gaux(usq, 3); // local+recycle, weighs 3
  d_0 = u_0*usq; // auto+recycle, weighs 1
  tmp0 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto, weighs 5
  nai_000_200_1 = tmp0 + nai_000_100_1*p_0 - u_0*(-d_0 + nai_000_000_2*p_0); // local, weighs 8
  tmp6 = pow(a_a,1.25)*pow(b_a,1.75); // auto, weighs 5
  out[0] = 2.34611633910158*tmp6*(nai_000_200_0*v_0 + (nai_000_100_0 - nai_000_100_1)/ab_a - nai_000_200_1*u_0); // final, weighs 11
  tmp15 = -tmp15 + d_1*v_0; // local+recycle, weighs 3
  nai_000_100_0 = -tmp12 + d_2*v_0; // local+recycle, weighs 3
  tmp12 = -d_0 + nai_000_000_2*v_0; // local+recycle, weighs 3
  nai_000_100_1 = d_2*u_1; // auto+recycle, weighs 1
  d_0 = -nai_000_100_1 + d_1*p_1; // local+recycle, weighs 3
  tmp11 = nai_000_000_2*u_1; // auto, weighs 1
  nai_000_010_1 = -tmp11 + d_2*p_1; // local, weighs 3
  out[1] = 4.06359269979142*tmp6*(p_0*(p_1*tmp15 - nai_000_100_0*u_1) + (0.5*d_0 - 0.5*nai_000_010_1)/ab_a - u_0*(nai_000_100_0*p_1 - tmp12*u_1)); // final, weighs 20
  tmp15 = p_2*tmp15 - nai_000_100_0*u_2; // local+recycle, weighs 4
  nai_000_100_0 = nai_000_100_0*p_2 - tmp12*u_2; // local+recycle, weighs 4
  tmp12 = d_2*u_2; // auto+recycle, weighs 1
  nai_000_001_0 = -tmp12 + d_1*p_2; // local, weighs 3
  tmp10 = nai_000_000_2*u_2; // auto, weighs 1
  nai_000_001_1 = -tmp10 + d_2*p_2; // local, weighs 3
  tmp2 = (0.5*nai_000_001_0 - 0.5*nai_000_001_1)/ab_a; // auto, weighs 5
  out[2] = 4.06359269979142*tmp6*(tmp2 + p_0*tmp15 - nai_000_100_0*u_0); // final, weighs 7
  nai_000_020_0 = tmp1 + d_0*p_1 - nai_000_010_1*u_1; // local, weighs 5
  tmp8 = u_1*usq; // auto, weighs 1
  nai_000_020_1 = tmp0 + nai_000_010_1*p_1 - u_1*(-tmp8 + nai_000_000_2*p_1); // local, weighs 8
  out[3] = 2.34611633910158*tmp6*(nai_000_020_0*v_0 - nai_000_020_1*u_0); // final, weighs 6
  out[4] = 4.06359269979142*tmp6*(p_1*tmp15 - nai_000_100_0*u_1); // final, weighs 6
  tmp15 = tmp1 + nai_000_001_0*p_2 - nai_000_001_1*u_2; // local+recycle, weighs 5
  nai_000_100_0 = u_2*usq; // auto+recycle, weighs 1
  usq = tmp0 + nai_000_001_1*p_2 - u_2*(-nai_000_100_0 + nai_000_000_2*p_2); // local+recycle, weighs 8
  out[5] = 2.34611633910158*tmp6*(tmp15*v_0 - u_0*usq); // final, weighs 6
  out[6] = 2.34611633910158*tmp6*(nai_000_200_0*v_1 - nai_000_200_1*u_1); // final, weighs 6
  v_0 = -nai_000_100_1 + d_1*v_1; // local+recycle, weighs 3
  tmp11 = -tmp11 + d_2*v_1; // local+recycle, weighs 3
  tmp8 = -tmp8 + nai_000_000_2*v_1; // local+recycle, weighs 3
  out[7] = 4.06359269979142*tmp6*(p_0*(tmp1 + p_1*v_0 - tmp11*u_1) - u_0*(tmp0 + p_1*tmp11 - tmp8*u_1)); // final, weighs 16
  nai_000_100_1 = p_2*v_0 - tmp11*u_2; // local+recycle, weighs 4
  v_0 = p_2*tmp11 - tmp8*u_2; // local+recycle, weighs 4
  out[8] = 4.06359269979142*tmp6*(nai_000_100_1*p_0 - u_0*v_0); // final, weighs 6
  out[9] = 2.34611633910158*tmp6*(nai_000_020_0*v_1 + (d_0 - nai_000_010_1)/ab_a - nai_000_020_1*u_1); // final, weighs 11
  out[10] = 4.06359269979142*tmp6*(tmp2 + nai_000_100_1*p_1 - u_1*v_0); // final, weighs 7
  out[11] = 2.34611633910158*tmp6*(tmp15*v_1 - u_1*usq); // final, weighs 6
  out[12] = 2.34611633910158*tmp6*(nai_000_200_0*v_2 - nai_000_200_1*u_2); // final, weighs 6
  v_1 = -tmp12 + d_1*v_2; // local+recycle, weighs 3
  d_0 = -tmp10 + d_2*v_2; // local+recycle, weighs 3
  tmp11 = -nai_000_100_0 + nai_000_000_2*v_2; // local+recycle, weighs 3
  out[13] = 4.06359269979142*tmp6*(p_0*(p_1*v_1 - d_0*u_1) - u_0*(d_0*p_1 - tmp11*u_1)); // final, weighs 14
  nai_000_010_1 = tmp1 + p_2*v_1 - d_0*u_2; // local+recycle, weighs 5
  tmp0 = tmp0 + d_0*p_2 - tmp11*u_2; // local+recycle, weighs 5
  out[14] = 4.06359269979142*tmp6*(nai_000_010_1*p_0 - tmp0*u_0); // final, weighs 6
  out[15] = 2.34611633910158*tmp6*(nai_000_020_0*v_2 - nai_000_020_1*u_2); // final, weighs 6
  out[16] = 4.06359269979142*tmp6*(nai_000_010_1*p_1 - tmp0*u_1); // final, weighs 6
  out[17] = 2.34611633910158*tmp6*(tmp15*v_2 + (nai_000_001_0 - nai_000_001_1)/ab_a - u_2*usq); // final, weighs 11
  // total weight = 367
}

void gint2_nai_cD_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 48
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_000_002_1, nai_000_002_2, nai_000_020_0, nai_000_020_1, nai_001_010_1, nai_010_000_2, nai_010_000_3, nai_010_001_1, nai_010_010_0, nai_010_010_1, nai_100_000_0, nai_100_000_1, nai_100_000_2, nai_100_001_1, nai_110_000_2, nai_200_000_0, nai_200_000_1, nai_200_000_2, p_0, p_1, p_2, tmp0, tmp1, tmp11, tmp14, tmp17, tmp19, tmp2, tmp21, tmp22, tmp24, tmp29, tmp33, tmp36, tmp39, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  usq = 6.28318530717959*d_0*gaux(usq, 4); // local+recycle, weighs 3
  d_0 = u_0*usq; // auto+recycle, weighs 1
  tmp2 = (0.5*nai_000_000_2 - 0.5*nai_000_000_3)/ab_a; // auto, weighs 5
  nai_200_000_2 = tmp2 + nai_100_000_2*v_0 - u_0*(-d_0 + nai_000_000_3*v_0); // local, weighs 8
  tmp14 = 0.5*nai_200_000_0 - 0.5*nai_200_000_1; // auto, weighs 3
  tmp21 = pow(a_a,1.75)*pow(b_a,1.75); // auto, weighs 5
  out[0] = 2.70906179986095*tmp21*(p_0*(nai_200_000_0*p_0 + (nai_100_000_0 - nai_100_000_1)/ab_a - nai_200_000_1*u_0) + (tmp14 + tmp17 - tmp11 + nai_100_000_0*p_0 - nai_100_000_1*p_0)/ab_a - u_0*(nai_200_000_1*p_0 + (nai_100_000_1 - nai_100_000_2)/ab_a - nai_200_000_2*u_0)); // final, weighs 35
  tmp17 = nai_200_000_0*p_1 - nai_200_000_1*u_1; // local+recycle, weighs 4
  tmp11 = nai_200_000_1*p_1 - nai_200_000_2*u_1; // local+recycle, weighs 4
  out[1] = 4.69223267820315*tmp21*(p_0*tmp17 + (nai_100_000_0*p_1 + nai_100_000_2*u_1 - nai_100_000_1*p_1 - nai_100_000_1*u_1)/ab_a - tmp11*u_0); // final, weighs 18
  nai_200_000_0 = nai_200_000_0*p_2 - nai_200_000_1*u_2; // local+recycle, weighs 4
  nai_200_000_1 = nai_200_000_1*p_2 - nai_200_000_2*u_2; // local+recycle, weighs 4
  nai_200_000_2 = nai_100_000_0*p_2 - nai_100_000_1*u_2; // local+recycle, weighs 4
  nai_100_001_1 = nai_100_000_1*p_2 - nai_100_000_2*u_2; // local, weighs 4
  out[2] = 4.69223267820315*tmp21*(nai_200_000_0*p_0 + (nai_200_000_2 - nai_100_001_1)/ab_a - nai_200_000_1*u_0); // final, weighs 11
  tmp14 = tmp14/ab_a; // auto+recycle, weighs 2
  out[3] = 2.70906179986095*tmp21*(tmp14 + p_1*tmp17 - tmp11*u_1); // final, weighs 7
  out[4] = 4.69223267820315*tmp21*(nai_200_000_0*p_1 - nai_200_000_1*u_1); // final, weighs 6
  out[5] = 2.70906179986095*tmp21*(tmp14 + nai_200_000_0*p_2 - nai_200_000_1*u_2); // final, weighs 7
  tmp36 = -tmp36 + d_2*p_0; // local+recycle, weighs 3
  tmp17 = tmp1 + p_0*(-tmp39 + d_1*p_0) - tmp36*u_0; // local+recycle, weighs 8
  nai_200_000_0 = -tmp33 + nai_000_000_2*p_0; // local+recycle, weighs 3
  tmp33 = tmp0 + p_0*tmp36 - nai_200_000_0*u_0; // local+recycle, weighs 5
  tmp11 = tmp2 + nai_200_000_0*p_0 - u_0*(-d_0 + nai_000_000_3*p_0); // local+recycle, weighs 8
  nai_200_000_1 = d_2*u_1; // auto+recycle, weighs 1
  d_0 = -nai_200_000_1 + d_1*v_1; // local+recycle, weighs 3
  tmp14 = nai_000_000_2*u_1; // auto+recycle, weighs 1
  tmp39 = -tmp14 + d_2*v_1; // local+recycle, weighs 3
  tmp36 = tmp39*u_0; // auto+recycle, weighs 1
  nai_200_000_0 = nai_000_000_3*u_1; // auto+recycle, weighs 1
  nai_010_000_2 = -nai_200_000_0 + nai_000_000_2*v_1; // local, weighs 3
  tmp22 = nai_010_000_2*u_0; // auto, weighs 1
  out[6] = 4.69223267820315*tmp21*(v_0*(tmp17*v_1 - tmp33*u_1) + (tmp22 - tmp36 + d_0*p_0 - p_0*tmp39)/ab_a - u_0*(tmp33*v_1 - tmp11*u_1)); // final, weighs 24
  tmp36 = -tmp36 + d_0*v_0; // local+recycle, weighs 3
  tmp22 = -tmp22 + tmp39*v_0; // local+recycle, weighs 3
  nai_100_000_0 = (0.5*nai_100_000_0 - 0.5*nai_100_000_1)/ab_a; // auto+recycle, weighs 5
  tmp29 = u_1*usq; // auto, weighs 1
  nai_010_000_3 = -tmp29 + nai_000_000_3*v_1; // local, weighs 3
  nai_110_000_2 = nai_010_000_2*v_0 - nai_010_000_3*u_0; // local, weighs 4
  nai_100_000_1 = (0.5*nai_100_000_1 - 0.5*nai_100_000_2)/ab_a; // auto+recycle, weighs 5
  nai_100_000_2 = tmp1 - tmp39*u_1; // auto+recycle, weighs 3
  nai_010_010_0 = nai_100_000_2 + d_0*p_1; // local, weighs 2
  tmp19 = tmp0 - nai_010_000_2*u_1; // auto, weighs 3
  nai_010_010_1 = tmp19 + p_1*tmp39; // local, weighs 2
  out[7] = 8.12718539958285*tmp21*(p_0*(nai_100_000_0 + p_1*tmp36 - tmp22*u_1) + (0.5*nai_010_010_0 - 0.5*nai_010_010_1)/ab_a - u_0*(nai_100_000_1 + p_1*tmp22 - nai_110_000_2*u_1)); // final, weighs 22
  tmp36 = p_2*tmp36 - tmp22*u_2; // local+recycle, weighs 4
  tmp22 = p_2*tmp22 - nai_110_000_2*u_2; // local+recycle, weighs 4
  nai_110_000_2 = d_0*p_2 - tmp39*u_2; // local+recycle, weighs 4
  nai_010_001_1 = p_2*tmp39 - nai_010_000_2*u_2; // local, weighs 4
  out[8] = 8.12718539958285*tmp21*(p_0*tmp36 + (0.5*nai_110_000_2 - 0.5*nai_010_001_1)/ab_a - tmp22*u_0); // final, weighs 12
  nai_200_000_1 = -nai_200_000_1 + d_1*p_1; // local+recycle, weighs 3
  tmp14 = -tmp14 + d_2*p_1; // local+recycle, weighs 3
  nai_000_020_0 = tmp1 + nai_200_000_1*p_1 - tmp14*u_1; // local, weighs 5
  nai_200_000_0 = -nai_200_000_0 + nai_000_000_2*p_1; // local+recycle, weighs 3
  nai_000_020_1 = tmp0 + p_1*tmp14 - nai_200_000_0*u_1; // local, weighs 5
  tmp29 = tmp2 + nai_200_000_0*p_1 - u_1*(-tmp29 + nai_000_000_3*p_1); // local+recycle, weighs 8
  out[9] = 4.69223267820315*tmp21*(v_0*(nai_000_020_0*v_1 + (nai_200_000_1 - tmp14)/ab_a - nai_000_020_1*u_1) - u_0*(nai_000_020_1*v_1 + (tmp14 - nai_200_000_0)/ab_a - tmp29*u_1)); // final, weighs 24
  out[10] = 8.12718539958285*tmp21*(p_1*tmp36 + (0.5*nai_200_000_2 - 0.5*nai_100_001_1)/ab_a - tmp22*u_1); // final, weighs 12
  tmp36 = d_2*u_2; // auto+recycle, weighs 1
  nai_200_000_0 = -tmp36 + d_1*p_2; // local+recycle, weighs 3
  tmp22 = nai_000_000_2*u_2; // auto+recycle, weighs 1
  nai_200_000_1 = -tmp22 + d_2*p_2; // local+recycle, weighs 3
  nai_200_000_2 = tmp1 + nai_200_000_0*p_2 - nai_200_000_1*u_2; // local+recycle, weighs 5
  nai_100_001_1 = nai_000_000_3*u_2; // auto+recycle, weighs 1
  tmp14 = -nai_100_001_1 + nai_000_000_2*p_2; // local+recycle, weighs 3
  nai_000_002_1 = tmp0 + nai_200_000_1*p_2 - tmp14*u_2; // local, weighs 5
  usq = u_2*usq; // auto+recycle, weighs 1
  nai_000_002_2 = tmp2 + p_2*tmp14 - u_2*(-usq + nai_000_000_3*p_2); // local, weighs 8
  out[11] = 4.69223267820315*tmp21*(v_0*(nai_200_000_2*v_1 - nai_000_002_1*u_1) - u_0*(nai_000_002_1*v_1 - nai_000_002_2*u_1)); // final, weighs 14
  tmp17 = tmp17*v_2 - tmp33*u_2; // local+recycle, weighs 4
  tmp33 = tmp33*v_2 - tmp11*u_2; // local+recycle, weighs 4
  tmp11 = -tmp36 + d_1*v_2; // local+recycle, weighs 3
  d_1 = -tmp22 + d_2*v_2; // local+recycle, weighs 3
  d_2 = d_1*u_0; // auto+recycle, weighs 1
  nai_000_000_2 = -nai_100_001_1 + nai_000_000_2*v_2; // local+recycle, weighs 3
  tmp36 = nai_000_000_2*u_0; // auto+recycle, weighs 1
  out[12] = 4.69223267820315*tmp21*(tmp17*v_0 + (tmp36 - d_2 + p_0*tmp11 - d_1*p_0)/ab_a - tmp33*u_0); // final, weighs 16
  tmp22 = -d_2 + tmp11*v_0; // local+recycle, weighs 3
  nai_100_001_1 = -tmp36 + d_1*v_0; // local+recycle, weighs 3
  d_2 = -usq + nai_000_000_3*v_2; // local+recycle, weighs 3
  tmp36 = nai_000_000_2*v_0 - d_2*u_0; // local+recycle, weighs 4
  nai_000_000_3 = d_1*u_1; // auto+recycle, weighs 1
  usq = -nai_000_000_3 + p_1*tmp11; // local+recycle, weighs 3
  tmp24 = nai_000_000_2*u_1; // auto, weighs 1
  nai_001_010_1 = -tmp24 + d_1*p_1; // local, weighs 3
  out[13] = 8.12718539958285*tmp21*(p_0*(p_1*tmp22 - nai_100_001_1*u_1) + (0.5*usq - 0.5*nai_001_010_1)/ab_a - u_0*(nai_100_001_1*p_1 - tmp36*u_1)); // final, weighs 20
  nai_100_000_0 = nai_100_000_0 + p_2*tmp22 - nai_100_001_1*u_2; // local+recycle, weighs 5
  tmp36 = nai_100_000_1 + nai_100_001_1*p_2 - tmp36*u_2; // local+recycle, weighs 5
  nai_100_000_1 = tmp1 - d_1*u_2; // auto+recycle, weighs 3
  tmp1 = nai_100_000_1 + p_2*tmp11; // local+recycle, weighs 2
  tmp22 = tmp0 - nai_000_000_2*u_2; // auto+recycle, weighs 3
  tmp0 = tmp22 + d_1*p_2; // local+recycle, weighs 2
  nai_100_001_1 = (0.5*tmp1 - 0.5*tmp0)/ab_a; // auto+recycle, weighs 5
  out[14] = 8.12718539958285*tmp21*(nai_100_001_1 + nai_100_000_0*p_0 - tmp36*u_0); // final, weighs 7
  nai_000_020_0 = nai_000_020_0*v_2 - nai_000_020_1*u_2; // local+recycle, weighs 4
  tmp29 = nai_000_020_1*v_2 - tmp29*u_2; // local+recycle, weighs 4
  out[15] = 4.69223267820315*tmp21*(nai_000_020_0*v_0 - tmp29*u_0); // final, weighs 6
  out[16] = 8.12718539958285*tmp21*(nai_100_000_0*p_1 - tmp36*u_1); // final, weighs 6
  tmp36 = nai_200_000_2*v_2 + (nai_200_000_0 - nai_200_000_1)/ab_a - nai_000_002_1*u_2; // local+recycle, weighs 9
  nai_200_000_0 = nai_000_002_1*v_2 + (nai_200_000_1 - tmp14)/ab_a - nai_000_002_2*u_2; // local+recycle, weighs 9
  out[17] = 4.69223267820315*tmp21*(tmp36*v_0 - nai_200_000_0*u_0); // final, weighs 6
  nai_000_020_1 = nai_100_000_2 + d_0*v_1; // local+recycle, weighs 2
  nai_100_000_2 = tmp19 + tmp39*v_1; // local+recycle, weighs 2
  tmp19 = tmp2 + nai_010_000_2*v_1 - nai_010_000_3*u_1; // local+recycle, weighs 5
  nai_200_000_1 = 0.5*nai_000_020_1 - 0.5*nai_100_000_2; // auto+recycle, weighs 3
  nai_000_002_2 = nai_200_000_1/ab_a; // auto+recycle, weighs 2
  out[18] = 2.70906179986095*tmp21*(nai_000_002_2 + p_0*(nai_000_020_1*p_0 - nai_100_000_2*u_0) - u_0*(nai_100_000_2*p_0 - tmp19*u_0)); // final, weighs 15
  v_0 = nai_000_020_1*p_1 + (d_0 - tmp39)/ab_a - nai_100_000_2*u_1; // local+recycle, weighs 9
  nai_200_000_2 = nai_100_000_2*p_1 + (tmp39 - nai_010_000_2)/ab_a - tmp19*u_1; // local+recycle, weighs 9
  out[19] = 4.69223267820315*tmp21*(p_0*v_0 - nai_200_000_2*u_0); // final, weighs 6
  tmp14 = nai_000_020_1*p_2 - nai_100_000_2*u_2; // local+recycle, weighs 4
  nai_000_002_1 = nai_100_000_2*p_2 - tmp19*u_2; // local+recycle, weighs 4
  out[20] = 4.69223267820315*tmp21*(p_0*tmp14 - nai_000_002_1*u_0); // final, weighs 6
  out[21] = 2.70906179986095*tmp21*(p_1*v_0 + (nai_010_010_0 + nai_200_000_1 - nai_010_010_1)/ab_a - nai_200_000_2*u_1); // final, weighs 12
  out[22] = 4.69223267820315*tmp21*(p_1*tmp14 + (nai_110_000_2 - nai_010_001_1)/ab_a - nai_000_002_1*u_1); // final, weighs 11
  out[23] = 2.70906179986095*tmp21*(nai_000_002_2 + p_2*tmp14 - nai_000_002_1*u_2); // final, weighs 7
  out[24] = 4.69223267820315*tmp21*(tmp17*v_1 - tmp33*u_1); // final, weighs 6
  nai_100_000_0 = -nai_000_000_3 + tmp11*v_1; // local+recycle, weighs 3
  nai_010_000_3 = -tmp24 + d_1*v_1; // local+recycle, weighs 3
  nai_110_000_2 = nai_000_000_2*v_1 - d_2*u_1; // local+recycle, weighs 4
  out[25] = 8.12718539958285*tmp21*(p_0*(nai_100_000_0*p_1 + (0.5*tmp11 - 0.5*d_1)/ab_a - nai_010_000_3*u_1) - u_0*(nai_010_000_3*p_1 + (0.5*d_1 - 0.5*nai_000_000_2)/ab_a - nai_110_000_2*u_1)); // final, weighs 26
  nai_010_001_1 = nai_100_000_0*p_2 + (0.5*d_0 - 0.5*tmp39)/ab_a - nai_010_000_3*u_2; // local+recycle, weighs 10
  tmp17 = nai_010_000_3*p_2 + (0.5*tmp39 - 0.5*nai_010_000_2)/ab_a - nai_110_000_2*u_2; // local+recycle, weighs 10
  out[26] = 8.12718539958285*tmp21*(nai_010_001_1*p_0 - tmp17*u_0); // final, weighs 6
  out[27] = 4.69223267820315*tmp21*(nai_000_020_0*v_1 + (usq - nai_001_010_1)/ab_a - tmp29*u_1); // final, weighs 11
  out[28] = 8.12718539958285*tmp21*(nai_100_001_1 + nai_010_001_1*p_1 - tmp17*u_1); // final, weighs 7
  out[29] = 4.69223267820315*tmp21*(tmp36*v_1 - nai_200_000_0*u_1); // final, weighs 6
  nai_200_000_0 = nai_100_000_1 + tmp11*v_2; // local+recycle, weighs 2
  nai_000_020_1 = tmp22 + d_1*v_2; // local+recycle, weighs 2
  nai_010_000_2 = tmp2 + nai_000_000_2*v_2 - d_2*u_2; // local+recycle, weighs 5
  tmp22 = 0.5*nai_200_000_0 - 0.5*nai_000_020_1; // auto+recycle, weighs 3
  nai_000_000_3 = tmp22/ab_a; // auto+recycle, weighs 2
  out[30] = 2.70906179986095*tmp21*(nai_000_000_3 + p_0*(nai_200_000_0*p_0 - nai_000_020_1*u_0) - u_0*(nai_000_020_1*p_0 - nai_010_000_2*u_0)); // final, weighs 15
  tmp33 = nai_200_000_0*p_1 - nai_000_020_1*u_1; // local+recycle, weighs 4
  nai_100_000_2 = nai_000_020_1*p_1 - nai_010_000_2*u_1; // local+recycle, weighs 4
  out[31] = 4.69223267820315*tmp21*(p_0*tmp33 - nai_100_000_2*u_0); // final, weighs 6
  nai_010_010_0 = nai_200_000_0*p_2 + (tmp11 - d_1)/ab_a - nai_000_020_1*u_2; // local+recycle, weighs 9
  tmp19 = nai_000_020_1*p_2 + (d_1 - nai_000_000_2)/ab_a - nai_010_000_2*u_2; // local+recycle, weighs 9
  out[32] = 4.69223267820315*tmp21*(nai_010_010_0*p_0 - tmp19*u_0); // final, weighs 6
  out[33] = 2.70906179986095*tmp21*(nai_000_000_3 + p_1*tmp33 - nai_100_000_2*u_1); // final, weighs 7
  out[34] = 4.69223267820315*tmp21*(nai_010_010_0*p_1 - tmp19*u_1); // final, weighs 6
  out[35] = 2.70906179986095*tmp21*(nai_010_010_0*p_2 + (tmp1 + tmp22 - tmp0)/ab_a - tmp19*u_2); // final, weighs 12
  // total weight = 928
}

void gint2_nai_cF_cD(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place permutation of the output of another routine.
  double tmp;
  gint2_nai_cD_cF(b, b_a, a, a_a, c, out);
  tmp = out[1];
  out[1] = out[10];
  out[10] = out[41];
  out[41] = out[56];
  out[56] = out[29];
  out[29] = out[54];
  out[54] = out[9];
  out[9] = out[31];
  out[31] = out[15];
  out[15] = out[32];
  out[32] = out[25];
  out[25] = out[14];
  out[14] = out[22];
  out[22] = out[43];
  out[43] = out[17];
  out[17] = out[52];
  out[52] = out[48];
  out[48] = out[8];
  out[8] = out[21];
  out[21] = out[33];
  out[33] = out[35];
  out[35] = out[55];
  out[55] = out[19];
  out[19] = out[13];
  out[13] = out[12];
  out[12] = out[2];
  out[2] = out[20];
  out[20] = out[23];
  out[23] = out[53];
  out[53] = out[58];
  out[58] = out[49];
  out[49] = out[18];
  out[18] = out[3];
  out[3] = out[30];
  out[30] = out[5];
  out[5] = out[50];
  out[50] = out[28];
  out[28] = out[44];
  out[44] = out[27];
  out[27] = out[34];
  out[34] = out[45];
  out[45] = out[37];
  out[37] = out[16];
  out[16] = out[42];
  out[42] = out[7];
  out[7] = out[11];
  out[11] = out[51];
  out[51] = out[38];
  out[38] = out[26];
  out[26] = out[24];
  out[24] = out[4];
  out[4] = out[40];
  out[40] = out[46];
  out[46] = out[47];
  out[47] = out[57];
  out[57] = out[39];
  out[39] = out[36];
  out[36] = out[6];
  out[6] = tmp;
}

void gint2_nai_pF_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 15
  double tmp0, tmp1, tmp10, tmp11, tmp12, tmp13, tmp14, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  gint2_nai_cF_cF(a, a_a, b, b_a, c, out);
  tmp0 = 1.09544511501033*out[50] - 0.273861278752583*out[30] - 0.612372435695794*out[0]; // stacked, weighs 5
  tmp1 = out[90] - 0.670820393249937*out[20] - 0.670820393249937*out[70]; // stacked, weighs 4
  tmp2 = 0.866025403784439*out[20] - 0.866025403784439*out[70]; // stacked, weighs 3
  out[20] = 1.09544511501033*out[80] - 0.273861278752583*out[10] - 0.612372435695794*out[60]; // final, weighs 5
  tmp3 = 1.09544511501033*out[52] - 0.273861278752583*out[32] - 0.612372435695794*out[2]; // stacked, weighs 5
  tmp4 = 0.866025403784439*out[22] - 0.866025403784439*out[72]; // stacked, weighs 3
  tmp5 = out[92] - 0.670820393249937*out[22] - 0.670820393249937*out[72]; // stacked, weighs 4
  out[22] = 1.09544511501033*out[82] - 0.273861278752583*out[12] - 0.612372435695794*out[62]; // final, weighs 5
  tmp6 = 1.09544511501033*out[54] - 0.273861278752583*out[34] - 0.612372435695794*out[4]; // stacked, weighs 5
  tmp7 = out[94] - 0.670820393249937*out[24] - 0.670820393249937*out[74]; // stacked, weighs 4
  tmp8 = 0.866025403784439*out[24] - 0.866025403784439*out[74]; // stacked, weighs 3
  out[24] = 1.09544511501033*out[84] - 0.273861278752583*out[14] - 0.612372435695794*out[64]; // final, weighs 5
  tmp9 = 1.09544511501033*out[56] - 0.273861278752583*out[36] - 0.612372435695794*out[6]; // stacked, weighs 5
  tmp10 = 0.866025403784439*out[26] - 0.866025403784439*out[76]; // stacked, weighs 3
  tmp11 = out[96] - 0.670820393249937*out[26] - 0.670820393249937*out[76]; // stacked, weighs 4
  out[26] = 1.09544511501033*out[86] - 0.273861278752583*out[16] - 0.612372435695794*out[66]; // final, weighs 5
  tmp12 = 1.09544511501033*out[58] - 0.273861278752583*out[38] - 0.612372435695794*out[8]; // stacked, weighs 5
  tmp13 = out[98] - 0.670820393249937*out[28] - 0.670820393249937*out[78]; // stacked, weighs 4
  tmp14 = 0.866025403784439*out[28] - 0.866025403784439*out[78]; // stacked, weighs 3
  out[28] = 1.09544511501033*out[88] - 0.273861278752583*out[18] - 0.612372435695794*out[68]; // final, weighs 5
  out[50] = 0.790569415042095*out[0] - 1.06066017177982*out[30]; // final, weighs 3
  out[30] = tmp2; // final, weighs 0
  out[0] = tmp1; // final, weighs 0
  tmp2 = 1.09544511501033*out[51] - 0.273861278752583*out[31] - 0.612372435695794*out[1]; // stacked, weighs 5
  out[51] = 0.790569415042095*out[1] - 1.06066017177982*out[31]; // final, weighs 3
  out[1] = out[91] - 0.670820393249937*out[21] - 0.670820393249937*out[71]; // final, weighs 4
  out[31] = 0.866025403784439*out[21] - 0.866025403784439*out[71]; // final, weighs 3
  out[21] = 1.09544511501033*out[81] - 0.273861278752583*out[11] - 0.612372435695794*out[61]; // final, weighs 5
  out[52] = 0.790569415042095*out[2] - 1.06066017177982*out[32]; // final, weighs 3
  out[2] = tmp5; // final, weighs 0
  out[32] = tmp4; // final, weighs 0
  tmp1 = 1.09544511501033*out[53] - 0.273861278752583*out[33] - 0.612372435695794*out[3]; // stacked, weighs 5
  out[53] = 0.790569415042095*out[3] - 1.06066017177982*out[33]; // final, weighs 3
  out[3] = out[93] - 0.670820393249937*out[23] - 0.670820393249937*out[73]; // final, weighs 4
  out[33] = 0.866025403784439*out[23] - 0.866025403784439*out[73]; // final, weighs 3
  out[23] = 1.09544511501033*out[83] - 0.273861278752583*out[13] - 0.612372435695794*out[63]; // final, weighs 5
  out[54] = 0.790569415042095*out[4] - 1.06066017177982*out[34]; // final, weighs 3
  out[34] = tmp8; // final, weighs 0
  out[4] = tmp7; // final, weighs 0
  tmp5 = 1.09544511501033*out[55] - 0.273861278752583*out[35] - 0.612372435695794*out[5]; // stacked, weighs 5
  out[55] = 0.790569415042095*out[5] - 1.06066017177982*out[35]; // final, weighs 3
  out[5] = out[95] - 0.670820393249937*out[25] - 0.670820393249937*out[75]; // final, weighs 4
  out[35] = 0.866025403784439*out[25] - 0.866025403784439*out[75]; // final, weighs 3
  out[25] = 1.09544511501033*out[85] - 0.273861278752583*out[15] - 0.612372435695794*out[65]; // final, weighs 5
  out[56] = 0.790569415042095*out[6] - 1.06066017177982*out[36]; // final, weighs 3
  out[6] = tmp11; // final, weighs 0
  out[36] = tmp10; // final, weighs 0
  tmp4 = 1.09544511501033*out[57] - 0.273861278752583*out[37] - 0.612372435695794*out[7]; // stacked, weighs 5
  out[57] = 0.790569415042095*out[7] - 1.06066017177982*out[37]; // final, weighs 3
  out[7] = out[97] - 0.670820393249937*out[27] - 0.670820393249937*out[77]; // final, weighs 4
  out[37] = 0.866025403784439*out[27] - 0.866025403784439*out[77]; // final, weighs 3
  out[27] = 1.09544511501033*out[87] - 0.273861278752583*out[17] - 0.612372435695794*out[67]; // final, weighs 5
  out[58] = 0.790569415042095*out[8] - 1.06066017177982*out[38]; // final, weighs 3
  out[38] = tmp14; // final, weighs 0
  out[8] = tmp13; // final, weighs 0
  tmp8 = 1.09544511501033*out[59] - 0.273861278752583*out[39] - 0.612372435695794*out[9]; // stacked, weighs 5
  out[59] = 0.790569415042095*out[9] - 1.06066017177982*out[39]; // final, weighs 3
  out[9] = out[99] - 0.670820393249937*out[29] - 0.670820393249937*out[79]; // final, weighs 4
  out[39] = 0.866025403784439*out[29] - 0.866025403784439*out[79]; // final, weighs 3
  out[29] = 1.09544511501033*out[89] - 0.273861278752583*out[19] - 0.612372435695794*out[69]; // final, weighs 5
  out[60] = 1.06066017177982*out[10] - 0.790569415042095*out[60]; // final, weighs 3
  out[10] = tmp0; // final, weighs 0
  out[61] = 1.06066017177982*out[11] - 0.790569415042095*out[61]; // final, weighs 3
  out[11] = tmp2; // final, weighs 0
  out[62] = 1.06066017177982*out[12] - 0.790569415042095*out[62]; // final, weighs 3
  out[12] = tmp3; // final, weighs 0
  out[63] = 1.06066017177982*out[13] - 0.790569415042095*out[63]; // final, weighs 3
  out[13] = tmp1; // final, weighs 0
  out[64] = 1.06066017177982*out[14] - 0.790569415042095*out[64]; // final, weighs 3
  out[14] = tmp6; // final, weighs 0
  out[65] = 1.06066017177982*out[15] - 0.790569415042095*out[65]; // final, weighs 3
  out[15] = tmp5; // final, weighs 0
  out[66] = 1.06066017177982*out[16] - 0.790569415042095*out[66]; // final, weighs 3
  out[16] = tmp9; // final, weighs 0
  out[67] = 1.06066017177982*out[17] - 0.790569415042095*out[67]; // final, weighs 3
  out[17] = tmp4; // final, weighs 0
  out[68] = 1.06066017177982*out[18] - 0.790569415042095*out[68]; // final, weighs 3
  out[18] = tmp12; // final, weighs 0
  out[69] = 1.06066017177982*out[19] - 0.790569415042095*out[69]; // final, weighs 3
  out[19] = tmp8; // final, weighs 0
  // total weight = 230
}

void gint2_nai_pD_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // In-place modification of the output of another routine.
  // Number of local variables: 11
  double tmp0, tmp1, tmp10, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  gint2_nai_cD_cF(a, a_a, b, b_a, c, out);
  tmp0 = out[20]; // stacked, weighs 0
  out[20] = out[40]; // final, weighs 0
  tmp1 = out[21]; // stacked, weighs 0
  out[21] = out[41]; // final, weighs 0
  tmp2 = out[22]; // stacked, weighs 0
  out[22] = out[42]; // final, weighs 0
  tmp3 = out[23]; // stacked, weighs 0
  out[23] = out[43]; // final, weighs 0
  tmp4 = out[24]; // stacked, weighs 0
  out[24] = out[44]; // final, weighs 0
  tmp5 = out[25]; // stacked, weighs 0
  out[25] = out[45]; // final, weighs 0
  tmp6 = out[26]; // stacked, weighs 0
  out[26] = out[46]; // final, weighs 0
  tmp7 = out[27]; // stacked, weighs 0
  out[27] = out[47]; // final, weighs 0
  tmp8 = out[28]; // stacked, weighs 0
  out[28] = out[48]; // final, weighs 0
  tmp9 = out[29]; // stacked, weighs 0
  out[29] = out[49]; // final, weighs 0
  tmp10 = out[50] - 0.5*out[0] - 0.5*out[30]; // stacked, weighs 4
  out[30] = 0.866025403784439*out[0] - 0.866025403784439*out[30]; // final, weighs 3
  out[0] = tmp10; // final, weighs 0
  tmp10 = out[51] - 0.5*out[1] - 0.5*out[31]; // stacked, weighs 4
  out[31] = 0.866025403784439*out[1] - 0.866025403784439*out[31]; // final, weighs 3
  out[1] = tmp10; // final, weighs 0
  tmp10 = out[52] - 0.5*out[2] - 0.5*out[32]; // stacked, weighs 4
  out[32] = 0.866025403784439*out[2] - 0.866025403784439*out[32]; // final, weighs 3
  out[2] = tmp10; // final, weighs 0
  tmp10 = out[53] - 0.5*out[3] - 0.5*out[33]; // stacked, weighs 4
  out[33] = 0.866025403784439*out[3] - 0.866025403784439*out[33]; // final, weighs 3
  out[3] = tmp10; // final, weighs 0
  tmp10 = out[54] - 0.5*out[34] - 0.5*out[4]; // stacked, weighs 4
  out[34] = 0.866025403784439*out[4] - 0.866025403784439*out[34]; // final, weighs 3
  out[4] = tmp10; // final, weighs 0
  tmp10 = out[55] - 0.5*out[35] - 0.5*out[5]; // stacked, weighs 4
  out[35] = 0.866025403784439*out[5] - 0.866025403784439*out[35]; // final, weighs 3
  out[5] = tmp10; // final, weighs 0
  tmp10 = out[56] - 0.5*out[36] - 0.5*out[6]; // stacked, weighs 4
  out[36] = 0.866025403784439*out[6] - 0.866025403784439*out[36]; // final, weighs 3
  out[6] = tmp10; // final, weighs 0
  tmp10 = out[57] - 0.5*out[37] - 0.5*out[7]; // stacked, weighs 4
  out[37] = 0.866025403784439*out[7] - 0.866025403784439*out[37]; // final, weighs 3
  out[7] = tmp10; // final, weighs 0
  tmp10 = out[58] - 0.5*out[38] - 0.5*out[8]; // stacked, weighs 4
  out[38] = 0.866025403784439*out[8] - 0.866025403784439*out[38]; // final, weighs 3
  out[8] = tmp10; // final, weighs 0
  tmp10 = out[59] - 0.5*out[39] - 0.5*out[9]; // stacked, weighs 4
  out[39] = 0.866025403784439*out[9] - 0.866025403784439*out[39]; // final, weighs 3
  out[9] = tmp10; // final, weighs 0
  out[40] = out[10]; // final, weighs 0
  out[10] = tmp0; // final, weighs 0
  out[41] = out[11]; // final, weighs 0
  out[11] = tmp1; // final, weighs 0
  out[42] = out[12]; // final, weighs 0
  out[12] = tmp2; // final, weighs 0
  out[43] = out[13]; // final, weighs 0
  out[13] = tmp3; // final, weighs 0
  out[44] = out[14]; // final, weighs 0
  out[14] = tmp4; // final, weighs 0
  out[45] = out[15]; // final, weighs 0
  out[15] = tmp5; // final, weighs 0
  out[46] = out[16]; // final, weighs 0
  out[16] = tmp6; // final, weighs 0
  out[47] = out[17]; // final, weighs 0
  out[17] = tmp7; // final, weighs 0
  out[48] = out[18]; // final, weighs 0
  out[18] = tmp8; // final, weighs 0
  out[49] = out[19]; // final, weighs 0
  out[19] = tmp9; // final, weighs 0
  // total weight = 70
}

void gint2_nai_SP_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 60
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_000_001_0, nai_000_001_1, nai_000_001_2, nai_000_002_0, nai_000_002_1, nai_000_003_0, nai_000_010_0, nai_000_010_1, nai_000_010_2, nai_000_011_0, nai_000_011_1, nai_000_020_0, nai_000_020_1, nai_000_030_0, nai_000_100_0, nai_000_100_1, nai_000_100_2, nai_000_200_0, nai_000_200_1, nai_000_200_2, nai_000_300_0, nai_100_001_1, nai_100_020_1, p_0, p_1, p_2, tmp1, tmp12, tmp19, tmp2, tmp23, tmp24, tmp25, tmp26, tmp27, tmp28, tmp29, tmp33, tmp34, tmp35, tmp36, tmp37, tmp38, tmp39, tmp40, tmp41, tmp5, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  tmp12 = pow(b_a,2.25); // auto, weighs 2
  tmp23 = tmp12*pow(a_a,0.75); // auto, weighs 3
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
  usq = 6.28318530717959*d_0*gaux(usq, 4); // local+recycle, weighs 3
  d_0 = u_0*usq; // auto+recycle, weighs 1
  tmp23 = (0.5*nai_000_000_2 - 0.5*nai_000_000_3)/ab_a; // auto+recycle, weighs 5
  nai_000_200_2 = tmp23 + nai_000_100_2*p_0 - u_0*(-d_0 + nai_000_000_3*p_0); // local, weighs 8
  nai_000_100_1 = (nai_000_100_1 - nai_000_100_2)/ab_a - nai_000_200_2*u_0; // auto+recycle, weighs 7
  nai_000_100_2 = nai_000_100_1 + nai_000_200_1*p_0; // local+recycle, weighs 2
  tmp12 = tmp12*pow(a_a,1.25); // auto+recycle, weighs 3
  out[10] = 2.09843024694163*tmp12*(nai_000_300_0*v_0 + (1.5*nai_000_200_0 - 1.5*nai_000_200_1)/ab_a - nai_000_100_2*u_0); // final, weighs 12
  nai_000_100_0 = nai_000_100_0 + nai_000_200_0*v_0; // local+recycle, weighs 2
  nai_000_100_1 = nai_000_100_1 + nai_000_200_1*v_0; // local+recycle, weighs 2
  out[11] = 4.69223267820315*tmp12*(nai_000_100_0*p_1 - nai_000_100_1*u_1); // final, weighs 6
  out[12] = 4.69223267820315*tmp12*(nai_000_100_0*p_2 - nai_000_100_1*u_2); // final, weighs 6
  nai_000_100_0 = -tmp27 + nai_000_020_0*v_0; // local+recycle, weighs 3
  tmp27 = u_1*usq; // auto+recycle, weighs 1
  nai_000_100_1 = tmp23 + nai_000_010_2*p_1 - u_1*(-tmp27 + nai_000_000_3*p_1); // local+recycle, weighs 8
  nai_100_020_1 = nai_000_020_1*v_0 - nai_000_100_1*u_0; // local, weighs 4
  tmp5 = (0.5*nai_000_020_0 - 0.5*nai_000_020_1)/ab_a; // auto, weighs 5
  out[13] = 4.69223267820315*tmp12*(tmp5 + nai_000_100_0*p_0 - nai_100_020_1*u_0); // final, weighs 7
  tmp38 = -tmp38 + d_2*v_0; // local+recycle, weighs 3
  tmp35 = -tmp35 + nai_000_000_2*v_0; // local+recycle, weighs 3
  nai_100_001_1 = p_2*tmp38 - tmp35*u_2; // local, weighs 4
  out[14] = 8.12718539958285*tmp12*(p_0*(p_1*(p_2*(-tmp41 + d_1*v_0) - tmp38*u_2) - nai_100_001_1*u_1) + (0.5*nai_000_011_0 - 0.5*nai_000_011_1)/ab_a - u_0*(nai_100_001_1*p_1 - u_1*(p_2*tmp35 - u_2*(-d_0 + nai_000_000_3*v_0)))); // final, weighs 34
  d_0 = -tmp29 + nai_000_002_0*v_0; // local+recycle, weighs 3
  tmp41 = u_2*usq; // auto+recycle, weighs 1
  nai_000_011_0 = tmp23 + nai_000_001_2*p_2 - u_2*(-tmp41 + nai_000_000_3*p_2); // local+recycle, weighs 8
  nai_000_011_1 = nai_000_002_1*v_0 - nai_000_011_0*u_0; // local+recycle, weighs 4
  tmp29 = (0.5*nai_000_002_0 - 0.5*nai_000_002_1)/ab_a; // auto+recycle, weighs 5
  out[15] = 4.69223267820315*tmp12*(tmp29 + d_0*p_0 - nai_000_011_1*u_0); // final, weighs 7
  tmp38 = (nai_000_010_1 - nai_000_010_2)/ab_a - nai_000_100_1*u_1; // auto+recycle, weighs 7
  tmp35 = tmp38 + nai_000_020_1*p_1; // local+recycle, weighs 2
  out[16] = 2.09843024694163*tmp12*(nai_000_030_0*v_0 - tmp35*u_0); // final, weighs 6
  out[17] = 4.69223267820315*tmp12*(nai_000_100_0*p_2 - nai_100_020_1*u_2); // final, weighs 6
  out[18] = 4.69223267820315*tmp12*(d_0*p_1 - nai_000_011_1*u_1); // final, weighs 6
  nai_100_001_1 = (nai_000_001_1 - nai_000_001_2)/ab_a - nai_000_011_0*u_2; // auto+recycle, weighs 7
  usq = nai_100_001_1 + nai_000_002_1*p_2; // local+recycle, weighs 2
  out[19] = 2.09843024694163*tmp12*(nai_000_003_0*v_0 - u_0*usq); // final, weighs 6
  out[20] = 2.09843024694163*tmp12*(nai_000_300_0*v_1 - nai_000_100_2*u_1); // final, weighs 6
  v_0 = -tmp25 + nai_000_200_0*v_1; // local+recycle, weighs 3
  d_0 = nai_000_200_1*v_1 - nai_000_200_2*u_1; // local+recycle, weighs 4
  nai_000_100_0 = (0.5*nai_000_200_0 - 0.5*nai_000_200_1)/ab_a; // auto+recycle, weighs 5
  out[21] = 4.69223267820315*tmp12*(nai_000_100_0 + p_1*v_0 - d_0*u_1); // final, weighs 7
  out[22] = 4.69223267820315*tmp12*(p_2*v_0 - d_0*u_2); // final, weighs 6
  tmp25 = nai_000_010_0 + nai_000_020_0*v_1; // local+recycle, weighs 2
  nai_000_010_0 = tmp38 + nai_000_020_1*v_1; // local+recycle, weighs 2
  out[23] = 4.69223267820315*tmp12*(p_0*tmp25 - nai_000_010_0*u_0); // final, weighs 6
  tmp37 = -tmp37 + d_2*v_1; // local+recycle, weighs 3
  nai_000_010_1 = -tmp34 + nai_000_000_2*v_1; // local+recycle, weighs 3
  tmp34 = p_2*tmp37 - nai_000_010_1*u_2; // local+recycle, weighs 4
  out[24] = 8.12718539958285*tmp12*(p_0*(p_1*(p_2*(-tmp40 + d_1*v_1) - tmp37*u_2) + (0.5*nai_000_001_0 - 0.5*nai_000_001_1)/ab_a - tmp34*u_1) - u_0*(p_1*tmp34 + (0.5*nai_000_001_1 - 0.5*nai_000_001_2)/ab_a - u_1*(nai_000_010_1*p_2 - u_2*(-tmp27 + nai_000_000_3*v_1)))); // final, weighs 40
  nai_000_010_2 = -tmp28 + nai_000_002_0*v_1; // local+recycle, weighs 3
  tmp27 = nai_000_002_1*v_1 - nai_000_011_0*u_1; // local+recycle, weighs 4
  out[25] = 4.69223267820315*tmp12*(nai_000_010_2*p_0 - tmp27*u_0); // final, weighs 6
  out[26] = 2.09843024694163*tmp12*(nai_000_030_0*v_1 + (1.5*nai_000_020_0 - 1.5*nai_000_020_1)/ab_a - tmp35*u_1); // final, weighs 12
  out[27] = 4.69223267820315*tmp12*(p_2*tmp25 - nai_000_010_0*u_2); // final, weighs 6
  out[28] = 4.69223267820315*tmp12*(tmp29 + nai_000_010_2*p_1 - tmp27*u_1); // final, weighs 7
  out[29] = 2.09843024694163*tmp12*(nai_000_003_0*v_1 - u_1*usq); // final, weighs 6
  out[30] = 2.09843024694163*tmp12*(nai_000_300_0*v_2 - nai_000_100_2*u_2); // final, weighs 6
  nai_000_001_0 = -tmp24 + nai_000_200_0*v_2; // local+recycle, weighs 3
  nai_000_001_1 = nai_000_200_1*v_2 - nai_000_200_2*u_2; // local+recycle, weighs 4
  out[31] = 4.69223267820315*tmp12*(nai_000_001_0*p_1 - nai_000_001_1*u_1); // final, weighs 6
  out[32] = 4.69223267820315*tmp12*(nai_000_100_0 + nai_000_001_0*p_2 - nai_000_001_1*u_2); // final, weighs 7
  nai_000_011_0 = -tmp26 + nai_000_020_0*v_2; // local+recycle, weighs 3
  nai_000_001_2 = nai_000_020_1*v_2 - nai_000_100_1*u_2; // local+recycle, weighs 4
  out[33] = 4.69223267820315*tmp12*(nai_000_011_0*p_0 - nai_000_001_2*u_0); // final, weighs 6
  nai_000_011_1 = -tmp36 + d_2*v_2; // local+recycle, weighs 3
  tmp29 = -tmp33 + nai_000_000_2*v_2; // local+recycle, weighs 3
  nai_000_000_2 = tmp1 + nai_000_011_1*p_2 - tmp29*u_2; // local+recycle, weighs 5
  out[34] = 8.12718539958285*tmp12*(p_0*(p_1*(tmp2 + p_2*(-tmp39 + d_1*v_2) - nai_000_011_1*u_2) - nai_000_000_2*u_1) - u_0*(nai_000_000_2*p_1 - u_1*(tmp23 + p_2*tmp29 - u_2*(-tmp41 + nai_000_000_3*v_2)))); // final, weighs 30
  tmp38 = tmp19 + nai_000_002_0*v_2; // local+recycle, weighs 2
  nai_000_100_1 = nai_100_001_1 + nai_000_002_1*v_2; // local+recycle, weighs 2
  out[35] = 4.69223267820315*tmp12*(p_0*tmp38 - nai_000_100_1*u_0); // final, weighs 6
  out[36] = 2.09843024694163*tmp12*(nai_000_030_0*v_2 - tmp35*u_2); // final, weighs 6
  out[37] = 4.69223267820315*tmp12*(tmp5 + nai_000_011_0*p_2 - nai_000_001_2*u_2); // final, weighs 7
  out[38] = 4.69223267820315*tmp12*(p_1*tmp38 - nai_000_100_1*u_1); // final, weighs 6
  out[39] = 2.09843024694163*tmp12*(nai_000_003_0*v_2 + (1.5*nai_000_002_0 - 1.5*nai_000_002_1)/ab_a - u_2*usq); // final, weighs 12
  // total weight = 690
}

void gint2_nai_S_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 18
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_100_0, nai_000_100_1, nai_000_200_0, nai_000_200_1, p_0, p_1, p_2, tmp1, tmp5, u_0, u_1, u_2, usq;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  usq = 6.28318530717959*d_0*gaux(usq, 3); // local+recycle, weighs 3
  d_0 = (0.5*d_2 - 0.5*nai_000_000_2)/ab_a; // auto+recycle, weighs 5
  nai_000_200_1 = d_0 + nai_000_100_1*p_0 - u_0*(nai_000_000_2*p_0 - u_0*usq); // local, weighs 9
  tmp5 = pow(a_a,0.75)*pow(b_a,2.25); // auto, weighs 5
  out[0] = 1.04921512347081*tmp5*(nai_000_200_0*p_0 + (nai_000_100_0 - nai_000_100_1)/ab_a - nai_000_200_1*u_0); // final, weighs 11
  out[1] = 2.34611633910158*tmp5*(nai_000_200_0*p_1 - nai_000_200_1*u_1); // final, weighs 6
  out[2] = 2.34611633910158*tmp5*(nai_000_200_0*p_2 - nai_000_200_1*u_2); // final, weighs 6
  nai_000_100_0 = d_1*p_1 - d_2*u_1; // local+recycle, weighs 4
  nai_000_100_1 = d_2*p_1 - nai_000_000_2*u_1; // local+recycle, weighs 4
  nai_000_200_0 = tmp1 + nai_000_100_0*p_1 - nai_000_100_1*u_1; // local+recycle, weighs 5
  nai_000_200_1 = d_0 + nai_000_100_1*p_1 - u_1*(nai_000_000_2*p_1 - u_1*usq); // local+recycle, weighs 9
  out[3] = 2.34611633910158*tmp5*(nai_000_200_0*p_0 - nai_000_200_1*u_0); // final, weighs 6
  d_1 = d_1*p_2 - d_2*u_2; // local+recycle, weighs 4
  d_2 = d_2*p_2 - nai_000_000_2*u_2; // local+recycle, weighs 4
  nai_000_000_2 = nai_000_000_2*p_2 - u_2*usq; // local+recycle, weighs 4
  out[4] = 4.06359269979142*tmp5*(p_0*(d_1*p_1 - d_2*u_1) - u_0*(d_2*p_1 - nai_000_000_2*u_1)); // final, weighs 14
  tmp1 = tmp1 + d_1*p_2 - d_2*u_2; // local+recycle, weighs 5
  usq = d_0 + d_2*p_2 - nai_000_000_2*u_2; // local+recycle, weighs 5
  out[5] = 2.34611633910158*tmp5*(p_0*tmp1 - u_0*usq); // final, weighs 6
  out[6] = 1.04921512347081*tmp5*(nai_000_200_0*p_1 + (nai_000_100_0 - nai_000_100_1)/ab_a - nai_000_200_1*u_1); // final, weighs 11
  out[7] = 2.34611633910158*tmp5*(nai_000_200_0*p_2 - nai_000_200_1*u_2); // final, weighs 6
  out[8] = 2.34611633910158*tmp5*(p_1*tmp1 - u_1*usq); // final, weighs 6
  out[9] = 1.04921512347081*tmp5*(p_2*tmp1 + (d_1 - d_2)/ab_a - u_2*usq); // final, weighs 11
  // total weight = 229
}

void gint2_nai_P_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 52
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_000_001_0, nai_000_001_1, nai_000_001_2, nai_000_010_1, nai_000_010_2, nai_000_020_0, nai_000_020_1, nai_000_020_2, nai_000_030_0, nai_000_100_0, nai_000_100_1, nai_000_100_2, nai_000_200_0, nai_000_200_1, nai_000_200_2, nai_000_300_0, nai_010_200_1, nai_100_001_1, nai_100_020_0, nai_100_020_1, p_0, p_1, p_2, tmp0, tmp1, tmp2, tmp21, tmp23, tmp25, tmp26, tmp27, tmp28, tmp29, tmp3, tmp30, tmp31, tmp33, tmp4, tmp5, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  usq = 6.28318530717959*d_0*gaux(usq, 4); // local+recycle, weighs 3
  d_0 = u_0*usq; // auto+recycle, weighs 1
  tmp0 = (0.5*nai_000_000_2 - 0.5*nai_000_000_3)/ab_a; // auto, weighs 5
  nai_000_200_2 = tmp0 + nai_000_100_2*p_0 - u_0*(-d_0 + nai_000_000_3*p_0); // local, weighs 8
  nai_000_100_1 = (nai_000_100_1 - nai_000_100_2)/ab_a - nai_000_200_2*u_0; // auto+recycle, weighs 7
  nai_000_100_2 = nai_000_100_1 + nai_000_200_1*p_0; // local+recycle, weighs 2
  tmp21 = pow(a_a,1.25)*pow(b_a,2.25); // auto, weighs 5
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
  tmp23 = u_1*usq; // auto, weighs 1
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
  out[4] = 8.12718539958285*tmp21*(p_0*(p_1*(p_2*(-tmp33 + d_1*v_0) - tmp30*u_2) - nai_100_001_1*u_1) + (0.5*nai_000_001_0*p_1 + 0.5*nai_000_001_2*u_1 - 0.5*nai_000_001_1*p_1 - 0.5*nai_000_001_1*u_1)/ab_a - u_0*(nai_100_001_1*p_1 - u_1*(p_2*tmp27 - u_2*(-d_0 + nai_000_000_3*v_0)))); // final, weighs 42
  d_0 = tmp2 + nai_000_001_0*p_2 - nai_000_001_1*u_2; // local+recycle, weighs 5
  tmp33 = tmp1 + nai_000_001_1*p_2 - nai_000_001_2*u_2; // local+recycle, weighs 5
  tmp30 = d_0*v_0 - tmp33*u_0; // local+recycle, weighs 4
  tmp27 = u_2*usq; // auto+recycle, weighs 1
  nai_100_001_1 = tmp0 + nai_000_001_2*p_2 - u_2*(-tmp27 + nai_000_000_3*p_2); // local+recycle, weighs 8
  usq = tmp33*v_0 - nai_100_001_1*u_0; // local+recycle, weighs 4
  tmp3 = (0.5*d_0 - 0.5*tmp33)/ab_a; // auto, weighs 5
  out[5] = 4.69223267820315*tmp21*(tmp3 + p_0*tmp30 - u_0*usq); // final, weighs 7
  nai_000_100_1 = (nai_000_100_1 - nai_000_010_1)/ab_a - nai_000_020_1*u_1; // auto+recycle, weighs 7
  nai_000_030_0 = nai_000_100_1 + nai_000_020_0*p_1; // local, weighs 2
  nai_000_010_1 = (nai_000_010_1 - nai_000_010_2)/ab_a - nai_000_020_2*u_1; // auto+recycle, weighs 7
  nai_000_010_2 = nai_000_010_1 + nai_000_020_1*p_1; // local+recycle, weighs 2
  out[6] = 2.09843024694163*tmp21*(nai_000_030_0*v_0 - nai_000_010_2*u_0); // final, weighs 6
  out[7] = 4.69223267820315*tmp21*(nai_100_020_0*p_2 - nai_100_020_1*u_2); // final, weighs 6
  out[8] = 4.69223267820315*tmp21*(p_1*tmp30 - u_1*usq); // final, weighs 6
  nai_100_020_0 = (nai_000_001_0 - nai_000_001_1)/ab_a - tmp33*u_2; // auto+recycle, weighs 7
  nai_100_020_1 = nai_100_020_0 + d_0*p_2; // local+recycle, weighs 2
  usq = (nai_000_001_1 - nai_000_001_2)/ab_a - nai_100_001_1*u_2; // auto+recycle, weighs 7
  tmp30 = usq + p_2*tmp33; // local+recycle, weighs 2
  out[9] = 2.09843024694163*tmp21*(nai_100_020_1*v_0 - tmp30*u_0); // final, weighs 6
  out[10] = 2.09843024694163*tmp21*(nai_000_300_0*v_1 - nai_000_100_2*u_1); // final, weighs 6
  v_0 = nai_000_200_0*v_1 - nai_000_200_1*u_1; // local+recycle, weighs 4
  nai_010_200_1 = nai_000_200_1*v_1 - nai_000_200_2*u_1; // local, weighs 4
  tmp4 = (0.5*nai_000_200_0 - 0.5*nai_000_200_1)/ab_a; // auto, weighs 5
  out[11] = 4.69223267820315*tmp21*(tmp4 + p_1*v_0 - nai_010_200_1*u_1); // final, weighs 7
  out[12] = 4.69223267820315*tmp21*(p_2*v_0 - nai_010_200_1*u_2); // final, weighs 6
  nai_000_100_1 = nai_000_100_1 + nai_000_020_0*v_1; // local+recycle, weighs 2
  nai_000_010_1 = nai_000_010_1 + nai_000_020_1*v_1; // local+recycle, weighs 2
  out[13] = 4.69223267820315*tmp21*(nai_000_100_1*p_0 - nai_000_010_1*u_0); // final, weighs 6
  v_0 = -tmp29 + d_2*v_1; // local+recycle, weighs 3
  nai_010_200_1 = -tmp26 + nai_000_000_2*v_1; // local+recycle, weighs 3
  tmp29 = p_2*v_0 - nai_010_200_1*u_2; // local+recycle, weighs 4
  out[14] = 8.12718539958285*tmp21*(p_0*(p_1*(p_2*(-nai_000_100_0 + d_1*v_1) - u_2*v_0) + (0.5*nai_000_001_0 - 0.5*nai_000_001_1)/ab_a - tmp29*u_1) - u_0*(p_1*tmp29 + (0.5*nai_000_001_1 - 0.5*nai_000_001_2)/ab_a - u_1*(nai_010_200_1*p_2 - u_2*(-tmp23 + nai_000_000_3*v_1)))); // final, weighs 40
  tmp26 = d_0*v_1 - tmp33*u_1; // local+recycle, weighs 4
  tmp23 = tmp33*v_1 - nai_100_001_1*u_1; // local+recycle, weighs 4
  out[15] = 4.69223267820315*tmp21*(p_0*tmp26 - tmp23*u_0); // final, weighs 6
  out[16] = 2.09843024694163*tmp21*(nai_000_030_0*v_1 + (1.5*nai_000_020_0 - 1.5*nai_000_020_1)/ab_a - nai_000_010_2*u_1); // final, weighs 12
  out[17] = 4.69223267820315*tmp21*(nai_000_100_1*p_2 - nai_000_010_1*u_2); // final, weighs 6
  out[18] = 4.69223267820315*tmp21*(tmp3 + p_1*tmp26 - tmp23*u_1); // final, weighs 7
  out[19] = 2.09843024694163*tmp21*(nai_100_020_1*v_1 - tmp30*u_1); // final, weighs 6
  out[20] = 2.09843024694163*tmp21*(nai_000_300_0*v_2 - nai_000_100_2*u_2); // final, weighs 6
  nai_000_200_0 = nai_000_200_0*v_2 - nai_000_200_1*u_2; // local+recycle, weighs 4
  nai_100_001_1 = nai_000_200_1*v_2 - nai_000_200_2*u_2; // local+recycle, weighs 4
  out[21] = 4.69223267820315*tmp21*(nai_000_200_0*p_1 - nai_100_001_1*u_1); // final, weighs 6
  out[22] = 4.69223267820315*tmp21*(tmp4 + nai_000_200_0*p_2 - nai_100_001_1*u_2); // final, weighs 7
  nai_000_001_0 = nai_000_020_0*v_2 - nai_000_020_1*u_2; // local+recycle, weighs 4
  nai_000_001_1 = nai_000_020_1*v_2 - nai_000_020_2*u_2; // local+recycle, weighs 4
  out[23] = 4.69223267820315*tmp21*(nai_000_001_0*p_0 - nai_000_001_1*u_0); // final, weighs 6
  nai_000_001_2 = -tmp28 + d_2*v_2; // local+recycle, weighs 3
  nai_000_100_2 = -tmp25 + nai_000_000_2*v_2; // local+recycle, weighs 3
  tmp1 = tmp1 + nai_000_001_2*p_2 - nai_000_100_2*u_2; // local+recycle, weighs 5
  out[24] = 8.12718539958285*tmp21*(p_0*(p_1*(tmp2 + p_2*(-tmp31 + d_1*v_2) - nai_000_001_2*u_2) - tmp1*u_1) - u_0*(p_1*tmp1 - u_1*(tmp0 + nai_000_100_2*p_2 - u_2*(-tmp27 + nai_000_000_3*v_2)))); // final, weighs 30
  nai_000_200_1 = nai_100_020_0 + d_0*v_2; // local+recycle, weighs 2
  usq = usq + tmp33*v_2; // local+recycle, weighs 2
  out[25] = 4.69223267820315*tmp21*(nai_000_200_1*p_0 - u_0*usq); // final, weighs 6
  out[26] = 2.09843024694163*tmp21*(nai_000_030_0*v_2 - nai_000_010_2*u_2); // final, weighs 6
  out[27] = 4.69223267820315*tmp21*(tmp5 + nai_000_001_0*p_2 - nai_000_001_1*u_2); // final, weighs 7
  out[28] = 4.69223267820315*tmp21*(nai_000_200_1*p_1 - u_1*usq); // final, weighs 6
  out[29] = 2.09843024694163*tmp21*(nai_100_020_1*v_2 + (1.5*d_0 - 1.5*tmp33)/ab_a - tmp30*u_2); // final, weighs 12
  // total weight = 645
}

void gint2_nai_cD_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 90
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_000_000_4, nai_000_001_1, nai_000_001_2, nai_000_001_3, nai_000_002_0, nai_000_002_1, nai_000_002_2, nai_000_002_3, nai_000_003_2, nai_000_010_1, nai_000_010_2, nai_000_010_3, nai_000_020_0, nai_000_020_1, nai_000_020_2, nai_000_020_3, nai_000_030_0, nai_000_030_1, nai_000_100_0, nai_000_100_1, nai_000_100_2, nai_000_100_3, nai_000_200_0, nai_000_200_1, nai_000_200_2, nai_000_200_3, nai_000_300_0, nai_000_300_1, nai_010_001_0, nai_010_001_1, nai_010_001_2, nai_010_100_1, nai_010_200_0, nai_010_200_1, nai_010_300_0, nai_010_300_1, nai_100_001_1, nai_100_001_2, nai_100_002_0, nai_100_002_1, nai_100_020_0, nai_100_020_1, nai_110_200_0, nai_200_000_0, nai_200_000_1, nai_200_000_2, nai_200_100_1, p_0, p_1, p_2, tmp0, tmp1, tmp15, tmp2, tmp26, tmp27, tmp3, tmp30, tmp31, tmp32, tmp38, tmp43, tmp5, tmp53, tmp54, tmp55, tmp56, tmp64, tmp66, tmp67, tmp68, tmp69, tmp70, tmp71, tmp73, tmp74, tmp77, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  usq = 6.28318530717959*d_0*gaux(usq, 5); // local+recycle, weighs 3
  d_0 = u_0*usq; // auto+recycle, weighs 1
  tmp3 = (0.5*nai_000_000_3 - 0.5*nai_000_000_4)/ab_a; // auto, weighs 5
  nai_000_200_3 = tmp3 + nai_000_100_3*p_0 - u_0*(-d_0 + nai_000_000_4*p_0); // local, weighs 8
  nai_000_100_3 = nai_000_200_2*p_0 + (nai_000_100_2 - nai_000_100_3)/ab_a - nai_000_200_3*u_0; // local+recycle, weighs 9
  nai_000_100_2 = nai_000_100_0 + nai_000_200_0*v_0; // local+recycle, weighs 2
  nai_000_100_0 = nai_000_100_1 + nai_000_200_1*v_0; // local+recycle, weighs 2
  nai_000_100_1 = 0.5*nai_000_300_0 - 0.5*nai_000_300_1; // auto+recycle, weighs 3
  tmp53 = pow(a_a,1.75)*pow(b_a,2.25); // auto, weighs 5
  out[0] = 2.4230585358948*tmp53*(v_0*(nai_000_300_0*v_0 + (1.5*nai_000_200_0 - 1.5*nai_000_200_1)/ab_a - nai_000_300_1*u_0) + (nai_000_100_1 + 1.5*nai_000_100_2 - 1.5*nai_000_100_0)/ab_a - u_0*(nai_000_300_1*v_0 + (1.5*nai_000_200_1 - 1.5*nai_000_200_2)/ab_a - nai_000_100_3*u_0)); // final, weighs 33
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
  d_0 = tmp3 + tmp68*v_0 - u_0*(-d_0 + nai_000_000_4*v_0); // local+recycle, weighs 8
  tmp31 = 0.5*nai_200_000_1 - 0.5*nai_200_000_2; // auto, weighs 3
  tmp26 = nai_200_100_1*p_0 + (tmp27 + tmp31 - tmp26 - p_0*tmp71)/ab_a - u_0*(nai_200_000_2*p_0 + (tmp71 - tmp68)/ab_a - d_0*u_0); // local+recycle, weighs 22
  out[1] = 5.4181235997219*tmp53*(p_1*tmp43 - tmp26*u_1); // final, weighs 6
  out[2] = 5.4181235997219*tmp53*(p_2*tmp43 - tmp26*u_2); // final, weighs 6
  nai_200_100_1 = nai_200_000_1*p_1 - nai_200_000_2*u_1; // local+recycle, weighs 4
  tmp43 = tmp32/ab_a; // auto+recycle, weighs 2
  tmp27 = tmp43 + p_1*(nai_200_000_0*p_1 - nai_200_000_1*u_1) - nai_200_100_1*u_1; // local+recycle, weighs 9
  tmp32 = tmp31/ab_a; // auto+recycle, weighs 2
  tmp31 = tmp32 + nai_200_100_1*p_1 - u_1*(nai_200_000_2*p_1 - d_0*u_1); // local+recycle, weighs 9
  tmp26 = d_2*u_1; // auto+recycle, weighs 1
  nai_200_100_1 = -tmp26 + d_1*p_1; // local+recycle, weighs 3
  tmp73 = nai_000_000_2*u_1; // auto, weighs 1
  nai_000_010_1 = -tmp73 + d_2*p_1; // local, weighs 3
  nai_000_020_0 = tmp2 + nai_200_100_1*p_1 - nai_000_010_1*u_1; // local, weighs 5
  tmp70 = nai_000_000_3*u_1; // auto, weighs 1
  nai_000_010_2 = -tmp70 + nai_000_000_2*p_1; // local, weighs 3
  nai_000_020_1 = tmp1 + nai_000_010_1*p_1 - nai_000_010_2*u_1; // local, weighs 5
  nai_100_020_0 = nai_000_020_0*v_0 - nai_000_020_1*u_0; // local, weighs 4
  tmp67 = nai_000_000_4*u_1; // auto, weighs 1
  nai_000_010_3 = -tmp67 + nai_000_000_3*p_1; // local, weighs 3
  nai_000_020_2 = tmp0 + nai_000_010_2*p_1 - nai_000_010_3*u_1; // local, weighs 5
  nai_100_020_1 = nai_000_020_1*v_0 - nai_000_020_2*u_0; // local, weighs 4
  out[3] = 5.4181235997219*tmp53*(p_0*tmp27 + (nai_100_020_0 - nai_100_020_1)/ab_a - tmp31*u_0); // final, weighs 11
  nai_200_000_0 = nai_200_000_0*p_2 - nai_200_000_1*u_2; // local+recycle, weighs 4
  nai_200_000_1 = nai_200_000_1*p_2 - nai_200_000_2*u_2; // local+recycle, weighs 4
  d_0 = nai_200_000_2*p_2 - d_0*u_2; // local+recycle, weighs 4
  nai_200_000_2 = p_2*tmp77 - tmp74*u_2; // local+recycle, weighs 4
  nai_100_001_1 = p_2*tmp74 - tmp71*u_2; // local, weighs 4
  nai_100_001_2 = p_2*tmp71 - tmp68*u_2; // local, weighs 4
  out[4] = 9.3844653564063*tmp53*(p_0*(nai_200_000_0*p_1 - nai_200_000_1*u_1) + (nai_100_001_2*u_1 + nai_200_000_2*p_1 - nai_100_001_1*p_1 - nai_100_001_1*u_1)/ab_a - u_0*(nai_200_000_1*p_1 - d_0*u_1)); // final, weighs 26
  tmp43 = tmp43 + nai_200_000_0*p_2 - nai_200_000_1*u_2; // local+recycle, weighs 5
  nai_200_000_0 = tmp32 + nai_200_000_1*p_2 - d_0*u_2; // local+recycle, weighs 5
  tmp32 = d_2*u_2; // auto+recycle, weighs 1
  nai_200_000_1 = -tmp32 + d_1*p_2; // local+recycle, weighs 3
  d_0 = nai_000_000_2*u_2; // auto+recycle, weighs 1
  nai_000_001_1 = -d_0 + d_2*p_2; // local, weighs 3
  nai_000_002_0 = tmp2 + nai_200_000_1*p_2 - nai_000_001_1*u_2; // local, weighs 5
  tmp69 = nai_000_000_3*u_2; // auto, weighs 1
  nai_000_001_2 = -tmp69 + nai_000_000_2*p_2; // local, weighs 3
  nai_000_002_1 = tmp1 + nai_000_001_1*p_2 - nai_000_001_2*u_2; // local, weighs 5
  nai_100_002_0 = nai_000_002_0*v_0 - nai_000_002_1*u_0; // local, weighs 4
  tmp66 = nai_000_000_4*u_2; // auto, weighs 1
  nai_000_001_3 = -tmp66 + nai_000_000_3*p_2; // local, weighs 3
  nai_000_002_2 = tmp0 + nai_000_001_2*p_2 - nai_000_001_3*u_2; // local, weighs 5
  nai_100_002_1 = nai_000_002_1*v_0 - nai_000_002_2*u_0; // local, weighs 4
  out[5] = 5.4181235997219*tmp53*(p_0*tmp43 + (nai_100_002_0 - nai_100_002_1)/ab_a - nai_200_000_0*u_0); // final, weighs 11
  nai_200_100_1 = (nai_200_100_1 - nai_000_010_1)/ab_a - nai_000_020_1*u_1; // auto+recycle, weighs 7
  nai_000_030_0 = nai_200_100_1 + nai_000_020_0*p_1; // local, weighs 2
  nai_000_010_1 = (nai_000_010_1 - nai_000_010_2)/ab_a - nai_000_020_2*u_1; // auto+recycle, weighs 7
  nai_000_030_1 = nai_000_010_1 + nai_000_020_1*p_1; // local, weighs 2
  tmp64 = u_1*usq; // auto, weighs 1
  nai_000_020_3 = tmp3 + nai_000_010_3*p_1 - u_1*(-tmp64 + nai_000_000_4*p_1); // local, weighs 8
  nai_000_010_2 = (nai_000_010_2 - nai_000_010_3)/ab_a - nai_000_020_3*u_1; // auto+recycle, weighs 7
  nai_000_010_3 = nai_000_010_2 + nai_000_020_2*p_1; // local+recycle, weighs 2
  tmp30 = 0.5*nai_000_030_0 - 0.5*nai_000_030_1; // auto, weighs 3
  tmp5 = tmp30/ab_a; // auto, weighs 2
  out[6] = 2.4230585358948*tmp53*(tmp5 + v_0*(nai_000_030_0*v_0 - nai_000_030_1*u_0) - u_0*(nai_000_030_1*v_0 - nai_000_010_3*u_0)); // final, weighs 15
  out[7] = 5.4181235997219*tmp53*(p_2*tmp27 - tmp31*u_2); // final, weighs 6
  out[8] = 5.4181235997219*tmp53*(p_1*tmp43 - nai_200_000_0*u_1); // final, weighs 6
  tmp43 = (nai_200_000_1 - nai_000_001_1)/ab_a - nai_000_002_1*u_2; // auto+recycle, weighs 7
  nai_200_000_0 = tmp43 + nai_000_002_0*p_2; // local+recycle, weighs 2
  tmp27 = (nai_000_001_1 - nai_000_001_2)/ab_a - nai_000_002_2*u_2; // auto+recycle, weighs 7
  tmp31 = tmp27 + nai_000_002_1*p_2; // local+recycle, weighs 2
  usq = u_2*usq; // auto+recycle, weighs 1
  nai_000_002_3 = tmp3 + nai_000_001_3*p_2 - u_2*(-usq + nai_000_000_4*p_2); // local, weighs 8
  nai_000_001_3 = (nai_000_001_2 - nai_000_001_3)/ab_a - nai_000_002_3*u_2; // auto+recycle, weighs 7
  nai_000_003_2 = nai_000_001_3 + nai_000_002_2*p_2; // local, weighs 2
  tmp38 = 0.5*nai_200_000_0 - 0.5*tmp31; // auto, weighs 3
  tmp15 = tmp38/ab_a; // auto, weighs 2
  out[9] = 2.4230585358948*tmp53*(tmp15 + v_0*(nai_200_000_0*v_0 - tmp31*u_0) - u_0*(tmp31*v_0 - nai_000_003_2*u_0)); // final, weighs 15
  nai_010_300_0 = nai_000_300_0*v_1 - nai_000_300_1*u_1; // local, weighs 4
  nai_010_300_1 = nai_000_300_1*v_1 - nai_000_100_3*u_1; // local, weighs 4
  nai_010_200_0 = nai_000_200_0*v_1 - nai_000_200_1*u_1; // local, weighs 4
  nai_010_200_1 = nai_000_200_1*v_1 - nai_000_200_2*u_1; // local, weighs 4
  out[10] = 4.19686049388326*tmp53*(nai_010_300_0*v_0 + (1.5*nai_010_200_0 - 1.5*nai_010_200_1)/ab_a - nai_010_300_1*u_0); // final, weighs 12
  tmp26 = -tmp26 + d_1*v_1; // local+recycle, weighs 3
  tmp73 = -tmp73 + d_2*v_1; // local+recycle, weighs 3
  tmp56 = tmp73*u_0; // auto, weighs 1
  tmp70 = -tmp70 + nai_000_000_2*v_1; // local+recycle, weighs 3
  tmp55 = tmp70*u_0; // auto, weighs 1
  nai_010_100_1 = -tmp55 + p_0*tmp73; // local, weighs 3
  nai_110_200_0 = nai_010_200_0*v_0 + (-nai_010_100_1 - tmp56 + p_0*tmp26)/ab_a - nai_010_200_1*u_0; // local, weighs 12
  tmp67 = -tmp67 + nai_000_000_3*v_1; // local+recycle, weighs 3
  tmp54 = tmp67*u_0; // auto, weighs 1
  nai_010_100_1 = nai_010_200_1*v_0 + (nai_010_100_1 + tmp54 - p_0*tmp70)/ab_a - u_0*(nai_000_200_2*v_1 - nai_000_200_3*u_1); // local+recycle, weighs 15
  nai_000_100_2 = (0.5*nai_000_100_2 - 0.5*nai_000_100_0)/ab_a; // auto+recycle, weighs 5
  out[11] = 9.3844653564063*tmp53*(nai_000_100_2 + nai_110_200_0*p_1 - nai_010_100_1*u_1); // final, weighs 7
  out[12] = 9.3844653564063*tmp53*(nai_110_200_0*p_2 - nai_010_100_1*u_2); // final, weighs 6
  nai_000_100_0 = nai_200_100_1 + nai_000_020_0*v_1; // local+recycle, weighs 2
  nai_200_100_1 = nai_000_010_1 + nai_000_020_1*v_1; // local+recycle, weighs 2
  nai_000_010_1 = nai_000_100_0*v_0 - nai_200_100_1*u_0; // local+recycle, weighs 4
  nai_010_100_1 = nai_200_100_1*v_0 - u_0*(nai_000_010_2 + nai_000_020_2*v_1); // local+recycle, weighs 6
  nai_110_200_0 = (0.5*nai_000_100_0 - 0.5*nai_200_100_1)/ab_a; // auto+recycle, weighs 5
  out[13] = 9.3844653564063*tmp53*(nai_110_200_0 + nai_000_010_1*p_0 - nai_010_100_1*u_0); // final, weighs 7
  nai_000_010_2 = -tmp55 + tmp73*v_0; // local+recycle, weighs 3
  tmp54 = -tmp54 + tmp70*v_0; // local+recycle, weighs 3
  tmp55 = nai_000_010_2*p_2 - tmp54*u_2; // local+recycle, weighs 4
  tmp64 = -tmp64 + nai_000_000_4*v_1; // local+recycle, weighs 3
  nai_010_001_0 = p_2*tmp26 - tmp73*u_2; // local, weighs 4
  nai_010_001_1 = p_2*tmp73 - tmp70*u_2; // local, weighs 4
  nai_010_001_2 = p_2*tmp70 - tmp67*u_2; // local, weighs 4
  out[14] = 16.2543707991657*tmp53*(p_0*(p_1*(p_2*(-tmp56 + tmp26*v_0) - nai_000_010_2*u_2) + (0.5*nai_200_000_2 - 0.5*nai_100_001_1)/ab_a - tmp55*u_1) + (0.5*nai_010_001_0*p_1 + 0.5*nai_010_001_2*u_1 + 0.5*(0.5*nai_200_000_1 - 0.5*nai_000_001_1)/ab_a - 0.5*nai_010_001_1*p_1 - 0.5*nai_010_001_1*u_1 - 0.5*(0.5*nai_000_001_1 - 0.5*nai_000_001_2)/ab_a)/ab_a - u_0*(p_1*tmp55 + (0.5*nai_100_001_1 - 0.5*nai_100_001_2)/ab_a - u_1*(p_2*tmp54 - u_2*(tmp67*v_0 - tmp64*u_0)))); // final, weighs 69
  tmp55 = nai_000_002_0*v_1 - nai_000_002_1*u_1; // local+recycle, weighs 4
  nai_000_010_2 = nai_000_002_1*v_1 - nai_000_002_2*u_1; // local+recycle, weighs 4
  tmp54 = tmp55*v_0 - nai_000_010_2*u_0; // local+recycle, weighs 4
  nai_200_000_1 = nai_000_010_2*v_0 - u_0*(nai_000_002_2*v_1 - nai_000_002_3*u_1); // local+recycle, weighs 8
  out[15] = 9.3844653564063*tmp53*(p_0*tmp54 + (0.5*tmp55 - 0.5*nai_000_010_2)/ab_a - nai_200_000_1*u_0); // final, weighs 12
  nai_000_002_3 = nai_000_030_0*v_1 + (1.5*nai_000_020_0 - 1.5*nai_000_020_1)/ab_a - nai_000_030_1*u_1; // local+recycle, weighs 10
  nai_000_001_1 = nai_000_030_1*v_1 + (1.5*nai_000_020_1 - 1.5*nai_000_020_2)/ab_a - nai_000_010_3*u_1; // local+recycle, weighs 10
  out[16] = 4.19686049388326*tmp53*(nai_000_002_3*v_0 - nai_000_001_1*u_0); // final, weighs 6
  out[17] = 9.3844653564063*tmp53*(nai_000_010_1*p_2 - nai_010_100_1*u_2); // final, weighs 6
  out[18] = 9.3844653564063*tmp53*(p_1*tmp54 + (0.5*nai_100_002_0 - 0.5*nai_100_002_1)/ab_a - nai_200_000_1*u_1); // final, weighs 12
  nai_000_001_2 = nai_200_000_0*v_1 - tmp31*u_1; // local+recycle, weighs 4
  nai_100_002_0 = tmp31*v_1 - nai_000_003_2*u_1; // local+recycle, weighs 4
  out[19] = 4.19686049388326*tmp53*(nai_000_001_2*v_0 - nai_100_002_0*u_0); // final, weighs 6
  nai_100_002_1 = nai_000_300_0*v_2 - nai_000_300_1*u_2; // local+recycle, weighs 4
  nai_000_300_0 = nai_000_300_1*v_2 - nai_000_100_3*u_2; // local+recycle, weighs 4
  nai_200_000_2 = nai_000_200_0*v_2 - nai_000_200_1*u_2; // local+recycle, weighs 4
  nai_100_001_1 = nai_000_200_1*v_2 - nai_000_200_2*u_2; // local+recycle, weighs 4
  out[20] = 4.19686049388326*tmp53*(nai_100_002_1*v_0 + (1.5*nai_200_000_2 - 1.5*nai_100_001_1)/ab_a - nai_000_300_0*u_0); // final, weighs 12
  nai_100_001_2 = -tmp32 + d_1*v_2; // local+recycle, weighs 3
  tmp56 = -d_0 + d_2*v_2; // local+recycle, weighs 3
  nai_000_010_1 = tmp56*u_0; // auto+recycle, weighs 1
  nai_010_100_1 = -tmp69 + nai_000_000_2*v_2; // local+recycle, weighs 3
  tmp54 = nai_010_100_1*u_0; // auto+recycle, weighs 1
  nai_000_100_3 = -tmp54 + p_0*tmp56; // local+recycle, weighs 3
  nai_000_000_2 = nai_200_000_2*v_0 + (-nai_000_010_1 - nai_000_100_3 + nai_100_001_2*p_0)/ab_a - nai_100_001_1*u_0; // local+recycle, weighs 12
  nai_000_300_1 = nai_000_200_2*v_2 - nai_000_200_3*u_2; // local+recycle, weighs 4
  nai_000_200_0 = -tmp66 + nai_000_000_3*v_2; // local+recycle, weighs 3
  nai_000_000_3 = nai_000_200_0*u_0; // auto+recycle, weighs 1
  tmp32 = nai_100_001_1*v_0 + (nai_000_000_3 + nai_000_100_3 - nai_010_100_1*p_0)/ab_a - nai_000_300_1*u_0; // local+recycle, weighs 11
  out[21] = 9.3844653564063*tmp53*(nai_000_000_2*p_1 - tmp32*u_1); // final, weighs 6
  out[22] = 9.3844653564063*tmp53*(nai_000_100_2 + nai_000_000_2*p_2 - tmp32*u_2); // final, weighs 7
  nai_200_000_1 = nai_000_020_0*v_2 - nai_000_020_1*u_2; // local+recycle, weighs 4
  nai_000_100_2 = nai_000_020_1*v_2 - nai_000_020_2*u_2; // local+recycle, weighs 4
  nai_000_200_1 = nai_200_000_1*v_0 - nai_000_100_2*u_0; // local+recycle, weighs 4
  d_0 = nai_000_020_2*v_2 - nai_000_020_3*u_2; // local+recycle, weighs 4
  tmp69 = nai_000_100_2*v_0 - d_0*u_0; // local+recycle, weighs 4
  out[23] = 9.3844653564063*tmp53*(nai_000_200_1*p_0 + (0.5*nai_200_000_1 - 0.5*nai_000_100_2)/ab_a - tmp69*u_0); // final, weighs 12
  tmp66 = -tmp54 + tmp56*v_0; // local+recycle, weighs 3
  nai_000_200_3 = -nai_000_000_3 + nai_010_100_1*v_0; // local+recycle, weighs 3
  d_1 = p_2*tmp66 + (0.5*tmp74 - 0.5*tmp71)/ab_a - nai_000_200_3*u_2; // local+recycle, weighs 10
  d_2 = -usq + nai_000_000_4*v_2; // local+recycle, weighs 3
  nai_000_000_4 = tmp2 - tmp56*u_2; // auto+recycle, weighs 3
  nai_000_020_3 = nai_000_000_4 + nai_100_001_2*p_2; // local+recycle, weighs 2
  nai_000_020_0 = tmp1 - nai_010_100_1*u_2; // auto+recycle, weighs 3
  nai_000_020_1 = nai_000_020_0 + p_2*tmp56; // local+recycle, weighs 2
  tmp54 = tmp0 - nai_000_200_0*u_2; // auto+recycle, weighs 3
  nai_000_020_2 = tmp54 + nai_010_100_1*p_2; // local+recycle, weighs 2
  out[24] = 16.2543707991657*tmp53*(p_0*(p_1*(p_2*(-nai_000_010_1 + nai_100_001_2*v_0) + (0.5*tmp77 - 0.5*tmp74)/ab_a - tmp66*u_2) - d_1*u_1) + (0.5*nai_000_020_2*u_1 + 0.5*nai_000_020_3*p_1 - 0.5*nai_000_020_1*p_1 - 0.5*nai_000_020_1*u_1)/ab_a - u_0*(d_1*p_1 - u_1*(nai_000_200_3*p_2 + (0.5*tmp71 - 0.5*tmp68)/ab_a - u_2*(nai_000_200_0*v_0 - d_2*u_0)))); // final, weighs 55
  nai_000_100_3 = tmp43 + nai_000_002_0*v_2; // local+recycle, weighs 2
  nai_000_200_2 = tmp27 + nai_000_002_1*v_2; // local+recycle, weighs 2
  nai_000_000_2 = nai_000_100_3*v_0 - nai_000_200_2*u_0; // local+recycle, weighs 4
  tmp74 = nai_000_001_3 + nai_000_002_2*v_2; // local+recycle, weighs 2
  tmp43 = nai_000_200_2*v_0 - tmp74*u_0; // local+recycle, weighs 4
  nai_000_000_3 = (0.5*nai_000_100_3 - 0.5*nai_000_200_2)/ab_a; // auto+recycle, weighs 5
  out[25] = 9.3844653564063*tmp53*(nai_000_000_3 + nai_000_000_2*p_0 - tmp43*u_0); // final, weighs 7
  tmp71 = nai_000_030_0*v_2 - nai_000_030_1*u_2; // local+recycle, weighs 4
  tmp27 = nai_000_030_1*v_2 - nai_000_010_3*u_2; // local+recycle, weighs 4
  out[26] = 4.19686049388326*tmp53*(tmp71*v_0 - tmp27*u_0); // final, weighs 6
  out[27] = 9.3844653564063*tmp53*(nai_000_200_1*p_2 + (0.5*nai_100_020_0 - 0.5*nai_100_020_1)/ab_a - tmp69*u_2); // final, weighs 12
  out[28] = 9.3844653564063*tmp53*(nai_000_000_2*p_1 - tmp43*u_1); // final, weighs 6
  tmp32 = nai_200_000_0*v_2 + (1.5*nai_000_002_0 - 1.5*nai_000_002_1)/ab_a - tmp31*u_2; // local+recycle, weighs 10
  nai_000_200_1 = tmp31*v_2 + (1.5*nai_000_002_1 - 1.5*nai_000_002_2)/ab_a - nai_000_003_2*u_2; // local+recycle, weighs 10
  out[29] = 4.19686049388326*tmp53*(tmp32*v_0 - nai_000_200_1*u_0); // final, weighs 6
  usq = nai_000_100_1/ab_a; // auto+recycle, weighs 2
  out[30] = 2.4230585358948*tmp53*(usq + nai_010_300_0*v_1 - nai_010_300_1*u_1); // final, weighs 7
  v_0 = tmp2 - tmp73*u_1; // auto+recycle, weighs 3
  nai_000_002_0 = v_0 + tmp26*v_1; // local+recycle, weighs 2
  tmp69 = tmp1 - tmp70*u_1; // auto+recycle, weighs 3
  nai_000_002_1 = tmp69 + tmp73*v_1; // local+recycle, weighs 2
  tmp66 = tmp0 - tmp67*u_1; // auto+recycle, weighs 3
  nai_000_001_3 = tmp66 + tmp70*v_1; // local+recycle, weighs 2
  nai_000_003_2 = nai_000_002_1*p_0 - nai_000_001_3*u_0; // local+recycle, weighs 4
  nai_010_300_0 = 0.5*nai_000_002_0 - 0.5*nai_000_002_1; // auto+recycle, weighs 3
  nai_010_300_1 = nai_010_300_0/ab_a; // auto+recycle, weighs 2
  nai_000_002_2 = nai_010_300_1 + p_0*(nai_000_002_0*p_0 - nai_000_002_1*u_0) - nai_000_003_2*u_0; // local+recycle, weighs 9
  tmp31 = tmp3 + tmp67*v_1 - tmp64*u_1; // local+recycle, weighs 5
  nai_000_200_3 = 0.5*nai_000_002_1 - 0.5*nai_000_001_3; // auto+recycle, weighs 3
  d_1 = nai_000_200_3/ab_a; // auto+recycle, weighs 2
  tmp77 = d_1 + nai_000_003_2*p_0 - u_0*(nai_000_001_3*p_0 - tmp31*u_0); // local+recycle, weighs 9
  out[31] = 5.4181235997219*tmp53*(nai_000_002_2*p_1 + (nai_010_200_0 - nai_010_200_1)/ab_a - tmp77*u_1); // final, weighs 11
  out[32] = 5.4181235997219*tmp53*(nai_000_002_2*p_2 - tmp77*u_2); // final, weighs 6
  tmp68 = nai_000_002_1*p_1 + (tmp73 - tmp70)/ab_a - nai_000_001_3*u_1; // local+recycle, weighs 9
  nai_000_030_0 = tmp69 + p_1*tmp73; // local+recycle, weighs 2
  nai_000_010_1 = p_1*(nai_000_002_0*p_1 + (tmp26 - tmp73)/ab_a - nai_000_002_1*u_1) + (nai_010_300_0 + v_0 - nai_000_030_0 + p_1*tmp26)/ab_a - tmp68*u_1; // local+recycle, weighs 21
  nai_000_030_1 = p_1*tmp68 + (nai_000_030_0 + nai_000_200_3 - tmp66 - p_1*tmp70)/ab_a - u_1*(nai_000_001_3*p_1 + (tmp70 - tmp67)/ab_a - tmp31*u_1); // local+recycle, weighs 22
  out[33] = 5.4181235997219*tmp53*(nai_000_010_1*p_0 - nai_000_030_1*u_0); // final, weighs 6
  tmp64 = nai_000_002_0*p_2 - nai_000_002_1*u_2; // local+recycle, weighs 4
  nai_100_020_0 = nai_000_002_1*p_2 - nai_000_001_3*u_2; // local+recycle, weighs 4
  nai_000_010_3 = nai_000_001_3*p_2 - tmp31*u_2; // local+recycle, weighs 4
  out[34] = 9.3844653564063*tmp53*(p_0*(p_1*tmp64 + (nai_010_001_0 - nai_010_001_1)/ab_a - nai_100_020_0*u_1) - u_0*(nai_100_020_0*p_1 + (nai_010_001_1 - nai_010_001_2)/ab_a - nai_000_010_3*u_1)); // final, weighs 24
  nai_100_020_1 = nai_010_300_1 + p_2*tmp64 - nai_100_020_0*u_2; // local+recycle, weighs 5
  tmp0 = d_1 + nai_100_020_0*p_2 - nai_000_010_3*u_2; // local+recycle, weighs 5
  out[35] = 5.4181235997219*tmp53*(nai_100_020_1*p_0 - tmp0*u_0); // final, weighs 6
  out[36] = 2.4230585358948*tmp53*(nai_000_002_3*v_1 + (tmp30 + 1.5*nai_000_100_0 - 1.5*nai_200_100_1)/ab_a - nai_000_001_1*u_1); // final, weighs 13
  out[37] = 5.4181235997219*tmp53*(nai_000_010_1*p_2 - nai_000_030_1*u_2); // final, weighs 6
  out[38] = 5.4181235997219*tmp53*(nai_100_020_1*p_1 + (tmp55 - nai_000_010_2)/ab_a - tmp0*u_1); // final, weighs 11
  out[39] = 2.4230585358948*tmp53*(tmp15 + nai_000_001_2*v_1 - nai_100_002_0*u_1); // final, weighs 7
  out[40] = 4.19686049388326*tmp53*(nai_100_002_1*v_1 - nai_000_300_0*u_1); // final, weighs 6
  nai_000_000_2 = nai_200_000_2*v_1 - nai_100_001_1*u_1; // local+recycle, weighs 4
  tmp43 = nai_100_001_1*v_1 - nai_000_300_1*u_1; // local+recycle, weighs 4
  out[41] = 9.3844653564063*tmp53*(nai_000_000_2*p_1 + (0.5*nai_200_000_2 - 0.5*nai_100_001_1)/ab_a - tmp43*u_1); // final, weighs 12
  out[42] = 9.3844653564063*tmp53*(nai_000_000_2*p_2 + (0.5*nai_010_200_0 - 0.5*nai_010_200_1)/ab_a - tmp43*u_2); // final, weighs 12
  nai_200_000_0 = tmp56*u_1; // auto+recycle, weighs 1
  nai_000_100_1 = nai_010_100_1*u_1; // auto+recycle, weighs 1
  nai_000_300_1 = -nai_000_100_1 + p_1*tmp56; // local+recycle, weighs 3
  tmp2 = nai_200_000_1*v_1 + (-nai_000_300_1 - nai_200_000_0 + nai_100_001_2*p_1)/ab_a - nai_000_100_2*u_1; // local+recycle, weighs 12
  tmp1 = nai_000_200_0*u_1; // auto+recycle, weighs 1
  nai_000_002_3 = nai_000_100_2*v_1 + (nai_000_300_1 + tmp1 - nai_010_100_1*p_1)/ab_a - d_0*u_1; // local+recycle, weighs 11
  out[43] = 9.3844653564063*tmp53*(p_0*tmp2 - nai_000_002_3*u_0); // final, weighs 6
  v_0 = -nai_000_100_1 + tmp56*v_1; // local+recycle, weighs 3
  d_0 = -tmp1 + nai_010_100_1*v_1; // local+recycle, weighs 3
  nai_000_001_1 = p_2*v_0 + (0.5*tmp73 - 0.5*tmp70)/ab_a - d_0*u_2; // local+recycle, weighs 10
  out[44] = 16.2543707991657*tmp53*(p_0*(p_1*(p_2*(-nai_200_000_0 + nai_100_001_2*v_1) + (0.5*tmp26 - 0.5*tmp73)/ab_a - u_2*v_0) + (0.5*nai_000_020_3 - 0.5*nai_000_020_1)/ab_a - nai_000_001_1*u_1) - u_0*(nai_000_001_1*p_1 + (0.5*nai_000_020_1 - 0.5*nai_000_020_2)/ab_a - u_1*(d_0*p_2 + (0.5*tmp70 - 0.5*tmp67)/ab_a - u_2*(nai_000_200_0*v_1 - d_2*u_1)))); // final, weighs 53
  nai_000_002_0 = nai_000_100_3*v_1 - nai_000_200_2*u_1; // local+recycle, weighs 4
  tmp69 = nai_000_200_2*v_1 - tmp74*u_1; // local+recycle, weighs 4
  out[45] = 9.3844653564063*tmp53*(nai_000_002_0*p_0 - tmp69*u_0); // final, weighs 6
  out[46] = 4.19686049388326*tmp53*(tmp71*v_1 + (1.5*nai_200_000_1 - 1.5*nai_000_100_2)/ab_a - tmp27*u_1); // final, weighs 12
  out[47] = 9.3844653564063*tmp53*(nai_110_200_0 + p_2*tmp2 - nai_000_002_3*u_2); // final, weighs 7
  out[48] = 9.3844653564063*tmp53*(nai_000_000_3 + nai_000_002_0*p_1 - tmp69*u_1); // final, weighs 7
  out[49] = 4.19686049388326*tmp53*(tmp32*v_1 - nai_000_200_1*u_1); // final, weighs 6
  out[50] = 2.4230585358948*tmp53*(usq + nai_100_002_1*v_2 - nai_000_300_0*u_2); // final, weighs 7
  nai_000_001_2 = nai_000_000_4 + nai_100_001_2*v_2; // local+recycle, weighs 2
  nai_000_002_1 = nai_000_020_0 + tmp56*v_2; // local+recycle, weighs 2
  nai_100_002_0 = tmp54 + nai_010_100_1*v_2; // local+recycle, weighs 2
  tmp66 = nai_000_002_1*p_0 - nai_100_002_0*u_0; // local+recycle, weighs 4
  nai_000_001_3 = 0.5*nai_000_001_2 - 0.5*nai_000_002_1; // auto+recycle, weighs 3
  nai_000_003_2 = nai_000_001_3/ab_a; // auto+recycle, weighs 2
  tmp15 = nai_000_003_2 + p_0*(nai_000_001_2*p_0 - nai_000_002_1*u_0) - tmp66*u_0; // local+recycle, weighs 9
  nai_010_300_0 = tmp3 + nai_000_200_0*v_2 - d_2*u_2; // local+recycle, weighs 5
  nai_010_300_1 = 0.5*nai_000_002_1 - 0.5*nai_100_002_0; // auto+recycle, weighs 3
  nai_010_200_0 = nai_010_300_1/ab_a; // auto+recycle, weighs 2
  nai_010_200_1 = nai_010_200_0 + p_0*tmp66 - u_0*(nai_100_002_0*p_0 - nai_010_300_0*u_0); // local+recycle, weighs 9
  out[51] = 5.4181235997219*tmp53*(p_1*tmp15 - nai_010_200_1*u_1); // final, weighs 6
  out[52] = 5.4181235997219*tmp53*(p_2*tmp15 + (nai_200_000_2 - nai_100_001_1)/ab_a - nai_010_200_1*u_2); // final, weighs 11
  nai_000_002_2 = nai_000_002_1*p_1 - nai_100_002_0*u_1; // local+recycle, weighs 4
  nai_100_002_1 = nai_000_003_2 + p_1*(nai_000_001_2*p_1 - nai_000_002_1*u_1) - nai_000_002_2*u_1; // local+recycle, weighs 9
  tmp31 = nai_010_200_0 + nai_000_002_2*p_1 - u_1*(nai_100_002_0*p_1 - nai_010_300_0*u_1); // local+recycle, weighs 9
  out[53] = 5.4181235997219*tmp53*(nai_100_002_1*p_0 - tmp31*u_0); // final, weighs 6
  tmp3 = nai_000_001_2*p_2 + (nai_100_001_2 - tmp56)/ab_a - nai_000_002_1*u_2; // local+recycle, weighs 9
  nai_000_200_3 = nai_000_002_1*p_2 + (tmp56 - nai_010_100_1)/ab_a - nai_100_002_0*u_2; // local+recycle, weighs 9
  d_1 = nai_100_002_0*p_2 + (nai_010_100_1 - nai_000_200_0)/ab_a - nai_010_300_0*u_2; // local+recycle, weighs 9
  out[54] = 9.3844653564063*tmp53*(p_0*(p_1*tmp3 - nai_000_200_3*u_1) - u_0*(nai_000_200_3*p_1 - d_1*u_1)); // final, weighs 14
  d_2 = p_2*tmp3 + (nai_000_001_3 + nai_000_020_3 - nai_000_020_1)/ab_a - nai_000_200_3*u_2; // local+recycle, weighs 10
  tmp77 = nai_000_200_3*p_2 + (nai_000_020_1 + nai_010_300_1 - nai_000_020_2)/ab_a - d_1*u_2; // local+recycle, weighs 10
  out[55] = 5.4181235997219*tmp53*(d_2*p_0 - tmp77*u_0); // final, weighs 6
  out[56] = 2.4230585358948*tmp53*(tmp5 + tmp71*v_2 - tmp27*u_2); // final, weighs 7
  out[57] = 5.4181235997219*tmp53*(nai_100_002_1*p_2 + (nai_200_000_1 - nai_000_100_2)/ab_a - tmp31*u_2); // final, weighs 11
  out[58] = 5.4181235997219*tmp53*(d_2*p_1 - tmp77*u_1); // final, weighs 6
  out[59] = 2.4230585358948*tmp53*(tmp32*v_2 + (tmp38 + 1.5*nai_000_100_3 - 1.5*nai_000_200_2)/ab_a - nai_000_200_1*u_2); // final, weighs 13
  // total weight = 1832
}

void gint2_nai_cF_cF(double* a, double a_a, double* b, double b_a, double* c, double* out)
{
  // Generated code based on 'Common SubExpression' analysis.
  // Number of local variables: 137
  double ab_a, d_0, d_1, d_2, nai_000_000_2, nai_000_000_3, nai_000_000_4, nai_000_000_5, nai_000_001_1, nai_000_001_2, nai_000_001_3, nai_000_001_4, nai_000_002_0, nai_000_002_1, nai_000_002_2, nai_000_002_3, nai_000_010_1, nai_000_010_2, nai_000_010_3, nai_000_010_4, nai_000_020_0, nai_000_020_1, nai_000_020_2, nai_000_020_3, nai_000_020_4, nai_010_001_1, nai_010_002_1, nai_010_020_0, nai_010_020_1, nai_010_020_2, nai_020_001_0, nai_020_001_1, nai_020_001_2, nai_020_010_2, nai_020_200_2, nai_100_000_0, nai_100_000_1, nai_100_000_2, nai_100_000_3, nai_100_000_4, nai_100_001_1, nai_100_030_1, nai_110_000_0, nai_110_000_1, nai_110_000_2, nai_110_000_3, nai_110_001_0, nai_110_001_1, nai_110_001_2, nai_110_020_0, nai_110_020_1, nai_120_001_1, nai_200_000_0, nai_200_000_1, nai_200_000_2, nai_200_000_3, nai_200_000_4, nai_200_001_0, nai_200_001_1, nai_200_001_2, nai_200_002_1, nai_200_003_0, nai_200_010_0, nai_200_010_1, nai_200_010_2, nai_200_020_0, nai_200_020_1, nai_200_030_0, nai_210_001_1, nai_300_000_0, nai_300_000_1, nai_300_000_2, nai_300_000_3, nai_300_100_0, nai_300_100_1, nai_300_200_0, nai_300_200_1, p_0, p_1, p_2, tmp0, tmp1, tmp106, tmp107, tmp108, tmp11, tmp110, tmp115, tmp12, tmp121, tmp122, tmp123, tmp124, tmp126, tmp129, tmp131, tmp133, tmp138, tmp139, tmp140, tmp157, tmp159, tmp160, tmp161, tmp162, tmp163, tmp164, tmp165, tmp166, tmp167, tmp169, tmp170, tmp173, tmp2, tmp20, tmp21, tmp3, tmp33, tmp34, tmp4, tmp62, tmp63, tmp69, tmp70, tmp74, tmp77, tmp8, tmp80, tmp87, tmp88, u_0, u_1, u_2, usq, v_0, v_1, v_2;
  ab_a = a_a + b_a; // local, weighs 1
  p_0 = (a[0]*a_a + b[0]*b_a)/ab_a; // local, weighs 5
  p_1 = (a[1]*a_a + b[1]*b_a)/ab_a; // local, weighs 5
  p_2 = (a[2]*a_a + b[2]*b_a)/ab_a; // local, weighs 5
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
  usq = 6.28318530717959*d_0*gaux(usq, 6); // local+recycle, weighs 3
  d_0 = u_0*usq; // auto+recycle, weighs 1
  tmp4 = (0.5*nai_000_000_4 - 0.5*nai_000_000_5)/ab_a; // auto, weighs 5
  nai_200_000_4 = tmp4 + nai_100_000_4*v_0 - u_0*(-d_0 + nai_000_000_5*v_0); // local, weighs 8
  nai_100_000_4 = (nai_100_000_3 - nai_100_000_4)/ab_a - nai_200_000_4*u_0; // auto+recycle, weighs 7
  nai_300_000_3 = nai_100_000_4 + nai_200_000_3*v_0; // local, weighs 2
  tmp108 = tmp108 + nai_200_000_2*p_0; // local+recycle, weighs 2
  tmp69 = 0.5*nai_300_000_1 - 0.5*nai_300_000_2; // auto, weighs 3
  nai_300_200_1 = nai_300_100_1*p_0 + (tmp69 + 1.5*tmp107 - 1.5*tmp108)/ab_a - u_0*(nai_300_000_2*p_0 + (1.5*nai_200_000_2 - 1.5*nai_200_000_3)/ab_a - nai_300_000_3*u_0); // local, weighs 21
  tmp63 = tmp63 + nai_100_000_1*p_0; // local+recycle, weighs 2
  tmp77 = 0.5*nai_200_000_0 - 0.5*nai_200_000_1; // auto, weighs 3
  tmp106 = p_0*tmp106 + (tmp110 + tmp77 - tmp63 + nai_100_000_0*p_0)/ab_a - tmp107*u_0; // local+recycle, weighs 12
  tmp62 = tmp62 + nai_100_000_2*p_0; // local+recycle, weighs 2
  tmp110 = 0.5*nai_200_000_1 - 0.5*nai_200_000_2; // auto+recycle, weighs 3
  tmp63 = p_0*tmp107 + (tmp110 + tmp63 - tmp62)/ab_a - tmp108*u_0; // local+recycle, weighs 10
  tmp107 = pow(a_a,2.25)*pow(b_a,2.25); // auto+recycle, weighs 5
  out[0] = 2.16724943988876*tmp107*(nai_300_200_0*p_0 + (nai_300_100_0 - nai_300_100_1 + 1.5*tmp106 - 1.5*tmp63)/ab_a - nai_300_200_1*u_0); // final, weighs 15
  out[1] = 4.84611707178961*tmp107*(nai_300_200_0*p_1 - nai_300_200_1*u_1); // final, weighs 6
  out[2] = 4.84611707178961*tmp107*(nai_300_200_0*p_2 - nai_300_200_1*u_2); // final, weighs 6
  nai_300_200_0 = nai_300_000_0*p_1 - nai_300_000_1*u_1; // local+recycle, weighs 4
  nai_300_100_0 = nai_300_000_1*p_1 - nai_300_000_2*u_1; // local+recycle, weighs 4
  nai_300_200_1 = tmp70/ab_a; // auto+recycle, weighs 2
  nai_300_100_1 = nai_300_200_1 + nai_300_200_0*p_1 - nai_300_100_0*u_1; // local+recycle, weighs 5
  tmp70 = tmp69/ab_a; // auto+recycle, weighs 2
  tmp69 = tmp70 + nai_300_100_0*p_1 - u_1*(nai_300_000_2*p_1 - nai_300_000_3*u_1); // local+recycle, weighs 9
  tmp126 = nai_200_000_1*u_1; // auto, weighs 1
  nai_200_010_0 = -tmp126 + nai_200_000_0*p_1; // local, weighs 3
  tmp124 = nai_200_000_2*u_1; // auto, weighs 1
  nai_200_010_1 = -tmp124 + nai_200_000_1*p_1; // local, weighs 3
  tmp77 = tmp77/ab_a; // auto+recycle, weighs 2
  nai_200_020_0 = tmp77 + nai_200_010_0*p_1 - nai_200_010_1*u_1; // local, weighs 5
  tmp122 = nai_200_000_3*u_1; // auto, weighs 1
  nai_200_010_2 = -tmp122 + nai_200_000_2*p_1; // local, weighs 3
  tmp110 = tmp110/ab_a; // auto+recycle, weighs 2
  nai_200_020_1 = tmp110 + nai_200_010_1*p_1 - nai_200_010_2*u_1; // local, weighs 5
  tmp11 = (1.5*nai_200_020_0 - 1.5*nai_200_020_1)/ab_a; // auto, weighs 5
  out[3] = 4.84611707178961*tmp107*(tmp11 + nai_300_100_1*p_0 - tmp69*u_0); // final, weighs 7
  nai_300_000_0 = nai_300_000_0*p_2 - nai_300_000_1*u_2; // local+recycle, weighs 4
  nai_300_000_1 = nai_300_000_1*p_2 - nai_300_000_2*u_2; // local+recycle, weighs 4
  nai_300_000_3 = nai_300_000_2*p_2 - nai_300_000_3*u_2; // local+recycle, weighs 4
  nai_300_000_2 = nai_200_000_1*u_2; // auto+recycle, weighs 1
  nai_200_001_0 = -nai_300_000_2 + nai_200_000_0*p_2; // local, weighs 3
  tmp123 = nai_200_000_2*u_2; // auto, weighs 1
  nai_200_001_1 = -tmp123 + nai_200_000_1*p_2; // local, weighs 3
  tmp121 = nai_200_000_3*u_2; // auto, weighs 1
  nai_200_001_2 = -tmp121 + nai_200_000_2*p_2; // local, weighs 3
  out[4] = 8.39372098776651*tmp107*(p_0*(nai_300_000_0*p_1 - nai_300_000_1*u_1) + (1.5*nai_200_001_0*p_1 + 1.5*nai_200_001_2*u_1 - 1.5*nai_200_001_1*p_1 - 1.5*nai_200_001_1*u_1)/ab_a - u_0*(nai_300_000_1*p_1 - nai_300_000_3*u_1)); // final, weighs 28
  nai_300_200_1 = nai_300_200_1 + nai_300_000_0*p_2 - nai_300_000_1*u_2; // local+recycle, weighs 5
  tmp70 = tmp70 + nai_300_000_1*p_2 - nai_300_000_3*u_2; // local+recycle, weighs 5
  nai_300_000_3 = tmp77 + nai_200_001_0*p_2 - nai_200_001_1*u_2; // local+recycle, weighs 5
  nai_200_002_1 = tmp110 + nai_200_001_1*p_2 - nai_200_001_2*u_2; // local, weighs 5
  tmp12 = (1.5*nai_300_000_3 - 1.5*nai_200_002_1)/ab_a; // auto, weighs 5
  out[5] = 4.84611707178961*tmp107*(tmp12 + nai_300_200_1*p_0 - tmp70*u_0); // final, weighs 7
  out[6] = 2.16724943988876*tmp107*(nai_300_100_1*p_1 + (nai_300_200_0 - nai_300_100_0)/ab_a - tmp69*u_1); // final, weighs 11
  out[7] = 4.84611707178961*tmp107*(nai_300_100_1*p_2 - tmp69*u_2); // final, weighs 6
  out[8] = 4.84611707178961*tmp107*(nai_300_200_1*p_1 - tmp70*u_1); // final, weighs 6
  out[9] = 2.16724943988876*tmp107*(nai_300_200_1*p_2 + (nai_300_000_0 - nai_300_000_1)/ab_a - tmp70*u_2); // final, weighs 11
  tmp69 = -tmp173 + d_1*p_0; // local+recycle, weighs 3
  nai_300_200_1 = -tmp170 + d_2*p_0; // local+recycle, weighs 3
  nai_300_100_1 = tmp3 + p_0*tmp69 - nai_300_200_1*u_0; // local+recycle, weighs 5
  tmp173 = -tmp167 + nai_000_000_2*p_0; // local+recycle, weighs 3
  tmp170 = tmp2 + nai_300_200_1*p_0 - tmp173*u_0; // local+recycle, weighs 5
  tmp167 = (tmp69 - nai_300_200_1)/ab_a - tmp170*u_0; // auto+recycle, weighs 7
  nai_300_000_0 = tmp167 + nai_300_100_1*p_0; // local+recycle, weighs 2
  tmp164 = -tmp164 + nai_000_000_3*p_0; // local+recycle, weighs 3
  tmp70 = tmp1 + p_0*tmp173 - tmp164*u_0; // local+recycle, weighs 5
  nai_300_200_0 = (nai_300_200_1 - tmp173)/ab_a - tmp70*u_0; // auto+recycle, weighs 7
  nai_300_000_1 = nai_300_200_0 + p_0*tmp170; // local+recycle, weighs 2
  nai_300_100_0 = -tmp161 + nai_000_000_4*p_0; // local+recycle, weighs 3
  tmp161 = tmp0 + p_0*tmp164 - nai_300_100_0*u_0; // local+recycle, weighs 5
  tmp69 = (tmp173 - tmp164)/ab_a - tmp161*u_0; // auto+recycle, weighs 7
  nai_300_200_1 = tmp69 + p_0*tmp70; // local+recycle, weighs 2
  tmp173 = nai_300_000_1*v_0 + (1.5*tmp170 - 1.5*tmp70)/ab_a - nai_300_200_1*u_0; // local+recycle, weighs 10
  nai_300_200_0 = nai_300_200_0 + tmp170*v_0; // local+recycle, weighs 2
  tmp87 = 0.5*nai_300_000_0 - 0.5*nai_300_000_1; // auto, weighs 3
  tmp167 = v_0*(nai_300_000_0*v_0 + (1.5*nai_300_100_1 - 1.5*tmp170)/ab_a - nai_300_000_1*u_0) + (tmp87 + 1.5*tmp167 - 1.5*nai_300_200_0 + 1.5*nai_300_100_1*v_0)/ab_a - tmp173*u_0; // local+recycle, weighs 24
  d_0 = tmp4 + nai_300_100_0*p_0 - u_0*(-d_0 + nai_000_000_5*p_0); // local+recycle, weighs 8
  tmp164 = p_0*tmp161 + (tmp164 - nai_300_100_0)/ab_a - d_0*u_0; // local+recycle, weighs 9
  nai_300_100_0 = 0.5*nai_300_000_1 - 0.5*nai_300_200_1; // auto+recycle, weighs 3
  tmp69 = tmp173*v_0 + (nai_300_100_0 + 1.5*nai_300_200_0 - 1.5*tmp69 - 1.5*tmp70*v_0)/ab_a - u_0*(nai_300_200_1*v_0 + (1.5*tmp70 - 1.5*tmp161)/ab_a - tmp164*u_0); // local+recycle, weighs 24
  out[10] = 4.84611707178961*tmp107*(tmp167*v_1 - tmp69*u_1); // final, weighs 6
  tmp173 = tmp106*v_1 - tmp63*u_1; // local+recycle, weighs 4
  nai_300_200_0 = 0.5*nai_200_000_2 - 0.5*nai_200_000_3; // auto+recycle, weighs 3
  nai_100_000_4 = p_0*tmp108 + (nai_300_200_0 + tmp62 - tmp115 - nai_100_000_3*p_0)/ab_a - u_0*(nai_100_000_4 + nai_200_000_3*p_0); // local+recycle, weighs 15
  tmp115 = tmp63*v_1 - nai_100_000_4*u_1; // local+recycle, weighs 4
  tmp108 = (0.5*tmp106 - 0.5*tmp63)/ab_a; // auto+recycle, weighs 5
  out[11] = 10.8362471994438*tmp107*(tmp108 + p_1*tmp173 - tmp115*u_1); // final, weighs 7
  out[12] = 10.8362471994438*tmp107*(p_2*tmp173 - tmp115*u_2); // final, weighs 6
  nai_200_010_0 = nai_200_020_0*v_1 + (nai_200_010_0 - nai_200_010_1)/ab_a - nai_200_020_1*u_1; // local+recycle, weighs 9
  tmp173 = nai_200_000_4*u_1; // auto+recycle, weighs 1
  tmp62 = nai_300_200_0/ab_a; // auto+recycle, weighs 2
  nai_300_200_0 = tmp62 + nai_200_010_2*p_1 - u_1*(-tmp173 + nai_200_000_3*p_1); // local+recycle, weighs 8
  tmp115 = nai_200_020_1*v_1 + (nai_200_010_1 - nai_200_010_2)/ab_a - nai_300_200_0*u_1; // local+recycle, weighs 9
  nai_200_010_1 = d_2*u_1; // auto+recycle, weighs 1
  nai_200_010_2 = -nai_200_010_1 + d_1*p_1; // local+recycle, weighs 3
  tmp169 = nai_000_000_2*u_1; // auto, weighs 1
  nai_000_010_1 = -tmp169 + d_2*p_1; // local, weighs 3
  nai_000_020_0 = tmp3 + nai_200_010_2*p_1 - nai_000_010_1*u_1; // local, weighs 5
  tmp166 = nai_000_000_3*u_1; // auto, weighs 1
  nai_000_010_2 = -tmp166 + nai_000_000_2*p_1; // local, weighs 3
  nai_000_020_1 = tmp2 + nai_000_010_1*p_1 - nai_000_010_2*u_1; // local, weighs 5
  nai_200_010_2 = (nai_200_010_2 - nai_000_010_1)/ab_a - nai_000_020_1*u_1; // auto+recycle, weighs 7
  nai_010_020_0 = nai_200_010_2 + nai_000_020_0*v_1; // local, weighs 2
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
  out[13] = 10.8362471994438*tmp107*(nai_200_010_0*p_0 + (nai_110_020_0 - nai_110_020_1)/ab_a - tmp115*u_0); // final, weighs 11
  tmp124 = -tmp124 + nai_200_000_1*v_1; // local+recycle, weighs 3
  tmp122 = -tmp122 + nai_200_000_2*v_1; // local+recycle, weighs 3
  nai_210_001_1 = p_2*tmp124 - tmp122*u_2; // local, weighs 4
  nai_200_010_1 = -nai_200_010_1 + d_1*v_1; // local+recycle, weighs 3
  tmp169 = -tmp169 + d_2*v_1; // local+recycle, weighs 3
  tmp140 = tmp169*u_0; // auto, weighs 1
  nai_110_000_0 = -tmp140 + nai_200_010_1*v_0; // local, weighs 3
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
  out[14] = 18.7689307128126*tmp107*(p_0*(p_1*(p_2*(-tmp126 + nai_200_000_0*v_1) - tmp124*u_2) + (0.5*nai_200_001_0 - 0.5*nai_200_001_1)/ab_a - nai_210_001_1*u_1) + (nai_110_001_0*p_1 + nai_110_001_2*u_1 + (-0.5*nai_100_001_1 + 0.5*nai_100_000_0*p_2 - 0.5*nai_100_000_1*u_2)/ab_a - nai_110_001_1*p_1 - nai_110_001_1*u_1 - (0.5*nai_100_001_1 + 0.5*nai_100_000_3*u_2 - 0.5*nai_100_000_2*p_2)/ab_a)/ab_a - u_0*(nai_210_001_1*p_1 + (0.5*nai_200_001_1 - 0.5*nai_200_001_2)/ab_a - u_1*(p_2*tmp122 - u_2*(-tmp173 + nai_200_000_3*v_1)))); // final, weighs 73
  tmp126 = nai_300_000_3*v_1 - nai_200_002_1*u_1; // local+recycle, weighs 4
  tmp124 = nai_200_000_4*u_2; // auto+recycle, weighs 1
  nai_200_000_4 = tmp62 + nai_200_001_2*p_2 - u_2*(-tmp124 + nai_200_000_3*p_2); // local+recycle, weighs 8
  tmp173 = nai_200_002_1*v_1 - nai_200_000_4*u_1; // local+recycle, weighs 4
  tmp122 = d_2*u_2; // auto+recycle, weighs 1
  nai_210_001_1 = -tmp122 + d_1*p_2; // local+recycle, weighs 3
  nai_100_001_1 = nai_000_000_2*u_2; // auto+recycle, weighs 1
  nai_000_001_1 = -nai_100_001_1 + d_2*p_2; // local, weighs 3
  nai_000_002_0 = tmp3 + nai_210_001_1*p_2 - nai_000_001_1*u_2; // local, weighs 5
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
  out[15] = 10.8362471994438*tmp107*(nai_010_002_1 + p_0*tmp126 - tmp173*u_0); // final, weighs 7
  nai_200_010_2 = nai_200_010_2 + nai_000_020_0*p_1; // local+recycle, weighs 2
  nai_000_010_1 = nai_000_010_1 + nai_000_020_1*p_1; // local+recycle, weighs 2
  nai_000_010_2 = nai_000_010_2 + nai_000_020_2*p_1; // local+recycle, weighs 2
  nai_100_030_1 = nai_000_010_1*v_0 - nai_000_010_2*u_0; // local, weighs 4
  tmp74 = 0.5*nai_200_010_2 - 0.5*nai_000_010_1; // auto, weighs 3
  tmp21 = tmp74/ab_a; // auto, weighs 2
  nai_200_030_0 = tmp21 + v_0*(nai_200_010_2*v_0 - nai_000_010_1*u_0) - nai_100_030_1*u_0; // local, weighs 9
  tmp157 = u_1*usq; // auto, weighs 1
  nai_000_020_4 = tmp4 + nai_000_010_4*p_1 - u_1*(-tmp157 + nai_000_000_5*p_1); // local, weighs 8
  nai_000_010_3 = nai_000_020_3*p_1 + (nai_000_010_3 - nai_000_010_4)/ab_a - nai_000_020_4*u_1; // local+recycle, weighs 9
  nai_000_010_4 = 0.5*nai_000_010_1 - 0.5*nai_000_010_2; // auto+recycle, weighs 3
  tmp20 = nai_000_010_4/ab_a; // auto, weighs 2
  nai_100_030_1 = tmp20 + nai_100_030_1*v_0 - u_0*(nai_000_010_2*v_0 - nai_000_010_3*u_0); // local+recycle, weighs 9
  out[16] = 4.84611707178961*tmp107*(tmp11 + nai_200_030_0*v_1 - nai_100_030_1*u_1); // final, weighs 7
  out[17] = 10.8362471994438*tmp107*(nai_200_010_0*p_2 - tmp115*u_2); // final, weighs 6
  out[18] = 10.8362471994438*tmp107*(p_1*tmp126 + (0.5*nai_300_000_3 - 0.5*nai_200_002_1)/ab_a - tmp173*u_1); // final, weighs 12
  tmp115 = (nai_210_001_1 - nai_000_001_1)/ab_a - nai_000_002_1*u_2; // auto+recycle, weighs 7
  tmp126 = tmp115 + nai_000_002_0*p_2; // local+recycle, weighs 2
  nai_200_010_0 = (nai_000_001_1 - nai_000_001_2)/ab_a - nai_000_002_2*u_2; // auto+recycle, weighs 7
  tmp173 = nai_200_010_0 + nai_000_002_1*p_2; // local+recycle, weighs 2
  tmp11 = (nai_000_001_2 - nai_000_001_3)/ab_a - nai_000_002_3*u_2; // auto+recycle, weighs 7
  nai_210_001_1 = tmp11 + nai_000_002_2*p_2; // local+recycle, weighs 2
  nai_000_001_1 = tmp173*v_0 - nai_210_001_1*u_0; // local+recycle, weighs 4
  nai_000_001_2 = 0.5*tmp126 - 0.5*tmp173; // auto+recycle, weighs 3
  tmp34 = nai_000_001_2/ab_a; // auto, weighs 2
  nai_200_003_0 = tmp34 + v_0*(tmp126*v_0 - tmp173*u_0) - nai_000_001_1*u_0; // local, weighs 9
  usq = u_2*usq; // auto+recycle, weighs 1
  nai_000_001_3 = (nai_000_001_3 - nai_000_001_4)/ab_a - u_2*(tmp4 + nai_000_001_4*p_2 - u_2*(-usq + nai_000_000_5*p_2)); // auto+recycle, weighs 15
  nai_000_001_4 = nai_000_001_3 + nai_000_002_3*p_2; // local+recycle, weighs 2
  tmp88 = 0.5*tmp173 - 0.5*nai_210_001_1; // auto, weighs 3
  tmp33 = tmp88/ab_a; // auto, weighs 2
  nai_000_001_1 = tmp33 + nai_000_001_1*v_0 - u_0*(nai_210_001_1*v_0 - nai_000_001_4*u_0); // local+recycle, weighs 9
  out[19] = 4.84611707178961*tmp107*(nai_200_003_0*v_1 - nai_000_001_1*u_1); // final, weighs 6
  out[20] = 4.84611707178961*tmp107*(tmp167*v_2 - tmp69*u_2); // final, weighs 6
  tmp106 = tmp106*v_2 - tmp63*u_2; // local+recycle, weighs 4
  nai_100_000_4 = tmp63*v_2 - nai_100_000_4*u_2; // local+recycle, weighs 4
  out[21] = 10.8362471994438*tmp107*(p_1*tmp106 - nai_100_000_4*u_1); // final, weighs 6
  out[22] = 10.8362471994438*tmp107*(tmp108 + p_2*tmp106 - nai_100_000_4*u_2); // final, weighs 7
  tmp108 = nai_200_020_0*v_2 - nai_200_020_1*u_2; // local+recycle, weighs 4
  tmp69 = nai_200_020_1*v_2 - nai_300_200_0*u_2; // local+recycle, weighs 4
  tmp167 = nai_000_020_0*v_2 - nai_000_020_1*u_2; // local+recycle, weighs 4
  tmp63 = nai_000_020_1*v_2 - nai_000_020_2*u_2; // local+recycle, weighs 4
  tmp106 = nai_000_020_2*v_2 - nai_000_020_3*u_2; // local+recycle, weighs 4
  nai_300_200_0 = (tmp106*u_0 + tmp167*v_0 - tmp63*u_0 - tmp63*v_0)/ab_a; // auto+recycle, weighs 11
  out[23] = 10.8362471994438*tmp107*(nai_300_200_0 + p_0*tmp108 - tmp69*u_0); // final, weighs 7
  nai_100_000_4 = -tmp123 + nai_200_000_1*v_2; // local+recycle, weighs 3
  tmp123 = -tmp121 + nai_200_000_2*v_2; // local+recycle, weighs 3
  tmp121 = tmp110 + nai_100_000_4*p_2 - tmp123*u_2; // local+recycle, weighs 5
  d_1 = -tmp122 + d_1*v_2; // local+recycle, weighs 3
  d_2 = -nai_100_001_1 + d_2*v_2; // local+recycle, weighs 3
  tmp110 = d_2*u_0; // auto+recycle, weighs 1
  tmp122 = -tmp110 + d_1*v_0; // local+recycle, weighs 3
  nai_100_001_1 = -tmp165 + nai_000_000_2*v_2; // local+recycle, weighs 3
  tmp165 = nai_100_001_1*u_0; // auto+recycle, weighs 1
  nai_200_000_1 = -tmp165 + d_2*v_0; // local+recycle, weighs 3
  nai_200_000_2 = p_2*tmp122 + (0.5*nai_100_000_0 - 0.5*nai_100_000_1)/ab_a - nai_200_000_1*u_2; // local+recycle, weighs 10
  nai_100_000_0 = -tmp162 + nai_000_000_3*v_2; // local+recycle, weighs 3
  nai_000_000_2 = nai_100_000_0*u_0; // auto+recycle, weighs 1
  nai_000_000_3 = -nai_000_000_2 + nai_100_001_1*v_0; // local+recycle, weighs 3
  tmp162 = nai_200_000_1*p_2 + (0.5*nai_100_000_1 - 0.5*nai_100_000_2)/ab_a - nai_000_000_3*u_2; // local+recycle, weighs 10
  tmp159 = -tmp159 + nai_000_000_4*v_2; // local+recycle, weighs 3
  nai_000_000_4 = nai_100_000_0*v_0 - tmp159*u_0; // local+recycle, weighs 4
  nai_100_000_3 = nai_000_000_3*p_2 + (0.5*nai_100_000_2 - 0.5*nai_100_000_3)/ab_a - nai_000_000_4*u_2; // local+recycle, weighs 10
  out[24] = 18.7689307128126*tmp107*(p_0*(p_1*(tmp77 + p_2*(-nai_300_000_2 + nai_200_000_0*v_2) - nai_100_000_4*u_2) - tmp121*u_1) + (nai_100_000_3*u_1 + nai_200_000_2*p_1 - p_1*tmp162 - tmp162*u_1)/ab_a - u_0*(p_1*tmp121 - u_1*(tmp62 + p_2*tmp123 - u_2*(-tmp124 + nai_200_000_3*v_2)))); // final, weighs 42
  tmp62 = nai_300_000_3*v_2 + (nai_200_001_0 - nai_200_001_1)/ab_a - nai_200_002_1*u_2; // local+recycle, weighs 9
  nai_100_000_4 = nai_200_002_1*v_2 + (nai_200_001_1 - nai_200_001_2)/ab_a - nai_200_000_4*u_2; // local+recycle, weighs 9
  nai_300_000_3 = tmp115 + nai_000_002_0*v_2; // local+recycle, weighs 2
  nai_200_002_1 = nai_200_010_0 + nai_000_002_1*v_2; // local+recycle, weighs 2
  tmp115 = nai_300_000_3*v_0 - nai_200_002_1*u_0; // local+recycle, weighs 4
  nai_200_000_3 = tmp11 + nai_000_002_2*v_2; // local+recycle, weighs 2
  nai_200_010_0 = nai_200_002_1*v_0 - nai_200_000_3*u_0; // local+recycle, weighs 4
  out[25] = 10.8362471994438*tmp107*(p_0*tmp62 + (tmp115 - nai_200_010_0)/ab_a - nai_100_000_4*u_0); // final, weighs 11
  out[26] = 4.84611707178961*tmp107*(nai_200_030_0*v_2 - nai_100_030_1*u_2); // final, weighs 6
  out[27] = 10.8362471994438*tmp107*(p_2*tmp108 + (0.5*nai_200_020_0 - 0.5*nai_200_020_1)/ab_a - tmp69*u_2); // final, weighs 12
  out[28] = 10.8362471994438*tmp107*(p_1*tmp62 - nai_100_000_4*u_1); // final, weighs 6
  out[29] = 4.84611707178961*tmp107*(tmp12 + nai_200_003_0*v_2 - nai_000_001_1*u_2); // final, weighs 7
  tmp124 = nai_300_000_1*v_1 - nai_300_200_1*u_1; // local+recycle, weighs 4
  nai_300_000_2 = tmp87/ab_a; // auto+recycle, weighs 2
  nai_200_001_0 = nai_300_000_2 + v_1*(nai_300_000_0*v_1 - nai_300_000_1*u_1) - tmp124*u_1; // local+recycle, weighs 9
  tmp123 = nai_300_100_0/ab_a; // auto+recycle, weighs 2
  nai_200_001_1 = tmp123 + tmp124*v_1 - u_1*(nai_300_200_1*v_1 - tmp164*u_1); // local+recycle, weighs 9
  tmp121 = tmp3 - tmp169*u_1; // auto+recycle, weighs 3
  nai_200_001_2 = tmp121 + nai_200_010_1*v_1; // local+recycle, weighs 2
  nai_200_000_4 = tmp2 - tmp166*u_1; // auto+recycle, weighs 3
  nai_100_000_1 = nai_200_000_4 + tmp169*v_1; // local+recycle, weighs 2
  nai_200_020_1 = nai_100_000_1*u_0; // auto+recycle, weighs 1
  tmp11 = -nai_200_020_1 + nai_200_001_2*p_0; // local+recycle, weighs 3
  nai_200_000_0 = tmp1 - tmp163*u_1; // auto+recycle, weighs 3
  nai_100_000_2 = nai_200_000_0 + tmp166*v_1; // local+recycle, weighs 2
  tmp77 = nai_100_000_2*u_0; // auto+recycle, weighs 1
  nai_200_020_0 = -tmp77 + nai_100_000_1*p_0; // local+recycle, weighs 3
  nai_000_001_1 = 0.5*nai_200_001_2 - 0.5*nai_100_000_1; // auto+recycle, weighs 3
  nai_200_003_0 = nai_000_001_1/ab_a; // auto+recycle, weighs 2
  nai_100_030_1 = nai_200_003_0 + p_0*tmp11 - nai_200_020_0*u_0; // local+recycle, weighs 5
  nai_200_030_0 = tmp0 - tmp160*u_1; // auto+recycle, weighs 3
  tmp62 = nai_200_030_0 + tmp163*v_1; // local+recycle, weighs 2
  tmp87 = tmp62*u_0; // auto+recycle, weighs 1
  nai_300_100_0 = -tmp87 + nai_100_000_2*p_0; // local+recycle, weighs 3
  nai_100_000_4 = 0.5*nai_100_000_1 - 0.5*nai_100_000_2; // auto+recycle, weighs 3
  tmp12 = nai_100_000_4/ab_a; // auto+recycle, weighs 2
  tmp108 = tmp12 + nai_200_020_0*p_0 - nai_300_100_0*u_0; // local+recycle, weighs 5
  tmp69 = (1.5*nai_100_030_1 - 1.5*tmp108)/ab_a; // auto+recycle, weighs 5
  out[30] = 4.84611707178961*tmp107*(tmp69 + nai_200_001_0*v_0 - nai_200_001_1*u_0); // final, weighs 7
  tmp124 = nai_100_030_1*v_0 + (tmp11 - nai_200_020_0)/ab_a - tmp108*u_0; // local+recycle, weighs 9
  tmp11 = tmp4 + tmp160*v_1 - u_1*(-tmp157 + nai_000_000_5*v_1); // local+recycle, weighs 8
  tmp157 = tmp11*u_0; // auto+recycle, weighs 1
  tmp80 = 0.5*nai_100_000_2 - 0.5*tmp62; // auto, weighs 3
  tmp8 = tmp80/ab_a; // auto, weighs 2
  nai_020_200_2 = tmp8 + nai_300_100_0*p_0 - u_0*(-tmp157 + p_0*tmp62); // local, weighs 8
  nai_300_100_0 = tmp108*v_0 + (nai_200_020_0 - nai_300_100_0)/ab_a - nai_020_200_2*u_0; // local+recycle, weighs 9
  nai_200_020_0 = tmp170*v_1 - tmp70*u_1; // local+recycle, weighs 4
  tmp139 = -tmp139 + p_0*tmp169; // local+recycle, weighs 3
  tmp140 = v_0*(nai_300_100_1*v_1 - tmp170*u_1) + (-tmp139 - tmp140 + nai_200_010_1*p_0)/ab_a - nai_200_020_0*u_0; // local+recycle, weighs 16
  tmp139 = nai_200_020_0*v_0 + (tmp138 + tmp139 - p_0*tmp166)/ab_a - u_0*(tmp70*v_1 - tmp161*u_1); // local+recycle, weighs 15
  out[31] = 10.8362471994438*tmp107*(p_1*tmp124 + (tmp140 - tmp139)/ab_a - nai_300_100_0*u_1); // final, weighs 11
  out[32] = 10.8362471994438*tmp107*(p_2*tmp124 - nai_300_100_0*u_2); // final, weighs 6
  nai_300_100_0 = (nai_200_010_1 - tmp169)/ab_a - nai_100_000_1*u_1; // auto+recycle, weighs 7
  tmp124 = nai_300_100_0 + nai_200_001_2*p_1; // local+recycle, weighs 2
  nai_200_020_0 = (tmp169 - tmp166)/ab_a - nai_100_000_2*u_1; // auto+recycle, weighs 7
  tmp138 = nai_200_020_0 + nai_100_000_1*p_1; // local+recycle, weighs 2
  nai_200_000_4 = nai_200_000_4 + p_1*tmp169; // local+recycle, weighs 2
  nai_000_001_1 = p_1*tmp124 + (nai_000_001_1 + tmp121 - nai_200_000_4 + nai_200_010_1*p_1)/ab_a - tmp138*u_1; // local+recycle, weighs 12
  tmp121 = (tmp166 - tmp163)/ab_a - tmp62*u_1; // auto+recycle, weighs 7
  nai_020_010_2 = tmp121 + nai_100_000_2*p_1; // local, weighs 2
  nai_200_000_0 = nai_200_000_0 + p_1*tmp166; // local+recycle, weighs 2
  nai_100_000_4 = p_1*tmp138 + (nai_100_000_4 + nai_200_000_4 - nai_200_000_0)/ab_a - nai_020_010_2*u_1; // local+recycle, weighs 10
  nai_200_000_4 = nai_000_001_1*v_0 - nai_100_000_4*u_0; // local+recycle, weighs 4
  tmp160 = (tmp163 - tmp160)/ab_a - tmp11*u_1; // auto+recycle, weighs 7
  nai_200_030_0 = nai_020_010_2*p_1 + (nai_200_000_0 + tmp80 - nai_200_030_0 - p_1*tmp163)/ab_a - u_1*(tmp160 + p_1*tmp62); // local+recycle, weighs 15
  tmp80 = nai_100_000_4*v_0 - nai_200_030_0*u_0; // local+recycle, weighs 4
  nai_200_000_0 = (0.5*nai_000_001_1 - 0.5*nai_100_000_4)/ab_a; // auto+recycle, weighs 5
  out[33] = 10.8362471994438*tmp107*(nai_200_000_0 + nai_200_000_4*p_0 - tmp80*u_0); // final, weighs 7
  tmp77 = -tmp77 + nai_100_000_1*v_0; // local+recycle, weighs 3
  tmp87 = -tmp87 + nai_100_000_2*v_0; // local+recycle, weighs 3
  nai_120_001_1 = p_2*tmp77 - tmp87*u_2; // local, weighs 4
  tmp133 = nai_100_000_1*u_2; // auto, weighs 1
  nai_020_001_0 = -tmp133 + nai_200_001_2*p_2; // local, weighs 3
  tmp131 = nai_100_000_2*u_2; // auto, weighs 1
  nai_020_001_1 = -tmp131 + nai_100_000_1*p_2; // local, weighs 3
  nai_010_001_1 = p_2*tmp169 - tmp166*u_2; // local, weighs 4
  tmp129 = tmp62*u_2; // auto, weighs 1
  nai_020_001_2 = -tmp129 + nai_100_000_2*p_2; // local, weighs 3
  out[34] = 18.7689307128126*tmp107*(p_0*(p_1*(p_2*(-nai_200_020_1 + nai_200_001_2*v_0) - tmp77*u_2) + (nai_110_001_0 - nai_110_001_1)/ab_a - nai_120_001_1*u_1) + (0.5*nai_020_001_0*p_1 + 0.5*nai_020_001_2*u_1 + 0.5*(-nai_010_001_1 + nai_200_010_1*p_2 - tmp169*u_2)/ab_a - 0.5*nai_020_001_1*p_1 - 0.5*nai_020_001_1*u_1 - 0.5*(nai_010_001_1 + tmp163*u_2 - p_2*tmp166)/ab_a)/ab_a - u_0*(nai_120_001_1*p_1 + (nai_110_001_1 - nai_110_001_2)/ab_a - u_1*(p_2*tmp87 - u_2*(-tmp157 + tmp62*v_0)))); // final, weighs 71
  nai_200_020_1 = nai_200_003_0 + nai_020_001_0*p_2 - nai_020_001_1*u_2; // local+recycle, weighs 5
  tmp77 = tmp12 + nai_020_001_1*p_2 - nai_020_001_2*u_2; // local+recycle, weighs 5
  nai_110_001_1 = nai_200_020_1*v_0 - tmp77*u_0; // local+recycle, weighs 4
  nai_110_001_2 = tmp11*u_2; // auto+recycle, weighs 1
  nai_110_001_0 = tmp8 + nai_020_001_2*p_2 - u_2*(-nai_110_001_2 + p_2*tmp62); // local+recycle, weighs 8
  tmp157 = tmp77*v_0 - nai_110_001_0*u_0; // local+recycle, weighs 4
  out[35] = 10.8362471994438*tmp107*(nai_110_001_1*p_0 + (0.5*nai_200_020_1 - 0.5*tmp77)/ab_a - tmp157*u_0); // final, weighs 12
  tmp87 = nai_000_010_1*v_1 + (1.5*nai_000_020_1 - 1.5*nai_000_020_2)/ab_a - nai_000_010_2*u_1; // local+recycle, weighs 10
  nai_120_001_1 = v_1*(nai_200_010_2*v_1 + (1.5*nai_000_020_0 - 1.5*nai_000_020_1)/ab_a - nai_000_010_1*u_1) + (tmp74 + 1.5*nai_010_020_0 - 1.5*nai_010_020_1)/ab_a - tmp87*u_1; // local+recycle, weighs 21
  nai_010_001_1 = tmp87*v_1 + (nai_000_010_4 + 1.5*nai_010_020_1 - 1.5*nai_010_020_2)/ab_a - u_1*(nai_000_010_2*v_1 + (1.5*nai_000_020_2 - 1.5*nai_000_020_3)/ab_a - nai_000_010_3*u_1); // local+recycle, weighs 21
  out[36] = 4.84611707178961*tmp107*(nai_120_001_1*v_0 - nai_010_001_1*u_0); // final, weighs 6
  out[37] = 10.8362471994438*tmp107*(nai_200_000_4*p_2 - tmp80*u_2); // final, weighs 6
  out[38] = 10.8362471994438*tmp107*(nai_010_002_1 + nai_110_001_1*p_1 - tmp157*u_1); // final, weighs 7
  nai_200_000_4 = tmp173*v_1 - nai_210_001_1*u_1; // local+recycle, weighs 4
  tmp11 = tmp34 + v_1*(tmp126*v_1 - tmp173*u_1) - nai_200_000_4*u_1; // local+recycle, weighs 9
  nai_010_020_0 = tmp33 + nai_200_000_4*v_1 - u_1*(nai_210_001_1*v_1 - nai_000_001_4*u_1); // local+recycle, weighs 9
  out[39] = 4.84611707178961*tmp107*(tmp11*v_0 - nai_010_020_0*u_0); // final, weighs 6
  nai_110_001_1 = nai_300_000_0*v_2 - nai_300_000_1*u_2; // local+recycle, weighs 4
  tmp34 = nai_300_000_1*v_2 - nai_300_200_1*u_2; // local+recycle, weighs 4
  nai_010_002_1 = nai_300_200_1*v_2 - tmp164*u_2; // local+recycle, weighs 4
  tmp33 = nai_300_100_1*v_2 - tmp170*u_2; // local+recycle, weighs 4
  nai_000_020_2 = tmp170*v_2 - tmp70*u_2; // local+recycle, weighs 4
  nai_010_020_1 = tmp33*v_1 - nai_000_020_2*u_1; // local+recycle, weighs 4
  nai_000_010_4 = tmp70*v_2 - tmp161*u_2; // local+recycle, weighs 4
  nai_000_020_0 = nai_000_020_2*v_1 - nai_000_010_4*u_1; // local+recycle, weighs 4
  out[40] = 8.39372098776651*tmp107*(v_0*(nai_110_001_1*v_1 - tmp34*u_1) + (1.5*nai_010_020_1 - 1.5*nai_000_020_0)/ab_a - u_0*(tmp34*v_1 - nai_010_002_1*u_1)); // final, weighs 20
  tmp74 = d_2*u_1; // auto+recycle, weighs 1
  tmp157 = -tmp74 + d_1*v_1; // local+recycle, weighs 3
  tmp80 = nai_100_001_1*u_1; // auto+recycle, weighs 1
  nai_010_020_2 = -tmp80 + d_2*v_1; // local+recycle, weighs 3
  nai_000_020_1 = nai_010_020_2*u_0; // auto+recycle, weighs 1
  nai_300_000_0 = nai_100_000_0*u_1; // auto+recycle, weighs 1
  tmp164 = -nai_300_000_0 + nai_100_001_1*v_1; // local+recycle, weighs 3
  tmp70 = tmp164*u_0; // auto+recycle, weighs 1
  tmp87 = -tmp70 + nai_010_020_2*p_0; // local+recycle, weighs 3
  nai_300_000_1 = nai_010_020_1*v_0 + (-nai_000_020_1 - tmp87 + p_0*tmp157)/ab_a - nai_000_020_0*u_0; // local+recycle, weighs 12
  nai_300_200_1 = tmp159*u_1; // auto+recycle, weighs 1
  nai_300_100_1 = -nai_300_200_1 + nai_100_000_0*v_1; // local+recycle, weighs 3
  nai_200_000_4 = nai_300_100_1*u_0; // auto+recycle, weighs 1
  tmp170 = nai_000_020_0*v_0 + (nai_200_000_4 + tmp87 - p_0*tmp164)/ab_a - u_0*(nai_000_010_4*v_1 - u_1*(tmp161*v_2 - d_0*u_2)); // local+recycle, weighs 19
  tmp165 = -tmp165 + d_2*p_0; // local+recycle, weighs 3
  tmp33 = tmp33*v_0 + (-tmp110 - tmp165 + d_1*p_0)/ab_a - nai_000_020_2*u_0; // local+recycle, weighs 12
  nai_000_020_2 = nai_000_020_2*v_0 + (nai_000_000_2 + tmp165 - nai_100_001_1*p_0)/ab_a - nai_000_010_4*u_0; // local+recycle, weighs 11
  out[41] = 18.7689307128126*tmp107*(nai_300_000_1*p_1 + (0.5*tmp33 - 0.5*nai_000_020_2)/ab_a - tmp170*u_1); // final, weighs 12
  out[42] = 18.7689307128126*tmp107*(nai_300_000_1*p_2 + (0.5*tmp140 - 0.5*tmp139)/ab_a - tmp170*u_2); // final, weighs 12
  tmp140 = -tmp80 + d_2*p_1; // local+recycle, weighs 3
  nai_000_010_4 = tmp167*v_1 + (-tmp140 - tmp74 + d_1*p_1)/ab_a - tmp63*u_1; // local+recycle, weighs 12
  tmp139 = -nai_300_000_0 + nai_100_001_1*p_1; // local+recycle, weighs 3
  tmp74 = tmp63*v_1 + (tmp140 - tmp139)/ab_a - tmp106*u_1; // local+recycle, weighs 9
  tmp80 = nai_000_010_4*v_0 - tmp74*u_0; // local+recycle, weighs 4
  nai_000_020_4 = tmp74*v_0 - u_0*(tmp106*v_1 + (nai_300_200_1 + tmp139 - nai_100_000_0*p_1)/ab_a - u_1*(nai_000_020_3*v_2 - nai_000_020_4*u_2)); // local+recycle, weighs 19
  out[43] = 18.7689307128126*tmp107*(p_0*tmp80 + (0.5*nai_000_010_4 - 0.5*tmp74)/ab_a - nai_000_020_4*u_0); // final, weighs 12
  nai_300_000_0 = -tmp70 + nai_010_020_2*v_0; // local+recycle, weighs 3
  tmp70 = -nai_200_000_4 + tmp164*v_0; // local+recycle, weighs 3
  tmp87 = nai_300_000_0*p_2 + (0.5*nai_110_000_1 - 0.5*nai_110_000_2)/ab_a - tmp70*u_2; // local+recycle, weighs 10
  nai_300_000_1 = -usq + nai_000_000_5*v_2; // local+recycle, weighs 3
  nai_000_000_5 = p_2*tmp157 + (0.5*nai_200_010_1 - 0.5*tmp169)/ab_a - nai_010_020_2*u_2; // local+recycle, weighs 10
  tmp161 = nai_010_020_2*p_2 + (0.5*tmp169 - 0.5*tmp166)/ab_a - tmp164*u_2; // local+recycle, weighs 10
  nai_200_010_1 = tmp3 - d_2*u_2; // auto+recycle, weighs 3
  nai_300_200_1 = nai_200_010_1 + d_1*p_2; // local+recycle, weighs 2
  usq = tmp2 - nai_100_001_1*u_2; // auto+recycle, weighs 3
  d_0 = usq + d_2*p_2; // local+recycle, weighs 2
  nai_200_000_4 = p_2*tmp164 + (0.5*tmp166 - 0.5*tmp163)/ab_a - nai_300_100_1*u_2; // local+recycle, weighs 10
  nai_000_000_2 = tmp1 - nai_100_000_0*u_2; // auto+recycle, weighs 3
  tmp170 = nai_000_000_2 + nai_100_001_1*p_2; // local+recycle, weighs 2
  out[44] = 32.5087415983314*tmp107*(p_0*(p_1*(p_2*(-nai_000_020_1 + tmp157*v_0) + (0.5*nai_110_000_0 - 0.5*nai_110_000_1)/ab_a - nai_300_000_0*u_2) + (0.5*nai_200_000_2 - 0.5*tmp162)/ab_a - tmp87*u_1) + (0.5*nai_000_000_5*p_1 + 0.5*nai_200_000_4*u_1 + 0.5*(0.5*nai_300_200_1 - 0.5*d_0)/ab_a - 0.5*p_1*tmp161 - 0.5*tmp161*u_1 - 0.5*(0.5*d_0 - 0.5*tmp170)/ab_a)/ab_a - u_0*(p_1*tmp87 + (0.5*tmp162 - 0.5*nai_100_000_3)/ab_a - u_1*(p_2*tmp70 + (0.5*nai_110_000_2 - 0.5*nai_110_000_3)/ab_a - u_2*(nai_300_100_1*v_0 - u_0*(tmp159*v_1 - nai_300_000_1*u_1))))); // final, weighs 85
  tmp3 = nai_300_000_3*v_1 - nai_200_002_1*u_1; // local+recycle, weighs 4
  tmp110 = nai_200_002_1*v_1 - nai_200_000_3*u_1; // local+recycle, weighs 4
  tmp2 = tmp3*v_0 - tmp110*u_0; // local+recycle, weighs 4
  tmp163 = tmp110*v_0 - u_0*(nai_200_000_3*v_1 - u_1*(nai_000_001_3 + nai_000_002_3*v_2)); // local+recycle, weighs 10
  out[45] = 18.7689307128126*tmp107*(p_0*tmp2 + (0.5*tmp3 - 0.5*tmp110)/ab_a - tmp163*u_0); // final, weighs 12
  nai_110_000_2 = nai_200_010_2*v_2 - nai_000_010_1*u_2; // local+recycle, weighs 4
  tmp165 = nai_000_010_1*v_2 - nai_000_010_2*u_2; // local+recycle, weighs 4
  tmp162 = nai_000_010_2*v_2 - nai_000_010_3*u_2; // local+recycle, weighs 4
  out[46] = 8.39372098776651*tmp107*(v_0*(nai_110_000_2*v_1 + (1.5*tmp167 - 1.5*tmp63)/ab_a - tmp165*u_1) - u_0*(tmp165*v_1 + (1.5*tmp63 - 1.5*tmp106)/ab_a - tmp162*u_1)); // final, weighs 26
  out[47] = 18.7689307128126*tmp107*(p_2*tmp80 + (0.5*nai_110_020_0 - 0.5*nai_110_020_1)/ab_a - nai_000_020_4*u_2); // final, weighs 12
  out[48] = 18.7689307128126*tmp107*(p_1*tmp2 + (0.5*tmp115 - 0.5*nai_200_010_0)/ab_a - tmp163*u_1); // final, weighs 12
  nai_000_001_3 = tmp126*v_2 + (1.5*nai_000_002_0 - 1.5*nai_000_002_1)/ab_a - tmp173*u_2; // local+recycle, weighs 10
  nai_000_010_3 = tmp173*v_2 + (1.5*nai_000_002_1 - 1.5*nai_000_002_2)/ab_a - nai_210_001_1*u_2; // local+recycle, weighs 10
  tmp169 = nai_210_001_1*v_2 + (1.5*nai_000_002_2 - 1.5*nai_000_002_3)/ab_a - nai_000_001_4*u_2; // local+recycle, weighs 10
  out[49] = 8.39372098776651*tmp107*(v_0*(nai_000_001_3*v_1 - nai_000_010_3*u_1) - u_0*(nai_000_010_3*v_1 - tmp169*u_1)); // final, weighs 14
  tmp140 = nai_300_000_2 + nai_110_001_1*v_2 - tmp34*u_2; // local+recycle, weighs 5
  nai_110_000_0 = tmp123 + tmp34*v_2 - nai_010_002_1*u_2; // local+recycle, weighs 5
  nai_000_010_1 = nai_200_010_1 + d_1*v_2; // local+recycle, weighs 2
  nai_110_020_0 = usq + d_2*v_2; // local+recycle, weighs 2
  nai_110_000_3 = nai_110_020_0*u_0; // auto+recycle, weighs 1
  nai_000_020_3 = -nai_110_000_3 + nai_000_010_1*p_0; // local+recycle, weighs 3
  tmp166 = nai_000_000_2 + nai_100_001_1*v_2; // local+recycle, weighs 2
  tmp139 = tmp166*u_0; // auto+recycle, weighs 1
  nai_110_000_1 = -tmp139 + nai_110_020_0*p_0; // local+recycle, weighs 3
  nai_000_010_2 = 0.5*nai_000_010_1 - 0.5*nai_110_020_0; // auto+recycle, weighs 3
  tmp80 = nai_000_010_2/ab_a; // auto+recycle, weighs 2
  nai_000_020_4 = tmp80 + nai_000_020_3*p_0 - nai_110_000_1*u_0; // local+recycle, weighs 5
  nai_110_020_1 = tmp0 - tmp159*u_2; // auto+recycle, weighs 3
  nai_000_020_1 = nai_110_020_1 + nai_100_000_0*v_2; // local+recycle, weighs 2
  tmp106 = nai_000_020_1*u_0; // auto+recycle, weighs 1
  nai_300_000_0 = -tmp106 + p_0*tmp166; // local+recycle, weighs 3
  nai_100_000_3 = 0.5*nai_110_020_0 - 0.5*tmp166; // auto+recycle, weighs 3
  tmp1 = nai_100_000_3/ab_a; // auto+recycle, weighs 2
  nai_200_000_2 = tmp1 + nai_110_000_1*p_0 - nai_300_000_0*u_0; // local+recycle, weighs 5
  tmp70 = (1.5*nai_000_020_4 - 1.5*nai_200_000_2)/ab_a; // auto+recycle, weighs 5
  out[50] = 4.84611707178961*tmp107*(tmp70 + tmp140*v_0 - nai_110_000_0*u_0); // final, weighs 7
  tmp87 = nai_000_020_4*v_0 + (nai_000_020_3 - nai_110_000_1)/ab_a - nai_200_000_2*u_0; // local+recycle, weighs 9
  nai_300_000_1 = tmp4 + tmp159*v_2 - nai_300_000_1*u_2; // local+recycle, weighs 5
  tmp0 = nai_300_000_1*u_0; // auto+recycle, weighs 1
  tmp115 = 0.5*tmp166 - 0.5*nai_000_020_1; // auto+recycle, weighs 3
  tmp126 = tmp115/ab_a; // auto+recycle, weighs 2
  nai_200_010_0 = tmp126 + nai_300_000_0*p_0 - u_0*(-tmp0 + nai_000_020_1*p_0); // local+recycle, weighs 8
  nai_200_010_1 = nai_200_000_2*v_0 + (nai_110_000_1 - nai_300_000_0)/ab_a - nai_200_010_0*u_0; // local+recycle, weighs 9
  out[51] = 10.8362471994438*tmp107*(p_1*tmp87 - nai_200_010_1*u_1); // final, weighs 6
  out[52] = 10.8362471994438*tmp107*(p_2*tmp87 + (tmp33 - nai_000_020_2)/ab_a - nai_200_010_1*u_2); // final, weighs 11
  nai_300_000_2 = nai_110_020_0*u_1; // auto+recycle, weighs 1
  tmp123 = -nai_300_000_2 + nai_000_010_1*p_1; // local+recycle, weighs 3
  usq = tmp166*u_1; // auto+recycle, weighs 1
  tmp4 = -usq + nai_110_020_0*p_1; // local+recycle, weighs 3
  tmp173 = tmp80 + p_1*tmp123 - tmp4*u_1; // local+recycle, weighs 5
  nai_000_000_2 = nai_000_020_1*u_1; // auto+recycle, weighs 1
  tmp167 = -nai_000_000_2 + p_1*tmp166; // local+recycle, weighs 3
  tmp2 = tmp1 + p_1*tmp4 - tmp167*u_1; // local+recycle, weighs 5
  tmp63 = tmp173*v_0 - tmp2*u_0; // local+recycle, weighs 4
  nai_210_001_1 = nai_300_000_1*u_1; // auto+recycle, weighs 1
  nai_200_010_2 = tmp126 + p_1*tmp167 - u_1*(-nai_210_001_1 + nai_000_020_1*p_1); // local+recycle, weighs 8
  tmp163 = tmp2*v_0 - nai_200_010_2*u_0; // local+recycle, weighs 4
  out[53] = 10.8362471994438*tmp107*(p_0*tmp63 + (0.5*tmp173 - 0.5*tmp2)/ab_a - tmp163*u_0); // final, weighs 12
  nai_110_001_1 = -tmp139 + nai_110_020_0*v_0; // local+recycle, weighs 3
  nai_000_002_0 = -tmp106 + tmp166*v_0; // local+recycle, weighs 3
  tmp34 = nai_110_001_1*p_2 + (nai_200_000_1 - nai_000_000_3)/ab_a - nai_000_002_0*u_2; // local+recycle, weighs 9
  nai_000_002_1 = (d_1 - d_2)/ab_a - nai_110_020_0*u_2; // auto+recycle, weighs 7
  nai_000_002_2 = nai_000_002_1 + nai_000_010_1*p_2; // local+recycle, weighs 2
  nai_010_002_1 = (d_2 - nai_100_001_1)/ab_a - tmp166*u_2; // auto+recycle, weighs 7
  nai_000_001_4 = nai_010_002_1 + nai_110_020_0*p_2; // local+recycle, weighs 2
  tmp33 = (nai_100_001_1 - nai_100_000_0)/ab_a - nai_000_020_1*u_2; // auto+recycle, weighs 7
  nai_000_002_3 = tmp33 + p_2*tmp166; // local+recycle, weighs 2
  out[54] = 18.7689307128126*tmp107*(p_0*(p_1*(p_2*(-nai_110_000_3 + nai_000_010_1*v_0) + (tmp122 - nai_200_000_1)/ab_a - nai_110_001_1*u_2) - tmp34*u_1) + (0.5*nai_000_002_2*p_1 + 0.5*nai_000_002_3*u_1 - 0.5*nai_000_001_4*p_1 - 0.5*nai_000_001_4*u_1)/ab_a - u_0*(p_1*tmp34 - u_1*(nai_000_002_0*p_2 + (nai_000_000_3 - nai_000_000_4)/ab_a - u_2*(-tmp0 + nai_000_020_1*v_0)))); // final, weighs 52
  nai_000_020_2 = nai_000_002_2*p_2 + (nai_000_010_2 + nai_300_200_1 - d_0)/ab_a - nai_000_001_4*u_2; // local+recycle, weighs 10
  nai_110_000_3 = nai_000_001_4*p_2 + (d_0 + nai_100_000_3 - tmp170)/ab_a - nai_000_002_3*u_2; // local+recycle, weighs 10
  nai_000_020_3 = nai_000_020_2*v_0 - nai_110_000_3*u_0; // local+recycle, weighs 4
  tmp139 = (nai_100_000_0 - tmp159)/ab_a - nai_300_000_1*u_2; // auto+recycle, weighs 7
  nai_110_000_1 = nai_000_002_3*p_2 + (tmp115 + tmp170 - nai_110_020_1 - nai_100_000_0*p_2)/ab_a - u_2*(tmp139 + nai_000_020_1*p_2); // local+recycle, weighs 15
  nai_000_010_2 = nai_110_000_3*v_0 - nai_110_000_1*u_0; // local+recycle, weighs 4
  tmp80 = (0.5*nai_000_020_2 - 0.5*nai_110_000_3)/ab_a; // auto+recycle, weighs 5
  out[55] = 10.8362471994438*tmp107*(tmp80 + nai_000_020_3*p_0 - nai_000_010_2*u_0); // final, weighs 7
  nai_110_020_1 = tmp21 + nai_110_000_2*v_2 - tmp165*u_2; // local+recycle, weighs 5
  nai_200_000_1 = tmp20 + tmp165*v_2 - tmp162*u_2; // local+recycle, weighs 5
  out[56] = 4.84611707178961*tmp107*(nai_110_020_1*v_0 - nai_200_000_1*u_0); // final, weighs 6
  out[57] = 10.8362471994438*tmp107*(nai_300_200_0 + p_2*tmp63 - tmp163*u_2); // final, weighs 7
  out[58] = 10.8362471994438*tmp107*(nai_000_020_3*p_1 - nai_000_010_2*u_1); // final, weighs 6
  tmp106 = nai_000_001_3*v_2 + (nai_000_001_2 + 1.5*nai_300_000_3 - 1.5*nai_200_002_1)/ab_a - nai_000_010_3*u_2; // local+recycle, weighs 11
  nai_300_000_0 = nai_000_010_3*v_2 + (tmp88 + 1.5*nai_200_002_1 - 1.5*nai_200_000_3)/ab_a - tmp169*u_2; // local+recycle, weighs 11
  out[59] = 4.84611707178961*tmp107*(tmp106*v_0 - nai_300_000_0*u_0); // final, weighs 6
  nai_000_000_4 = nai_300_100_0 + nai_200_001_2*v_1; // local+recycle, weighs 2
  nai_100_000_3 = nai_200_020_0 + nai_100_000_1*v_1; // local+recycle, weighs 2
  tmp1 = nai_000_000_4*p_0 - nai_100_000_3*u_0; // local+recycle, weighs 4
  nai_300_200_0 = tmp121 + nai_100_000_2*v_1; // local+recycle, weighs 2
  tmp87 = nai_100_000_3*p_0 - nai_300_200_0*u_0; // local+recycle, weighs 4
  nai_300_000_1 = 0.5*nai_000_000_4 - 0.5*nai_100_000_3; // auto+recycle, weighs 3
  nai_300_100_0 = nai_300_000_1/ab_a; // auto+recycle, weighs 2
  nai_300_000_3 = nai_300_100_0 + p_0*tmp1 - tmp87*u_0; // local+recycle, weighs 5
  nai_200_002_1 = tmp160 + tmp62*v_1; // local+recycle, weighs 2
  tmp0 = 0.5*nai_100_000_3 - 0.5*nai_300_200_0; // auto+recycle, weighs 3
  tmp115 = tmp0/ab_a; // auto+recycle, weighs 2
  nai_200_000_3 = tmp115 + p_0*tmp87 - u_0*(nai_300_200_0*p_0 - nai_200_002_1*u_0); // local+recycle, weighs 9
  out[60] = 2.16724943988876*tmp107*(nai_300_000_3*p_0 + (tmp1 - tmp87)/ab_a - nai_200_000_3*u_0); // final, weighs 11
  out[61] = 4.84611707178961*tmp107*(tmp69 + nai_300_000_3*p_1 - nai_200_000_3*u_1); // final, weighs 7
  out[62] = 4.84611707178961*tmp107*(nai_300_000_3*p_2 - nai_200_000_3*u_2); // final, weighs 6
  tmp69 = nai_000_000_4*p_1 + (1.5*nai_200_001_2 - 1.5*nai_100_000_1)/ab_a - nai_100_000_3*u_1; // local+recycle, weighs 10
  tmp126 = nai_100_000_3*p_1 + (1.5*nai_100_000_1 - 1.5*nai_100_000_2)/ab_a - nai_300_200_0*u_1; // local+recycle, weighs 10
  tmp124 = p_1*tmp69 + (nai_300_000_1 + 1.5*tmp124 - 1.5*tmp138)/ab_a - tmp126*u_1; // local+recycle, weighs 11
  nai_200_010_1 = p_1*tmp126 + (tmp0 + 1.5*tmp138 - 1.5*nai_020_010_2)/ab_a - u_1*(nai_300_200_0*p_1 + (1.5*nai_100_000_2 - 1.5*tmp62)/ab_a - nai_200_002_1*u_1); // local+recycle, weighs 21
  out[63] = 4.84611707178961*tmp107*(p_0*tmp124 - nai_200_010_1*u_0); // final, weighs 6
  nai_300_200_1 = nai_000_000_4*p_2 - nai_100_000_3*u_2; // local+recycle, weighs 4
  tmp121 = nai_100_000_3*p_2 - nai_300_200_0*u_2; // local+recycle, weighs 4
  nai_020_010_2 = nai_300_200_0*p_2 - nai_200_002_1*u_2; // local+recycle, weighs 4
  out[64] = 8.39372098776651*tmp107*(p_0*(nai_300_200_1*p_1 + (1.5*nai_020_001_0 - 1.5*nai_020_001_1)/ab_a - tmp121*u_1) - u_0*(p_1*tmp121 + (1.5*nai_020_001_1 - 1.5*nai_020_001_2)/ab_a - nai_020_010_2*u_1)); // final, weighs 26
  v_0 = nai_300_100_0 + nai_300_200_1*p_2 - tmp121*u_2; // local+recycle, weighs 5
  d_0 = tmp115 + p_2*tmp121 - nai_020_010_2*u_2; // local+recycle, weighs 5
  out[65] = 4.84611707178961*tmp107*(p_0*v_0 - d_0*u_0); // final, weighs 6
  out[66] = 2.16724943988876*tmp107*(p_1*tmp124 + (tmp69 - tmp126 + 1.5*nai_000_001_1 - 1.5*nai_100_000_4)/ab_a - nai_200_010_1*u_1); // final, weighs 15
  out[67] = 4.84611707178961*tmp107*(p_2*tmp124 - nai_200_010_1*u_2); // final, weighs 6
  d_1 = (1.5*nai_200_020_1 - 1.5*tmp77)/ab_a; // auto+recycle, weighs 5
  out[68] = 4.84611707178961*tmp107*(d_1 + p_1*v_0 - d_0*u_1); // final, weighs 7
  out[69] = 2.16724943988876*tmp107*(p_2*v_0 + (nai_300_200_1 - tmp121)/ab_a - d_0*u_2); // final, weighs 11
  out[70] = 4.84611707178961*tmp107*(nai_200_001_0*v_2 - nai_200_001_1*u_2); // final, weighs 6
  d_2 = nai_100_030_1*v_2 - tmp108*u_2; // local+recycle, weighs 4
  nai_100_000_0 = tmp108*v_2 - nai_020_200_2*u_2; // local+recycle, weighs 4
  tmp170 = (nai_010_020_1 - nai_000_020_0)/ab_a; // auto+recycle, weighs 4
  out[71] = 10.8362471994438*tmp107*(tmp170 + d_2*p_1 - nai_100_000_0*u_1); // final, weighs 7
  out[72] = 10.8362471994438*tmp107*(d_2*p_2 + (0.5*nai_100_030_1 - 0.5*tmp108)/ab_a - nai_100_000_0*u_2); // final, weighs 12
  nai_000_000_3 = nai_000_001_1*v_2 - nai_100_000_4*u_2; // local+recycle, weighs 4
  tmp63 = nai_100_000_4*v_2 - nai_200_030_0*u_2; // local+recycle, weighs 4
  out[73] = 10.8362471994438*tmp107*(nai_000_000_3*p_0 - tmp63*u_0); // final, weighs 6
  nai_200_020_0 = -tmp131 + nai_100_000_1*v_2; // local+recycle, weighs 3
  tmp122 = -tmp129 + nai_100_000_2*v_2; // local+recycle, weighs 3
  tmp163 = tmp12 + nai_200_020_0*p_2 - tmp122*u_2; // local+recycle, weighs 5
  out[74] = 18.7689307128126*tmp107*(p_0*(p_1*(nai_200_003_0 + p_2*(-tmp133 + nai_200_001_2*v_2) - nai_200_020_0*u_2) + (nai_000_000_5 - tmp161)/ab_a - tmp163*u_1) - u_0*(p_1*tmp163 + (tmp161 - nai_200_000_4)/ab_a - u_1*(tmp8 + p_2*tmp122 - u_2*(-nai_110_001_2 + tmp62*v_2)))); // final, weighs 40
  tmp138 = nai_200_020_1*v_2 + (nai_020_001_0 - nai_020_001_1)/ab_a - tmp77*u_2; // local+recycle, weighs 9
  nai_110_000_2 = tmp77*v_2 + (nai_020_001_1 - nai_020_001_2)/ab_a - nai_110_001_0*u_2; // local+recycle, weighs 9
  out[75] = 10.8362471994438*tmp107*(p_0*tmp138 - nai_110_000_2*u_0); // final, weighs 6
  out[76] = 4.84611707178961*tmp107*(nai_120_001_1*v_2 - nai_010_001_1*u_2); // final, weighs 6
  out[77] = 10.8362471994438*tmp107*(nai_200_000_0 + nai_000_000_3*p_2 - tmp63*u_2); // final, weighs 7
  out[78] = 10.8362471994438*tmp107*(p_1*tmp138 + (tmp3 - tmp110)/ab_a - nai_110_000_2*u_1); // final, weighs 11
  out[79] = 4.84611707178961*tmp107*(d_1 + tmp11*v_2 - nai_010_020_0*u_2); // final, weighs 7
  out[80] = 4.84611707178961*tmp107*(tmp140*v_1 - nai_110_000_0*u_1); // final, weighs 6
  nai_110_001_1 = nai_000_020_4*v_1 - nai_200_000_2*u_1; // local+recycle, weighs 4
  nai_100_001_1 = nai_200_000_2*v_1 - nai_200_010_0*u_1; // local+recycle, weighs 4
  out[81] = 10.8362471994438*tmp107*(nai_110_001_1*p_1 + (0.5*nai_000_020_4 - 0.5*nai_200_000_2)/ab_a - nai_100_001_1*u_1); // final, weighs 12
  out[82] = 10.8362471994438*tmp107*(tmp170 + nai_110_001_1*p_2 - nai_100_001_1*u_2); // final, weighs 7
  nai_000_001_1 = tmp173*v_1 + (tmp123 - tmp4)/ab_a - tmp2*u_1; // local+recycle, weighs 9
  nai_000_002_0 = tmp2*v_1 + (tmp4 - tmp167)/ab_a - nai_200_010_2*u_1; // local+recycle, weighs 9
  out[83] = 10.8362471994438*tmp107*(nai_000_001_1*p_0 - nai_000_002_0*u_0); // final, weighs 6
  tmp165 = -usq + nai_110_020_0*v_1; // local+recycle, weighs 3
  nai_000_001_2 = -nai_000_000_2 + tmp166*v_1; // local+recycle, weighs 3
  tmp34 = p_2*tmp165 + (nai_010_020_2 - tmp164)/ab_a - nai_000_001_2*u_2; // local+recycle, weighs 9
  out[84] = 18.7689307128126*tmp107*(p_0*(p_1*(p_2*(-nai_300_000_2 + nai_000_010_1*v_1) + (tmp157 - nai_010_020_2)/ab_a - tmp165*u_2) + (0.5*nai_000_002_2 - 0.5*nai_000_001_4)/ab_a - tmp34*u_1) - u_0*(p_1*tmp34 + (0.5*nai_000_001_4 - 0.5*nai_000_002_3)/ab_a - u_1*(nai_000_001_2*p_2 + (tmp164 - nai_300_100_1)/ab_a - u_2*(-nai_210_001_1 + nai_000_020_1*v_1)))); // final, weighs 50
  nai_200_003_0 = nai_000_020_2*v_1 - nai_110_000_3*u_1; // local+recycle, weighs 4
  tmp162 = nai_110_000_3*v_1 - nai_110_000_1*u_1; // local+recycle, weighs 4
  out[85] = 10.8362471994438*tmp107*(nai_200_003_0*p_0 - tmp162*u_0); // final, weighs 6
  nai_000_001_3 = (1.5*tmp173 - 1.5*tmp2)/ab_a; // auto+recycle, weighs 5
  out[86] = 4.84611707178961*tmp107*(nai_000_001_3 + nai_110_020_1*v_1 - nai_200_000_1*u_1); // final, weighs 7
  out[87] = 10.8362471994438*tmp107*(nai_000_001_1*p_2 + (nai_000_010_4 - tmp74)/ab_a - nai_000_002_0*u_2); // final, weighs 11
  out[88] = 10.8362471994438*tmp107*(tmp80 + nai_200_003_0*p_1 - tmp162*u_1); // final, weighs 7
  out[89] = 4.84611707178961*tmp107*(tmp106*v_1 - nai_300_000_0*u_1); // final, weighs 6
  tmp159 = nai_000_002_1 + nai_000_010_1*v_2; // local+recycle, weighs 2
  tmp88 = nai_010_002_1 + nai_110_020_0*v_2; // local+recycle, weighs 2
  nai_000_010_3 = p_0*tmp159 - tmp88*u_0; // local+recycle, weighs 4
  tmp169 = tmp33 + tmp166*v_2; // local+recycle, weighs 2
  tmp140 = p_0*tmp88 - tmp169*u_0; // local+recycle, weighs 4
  nai_110_000_0 = 0.5*tmp159 - 0.5*tmp88; // auto+recycle, weighs 3
  nai_010_020_1 = nai_110_000_0/ab_a; // auto+recycle, weighs 2
  tmp160 = nai_010_020_1 + nai_000_010_3*p_0 - tmp140*u_0; // local+recycle, weighs 5
  nai_110_001_2 = tmp139 + nai_000_020_1*v_2; // local+recycle, weighs 2
  nai_000_010_4 = 0.5*tmp88 - 0.5*tmp169; // auto+recycle, weighs 3
  tmp20 = nai_000_010_4/ab_a; // auto+recycle, weighs 2
  nai_000_020_3 = tmp20 + p_0*tmp140 - u_0*(p_0*tmp169 - nai_110_001_2*u_0); // local+recycle, weighs 9
  out[90] = 2.16724943988876*tmp107*(p_0*tmp160 + (nai_000_010_3 - tmp140)/ab_a - nai_000_020_3*u_0); // final, weighs 11
  out[91] = 4.84611707178961*tmp107*(p_1*tmp160 - nai_000_020_3*u_1); // final, weighs 6
  out[92] = 4.84611707178961*tmp107*(tmp70 + p_2*tmp160 - nai_000_020_3*u_2); // final, weighs 7
  nai_000_020_0 = p_1*tmp159 - tmp88*u_1; // local+recycle, weighs 4
  tmp139 = p_1*tmp88 - tmp169*u_1; // local+recycle, weighs 4
  nai_110_000_1 = nai_010_020_1 + nai_000_020_0*p_1 - tmp139*u_1; // local+recycle, weighs 5
  nai_110_001_0 = tmp20 + p_1*tmp139 - u_1*(p_1*tmp169 - nai_110_001_2*u_1); // local+recycle, weighs 9
  out[93] = 4.84611707178961*tmp107*(nai_110_000_1*p_0 - nai_110_001_0*u_0); // final, weighs 6
  nai_000_010_2 = p_2*tmp159 + (1.5*nai_000_010_1 - 1.5*nai_110_020_0)/ab_a - tmp88*u_2; // local+recycle, weighs 10
  nai_100_030_1 = p_2*tmp88 + (1.5*nai_110_020_0 - 1.5*tmp166)/ab_a - tmp169*u_2; // local+recycle, weighs 10
  tmp74 = p_2*tmp169 + (1.5*tmp166 - 1.5*nai_000_020_1)/ab_a - nai_110_001_2*u_2; // local+recycle, weighs 10
  out[94] = 8.39372098776651*tmp107*(p_0*(nai_000_010_2*p_1 - nai_100_030_1*u_1) - u_0*(nai_100_030_1*p_1 - tmp74*u_1)); // final, weighs 14
  tmp21 = nai_000_010_2*p_2 + (nai_110_000_0 + 1.5*nai_000_002_2 - 1.5*nai_000_001_4)/ab_a - nai_100_030_1*u_2; // local+recycle, weighs 11
  nai_200_030_0 = nai_100_030_1*p_2 + (nai_000_010_4 + 1.5*nai_000_001_4 - 1.5*nai_000_002_3)/ab_a - tmp74*u_2; // local+recycle, weighs 11
  out[95] = 4.84611707178961*tmp107*(p_0*tmp21 - nai_200_030_0*u_0); // final, weighs 6
  out[96] = 2.16724943988876*tmp107*(nai_110_000_1*p_1 + (nai_000_020_0 - tmp139)/ab_a - nai_110_001_0*u_1); // final, weighs 11
  out[97] = 4.84611707178961*tmp107*(nai_000_001_3 + nai_110_000_1*p_2 - nai_110_001_0*u_2); // final, weighs 7
  out[98] = 4.84611707178961*tmp107*(p_1*tmp21 - nai_200_030_0*u_1); // final, weighs 6
  out[99] = 2.16724943988876*tmp107*(p_2*tmp21 + (nai_000_010_2 - nai_100_030_1 + 1.5*nai_000_020_2 - 1.5*nai_110_000_3)/ab_a - nai_200_030_0*u_2); // final, weighs 15
  // total weight = 3600
}

typedef void (*fntype)(double*, double, double*, double, double*, double*);
static const fntype fns[49] = {gint2_nai_pF_pF, gint2_nai_pD_pF, gint2_nai_SP_pF, gint2_nai_S_pF, gint2_nai_P_pF, gint2_nai_cD_pF, gint2_nai_cF_pF, gint2_nai_pF_pD, gint2_nai_pD_pD, gint2_nai_SP_pD, gint2_nai_S_pD, gint2_nai_P_pD, gint2_nai_cD_pD, gint2_nai_cF_pD, gint2_nai_pF_SP, gint2_nai_pD_SP, gint2_nai_SP_SP, gint2_nai_S_SP, gint2_nai_P_SP, gint2_nai_cD_SP, gint2_nai_cF_SP, gint2_nai_pF_S, gint2_nai_pD_S, gint2_nai_SP_S, gint2_nai_S_S, gint2_nai_P_S, gint2_nai_cD_S, gint2_nai_cF_S, gint2_nai_pF_P, gint2_nai_pD_P, gint2_nai_SP_P, gint2_nai_S_P, gint2_nai_P_P, gint2_nai_cD_P, gint2_nai_cF_P, gint2_nai_pF_cD, gint2_nai_pD_cD, gint2_nai_SP_cD, gint2_nai_S_cD, gint2_nai_P_cD, gint2_nai_cD_cD, gint2_nai_cF_cD, gint2_nai_pF_cF, gint2_nai_pD_cF, gint2_nai_SP_cF, gint2_nai_S_cF, gint2_nai_P_cF, gint2_nai_cD_cF, gint2_nai_cF_cF};

void gint2_nai_dispatch(int a_s, double* a, double a_a, int b_s, double* b, double b_a, double* c, double* out)
{
  fns[3+a_s+7*(3+b_s)](a, a_a, b, b_a, c, out);
}
