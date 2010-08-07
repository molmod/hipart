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


#ifndef GINT2_NAI_AUTO_H
#define GINT2_NAI_AUTO_H

#define MAX_SHELL 3
#define NUM_SHELL_TYPES 7
#define MAX_SHELL_DOF 10
#define CHECK_ALLOC(pointer) if (pointer==NULL) {result = -1; goto EXIT; }
#define CHECK_SHELL(shell_type) if (abs(shell_type) > MAX_SHELL) { result = -2; goto EXIT; }
#define GET_SHELL_DOF(shell_type) ((shell_type<-1)?(-2*shell_type+1):((shell_type==-1)?(4):(((shell_type+1)*(shell_type+2))/2)))

void gint2_nai_pF_pF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pD_pF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_SP_pF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_S_pF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_P_pF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cD_pF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cF_pF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pF_pD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pD_pD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_SP_pD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_S_pD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_P_pD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cD_pD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cF_pD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pF_SP(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pD_SP(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_SP_SP(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_S_SP(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_P_SP(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cD_SP(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cF_SP(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pF_S(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pD_S(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_SP_S(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_S_S(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_P_S(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cD_S(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cF_S(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pF_P(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pD_P(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_SP_P(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_S_P(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_P_P(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cD_P(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cF_P(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pF_cD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pD_cD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_SP_cD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_S_cD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_P_cD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cD_cD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cF_cD(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pF_cF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_pD_cF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_SP_cF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_S_cF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_P_cF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cD_cF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_cF_cF(double* a, double a_a, double* b, double b_a, double* c, double* out);
void gint2_nai_dispatch(int a_s, double* a, double a_a, int b_s, double* b, double b_a, double* c, double* out);

#endif
