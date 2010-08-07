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


#ifndef GINT1_FN_H
#define GINT1_FN_H

int gint1_fns_basis(double* weights, double* orbs, double* points,
  double* centers, int* shell_types, int* shell_map,  int* num_primitives,
  double* ccoeffs, double* exponents, int num_orbs, int num_dof, int num_points,
  int num_centers, int num_shells, int num_ccoeffs, int num_exponents);

int gint1_fn_dmat(double* dmat, double* density, double* points,
  double* centers, int* shell_types, int* shell_map,  int* num_primitives,
  double* ccoeffs, double* exponents, int num_dmat, int num_points,
  int num_centers, int num_shells, int num_ccoeffs, int num_exponents);

int gint1_fn_overlap(double* overlap, double* points, double* weights,
  double* centers, int* shell_types, int* shell_map, int* num_primitives,
  double* ccoeffs, double* exponents, int num_overlap, int num_points,
  int num_centers, int num_shells, int num_ccoeffs, int num_exponents);

#endif
