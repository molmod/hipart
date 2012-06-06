// HiPart is a program to analyze the electronic structure of molecules with
// fuzzy-atom partitioning methods.
// Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
//--


#include <stdlib.h>
#include <string.h>


int reorder_dmat_c(double* dmat, int* permutation, int num_dof) {
  int result, size, i, j, ip, jp, k, kp;
  double *tmp;

  // The density matrix is a symetric matrix and we suppose that it is
  // represented by its lower triangular part in compact row-major storage.

  // The permutation vector contains the old positions of the rows and columns.

  result = 0;

  size = ((num_dof+1)*num_dof)/2;
  tmp = malloc(size*sizeof(double));
  if (tmp==NULL) { result = -1; goto EXIT; }

  // Copy and reorder
  k = 0;
  for (i=0; i<num_dof; i++) {
    ip = permutation[i];
    for (j=0; j<=i; j++) {
      jp = permutation[j];
      if (ip >= jp) {
        kp = (ip*(ip+1))/2+jp;
      } else {
        kp = (jp*(jp+1))/2+ip;
      }
      // kp is the position in the old array.
      // k is the position in the new array.
      // printf("i=%d  j=%d  k=%d    ip=%d  jp=%d  kp=%d\n", i, j, k, ip, jp, kp);
      tmp[k] = dmat[kp];
      k++;
    }
  }

  // fast bitwise copy;
  memcpy(dmat, tmp, size*sizeof(double));

EXIT:
  free(tmp);
  return result;
}


void add_orbitals_to_dmat(double* orbitals, double* dmat, int num_orbitals, int num_dof) {
  int i, j, k, o;
  double *orbital;

  orbital = orbitals;
  for (o=0; o<num_orbitals; o++) {
    k = 0;
    for (i=0; i<num_dof; i++) {
      for (j=0; j<=i; j++) {
        dmat[k] += orbital[i]*orbital[j];
        k++;
      }
    }
    orbital += num_dof;
  }
}


void dmat_to_full(double* dmat, double* full, int num_dof) {
  int i, j;

  for (i=0; i<num_dof; i++) {
    for (j=0; j<=i; j++) {
      full[i*num_dof+j] = *dmat;
      full[j*num_dof+i] = *dmat;
      dmat++;
    }
  }
}
