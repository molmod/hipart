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

#include <math.h>
#include "gaux.h"

#ifdef BOUNDS_CHECKING
#include <stdio.h>
#endif

#include "gaux.inc"
#define SQRT_PI_d2      8.86226925452758013649e-1 // sqrt(pi)/2

/**********************************************************
 *
 *   auxiliary function for the evaluation of coulomb integrals
 *
 **********************************************************/

int gaux_get_maxn(void) {
  return gaux_maxn-6;
}

double gaux_fn_limit(double T, int n) {
  if (T<0 || n<0) {
    return -1;
  } else {
    int i;
    double result;
    result = 1;
    for (i=1; i<=n; i++) {
      result *= 0.5*(2*i-1)/T;
    }
    return result*SQRT_PI_d2/sqrt(T);
  }
}

double gaux(double T, int n) {
  if (T<0 || n<0 || n>gaux_get_maxn()) {
#ifdef BOUNDS_CHECKING
    printf("BOUNDS CHECKING ERROR in gaux_fn: T=%e n=%i n_max=%i\n", T, n, gint_gaux__get_maxn());
#endif
    return -1;
  } else {
    if (isnan(T)) {
      return 0;
    } else if (round(T*gaux_resolution)>=(gaux_Fsizes[n]-1)) {
      return gaux_fn_limit(T, n);
    } else {
      double Tdelta, tmp, result;
      int i;
      i = (int)(round(T*gaux_resolution));
      //printf("i=%i  n=%i  T=%f\n", i, n, T);
      // It used to be like this, but it gives silly roundoff errors:
      //Tdelta = ((float)i)/gaux_resolution - T;
      // This is slightly better:
      Tdelta = (i - T*gaux_resolution)/gaux_resolution;
      result = 0.0;
      tmp = 1;    result += gaux_Fdata[n][i];
      tmp *= Tdelta;   result += gaux_Fdata[n+1][i]*tmp;
      tmp *= Tdelta/2; result += gaux_Fdata[n+2][i]*tmp;
      tmp *= Tdelta/3; result += gaux_Fdata[n+3][i]*tmp;
      tmp *= Tdelta/4; result += gaux_Fdata[n+4][i]*tmp;
      tmp *= Tdelta/5; result += gaux_Fdata[n+5][i]*tmp;
      tmp *= Tdelta/6; result += gaux_Fdata[n+6][i]*tmp;
      return result;
    }
  }
}
