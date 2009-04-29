// HiPart is a tool to analyse molecular densities with the hirshfeld partitioning scheme
// Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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



#include <stdlib.h>
#include <stdio.h>

void spline_construct(double *x, double *y, double *d, int n) {
  // Constructs the derivatives (*d) based on the *x and *y values
  int i;
  double tmp;
  double *diag_low;
  double *diag_mid;
  double *diag_up;
  double *right;

  diag_low = malloc((n-1)*sizeof(double));
  diag_mid = malloc(n*sizeof(double));
  diag_up = malloc((n-1)*sizeof(double));
  right =  malloc(n*sizeof(double));

  // A) Setup the tri-diagonal system
  diag_mid[0] = 2*(x[1]-x[0]);
  diag_up[0] = (x[1]-x[0]);
  right[0] = 3*(y[1]-y[0]);
  for (i=1; i<n-1; i++) {
    diag_low[i-1] = (x[i]-x[i-1]);
    diag_mid[i] = 2*(x[i]-x[i-1]) + 2*(x[i+1]-x[i]);
    diag_up[i] = (x[i+1]-x[i]);
    right[i] = 3*(y[i+1]-y[i-1]);
  }
  diag_low[n-2] = (x[n-1]-x[n-2]);
  diag_mid[n-1] = 2*(x[n-1]-x[n-2]);
  right[n-1] = 3*(y[n-1]-y[n-2]);

  // B) Solve the tri-diagonal system
  diag_up[0] /= diag_mid[0];
  right[0] /= diag_mid[0];
  for(i = 1; i < n; i++){
    tmp = 1.0/(diag_mid[i] - diag_up[i - 1]*diag_low[i-1]);
    if (i < n-1) diag_up[i] *= tmp;
    right[i] = (right[i] - right[i - 1]*diag_low[i-1])*tmp;
  }

  d[n - 1] = right[n - 1];
  for(i = n - 2; i >= 0; i--) {
    d[i] = right[i] - diag_up[i]*d[i + 1];
  }
  free(diag_low);
  free(diag_mid);
  free(diag_up);
  free(right);
}

void spline_eval(double *x, double *y, double *d, int n, double* new_x, double* new_y, int new_n) {
    int i;
    for (i=0; i<new_n; i++) {
        double x_cur;
        int j;
        x_cur = new_x[i];
        if (x_cur <= x[0]) {
          new_y[i] = y[0];
          continue;
        } else if (x_cur >= x[n-1]) {
          new_y[i] = y[n-1];
          continue;
        } else {
          //printf("x_cur = %f\n", x_cur);
          // find the index of the interval in which new_x[i] lies.
          int j_low, j_high;
          j = n/2;
          j_low = 0;
          j_high = n-1;
          while (1) {
            //printf("j_low=%i    j=%i    j_high=%i    x[j]=%f    x[j+1]=%f    oops=%i\n", j_low, j, j_high, x[j], x[j+1], (j_low > j_high));
            if (j_low > j_high) break;
            if (x_cur < x[j]) {
              j_high = j;
              j = j_low + (j_high - j_low)/2;
            } else if (x_cur <= x[j+1]) {
              break;
            } else {
              j_low = j+1;
              j = j_low + (j_high - j_low)/2;
            }
          }
        }
        //printf("final j=%i\n", j);
        // do the interpolation
        {
          double u, h, z;
          u = x_cur - x[j];
          h = x[j+1]-x[j];
          z = y[j+1]-y[j];
          new_y[i] = y[j] + u*(d[j] +
            u/h*(3*z/h - 2*d[j] - d[j+1] +
              u/h*(-2*z/h + d[j] + d[j+1])
            )
          );
        }
    }
}

double spline_int(double *x, double *y, double *d, int n) {
  int i;
  double result, h, z;
  result = 0;
  for (i=0; i<n-1; i++) {
    h = x[i+1]-x[i];
    z = y[i+1]-y[i];
    result += h*(y[i] + d[i]*h/2.0);
    result += h*(3*z-h*(2*d[i]+d[i+1]))/3.0;
    result += h*(-2*z+h*(d[i]+d[i+1]))/4.0;
  }
  return result;
}

void spline_cumul_int(double *x, double *y, double *d, double* yint, int n) {
  int i;
  double result, h, z;
  result = 0;
  for (i=0; i<n-1; i++) {
    yint[i] = result;
    h = x[i+1]-x[i];
    z = y[i+1]-y[i];
    result += h*(y[i] + d[i]*h/2.0);
    result += h*(3*z-h*(2*d[i]+d[i+1]))/3.0;
    result += h*(-2*z+h*(d[i]+d[i+1]))/4.0;
  }
  yint[n-1] = result;
}

