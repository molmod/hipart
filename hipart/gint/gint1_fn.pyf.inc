! HiPart is a software toolkit to analyse molecular densities with the hirshfeld partitioning scheme.
! Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>
!
! This file is part of HiPart.
!
! HiPart is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 3
! of the License, or (at your option) any later version.
!
! HiPart is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, see <http://www.gnu.org/licenses/>
!
! --


  subroutine gint1_fn_pF(a, a_a, p, out)
    intent(c) gint1_fn_pF
    intent(c)
    double precision intent(in) :: a(3), a_a, p(3)
    double precision intent(inout) :: out(7)
  end subroutine gint1_fn_pF

  subroutine gint1_fn_pD(a, a_a, p, out)
    intent(c) gint1_fn_pD
    intent(c)
    double precision intent(in) :: a(3), a_a, p(3)
    double precision intent(inout) :: out(5)
  end subroutine gint1_fn_pD

  subroutine gint1_fn_SP(a, a_a, p, out)
    intent(c) gint1_fn_SP
    intent(c)
    double precision intent(in) :: a(3), a_a, p(3)
    double precision intent(inout) :: out(4)
  end subroutine gint1_fn_SP

  subroutine gint1_fn_S(a, a_a, p, out)
    intent(c) gint1_fn_S
    intent(c)
    double precision intent(in) :: a(3), a_a, p(3)
    double precision intent(inout) :: out(1)
  end subroutine gint1_fn_S

  subroutine gint1_fn_P(a, a_a, p, out)
    intent(c) gint1_fn_P
    intent(c)
    double precision intent(in) :: a(3), a_a, p(3)
    double precision intent(inout) :: out(3)
  end subroutine gint1_fn_P

  subroutine gint1_fn_cD(a, a_a, p, out)
    intent(c) gint1_fn_cD
    intent(c)
    double precision intent(in) :: a(3), a_a, p(3)
    double precision intent(inout) :: out(6)
  end subroutine gint1_fn_cD

  subroutine gint1_fn_cF(a, a_a, p, out)
    intent(c) gint1_fn_cF
    intent(c)
    double precision intent(in) :: a(3), a_a, p(3)
    double precision intent(inout) :: out(10)
  end subroutine gint1_fn_cF

  integer function gint1_fn_basis(weights, fns, points, centers, shell_types, shell_map, num_primitives, ccoeffs, exponents, num_weights, num_points, num_centers, num_shells, num_ccoeffs, num_exponents)
    intent(c) gint1_fn_basis
    intent(c)
    double precision intent(in) :: weights(num_weights)
    double precision intent(inout) :: fns(num_points)
    double precision intent(in) :: points(num_points,3)
    double precision intent(in)  :: centers(num_centers,3)
    integer intent(int) :: shell_types(num_shells)
    integer intent(int) :: shell_map(num_shells)
    integer intent(int) :: num_primitives(num_shells)
    double precision intent(in) :: ccoeffs(num_ccoeffs)
    double precision intent(in) :: exponents(num_exponents)
    integer intent(hide), depend(weights) :: num_weights=len(weights)
    integer intent(hide), depend(points) :: num_points=len(points)
    integer intent(hide), depend(centers) :: num_centers=len(centers)
    integer intent(hide), depend(shell_types) :: num_shells=len(shell_types)
    integer intent(hide), depend(ccoeffs) :: num_ccoeffs=len(ccoeffs)
    integer intent(hide), depend(exponents) :: num_exponents=len(exponents)
  end function gint1_fn_basis

  integer function gint1_fn_dmat(dmat, density, points, centers, shell_types, shell_map, num_primitives, ccoeffs, exponents, num_dmat, num_points, num_centers, num_shells, num_ccoeffs, num_exponents)
    intent(c) gint1_fn_dmat
    intent(c)
    double precision intent(in) :: dmat(num_dmat)
    double precision intent(inout) :: density(num_points)
    double precision intent(in) :: points(num_points,3)
    double precision intent(in)  :: centers(num_centers,3)
    integer intent(int) :: shell_types(num_shells)
    integer intent(int) :: shell_map(num_shells)
    integer intent(int) :: num_primitives(num_shells)
    double precision intent(in) :: ccoeffs(num_ccoeffs)
    double precision intent(in) :: exponents(num_exponents)
    integer intent(hide), depend(dmat) :: num_dmat=len(dmat)
    integer intent(hide), depend(points) :: num_points=len(points)
    integer intent(hide), depend(centers) :: num_centers=len(centers)
    integer intent(hide), depend(shell_types) :: num_shells=len(shell_types)
    integer intent(hide), depend(ccoeffs) :: num_ccoeffs=len(ccoeffs)
    integer intent(hide), depend(exponents) :: num_exponents=len(exponents)
  end function gint1_fn_dmat

  integer function gint1_fn_overlap(overlap, points, weights, centers, shell_types, shell_map, num_primitives, ccoeffs, exponents, num_overlap, num_points, num_centers, num_shells, num_ccoeffs, num_exponents)
    intent(c) gint1_fn_overlap
    intent(c)
    double precision intent(inout) :: overlap(num_overlap)
    double precision intent(in) :: points(num_points,3)
    double precision intent(in) :: weights(num_points)
    double precision intent(in)  :: centers(num_centers,3)
    integer intent(int) :: shell_types(num_shells)
    integer intent(int) :: shell_map(num_shells)
    integer intent(int) :: num_primitives(num_shells)
    double precision intent(in) :: ccoeffs(num_ccoeffs)
    double precision intent(in) :: exponents(num_exponents)
    integer intent(hide), depend(overlap) :: num_overlap=len(overlap)
    integer intent(hide), depend(points) :: num_points=len(points)
    integer intent(hide), depend(centers) :: num_centers=len(centers)
    integer intent(hide), depend(shell_types) :: num_shells=len(shell_types)
    integer intent(hide), depend(ccoeffs) :: num_ccoeffs=len(ccoeffs)
    integer intent(hide), depend(exponents) :: num_exponents=len(exponents)
  end function gint1_fn_overlap
