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


  subroutine gint2_nai_pF_pF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pF_pF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(49)
  end subroutine gint2_nai_pF_pF

  subroutine gint2_nai_pD_pF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pD_pF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(35)
  end subroutine gint2_nai_pD_pF

  subroutine gint2_nai_SP_pF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_SP_pF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(28)
  end subroutine gint2_nai_SP_pF

  subroutine gint2_nai_S_pF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_S_pF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(7)
  end subroutine gint2_nai_S_pF

  subroutine gint2_nai_P_pF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_P_pF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(21)
  end subroutine gint2_nai_P_pF

  subroutine gint2_nai_cD_pF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cD_pF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(42)
  end subroutine gint2_nai_cD_pF

  subroutine gint2_nai_cF_pF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cF_pF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(70)
  end subroutine gint2_nai_cF_pF

  subroutine gint2_nai_pF_pD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pF_pD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(35)
  end subroutine gint2_nai_pF_pD

  subroutine gint2_nai_pD_pD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pD_pD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(25)
  end subroutine gint2_nai_pD_pD

  subroutine gint2_nai_SP_pD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_SP_pD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(20)
  end subroutine gint2_nai_SP_pD

  subroutine gint2_nai_S_pD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_S_pD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(5)
  end subroutine gint2_nai_S_pD

  subroutine gint2_nai_P_pD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_P_pD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(15)
  end subroutine gint2_nai_P_pD

  subroutine gint2_nai_cD_pD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cD_pD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(30)
  end subroutine gint2_nai_cD_pD

  subroutine gint2_nai_cF_pD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cF_pD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(50)
  end subroutine gint2_nai_cF_pD

  subroutine gint2_nai_pF_SP(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pF_SP
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(28)
  end subroutine gint2_nai_pF_SP

  subroutine gint2_nai_pD_SP(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pD_SP
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(20)
  end subroutine gint2_nai_pD_SP

  subroutine gint2_nai_SP_SP(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_SP_SP
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(16)
  end subroutine gint2_nai_SP_SP

  subroutine gint2_nai_S_SP(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_S_SP
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(4)
  end subroutine gint2_nai_S_SP

  subroutine gint2_nai_P_SP(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_P_SP
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(12)
  end subroutine gint2_nai_P_SP

  subroutine gint2_nai_cD_SP(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cD_SP
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(24)
  end subroutine gint2_nai_cD_SP

  subroutine gint2_nai_cF_SP(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cF_SP
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(40)
  end subroutine gint2_nai_cF_SP

  subroutine gint2_nai_pF_S(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pF_S
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(7)
  end subroutine gint2_nai_pF_S

  subroutine gint2_nai_pD_S(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pD_S
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(5)
  end subroutine gint2_nai_pD_S

  subroutine gint2_nai_SP_S(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_SP_S
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(4)
  end subroutine gint2_nai_SP_S

  subroutine gint2_nai_S_S(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_S_S
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(1)
  end subroutine gint2_nai_S_S

  subroutine gint2_nai_P_S(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_P_S
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(3)
  end subroutine gint2_nai_P_S

  subroutine gint2_nai_cD_S(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cD_S
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(6)
  end subroutine gint2_nai_cD_S

  subroutine gint2_nai_cF_S(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cF_S
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(10)
  end subroutine gint2_nai_cF_S

  subroutine gint2_nai_pF_P(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pF_P
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(21)
  end subroutine gint2_nai_pF_P

  subroutine gint2_nai_pD_P(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pD_P
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(15)
  end subroutine gint2_nai_pD_P

  subroutine gint2_nai_SP_P(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_SP_P
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(12)
  end subroutine gint2_nai_SP_P

  subroutine gint2_nai_S_P(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_S_P
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(3)
  end subroutine gint2_nai_S_P

  subroutine gint2_nai_P_P(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_P_P
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(9)
  end subroutine gint2_nai_P_P

  subroutine gint2_nai_cD_P(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cD_P
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(18)
  end subroutine gint2_nai_cD_P

  subroutine gint2_nai_cF_P(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cF_P
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(30)
  end subroutine gint2_nai_cF_P

  subroutine gint2_nai_pF_cD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pF_cD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(42)
  end subroutine gint2_nai_pF_cD

  subroutine gint2_nai_pD_cD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pD_cD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(30)
  end subroutine gint2_nai_pD_cD

  subroutine gint2_nai_SP_cD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_SP_cD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(24)
  end subroutine gint2_nai_SP_cD

  subroutine gint2_nai_S_cD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_S_cD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(6)
  end subroutine gint2_nai_S_cD

  subroutine gint2_nai_P_cD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_P_cD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(18)
  end subroutine gint2_nai_P_cD

  subroutine gint2_nai_cD_cD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cD_cD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(36)
  end subroutine gint2_nai_cD_cD

  subroutine gint2_nai_cF_cD(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cF_cD
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(60)
  end subroutine gint2_nai_cF_cD

  subroutine gint2_nai_pF_cF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pF_cF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(70)
  end subroutine gint2_nai_pF_cF

  subroutine gint2_nai_pD_cF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_pD_cF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(50)
  end subroutine gint2_nai_pD_cF

  subroutine gint2_nai_SP_cF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_SP_cF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(40)
  end subroutine gint2_nai_SP_cF

  subroutine gint2_nai_S_cF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_S_cF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(10)
  end subroutine gint2_nai_S_cF

  subroutine gint2_nai_P_cF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_P_cF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(30)
  end subroutine gint2_nai_P_cF

  subroutine gint2_nai_cD_cF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cD_cF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(60)
  end subroutine gint2_nai_cD_cF

  subroutine gint2_nai_cF_cF(a, a_a, b, b_a, c, out)
    intent(c) gint2_nai_cF_cF
    intent(c)
    double precision intent(in) :: a(3), a_a, b(3), b_a, c(3)
    double precision intent(inout) :: out(100)
  end subroutine gint2_nai_cF_cF
