! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/scf/mixer/broyden.f90
!> Provides an electronic mixer implementation

!> Implementation of a modified Broyden mixing
module tblite_scf_mixer_broyden
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_lapack, only : getrf, getrs
   use tblite_scf_mixer_type, only : mixer_type
   use iso_c_binding
   implicit none
   private

   public :: new_broyden


   !> Configuration for the Broyden mixer
   type, public :: broyden_input
      !> Number of steps to keep in memory
      integer :: memory
      !> Damping parameter
      real(wp) :: damp
   end type broyden_input

   !> Electronic mixer using modified Broyden scheme
   type, public, extends(mixer_type) :: broyden_mixer
      integer :: ndim
      integer :: memory
      integer :: iter
      integer :: iset
      integer :: idif
      integer :: iget
      real(wp) :: damp
      real(wp), allocatable :: df(:, :)
      real(wp), allocatable :: u(:, :)
      real(wp), allocatable :: a(:, :)
      real(wp), allocatable :: dq(:)
      real(wp), allocatable :: dqlast(:)
      real(wp), allocatable :: qlast_in(:)
      real(wp), allocatable :: omega(:)
      real(wp), allocatable :: q_in(:)
   contains
      !> Apply mixing to the density
      procedure :: next
      !> Set new density from 1D array
      procedure :: set_1d
      !> Set difference between new and old density from 1D array
      procedure :: diff_1d
      !> Get density as 1D array
      procedure :: get_1d
      !> Get error metric from mixing
      procedure :: get_error
   end type broyden_mixer

   interface
   subroutine c_broyden(n,iter,memory,q,qlast,dq,dqlast,df,u,a,omega,alpha) bind(C, name="c_broyden")
      use iso_c_binding
      implicit none
      integer(c_int), value :: n
      integer(c_int), value :: iter
      integer(c_int), value :: memory
      real(c_double) :: q(*)
      real(c_double) :: qlast(*)
      real(c_double) :: dq(*)
      real(c_double)  :: dqlast(*)
      real(c_double) :: df(*)
      real(c_double) :: u(*)
      real(c_double) :: a(*)
      real(c_double) :: omega(*)
      real(c_double), value :: alpha
   end subroutine
   pure subroutine c_get_error(n,dq,error) bind(C, name="c_get_broyden_error")
      use iso_c_binding
      implicit none
      integer(c_int), value, intent(in) :: n
      real(c_double), intent(in) :: dq(*)
      real(c_double), intent(out) :: error
   end subroutine

end interface

contains

!> Create new instance of electronic mixer
subroutine new_broyden(self, ndim, input)
   !> Instance of the mixer
   type(broyden_mixer), intent(out) :: self
   !> Number of variables to consider
   integer, intent(in) :: ndim
   !> Configuration of the Broyden mixer
   type(broyden_input), intent(in) :: input

   self%ndim = ndim
   self%memory = input%memory
   self%iter = 0
   self%iset = 0
   self%idif = 0
   self%iget = 0
   self%damp = input%damp
   allocate(self%df(ndim, input%memory))
   allocate(self%u(ndim, input%memory))
   allocate(self%a(input%memory, input%memory))
   allocate(self%dq(ndim))
   allocate(self%dqlast(ndim))
   allocate(self%qlast_in(ndim))
   allocate(self%omega(input%memory))
   allocate(self%q_in(ndim))
end subroutine new_broyden

!> Set new density from 1D array
subroutine set_1d(self, qvec)
   !> Instance of the mixer
   class(broyden_mixer), intent(inout) :: self
   !> Density vector
   real(wp), intent(in) :: qvec(:)
   self%q_in(self%iset+1:self%iset+size(qvec)) = qvec
   self%iset = self%iset + size(qvec)
end subroutine set_1d

!> Set difference between new and old density from 1D array
subroutine diff_1d(self, qvec)
   !> Instance of the mixer
   class(broyden_mixer), intent(inout) :: self
   !> Density vector
   real(wp), intent(in) :: qvec(:)
   self%dq(self%idif+1:self%idif+size(qvec)) = qvec &
      & - self%q_in(self%idif+1:self%idif+size(qvec))
   self%idif = self%idif + size(qvec)
end subroutine diff_1d

!> Apply mixing to the density
subroutine next(self, error)
   !> Instance of the mixer
   class(broyden_mixer), intent(inout) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: info

   self%iset = 0
   self%idif = 0
   self%iget = 0
   self%iter = self%iter + 1

   if (self%iter == 1) then
      self%dqlast(:) = self%dq
      self%qlast_in(:) = self%q_in
      self%q_in(:) = self%q_in + self%damp * self%dq
      return
   end if

   call c_broyden(self%ndim,self%iter,self%memory,self%q_in,self%qlast_in,self%dq,self%dqlast,self%df,self%u,self%a,self%omega,self%damp)

end subroutine next

!> Get density as 1D array
subroutine get_1d(self, qvec)
   !> Instance of the mixer
   class(broyden_mixer), intent(inout) :: self
   !> Density vector
   real(wp), intent(out) :: qvec(:)
   qvec(:) = self%q_in(self%iget+1:self%iget+size(qvec))
   self%iget = self%iget + size(qvec)
end subroutine get_1d

pure function get_error(self) result(error)
   class(broyden_mixer), intent(in) :: self
   real(wp) :: error
   ! integer :: i
   ! error = 0.0_wp
   ! do i = 1, size(self%dq)
   !    error = error + self%dq(i)**2 / size(self%dq)
   ! end do
   ! error = sqrt(error)
   call c_get_error(self%ndim,self%dq,error)
end function get_error

end module tblite_scf_mixer_broyden
