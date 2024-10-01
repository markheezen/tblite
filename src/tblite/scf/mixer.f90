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

!> @dir tblite/scf/mixer
!> Routines for implementing electronic mixing

!> @file tblite/scf/mixer.f90
!> Proxy module for electronic mixing routines

!> Provides an electronic mixer implementation
module tblite_scf_mixer
   use tblite_xtb_calculator, only : xtb_calculator
   use mctc_io, only : structure_type
   use tblite_scf_info, only : scf_info
   use tblite_integral_type, only : integral_type
   use tblite_basis_type, only : basis_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_wavefunction, only : wavefunction_type
   use iso_c_binding
   implicit none



   !> Electronic mixer
   type, public, abstract :: mixer_type
      integer :: ndim
      integer :: memory
      type(c_ptr) :: ptr
   contains
      !> Apply mixing to the density
      procedure :: next
      !> Set new object to mix
      procedure(set), deferred :: set
      !> Set difference between two consecutive objects to mix
      procedure(diff), deferred :: diff
      !> Get mixed object
      procedure(get), deferred :: get
   end type mixer_type

   abstract interface
   subroutine set(self, iscf, wfn, info, ints)
      import :: mixer_type, wavefunction_type, scf_info, integral_type
      !> Instance of the mixer
      class(mixer_type), intent(inout) :: self
      !> Current iteration count
      integer, intent(inout) :: iscf
      !> Wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Info data
      type(scf_info) :: info
      !> Integral container
      type(integral_type), intent(in) :: ints
   end subroutine set

   subroutine diff(self, wfn, info)
      import :: mixer_type, wavefunction_type, scf_info
      !> Instance of the mixer
      class(mixer_type), intent(inout) :: self
      !> Wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Info data
      type(scf_info) :: info
   end subroutine diff

   subroutine get(self, bas, wfn, info)
      import :: mixer_type, basis_type, wavefunction_type, scf_info
      !> Instance of the mixer
      class(mixer_type), intent(inout) :: self
      !> Basis functions data
      type(basis_type), intent(in) :: bas      
      !> Wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Info data
      type(scf_info) :: info
   end subroutine get
   end interface

   interface
   subroutine next_mixer(mixer,iter) bind(C,name="Next")
      use iso_c_binding
      type(c_ptr), value, intent(in) :: mixer
      integer(c_int), value, intent(in) :: iter
   end subroutine next_mixer
   end interface
 
   contains

   subroutine next(self, iter)
      !> Mixer object
      class(mixer_type), intent(inout) :: self
      !> SCF Iteration
      integer, intent(in) :: iter

      call next_mixer(self%ptr, iter)
   end subroutine next

end module tblite_scf_mixer
