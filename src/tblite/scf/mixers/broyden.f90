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

!> @file tblite/scf/mixers/broyden.f90
!> Implementing Broyden mixing

!> Provides an electronic mixer implementation
module tblite_scf_mixer_broyden
   use mctc_env, only : wp
   use tblite_scf_mixer
   use tblite_basis_type, only : basis_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_wavefunction, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_scf_info, only : scf_info
   use tblite_integral_type, only : integral_type
   use tblite_scf_utils, only: get_qat_from_qsh
   use iso_c_binding
   implicit none

!> Broyden mixer
   type, extends(mixer_type) :: broyden_type

   contains
      !> Get dimension of object to mix
      procedure :: get_dimension
      !> Set new object to mix
      procedure :: set => set_broyden
      !> Set difference between two consecutive objects to mix
      procedure :: diff => diff_broyden
      !> Get mixed object
      procedure :: get => get_broyden
      !> Get error metric from mixing
    !   procedure :: get_error
   end type broyden_type

   interface
      type(c_ptr) function c_new_broyden(ndim, memory, alpha, nao) bind(C,name="SetupBroyden")
         use iso_c_binding
         integer(c_int), value :: ndim
         integer(c_int), value :: memory
         real(c_double), value :: alpha
         integer(c_int), value :: nao
      end function c_new_broyden

      subroutine set_mixer_data(mixer,target,size) bind(C,name="SetData")
         use iso_c_binding
         use mctc_env, only : wp
         type(c_ptr), value, intent(in) :: mixer
         real(wp), intent(in) :: target(*)
         integer(c_int), value, intent(in) :: size
      end subroutine set_mixer_data

      subroutine get_mixer_data(mixer,target,size) bind(C,name="GetData")
        use iso_c_binding
        use mctc_env, only : wp
        type(c_ptr), value, intent(in) :: mixer
        real(wp), intent(inout) :: target(*)
        integer(c_int), value, intent(in) :: size
     end subroutine get_mixer_data

     subroutine diff_mixer_data(mixer,target,size) bind(C,name="DiffData")
        use iso_c_binding
        use mctc_env, only : wp
        type(c_ptr), value, intent(in) :: mixer
        real(wp), intent(in) :: target(*)
        integer(c_int), value, intent(in) :: size
     end subroutine diff_mixer_data
   end interface

contains

   !> Create a new instance of the Broyden mixer
   subroutine new_broyden(self, mol, calc, info, roks)
      !> Broyden object
      class(broyden_type), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Single-point calculator
      type(xtb_calculator), intent(in) :: calc
      !> Info data
      type(scf_info) :: info
      !> ROKS method
      logical, optional :: roks

      self%ndim = self%get_dimension(mol,calc%bas,info)
      if (present(roks)) self%ndim = self%ndim * 4
      self%memory = calc%mixer_mem
      self%ptr = c_new_broyden(self%ndim, self%memory, calc%mixer_damping, calc%bas%nao)
   end subroutine new_broyden

   !> Get the dimensionality of the mixer
   function get_dimension(self, mol, bas, info) result(ndim)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      !> Instance of the Broyden mixer
      class(broyden_type), intent(inout) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis functions
      type(basis_type), intent(in) :: bas
      !> Info data
      type(scf_info), intent(in) :: info
      integer :: ndim

      ndim = 0

      select case(info%charge)
       case(atom_resolved)
         ndim = ndim + mol%nat
       case(shell_resolved)
         ndim = ndim + bas%nsh
      end select

      select case(info%dipole)
       case(atom_resolved)
         ndim = ndim + 3*mol%nat
      end select

      select case(info%quadrupole)
       case(atom_resolved)
         ndim = ndim + 6*mol%nat
      end select
   end function get_dimension


   !> Set the vector to mix
   subroutine set_broyden(self, iscf, wfn, info, ints)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      !> Instance of the Broyden mixer
      class(broyden_type), intent(inout) :: self
      !> Current iteration count
      integer, intent(inout) :: iscf
      !> Wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Info data
      type(scf_info) :: info
      !> Integral container
      type(integral_type), intent(in) :: ints

      select case(info%charge)
       case(atom_resolved)
         call set_mixer_data(self%ptr, wfn%qat, size(wfn%qat))
       case(shell_resolved)
         call set_mixer_data(self%ptr, wfn%qsh, size(wfn%qsh))
      end select

      select case(info%dipole)
       case(atom_resolved)
         call set_mixer_data(self%ptr, wfn%dpat, size(wfn%dpat))
      end select

      select case(info%quadrupole)
       case(atom_resolved)
         call set_mixer_data(self%ptr, wfn%qpat, size(wfn%qpat))
      end select
   end subroutine set_broyden


!> Get the differences of the mixed vector
   subroutine diff_broyden(self, wfn, info)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      !> Instance of the Broyden mixer
      class(broyden_type), intent(inout) :: self
      !> Wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Info data
      type(scf_info) :: info

      select case(info%charge)
       case(atom_resolved)
         call diff_mixer_data(self%ptr, wfn%qat, size(wfn%qat))
       case(shell_resolved)
         call diff_mixer_data(self%ptr, wfn%qsh, size(wfn%qsh))
      end select

      select case(info%dipole)
       case(atom_resolved)
         call diff_mixer_data(self%ptr, wfn%dpat, size(wfn%dpat))
      end select

      select case(info%quadrupole)
       case(atom_resolved)
         call diff_mixer_data(self%ptr, wfn%qpat, size(wfn%qpat))
      end select
   end subroutine diff_broyden

!> Get the mixed vector
   subroutine get_broyden(self, bas, wfn, info)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      !> Instance of the Broyden mixer
      class(broyden_type), intent(inout) :: self
      !> Basis functions data
      type(basis_type), intent(in) :: bas
      !> Wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Info data
      type(scf_info) :: info


      select case(info%charge)
       case(atom_resolved)
         call get_mixer_data(self%ptr, wfn%qat, size(wfn%qat))
       case(shell_resolved)
         call get_mixer_data(self%ptr, wfn%qsh, size(wfn%qsh))
         call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
      end select

      select case(info%dipole)
       case(atom_resolved)
         call get_mixer_data(self%ptr, wfn%dpat, size(wfn%dpat))
      end select

      select case(info%quadrupole)
       case(atom_resolved)
         call get_mixer_data(self%ptr, wfn%qpat, size(wfn%qpat))
      end select
   end subroutine get_broyden

end module tblite_scf_mixer_broyden
