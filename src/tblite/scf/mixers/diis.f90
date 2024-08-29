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

!> @file tblite/scf/mixers/diis.f90
!> Implementing DIIS (direct inversion in the iterative subspace) mixing

module tblite_scf_mixer_diis
   use mctc_env, only : wp
   use tblite_scf_mixer
   use tblite_basis_type, only : basis_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_wavefunction, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_scf_info, only : scf_info
   use tblite_integral_type, only : integral_type
   use iso_c_binding
   implicit none

   !> DIIS mixer
   type, extends(mixer_type) :: diis_type

   contains
      !> Set new object to mix
      procedure :: set => set_diis
        !> Set difference between two consecutive objects to mix
      procedure :: diff => diff_diis
      !> Get mixed object
      procedure :: get => get_diis
      !> Get error metric from mixing
    !   procedure :: get_error
   end type diis_type

   interface
      type(c_ptr) function c_new_diis(ndim, memory, alpha, nao) bind(C,name="SetupDIIS")
         use iso_c_binding
         integer(c_int), value :: ndim
         integer(c_int), value :: memory
         real(c_double), value :: alpha
         integer(c_int), value :: nao
      end function c_new_diis
      subroutine copy_matrices(mixer,iscf,fock,density,overlap) bind(C,name="CopyMatrices")
         use iso_c_binding
         use mctc_env, only : wp
         type(c_ptr), value, intent(in) :: mixer
         integer, intent(in), value :: iscf
         real(wp), intent(in) :: fock(*)
         real(wp), intent(in) :: density(*)
         real(wp), intent(in) :: overlap(*)
      end subroutine copy_matrices

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
   !> Create a new instance of the DIIS mixer
   subroutine new_diis(self, mol, calc, info)
      !> Broyden object
      class(diis_type), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Single-point calculator
      type(xtb_calculator), intent(in) :: calc
      !> Info data
      type(scf_info) :: info

      self%ndim = calc%bas%nao**2
      self%memory = 10
      self%ptr = c_new_diis(self%ndim, self%memory, calc%mixer_damping, calc%bas%nao)
   end subroutine new_diis

   !> Set the vector to mix
   subroutine set_diis(self, iscf, wfn, info, ints)
    use tblite_scf_info, only : atom_resolved, shell_resolved
    !> Instance of the Broyden mixer
    class(diis_type), intent(inout) :: self
    !> Current iteration count
    integer, intent(inout) :: iscf
    !> Wavefunction data
    type(wavefunction_type), intent(inout) :: wfn
    !> Info data
    type(scf_info) :: info
    !> Integral container
    type(integral_type), intent(in) :: ints

      call copy_matrices(self%ptr,iscf,wfn%coeff(:,:,1),wfn%density(:,:,1),ints%overlap)
   end subroutine set_diis

   !> Get the differences of the mixed vector
   subroutine diff_diis(self, wfn, info)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Info data
      type(scf_info) :: info

      call diff_mixer_data(self%ptr, wfn%coeff(:,:,1), size(wfn%coeff(:,:,1)))
   end subroutine diff_diis

!> Get the mixed vector
   subroutine get_diis(self, bas, wfn, info)
      use tblite_scf_info, only : atom_resolved, shell_resolved
      !> Instance of the DIIS mixer
      class(diis_type), intent(inout) :: self
      !> Basis functions data
      type(basis_type), intent(in) :: bas
      !> Wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Info data
      type(scf_info) :: info

      call get_mixer_data(self%ptr,wfn%coeff(:,:,1),size(wfn%coeff(:,:,1)))
   end subroutine get_diis

end module tblite_scf_mixer_diis
