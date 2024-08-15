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
   use mctc_env, only : wp
   use tblite_basis_type, only : basis_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_wavefunction, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_scf_info, only : scf_info
   use iso_c_binding
   implicit none

   interface
      type(c_ptr) function new_broyden(ndim, memory, alpha) bind(C,name="SetupBroyden")
         use iso_c_binding
         integer(c_int), value :: ndim
         integer(c_int), value :: memory
         real(c_double), value :: alpha
      end function new_broyden
   end interface

contains

!> Create a new instance of the mixer
   subroutine new_mixer(mixer, type, mol, calc, wfn, info)
      !> Pointer to the mixer on exit
      type(c_ptr), intent(out) :: mixer
      !> Type of mixer to use (Broyden=0, DIIS=1)
      integer, intent(in) :: type
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Single-point calculator
      type(xtb_calculator), intent(in) :: calc
      !> Wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Info data
      type(scf_info) :: info

      integer :: ndim

      select case(type)
       case(0)
         mixer = new_broyden(get_mixer_dimension(mol,calc%bas,info), calc%max_iter, calc%mixer_damping)
       case(1)
         write(*,*) "DIIS mixer"
         ! mixer = new_diis(10, calc%max_iter, calc%mixer_damping)
         stop
      end select
   end subroutine new_mixer

!> Get the dimensionality of the mixer
   function get_mixer_dimension(mol, bas, info) result(ndim)
      use tblite_scf_info, only : atom_resolved, shell_resolved
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
   end function get_mixer_dimension

end module tblite_scf_mixer
