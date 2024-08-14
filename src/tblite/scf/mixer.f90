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

   integer(c_int) function get_broyden_dimension(charge, dipole, quadrupole, nat, nsh) bind(C,name="GetBroydenDimension")
      use iso_c_binding
      integer(c_int), value :: charge
      integer(c_int), value :: dipole
      integer(c_int), value :: quadrupole
      integer(c_int), value :: nat
      integer(c_int), value :: nsh
   end function get_broyden_dimension

   subroutine next_broyden(mixer,iter) bind(C,name="NextBroyden")
      use iso_c_binding
      type(c_ptr), value, intent(in) :: mixer
      integer(c_int), value, intent(in) :: iter
  end subroutine next_broyden

   real(c_double) function get_broyden_error(mixer) bind(C,name="GetBroydenError")
      use iso_c_binding
      type(c_ptr), value :: mixer
   end function get_broyden_error

   subroutine destroy_broyden(mixer) bind(C,name="DestroyBroyden")
      use iso_c_binding
      type(c_ptr), value :: mixer
   end subroutine destroy_broyden
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
         ndim = wfn%nspin*get_broyden_dimension(info%charge, info%dipole, info%quadrupole, mol%nat, calc%bas%nsh)
         mixer = new_broyden(ndim, calc%max_iter, calc%mixer_damping)
      case(1)
         write(*,*) "DIIS mixer"
         stop
   end select
end subroutine new_mixer

!> Destruct the mixer
subroutine next_mixer(mixer, type, iter)
   !> Pointer to the mixer
   type(c_ptr), intent(in) :: mixer
   !> Type of mixer to use (Broyden=0, DIIS=1)
   integer, intent(in) :: type
   !> Iteration in the SCF process
   integer, intent(in) :: iter

   select case(type)
      case(0)
         call next_broyden(mixer,iter)
      case(1)
         write(*,*) "DIIS mixer"
         stop
   end select
end subroutine next_mixer

!> Get the density error
function get_error(mixer, type) result(error)
   !> Pointer to the mixer
   type(c_ptr), intent(out) :: mixer
   !> Type of mixer to use (Broyden=0, DIIS=1)
   integer, intent(in) :: type
   !> Density error
   double precision :: error

   select case(type)
      case(0)
         error=get_broyden_error(mixer)
      case(1)
         write(*,*) "DIIS mixer"
         stop
   end select
end function get_error

!> Destruct the mixer
subroutine destroy_mixer(mixer, type)
   !> Pointer to the mixer
   type(c_ptr), intent(out) :: mixer
   !> Type of mixer to use (Broyden=0, DIIS=1)
   integer, intent(in) :: type

   select case(type)
      case(0)
         call destroy_broyden(mixer)
      case(1)
         write(*,*) "DIIS mixer"
         stop
   end select
end subroutine destroy_mixer


end module tblite_scf_mixer