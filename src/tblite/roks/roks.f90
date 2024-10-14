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

!> @file tblite/roks/roks.f90
!> Provides main entry point for performing a Restricted Open Shell Kohn
!> Sham calculation.

!> Implementation of the single point calculation for a self-consistent
!> extended tight-binding Hamiltonian.

module tblite_roks_singlepoint
   use mctc_env, only : error_type, fatal_error
   use tblite_context, only : context_type
   implicit none
   private

   public :: roks_singlepoint

contains


!> Entry point for performing the ROKS calculation
   subroutine roks_singlepoint(ctx)
      !> Calculation context
      type(context_type), intent(inout) :: ctx
      
      type(error_type), allocatable :: error

      write(*,*) 
      call fatal_error(error, "This functionality is not implemented yet")
      call ctx%set_error(error)

   end subroutine roks_singlepoint

end module tblite_roks_singlepoint