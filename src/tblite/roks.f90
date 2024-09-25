! ! This file is part of tblite.
! ! SPDX-Identifier: LGPL-3.0-or-later
! !
! ! tblite is free software: you can redistribute it and/or modify it under
! ! the terms of the GNU Lesser General Public License as published by
! ! the Free Software Foundation, either version 3 of the License, or
! ! (at your option) any later version.
! !
! ! tblite is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! ! GNU Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public License
! ! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

! !> @dir tblite/roks
! !> Contains the implementation for the Restricted Open Shell Kohn Sham method

! !> @file tblite/roks.f90
! !> Reexports the ROKS related functionality.

! !> Proxy module for reexports from ROKS related modules
! module tblite_roks
!     use tblite_wavefunction_guess, only : sad_guess, eeq_guess, shell_partition
!     use tblite_wavefunction_mulliken, only : get_molecular_dipole_moment, &
!        & get_molecular_quadrupole_moment, get_mayer_bond_orders
!     use tblite_wavefunction_spin, only : magnet_to_updown, updown_to_magnet
!     use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction, &
!        & get_density_matrix, get_alpha_beta_occupation
!     implicit none
!     private

!     public :: roks_singlepoint
 
!  end module tblite_roks

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

!> @dir tblite/xtb
!> Routines related to the extend tight-binding Hamiltonian

!> @file tblite/xtb.f90
!> Proxy module for extended tight-binding Hamiltonian related settings

!> This module contains reexports for the #tblite_roks_calculator::roks_calculator class
!> and the possible parametrizations via #tblite_roks_gfn1, #tblite_roks_gfn2, and
!> #tblite_roks_ipea1. The main entry point for performing calculations is
!> provided with #tblite_roks_singlepoint::roks_singlepoint.
module tblite_roks
    use tblite_roks_singlepoint, only : roks_singlepoint
    implicit none
 
    public :: roks_singlepoint
 end module tblite_roks