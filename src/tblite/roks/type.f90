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

!> @file tblite/roks/type.f90
!> Provides a base class for the ROKS method

module tblite_roks_context
    use mctc_env, only : wp
    use tblite_wavefunction_type, only : wavefunction_type
    use tblite_container, only : container_cache, container_list
    use tblite_disp, only : dispersion_type
    use tblite_scf_potential, only: potential_type
    use tblite_xtb_calculator, only: xtb_calculator
    use tblite_xtb_coulomb, only : tb_coulomb
    use tblite_integral_type, only : integral_type
    implicit none
    private
 
    public :: new_roks

    type, public :: roks_type

 end type roks_type

contains


subroutine new_roks(self, spin, wfn, ccache, dcache, icache, hcache, rcache, ecache, pot, calc, ints)
    type(roks_type), intent(out) :: self
    integer, intent(in) :: spin ! 1=alpha, 2=beta
    type(wavefunction_type), intent(in) :: wfn
    type(container_cache), intent(in) :: ccache
    type(container_cache), intent(in) :: dcache
    type(container_cache), intent(in) :: icache
    type(container_cache), intent(in) :: hcache
    type(container_cache), intent(in) :: rcache
    type(container_cache), intent(in) :: ecache
    type(potential_type), intent(in) :: pot
    type(xtb_calculator), intent(in) :: calc
    type(integral_type), intent(in) :: ints

end subroutine new_roks


end module tblite_roks_context