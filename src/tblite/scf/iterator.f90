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

!> @file tblite/scf/iterator.f90
!> Provides the implementation of the actual self-consistent field iteractions

!> Iterator for evaluating the Hamiltonian self-consistently
module tblite_scf_iterator
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_cache, container_list
   use tblite_disp, only : dispersion_type
   use tblite_integral_type, only : integral_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges, &
      & get_mulliken_atomic_multipoles
   use tblite_xtb_coulomb, only : tb_coulomb
   use tblite_scf_mixers, only: mixers_type
   use tblite_scf_info, only : scf_info
   use tblite_scf_potential, only : potential_type, add_pot_to_h1
   use tblite_scf_solver, only : solver_type
   use tblite_scf_utils, only : get_electronic_energy, reduce, get_qat_from_qsh, &
      & get_density
   use tblite_exchange_type, only : exchange_type
   implicit none
   private

   public :: next_scf

contains

!> Evaluate self-consistent iteration for the density-dependent Hamiltonian
subroutine next_scf(iscf, mol, bas, wfn, solver, mixer, info, coulomb, dispersion, &
      & interactions, exchange, ints, pot, cache, dcache, icache, ecache, &
      & energies, error)
   !> Current iteration count
   integer, intent(inout) :: iscf
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Solver for the general eigenvalue problem
   class(solver_type), intent(inout) :: solver
   !> Convergence accelerator
   class(mixers_type), intent(inout) :: mixer
   !> Information on wavefunction data used to construct Hamiltonian
   type(scf_info), intent(in) :: info
   !> Container for coulombic interactions
   type(tb_coulomb), intent(in), optional :: coulomb
   !> Container for dispersion interactions
   class(dispersion_type), intent(in), optional :: dispersion
   !> Container for general interactions
   type(container_list), intent(in), optional :: interactions
   !> Container for exchange interactions
   class(exchange_type), intent(in), optional :: exchange

   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Density dependent potential shifts
   type(potential_type), intent(inout) :: pot
   !> Restart data for coulombic interactions
   type(container_cache), intent(inout) :: cache
   !> Restart data for dispersion interactions
   type(container_cache), intent(inout) :: dcache
   !> Restart data for interaction containers
   type(container_cache), intent(inout) :: icache
   !> Restart data for Mulliken exchange
   type(container_cache), intent(inout) :: ecache

   !> Self-consistent energy
   real(wp), intent(inout) :: energies(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: eao(:)
   real(wp) :: ts

   !we do this to check whether a density-based guess was used
   if (wfn%density(1, 1, 1) .ne. 0.0_wp) then
      call get_mulliken_shell_charges(bas, ints%overlap, wfn%density, wfn%n0sh, &
      & wfn%qsh)
      call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
      call get_mulliken_atomic_multipoles(bas, ints%dipole, wfn%density, &
      & wfn%dpat)
      call get_mulliken_atomic_multipoles(bas, ints%quadrupole, wfn%density, &
      & wfn%qpat)
   end if
   
   iscf = iscf + 1

   if (allocated(mixer%broyden) .and. iscf > 1) then
      call mixer%broyden%next(iscf)
      call mixer%broyden%get(wfn%qat, wfn%qsh, wfn%dpat, wfn%qpat, info)
   end if 
   
   call pot%reset
   if (present(coulomb)) then
      call coulomb%get_potential(mol, cache, wfn, pot)
   end if
   if (present(dispersion)) then
      call dispersion%get_potential(mol, dcache, wfn, pot)
   end if
   if (present(interactions)) then
      call interactions%get_potential(mol, icache, wfn, pot)
   end if

   if (present(exchange) .and. iscf > 1) then 
      call exchange%get_potential_w_overlap(mol, ecache, wfn, pot, ints%overlap)
   end if
   call add_pot_to_h1(bas, ints, pot, wfn%coeff)
   
   if (allocated(mixer%broyden) .and. iscf > 1) then
      call mixer%broyden%set(wfn%qat, wfn%qsh, wfn%dpat, wfn%qpat, info)
   end if

   if (allocated(mixer%diis)) call mixer%diis%set(wfn%coeff(:,:,1))

   if (allocated(mixer%diis) .and. iscf > 1) then
      call mixer%diis%construct_error(wfn%coeff(:,:,1), wfn%density(:,:,1), ints%overlap)
      call mixer%diis%next(iscf)
      call mixer%diis%get(wfn%coeff(:,:,1))
   endif

   call get_density(wfn, solver, ints, ts, error)
   if (allocated(error)) return
   
   call get_mulliken_shell_charges(bas, ints%overlap, wfn%density, wfn%n0sh, &
      & wfn%qsh)
   call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
   
   call get_mulliken_atomic_multipoles(bas, ints%dipole, wfn%density, &
      & wfn%dpat)
   call get_mulliken_atomic_multipoles(bas, ints%quadrupole, wfn%density, &
      & wfn%qpat)

   if (allocated(mixer%broyden)) then 
      call mixer%broyden%diff(wfn%qat, wfn%qsh, wfn%dpat, wfn%qpat, info)
   end if

   if (allocated(mixer%diis) .and. iscf == 1) then
      call mixer%diis%diff(wfn%coeff(:,:,1))
   end if
   
   allocate(eao(bas%nao), source=0.0_wp)
   call get_electronic_energy(ints%hamiltonian, wfn%density, eao)

   energies(:) = ts / size(energies)
   call reduce(energies, eao, bas%ao2at)
   if (present(coulomb)) then
      call coulomb%get_energy(mol, cache, wfn, energies)
   end if
   if (present(dispersion)) then
      call dispersion%get_energy(mol, dcache, wfn, energies)
   end if
   if (present(interactions)) then
      call interactions%get_energy(mol, icache, wfn, energies)
   end if
   if (present(exchange) .and. iscf > 1) then 
      call exchange%get_energy(mol, ecache, wfn, energies)
   end if

end subroutine next_scf

end module tblite_scf_iterator
