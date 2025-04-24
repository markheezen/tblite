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
!> Provides main entry point for performing a Restricted Open Shell Kohn-
!> Sham calculation.

!> Implementation of the ROKS calculation using a self-consistent
!> extended tight-binding Hamiltonian.

module tblite_roks_singlepoint
    use mctc_env, only : wp, error_type, fatal_error
    use mctc_io, only : structure_type
    use tblite_container, only : container_cache
    use tblite_integral_type, only : integral_type
    use tblite_output_format, only : format_string
    use tblite_context, only : context_type
    use tblite_output_format, only : format_string
    use tblite_results, only : results_type
    use tblite_scf, only : potential_type, new_potential
    use tblite_scf_solver, only : solver_type
    use tblite_scf_mixer, only: mixer_type
    use tblite_scf_mixers, only: mixers_type, destroy, setup
    use tblite_timer, only : timer_type, format_time
    use tblite_wavefunction, only : wavefunction_type
    use tblite_xtb_calculator, only : xtb_calculator
    use tblite_post_processing_list, only : post_processing_list
    use tblite_xtb_singlepoint, only : xtb_singlepoint
    use tblite_blas
 
    use tblite_basis_type, only : get_cutoff
    use tblite_cutoff, only : get_lattice_points
    use tblite_adjlist, only : adjacency_list, new_adjacency_list
    use tblite_integral_type, only : new_integral
    use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, get_occupation, &
       & get_hamiltonian_gradient
 
    use tblite_debug
    implicit none
    private
 
    public :: roks_singlepoint
 
 contains
 
    !> Entry point for performing single point calculation using the ROKS methodology
    subroutine roks_singlepoint(ctx, mol, calc, wfn, accuracy, energy, gradient, sigma, &
    & verbosity, results, post_process, roks_start)
       !> Calculation context
       type(context_type), intent(inout) :: ctx
       !> Molecular structure data
       type(structure_type), intent(inout) :: mol
       !> Single-point calculator
       type(xtb_calculator), intent(inout) :: calc
       !> Wavefunction data
       type(wavefunction_type), intent(inout) :: wfn
       !> Accuracy for computation
       real(wp), intent(in) :: accuracy
       !> Total energy
       real(wp), intent(out) :: energy
       !> Gradient with respect to cartesian coordinates
       real(wp), contiguous, intent(out), optional :: gradient(:, :)
       !> Strain derivatives with respect to strain deformations
       real(wp), contiguous, intent(out), optional :: sigma(:, :)
       !> Verbosity level of output
       integer, intent(in), optional :: verbosity
       !> Container for storing additional results
       type(results_type), intent(out), optional :: results
       type(post_processing_list), intent(inout), optional :: post_process
       !> Initial guess for ROKS calculation
       integer, intent(in) :: roks_start
 
       type(potential_type) :: pot
       logical :: econverged = .false., pconverged = .false., converged = .false.
       integer :: prlevel
       real(wp) :: pconv, econv
       class(solver_type), allocatable :: solver
       type(mixers_type) :: mixer
       type(error_type), allocatable :: error
       real(wp), allocatable :: eelec(:), eao(:)
       !> Container for storing results from first iteration
       type(results_type) :: results_first
 
       integer :: iscf
       real(wp) :: err
       type(timer_type) :: timer
       type(integral_type) :: ints
 
       real(wp), dimension(calc%bas%nao, calc%bas%nao) :: temp1, temp2, temp3
 
       write(*,*) "ROKS keyword succesfully used"
 
    end subroutine roks_singlepoint
 
 end module tblite_roks_singlepoint
 