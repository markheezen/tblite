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

!> @file tblite/xtb/singlepoint.f90
!> Provides main entry point for performing single point calculations with
!> a #roks_calculator instance.

!> Implementation of the single point calculation for a self-consistent
!> extended tight-binding Hamiltonian.

module tblite_roks_singlepoint
   use tblite_roks_context, only: roks_type, new_roks
   use mctc_env, only : wp, error_type, fatal_error, get_variable
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type, only : get_cutoff, basis_type
   use tblite_blas, only : gemv
   use tblite_container, only : container_cache
   use tblite_context, only : context_type, escape
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_lapack_sygvr, only : sygvr_solver
   use tblite_output_format, only : format_string
   use tblite_results, only : results_type
   use tblite_scf, only : scf_info, next_scf, &
   & potential_type, new_potential
   use tblite_scf_solver, only : solver_type
   use tblite_scf_mixer, only: mixer_type
   use tblite_scf_mixers, only: mixers_type, destroy, setup
   use tblite_scf_mixer_broyden, only: broyden_type, new_broyden
   use tblite_scf_mixer_diis, only: diis_type, new_diis
   use tblite_scf_potential, only : add_pot_to_h1
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction, only : wavefunction_type, get_density_matrix, &
   & get_alpha_beta_occupation, &
   & magnet_to_updown, updown_to_magnet
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, get_occupation, &
   & get_hamiltonian_gradient
   use tblite_post_processing_type, only : collect_containers_caches
   use tblite_post_processing_list, only : post_processing_list
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges, &
   & get_mulliken_atomic_multipoles
   use tblite_scf_iterator, only: get_qat_from_qsh
   use iso_c_binding
   implicit none
   private

   public :: roks_singlepoint

contains


!> Entry point for performing single point calculation using the xTB calculator
   subroutine roks_singlepoint(ctx, mol, calc, wfn, accuracy, energy, gradient, sigma, &
   & verbosity, results, post_process)
      !> Calculation context
      type(context_type), intent(inout) :: ctx
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Single-point calculator
      type(xtb_calculator), intent(in) :: calc
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
      
      type(error_type), allocatable :: error

      write(*,*) 
      call fatal_error(error, "This functionality is not implemented yet")
      call ctx%set_error(error)

   end subroutine roks_singlepoint

end module tblite_roks_singlepoint