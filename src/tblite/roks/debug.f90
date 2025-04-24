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

!> @file tblite/roks/debug.f90
!> Provides routines to easily print matrices and vectors for debugging
!> purposes.

module tblite_debug
    use mctc_env, only: wp

    public

    integer, parameter :: matrix_size = 6

contains

    subroutine print_mat(mat)
        real(wp), intent(in) :: mat(matrix_size,matrix_size)

        integer :: i

        do i = 1, matrix_size
        write(*,'(16F10.6)') mat(i,:)
        end do

    end subroutine print_mat

    subroutine print_vec(vec)
        real(wp), intent(in) :: vec(matrix_size)

        write(*,'(16F10.6)') vec(:)

    end subroutine print_vec

end module tblite_debug
