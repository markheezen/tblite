module tblite_purification_solver
   use mctc_env, only : sp, dp, error_type, wp
   use tblite_scf_solver, only : solver_type
   use tblite_lapack_solver, only : lapack_solver
   use tblite_cusolver_sygvd, only : new_sygvd_gpu, sygvd_cusolver
   use tblite_lapack_sygvd, only : new_sygvd, sygvd_solver
   use tblite_lapack_solver, only : solver_type
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction_type, only : wavefunction_type
   use iso_c_binding
   implicit none
   private

   public :: new_purification

   public :: purification_precision, purification_type, purification_runmode

   type enum_purification_type
       !> no puridfication
       integer(c_size_t) :: none = 0
       !> TC2 purification
       integer(c_size_t) :: tc2 = 2
       !> accelerated SP2
       integer(c_size_t) :: tc2accel = 42
       !> Mc Weeney Purification
       integer(c_size_t) :: mcweeney = 3
       !> TRS4 purification
       integer(c_size_t) :: trs4 = 4
   end type

   type enum_precision
       !> mixed precision
       integer(c_size_t) :: mixed = 0
       !> full double precision
       integer(c_size_t) :: double = 1
       !> full single precision
       integer(c_size_t) :: single = 2
       !> approx FP64
       integer(c_size_t) :: approx = 3
   end type

   type enum_runmode
       !> default switching CPU and GPU
       integer(c_size_t) :: default = 0
       !> all on the CPU
       integer(c_size_t) :: cpu = 1
       !> all GPU
       integer(c_size_t) :: gpu = 2
   end type
   
  type(enum_precision), parameter :: purification_precision = enum_precision()
  type(enum_purification_type), parameter :: purification_type = enum_purification_type()
  type(enum_runmode), parameter :: purification_runmode = enum_runmode()


   type, public, extends(solver_type) :: purification_solver
     integer(c_size_t) :: precision = purification_precision%mixed
     integer(c_size_t) :: type = purification_type%tc2
     integer(c_size_t) :: runmode = purification_runmode%default
     real(c_double) :: thresh = 1.0e-07
     type(c_ptr) :: solver_ptr
     integer(c_size_t) :: maxiter = 200
     integer :: iscf = 0
     class(solver_type), allocatable :: lapack_solv
     type(timer_type) :: timer
     logical :: trans = .false.
  contains
     procedure :: solve_sp
     procedure :: solve_dp
     procedure :: purify_dp
     procedure :: got_transform
     procedure :: delete
     procedure :: get_density_matrix
     procedure :: get_energy_w_density_matrix
  end type

  interface
     subroutine GetBoundsGershgorin(n, H, cpot, hmax, hmin) bind(C, name="GetBoundsGershgorin")
        use iso_c_binding
        integer(c_size_t), value, intent(in) :: n
        type(c_ptr), value, intent(in) :: H
        real(c_double), intent(out) :: cpot, hmax, hmin
     end subroutine
     subroutine GetDensityGuessAO(nao, nel, Fock, Overlap, Dens) bind(C, name="GetDensityGuessAO")
        use iso_c_binding
        integer(c_size_t), value, intent(in) :: nao
        real(c_double), value, intent(in) :: nel
        type(c_ptr), value, intent(in) :: Fock
        type(c_ptr), value, intent(in) :: Overlap
        type(c_ptr), value :: Dens
     end subroutine
     subroutine GetDensityGuess(nao, efermi, hmax, hmin, nel, Fock, Dens) bind(C, name="GetDensityGuess")
        use iso_c_binding
        integer(c_size_t), value, intent(in) :: nao
        real(c_double), value, intent(in) :: nel, efermi, hmax, hmin
        type(c_ptr), value, intent(in) :: Fock
        type(c_ptr), value :: Dens
     end subroutine
     subroutine PurifyDensity(nmo, nel, threshPP, maxit, Dens, type, run, prec) bind(C, name="PurifyDensity")
        use iso_c_binding
        integer(c_size_t), value, intent(in) :: nmo
        real(c_double), value, intent(in) :: nel
        real(c_double), value, intent(in) :: threshPP
        integer(c_size_t), value, intent(in) :: maxit
        type(c_ptr), value :: Dens
        integer(c_size_t), value, intent(in) :: type, run, prec
     end subroutine
     type(c_ptr) function PurifyFockSetup( nmo, threshPP, maxit, type, run, prec, prlevel) bind(C, name="PurifyFockSetup")
        use iso_c_binding
        integer(c_size_t), value, intent(in) :: nmo
        real(c_double), value, intent(in) :: threshPP
        integer(c_size_t), value, intent(in) :: maxit, prlevel
        integer(c_size_t), value, intent(in) :: type, run, prec
     end function
     subroutine SetTransformationMatrix(ptr, X) bind(C, name="SetTransformationMatrix")
        use iso_c_binding
        type(c_ptr), value :: ptr
        real(c_double),  intent(in) :: X(*)
     end subroutine
     subroutine SetOrbitalEnergiesMatrix(ptr, oen) bind(C, name="SetOrbitalEnergiesMatrix")
        use iso_c_binding
        type(c_ptr), value :: ptr
        real(c_double),  intent(in) :: oen(*)
     end subroutine
     subroutine GetDensityAO(ptr, Fock, Dens, elnum) bind(C, name="GetDensityAO")
        use iso_c_binding
        type(c_ptr), value :: ptr
        real(c_double) :: Fock(*)
        real(c_double) :: Dens(*)
        real(c_double), value :: elnum
     end subroutine
     subroutine GetEnergyWDensityAO(ptr, Fock, Dens, WDens) bind(C, name="GetEnergyWDensityAO")
        use iso_c_binding
        type(c_ptr), value :: ptr
        real(c_double) :: Fock(*)
        real(c_double) :: Dens(*)
        real(c_double) :: WDens(*)
     end subroutine
     subroutine DeletePointer(ptr) bind(C, name="DeletePointer")
        use iso_c_binding
        type(c_ptr), value :: ptr
     end subroutine
  end interface

  contains

  function got_transform(self) result(trans)
     class(purification_solver) :: self
     logical :: trans
     trans = self%trans;
  end function

  subroutine new_purification(self, type_, run_, prec_, ndim, maxiter, thresh)
     type(purification_solver) :: self
     integer(c_size_t) :: type_, run_, prec_
     integer :: ndim
     integer(c_size_t), optional :: maxiter
     real(c_double), optional :: thresh
     self%type = type_
     self%precision = prec_
     self%runmode = run_

     call self%timer%push("Setup LAPACK")
     !use LAPACK for molecules smaller than 750 basis functions
     select case(run_)
     case(purification_runmode%cpu)
        block
           type(sygvd_solver), allocatable :: tmp
           allocate(tmp) 
           call new_sygvd(tmp, ndim)
           call move_alloc(tmp, self%lapack_solv)
        end block
     case(purification_runmode%gpu)
        block
           type(sygvd_cusolver), allocatable :: tmp
           allocate(tmp) 
           call new_sygvd_gpu(tmp, ndim, .true.)
           call move_alloc(tmp, self%lapack_solv)
        end block
     case(purification_runmode%default)
        if (ndim < 750) then
           block
              type(sygvd_solver), allocatable :: tmp
              allocate(tmp) 
              call new_sygvd(tmp, ndim)
              call move_alloc(tmp, self%lapack_solv)
           end block
        else
           block
              type(sygvd_cusolver), allocatable :: tmp
              allocate(tmp) 
              call new_sygvd_gpu(tmp, ndim, .true.)
              call move_alloc(tmp, self%lapack_solv)
           end block
           end if
     end select


     call self%timer%pop()
     call self%timer%push("Setup TC")
     
     if (present(maxiter)) self%maxiter = maxiter
     if (present(thresh)) self%thresh = thresh

     self%solver_ptr =  PurifyFockSetup(int(ndim, kind=c_size_t), self%thresh, self%maxiter, self%type, self%runmode, self%precision, int(1, kind=c_size_t))
     call self%timer%pop()
  end subroutine

  subroutine solve_sp(self, hmat, smat, eval, error)
     class(purification_solver), intent(inout) :: self
     real(sp), contiguous, intent(inout) :: hmat(:, :)
     real(sp), contiguous, intent(in) :: smat(:, :)
     real(sp), contiguous, intent(inout) :: eval(:)
     type(error_type), allocatable, intent(out) :: error
  end subroutine



  subroutine solve_dp(self, hmat, smat, eval, error)
     class(purification_solver), intent(inout) :: self
     real(dp), contiguous, intent(inout) :: hmat(:, :)
     real(dp), contiguous, intent(in) :: smat(:, :)
     real(dp), contiguous, intent(inout) :: eval(:)
     real(c_double), allocatable, target :: tmp(:, :)
     type(error_type), allocatable, intent(out) :: error
     call self%timer%push("Diagonalize")
     call self%lapack_solv%solve(hmat, smat, eval, error)
     allocate(tmp(size(hmat, dim=1), size(hmat, dim=2)), source=0.0_c_double)
     tmp = hmat
     call SetTransformationMatrix(self%solver_ptr, tmp)
     call SetOrbitalEnergiesMatrix(self%solver_ptr, eval)
     call self%timer%pop()
     self%trans = .true.
  end subroutine solve_dp

  subroutine purify_dp(self, hmat, smat, dens, nel, error)
       class(purification_solver), intent(inout) :: self
       real(c_double), contiguous, target, intent(inout)  :: hmat(:, :)
       real(c_double), contiguous, target, intent(in) :: smat(:, :)
       real(c_double), contiguous, target, intent(inout) :: dens(:, :)
       real(c_double) :: nel
       type(error_type), allocatable, intent(out) :: error

        self%iscf = self%iscf +1
        call self%timer%push("Purification") 
        call GetDensityAO(self%solver_ptr, hmat, dens, nel/2)
        call self%timer%pop()
        
   end subroutine purify_dp

   subroutine delete(self)
     class(purification_solver) :: self
     call self%timer%push("Delete")
     call self%lapack_solv%delete() 
     deallocate(self%lapack_solv)
     call DeletePointer(self%solver_ptr)
     self%solver_ptr = c_null_ptr
     call self%timer%pop()
     block
        integer :: it
        real(wp) :: stime
        real(wp) :: ttime = 0.0_wp
        character(len=*), parameter :: label(*) = [character(len=20):: &
        & "Setup LAPACK", "Setup TC", "Purification", "TransformD", "Diagonalize", "Delete"]
        do it = 1, size(label)
           stime = self%timer%get(label(it))
           ttime = ttime + stime
           if (stime <= epsilon(0.0_wp)) cycle
           write(*,*) " - "//label(it)//format_time(stime)
        end do
        write(*,*) "___________________________________"             
        write(*,*) "   "//"Total: "//format_time(ttime)
     end block
   end subroutine

   subroutine get_density_matrix(self, focc, coeff, pmat)
     class(purification_solver) :: self
     real(wp), intent(in) :: focc(:)
     real(wp), contiguous, intent(in) :: coeff(:, :)
     real(wp), contiguous, intent(out) :: pmat(:, :)
  
     real(wp), allocatable :: scratch(:, :)
     integer :: iao, jao
     call self%lapack_solv%get_density_matrix(focc, coeff, pmat)
  end subroutine get_density_matrix


  subroutine get_energy_w_density_matrix(self, wfn, wdensity)
     class(purification_solver) :: self
     type(wavefunction_type), intent(inout) :: wfn
     real(wp) :: wdensity(:,:,:)
     real(wp), allocatable :: tmp(:)
     integer :: spin
     do spin = 1, wfn%nspin
        call GetEnergyWDensityAO(self%solver_ptr, wfn%coeff(:, :, spin), wfn%density(:, :, spin), wdensity(:, :, spin))         
     end do
  end subroutine


end module