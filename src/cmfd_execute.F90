module cmfd_execute 

!==============================================================================
! CMFD_EXECUTE -- This module is the highest level cmfd module that controls the
! cross section generation, diffusion calculation, and source re-weighting
!==============================================================================

  implicit none
  private
  public :: execute_cmfd, cmfd_init_batch

# ifdef PETSC
#   include <finclude/petsc.h90>
# endif

contains

!==============================================================================
! EXECUTE_CMFD
!==============================================================================

  subroutine execute_cmfd()

# ifdef PETSC

    use cmfd_data,              only: set_up_cmfd
    use cmfd_message_passing,   only: petsc_init_mpi, cmfd_bcast
    use cmfd_power_solver,      only: cmfd_power_execute
    use cmfd_snes_solver,       only: cmfd_snes_execute
    use constants,              only: CMFD_EIG, MAX_FILE_LEN
    use error,                  only: warning, fatal_error 
    use global,                 only: n_procs_cmfd, cmfd,                       &
                                      cmfd_solver_type, time_cmfd,              &
                                      cmfd_run_adjoint, cmfd_write_hdf5,        &
                                      cmfd_feedback, cmfd_hold_weights,         &
                                      cmfd_inact_flush, cmfd_keff_tol,          &
                                      cmfd_act_flush, cmfd_multiset,            & 
                                      cmfd_set_size, current_batch, keff,       &
                                      n_batches, message, master, mpi_err,      &
                                      rank, cmfd_accum, n_inactive, path_output, &
                                      spec_rad_on
    ! BEGIN VARIABLES ADDED BY K.KEADY ON 12/10/2013 
    use string,                 only: to_str
    integer :: i                          ! iteration counter for x
    integer :: j                          ! iteration counter for y
    integer :: k                          ! iteration counter for z
    integer :: g                          ! iteration counter for groups
    integer :: nx                         ! number of mesh cells in x
    integer :: ny                         ! number of mesh cells in y
    integer :: nz                         ! number of mesh cells in z
    integer :: ng                         ! number of groups
    character(MAX_FILE_LEN) :: filename   ! output filename
    ! END VARIABLES ADDED
    
    ! stop cmfd timer
    if (master) then
      call time_cmfd % start()
    end if

    ! filter processors (lowest PETSc group)
    if (rank < n_procs_cmfd) then

      ! set up cmfd data (master only)
      if (master) call set_up_cmfd()

      ! broadcast cmfd to all petsc procs
      call cmfd_bcast()

      ! process solver options
      call process_cmfd_options()

    end if

    ! filter processors (lowest PETSc group)
    if (rank < n_procs_cmfd) then

      ! call solver
      if (trim(cmfd_solver_type) == 'power') then
        call cmfd_power_execute()
      elseif (trim(cmfd_solver_type) == 'jfnk') then
        call cmfd_snes_execute()
      else
        message = 'solver type became invalid after input processing'
        call fatal_error() 
      end if

      ! perform any last batch tasks 
      if (current_batch == n_batches) then

        ! check for adjoint run
        if (cmfd_run_adjoint) then
          if (trim(cmfd_solver_type) == 'power') then
            call cmfd_power_execute(adjoint = .true.)
          elseif (trim(cmfd_solver_type) == 'jfnk') then
            call cmfd_snes_execute(adjoint = .true.)
          end if
        end if
        
        !If accumulating, we need to store the final CMFD flux here...
        if(cmfd_accum) then
          cmfd % phi_final = cmfd % phi
        end if
      end if
      
      ! If calcing spectral radius, do it now!
      if (spec_rad_on) call calc_spec_rad()
      
      ! if not accumulating tallies, accum sum and sum sq
      ! for this batch
      if(current_batch > n_inactive .and. .not. cmfd_accum) then
        cmfd % phi_sum = cmfd % phi_sum + cmfd % phi
        cmfd % phi_sum_sq = cmfd % phi_sum_sq + (cmfd % phi * cmfd % phi)
      ! if using multiset method (prescribed flushes), do the same
      elseif(current_batch > n_inactive .and. cmfd_multiset) then
        cmfd % phi_sum = cmfd % phi_sum + cmfd % phi
        cmfd % phi_sum_sq = cmfd % phi_sum_sq + (cmfd % phi * cmfd % phi)
      end if
      
      ! BEGIN SECTION ADDED BY K.KEADY ON 12/10/2013 TO LOOK AT CYCLE-TO-CYCLE
      ! VARIATION IN LOW-ORDER EIGENFUNCTIONS
      ! open output file for this cycle
      !filename = trim(path_output) // trim(to_str(current_batch)) //&
      !           '_cyc_flux.out'
      !open(UNIT=CMFD_EIG, FILE=filename, ACTION='write')
      ! extract spatial and energy indices from object
      !nx = cmfd % indices(1)
      !ny = cmfd % indices(2)
      !nz = cmfd % indices(3)
      !ng = cmfd % indices(4)
      ! begin loop around space and energy groups
      !ZLOOP: do k = 1, nz
      !  XLOOP: do i = 1, nx
      !    YLOOP: do j = 1, ny
      !      GROUPG: do g = 1, ng
      !        write(CMFD_EIG,'(1X,I0,1X,I0,1X,I0,1X,I0,1X,E20.7)') &
      !              i,j,k,g,cmfd % phi(g + ng*(i-1) + ng*nx*(j-1) &
      !              + ng*nx*ny*(k-1))
      !    end do GROUPG
      !  end do YLOOP
      !end do XLOOP
    !end do ZLOOP
    ! close file
    !close(CMFD_EIG)
    ! END SECTION ADDED BY K.KEADY
    
    end if

    ! calculate fission source
    call calc_fission_source()

    ! calculate weight factors if using std feedback
    if (cmfd_feedback .and. .not. cmfd_multiset) then
      call cmfd_reweight(.true.)
    ! also calc wt factors if we're using the multiset method and this is NOT
    ! the first cycle in a set 
    elseif (cmfd_multiset) then
      if(mod(current_batch,cmfd_set_size)/= 1) then
        call cmfd_reweight(.true.)
      elseif(mod(current_batch,cmfd_set_size)== 1) then
        message = "Batch " // trim(to_str(current_batch)) // " feedback OFF" 
        call warning()
      end if
    end if
    
    ! stop cmfd timer
    if (master) then
      call time_cmfd % stop()
    end if

    ! wait here for all procs
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

# endif

  end subroutine execute_cmfd

!==============================================================================
! CMFD_INIT_BATCH
!==============================================================================

  subroutine cmfd_init_batch()

    use global,            only: cmfd_begin, cmfd_on, cmfd_tally_on,         &
                                 cmfd_inact_flush, cmfd_act_flush, cmfd_run, &
                                 current_batch, cmfd_hold_weights,           &
                                 cmfd_inactive, n_inactive, cmfd_cfl,        &
                                 cmfd_accum, n_batches, gen_per_batch, cmfd, &
                                 cmfd_multiset, cmfd_set_size, message
    use error,             only: warning

    ! initially set tally flush to F for false; will be changed if the
    ! tally reset function is called
    cmfd_cfl = 'F'
    
    ! Allocate k storage (if it's the first batch)
    if (.not. allocated(cmfd % k_gen_cmfd)) allocate(cmfd % k_gen_cmfd(n_batches*gen_per_batch))
    ! check to activate CMFD diffusion and possible feedback
    ! this guarantees that when cmfd begins at least one batch of tallies are
    ! accumulated
    if (cmfd_run .and. cmfd_begin == current_batch) then
      cmfd_on = .true.
      cmfd_tally_on = .true.
      ! allocate array for cmfd batch keff
    end if

    ! check to flush cmfd tallies for active batches, no more inactive flush
    if (cmfd_run .and. cmfd_act_flush == current_batch) then
      call cmfd_tally_reset()
      cmfd_tally_on = .true.
      cmfd_inact_flush(2) = -1
    end if

    ! check to flush cmfd tallies (>= on number of
    ! flushes important as the code will flush on the first batch which we
    ! dont want to count)
    
    ! check for different cases on which to flush tallies:
    if(cmfd_run .and. cmfd_begin < current_batch) then
      ! If accum=false, flush tallies every cycle (active or not)
      if (.not. cmfd_accum) then
        cmfd_hold_weights = .true.
        call cmfd_tally_reset()
      ! If using the multiset method... 
      elseif (cmfd_multiset) then
        ! flush every (cmfd_set_size) cycles
        if(mod(current_batch,cmfd_set_size)== 1) then
          cmfd_hold_weights = .true.
          call cmfd_tally_reset()
        end if
      ! If inactive=false, flush tallies every cycle during inactive batches 
      elseif (.not. cmfd_inactive .and. current_batch < n_inactive) then
        cmfd_hold_weights = .true.
        call cmfd_tally_reset()
      ! If this is a prescribed flush cycle, flush...
      elseif (mod(current_batch,cmfd_inact_flush(1))== 0 .and. &
          cmfd_inact_flush(2) > 0) then
        cmfd_hold_weights = .true.
        call cmfd_tally_reset()
        cmfd_inact_flush(2) = cmfd_inact_flush(2) - 1
      end if  
    end if


  end subroutine cmfd_init_batch

# ifdef PETSC

!==============================================================================
! PROCESS_CMFD_OPTIONS 
!==============================================================================

  subroutine process_cmfd_options()

    use global,       only: cmfd_snes_monitor, cmfd_ksp_monitor, mpi_err

    ! check for snes monitor
    if (cmfd_snes_monitor) call PetscOptionsSetValue("-snes_monitor", &
         "stdout", mpi_err)

    ! check for ksp monitor
    if (cmfd_ksp_monitor) call PetscOptionsSetValue("-ksp_monitor", &
         "stdout", mpi_err)

    end subroutine process_cmfd_options

!===============================================================================
! CALC_FISSION_SOURCE calculates the cmfd fission source
!===============================================================================

  subroutine calc_fission_source()

    use constants,  only: CMFD_NOACCEL, ZERO, TWO
    use global,     only: cmfd, cmfd_coremap, master, mpi_err, entropy_on, &
                                  fcpi_active

    integer :: nx ! maximum number of cells in x direction
    integer :: ny ! maximum number of cells in y direction
    integer :: nz ! maximum number of cells in z direction
    integer :: ng ! maximum number of energy groups
    integer :: n  ! total size
    integer :: i ! iteration counter for x
    integer :: j ! iteration counter for y
    integer :: k ! iteration counter for z
    integer :: g ! iteration counter for groups
    integer :: idx ! index in vector
    real(8) :: hxyz(3) ! cell dimensions of current ijk cell
    real(8) :: vol     ! volume of cell
    real(8),allocatable :: source(:,:,:,:)  ! tmp source array for entropy

    ! get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)
    n  = ng*nx*ny*nz

    ! allocate cmfd source if not already allocated and allocate buffer
    if (.not. allocated(cmfd%cmfd_src)) allocate(cmfd%cmfd_src(ng,nx,ny,nz))
    ! allocate cmfd scattering src if fcpi run and not allocated
    if(fcpi_active) then
      if (.not. allocated(cmfd%cmfd_scatsrc)) allocate(cmfd%cmfd_scatsrc(ng,nx,ny,nz))
    end if

    ! reset cmfd source to 0
    cmfd%cmfd_src = ZERO

    ! only perform for master
    if (master) then

      ! loop around indices to map to cmfd object
      ZLOOP: do k = 1, nz

        YLOOP: do j = 1, ny

          XLOOP: do i = 1, nx

            GROUP: do g = 1, ng

              ! check for core map
              if (cmfd_coremap) then
                if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) then
                  cycle
                end if
              end if

              ! get dimensions of cell
              hxyz = cmfd%hxyz(:,i,j,k)

              ! calculate volume
              vol = hxyz(1)*hxyz(2)*hxyz(3)

              ! get first index
              idx = get_matrix_idx(1,i,j,k,ng,nx,ny)

              ! compute fission source
              ! UPDATE K.Keady 10/12 -- we're flipping this around :D
              !cmfd%cmfd_src(g,i,j,k) = sum(cmfd%nfissxs(:,g,i,j,k) * &
              !     cmfd%phi(idx:idx+(ng-1)))*vol
              cmfd%cmfd_src(ng-g+1,i,j,k) = sum(cmfd%nfissxs(:,g,i,j,k) * &
                   cmfd%phi(idx:idx+(ng-1)))*vol
              if (fcpi_active) then
                cmfd %cmfd_scatsrc(ng-g+1,i,j,k) = sum(cmfd%scattxs(:,g,i,j,k) * &
                    cmfd%phi(idx:idx+(ng-1)))*vol
              end if
            end do GROUP

          end do XLOOP

        end do YLOOP

      end do ZLOOP

      ! normalize source such that it sums to 1.0
      cmfd%cmfd_src = cmfd%cmfd_src/sum(cmfd%cmfd_src)
      ! same with scattering source
      if(fcpi_active) then
        cmfd%cmfd_scatsrc = cmfd%cmfd_scatsrc/sum(cmfd%cmfd_scatsrc)      
      end if
      
      ! compute entropy
      if (entropy_on) then

        ! allocate tmp array
        if (.not.allocated(source)) allocate(source(ng,nx,ny,nz))

        ! initialize the source
        source = ZERO

        ! compute log
        where (cmfd%cmfd_src > ZERO)
          source = cmfd%cmfd_src*log(cmfd%cmfd_src)/log(TWO)
        end where

        ! sum that source
        cmfd%entropy = -sum(source)

        ! deallocate tmp array
        if (allocated(source)) deallocate(source)

      end if

      ! normalize source so average is 1.0
      cmfd%cmfd_src = cmfd%cmfd_src/sum(cmfd%cmfd_src)*cmfd%norm
      if(fcpi_active) then
        cmfd%cmfd_scatsrc = cmfd%cmfd_scatsrc/sum(cmfd%cmfd_scatsrc)*cmfd%norm
      end if
    end if

    ! broadcast full source to all procs
    call MPI_BCAST(cmfd%cmfd_src, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)

  end subroutine calc_fission_source

!===============================================================================
! CMFD_REWEIGHT
!===============================================================================

  subroutine cmfd_reweight(new_weights)

    use constants,   only: ZERO, ONE, FISS_P, SCAT_P
    use error,       only: warning, fatal_error
    use global,      only: n_particles, meshes, source_bank, work,             &
                           n_user_meshes, message, cmfd, master, mpi_err,       &
                           fcpi_active, current_batch
    use mesh_header, only: StructuredMesh
    use mesh,        only: count_bank_sites, get_mesh_indices   
!#ifndef MULTIGROUP  
    use search,      only: binary_search
!#endif

    ! local variables
    integer :: j
    integer :: k
    integer :: s ! iterator to transpose cmfd src matrix
    integer :: nx ! maximum number of cells in x direction
    integer :: ny ! maximum number of cells in y direction
    integer :: nz ! maximum number of cells in z direction
    integer :: ng ! maximum number of energy groups
    integer :: i ! iteration counter
    integer :: ijk(3) ! spatial bin location
    integer :: e_bin ! energy bin of source particle
    integer :: n_grps ! number of energy groups
    logical :: in_mesh ! source site is inside mesh
    logical :: new_weights ! calcualte new weights
    logical :: outside ! any source sites outside mesh
    type(StructuredMesh), pointer :: m ! point to mesh
    real(8), allocatable :: egrid(:)
    real(8), allocatable :: temp_cmfd_src(:,:,:,:)    
    
    ! associate pointer
    m => meshes(n_user_meshes + 1)

    ! get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)


    ! allocate arrays in cmfd object (can take out later extend to multigroup)
    if (.not.allocated(cmfd%sourcecounts)) then 
      allocate(cmfd%sourcecounts(ng,nx,ny,nz))
      cmfd % sourcecounts = 0
    end if
    if (.not.allocated(cmfd%weightfactors)) then 
      allocate(cmfd%weightfactors(ng,nx,ny,nz))
      cmfd % weightfactors = ONE
    end if
    
    if(fcpi_active) then
      if(.not.allocated(cmfd%scat_sourcecounts)) then
        allocate(cmfd%scat_sourcecounts(ng,nx,ny,nz))
        cmfd % scat_sourcecounts = 0
      end if
      if (.not.allocated(cmfd%scat_weightfactors)) then 
        allocate(cmfd%scat_weightfactors(ng,nx,ny,nz))
        cmfd % scat_weightfactors = ONE
      end if
    end if
    
    
    ! allocate energy grid and reverse cmfd energy grid
    if (.not. allocated(egrid)) allocate(egrid(ng+1))
    egrid = (/(cmfd%egrid(ng-i+2),i = 1,ng+1)/)
    
    ! compute new weight factors
    if (new_weights) then

      ! zero out weights
      cmfd%weightfactors = ZERO

      ! count bank sites in mesh
      call count_bank_sites(m, source_bank, cmfd%sourcecounts, egrid, &
           sites_outside=outside, part_type=FISS_P, size_bank=work)

      ! if fcpi is active, count scattering source particles too
      if(fcpi_active) then
           call count_bank_sites(m, source_bank, cmfd%scat_sourcecounts, egrid, &
           sites_outside=outside, part_type=SCAT_P, size_bank=work)            
      end if
           
     !print *,"mc sites in fiss src: ", cmfd%sourcecounts
     !print *,"cmfd fiss src shape: ", cmfd%cmfd_src    
     !print *,"mc sites in scat src: ", cmfd%scat_sourcecounts
     !print *,"cmfd scat src shape: ",cmfd%cmfd_scatsrc

      ! check for sites outside of the mesh
      if (master .and. outside) then
        message = "Source sites outside of the CMFD mesh!"
        call fatal_error()
      end if

      ! have master compute weight factors
      if (master) then
        where(cmfd%cmfd_src > ZERO .and. cmfd%sourcecounts > ZERO)
          cmfd%weightfactors = cmfd%cmfd_src/sum(cmfd%cmfd_src)* &
                               sum(cmfd%sourcecounts) / cmfd%sourcecounts
        end where
        ! do this for scattering source too, if fcpi is active
        if  (fcpi_active) then
          where(cmfd%cmfd_scatsrc > ZERO .and. cmfd%scat_sourcecounts > ZERO)
            cmfd%scat_weightfactors = cmfd%cmfd_scatsrc/sum(cmfd%cmfd_scatsrc)* &
                                 sum(cmfd%scat_sourcecounts) / cmfd%scat_sourcecounts 
         end where
        end if
      end if

      !print *,"WEIGHT FACTORS:    ",cmfd%weightfactors 

      ! broadcast weight factors to all procs
      call MPI_BCAST(cmfd%weightfactors, ng*nx*ny*nz, MPI_REAL8, 0, &
           MPI_COMM_WORLD, mpi_err)

   end if

    ! begin loop over source bank
    do i = 1, int(work, 4) 

      ! determine spatial bin
      call get_mesh_indices(m, source_bank(i)%xyz, ijk, in_mesh)

      ! determine energy bin
      n_grps = size(cmfd%egrid) - 1
      if (source_bank(i) % E < cmfd%egrid(1)) then
        e_bin = 1
        message = 'source pt below energy grid'
        call warning()
      elseif (source_bank(i) % E > cmfd%egrid(n_grps+1)) then
        e_bin = n_grps
        message = 'source pt above energy grid'
        call warning()
      else
        e_bin = binary_search(cmfd%egrid, n_grps + 1, source_bank(i) % E)

      end if

#ifndef MULTIGROUP
      ! reverse order (lowest energy GROUP is highest ENERGY,
      ! for continuous energy MC case)
      e_bin = n_grps - e_bin + 1
#endif

      
      ! check for outside of mesh
      if (.not. in_mesh) then
        message = 'Source site found outside of CMFD mesh!'
        call fatal_error()
      end if

      ! reweight particle
      if(source_bank(i)%type == FISS_P) then
        source_bank(i)%wgt = source_bank(i)%wgt * &
             cmfd%weightfactors(e_bin,ijk(1),ijk(2),ijk(3))
      else
        ! it's a scatter particle!
        ! this should ONLY happen if fcpi is active-- add a check!!! 
        if(.not. fcpi_active) then
          message = "WHOA NO NO NO"
          call fatal_error()
        end if
        source_bank(i)%wgt = source_bank(i)%wgt * &
             cmfd%scat_weightfactors(e_bin,ijk(1),ijk(2),ijk(3))        
           !print *,source_bank(i)%wgt ,cmfd%weightfactors(e_bin,ijk(1),ijk(2),ijk(3))
             end if
    end do

    ! deallocate
    if (allocated(egrid)) deallocate(egrid)

  end subroutine cmfd_reweight

!===============================================================================
! CALC_SPEC_RAD estimates the spectral radius of the iteration scheme
!===============================================================================
  subroutine calc_spec_rad()
    
    use global, only: cmfd, overall_gen
    
    real(8) :: spec_rad  ! The spectral radius estimate
    real(8) :: norm_num
    real(8) :: norm_denom
    
    norm_num = 0.0
    norm_denom = 0.0
    
    if(overall_gen > 3) then
      norm_num = sum((cmfd % flux - cmfd % flux_old)*(cmfd % flux - cmfd % flux_old))
      norm_denom = sum((cmfd % flux_old - cmfd % flux_old_old)*(cmfd % flux_old - cmfd % flux_old_old))
      
      spec_rad = sqrt(norm_num)/sqrt(norm_denom)
      print *,"Spectral radius estimate: ", spec_rad
    end if
    
  end subroutine calc_spec_rad
  
!===============================================================================
! GET_MATRIX_IDX takes (x,y,z,g) indices and computes location in matrix
!===============================================================================

  function get_matrix_idx(g, i, j, k, ng, nx, ny) result (matidx)

    use global, only: cmfd, cmfd_coremap

    integer :: matidx ! the index location in matrix
    integer :: i      ! current x index
    integer :: j      ! current y index
    integer :: k      ! current z index
    integer :: g      ! current group index
    integer :: nx     ! maximum number of cells in x direction
    integer :: ny     ! maximum number of cells in y direction
    integer :: ng     ! maximum number of energy groups

    ! check if coremap is used
    if (cmfd_coremap) then

      ! get idx from core map
      matidx = ng*(cmfd % coremap(i,j,k)) - (ng - g)

    else

      ! compute index
      matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    end if

  end function get_matrix_idx

# endif

!===============================================================================
! CMFD_TALLY_RESET
!===============================================================================

  subroutine cmfd_tally_reset()

    use global,  only: n_cmfd_tallies, cmfd_tallies, message, &
                       cmfd_cfl
    use output,  only: write_message
    use tally,   only: reset_result

    integer :: i ! loop counter

    ! print message
    !message = "CMFD tallies reset"
    !call write_message(7)
    cmfd_cfl = 'T'

    ! begin loop around CMFD tallies
    do i = 1, n_cmfd_tallies

      ! reset that tally
      cmfd_tallies(i) % n_realizations = 0
      call reset_result(cmfd_tallies(i) % results)

    end do

  end subroutine cmfd_tally_reset

end module cmfd_execute
