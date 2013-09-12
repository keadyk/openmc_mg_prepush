module cmfd_output

#ifdef PETSC

! This modules cleans up cmfd objects and echos the results

  implicit none
  private
  public :: finalize_cmfd

contains

!===============================================================================
! FINALIZE_CMFD
!===============================================================================

  subroutine finalize_cmfd() 

    use global,      only: cmfd, cmfd_write_balance, cmfd_write_hdf5, &
                           master, mpi_err
    use cmfd_header, only: deallocate_cmfd

    ! finalize petsc
    call PetscFinalize(mpi_err)

    ! write out final neutron balance
    if (master .and. cmfd_write_balance) call write_neutron_balance()
    
    ! write out cmfd eigenfunction
    if (master) call write_cmfd_eigenfunction(cmfd)

    ! deallocate cmfd object
    call deallocate_cmfd(cmfd)

  end subroutine finalize_cmfd 

!===============================================================================
! WRITE_NEUTRON_BALANCE 
!===============================================================================

  subroutine write_neutron_balance()

    use cmfd_data,    only: neutron_balance
    use constants,    only: CMFD_BALANCE, MAX_FILE_LEN
    use global,       only: path_output

    character(MAX_FILE_LEN) :: filename

    filename = trim(path_output) // 'neutron_balance.out'

    ! open file for output
    open(UNIT=CMFD_BALANCE, FILE=filename, ACTION='write')

    ! write out the tally
    call neutron_balance(CMFD_BALANCE) 

    ! close file
    close(CMFD_BALANCE)

  end subroutine write_neutron_balance
  
!===============================================================================
! WRITE_CMFD_EIGENFUNCTION writes out the the CMFD eigenfunction to a file 
! with name '#_CMFD.out' or '#_HCMFD.out', where # is the initial random number
!  seed used for the MC simulation 
!===============================================================================

  subroutine write_cmfd_eigenfunction(this)

    use cmfd_header,  only: cmfd_type
    use constants,    only: CMFD_EIG, CMFD_NOACCEL, MAX_FILE_LEN
    use global,       only: path_output, seed, cmfd, use_functs, cmfd_coremap
    use string,       only: to_str
    
    integer :: i        ! x-index
    integer :: j        ! y-index
    integer :: k        ! z-index 
    integer :: g        ! group index
    integer :: nx       ! number of mesh cells in x
    integer :: ny       ! number of mesh cells in y
    integer :: nz       ! number of mesh cells in z
    integer :: ng       ! number of groups
    integer :: matidx   ! index of eigenfunction entry
    character(MAX_FILE_LEN) :: filename   ! output filename
    
    type(cmfd_type) :: this   ! cmfd object

    ! extract spatial and energy indices from object
    nx = this % indices(1)
    ny = this % indices(2)
    nz = this % indices(3)
    ng = this % indices(4)

    if(use_functs) then
      filename = trim(path_output) // trim(to_str(seed)) // '_HCMFD.out'
    else
      filename = trim(path_output) // trim(to_str(seed)) // '_CMFD.out'
    end if
    ! open file for output
    open(UNIT=CMFD_EIG, FILE=filename, ACTION='write')

    ! write out the tally
    write(CMFD_EIG,'(A)') "#Location (i,j,k), Group, Eigenfunction"    
    
    ! begin loop around space and energy groups
    ZLOOP: do k = 1, nz

      XLOOP: do i = 1, nx

        YLOOP: do j = 1, ny

          GROUPG: do g = 1, ng
              ! check if coremap is used
            if (cmfd_coremap) then
              ! get idx from core map
              matidx = ng*(this % coremap(i,j,k)) - (ng - g)
              
              if (this % coremap(i,j,k) == CMFD_NOACCEL) then
                ! (chose a non-zero value so the rsd calc wouldn't get f-ed up)
                this % phi_final(matidx) = 1.00
              end if
              
              ! write out info to file
              write(CMFD_EIG,'(1X,I0,1X,I0,1X,I0,1X,I0,1X,E20.7)') &
                          i,j,k,g,this % phi_final(matidx)
            else
              ! compute the location of this entry
              matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)
              
              ! write out info to file
              write(CMFD_EIG,'(1X,I0,1X,I0,1X,I0,1X,I0,1X,E20.7)') &
                          i,j,k,g,this % phi_final(matidx)
            end if
            
          end do GROUPG
          
        end do YLOOP
        
      end do XLOOP
      
    end do ZLOOP

    ! close file
    close(CMFD_EIG)

  end subroutine write_cmfd_eigenfunction


#endif

end module cmfd_output
