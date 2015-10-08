module bank_header

  implicit none

!===============================================================================
! BANK is used for storing fission sites in eigenvalue calculations. Since all
! the state information of a neutron is not needed, this type allows sites to be
! stored with less memory
!===============================================================================

  type Bank
    ! The 'sequence' attribute is used here to ensure that the data listed
    ! appears in the given order. This is important for MPI purposes when bank
    ! sites are sent from one processor to another.
    sequence

    real(8)    :: wgt    ! weight of bank site
    real(8)    :: xyz(3) ! location of bank particle
    real(8)    :: uvw(3) ! directional cosines
#ifdef MULTIGROUP
    integer    :: E      ! energy group
#else
    real(8)    :: E      ! energy
#endif
  end type Bank

!===============================================================================
! FCPI_BANK is used for storing intermediate fission sites in fcpi eigenvalue calculations. 
! It's essentially just the usual bank, but with an additional variable for the # of collisions
!===============================================================================
  ! Added by K. Keady on 10/8/15. Was going to just
  ! extend Bank type, but apparently you can't extend
  ! a sequence type...
  type FCPI_Bank
    ! The 'sequence' attribute is used here to ensure that the data listed
    ! appears in the given order. This is important for MPI purposes when bank
    ! sites are sent from one processor to another.
    sequence

    real(8)    :: wgt    ! weight of bank site
    real(8)    :: xyz(3) ! location of bank particle
    real(8)    :: uvw(3) ! directional cosines
#ifdef MULTIGROUP
    integer    :: E      ! energy group
#else
    real(8)    :: E      ! energy
#endif
    integer :: n_collision
  end type  FCPI_Bank
  
end module bank_header
