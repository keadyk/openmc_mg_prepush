module roi

  use constants
  use error,             only: fatal_error, warning
  use global
  use output,            only: write_message
  use particle_header,   only: LocalCoord
  use random_lcg,        only: prn
  use string,            only: to_str
  use tally,             only: score_surface_current

  implicit none

contains

!===============================================================================
! CHECK_CELL is called when the particle crosses into a new cell, and the roi
! method is enabled. It determines if the particle has crossed into or out of 
! a buffer or roi region, and calls split/roulette methods as needed.
!===============================================================================

  subroutine check_cell(last_cell)
    
    integer, intent(in)       :: last_cell  ! last cell particle was in

    real(8)                   :: ratio           ! density ratio
    type(Cell),       pointer :: c_new           ! pointer to new cell
    type(Cell),       pointer :: c_new0           ! pointer to new cell
    type(Cell),       pointer :: c_old           ! pointer to old
    ! type(LocalCoord), pointer :: coord => null() 

    c_new => cells(p % coord % cell)             ! set to point to current cell
    c_new0 => cells(p % coord0 % cell)             ! set to point to current cell
    c_old => cells(last_cell)                    ! set to point to last cell

    
    
    ratio = real(c_new % n_split, 8)/real(c_old % n_split, 8)
    
    if(ratio /= ONE) then
      message = "NEW: " // trim(to_str(c_new % id)) // ", OLD: " // trim(to_str(c_old % id)) 
      call write_message(5)
      message = "Ratio is: " // trim(to_str(ratio))
      call write_message(5)
    end if
    
    if (ratio < ONE) then
      ! Need to first bring current tallies 'up to date' before wt chg
      if (active_current_tallies % size() > 0) call score_surface_current()
      ! Update coordinates for tallying purposes
      p % last_xyz = p % coord0 % xyz
      call roulette(ratio)
    elseif (ratio > ONE) then
      ! Need to first bring current tallies 'up to date' before wt chg
      if (active_current_tallies % size() > 0) call score_surface_current()
      ! Update coordinates for tallying purposes
      p % last_xyz = p % coord0 % xyz

      call split(ratio)
    else
      !message = "Doing nothing!"
      !call write_message(5)
      ! Do nothing!
    end if
    
  end subroutine check_cell
  
!===============================================================================
! CHECK_SRC is called when the roi method is enabled, and a fission particle is
! sampled. If the fission particle is born in a region where higher particle 
! density is desired, the fission particle is automatically split 
!===============================================================================

  subroutine check_src(birth_cell)
    
    integer, intent(in)       :: birth_cell  ! cell particle was born into

    real(8)                   :: split_factor    ! density ratio
    type(Cell),       pointer :: c               ! pointer to new cell

    c => cells(birth_cell)                    ! set to point to birth cell
    
    split_factor = real(c % n_split)

    ! message = "Ratio is: " // trim(to_str(ratio))
    !call write_message(5)
    
    if (split_factor == ONE) then
      ! Do nothing!
      return
    elseif (split_factor > ONE) then
      call split(split_factor)
    end if
    
  end subroutine check_src
  
!===============================================================================
! ROULETTE selectively kills particles moving from a cell with higher preset 
! particle density to a cell with lower preset particle density. If the particle
! survives, its weight is adjusted to maintain a fair simulation
!===============================================================================
  subroutine roulette(ratio)

    real(8), intent(in) :: ratio           ! splitting ratio
    real(8)             :: rand            ! random number
    
    rand = prn();
    
    !message = "weight " // trim(to_str(p % wgt)) // " ratio/rand " // trim(to_str(ratio)) // " " // trim(to_str(rand))
    !call write_message(5)
        
    if (rand < ratio) then                ! bump up the particle's weight
      p % wgt = p % wgt / ratio
      !if (current_batch == 10 .or. current_batch == 11) then
      !if(current_batch == 10) then
      !message = "Survived w wgt " // trim(to_str(p % wgt))
      !call write_message(5)
      !end if
      return
    else                                   ! kill the particle
      p % alive = .false.
      p % wgt = 0.0
      !if (current_batch == 10 .or. current_batch == 11) then
        !if(current_batch == 10) then
      !message = "Killed!" 
      !call write_message(5)
      !end if
      return
    end if
  end subroutine roulette
  
!===============================================================================
! SPLIT breaks a particle into multiple daughters, with the exact number 
! dependent on the desired particle density ratio between the last and current 
! cells. After a split, the original particle keeps running (with appropriately 
! reduced weight)
!===============================================================================
  subroutine split(ratio)
    
    real(8), intent(in) :: ratio           ! splitting ratio

    integer    :: i                        ! loop index
    integer(4) :: tot_split                ! total number of splits
    real(8)    :: split_w                  ! weight of each split 
    
    ! Determine total number of split particles
    tot_split = int(ratio, 4)
    
    ! Determine weight of each split particle
    split_w = p % wgt / ratio
    
    ! Skip if split bank is full already
    if (n_sbank == n_split*n_split*(work/10)) then
      message = "Split bank full!"
      call write_message(5)
      return
    end if
    
    !if (current_batch == 10 .or. current_batch == 11) then
    !message = "Split into " // trim(to_str(tot_split)) // &
    !     " daughters of weight " // trim(to_str(split_w))
    !     call write_message(5)
    !end if     
    
    !message = "Split bank, max size: " // trim(to_str(n_sbank)) // " " // trim(to_str(3*work))
    !call write_message(5)
    
    ! Bank split neutrons, less the one that will continue running
    do i = int(n_sbank,4) + 1, int(min(n_sbank + (tot_split-1), n_split*n_split*(work/10)),4)

      ! Bank split neutrons by copying particle data
      split_bank(i) % xyz = p % coord0 % xyz
      
      ! Set direction of split bank site
      split_bank(i) % uvw = p % coord0 % uvw
      
      ! Set weight of split bank site
      split_bank(i) % wgt = split_w
      
      ! Set energy/energy group of split bank site
      split_bank(i) % E = p % E
      
    end do
    
    ! Decrease weight of currently running particle accordingly
    p % wgt = split_w
    
    ! Add new splits to bank (less the one that's still running)
    n_sbank = n_sbank + (tot_split - 1)
 
  end subroutine split
  
  subroutine adjust_weight_cutoffs()

    integer :: i  ! Loop index
    
    do i = 1, n_cells
      if(cells(i) % n_split > ONE) then
        weight_cutoffs(i) = weight_cutoffs(i) / real(cells(i) % n_split)
        weight_survives(i) = weight_survives(i) / real(cells(i) % n_split) 
      end if
    end do
    
  end subroutine adjust_weight_cutoffs
  
end module roi