module multigroup_physics

  use multigroup_header, only: Nuclide, Reaction
  use constants
  use multigroup_cross_section,     only: calculate_xs
  use endf,              only: reaction_name
  use error,             only: fatal_error, warning
  use fission,           only: nu_total
  use geometry,          only: find_cell, distance_to_boundary, cross_surface, &
                               cross_lattice
  use geometry_header,   only: Universe, BASE_UNIVERSE
  use global
  use material_header,   only: Material
  use mesh,              only: get_mesh_indices
  use output,            only: write_message
  use particle_header,   only: LocalCoord
  use random_lcg,        only: prn
  use roi,               only: check_cell, check_src
  use search,            only: binary_search
  use string,            only: to_str
  use tally,             only: score_analog_tally, score_tracklength_tally, &
                               score_surface_current

  implicit none
  
contains

!===============================================================================
! TRANSPORT encompasses the main logic for moving a particle through geometry.
!===============================================================================

  subroutine transport()

    integer :: surface_crossed ! surface which particle is on
    integer :: lattice_crossed ! lattice boundary which particle crossed
    integer :: n_event         ! number of collisions/crossings
    real(8) :: d_boundary      ! distance to nearest boundary
    real(8) :: d_collision     ! sampled distance to collision
    real(8) :: distance        ! distance particle travels
    logical :: found_cell      ! found cell which particle is in?
    type(LocalCoord), pointer :: coord => null()

    ! Display message if high verbosity or trace is on
    if (verbosity >= 9 .or. trace) then
      message = "Simulating Particle " // trim(to_str(p % id))
      call write_message()
    end if

    ! If the cell hasn't been determined based on the particle's location,
    ! initiate a search for the current cell
    if (p % coord % cell == NONE) then
      call find_cell(found_cell)

      ! Particle couldn't be located
      if (.not. found_cell) then
        message = "Could not locate particle " // trim(to_str(p % id))
        call fatal_error()
      end if

      ! set birth cell attribute
      p % cell_born = p % coord % cell
    end if
    
    ! Increment track dist in birth cell
    if(roi_on .and. active_batches) then
      track_dist(p % coord % cell) = track_dist(p % coord % cell) + 1
    end if
    
    ! Initialize number of events to zero
    n_event = 0

    ! Add particle's starting weight to count for normalizing tallies later
    ! ONLY IF IT IS NOT A SPLIT OR FCPI PARTICLE! The "parent" of a split particle
    ! has already been counted, so we shouldn't count it again
    if (.not. split_particle .and. .not. fcpi_particle) then
      total_weight = total_weight + p % wgt
      ! In the roi method, we also have to 'split' fission particles born
      ! in the active region...
      if(active_batches .and. roi_on) then
        call check_src(p % cell_born)
      end if
    end if
      
    ! Force calculation of cross-sections by setting last energy to zero 
    micro_xs % last_E = 0

    do while (p % alive)
      ! Calculate microscopic and macroscopic cross sections -- note: if the
      ! material is the same as the last material and the energy of the
      ! particle hasn't changed, we don't need to lookup cross sections again.
      
      if (p % material /= p % last_material) call calculate_xs()
      
      ! Find the distance to the nearest boundary
      call distance_to_boundary(d_boundary, surface_crossed, lattice_crossed)
    
      ! Sample a distance to collision
      if (material_xs % total == ZERO) then
        d_collision = INFINITY
      else
        d_collision = -log(prn()) / material_xs % total
      end if

      ! Select smaller of the two distances
      distance = min(d_boundary, d_collision)
      
      ! Advance particle
      coord => p % coord0
      do while (associated(coord))
        coord % xyz = coord % xyz + distance * coord % uvw
        ! Update last cell id so it's guaranteed to be lowest-level
        last_cell = coord % cell
        coord => coord % next
      end do
      
      ! Score track-length tallies
      if (active_tracklength_tallies % size() > 0) &
           call score_tracklength_tally(distance)

!      write(*, '(A10,E14.7)') "nusigf:   ", material_xs % nu_fission     
      ! Score track-length estimate of k-eff
      global_tallies(K_TRACKLENGTH) % value = &
           global_tallies(K_TRACKLENGTH) % value + p % wgt * distance * &
           material_xs % nu_fission

      if (d_collision > d_boundary) then
        ! ====================================================================
        ! PARTICLE CROSSES SURFACE
        p % coord % cell = NONE
        if (lattice_crossed /= NONE) then
          ! Particle crosses lattice boundary
          p % surface = NONE
                    p % event = EVENT_LATTICE
          call cross_lattice(lattice_crossed)

        else
          ! Particle crosses surface
          p % surface = surface_crossed
                    p % event = EVENT_SURFACE 
          call cross_surface(last_cell)
       
        end if
        
        if (roi_on .and. active_batches) then
          if (p % alive) call check_cell(last_cell)
          ! If particle survives, tally its track in the new cell 
          if (p % alive) track_dist(p % coord % cell) = track_dist(p % coord % cell) + 1
        end if   
      else     
        ! ====================================================================
        ! PARTICLE HAS COLLISION
        ! Score collision estimate of keff
        global_tallies(K_COLLISION) % value = &
             global_tallies(K_COLLISION) % value + p % wgt * &
             material_xs % nu_fission / material_xs % total

        p % surface = NONE
        p % event = NONE
        
        call collision()

        ! Save coordinates for tallying purposes
        p % last_xyz = p % coord0 % xyz

        ! Set last material to none since cross sections will need to be
        ! re-evaluated
        p % last_material = NONE

        ! Set all uvws to base level -- right now, after a collision, only the
        ! base level uvws are changed
        coord => p % coord0
        do while(associated(coord % next))
          if (coord % next % rotated) then
            ! If next level is rotated, apply rotation matrix
            coord % next % uvw = matmul(cells(coord % cell) % &
                 rotation, coord % uvw)
          else
            ! Otherwise, copy this level's direction
            coord % next % uvw = coord % uvw
          end if

          ! Advance coordinate level
          coord => coord % next
        end do
      end if

      ! If particle has too many events, display warning and kill it
      n_event = n_event + 1
      if (n_event == MAX_EVENTS) then
        message = "Particle " // trim(to_str(p%id)) // " underwent maximum &
             &number of events."
        call warning()
        p % alive = .false.
      end if

    end do

  end subroutine transport      
  
!===============================================================================
! COLLISION samples a nuclide and reaction and then calls the appropriate
! routine for that reaction
!===============================================================================

  subroutine collision()

    ! Store pre-collision particle properties
    p % last_wgt = p % wgt
    p % last_E   = p % E
    p % last_uvw = p % coord0 % uvw
    
    ! Add to collision counter for particle
    p % n_collision = p % n_collision + 1

    ! score surface current tallies -- this has to be done before the collision
    ! since the direction of the particle will change and we need to use the
    ! pre-collision direction to figure out what mesh surfaces were crossed

    if (active_current_tallies % size() > 0) then
       call score_surface_current()
    end if

    ! Sample nuclide/reaction for the material the particle is in
    call sample_reaction()

    ! Display information about collision
    if (verbosity >= 10 .or. trace) then
      message = "    " // trim(reaction_name(p % event_MT)) // " with " // &
           trim(adjustl(nuclides(p % event_nuclide) % name)) // &
           ". Group # = " // trim(to_str(p % E)) // "." 
      call write_message()
    end if

    ! Score collision estimator tallies -- this is done after a collision has
    ! occurred rather than before because we need information on the outgoing
    ! energy for any tallies with an outgoing energy filter

    if (active_analog_tallies % size() > 0) call score_analog_tally()

    ! Reset banked weight during collision
    p % n_bank   = 0
    p % wgt_bank = ZERO

  end subroutine collision
  
!===============================================================================
! SAMPLE_REACTION samples a nuclide based on the macroscopic cross sections for
! each nuclide within a material and calls the appropriate routine to process
!  the physics. Note that there is special logic when survival biasing is 
!  turned on since fission and disappearance are treated implicitly.
!===============================================================================

  subroutine sample_reaction()

    integer :: i            ! index over nuclides in a material
    integer :: i_nuclide    ! index in nuclides array
    integer :: i_group      ! index of energy group
    real(8) :: sigma        ! microscopic total xs for nuclide
    real(8) :: prob         ! cumulative probability
    real(8) :: cutoff       ! random number
    real(8) :: atom_density ! atom density of nuclide in atom/b-cm
    type(Material), pointer :: mat => null()
    type(Nuclide),  pointer :: nuc => null()
    type(Reaction), pointer :: rxn => null()

    ! Get pointer to current material
    mat => materials(p % material)

    ! ==========================================================================
    ! SAMPLE NUCLIDE WITHIN THE MATERIAL

    cutoff = prn() * material_xs % total
    prob = ZERO

    i = 0
    do while (prob < cutoff)
      i = i + 1

      ! Check to make sure that a nuclide was sampled
      if (i > mat % n_nuclides) then
        message = "Did not sample any nuclide during collision."
        call fatal_error()
      end if

      ! Find atom density and microscopic total cross section
      i_nuclide    = mat % nuclide(i)
      atom_density = mat % atom_density(i)
      sigma        = atom_density * micro_xs(i_nuclide) % total

      ! Increment probability to compare to cutoff
      prob = prob + sigma
    end do

    ! Get pointer to table, nuclide grid index
    nuc     => nuclides(i_nuclide)
    i_group =  micro_xs(i_nuclide) % index_grid

    ! Save which nuclide particle had collision with
    p % event_nuclide = i_nuclide

    ! ==========================================================================
    ! DISAPPEARANCE REACTIONS (ANALOG) OR IMPLICIT CAPTURE (SURVIVAL BIASING)

    if (survival_biasing) then
      ! Determine weight absorbed in survival biasing
      p % absorb_wgt = p % wgt * (micro_xs(i_nuclide) % absorption + micro_xs(i_nuclide) % fission)/ &
           micro_xs(i_nuclide) % total

      ! Adjust weight of particle by probability of absorption
      p % wgt = p % wgt - p % absorb_wgt
      p % last_wgt = p % wgt

      ! Score implicit absorption estimate of keff. Unlike the analog absorption
      ! estimate, this only needs to be scored to in one place.
      global_tallies(K_ABSORPTION) % value = &
           global_tallies(K_ABSORPTION) % value + p % absorb_wgt * &
           material_xs % nu_fission / (material_xs % absorption + material_xs % fission)

    else
      ! set cutoff variable for analog cases
      cutoff = prn() * micro_xs(i_nuclide) % total
      prob = ZERO

      ! Add disappearance cross-section to prob
      prob = prob + micro_xs(i_nuclide) % absorption 

          !message = "abs: " // trim(to_str(micro_xs(i_nuclide) % absorption)) // " and fission: " //  &
          ! trim(to_str(micro_xs(i_nuclide) % fission)) // "."
          ! call warning()

      ! See if disappearance reaction happens
      if (prob > cutoff) then
        ! Score absorption estimate of keff. Note that this appears in two
        ! places -- absorption reactions and fission reactions
        global_tallies(K_ABSORPTION) % value = &
             global_tallies(K_ABSORPTION) % value + p % wgt * &
             material_xs % nu_fission / (material_xs % absorption + material_xs % fission)

        global_tallies(K_GLOBAL_DENOM) % value = &
             global_tallies(K_GLOBAL_DENOM) % value + p % wgt            
             
             !print *,p%id, " is DEEEEAAAAD"
        p % alive = .false.
        p % event = EVENT_ABSORB
        p % event_MT = N_DISAPPEAR
        return
      end if
    end if

    ! ==========================================================================
    ! FISSION EVENTS (ANALOG) OR BANK EXPECTED FISSION SITES (IMPLICIT)

    if (nuc % fissionable) then
      ! If survival biasing is turned on, then no fission events actually occur
      ! since absorption is treated implicitly. However, we still need to bank
      ! sites so we sample a fission reaction (if there are multiple) and bank
      ! the expected number of fission neutrons created.
      if (survival_biasing) then
        cutoff = prn() * micro_xs(i_nuclide) % fission
        prob = ZERO
      end if

      FISSION_REACTION_LOOP: do i = 1, nuc % n_fission
        rxn => nuc % reactions(nuc % index_fission(i))
        
        ! add to cumulative probability (no partial fissions)
        prob = prob + micro_xs(i_nuclide) % fission

        ! Create fission bank sites if fission occurs
        if (prob > cutoff) then
          call create_fission_sites(i_nuclide, rxn)

          if (survival_biasing) then
            ! Since a fission reaction has been sampled, we can exit this
            ! loop
            exit FISSION_REACTION_LOOP
          else
            ! Score absorption estimate of keff. Note that this appears in
            ! two places -- absorption reactions and fission reactions
            global_tallies(K_ABSORPTION) % value = &
                 global_tallies(K_ABSORPTION) % value + p % wgt * &
                 material_xs % nu_fission / (material_xs % absorption + material_xs % fission)
            ! Score denom of global keff. Note that this appears in
            ! two places -- absorption reactions and fission reactions
            global_tallies(K_GLOBAL_DENOM) % value = &
                 global_tallies(K_GLOBAL_DENOM) % value + p % wgt 
            ! With no survival biasing, the particle is absorbed and so
            ! its life is over :( poor little fella
            p % alive = .false.
            p % event = EVENT_FISSION
            p % event_MT = rxn % MT
            return
          end if
        end if
      end do FISSION_REACTION_LOOP
    end if

    ! ==========================================================================
    ! WEIGHT CUTOFF (SURVIVAL BIASING ONLY)

    if (survival_biasing) then
      if (p % wgt < weight_cutoffs(last_cell)) then
        if (prn() < p % wgt / weight_survives(last_cell)) then
          p % wgt = weight_survives(last_cell)
          p % last_wgt = p % wgt
        else
          p % wgt = ZERO
          p % alive = .false.
          return
        end if
      end if

      ! At this point, we also need to set the cutoff variable for cases with
      ! survival biasing. The cutoff will be a random number times the
      ! scattering cross section

      cutoff = prn() * (micro_xs(i_nuclide) % total - &
           micro_xs(i_nuclide) % absorption)
      prob = ZERO
    end if

    ! ==========================================================================
    ! SCATTERING REACTION

    ! get pointer to scattering reaction
    rxn => nuc % reactions(1)

    ! Perform collision physics for multigroup scattering
    call multigroup_scatter(rxn)

    p % event_MT = N_LEVEL

    ! If we made it this far, it means that a scattering reaction took place
    ! since the absorption and fission blocks had return statements.
    p % event = EVENT_SCATTER

  end subroutine sample_reaction
 
!===============================================================================
! CREATE_FISSION_SITES determines the average total (prompt)
! neutrons produced from fission and creates appropriate bank sites.
!===============================================================================

  subroutine create_fission_sites(i_nuclide, rxn)

    integer, intent(in)     :: i_nuclide
    type(Reaction), pointer :: rxn

    integer :: i            ! loop index
    integer :: j            ! index on nu energy grid / precursor group
    integer :: lc           ! index before start of energies/nu values
    integer :: NR           ! number of interpolation regions
    integer :: NE           ! number of energies tabulated
    integer :: nu           ! actual number of neutrons produced
    integer :: ijk(3)       ! indices in ufs mesh
    integer :: E            ! incoming energy group of neutron
    integer :: E_out        ! outgoing energy group of fission neutron
    real(8) :: nu_t         ! total nu
    real(8) :: mu           ! fission neutron angular cosine
    real(8) :: phi          ! fission neutron azimuthal angle
    real(8) :: beta         ! delayed neutron fraction
    real(8) :: xi           ! random number
    real(8) :: prob         ! cumulative probability
    real(8) :: weight       ! weight adjustment for ufs method
    logical :: in_mesh      ! source site in ufs mesh?
    type(Nuclide),    pointer :: nuc

    ! Initialize mu
    mu = ZERO
    
    ! Get pointer to nuclide
    nuc => nuclides(i_nuclide)
    
    ! copy energy of neutron
    E = p % E

    ! Determine total nu
    nu_t = nu_total(nuc, E)
    
    ! If uniform fission source weighting is turned on, we increase or decrease
    ! the expected number of fission sites produced

    if (ufs) then
      ! Determine indices on ufs mesh for current location
      call get_mesh_indices(ufs_mesh, p % coord0 % xyz, ijk, in_mesh)
      if (.not. in_mesh) then
        message = "Source site outside UFS mesh!"
        call fatal_error()
      end if

      if (source_frac(1,ijk(1),ijk(2),ijk(3)) /= ZERO) then
        weight = ufs_mesh % volume_frac / source_frac(1,ijk(1),ijk(2),ijk(3))
      else
        weight = ONE
      end if
    else
      weight = ONE
    end if
    
    ! Sample number of neutrons produced
    if (survival_biasing) then
      ! Need to use the weight before survival biasing
      nu_t = (p % wgt + p % absorb_wgt) * micro_xs(i_nuclide) % fission / &
           (sample_keff * micro_xs(i_nuclide) % total) * nu_t * weight
    else 
      nu_t = p % wgt / sample_keff * nu_t * weight
      !print *,"sample k ", sample_keff, p % wgt, nu_t, weight
    end if
    if (prn() > nu_t - int(nu_t)) then
      nu = int(nu_t)
    else
      nu = int(nu_t) + 1
    end if

    ! tally global k numerator (neutron production)
    global_tallies(K_GLOBAL_NUM) % value = &
         global_tallies(K_GLOBAL_NUM) % value + nu * sample_keff     
    
    ! Edited by K.Keady on 10/10-- I TRIED TO MAKE THIS ELEGANT
    ! BUT NOTHING WORKS RIGHT WHEN EVERYTHING IS GLOBAL!!!!!
    ! Bank source neutrons
    ! If you're under max_coll and fcpi is activated, bank in intermed. bank
    if(fcpi_active .and. p % n_collision < max_coll) then
      if (nu == 0 .or. n_fbank == work) return
      
      !print *, p%id, " going into fiss bank w mult ", nu, " and coll ", p%n_collision
      do i = int(n_fbank,4) + 1, int(min(n_fbank + nu, work),4)
        if(nu < 0) then
          call fatal_error()
        end if
        if(nu > 100) then
          call fatal_error()
        end if
        ! Bank source neutrons by copying particle data
        int_fbank(i) % xyz = p % coord0 % xyz
        ! intermed. bank needs to know collision number
        int_fbank(i) % n_collision = p % n_collision
        ! Set weight of fission bank site
        int_fbank(i) % wgt = ONE/weight

        ! Sample cosine of angle 
        mu = TWO * prn() - ONE

        ! sample from prompt neutron chi distribution for the 
        ! fissioning nuclide
        E_out = sample_energy(nuc)

        ! Sample azimuthal angle uniformly in [0,2*pi)
        phi = TWO*PI*prn()
        int_fbank(i) % uvw(1) = mu
        int_fbank(i) % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
        int_fbank(i) % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

        ! set energy of fission neutron
        int_fbank(i) % E = E_out
      end do
      n_fbank = min(n_fbank + nu, 3*work)
      
    ! else, bank in final "end-of-cycle" bank
    else
      if (fcpi_active) then
        ! print *, "banking fiss with ncoll = ", p % n_collision
      end if
      if (nu == 0 .or. n_bank == 3*work) return
      do i = int(n_bank,4) + 1, int(min(n_bank + nu, 3*work),4)
        ! Bank source neutrons by copying particle data
        fission_bank(i) % xyz = p % coord0 % xyz
        ! Final bank needs to know particle type 
        fission_bank(i) % type = FISS_P
        ! Set weight of fission bank site
        fission_bank(i) % wgt = ONE/weight

        ! Sample cosine of angle -- fission neutrons are always emitted
        ! isotropically. Sometimes in ACE data, fission reactions actually have
        ! an angular distribution listed, but for those that do, it's simply just
        ! a uniform distribution in mu
        mu = TWO * prn() - ONE

        ! sample from prompt neutron chi distribution for the 
        ! fissioning nuclide
        E_out = sample_energy(nuc)

        ! Sample azimuthal angle uniformly in [0,2*pi)
        phi = TWO*PI*prn()
        fission_bank(i) % uvw(1) = mu
        fission_bank(i) % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
        fission_bank(i) % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

        ! set energy of fission neutron
        fission_bank(i) % E = E_out
      end do
      n_bank = min(n_bank + nu, 3*work)
      
    end if 

    ! Store total weight banked for analog fission tallies
    p % n_bank   = nu
    p % wgt_bank = nu/weight

  end subroutine create_fission_sites  
  
!===============================================================================
! MULTIGROUP_SCATTER treats the scattering of a neutron with a
! target using multigroup cross sections.
!===============================================================================

  subroutine multigroup_scatter(rxn)

    type(Reaction), pointer :: rxn

    real(8) :: mu      ! cosine of polar angle
    real(8) :: u       ! x-direction
    real(8) :: v       ! y-direction
    real(8) :: w       ! z-direction
    integer :: E       ! energy GROUP
    integer :: E_new   ! outgoing energy GROUP

    ! sample outgoing group
    E_new = sample_group(rxn, p % E)

    ! Sample scattering angle, given incoming, outgoing group
    mu = sample_angle(rxn, p % E, E_new)
!    write(*,'(E14.7)') mu
    
    ! copy directional cosines
    u = p % coord0 % uvw(1)
    v = p % coord0 % uvw(2)
    w = p % coord0 % uvw(3)

    ! change direction of particle
    call rotate_angle(u, v, w, mu)
    
    ! Set energy group and direction of particle
    p % E = E_new
    p % coord0 % uvw = (/ u, v, w /)

    ! Copy scattering cosine for tallies
    p % mu = mu
    
    if (fcpi_active .and. p % n_collision == max_coll) then
      ! Store it, then kill it
      call bank_fcpi_scatter()
      !print *, p%id, "banked"
    end if

  end subroutine multigroup_scatter

!===============================================================================
! BANK_FCPI_SCATTER banks a collided particle, then kills it in the current cycle
!===============================================================================
  subroutine bank_fcpi_scatter()
    integer :: i  ! index of this particle in bank
    ! YES, I KNOW IT'S CALLED THE FISSION BANK AND WE'RE PUTTING SCATTERS IN IT
    ! DEAL WITH IT
    ! print *, "banking scatter with ncoll = ", p % id, p % n_collision, (int(n_bank,4) + 1)
    !print *, p%id, " going into scatter bank"
    i = int(n_bank,4) + 1
    if(i <= 3*work) then
      ! Bank scattered neutron by copying particle data
      fission_bank(i) % xyz = p % coord0 % xyz
      ! Final bank needs to know particle type for cmfd sorting
      fission_bank(i) % type = SCAT_P
      fission_bank(i) % uvw = p % coord0 % uvw 
      fission_bank(i) % E = p % E 
      fission_bank(i) % wgt = p % wgt
    end if
    
    ! Update size of bank
    n_bank = min(n_bank + 1, 3*work)
    
    p % alive = .false.
    
  end subroutine bank_fcpi_scatter
  
  
!===============================================================================
! SAMPLE_GROUP samples the outgoing energy group, using the incident group and
! its associated scattering cross sections (stored in rxn).
!===============================================================================
  function sample_group(rxn, E) result(E_new)
  
  type(Reaction), pointer    :: rxn        ! sampled reaction (scattering)
  integer,        intent(in) :: E          ! incoming energy group
  integer         :: E_new      ! outgoing energy group
  real(8)         :: xi         ! random number
  real(8)         :: k          ! cumulative probability
  real(8)         :: cutoff     ! random number
  integer         :: offset     ! group search offset
  
  
  ! sample cutoff as fraction of total sigma_s for this energy group
  xi = prn()
  cutoff = xi * rxn % total_scatter(E)
  
  ! set offset and cumulative prob to zero
  k = 0.0
  offset = 0
  
  ! determine outgoing energy group
  E_new = rxn % max_scatter(E) - 1
  do while(k < cutoff) 
!    call warning()
    ! add this group's cross section
    k = k + rxn % sigma(rxn % group_index(E) + offset)
    offset = offset + 1
    E_new = E_new + 1
  end do
  
  end function sample_group
  
!===============================================================================
! SAMPLE_ANGLE samples the cosine of the angle between incident and exiting
! particle directions either from 32 equiprobable bins or from a tabular
! distribution.
!===============================================================================

  function sample_angle(rxn, E, E_new) result(mu)

    type(Reaction), pointer    :: rxn   ! reaction
    integer,        intent(in) :: E     ! incoming energy group
    integer,        intent(in) :: E_new ! outgoing energy group
    real(8)        :: xi      ! random number on [0,1)
    integer        :: index   ! index of distribution type for this g -> g' combo
    integer        :: type    ! angular distribution type
    integer        :: NT       ! number of total data points in ang distribution
    integer        :: lc      ! location in data array
    integer        :: NP      ! number of prob. data points in ang distribution
    integer        :: k       ! index on cosine grid
    real(8)        :: mu0     ! cosine in bin k
    real(8)        :: mu1     ! cosine in bin k+1
    real(8)        :: mu      ! final cosine sampled
    real(8)        :: c_k     ! cumulative frequency at k
    real(8)        :: c_k1    ! cumulative frequency at k+1
    real(8)        :: p0,p1   ! probability distribution

    ! check if reaction has angular distribution -- if not, sample outgoing
    ! angle isotropically
    if (.not. rxn % has_angle_dist) then
      mu = TWO * prn() - ONE
      return
    end if

    ! determine number of distribution data points
    NT = rxn % adist % n_data

    ! set index for this initial/final group combo (already determined)
    ! group_index(E) is the index for scattering OUT of group E into the
    ! highest possible energy group.
    ! (E_new - max_scatter(E)) gives the offset between the actual outgoing 
    ! energy group sampled and the highest possible energy group.
    index = rxn % group_index(E) + E_new - rxn % max_scatter(E)

!    write(*, '(I3,2X,I3,2X,I3)') E, E_new, index
    
    ! admittedly, this is kind of confusing... the index above gives the 
    ! location of the particular group-to-group scattering combination
    ! that we have sampled.  each group-to-group combo has it's own angular
    ! distribution, of length n_data.
    ! lc  = (index - 1) * rxn % adist % n_data + 1
    lc = (index - 1) * rxn % adist % n_data  ! Maybe this instead?
    
    type = rxn % adist % type(index)
      
    ! check whether this is an equiprobable bin or a tabular distribution
    if (type == ANGLE_ISOTROPIC) then
      mu = TWO * prn() - ONE
      
    elseif (type == ANGLE_NLEG_EQUI) then
      ! sample cosine bin
      xi = prn()
      k = 1 + int((NT-1) *xi)

!      message = "chosen bin: "
!      call warning()
!      write(*,'(I2)') k 
      
      ! calculate cosine
      mu0 = rxn % adist % data(lc + k)
      mu1 = rxn % adist % data(lc + k+1)
      
!      message = "bin bounds: "
!      call warning()
!      write(*,'(2E14.7)') mu0, mu1
      
!      FIX THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      mu = mu0 + ((NT-1) * xi - k + 1) * (mu1 - mu0)

    elseif (type == ANGLE_TABULAR) then
      ! NP = number of cum. probabilities listed
      NP = (NT - 1)/2

      ! determine outgoing cosine bin
      xi = prn()
!      message = "cum prob: "
!      call warning()
!      write(*,'(E14.7)') xi 
      c_k = rxn % adist % data(lc + 1)
!      message = "initial cum prob, assoc mu: "
!      call warning()
!      write(*, '(2E14.7)') c_k, rxn % adist % data(lc + NP + 1)
      do k = 1, NP
        c_k1 = rxn % adist % data(lc + k+1)
!        message = "initial cum prob, assoc mu: "
!        call warning()
!        write(*, '(2E14.7)') c_k1, rxn % adist % data(lc + k+1 + NP)
        if (xi < c_k1) exit
        c_k = c_k1
      end do

      ! check to make sure k is <= NP
      k = min(k, NP)

      p0  = rxn % adist % data(lc + NP + k)
      mu0 = rxn % adist % data(lc + k)
      ! Histogram interpolation
      if (p0 > ZERO) then
        mu = mu0 + (xi - c_k)/p0
      else
        mu = mu0
      end if

      ! Because of floating-point roundoff, it may be possible for mu to be
      ! outside of the range [-1,1). In these cases, we just set mu to exactly
      ! -1 or 1

      if (abs(mu) > ONE) mu = sign(ONE,mu)

    else
      message = "Unknown angular distribution type: " // trim(to_str(type))
      call fatal_error()
    end if

  end function sample_angle

!===============================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
! mu and through an azimuthal angle sampled uniformly. Note that this is done
! with direct sampling rather than rejection as is done in MCNP and SERPENT.
!===============================================================================

  subroutine rotate_angle(u, v, w, mu)

    real(8), intent(inout) :: u
    real(8), intent(inout) :: v
    real(8), intent(inout) :: w
    real(8), intent(in)    :: mu ! cosine of angle in lab or CM

    real(8) :: phi    ! azimuthal angle
    real(8) :: sinphi ! sine of azimuthal angle
    real(8) :: cosphi ! cosine of azimuthal angle
    real(8) :: a      ! sqrt(1 - mu^2)
    real(8) :: b      ! sqrt(1 - w^2)
    real(8) :: u0     ! original cosine in x direction
    real(8) :: v0     ! original cosine in y direction
    real(8) :: w0     ! original cosine in z direction

    ! Copy original directional cosines
    u0 = u
    v0 = v
    w0 = w

    ! Sample azimuthal angle in [0,2pi)
    phi = TWO * PI * prn()

    ! Precompute factors to save flops
    sinphi = sin(phi)
    cosphi = cos(phi)
    a = sqrt(max(ZERO, ONE - mu*mu))
    b = sqrt(max(ZERO, ONE - w0*w0))

    ! Need to treat special case where sqrt(1 - w**2) is close to zero by
    ! expanding about the v component rather than the w component
    if (b > 1e-10) then
      u = mu*u0 + a*(u0*w0*cosphi - v0*sinphi)/b
      v = mu*v0 + a*(v0*w0*cosphi + u0*sinphi)/b
      w = mu*w0 - a*b*cosphi
    else
      b = sqrt(ONE - v0*v0)
      u = mu*u0 + a*(u0*v0*cosphi + w0*sinphi)/b
      v = mu*v0 - a*b*cosphi
      w = mu*w0 + a*(v0*w0*cosphi - u0*sinphi)/b
    end if

  end subroutine rotate_angle

!===============================================================================
! SAMPLE_ENERGY samples the multi-group chi distribution 
! for a prompt fission neutron (dependent on the fissioned nuclide)
!===============================================================================
  
  function sample_energy(nuc) result(E)
        type(Nuclide), pointer   :: nuc  ! nuclide
        integer        ::   E        ! outgoing energy group
        real(8)        ::   sum      ! cumulative sum of chi values 
        real(8)        ::   cutoff   ! random number cutoff for chi
        
        cutoff = prn()  ! sample cutoff between 0 and 1
   
        sum = 0
        E = 0
        
        do while (sum < cutoff)
          E = E + 1
          sum = sum + nuc % chi_data(E)
        end do
        
  end function sample_energy
  
end module multigroup_physics