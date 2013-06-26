module multigroup_physics

  use multigroup_header, only: Nuclide, Reaction, DistEnergy
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
    integer :: last_cell       ! most recent cell particle was in
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

    ! Initialize number of events to zero
    n_event = 0

    ! Add paricle's starting weight to count for normalizing tallies later
    total_weight = total_weight + p % wgt

    ! Force calculation of cross-sections by setting last energy to zero 
    micro_xs % last_E = ZERO

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
        coord => coord % next
      end do

      ! Score track-length tallies
      if (active_tracklength_tallies % size() > 0) &
           call score_tracklength_tally(distance)

      ! Score track-length estimate of k-eff
      global_tallies(K_TRACKLENGTH) % value = &
           global_tallies(K_TRACKLENGTH) % value + p % wgt * distance * &
           material_xs % nu_fission

      if (d_collision > d_boundary) then
        ! ====================================================================
        ! PARTICLE CROSSES SURFACE

        last_cell = p % coord % cell
        p % coord % cell = NONE
        if (lattice_crossed /= NONE) then
          ! Particle crosses lattice boundary
          p % surface = NONE
          call cross_lattice(lattice_crossed)
          p % event = EVENT_LATTICE
        else
          ! Particle crosses surface
          p % surface = surface_crossed
          call cross_surface(last_cell)
          p % event = EVENT_SURFACE
        end if
      else
        ! ====================================================================
        ! PARTICLE HAS COLLISION

        ! Score collision estimate of keff
        global_tallies(K_COLLISION) % value = &
             global_tallies(K_COLLISION) % value + p % wgt * &
             material_xs % nu_fission / material_xs % total

        p % surface = NONE
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

    ! Add to collision counter for particle
    p % n_collision = p % n_collision + 1

    ! score surface current tallies -- this has to be done before the collision
    ! since the direction of the particle will change and we need to use the
    ! pre-collision direction to figure out what mesh surfaces were crossed

    if (active_current_tallies % size() > 0) call score_surface_current()

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
! each nuclide within a material and then samples a reaction for that nuclide
! and calls the appropriate routine to process the physics. Note that there is
! special logic when suvival biasing is turned on since fission and
! disappearance are treated implicitly.
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
      p % absorb_wgt = p % wgt * micro_xs(i_nuclide) % absorption / &
           micro_xs(i_nuclide) % total

      ! Adjust weight of particle by probability of absorption
      p % wgt = p % wgt - p % absorb_wgt
      p % last_wgt = p % wgt

      ! Score implicit absorption estimate of keff. Unlike the analog absorption
      ! estimate, this only needs to be scored to in one place.
      global_tallies(K_ABSORPTION) % value = &
           global_tallies(K_ABSORPTION) % value + p % absorb_wgt * &
           material_xs % nu_fission / material_xs % absorption

    else
      ! set cutoff variable for analog cases
      cutoff = prn() * micro_xs(i_nuclide) % total
      prob = ZERO

      ! Add disappearance cross-section to prob
      prob = prob + micro_xs(i_nuclide) % absorption - &
           micro_xs(i_nuclide) % fission

      ! See if disappearance reaction happens
      if (prob > cutoff) then
        ! Score absorption estimate of keff. Note that this appears in three
        ! places -- absorption reactions, total fission reactions, and
        ! first/second/etc chance fission reactions
        global_tallies(K_ABSORPTION) % value = &
             global_tallies(K_ABSORPTION) % value + p % wgt * &
             material_xs % nu_fission / material_xs % absorption

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

        ! Create fission bank sites if fission occus
        if (prob > cutoff) then
          call create_fission_sites(i_nuclide, rxn)

          if (survival_biasing) then
            ! Since a fission reaction has been sampled, we can exit this
            ! loop
            exit FISSION_REACTION_LOOP
          else
            ! Score absorption estimate of keff. Note that this appears in
            ! three places -- absorption reactions, total fission reactions,
            ! and first/second/etc chance fission reactions
            global_tallies(K_ABSORPTION) % value = &
                 global_tallies(K_ABSORPTION) % value + p % wgt * &
                 material_xs % nu_fission / material_xs % absorption

            ! With no survival biasing, the particle is absorbed and so
            ! its life is over :( poor little fella
            p % alive = .false.
            p % event = EVENT_FISSION
            p % event_MT = rxn % MT
            return
          end if
        end do FISSION_REACTION_LOOP

      end if
    end if

    ! ==========================================================================
    ! WEIGHT CUTOFF (SURVIVAL BIASING ONLY)

    if (survival_biasing) then
      if (p % wgt < weight_cutoff) then
        if (prn() < p % wgt / weight_survive) then
          p % wgt = weight_survive
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
    ! SCATTERING REACTIONS

    ! Shouldn't need these; if it ain't fission or absorp., it better be 
    ! scattering!!
    ! prob = prob + micro_xs(i_nuclide) % scattering
    ! if (prob < cutoff) then
      

    ! get pointer to scattering reaction
    rxn => nuc % reactions(1)

    ! Perform collision physics for multigroup scattering
    call multigroup_scatter(i_nuclide, rxn)

    p % event_MT = N_LEVEL

    ! If we made it this far, it means that a scattering reaction took place
    ! since the absorption and fission blocks had return statements.
    p % event = EVENT_SCATTER

  end subroutine sample_reaction
 
!===============================================================================
! MULTIGROUP_SCATTER treats the scattering of a neutron with a
! target using multigroup cross sections.
!===============================================================================

  subroutine multigroup_scatter(i_nuclide, rxn)

    integer, intent(in)     :: i_nuclide
    type(Reaction), pointer :: rxn

    real(8) :: awr     ! atomic weight ratio of target
    real(8) :: mu      ! cosine of polar angle
    real(8) :: vel     ! magnitude of velocity
    real(8) :: v_n(3)  ! velocity of neutron
    real(8) :: v_cm(3) ! velocity of center-of-mass
    real(8) :: v_t(3)  ! velocity of target nucleus
    real(8) :: u       ! x-direction
    real(8) :: v       ! y-direction
    real(8) :: w       ! z-direction
    real(8) :: E       ! energy GROUP
    real(8) :: E_new   ! outgoing energy GROUP
    type(Nuclide), pointer :: nuc => null()

    ! get pointer to nuclide
    nuc => nuclides(i_nuclide)

    vel = sqrt(p % E)
    awr = nuc % awr

    ! sample outgoing group
    E_new = sample_group(rxn, p % E)

    ! Sample scattering angle, given incoming, outgoing group
    mu = sample_angle(rxn, p % E, E_new)

    ! Change direction cosines according to mu
    call rotate_angle(u, v, w, mu)

    ! Set energy and direction of particle in LAB frame
    p % E = E_new
    p % coord0 % uvw = v_n / vel

    ! Copy scattering cosine for tallies
    p % mu = mu

  end subroutine multigroup_scatter

!===============================================================================
! SAMPLE_GROUP samples the outgoing energy group, using the incident group and
! its associated scattering cross sections (stored in rxn).
!===============================================================================
  function sample_group(rxn, E) result(E_new)
  
  type(Reaction), pointer    :: rxn        ! sampled reaction (scattering)
  real(8),        intent(in) :: E          ! incoming energy bin
  real(8)         :: k = 0.0    ! cumulative probability
  real(8)         :: cutoff     ! random number
  integer         :: offset = 0 ! group search offset
  
  ! sample cutoff as fraction of total sigma_s for this energy group
  cutoff = prn() * rxn % adist % total_scatter(E)
  
  ! determine outgoing energy group
  E_new = max_scatter(E) - 1 
  do while(k < cutoff) 
    E_new = E_new + 1
    ! add this group's cross section
    k = k + rxn % sigma(rxn % group_index(E) + offset)
    offset = offset + 1
  end do
  
  end function sample_group
  
!===============================================================================
! SAMPLE_ANGLE samples the cosine of the angle between incident and exiting
! particle directions either from 32 equiprobable bins or from a tabular
! distribution.
!===============================================================================

  function sample_angle(rxn, E, E_new) result(mu)

    type(Reaction), pointer    :: rxn   ! reaction
    real(8),        intent(in) :: E     ! incoming energy group
    real(8),        intent(in) :: E_new ! outgoing energy group
    real(8)        :: xi      ! random number on [0,1)
    integer        :: index   ! index of distribution type for this g -> g' combo
    integer        :: interp  ! type of interpolation
    integer        :: type    ! angular distribution type
    integer        :: n       ! number of incoming energy bins
    integer        :: lc      ! location in data array
    integer        :: NP      ! number of points in cos distribution
    integer        :: k       ! index on cosine grid
    real(8)        :: r       ! interpolation factor on incoming energy
    real(8)        :: frac    ! interpolation fraction on cosine
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
    n = rxn % adist % n_data

    ! set index for this initial/final group combo (already determined)
    ! group_index(E) is the index for scattering OUT of group E into the
    ! highest possible energy group.
    ! (E_new - max_scatter(E)) gives the offset between the actual outgoing 
    ! energy group sampled and the highest possible energy group.
    index = rxn % group_index(E) + E_new - rxn % max_scatter(E)
    
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
      k = 1 + int(n *xi)

      ! calculate cosine
      mu0 = rxn % adist % data(lc + k)
      mu1 = rxn % adist % data(lc + k+1)
      mu = mu0 + (n * xi - k) * (mu1 - mu0)

    elseif (type == ANGLE_TABULAR) then
      ! NP = number of cum. probabilities listed
      NP = (n - 1)/2

      ! determine outgoing cosine bin
      xi = prn()
      c_k = rxn % adist % data(lc + 1)
      do k = 1, NP
        c_k1 = rxn % adist % data(lc + k+1)
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
  
end module multigroup_physics