module multigroup_header

  use constants,   only: MAX_FILE_LEN
  ! use endf_header, only: Tab1

  implicit none

!===============================================================================
! DISTANGLE contains data for a tabular secondary angle distribution whether it
! be tabular or 32 equiprobable cosine bins
!===============================================================================

  type DistAngle
    integer              :: n_data      ! # of cosine values per group
    real(8), allocatable :: energy(:)   ! incoming energy grid
    integer, allocatable :: type(:)     ! type of distribution
    integer, allocatable :: location(:) ! location of each table
    real(8), allocatable :: data(:)     ! angular distribution data
    ! angular dist data will be n_data cosine values for ANGLE_NLEG_EQUI type,
    ! or n_data cosine values + (n_data - 1) cum. probs for ANGLE_TABULAR type  
  end type DistAngle

!===============================================================================
! DISTENERGY contains data for a secondary energy distribution for all
! scattering laws
!===============================================================================

  type DistEnergy
    integer    :: law                 ! secondary distribution law
    ! type(Tab1) :: p_valid             ! probability of law validity
    real(8), allocatable :: data(:)   ! energy distribution data

    ! For reactions that may have multiple energy distributions such as (n.2n),
    ! this pointer allows multiple laws to be stored
    type(DistEnergy), pointer :: next => null()
  end type DistEnergy

!===============================================================================
! REACTION contains the multigroup cross-sections for a single reaction
!===============================================================================

  type Reaction
    integer :: MT                            ! ENDF MT value
    integer :: multiplicity                  ! Number of secondary particles released
    logical :: scatter_in_cm                 ! scattering system in center-of-mass?
    real(8), allocatable :: total_scatter(:) ! if scatter rxn, total for each grp 
    real(8), allocatable :: sigma(:)         ! Cross section values
    integer, allocatable :: group_index(:)   ! Location in P0 xs where group data begin
    integer, allocatable :: max_scatter(:)   ! max outscatter group
    integer, allocatable :: min_scatter(:)   ! min outscatter group
    logical :: has_angle_dist                ! Angle distribution present?
    logical :: has_energy_dist               ! Energy distribution present?
    type(DistAngle)           :: adist       ! Secondary angular distribution
    type(DistEnergy), pointer :: edist       ! Secondary energy distribution
  end type Reaction

!===============================================================================
! NUCLIDE contains all the data for multi-group cross sections.
!===============================================================================

  type Nuclide
    character(10) :: name    ! name of nuclide, e.g. 92235.03c
    integer       :: zaid    ! Z and A identifier, e.g. 92235
    integer       :: listing ! index in xs_listings
    real(8)       :: awr     ! weight of nucleus in neutron masses
    real(8)       :: kT      ! temperature in MeV (k*T)

    ! Energy group information
    integer :: n_group                      ! # of energy groups
    integer, allocatable :: grid_index(:)   ! pointers to union grid
    real(8), allocatable :: group_energy(:) ! energy values corresponding to xs
    real(8), allocatable :: group_width(:)  ! width values corresponding to xs
    real(8), allocatable :: group_mass(:)   ! mass values corresponding to xs

    ! Microscopic cross sections
    real(8), allocatable :: total(:)      ! total cross section
    real(8), allocatable :: scattering(:) ! P0 scattering
    real(8), allocatable :: fission(:)    ! fission
    real(8), allocatable :: nu_fission(:) ! neutron production
    real(8), allocatable :: absorption(:) ! absorption (MT > 100)
    real(8), allocatable :: heating(:)    ! heating

    ! Fission information
    logical :: fissionable         ! nuclide is fissionable?
    logical :: has_partial_fission ! nuclide has partial fission reactions?
    integer :: n_fission           ! # of fission reactions
    integer, allocatable :: index_fission(:) ! indices in reactions

    ! Total fission neutron emission
    integer :: nu_t_type
    real(8), allocatable :: nu_t_data(:)

    ! Prompt fission neutron emission
    integer :: nu_p_type
    real(8), allocatable :: nu_p_data(:)
    
    ! Chi (fission fraction) information
    integer :: chi_type
    real(8), allocatable :: chi_data(:)

    ! Unresolved resonance data-- not valid for multigroup
    logical                :: urr_present = .false.

    ! Reactions
    integer :: n_reaction ! # of reactions
    type(Reaction), pointer :: reactions(:) => null()

  end type Nuclide

!===============================================================================
! XSLISTING contains data read from a cross_sections.xml file
!===============================================================================

  type XsListing
    character(12) :: name       ! table name, e.g. 92235.70c
    character(12) :: alias      ! table alias, e.g. U-235.70c
    integer       :: type       ! type of table (cont-E neutron, S(A,b), etc)
    integer       :: zaid       ! ZAID identifier = 1000*Z + A
    integer       :: filetype   ! ASCII or BINARY
    integer       :: location   ! location of table within library
    integer       :: recl       ! record length for library
    integer       :: entries    ! number of entries per record
    real(8)       :: awr        ! atomic weight ratio (# of neutron masses)
    real(8)       :: kT         ! Boltzmann constant * temperature (MeV)
    logical       :: metastable ! is this nuclide metastable?
    character(MAX_FILE_LEN) :: path ! path to library containing table
  end type XsListing

!===============================================================================
! NUCLIDEMICROXS contains cached microscopic cross sections for a
! particular nuclide in the current group
!===============================================================================

  type NuclideMicroXS
    integer :: index_grid      ! index on nuclide energy grid
    integer :: last_E          ! last evaluated energy GROUP
    real(8) :: total           ! microscropic total xs
    real(8) :: scattering      ! microscopic P0 scattering xs
    real(8) :: absorption      ! microscopic absorption xs
    real(8) :: fission         ! microscopic fission xs
    real(8) :: nu_fission      ! microscopic production xs

  end type NuclideMicroXS

!===============================================================================
! MATERIALMACROXS contains cached macroscopic cross sections for the material a
! particle is traveling through
!===============================================================================

  type MaterialMacroXS
    real(8) :: total         ! macroscopic total xs
    real(8) :: scattering    ! macroscopic P0 scattering xs
    real(8) :: absorption    ! macroscopic absorption xs
    real(8) :: fission       ! macroscopic fission xs
    real(8) :: nu_fission    ! macroscopic production xs
  end type MaterialMacroXS

end module multigroup_header

