module multigroup

  use multigroup_header, only: Nuclide, Reaction, XsListing, &
                              DistEnergy
  use constants
  use endf,              only: reaction_name
  use error,             only: fatal_error, warning
  use fission,           only: nu_total
  use global
  use material_header,   only: Material
  use output,            only: write_message
  use set_header,        only: SetChar
  use string,            only: str_to_int, str_to_real, lower_case, to_str

  implicit none

  integer :: NXS(16)             ! Descriptors for mg XSS tables
  integer :: JXS(32)             ! Pointers into mg XSS tables
  real(8), allocatable :: XSS(:) ! Cross section data
  integer :: XSS_index           ! current index in mg XSS data

  private :: NXS
  private :: JXS
  private :: XSS

contains

!===============================================================================
! READ_MG_XS reads all the cross sections for the problem (if multi-group  
! treatment is specified) and stores them in arrays. S(a,b) data is assumed
! to be taken into account in the multi-group cross section set.
!===============================================================================

  subroutine read_mg_xs()
      
    integer :: i            ! index in materials array
    integer :: j            ! index over nuclides in material
    integer :: i_listing    ! index in xs_listings array
    integer :: i_nuclide    ! index in nuclides
    integer :: m            ! position for sorting
    integer :: temp_nuclide ! temporary value for sorting
    integer :: temp_table   ! temporary value for sorting
    character(12)  :: name  ! name of isotope, e.g. 92235.70m
    character(12)  :: alias ! alias of nuclide, e.g. U-235.70m
    type(Material),   pointer :: mat => null()
    type(Nuclide),    pointer :: nuc => null()
    type(SetChar) :: already_read

    ! allocate arrays for ACE table storage and cross section cache
    allocate(nuclides(n_nuclides_total))
    allocate(micro_xs(n_nuclides_total))

    ! ==========================================================================
    ! READ ALL ACE CROSS SECTION TABLES
    ! Loop over all files
    MATERIAL_LOOP: do i = 1, n_materials
      mat => materials(i)

      NUCLIDE_LOOP: do j = 1, mat % n_nuclides
        name = mat % names(j)

        if (.not. already_read % contains(name)) then
          i_listing = xs_listing_dict % get_key(name)
          i_nuclide = nuclide_dict % get_key(name)
          name  = xs_listings(i_listing) % name
          alias = xs_listings(i_listing) % alias

          ! Keep track of what listing is associated with this nuclide
          nuc => nuclides(i_nuclide)
          nuc % listing = i_listing

          ! Read the mg table into the appropriate entry on the nuclides
          ! array
          call read_mg_table(i_nuclide, i_listing)

          ! Add name and alias to dictionary
          call already_read % add(name)
          call already_read % add(alias)
        end if
      end do NUCLIDE_LOOP
  end subroutine read_mg_xs   
  
!===============================================================================
! READ_MG_TABLE reads a single cross section table in ASCII format.
! This routine reads the header data for each table and then calls
! appropriate subroutines to parse the actual data.
!===============================================================================

   subroutine read_mg_table(i_table, i_listing)

    integer, intent(in) :: i_table   ! index in nuclides/sab_tables
    integer, intent(in) :: i_listing ! index in xs_listings

    integer       :: i             ! loop index for XSS records
    integer       :: j, j1, j2     ! indices in XSS
    integer       :: record_length ! Fortran record length
    integer       :: location      ! location of mg table
    integer       :: entries       ! number of entries on each record
    integer       :: length        ! length of mg table
    integer       :: in = 7        ! file unit
    integer       :: zaids(16)     ! list of ZAIDs (only used for S(a,b))
    real(8)       :: kT            ! temperature of table
    real(8)       :: awrs(16)      ! list of atomic weight ratios (not used)
    real(8)       :: awr           ! atomic weight ratio for table
    logical       :: file_exists   ! does mg library exist?
    character(7)  :: readable      ! is mg library readable?
    character(10) :: name          ! name of mg table
    character(10) :: date_         ! date mg library was processed
    character(10) :: mat           ! material identifier
    character(70) :: comment       ! comment for mg table
    character(MAX_FILE_LEN) :: filename ! path to mg cross section library
    type(Nuclide),   pointer :: nuc => null()
    type(XsListing), pointer :: listing => null()

    ! determine path, record length, and location of table
    listing => xs_listings(i_listing)
    filename      = listing % path
    record_length = listing % recl
    location      = listing % location
    entries       = listing % entries

    ! Check if multigroup library exists and is readable
    inquire(FILE=filename, EXIST=file_exists, READ=readable)
    if (.not. file_exists) then
      message = "Multigroup library '" // trim(filename) // "' does not exist!"
      call fatal_error()
    elseif (readable(1:3) == 'NO') then
      message = "Multigroup library '" // trim(filename) // "' is not readable! &
           &Change file permissions with chmod command."
      call fatal_error()
    end if

    ! display message
    message = "Loading multigroup cross section table: " // listing % name
    call write_message(6)

    ! =======================================================================
    ! READ ACE TABLE IN ASCII FORMAT
    ! Find location of table
    open(UNIT=in, FILE=filename, STATUS='old', ACTION='read')
    rewind(UNIT=in)
    do i = 1, location - 1
      read(UNIT=in, FMT=*)
    end do

    ! Read first line of header
    read(UNIT=in, FMT='(A10,2G12.0,1X,A10)') name, awr, kT, date_

    ! Read more header and NXS and JXS
    read(UNIT=in, FMT=100) comment, mat, & 
         (zaids(i), awrs(i), i=1,16), NXS, JXS
100 format(A70,A10/4(I7,F11.0)/4(I7,F11.0)/4(I7,F11.0)/4(I7,F11.0)/&
        ,8I9/8I9/8I9/8I9/8I9/8I9)

    ! determine table length
    length = NXS(1)
    allocate(XSS(length))

    ! Read XSS array
    read(UNIT=in, FMT='(4G20.0)') XSS

    ! Close ACE file
    close(UNIT=in)

    ! ==========================================================================
    ! PARSE DATA BASED ON NXS, JXS, AND XSS ARRAYS

    nuc => nuclides(i_table)
    nuc % name = name
    nuc % n_group = NXS(5)
    nuc % awr  = awr
    nuc % kT   = kT
    nuc % zaid = NXS(2)

    ! read all blocks
    call read_erg(nuc)
    call read_nu_data(nuc)
    call read_chi_data(nuc)
    call read_abs_data(nuc)
    call read_reactions(nuc)
    call read_angular_dist(nuc)
    call read_energy_dist(nuc)

    ! Currently subcritical fixed source calculations are not allowed. Thus,
    ! if any fissionable material is found in a fixed source calculation,
    ! abort the run.
    if (run_mode == MODE_FIXEDSOURCE .and. nuc % fissionable) then
      message = "Cannot have fissionable material in a fixed source run."
      call fatal_error()
    end if

    ! for fissionable nuclides, precalculate microscopic nu-fission cross
    ! sections so that we don't need to call the nu_total function during
    ! cross section lookups

    if (nuc % fissionable) call generate_nu_fission(nuc)

    deallocate(XSS)
    if(associated(nuc)) nullify(nuc)

  end subroutine read_ace_table
  
!===============================================================================
! READ_ERG - reads through the ERG block. This block contains the energy grid,
! with group centers and masses.
!===============================================================================

  subroutine read_erg(nuc)

    type(Nuclide), pointer :: nuc

    integer :: NE ! number of energy groups

    ! determine number of energy groups
    NE = NXS(5)
    nuc % n_group = NE

    ! allocate storage for energy grid and cross section arrays
    allocate(nuc % group_energy(NE))
    allocate(nuc % group_width(NE))
    allocate(nuc % group_mass(NE))
    allocate(nuc % total(NE))
    allocate(nuc % scattering(NE))
    allocate(nuc % fission(NE))
    allocate(nuc % nu_fission(NE))
    allocate(nuc % absorption(NE))

    ! initialize cross sections
    nuc % total      = ZERO
    nuc % scattering = ZERO
    nuc % fission    = ZERO
    nuc % nu_fission = ZERO
    nuc % absorption = ZERO

    ! Read data from XSS -- energy grid first, then  then group masses
    XSS_index = 1
    nuc % group_energy = get_real(NE)
    ! Continue reading energy group widths
    nuc % group_width = get_real(NE)
    ! Read energy group masses
    nuc % group_mass = get_real(NE)

  end subroutine read_erg
  
!===============================================================================
! READ_NU_DATA reads data given on the number of neutrons emitted from fission
! as a function of the incoming group of a neutron. This data may be broken
! down into prompt and delayed neutrons based on the value of NXS(10)
!===============================================================================

  subroutine read_nu_data(nuc)

    type(Nuclide), pointer :: nuc

    integer :: i      ! loop index
    integer :: JXS3   ! location of fission cross sections (0 if not fissionable)
    integer :: JXS4   ! location for fission nu data
    integer :: NXS10  ! number of nu types (i.e. total, prompt) present
    integer :: length ! length of data to allocate

    JXS3  = JXS(3)
    JXS4  = JXS(4)
    NXS10 = NXS(10)
    length = nuc % n_group ! Number of energy groups
    nuc % nu_d_type = NU_NONE ! No delayed nu data for multigroup (only total and/or prompt)
    
    if (JXS3 == 0) then
      ! =======================================================================
      ! NO PROMPT/TOTAL NU DATA
      nuc % nu_t_type = NU_NONE
      nuc % nu_p_type = NU_NONE

    else if (JXS3 /= 0 .and. NXS10 == 1) then
      ! =======================================================================
      ! PROMPT OR TOTAL NU DATA
      ! Set data types
      nuc % nu_t_type = NU_TABULAR
      nuc % nu_p_type = NU_NONE
      ! allocate space for total nu data
      allocate(nuc % nu_t_data(length))
      ! Set index to read nu data
      XSS_index = JXS4
      ! Read total nu data
      nuc % nu_t_data = get_real(length)

    else if (JXS3 /= 0 .and. NXS10 == 2) then
      ! =======================================================================
      ! PROMPT AND TOTAL NU DATA
      ! set data types
      nuc % nu_t_type = NU_TABULAR
      nuc % nu_p_type = NU_TABULAR
      ! allocate space for total and prompt nu data
      allocate(nuc % nu_t_data(length))
      allocate(nuc % nu_p_data(length))
      ! Set index to read nu data
      XSS_index = JXS4
      ! Read prompt nu data
      nuc % nu_p_data = get_real(length)
      ! Read total nu data
      nuc % nu_t_data = get_real(length)
    end

  end subroutine read_nu_data
 
!===============================================================================
! READ_CHI_DATA reads multigroup fission fraction data if the nuclide is 
! fissionable.  The data should be pre-normalized to sum to unity.
!===============================================================================
 
  subroutine read_chi_data(nuc)
      
    type(Nuclide), pointer :: nuc

    integer :: JXS3   ! location of fission cross sections (0 if not fissionable)
    integer :: JXS5   ! location for fission nu data
    integer :: length ! length of data to allocate
    
    JXS5 = JXS(5) ! location of chi data
    length = nuc % n_group
    
    if (JXS3 == 0) then
      ! =======================================================================
      ! NO CHI DATA (nuclide not fissionable)
      nuc % chi_type = CHI_NONE
    else
      ! =======================================================================
      ! TABULAR CHI DATA
      ! set type of chi value
      nuc % chi_type = CHI_TABULAR
      ! Allocate space for chi values
      allocate(nuc % chi_data(length))
      ! set index to read chi values
      XSS_index = JXS5
      ! read chi data
      nuc % chi_data = get_real(length)
    end if  
  end subroutine
  
!===============================================================================
! READ_REACTIONS currently reads the fission, absorption, and total cross
! sections for the nuclide.  It may be expanded to include other edit reactions.
!===============================================================================
  
  subroutine read_reactions(nuc)
      
      type(Nuclide), pointer :: nuc
      
      integer :: JXS2    ! location of total cross sections
      integer :: JXS3    ! location of fission cross sections
      integer :: JXS6    ! location of absorption cross section data
      integer :: length  ! length of data to allocate
      type(Reaction), pointer :: rxn => null()
  
      JXS2 = JXS(2)
      JXS3 = JXS(3)
      JXS6 = JXS(6)
      length = nuc % n_group
      
      ! allocate 3 reactions: fission, absorption, scattering
      nuc % n_reaction = 3
      allocate(nuc % reactions(3)) 
      
      nuc % fissionable = .false.

      if (JXS2 == 0) then 
        ! =======================================================================
        ! NO TOTAL XS DATA
        ! what do we do here?
      else 
        ! =======================================================================
        ! TABULAR TOTAL XS DATA
        rxn => nuc % reactions(1)

        ! set defaults
        rxn % has_angle_dist  = .false.
        rxn % MT = TOTAL_XS
        
        ! allocate space
        allocate(nuc % total(length))
        ! set index
        rxn
        XSS_index = JXS2
        ! read total x section data
        nuc % total = get_real(length)
      endif
      
      if (JXS3 == 0) then 
        ! =======================================================================
        ! NO FISSION XS DATA
      else 
        ! =======================================================================
        ! TABULAR FISSION XS DATA
        
        rxn => nuc % reactions(2)
        rxn % has_angle_dist  = .false.
        rxn % MT = N_FISSION
        
        nuc % fissionable = .true.
        ! allocate space
        allocate(nuc % fission(length))  
        ! set index
        XSS_index = JXS3
        ! read fission x section data
        nuc % fission = get_real(length)
      endif
      
      if (JXS6 == 0) then 
        ! =======================================================================
        ! NO ABS DATA
          ! what do we do here?
      else 
        ! =======================================================================
        ! TABULAR ABS DATA
        
        rxn => nuc % reactions(3)
        rxn % has_angle_dist  = .false.
        rxn % MT = TOTAL_XS
        ! allocate space
        allocate(nuc % absorption(length))
        ! set index
        XSS_index = JXS6
        ! read abs x section data
        nuc % absorption = get_real(length)
      endif

  end subroutine read_reactions(nuc)
  
!===============================================================================
! READ_SCATTERING reads in the P0 scattering matrix for each nuclide, as well as
! the angular distributions for exiting neutrons.  
!===============================================================================
  
  subroutine read_scattering(nuc)
    integer :: NXS3                 ! number of angular distribution variables
    integer :: NXS5                 ! total number of energy groups
    integer :: NXS6                 ! number of upscatter groups
    integer :: NXS7                 ! number of downscatter groups
    integer :: NXS8                 ! number of secondary particles
    integer :: NXS9                 ! type of angular distribution (EPB/discrete cosine)
    integer :: NXS12                ! incident particle identifier (used for check only)
    integer :: JXS11                ! location of secondary particle types
    integer :: JXS12                ! location of secondary group structure locators
    integer :: JXS13                ! location of P0 scattering table locators 
    integer :: JXS14                ! location of ang distribution types
    integer :: length               ! length of records to grab
    integer :: P0_index             ! location of incident particle P0 block
    integer :: PN_index             ! location of incident particle PN block
    integer, allocatable :: XPN(:)  ! locations for angular dist. blocks
    
    NXS5 = NXS(5)
    NXS6 = NXS(6)
    NXS7 = NXS(7)
    NXS9 = NXS(9)
    JXS13 = JXS(13)
    JXS16 = JXS(16)
    
    ! start by only reading/considering initial (incident) particle P0 scattering
    if(JXS13 == 0) then
      !=========================================================================
      ! NO P0 SCATTERING DATA
    else
      !=========================================================================
      ! TABULAR P0 SCATTERING DATA
      ! set length of table -- see LANL report LA-12704, pg 119
      length = int(NXS5 * (1 + NXS7 + NXS6) - 0.5 * (NXS7 * (NXS7 + 1) + NXS6 * (NXS6 + 1)))
      allocate(nuc % scattering(length))
      ! set XSS_index to location of incident particle PO_block
      XSS_index = JXS13
      P0_index = get_int(1)
      XSS_index = P0_index
      ! read correct number of values
      nuc % scattering = get_real(length)
      
      ! get location for incident particle XPN block (if one exists)
      XSS_index = JXS16
      PN_index = get_int(1)
      if(PN_index > 0) then
        XSS_index = PN_index
        allocate(nuc % )
      else
        !=======================================================================
        ! ALL SCATTERING IS ISOTROPIC-- DO NOTHING
      end if
    end if
    
  end subroutine read_scattering(nuc)
  
end module multigroup