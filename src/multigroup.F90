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
! READ_XS reads all the cross sections for the problem (if multi-group  
! treatment is specified) and stores them in arrays. S(a,b) data is assumed
! to be taken into account in the multi-group cross section set.
!===============================================================================

  subroutine read_xs()
      
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
    end do MATERIAL_LOOP
  end subroutine read_xs   
  
!===============================================================================
! READ_MG_TABLE reads a single cross section table in ASCII format.
! This routine reads the header data for each table and then calls
! appropriate subroutines to parse the actual data.
!===============================================================================

   subroutine read_mg_table(i_table, i_listing)

    integer, intent(in) :: i_table   ! index in nuclides
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
    call read_reactions(nuc)
    call read_scattering(nuc)

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

  end subroutine read_mg_table
  
!===============================================================================
! READ_ERG - reads through the ERG block. This block contains the energy grid,
! with group centers and masses.
!===============================================================================

  subroutine read_erg(nuc)

    type(Nuclide), pointer :: nuc

    integer :: NG ! number of energy groups

    ! determine number of energy groups
    NG = NXS(5)
    nuc % n_group = NG

    ! allocate storage for energy grid and cross section arrays
    ! skip scattering for now!
    allocate(nuc % group_energy(NG))
    allocate(nuc % group_width(NG))
    allocate(nuc % group_mass(NG))
    
    allocate(nuc % total(NG))
    allocate(nuc % fission(NG))
    allocate(nuc % nu_fission(NG))
    allocate(nuc % absorption(NG))

    ! initialize cross sections
    nuc % total      = ZERO
    nuc % fission    = ZERO
    nuc % nu_fission = ZERO
    nuc % absorption = ZERO

    ! Read data from XSS -- energy grid first, then  then group masses
    XSS_index = 1
    nuc % group_energy = get_real(NG)
    ! Continue reading energy group widths
    nuc % group_width = get_real(NG)
    ! Read energy group masses
    nuc % group_mass = get_real(NG)

  end subroutine read_erg
  
!===============================================================================
! READ_NU_DATA reads data given on the number of neutrons emitted from fission
! as a function of the incoming group of a neutron. This data may be broken
! down into prompt and total data depending on value of NXS10
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
    end if

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
  end subroutine read_chi_data
  
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
      ! remember, 'total' isn't a reaction, dummy! :)
      nuc % n_reaction = 3
      allocate(nuc % reactions(3)) 
      
      ! default value
      nuc % fissionable = .false.

      if (JXS2 == 0) then 
        ! =======================================================================
        ! NO TOTAL XS DATA
      else 
        ! =======================================================================
        ! TABULAR TOTAL XS DATA
        ! set index
        XSS_index = JXS2
        ! read total x section data (space already allocated)
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
        nuc % has_partial_fission = .false.
        nuc % n_fission = 1
        ! allocate space
        allocate(rxn % sigma(length)) 
        allocate(nuc % index_fission(1))
        ! set index
        XSS_index = JXS3
        ! read fission x section data
        rxn % sigma = get_real(length)
        nuc % fission = rxn % sigma
      endif
      
      if (JXS6 == 0) then 
        ! =======================================================================
        ! NO ABS DATA
      else 
        ! =======================================================================
        ! TABULAR ABS DATA
        
        rxn => nuc % reactions(3)
        rxn % has_angle_dist  = .false.
        rxn % MT = N_DISAPPEAR
        ! allocate space
        allocate(rxn % sigma(length))
        ! set index
        XSS_index = JXS6
        ! read abs x section data
        rxn % sigma = get_real(length)
        nuc % absorption = rxn % sigma
      endif
      
      ! skip reading secondary particle reactions for now! :)

  end subroutine read_reactions
  
!===============================================================================
! READ_SCATTERING reads in the P0 scattering matrix for each nuclide, as well as
! the angular distributions for exiting neutrons.  
!===============================================================================
  
  subroutine read_scattering(nuc)
    
    type(Nuclide), pointer :: nuc
    
    integer :: NLEG                 ! number of angular distribution variables
    integer :: NXS5                 ! total number of energy groups
    integer :: NXS6                 ! number of upscatter groups
    integer :: NXS7                 ! number of downscatter groups
    ! integer :: NXS8                 ! number of secondary particles
    integer :: ISANG                ! type of angular distribution (EPB/discrete cosine)
    integer :: NXS12                ! incident particle identifier (used for check only)
    ! integer :: JXS11                ! location of secondary particle types
    ! integer :: JXS12                ! location of secondary group structure locators
    integer :: JXS13                ! location of P0 scattering table locators 
    integer :: JXS14                ! location of ang distribution types
    integer :: JXS16                ! location of XPN block locators
    integer :: JXS17                ! location of PN block locators
    integer :: length               ! length of records to grab
    integer :: i                    ! scattering group iterator
    integer :: j                    ! iterator
    integer :: k                    ! scattering distribution iterator
    integer :: LP0                  ! location of incident particle P0 block
    integer :: LPN                  ! location of incident particle PN block
    integer :: LXPN                 ! location of incident particle XPN block
    integer :: base_type            ! base ang. dist. type for this nuclide
    integer :: this_index           ! index where data for this grp begins
    integer :: I1                   ! highest possible group for outscatter
    integer :: I2                   ! lowest possible group for outscatter
    integer, allocatable :: LPND(:)  ! locations of PND angular dist. blocks
    type(Reaction), pointer :: rxn => null()
          
    NLEG = NXS(3)
    NXS5 = NXS(5)
    NXS6 = NXS(6)
    NXS7 = NXS(7)
    ISANG = NXS(9)
    JXS13 = JXS(13)
    JXS16 = JXS(16)
    JXS17 = JXS(17)
    
    ! associate scattering with reaction 1
    rxn => nuc % reactions(1)
    
    ! MG xsections are in lab frame
    rxn % scatter_in_cm = .false.
    
    ! start by only reading/considering initial (incident) particle scattering
    if(JXS13 == 0) then
      !=========================================================================
      ! NO P0 SCATTERING DATA AT ALL -- MUST BE PURE ABSORBER
      rxn % has_angle_dist = .false.
      
    else
      !=========================================================================
      ! TABULAR P0 SCATTERING DATA
      ! allocate group indices, total scattering, min/max scatter groups  
      allocate(rxn % group_index(nuc % n_group))
      allocate(rxn % total_scatter(nuc % n_group))
      allocate(rxn % min_scatter(nuc % n_group))
      allocate(rxn % max_scatter(nuc % n_group))
      
      ! initialize to zeros
      rxn % group_index = ZERO
      rxn % total_scatter = ZERO
      rxn % min_scatter = ZERO
      rxn % max_scatter = ZERO
      
      ! set length of table -- see LANL report LA-12704, pg 119
      length = int(NXS5 * (1 + NXS7 + NXS6) - 0.5 * (NXS7 * (NXS7 + 1) + NXS6 &
               * (NXS6 + 1)))
      
      ! allocate scattering cross sections, reaction cross sections
      allocate(nuc % scattering(length))
      allocate(rxn % sigma(length))
      
      ! set XSS_index to location of incident particle PO xsects
      XSS_index = JXS13
      LP0 = int(XSS(XSS_index))
      XSS_index = LP0
      
      ! read correct number of P0 values, copy to reaction
      rxn % sigma = get_real(length)
      nuc % scattering = rxn % sigma
      
      ! sum up total P0 values (total scattering for each group)
      ! also store array index where each group's outscatter xs start
      this_index = 1
      do k=1,nuc % n_group
        ! store starting index for group
        rxn % group_index(k) = this_index
        
        ! determine min, max groups for outscatter
        rxn % max_scatter(k) = max(1, k - NXS6)
        rxn % min_scatter(k) = min(NXS5, k + NXS7)
        
        ! sum up outscatter sigmas for this energy group 
        do i=this_index,(this_index + rxn % min_scatter(k) - rxn % max_scatter(k))
          rxn % total_scatter(k) = rxn % total_scatter(k) + rxn % sigma(i)
        end do
        
        ! update starting index for next group
        this_index = this_index + (rxn % min_scatter(k) - rxn % max_scatter(k))
      end do
      
      ! get location for incident particle XPN block (if one exists)
      XSS_index = JXS16
      LXPN = int(XSS(XSS_index))
      
      ! get location for incident particle PN block (if one exists)
      XSS_index = JXS17
      LPN = int(XSS(XSS_index))

      if(LXPN > 0 .and. LPN > 0) then
        ! at least one ang distribution exists
        rxn % has_angle_dist = .true.
        allocate(rxn % adist % location(length))
        allocate(rxn % adist % type(length))
        
        ! allocate, then read in locations for PND blocks
        XSS_index = LXPN
        allocate(LPND(length))
        LPND = get_real(length)
        
        ! determine base distribution type
        if(ISANG == 1) then
          base_type = ANGLE_TABULAR
        else if(ISANG == 0) then
          base_type = ANGLE_NLEG_EQUI
        else
          message = "Invalid angular distribution (ISANG = " // trim(to_str   &
                    (ISANG)) // ")."
          call fatal_error()
        end if
        
        ! set # of data points per angular distribution
        rxn % adist % n_data = NLEG

        ! allocate space for all possible g->g' combo angular distributions
        allocate(rxn % adist % data(length * NLEG))
        
        do j = 1,length

          ! determine location (offset) where group data is located
          rxn % adist % location(j) = LPND(j)
          
          ! If this LPND = -1 or 0, no data in this block
          ! Set scattering type to isotropic, then cycle
          if(LPND(j) <= 0) then
            rxn % adist % type(j) = ANGLE_ISOTROPIC
            cycle
          end if
          
          ! We made it here-- there must be data!
          rxn % adist % type(j) = base_type
          XSS_index = LPN + LPND(j)
          rxn % adist % data((j - 1) * NLEG + 1:j * NLEG) = get_real(NLEG)
          
        end do
      else
        !=======================================================================
        ! ALL SCATTERING IS ISOTROPIC IN LAB SYSTEM-- DO NOTHING
        rxn % has_angle_dist = .false.
      end if
    end if
    
  end subroutine read_scattering
  
!===============================================================================
! GENERATE_NU_FISSION precalculates the microscopic nu-fission cross section for
! a given nuclide. This is done so that the nu_total function does not need to
! be called during cross section lookups.
!===============================================================================

  subroutine generate_nu_fission(nuc)

    type(Nuclide), pointer :: nuc

    integer :: i  ! nuclide energy group index
    real(8) :: nu ! # of neutrons per fission

    do i = 1, nuc % n_group

      ! determine total nu at given energy
      nu = nu_total(nuc, i)

      ! determine nu-fission microscopic cross section
      nuc % nu_fission(i) = nu * nuc % fission(i)
    end do

  end subroutine generate_nu_fission
  
!===============================================================================
! GET_INT returns an array of integers read from the current position in the XSS
! array
!===============================================================================

  function get_int(n_values) result(array)

    integer, intent(in) :: n_values        ! number of values to read
    integer             :: array(n_values) ! array of values

    array = int(XSS(XSS_index:XSS_index + n_values - 1))
    XSS_index = XSS_index + n_values

  end function get_int

!===============================================================================
! GET_REAL returns an array of real(8)s read from the current position in the
! XSS array
!===============================================================================

  function get_real(n_values) result(array)

    integer, intent(in) :: n_values        ! number of values to read
    real(8)             :: array(n_values) ! array of values

    array = XSS(XSS_index:XSS_index + n_values - 1)
    XSS_index = XSS_index + n_values

  end function get_real
  
end module multigroup