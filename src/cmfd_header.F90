module cmfd_header 

  use constants,  only: CMFD_NOACCEL, ZERO, ONE

  implicit none
  private
  public :: allocate_cmfd, allocate_funct, allocate_no_accum, deallocate_cmfd,&
            deallocate_funct, deallocate_no_accum

  type, public :: cmfd_type

    ! indices for problem
    integer :: indices(4)

    ! albedo boundary condition
    real(8) :: albedo(6)

    ! core overlay map
    integer, allocatable :: coremap(:,:,:)
    integer, allocatable :: indexmap(:,:)
    integer :: mat_dim = CMFD_NOACCEL 

    ! energy grid
    real(8), allocatable :: egrid(:)

    ! cross sections
    real(8), allocatable :: totalxs(:,:,:,:)
    real(8), allocatable :: p1scattxs(:,:,:,:)
    real(8), allocatable :: scattxs(:,:,:,:,:)
    real(8), allocatable :: nfissxs(:,:,:,:,:)

    ! diffusion coefficient
    real(8), allocatable :: diffcof(:,:,:,:)

    ! current
    real(8), allocatable :: current(:,:,:,:,:)
    
    ! edge mu-sq weighted flux
    real(8), allocatable :: mu_sq(:,:,:,:,:)

    ! flux
    real(8), allocatable :: flux(:,:,:,:)
    
    ! cell-average weighted current (if use_functs)
    real(8), allocatable :: vol_curr(:,:,:,:,:)

    ! coupling coefficients and equivalence parameters
    real(8), allocatable :: dtilde(:,:,:,:,:)
    real(8), allocatable :: dhat(:,:,:,:,:)

    ! dimensions of mesh cells ([hu,hv,hw],xloc,yloc,zloc)
    real(8), allocatable :: hxyz(:,:,:,:)

    ! source distributions
    real(8), allocatable :: cmfd_src(:,:,:,:)
    real(8), allocatable :: openmc_src(:,:,:,:)

    ! source sites in each mesh box
    real(8), allocatable :: sourcecounts(:,:,:,:)

    ! weight adjustment factors 
    real(8), allocatable :: weightfactors(:,:,:,:)

    ! eigenvector/eigenvalue from cmfd run
    real(8), allocatable :: phi(:)
    real(8) :: keff = ZERO
    
    ! accumulated eigenvector/eigenvector_sq from cmfd run
    ! (only used if .not. cmfd_accum)
    real(8), allocatable :: phi_sum(:)
    real(8), allocatable :: phi_sum_sq(:)

    ! FINAL eigenfunction from cmfd run
    real(8), allocatable :: phi_final(:)

    ! eigenvector/eigenvalue from adjoint run
    real(8), allocatable :: adj_phi(:)
    real(8) :: adj_keff = ZERO

    ! residual for neutron balance
    real(8), allocatable :: resnb(:,:,:,:)

    ! openmc source normalization factor
    real(8) :: norm = ONE

    ! Shannon entropy from cmfd fission source
    real(8) :: entropy 

  end type cmfd_type

contains

!==============================================================================
! ALLOCATE_CMFD
!==============================================================================

  subroutine allocate_cmfd(this)  
    
    type(cmfd_type) :: this 

    integer :: nx  ! number of mesh cells in x direction
    integer :: ny  ! number of mesh cells in y direction
    integer :: nz  ! number of mesh cells in z direction
    integer :: ng  ! number of energy groups

   ! extract spatial and energy indices from object
    nx = this % indices(1)
    ny = this % indices(2)
    nz = this % indices(3)
    ng = this % indices(4)

    ! allocate flux, cross sections and diffusion coefficient
    if (.not. allocated(this % flux))       allocate(this % flux(ng,nx,ny,nz))
    if (.not. allocated(this % totalxs))    allocate(this % totalxs(ng,nx,ny,nz))
    if (.not. allocated(this % p1scattxs))  allocate(this % p1scattxs(ng,nx,ny,nz))
    if (.not. allocated(this % scattxs))    allocate(this % scattxs(ng,ng,nx,ny,nz))
    if (.not. allocated(this % nfissxs))    allocate(this % nfissxs(ng,ng,nx,ny,nz))
    if (.not. allocated(this % diffcof))    allocate(this % diffcof(ng,nx,ny,nz))

    ! allocate dtilde and dhat 
    if (.not. allocated(this % dtilde))     allocate(this % dtilde(6,ng,nx,ny,nz))
    if (.not. allocated(this % dhat))       allocate(this % dhat(6,ng,nx,ny,nz))

    ! allocate dimensions for each box (here for general case)
    if (.not. allocated(this % hxyz))       allocate(this % hxyz(3,nx,ny,nz))

    ! allocate surface currents
    if (.not. allocated(this % current))    allocate(this % current(12,ng,nx,ny,nz))

    ! allocate source distributions
    if (.not. allocated(this % cmfd_src)) allocate(this % cmfd_src(ng,nx,ny,nz))
    if (.not. allocated(this % openmc_src)) allocate(this % openmc_src(ng,nx,ny,nz))
    
    ! allocate final flux distribution
    if (.not. allocated(this % phi_final)) allocate(this % phi_final(ng*nx*ny*nz))

    ! allocate source weight modification vars
    if (.not. allocated(this % sourcecounts)) allocate(this % sourcecounts(ng,nx,ny,nz))
    if (.not. allocated(this % weightfactors)) allocate(this % weightfactors(ng,nx,ny,nz))

    ! set everthing to 0 except weight multiply factors if feedback isnt on
    this % flux          = ZERO
    this % totalxs       = ZERO
    this % p1scattxs     = ZERO
    this % scattxs       = ZERO
    this % nfissxs       = ZERO
    this % diffcof       = ZERO
    this % dtilde        = ZERO
    this % dhat          = ZERO
    this % hxyz          = ZERO
    this % current       = ZERO
    this % cmfd_src      = ZERO
    this % openmc_src    = ZERO
    this % sourcecounts  = ZERO
    this % weightfactors = ONE
    this % phi_final     = ZERO

  end subroutine allocate_cmfd

!==============================================================================
! ALLOCATE_FUNCT allocates special variables used only for functional calcs
!==============================================================================

  subroutine allocate_funct(this)
    
    type(cmfd_type) :: this 

    integer :: nx  ! number of mesh cells in x direction
    integer :: ny  ! number of mesh cells in y direction
    integer :: nz  ! number of mesh cells in z direction
    integer :: ng  ! number of energy groups

   ! extract spatial and energy indices from object
    nx = this % indices(1)
    ny = this % indices(2)
    nz = this % indices(3)
    ng = this % indices(4)
    
    ! allocate variables for functional calculation
    if (.not. allocated(this % vol_curr)) allocate(this % vol_curr(3,ng,nx,ny,nz)) 
    ! All six faces of a 2-D cell; 3 direction cosines each face :)
    if (.not. allocated(this % mu_sq))    allocate(this % mu_sq(18,ng,nx,ny,nz)) 
    
    ! finally, initialize!
    this % vol_curr = ZERO
    this % mu_sq = ZERO
    
  end subroutine allocate_funct

!==============================================================================
! ALLOCATE_NO_ACCUM allocates arrays used only when tally accumulation is 
! turned OFF (these arrays are used to calc the CMFD RSDs)
!==============================================================================
  
  subroutine allocate_no_accum(this)
    
    type(cmfd_type) :: this 

    integer :: nx  ! number of mesh cells in x direction
    integer :: ny  ! number of mesh cells in y direction
    integer :: nz  ! number of mesh cells in z direction
    integer :: ng  ! number of energy groups

   ! extract spatial and energy indices from object
    nx = this % indices(1)
    ny = this % indices(2)
    nz = this % indices(3)
    ng = this % indices(4)
    
    ! These hold running phi sums
    if (.not. allocated(this % phi_sum)) allocate(this % phi_sum(ng*nx*ny*nz))
    ! These hold running phi sums squared
    if (.not. allocated(this % phi_sum_sq)) allocate(this % phi_sum_sq(ng*nx*ny*nz))    
    
    ! finally, initialize!
    this % phi_sum = ZERO
    this % phi_sum_sq = ZERO
    
  end subroutine allocate_no_accum
  
!===============================================================================
! DEALLOCATE_CMFD 
!===============================================================================

  subroutine deallocate_cmfd(this)

    type(cmfd_type) :: this

    if (allocated(this % egrid))         deallocate(this % egrid)
    if (allocated(this % totalxs))       deallocate(this % totalxs)
    if (allocated(this % p1scattxs))     deallocate(this % p1scattxs)
    if (allocated(this % scattxs))       deallocate(this % scattxs)
    if (allocated(this % nfissxs))       deallocate(this % nfissxs)
    if (allocated(this % diffcof))       deallocate(this % diffcof)
    if (allocated(this % current))       deallocate(this % current)
    if (allocated(this % flux))          deallocate(this % flux)
    if (allocated(this % dtilde))        deallocate(this % dtilde)
    if (allocated(this % dhat))          deallocate(this % dhat)
    if (allocated(this % hxyz))          deallocate(this % hxyz)
    if (allocated(this % coremap))       deallocate(this % coremap)
    if (allocated(this % indexmap))      deallocate(this % indexmap)
    if (allocated(this % phi))           deallocate(this % phi)
    if (allocated(this % phi_final))     deallocate(this % phi_final)
    if (allocated(this % sourcecounts))  deallocate(this % sourcecounts)
    if (allocated(this % weightfactors)) deallocate(this % weightfactors)
    if (allocated(this % cmfd_src))      deallocate(this % cmfd_src)
    if (allocated(this % openmc_src))    deallocate(this % openmc_src)
    
  end subroutine deallocate_cmfd

!===============================================================================
! DEALLOCATE_CMFD 
!===============================================================================

  subroutine deallocate_funct(this)

    type(cmfd_type) :: this    
    
    if (allocated(this % vol_curr))      deallocate(this % vol_curr) 
    if (allocated(this % mu_sq))         deallocate(this % mu_sq) 
    
  end subroutine deallocate_funct
  
!===============================================================================
! DEALLOCATE_NO_ACCUM 
!===============================================================================

  subroutine deallocate_no_accum(this)

    type(cmfd_type) :: this    
    
    if (allocated(this % phi_sum))      deallocate(this % phi_sum) 
    if (allocated(this % phi_sum_sq))   deallocate(this % phi_sum_sq) 
    
  end subroutine deallocate_no_accum
  
end module cmfd_header
