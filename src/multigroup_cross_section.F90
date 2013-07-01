module multigroup_cross_section

  use multigroup_header, only: Nuclide, Reaction
  use constants
  use error,             only: fatal_error
  use fission,           only: nu_total
  use global
  use material_header,   only: Material
  use random_lcg,        only: prn
  use search,            only: binary_search

  implicit none

contains

!===============================================================================
! CALCULATE_XS determines the macroscopic cross sections for the material the
! particle is currently traveling through.
!===============================================================================

  subroutine calculate_xs()

    integer :: i             ! loop index over nuclides
    integer :: i_nuclide     ! index into nuclides array
    real(8) :: atom_density  ! atom density of a nuclide
    type(Material),  pointer :: mat => null() ! current material

    ! Set all material macroscopic cross sections to zero
    material_xs % total      = ZERO
    material_xs % scattering = ZERO
    material_xs % absorption = ZERO
    material_xs % fission    = ZERO
    material_xs % nu_fission = ZERO

    ! Exit subroutine if material is void
    if (p % material == MATERIAL_VOID) return

    mat => materials(p % material)

    ! Add contribution from each nuclide in material
    do i = 1, mat % n_nuclides
      
      ! ========================================================================
      ! CALCULATE MICROSCOPIC CROSS SECTION

      ! Determine microscopic cross sections for this nuclide
      i_nuclide = mat % nuclide(i)

      ! Calculate microscopic cross section for this nuclide
      if (p % E /= micro_xs(i_nuclide) % last_E) then
        call calculate_nuclide_xs(i_nuclide)
      end if

      ! ========================================================================
      ! ADD TO MACROSCOPIC CROSS SECTION

      ! Copy atom density of nuclide in material
      atom_density = mat % atom_density(i)

      ! Add contributions to material macroscopic total cross section
      material_xs % total = material_xs % total + &
           atom_density * micro_xs(i_nuclide) % total

      ! Add contributions to material macroscopic scattering cross section
      material_xs % scattering = material_xs % scattering + &
           atom_density * micro_xs(i_nuclide) % scattering
           
      ! Add contributions to material macroscopic absorption cross section
      material_xs % absorption = material_xs % absorption + & 
           atom_density * micro_xs(i_nuclide) % absorption

      ! Add contributions to material macroscopic fission cross section
      material_xs % fission = material_xs % fission + &
           atom_density * micro_xs(i_nuclide) % fission

      ! Add contributions to material macroscopic nu-fission cross section
      material_xs % nu_fission = material_xs % nu_fission + &
           atom_density * micro_xs(i_nuclide) % nu_fission
 
    end do

  end subroutine calculate_xs
  
!===============================================================================
! CALCULATE_NUCLIDE_XS determines microscopic cross sections for a nuclide of a
! given index in the nuclides array in a particular energy group
!===============================================================================

  subroutine calculate_nuclide_xs(i_nuclide)

    integer, intent(in) :: i_nuclide ! index into nuclides array

!    integer :: i_group ! index on nuclide energy grid
    type(Nuclide),   pointer :: nuc => null()

    ! Set pointer to nuclide
    nuc => nuclides(i_nuclide)

    micro_xs(i_nuclide) % index_grid    = p % E

    ! Initialize nuclide cross-sections to zero
    micro_xs(i_nuclide) % fission    = ZERO
    micro_xs(i_nuclide) % nu_fission = ZERO

    ! Calculate microscopic nuclide total cross section
    micro_xs(i_nuclide) % total = nuc % total(p % E) 

    ! Calculate microscopic nuclide total cross section
    micro_xs(i_nuclide) % scattering = nuc % scattering(p % E) 

    ! Calculate microscopic nuclide absorption cross section
    micro_xs(i_nuclide) % absorption = nuc % absorption( &
         p % E) 

    if (nuc % fissionable) then
      ! Calculate microscopic nuclide total cross section
      micro_xs(i_nuclide) % fission = nuc % fission(p % E) 

      ! Calculate microscopic nuclide nu-fission cross section
      micro_xs(i_nuclide) % nu_fission =  nuc % nu_fission(p % E) 
    end if

      micro_xs(i_nuclide) % last_E = p % E
    
  end subroutine calculate_nuclide_xs

end module multigroup_cross_section
