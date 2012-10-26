! added by S. Karra 07/11/12

module Secondary_Continuum_module

  implicit none

  private

#include "definitions.h"

  type, public :: slab_type
    PetscReal :: length                       ! input - length of slab
    PetscReal :: area                         ! input - surface area
  end type slab_type
  
  type, public :: nested_cube_type
    PetscReal :: matrix_block_size            ! input - side of cube
    PetscReal :: fracture_spacing             ! input - fracture spacing
  end type nested_cube_type
  
  type, public :: nested_sphere_type
    PetscReal :: radius                       ! input - radius of sphere
  end type nested_sphere_type
  
  type, public :: sec_continuum_type
    PetscInt :: itype                         ! input - type of sec. continuum (slab, nested_cube, nested_sphere,....) 
    type(slab_type) :: slab
    type(nested_cube_type) :: nested_cube
    type(nested_sphere_type) :: nested_sphere 
  end type sec_continuum_type

  type, public :: sec_heat_type  
    PetscBool :: sec_temp_update               ! flag to check if the temp is updated
    PetscInt :: ncells                         ! number of secondary grid cells
    PetscReal :: aperture                      ! fracture aperture
    PetscReal :: epsilon                       ! vol. frac. of primary continuum
    type(sec_continuum_type) :: sec_continuum
    PetscReal, pointer :: sec_temp(:)          ! array of temp. at secondary grid cells
    PetscReal, pointer :: area(:)              ! surface area
    PetscReal, pointer :: vol(:)               ! volume     face      node       face
    PetscReal, pointer :: dm_plus(:)           ! see fig.    |----------o----------|
    PetscReal, pointer :: dm_minus(:)          ! see fig.      <dm_minus> <dm_plus>
    PetscReal :: interfacial_area              ! interfacial area between prim. and sec. per unit volume of prim.+sec.
    PetscBool :: log_spacing                   ! flag to check if log spacing is set
    PetscReal :: outer_spacing                 ! value of the outer most grid cell spacing
  end type sec_heat_type  
  
  type, public :: sec_transport_type  
    PetscBool :: sec_conc_update               ! flag to check if the temp is updated
    PetscInt :: ncells                         ! number of secondary grid cells
    PetscReal :: aperture                      ! fracture aperture
    PetscReal :: epsilon                       ! vol. frac. of primary continuum
    type(sec_continuum_type) :: sec_continuum
    PetscReal, pointer :: sec_conc(:)          ! array of aqueous species conc. at secondary grid cells
    PetscReal, pointer :: sec_mnrl_volfrac(:)  ! array of mineral vol fraction at secondary grid cells
    PetscInt, pointer :: sec_zeta(:)           ! array of zetas at secondary grid cells
    PetscReal, pointer :: area(:)              ! surface area
    PetscReal, pointer :: vol(:)               ! volume     face      node       face
    PetscReal, pointer :: dm_plus(:)           ! see fig.    |----------o----------|
    PetscReal, pointer :: dm_minus(:)          ! see fig.      <dm_minus> <dm_plus>
    PetscReal :: interfacial_area              ! interfacial area between prim. and sec. per unit volume of prim.+sec.
    PetscReal :: sec_mnrl_area                 ! secondary mineral surface area
    PetscBool :: log_spacing                   ! flag to check if log spacing is set
    PetscReal :: outer_spacing                 ! value of the outer most grid cell spacing
  end type sec_transport_type  


  public :: SecondaryContinuumType, &
            SecondaryContinuumSetProperties
            
  contains

! ************************************************************************** !
!
! SecondaryContinuumType: The area, volume, grid sizes for secondary continuum
! are calculated based on the input dimensions and geometry
! author: Satish Karra
! date: 07/11/12
!
! ************************************************************************** !
subroutine SecondaryContinuumType(sec_continuum,nmat,aream, &
            volm,dm1,dm2,aperture,epsilon,log_spacing,outer_spacing, &
            interfacial_area,option)

  use Option_module

  implicit none
  
  type(sec_continuum_type) :: sec_continuum

  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: igeom, nmat, m
  PetscReal :: aream(nmat), volm(nmat), dm1(nmat), dm2(nmat)
  PetscReal :: dy, r0, r1, aream0, am0, vm0, interfacial_area
  PetscReal :: num_density, aperture, epsilon, fracture_spacing
  PetscReal :: outer_spacing, matrix_block_size
  PetscReal :: grid_spacing(nmat)
  PetscBool :: log_spacing

  PetscInt, save :: icall

  data icall/0/

  igeom = sec_continuum%itype
  option%nsec_cells = nmat
    
  select case (igeom)      
    case(0) ! 1D Slab
    
      dy = sec_continuum%slab%length/nmat
      aream0 = sec_continuum%slab%area
      do m = 1, nmat
        volm(m) = dy*aream0
      enddo
      am0 = aream0
      vm0 = nmat*dy*aream0
      interfacial_area = am0/vm0
     
       do m = 1, nmat
        aream(m) = aream0
        dm1(m) = 0.5d0*dy
        dm2(m) = 0.5d0*dy
      enddo

      if (icall == 0 .and. OptionPrintToFile(option)) then
        icall = 1
        string = 'DCDM Multiple Continuum Model'
        write(option%fid_out,'(/,2x,a,/)') trim(string)
        string = 'Slab'
        write(option%fid_out,'(2x,a,/)') trim(string)
        num_density = (1.d0-epsilon)/vm0
        write(option%fid_out,'(2x,"number density: ",11x,1pe12.4," m^(-3)")') num_density
        write(option%fid_out,'(2x,"matrix block size: ",8x,1pe12.4," m")') sec_continuum%slab%length
        write(option%fid_out,'(2x,"epsilon: ",18x,1pe12.4)') epsilon
        write(option%fid_out,'(2x,"specific interfacial area: ",1pe12.4," m^(-1)")') interfacial_area
        do m = 1, nmat
          if (m == 1) write(option%fid_out,'(/,2x,"node matrix volume fraction")') 
          write(option%fid_out,'(2x,i3,3x,1pe12.4)') m,volm(m)/vm0*(1.d0 - epsilon)
        enddo
!       aperture = r0*(1.d0/(1.d0-epsilon)**(1.d0/3.d0)-1.d0)
!       write(option%fid_out,'(2x,"aperture: ",17x,1pe12.4," m")') aperture
      endif
          
    case(1) ! nested cubes

      if (sec_continuum%nested_cube%fracture_spacing > 0.d0) then

        fracture_spacing = sec_continuum%nested_cube%fracture_spacing
!        override epsilon if aperture defined
        if (aperture > 0.d0) then
          r0 = fracture_spacing - aperture
          epsilon = 1.d0 - (1.d0 + aperture/r0)**(-3.0)
        else if (epsilon > 0.d0) then
          r0 = fracture_spacing*(1.d0-epsilon)**(1.d0/3.d0)
          aperture = r0*((1.d0-epsilon)**(-1.d0/3.d0)-1.d0)
        endif
                                            
      else if (sec_continuum%nested_cube%matrix_block_size > 0.d0) then

        r0 = sec_continuum%nested_cube%matrix_block_size

!        override epsilon if aperture defined
        if (aperture > 0.d0) then
          fracture_spacing = r0 + aperture
          epsilon = 1.d0 - (1.d0 + aperture/r0)**(-3.0)
        else if (epsilon > 0.d0) then
          fracture_spacing = r0*(1.d0-epsilon)**(-1.d0/3.d0)
          aperture = fracture_spacing - r0
        endif
      endif
      
      if (log_spacing) then 
        
        matrix_block_size = r0
        call SecondaryContinuumCalcLogSpacing(matrix_block_size,outer_spacing, &
                                              nmat,grid_spacing,option)
        
        
        r0 = 2*grid_spacing(1)                                                                                                                                  
        dm1(1) = 0.5*grid_spacing(1)
        dm2(1) = 0.5*grid_spacing(1)
        volm(1) = r0**3
        aream(1) = 6.d0*r0**2         
        do m = 2, nmat
          dm1(m) = 0.5*grid_spacing(m)
          dm2(m) = 0.5*grid_spacing(m)
          r1 = r0 + 2*(dm1(m) + dm2(m))
          volm(m) = r1**3 - r0**3
          aream(m) = 6.d0*r1**2
          r0 = r1
        enddo
        r0 = matrix_block_size
        am0 = 6.d0*r0**2
        vm0 = r0**3
        interfacial_area = am0/vm0

      else
        dy = r0/nmat/2.d0
     
        r0 = 2.d0*dy
        volm(1) = r0**3
        do m = 2, nmat
          r1 = r0 + 2.d0*dy
          volm(m) = r1**3 - r0**3
          r0 = r1
        enddo
      
        r0 = 2.d0*dy
        aream(1) = 6.d0*r0**2
        dm1(1) = 0.5d0*dy
        dm2(1) = 0.5d0*dy
        do m = 2, nmat
          dm1(m) = 0.5d0*dy
          dm2(m) = 0.5d0*dy
          r0 = r0 + 2.d0*dy
          aream(m) = 6.d0*r0**2
        enddo
        r0 = real(2*nmat)*dy
        am0 = 6.d0*r0**2
        vm0 = r0**3
        interfacial_area = am0/vm0
      endif

      if (icall == 0 .and. OptionPrintToFile(option)) then
        icall = 1
        string = 'DCDM Multiple Continuum Model'
        write(option%fid_out,'(/,2x,a,/)') trim(string)
        string = 'Nested Cubes'
        write(option%fid_out,'(2x,a,/)') trim(string)
        num_density = (1.d0-epsilon)/vm0
        write(option%fid_out,'(2x,"number density: ",11x,1pe12.4," m^(-3)")') num_density
        write(option%fid_out,'(2x,"matrix block size: ",8x,1pe12.4," m")') r0
        write(option%fid_out,'(2x,"epsilon: ",18x,1pe12.4)') epsilon
        write(option%fid_out,'(2x,"specific interfacial area: ",1pe12.4," m^(-1)")') interfacial_area
        write(option%fid_out,'(2x,"fracture aperture: ",8x,1pe12.4," m")') aperture
        write(option%fid_out,'(2x,"fracture spacing: ",9x,1pe12.4," m")') fracture_spacing
        do m = 1, nmat
          if (m == 1) write(option%fid_out,'(/,2x,"node matrix volume fraction")') 
          write(option%fid_out,'(2x,i3,3x,1pe12.4)') m,volm(m)/vm0*(1.d0 - epsilon)
        enddo
      endif

    case(2) ! nested spheres
    
      dy = sec_continuum%nested_sphere%radius/nmat
      r0 = dy

      volm(1) = 4.d0/3.d0*pi*r0**3
      do m = 2, nmat
        r1 = r0 + dy
        volm(m) = 4.d0/3.d0*pi*(r1**3 - r0**3)
        r0 = r1
      enddo
      
      r0 = dy
      aream(1) = 4.d0*pi*r0**2
      dm1(1) = 0.5d0*dy
      dm2(1) = 0.5d0*dy
      do m = 2, nmat
        r0 = r0 + dy
        dm1(m) = 0.5d0*dy
        dm2(m) = 0.5d0*dy
        aream(m) = 4.d0*pi*r0**2
      enddo
      r0 = 0.5d0*real(2*nmat)*dy
      am0 = 4.d0*pi*r0**2
      vm0 = am0*r0/3.d0
      interfacial_area = am0/vm0

      if (icall == 0 .and. OptionPrintToFile(option)) then
        icall = 1
        string = 'DCDM Multiple Continuum Model'
        write(option%fid_out,'(/,2x,a,/)') trim(string)
        string = 'Nested Spheres'
        write(option%fid_out,'(2x,a,/)') trim(string)
        num_density = (1.d0-epsilon)/vm0
        write(option%fid_out,'(2x,"number density: ",11x,1pe12.4," m^(-3)")') num_density
        write(option%fid_out,'(2x,"sphere radius: ",8x,1pe12.4," m")') sec_continuum%nested_sphere%radius
        write(option%fid_out,'(2x,"epsilon: ",18x,1pe12.4)') epsilon
        write(option%fid_out,'(2x,"specific interfacial area: ",1pe12.4," m^(-1)")') interfacial_area
        do m = 1, nmat
          if (m == 1) write(option%fid_out,'(/,2x,"node matrix volume fraction")') 
          write(option%fid_out,'(2x,i3,3x,1pe12.4)') m,volm(m)/vm0*(1.d0 - epsilon)
        enddo

!       aperture = r0*(1.d0/(1.d0-epsilon)**(1.d0/3.d0)-1.d0)
!       write(option%fid_out,'(2x,"aperture: ",17x,1pe12.4," m")') aperture
      endif
                        
  end select
  
end subroutine SecondaryContinuumType

! ************************************************************************** !
!
! SecondaryContinuumSetProperties: The type, dimensions of the secondary
! continuum are set 
! author: Satish Karra
! date: 07/17/12
!
! ************************************************************************** !

subroutine SecondaryContinuumSetProperties(sec_continuum, &
                                           sec_continuum_name, & 
                                           sec_continuum_length, &
                                           sec_continuum_matrix_block_size, &
                                           sec_continuum_fracture_spacing, &
                                           sec_continuum_radius, &
                                           sec_continuum_area, &
                                           option)
                                    
  use Option_module
  use String_module
  
  implicit none
  
  type(sec_continuum_type) :: sec_continuum
  type(option_type) :: option
  PetscReal :: sec_continuum_matrix_block_size
  PetscReal :: sec_continuum_fracture_spacing
  PetscReal :: sec_continuum_length
  PetscReal :: sec_continuum_area
  PetscReal :: sec_continuum_radius
  character(len=MAXWORDLENGTH) :: sec_continuum_name

  call StringToUpper(sec_continuum_name)
  
  select case(trim(sec_continuum_name))
    case("SLAB")
      sec_continuum%itype = 0
      sec_continuum%slab%length = sec_continuum_length
      if(sec_continuum_area == 0.d0) then
        option%io_buffer = 'Keyword "AREA" not specified for SLAB type ' // &
                           'under SECONDARY_CONTINUUM'
        call printErrMsg(option)
      endif
      sec_continuum%slab%area = sec_continuum_area
    case("NESTED_CUBES")
      sec_continuum%itype = 1
      sec_continuum%nested_cube%matrix_block_size = sec_continuum_matrix_block_size
      sec_continuum%nested_cube%fracture_spacing = sec_continuum_fracture_spacing
    case("NESTED_SPHERES")
      sec_continuum%itype = 2
      sec_continuum%nested_sphere%radius = sec_continuum_radius
    case default
      option%io_buffer = 'Keyword "' // trim(sec_continuum_name) // '" not ' // &
                         'recognized in SecondaryContinuumSetProperties()'
      call printErrMsg(option)  
  end select
      
end subroutine SecondaryContinuumSetProperties  

! ************************************************************************** !
!
! SecondaryContinuumCalcLogSpacing: Given the matrix block size and the 
! grid spacing of the outer mode secondary continuum cell, a geometric
! series is assumed and the grid spacing of the rest of the cells is 
! calculated
! author: Satish Karra, LANL
! date: 07/17/12
!
! ************************************************************************** !

subroutine SecondaryContinuumCalcLogSpacing(matrix_size,outer_grid_size, &
                                            sec_num_cells,grid_spacing,option)
                                              
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: matrix_size, outer_grid_size
  PetscInt :: sec_num_cells
  PetscReal :: grid_spacing(sec_num_cells)
  PetscReal :: delta, delta_new, inner_grid_size
  PetscReal :: F, dF
  PetscReal, parameter :: tol = 1.d-12
  PetscInt, parameter :: maxit = 50
  PetscInt :: i 
  
  
  if (mod(sec_num_cells,2) /= 0) then
     option%io_buffer = 'NUM_CELLS under SECONDARY_CONTINUUM has to be' // &
                        ' even for logarithmic grid spacing'
      call printErrMsg(option)
  endif
  
  delta = 0.85d0
  
  do i = 1, maxit
    F = (1.d0 - delta)/(1.d0 - delta**sec_num_cells)*delta**(sec_num_cells-1) - &
        2.d0*outer_grid_size/matrix_size
    dF = (1.d0 + sec_num_cells*(delta - 1.d0) - delta**sec_num_cells)/ &
         (delta**sec_num_cells - 1.d0)**2*delta**(sec_num_cells - 2) 
    delta_new = delta + F/dF
    if ((abs(F) < tol)) exit
    delta = delta_new
    if (delta < 0.d0) delta = 0.5d0
    if (delta > 1.d0) delta = 0.9d0
  enddo

  inner_grid_size = outer_grid_size/delta**(sec_num_cells - 1)
  
  do i = 1, sec_num_cells
    grid_spacing(i) = inner_grid_size*delta**(i-1)
  enddo
    
end subroutine SecondaryContinuumCalcLogSpacing

end module Secondary_Continuum_module
            