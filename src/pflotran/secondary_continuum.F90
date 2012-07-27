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
    PetscReal :: length                       ! input - side of cube
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
    PetscBool :: sec_temp_update              ! flag to check if the temp is updated
    PetscInt :: ncells                        ! number of secondary grid cells
    PetscReal :: epsilon                      ! vol. frac. of primary continuum
    type(sec_continuum_type) :: sec_continuum
    PetscReal, pointer :: sec_temp(:)         ! array of temp. at secondary grid cells
    PetscReal, pointer :: area(:)             ! surface area
    PetscReal, pointer :: vol(:)              ! volume     face      node       face
    PetscReal, pointer :: dm_plus(:)          ! see fig.    |----------o----------|
    PetscReal, pointer :: dm_minus(:)         ! see fig.      <dm_minus> <dm_plus>
    PetscReal :: interfacial_area             ! interfacial area between prim. and sec. per unit volume of prim.+sec.
  end type sec_heat_type  

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
                                  volm,dm1,dm2,epsilon,interfacial_area,option)
  use option_module
  implicit none
  
  type(sec_continuum_type) :: sec_continuum

  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: igeom, nmat, m
  PetscReal :: aream(nmat), volm(nmat), dm1(nmat), dm2(nmat)
  PetscReal :: dy, r0, r1, aream0, am0, vm0, interfacial_area
  PetscReal :: num_density, aperture, epsilon

  PetscInt, save :: icall

  data icall/0/

  igeom = sec_continuum%itype
    
  select case (igeom)      
    case(0) ! 1D
    
      dy = sec_continuum%slab%length/nmat
      aream0 = sec_continuum%slab%area
      do m = 1, nmat
        volm(m) = dy*aream0
      enddo
      am0 = 1.d0*aream0
      vm0 = nmat*dy*aream0
      interfacial_area = am0/vm0
     
       do m = 1, nmat
        aream(m) = aream0
        dm1(m) = 0.5d0*dy
        dm2(m) = 0.5d0*dy
      enddo
          
    case(1) ! nested cubes
    
      dy = sec_continuum%nested_cube%length/nmat/2.d0    
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

      if (icall == 0) then
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

        aperture = r0*(1.d0/(1.d0-epsilon)**(1.d0/3.d0)-1.d0)
        write(option%fid_out,'(2x,"aperture: ",17x,1pe12.4," m")') aperture
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
                                           sec_continuum_area, &
                                           option)
                                    
  use Option_module
  use String_module
  
  implicit none
  
  type(sec_continuum_type) :: sec_continuum
  type(option_type) :: option
  PetscReal :: sec_continuum_length
  PetscReal :: sec_continuum_area
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
      sec_continuum%nested_cube%length = sec_continuum_length
    case("NESTED_SPHERES")
      sec_continuum%itype = 2
      sec_continuum%nested_sphere%radius = sec_continuum_length
    case default
      option%io_buffer = 'Keyword "' // trim(sec_continuum_name) // '" not ' // &
                         'recognized in SecondaryContinuumSetProperties()'
      call printErrMsg(option)  
  end select
      
end subroutine SecondaryContinuumSetProperties  

end module Secondary_Continuum_module
            