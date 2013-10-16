module Geometry_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
 
  type, public :: point3d_type
    PetscReal :: x
    PetscReal :: y
    PetscReal :: z
  end type point3d_type

  type, public :: polygonal_volume_type
    type(point3d_type), pointer :: xy_coordinates(:)
    type(point3d_type), pointer :: xz_coordinates(:)
    type(point3d_type), pointer :: yz_coordinates(:)
  end type polygonal_volume_type

  interface GeometryPointInPolygon
    module procedure GeometryPointInPolygon1
    module procedure GeometryPointInPolygon2
  end interface GeometryPointInPolygon
  
  public :: GeometryCreatePolygonalVolume, &
            GeometryReadCoordinates, &
            GeometryPointInPolygonalVolume, &
            GeometryCopyCoordinates, &
            GeometryDestroyPolygonalVolume

contains

! ************************************************************************** !
!
! GeometryCreatePolygonalVolume: Creates a polygonal volume.  I.e. a volume
!                              defined by 3 polygons, on each plane of the
!                              principle coordinate system x,y,z
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
function GeometryCreatePolygonalVolume()

  implicit none
  
  type(polygonal_volume_type), pointer :: GeometryCreatePolygonalVolume
  
  type(polygonal_volume_type), pointer :: polygonal_volume
  
  allocate(polygonal_volume)
  nullify(polygonal_volume%xy_coordinates)
  nullify(polygonal_volume%xz_coordinates)
  nullify(polygonal_volume%yz_coordinates)
  
  GeometryCreatePolygonalVolume => polygonal_volume

end function GeometryCreatePolygonalVolume

! ************************************************************************** !
!
! GeometryReadCoordinates: Reads a list of coordinates
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
subroutine GeometryReadCoordinates(input,option,region_name,coordinates)

  use Input_Aux_module
  use Option_module

  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: region_name
  type(point3d_type), pointer :: coordinates(:)
  
  PetscInt :: icount
  PetscInt, parameter :: max_num_coordinates = 100
  type(point3d_type) :: temp_coordinates(max_num_coordinates)

  icount = 0
  do
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'REGION')
    if (InputCheckExit(input,option)) exit              
    icount = icount + 1
    if (icount > max_num_coordinates) then
      write(option%io_buffer, &
            '(''Number of coordinates in region '',a, &
          &'' exceeds limit of '',i3)') region_name, max_num_coordinates
      option%io_buffer = trim(option%io_buffer) // &
        ' Increase size of PetscInt, parameter :: max_num_coordinates ' // &
        ' in RegionReadCoordinates()'
      call printErrMsg(option)
    endif
    call InputReadDouble(input,option,temp_coordinates(icount)%x) 
    call InputErrorMsg(input,option,'x-coordinate','REGION')
    call InputReadDouble(input,option,temp_coordinates(icount)%y)
    call InputErrorMsg(input,option,'y-coordinate','REGION')
    call InputReadDouble(input,option,temp_coordinates(icount)%z)
    call InputErrorMsg(input,option,'z-coordinate','REGION')
  enddo
  allocate(coordinates(icount))
  do icount = 1, size(coordinates)
    coordinates(icount)%x = temp_coordinates(icount)%x
    coordinates(icount)%y = temp_coordinates(icount)%y
    coordinates(icount)%z = temp_coordinates(icount)%z
  enddo

end subroutine GeometryReadCoordinates

! ************************************************************************** !
!
! GeometryPointInPolygonalVolume: Determines whether a point in xyz space is 
!                               within a polygonal volume defined by polygons
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
function GeometryPointInPolygonalVolume(x,y,z,polygonal_volume,option)
 
  use Option_module

  implicit none
  
  PetscReal :: x, y, z
  type(polygonal_volume_type) :: polygonal_volume
  type(option_type) :: option
  
  PetscBool :: xy, xz, yz
  PetscBool :: GeometryPointInPolygonalVolume
    
  GeometryPointInPolygonalVolume = PETSC_FALSE
  
  ! if a polygon is not defined in a particular direction, it is assumed that
  ! the point is within the "infinite" polygon
  
  ! XY plane
  if (associated(polygonal_volume%xy_coordinates)) then
    if (.not.GeometryPointInPolygon(x,y,z,Z_DIRECTION, &
                                  polygonal_volume%xy_coordinates)) then
      return
    endif
  endif

  ! XZ plane
  if (associated(polygonal_volume%xz_coordinates)) then
    if (.not.GeometryPointInPolygon(x,y,z,Y_DIRECTION, &
                                  polygonal_volume%xz_coordinates)) then
      return
    endif
  endif

  ! YZ plane
  if (associated(polygonal_volume%yz_coordinates)) then
    if (.not.GeometryPointInPolygon(x,y,z,X_DIRECTION, &
                                  polygonal_volume%yz_coordinates)) then
      return
    endif
  endif
  
  ! if point is not within any of polygons above, the function will return
  ! prior to this point.
  GeometryPointInPolygonalVolume = PETSC_TRUE
  
end function GeometryPointInPolygonalVolume

! ************************************************************************** !
!
! GeometryPointInPolygon1: Determines whether a point in xyz space is within
!                        a 2d polygon based on coordinate object
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
function GeometryPointInPolygon1(x,y,z,axis,coordinates)

  implicit none
  
  PetscReal :: x, y, z
  PetscInt :: axis
  type(point3d_type) :: coordinates(:)
  
  PetscInt, parameter :: max_num_coordinates = 100
  PetscReal :: xx, yy
  PetscInt :: i, num_coordinates
  PetscReal :: xx_array(max_num_coordinates), yy_array(max_num_coordinates)
  
  PetscBool :: GeometryPointInPolygon1
  
  GeometryPointInPolygon1 = PETSC_FALSE
  
  num_coordinates = size(coordinates)
  select case(axis)
    case(Z_DIRECTION)
      xx = x
      yy = y
      do i = 1, num_coordinates
        xx_array(i) = coordinates(i)%x
        yy_array(i) = coordinates(i)%y
      enddo
    case(Y_DIRECTION)
      xx = x
      yy = z
      do i = 1, num_coordinates
        xx_array(i) = coordinates(i)%x
        yy_array(i) = coordinates(i)%z
      enddo
    case(X_DIRECTION)
      xx = y
      yy = z
      do i = 1, num_coordinates
        xx_array(i) = coordinates(i)%y
        yy_array(i) = coordinates(i)%z
      enddo
  end select
  if (num_coordinates == 2) then
    GeometryPointInPolygon1 = &
      GeometryPointInRectangle(xx,yy,xx_array,yy_array)
  else
    GeometryPointInPolygon1 = &
      GeometryPointInPolygon2(xx,yy,xx_array,yy_array,num_coordinates)
  endif
  
end function GeometryPointInPolygon1

! ************************************************************************** !
!
! GeometryPointInPolygon2: Determines whether a point in xy space is within
!                       a polygon
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
function GeometryPointInPolygon2(x,y,x_array,y_array,num_coordinates)

  implicit none
  
  PetscReal :: x, y
  PetscReal :: x_array(:), y_array(:)
  PetscInt :: num_coordinates
  
  PetscBool :: GeometryPointInPolygon2
  PetscInt :: i, j
  
  GeometryPointInPolygon2 = PETSC_FALSE

  j = 1
  do i = 1, num_coordinates
    j = i + 1
    if (i == num_coordinates) j = 1
    if ((y_array(i) < y .and. y_array(j) >= y) .or. &
        (y_array(j) < y .and. y_array(i) >= y)) then
      if ((x_array(i) + &
           (y-y_array(i))/(y_array(j)-y_array(i))*(x_array(j)-x_array(i))) &
          < x) then
        GeometryPointInPolygon2 = .not.GeometryPointInPolygon2
      endif
    endif    
  enddo

end function GeometryPointInPolygon2

! ************************************************************************** !
!
! GeometryPointInRectangle: Determines whether a point in xy space is within
!                           a box
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
function GeometryPointInRectangle(x,y,x_array,y_array)

  implicit none
  
  PetscReal :: x, y
  PetscReal :: x_array(:), y_array(:)

  PetscBool :: GeometryPointInRectangle
  
  PetscReal :: x1, y1
  PetscReal :: x2, y2
  
  GeometryPointInRectangle = PETSC_FALSE
  
  ! only using first 2 values in each array
  
  x1 = minval(x_array(1:2))
  y1 = minval(y_array(1:2))
  x2 = maxval(x_array(1:2))
  y2 = maxval(y_array(1:2))
  
  
  if (x >= x1 .and. x <= x2 .and. y >= y1 .and. y <= y2) &
    GeometryPointInRectangle = PETSC_TRUE

end function GeometryPointInRectangle

! ************************************************************************** !
!
! GeometryCopyCoordinates: Deallocates a polygonal volume object
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
subroutine GeometryCopyCoordinates(coordinates_in,coordinates_out)

  implicit none
  
  type(point3d_type) :: coordinates_in(:)
  type(point3d_type), pointer :: coordinates_out(:)
  
  PetscInt :: num_coordinates
  PetscInt :: i
  
  num_coordinates = size(coordinates_in)
  allocate(coordinates_out(num_coordinates))
  do i = 1, num_coordinates
   coordinates_out(i)%x = coordinates_in(i)%x 
   coordinates_out(i)%y = coordinates_in(i)%y
   coordinates_out(i)%z = coordinates_in(i)%z 
  enddo
 
end subroutine GeometryCopyCoordinates

! ************************************************************************** !
!
! GeometryDestroyPolygonalVolume: Deallocates a polygonal volume object
! author: Glenn Hammond
! date: 11/01/09
!
! ************************************************************************** !
subroutine GeometryDestroyPolygonalVolume(polygonal_volume)

  implicit none
  
  type(polygonal_volume_type), pointer :: polygonal_volume
  
  if (.not.associated(polygonal_volume)) return
  
  if (associated(polygonal_volume%xy_coordinates)) &
    deallocate(polygonal_volume%xy_coordinates)
  nullify(polygonal_volume%xy_coordinates)
  if (associated(polygonal_volume%xz_coordinates)) &
    deallocate(polygonal_volume%xz_coordinates)
  nullify(polygonal_volume%xz_coordinates)
  if (associated(polygonal_volume%yz_coordinates)) &
    deallocate(polygonal_volume%yz_coordinates)
  nullify(polygonal_volume%yz_coordinates)
  
  deallocate(polygonal_volume)
  nullify(polygonal_volume)
  
end subroutine GeometryDestroyPolygonalVolume

end module Geometry_module
