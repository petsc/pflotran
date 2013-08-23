module Unstructured_Cell_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  ! grid cell type
  PetscInt, parameter, public :: HEX_TYPE          = 1
  PetscInt, parameter, public :: TET_TYPE          = 2
  PetscInt, parameter, public :: WEDGE_TYPE        = 3
  PetscInt, parameter, public :: PYR_TYPE          = 4
  ! 2D cell types:
  PetscInt, parameter, public :: TRI_TYPE          = 5
  PetscInt, parameter, public :: QUAD_TYPE         = 6

  ! grid cell properties
  PetscInt, parameter, public :: LINE_FACE_TYPE    = 1
  PetscInt, parameter, public :: TRI_FACE_TYPE     = 2
  PetscInt, parameter, public :: QUAD_FACE_TYPE    = 3
  PetscInt, parameter, public :: MAX_VERT_PER_FACE = 4
  PetscInt, parameter, public :: MAX_FACE_PER_CELL = 6
  PetscInt, parameter, public :: MAX_FACE_PER_CELL_SURF = 4

  type, public :: point_type
    PetscInt :: id
    PetscReal :: x
    PetscReal :: y
    PetscReal :: z
  end type point_type
  
  type, public :: plane_type
    PetscReal :: A
    PetscReal :: B
    PetscReal :: C
    PetscReal :: D
  end type plane_type
  
  public :: UCellComputeCentroid, &
            UCellComputeVolume, &
            UCellComputeArea, &
            UCellComputePlane, &
            UCellGetPlaneIntercept, &
            UCellProjectPointOntoPlane, &
            UCellComputeDistanceFromPlane, &
            UCellGetLineIntercept, &
            UCellGetNVertices, &
            UCellGetNFaces, &
            UCellGetNFaceVertices, &
            UCellGetFaceType, &
            UCellGetFaceVertices, &
            UCellGetNFaceVertsandVerts, &
            UCellTypeToWord, &
            UCellFaceTypeToWord, &
            UCellQuality
            
contains

! ************************************************************************** !
!
! UCellComputeCentroid: Computes the centroid a grid cell
! author: Glenn Hammond
! date: 10/30/09
!
! ************************************************************************** !
function UCellComputeCentroid(cell_type,vertices,option)

  use Option_module

  implicit none
  
  PetscInt :: cell_type
  type(point_type) :: vertices(*)
  type(option_type) :: option
  
  PetscReal :: UCellComputeCentroid(3)
  PetscInt :: ivertex
  
  UCellComputeCentroid = 0.d0
  select case(cell_type)
    case(HEX_TYPE)
      ! need something more sophisticated, but for now, just use average
      do ivertex = 1, 8
        UCellComputeCentroid(1) = UCellComputeCentroid(1) + vertices(ivertex)%x
        UCellComputeCentroid(2) = UCellComputeCentroid(2) + vertices(ivertex)%y
        UCellComputeCentroid(3) = UCellComputeCentroid(3) + vertices(ivertex)%z
      enddo
      UCellComputeCentroid = UCellComputeCentroid / 8.d0
    case(WEDGE_TYPE)
      ! need something more sophisticated, but for now, just use average
      do ivertex = 1, 6
        UCellComputeCentroid(1) = UCellComputeCentroid(1) + vertices(ivertex)%x
        UCellComputeCentroid(2) = UCellComputeCentroid(2) + vertices(ivertex)%y
        UCellComputeCentroid(3) = UCellComputeCentroid(3) + vertices(ivertex)%z
      enddo
      UCellComputeCentroid = UCellComputeCentroid / 6.d0
    case(PYR_TYPE)
      ! need something more sophisticated, but for now, just use average
      do ivertex = 1, 5
        UCellComputeCentroid(1) = UCellComputeCentroid(1) + vertices(ivertex)%x
        UCellComputeCentroid(2) = UCellComputeCentroid(2) + vertices(ivertex)%y
        UCellComputeCentroid(3) = UCellComputeCentroid(3) + vertices(ivertex)%z
      enddo
      UCellComputeCentroid = UCellComputeCentroid / 5.d0
    case(TET_TYPE, QUAD_TYPE)
      do ivertex = 1, 4
        UCellComputeCentroid(1) = UCellComputeCentroid(1) + vertices(ivertex)%x
        UCellComputeCentroid(2) = UCellComputeCentroid(2) + vertices(ivertex)%y
        UCellComputeCentroid(3) = UCellComputeCentroid(3) + vertices(ivertex)%z
      enddo
      UCellComputeCentroid = UCellComputeCentroid / 4.d0
    case(TRI_TYPE)
      do ivertex = 1, 3
        UCellComputeCentroid(1) = UCellComputeCentroid(1) + vertices(ivertex)%x
        UCellComputeCentroid(2) = UCellComputeCentroid(2) + vertices(ivertex)%y
        UCellComputeCentroid(3) = UCellComputeCentroid(3) + vertices(ivertex)%z
      enddo
      UCellComputeCentroid = UCellComputeCentroid / 3.d0
    case default
      option%io_buffer = 'Cell type not recognized'
      call printErrMsg(option)
  end select

end function UCellComputeCentroid

! ************************************************************************** !
!
! UCellComputeVolume: Computes the volume a grid cell
! author: Glenn Hammond
! date: 11/06/09
!
! ************************************************************************** !
function UCellComputeVolume(cell_type,vertices,option)

  use Utility_module, only : DotProduct, CrossProduct
  use Option_module

  implicit none
  
  PetscInt :: cell_type
  type(point_type) :: vertices(*)
  type(option_type) :: option
  
  PetscReal :: UCellComputeVolume
  PetscReal :: v(3)
  PetscReal :: l1, l2, l3
  PetscReal :: n1(3), n2(3), v1(3), v2(3)
  PetscReal :: area1, area2, dz
  PetscReal :: vv(3,8)
  PetscInt :: i, j
  
  UCellComputeVolume = 0.d0
  select case(cell_type)
    case(HEX_TYPE)
      ! split into 5 tetrahedron
      UCellComputeVolume = &
        UCellComputeVolumeOfTetrahedron(vertices(1),vertices(2), &
                                        vertices(3),vertices(6))
      UCellComputeVolume = &
        UCellComputeVolume + &
        UCellComputeVolumeOfTetrahedron(vertices(1),vertices(6), &
                                        vertices(8),vertices(5))
      UCellComputeVolume = &
        UCellComputeVolume + &
        UCellComputeVolumeOfTetrahedron(vertices(1),vertices(6), &
                                        vertices(3),vertices(8))
      UCellComputeVolume = &
        UCellComputeVolume + &
        UCellComputeVolumeOfTetrahedron(vertices(8),vertices(6), &
                                        vertices(3),vertices(7))
      UCellComputeVolume = &
        UCellComputeVolume + &
        UCellComputeVolumeOfTetrahedron(vertices(1),vertices(3), &
                                        vertices(4),vertices(8))
    case(WEDGE_TYPE)
      ! split into 3 tetrahedrons: v1v2v3v4, v3v4v5v6, v2v3v4v5
      UCellComputeVolume = &
        UCellComputeVolumeOfTetrahedron(vertices(1),vertices(2),vertices(3), &
                                        vertices(4))
      UCellComputeVolume = &
        UCellComputeVolume + &
        UCellComputeVolumeOfTetrahedron(vertices(3),vertices(4),vertices(5), &
                                        vertices(6))
      UCellComputeVolume = &
        UCellComputeVolume + &
        UCellComputeVolumeOfTetrahedron(vertices(2),vertices(3),vertices(4), &
                                        vertices(5))
    case(PYR_TYPE)
      ! split pyramid into two tets and compute
      UCellComputeVolume = &
        UCellComputeVolumeOfTetrahedron(vertices(1),vertices(2), &
                                        vertices(3),vertices(5))
      UCellComputeVolume = &
        UCellComputeVolume + &
        UCellComputeVolumeOfTetrahedron(vertices(3),vertices(4), &
                                        vertices(1),vertices(5))
    case(TET_TYPE)
      UCellComputeVolume = &
        UCellComputeVolumeOfTetrahedron(vertices(1),vertices(2),vertices(3), &
                                        vertices(4))
    case(QUAD_TYPE, TRI_TYPE)
      option%io_buffer = 'Cell type is QUAD or TRI, thus one should call '// &
        'UCellComputeArea instead of UCellComputeVolume'
      call printErrMsg(option)
    case default
      option%io_buffer = 'Cell type not recognized'
      call printErrMsg(option)
  end select

end function UCellComputeVolume

! ************************************************************************** !
!
! UCellComputeArea: Computes the area a 2D grid cell
! author: Gautam Bisht
! date: 03/17/12
!
! ************************************************************************** !
function UCellComputeArea(cell_type,vertices,option)

  use Utility_module, only : DotProduct, CrossProduct
  use Option_module

  implicit none
  
  PetscInt :: cell_type
  type(point_type) :: vertices(*)
  type(option_type) :: option
  
  PetscReal :: UCellComputeArea
  PetscReal :: v(3)
  PetscReal :: l1, l2, l3
  PetscReal :: n1(3), n2(3), v1(3), v2(3)
  PetscReal :: area1, area2, dz
  PetscReal :: vv(3,8)
  PetscInt :: i, j
  
  UCellComputeArea = 0.d0
  select case(cell_type)
    case(QUAD_TYPE)
      v1(1) = vertices(3)%x-vertices(2)%x
      v1(2) = vertices(3)%y-vertices(2)%y
      v1(3) = vertices(3)%z-vertices(2)%z
      v2(1) = vertices(1)%x-vertices(2)%x
      v2(2) = vertices(1)%y-vertices(2)%y
      v2(3) = vertices(1)%z-vertices(2)%z
      n1 = CrossProduct(v1,v2)
      area1 = 0.5d0*sqrt(DotProduct(n1,n1))
      
      v1(1) = vertices(3)%x-vertices(4)%x
      v1(2) = vertices(3)%y-vertices(4)%y
      v1(3) = vertices(3)%z-vertices(4)%z
      v2(1) = vertices(1)%x-vertices(4)%x
      v2(2) = vertices(1)%y-vertices(4)%y
      v2(3) = vertices(1)%z-vertices(4)%z
      n2 = CrossProduct(v1,v2)
      area2 = 0.5d0*sqrt(DotProduct(n2,n2))
      
      UCellComputeArea = (area1 + area2)
    case(TRI_TYPE)
      v1(1) = vertices(3)%x-vertices(2)%x
      v1(2) = vertices(3)%y-vertices(2)%y
      v1(3) = vertices(3)%z-vertices(2)%z
      v2(1) = vertices(1)%x-vertices(2)%x
      v2(2) = vertices(1)%y-vertices(2)%y
      v2(3) = vertices(1)%z-vertices(2)%z
      n1 = CrossProduct(v1,v2)
      area1 = 0.5d0*sqrt(DotProduct(n1,n1))

      UCellComputeArea = area1
    case default
      option%io_buffer = 'Cell type not recognized'
      call printErrMsg(option)
  end select

end function UCellComputeArea

! ************************************************************************** !
!
! UCellComputeVolumeOfTetrahedron: Computes the voluem of a tetrahedron
!                                  given four points
! author: Glenn Hammond
! date: 12/06/11
!
! ************************************************************************** !
function UCellComputeVolumeOfTetrahedron(point1,point2,point3,point4)

  use Utility_module, only : DotProduct, CrossProduct

  implicit none
  
  type(point_type) :: point1, point2, point3, point4
  
  PetscReal :: vv(3,4)
  PetscReal :: vv1_minus_vv4(3)
  PetscReal :: vv2_minus_vv4(3)
  PetscReal :: vv3_minus_vv4(3)
  PetscReal :: cross_2_minus_4_X_3_minus_4(3)
  PetscReal :: UCellComputeVolumeOfTetrahedron
  PetscInt :: i

  vv(1,1) = point1%x
  vv(2,1) = point1%y
  vv(3,1) = point1%z
  vv(1,2) = point2%x
  vv(2,2) = point2%y
  vv(3,2) = point2%z
  vv(1,3) = point3%x
  vv(2,3) = point3%y
  vv(3,3) = point3%z
  vv(1,4) = point4%x
  vv(2,4) = point4%y
  vv(3,4) = point4%z

  ! V = |(a-d).((b-d)x(c-d))| / 6

  !geh: Intel Visual Fortran creates temporary arrays and reports warnings
  !     to the screen.  Therefore, I will use temporary variables below

  !UCellComputeVolumeOfTetrahedron = dabs(DotProduct(vv(:,1)-vv(:,4), &
  !                                       CrossProduct(vv(:,2)-vv(:,4), &
  !                                                    vv(:,3)-vv(:,4)))) / &
  !                                  6.d0

  vv1_minus_vv4 = vv(:,1)-vv(:,4)
  vv2_minus_vv4 = vv(:,2)-vv(:,4)
  vv3_minus_vv4 = vv(:,3)-vv(:,4)
  cross_2_minus_4_X_3_minus_4 = CrossProduct(vv2_minus_vv4, &
                                             vv3_minus_vv4)
  UCellComputeVolumeOfTetrahedron = dabs(DotProduct(vv1_minus_vv4, &
                                         cross_2_minus_4_X_3_minus_4)) / &
                                    6.d0

end function UCellComputeVolumeOfTetrahedron

! ************************************************************************** !
!
! UCellComputePlane: Computes the plane intersected by 3 points
! author: Glenn Hammond
! date: 10/30/09
!
! ************************************************************************** !
subroutine UCellComputePlane(plane,point1,point2,point3)

  implicit none
  
  type(plane_type) :: plane
  type(point_type) :: point1, point2, point3
  
  PetscReal :: x1,y1,z1
  PetscReal :: x2,y2,z2
  PetscReal :: x3,y3,z3
  x1 = point1%x
  y1 = point1%y
  z1 = point1%z
  x2 = point2%x
  y2 = point2%y
  z2 = point2%z
  x3 = point3%x
  y3 = point3%y
  z3 = point3%z
  
  ! this grabbed from python script
  plane%A = y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2)
  plane%B = z1*(x2-x3)+z2*(x3-x1)+z3*(x1-x2)
  plane%C = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)
  plane%D = -1.*(x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1))

end subroutine UCellComputePlane

! ************************************************************************** !
!
! UCellProjectPointOntoPlane: Computes the intercept of a point with a plane
! author: Glenn Hammond
! date: 11/22/11
!
! ************************************************************************** !
subroutine UCellProjectPointOntoPlane(plane,point,intercept)

  implicit none
  
  type(plane_type) :: plane
  type(point_type) :: point
  type(point_type) :: intercept

  PetscReal :: scalar
  
  ! plane equation:
  !   A*x + B*y + C*z + D = 0

  scalar = (plane%A*point%x + plane%B*point%y + plane%C*point%z + plane%D) / &
           (plane%A*plane%A + plane%B*plane%B + plane%C*plane%C)
  
  intercept%x = point%x - plane%A * scalar
  intercept%y = point%y - plane%B * scalar
  intercept%z = point%z - plane%C * scalar

end subroutine UCellProjectPointOntoPlane

! ************************************************************************** !
!
! UCellGetPlaneIntercept: Computes the intercept of a line with a plane
! author: Glenn Hammond
! date: 10/30/09
!
! ************************************************************************** !
subroutine UCellGetPlaneIntercept(plane,point1,point2,intercept)

  implicit none
  
  type(plane_type) :: plane
  type(point_type) :: point1, point2
  type(point_type) :: intercept

  PetscReal :: x1,y1,z1
  PetscReal :: x2,y2,z2
  PetscReal :: u
    
  x1 = point1%x
  y1 = point1%y
  z1 = point1%z
  x2 = point2%x
  y2 = point2%y
  z2 = point2%z
 
 
  u = (plane%A*x1 + plane%B*y1 + plane%C*z1 + plane%D) / &
      (plane%A*(x1-x2) + plane%B*(y1-y2) + plane%C*(z1-z2))

  intercept%x = point1%x + u*(point2%x-point1%x)
  intercept%y = point1%y + u*(point2%y-point1%y)
  intercept%z = point1%z + u*(point2%z-point1%z)

end subroutine UCellGetPlaneIntercept

! ************************************************************************** !
!
! UCellComputeDistanceFromPlane: Computes the intercept of a line with a plane
! author: Glenn Hammond
! date: 10/24/11
!
! ************************************************************************** !
function UCellComputeDistanceFromPlane(plane,point)

  implicit none
  
  type(plane_type) :: plane
  type(point_type) :: point
  PetscReal :: UCellComputeDistanceFromPlane

  UCellComputeDistanceFromPlane = &
    (plane%A*point%x + plane%B*point%y + plane%C*point%z + plane%D) / &
    sqrt(plane%A*plane%A+plane%B*plane%B+plane%C*plane%C)
  
end function UCellComputeDistanceFromPlane

! ************************************************************************** !
!
! UCellGetNVertices: Returns number of vertices in a cell
! author: Glenn Hammond
! date: 10/24/11
!
! ************************************************************************** !
function UCellGetNVertices(cell_type,option)

  use Option_module
  implicit none
  
  PetscInt :: cell_type
  PetscInt :: UCellGetNVertices
  type(option_type) :: option
  
  select case(cell_type)
    case(HEX_TYPE)
      UCellGetNVertices = 8
    case(WEDGE_TYPE)
      UCellGetNVertices = 6
    case(PYR_TYPE)
      UCellGetNVertices = 5
    case(TET_TYPE, QUAD_TYPE)
      UCellGetNVertices = 4
    case(TRI_TYPE)
      UCellGetNVertices = 3
    case default
      option%io_buffer = 'Cell type not recognized'
      call printErrMsg(option)
  end select  
  
end function UCellGetNVertices

! ************************************************************************** !
!
! UCellGetNEdges: Returns number of edges in a cell
! author: Glenn Hammond
! date: 01/17/12
!
! ************************************************************************** !
function UCellGetNEdges(cell_type)

  implicit none
  
  PetscInt :: cell_type
  PetscInt :: UCellGetNEdges
  
  select case(cell_type)
    case(HEX_TYPE)
      UCellGetNEdges = 12
    case(WEDGE_TYPE)
      UCellGetNEdges = 9
    case(PYR_TYPE)
      UCellGetNEdges = 8
    case(TET_TYPE)
      UCellGetNEdges = 6
    case(QUAD_TYPE)
      UCellGetNEdges = 4
    case(TRI_TYPE)
      UCellGetNEdges = 3
  end select  
  
end function UCellGetNEdges

! ************************************************************************** !
!
! UCellGetNFaces: Returns number of faces in a cell
! author: Glenn Hammond
! date: 10/24/11
!
! ************************************************************************** !
function UCellGetNFaces(cell_type,option)

  use Option_module
  implicit none
  
  PetscInt :: cell_type
  PetscInt :: UCellGetNFaces
  type(option_type) :: option
  
  select case(cell_type)
    case(HEX_TYPE)
      UCellGetNFaces = 6
    case(WEDGE_TYPE,PYR_TYPE)
      UCellGetNFaces = 5
    case(TET_TYPE)
      UCellGetNFaces = 4
    case(QUAD_TYPE)
      UCellGetNFaces = 4
    case(TRI_TYPE)
      UCellGetNFaces = 3
    case default
      option%io_buffer = 'Cell type not recognized'
      call printErrMsg(option)
  end select  
  
end function UCellGetNFaces

! ************************************************************************** !
!
! UCellGetNFaceVertices: Returns number of vertices in a cell face
! author: Glenn Hammond
! date: 10/24/11
!
! ************************************************************************** !
function UCellGetNFaceVertices(cell_type,iface,option)

  use Option_module
  implicit none
  
  PetscInt :: cell_type
  PetscInt :: iface
  PetscInt :: UCellGetNFaceVertices
  type(option_type) :: option
  
  select case(cell_type)
    case(HEX_TYPE)
      UCellGetNFaceVertices = 4
    case(WEDGE_TYPE)
      if (iface > 3) then
        UCellGetNFaceVertices = 3
      else 
        UCellGetNFaceVertices = 4
      endif
    case(PYR_TYPE)
      if (iface > 4) then
        UCellGetNFaceVertices = 4
      else 
        UCellGetNFaceVertices = 3
      endif
    case(TET_TYPE)
      UCellGetNFaceVertices = 3
    case(QUAD_TYPE)
      UCellGetNFaceVertices = 2
    case(TRI_TYPE)
      UCellGetNFaceVertices = 2
    case default
      option%io_buffer = 'Cell type not recognized'
      call printErrMsg(option)
  end select
      
end function UCellGetNFaceVertices

! ************************************************************************** !
!
! UCellGetFaceType: Returns type of cell face
! author: Glenn Hammond
! date: 10/24/11
!
! ************************************************************************** !
function UCellGetFaceType(cell_type,iface,option)

  use Option_module
  implicit none
  
  PetscInt :: cell_type
  PetscInt :: iface
  PetscInt :: UCellGetFaceType
  type(option_type) :: option
  
  select case(cell_type)
    case(HEX_TYPE)
      UCellGetFaceType = QUAD_FACE_TYPE
    case(WEDGE_TYPE)
      if (iface > 3) then
        UCellGetFaceType = TRI_FACE_TYPE
      else 
        UCellGetFaceType = QUAD_FACE_TYPE
      endif
    case(PYR_TYPE)
      if (iface > 4) then
        UCellGetFaceType = QUAD_FACE_TYPE
      else 
        UCellGetFaceType = TRI_FACE_TYPE
      endif
    case(TET_TYPE)
      UCellGetFaceType = TRI_FACE_TYPE
    case(QUAD_TYPE)
      UCellGetFaceType = LINE_FACE_TYPE
    case(TRI_TYPE)
      UCellGetFaceType = LINE_FACE_TYPE
    case default
      option%io_buffer = 'Cell type not recognized'
      call printErrMsg(option)
  end select
  
end function UCellGetFaceType

! ************************************************************************** !
!
! UCellTypeToWord: Returns type of cell as a string
! author: Glenn Hammond
! date: 12/09/11
!
! ************************************************************************** !
function UCellTypeToWord(cell_type,option)

  use Option_module
  implicit none
  
  PetscInt :: cell_type
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: UCellTypeToWord
  
  select case(cell_type)
    case(HEX_TYPE)
      UCellTypeToWord = 'hexahedron'
    case(WEDGE_TYPE)
      UCellTypeToWord = 'wedge'
    case(PYR_TYPE)
      UCellTypeToWord = 'pyramid'
    case(TET_TYPE)
      UCellTypeToWord = 'tetrahedron'
    case default
      option%io_buffer = 'Cell type not recognized'
      call printErrMsg(option)
  end select
  
end function UCellTypeToWord

! ************************************************************************** !
!
! UCellFaceTypeToWord: Returns type of cell face as a string
! author: Glenn Hammond
! date: 12/09/11
!
! ************************************************************************** !
function UCellFaceTypeToWord(face_type,option)

  use Option_module
  implicit none
  
  PetscInt :: face_type
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: UCellFaceTypeToWord
  
  select case(face_type)
    case(TRI_FACE_TYPE)
      UCellFaceTypeToWord = 'triangle'
    case(QUAD_FACE_TYPE)
      UCellFaceTypeToWord = 'quadrilateral'
    case default
      option%io_buffer = 'Face type not recognized'
      call printErrMsg(option)
  end select
  
end function UCellFaceTypeToWord

! ************************************************************************** !
!
! UCellGetNFaceVertsandVerts: returns the numbber of vertices for a face and
!                             the vertices
! author: Glenn Hammond
! date: 12/06/11
!
! ************************************************************************** !
subroutine UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                      vertex_ids)
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscInt :: cell_type
  PetscInt :: iface
  PetscInt :: nvertices
  PetscInt :: vertex_ids(*)
  
  nvertices = UCellGetNFaceVertices(cell_type,iface,option)
  call UCellGetFaceVertices(option,cell_type,iface,vertex_ids)

end subroutine UCellGetNFaceVertsandVerts

! ************************************************************************** !
!
! UCellGetFaceVertices: returns vertex ids of a face
! author: Glenn Hammond
! date: 11/24/11
!
! ************************************************************************** !
subroutine UCellGetFaceVertices(option,cell_type,iface,vertex_ids)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscInt :: cell_type
  PetscInt :: iface
  PetscInt :: vertex_ids(*)
  
  select case(cell_type)
    case(HEX_TYPE)
      select case(iface)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
          vertex_ids(3) = 6
          vertex_ids(4) = 5
        case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
          vertex_ids(3) = 7
          vertex_ids(4) = 6
        case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 4
          vertex_ids(3) = 8
          vertex_ids(4) = 7
        case(4)
          vertex_ids(1) = 4
          vertex_ids(2) = 1
          vertex_ids(3) = 5
          vertex_ids(4) = 8
        case(5)
          vertex_ids(1) = 1
          vertex_ids(2) = 4
          vertex_ids(3) = 3
          vertex_ids(4) = 2
        case(6)
          vertex_ids(1) = 5
          vertex_ids(2) = 6
          vertex_ids(3) = 7
          vertex_ids(4) = 8
      end select
    case(WEDGE_TYPE)
      select case(iface)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
          vertex_ids(3) = 5
          vertex_ids(4) = 4
       case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
          vertex_ids(3) = 6
          vertex_ids(4) = 5
       case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 1
          vertex_ids(3) = 4
          vertex_ids(4) = 6
       case(4)
          vertex_ids(1) = 1
          vertex_ids(2) = 3
          vertex_ids(3) = 2
       case(5)
          vertex_ids(1) = 4
          vertex_ids(2) = 5
          vertex_ids(3) = 6
       case default
          option%io_buffer='Cell WEDGE_TYPE has only 5 faces'
          call printErrMsg(option)
       end select
    case(PYR_TYPE)
      select case(iface)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
          vertex_ids(3) = 5
       case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
          vertex_ids(3) = 5
       case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 4
          vertex_ids(3) = 5
       case(4)
          vertex_ids(1) = 4
          vertex_ids(2) = 1
          vertex_ids(3) = 5
       case(5)
          vertex_ids(1) = 1
          vertex_ids(2) = 4
          vertex_ids(3) = 3
          vertex_ids(4) = 2
       case default
          option%io_buffer='Cell PYR_TYPE has only 5 faces'
          call printErrMsg(option)
       end select
    case(TET_TYPE)
      select case(iface)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
          vertex_ids(3) = 4
        case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
          vertex_ids(3) = 4
        case(3)
          vertex_ids(1) = 1
          vertex_ids(2) = 4
          vertex_ids(3) = 3
        case(4)
          vertex_ids(1) = 1
          vertex_ids(2) = 3
          vertex_ids(3) = 2
        case default
          option%io_buffer='Cell TET_TYPE has only 4 faces'
          call printErrMsg(option)
      end select       
    case(QUAD_TYPE)
      select case(iface)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
        case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
        case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 4
        case(4)
          vertex_ids(1) = 4
          vertex_ids(2) = 1
        case default
          option%io_buffer='Cell QUAD_TYPE has only 4 faces'
          call printErrMsg(option)
      end select
    case(TRI_TYPE)
      select case(iface)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
        case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
        case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 1
        case default
          option%io_buffer='Cell TRI_TYPE has only 3 faces'
          call printErrMsg(option)
      end select
    case default
      option%io_buffer = 'Cell type not recognized'
      call printErrMsg(option)
  end select

end subroutine UCellGetFaceVertices

! ************************************************************************** !
!
! UCellGetEdgeVertices: returns vertex ids of an edge
! author: Glenn Hammond
! date: 01/17/12
!
! ************************************************************************** !
subroutine UCellGetEdgeVertices(cell_type,iedge,vertex_ids)

  implicit none
  
  PetscInt :: cell_type
  PetscInt :: iedge
  PetscInt :: vertex_ids(2)
  
  vertex_ids = -999
  
  select case(cell_type)
    case(HEX_TYPE)
      select case(iedge)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
        case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
        case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 4
        case(4)
          vertex_ids(1) = 4
          vertex_ids(2) = 1
        case(5)
          vertex_ids(1) = 5
          vertex_ids(2) = 6
        case(6)
          vertex_ids(1) = 6
          vertex_ids(2) = 7
        case(7)
          vertex_ids(1) = 7
          vertex_ids(2) = 8
        case(8)
          vertex_ids(1) = 8
          vertex_ids(2) = 5
        case(9)
          vertex_ids(1) = 1
          vertex_ids(2) = 5
        case(10)
          vertex_ids(1) = 2
          vertex_ids(2) = 6
        case(11)
          vertex_ids(1) = 3
          vertex_ids(2) = 7
        case(12)
          vertex_ids(1) = 4
          vertex_ids(2) = 8
      end select
    case(WEDGE_TYPE)
      select case(iedge)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
        case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
        case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 1
        case(4)
          vertex_ids(1) = 4
          vertex_ids(2) = 5
        case(5)
          vertex_ids(1) = 5
          vertex_ids(2) = 6
        case(6)
          vertex_ids(1) = 6
          vertex_ids(2) = 4
        case(7)
          vertex_ids(1) = 1
          vertex_ids(2) = 4
        case(8)
          vertex_ids(1) = 2
          vertex_ids(2) = 5
        case(9)
          vertex_ids(1) = 3
          vertex_ids(2) = 6
       end select
    case(PYR_TYPE)
      select case(iedge)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
        case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
        case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 4
        case(4)
          vertex_ids(1) = 4
          vertex_ids(2) = 1
        case(5)
          vertex_ids(1) = 1
          vertex_ids(2) = 5
        case(6)
          vertex_ids(1) = 2
          vertex_ids(2) = 5
        case(7)
          vertex_ids(1) = 3
          vertex_ids(2) = 5
        case(8)
          vertex_ids(1) = 4
          vertex_ids(2) = 5
       end select
    case(TET_TYPE)
      select case(iedge)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
        case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
        case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 1
        case(4)
          vertex_ids(1) = 1
          vertex_ids(2) = 4
        case(5)
          vertex_ids(1) = 2
          vertex_ids(2) = 4
        case(6)
          vertex_ids(1) = 3
          vertex_ids(2) = 4
      end select       
    case(QUAD_TYPE)
      select case(iedge)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
        case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
        case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 4
        case(4)
          vertex_ids(1) = 4
          vertex_ids(2) = 1
      end select
    case(TRI_TYPE)
      select case(iedge)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
        case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
        case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 1
      end select
  end select

end subroutine UCellGetEdgeVertices

! ************************************************************************** !
!
! UCellGetEdgeLength: returns length of an edge
! author: Glenn Hammond
! date: 01/17/12
!
! ************************************************************************** !
function UCellGetEdgeLength(cell_type,iedge,vertices)

  implicit none
  
  PetscInt :: cell_type
  PetscInt :: iedge
  type(point_type) :: vertices(*)
  
  PetscReal :: UCellGetEdgeLength
  
  PetscReal :: len_x, len_y, len_z
  PetscInt :: vertex_ids(2)
  
  call UCellGetEdgeVertices(cell_type,iedge,vertex_ids)
  
  len_x = vertices(vertex_ids(1))%x-vertices(vertex_ids(2))%x
  len_y = vertices(vertex_ids(1))%y-vertices(vertex_ids(2))%y
  len_z = vertices(vertex_ids(1))%z-vertices(vertex_ids(2))%z
  UCellGetEdgeLength = sqrt(len_x*len_x+len_y*len_y+len_z*len_z)

end function UCellGetEdgeLength

! ************************************************************************** !
!
! UCellQuality: returns vertex ids of a face
! author: Glenn Hammond
! date: 01/17/12
!
! ************************************************************************** !
function UCellQuality(cell_type,vertices,option)

  use Option_module

  implicit none
  
  PetscInt :: cell_type
  type(point_type) :: vertices(*)
  type(option_type) :: option
  
  PetscReal :: UCellQuality
  
  PetscInt :: iedge
  PetscReal :: max_side, min_side
  PetscReal :: len

  max_side = -1.d20
  min_side = 1.d20
  do iedge = 1, UCellGetNVertices(cell_type,option)
    len = UCellGetEdgeLength(cell_type,iedge,vertices)
    if (len > max_side) max_side = len
    if (len < min_side) min_side = len
  enddo
  
  UCellQuality = max_side / min_side

end function UCellQuality

! ************************************************************************** !
!
! UCellGetLineIntercept: Computes the intercept of a point with a line
! author: Gautam Bisht
! date: 02/26/12
!
! ************************************************************************** !
subroutine UCellGetLineIntercept(line_start,line_end,point,intercept)

  implicit none
  type(point_type) :: line_start
  type(point_type) :: line_end
  type(point_type) :: point
  type(point_type) :: intercept

  PetscReal :: dx,dy,dz
  PetscReal :: u, line_mag
  
  dx = (line_end%x - line_start%x)
  dy = (line_end%y - line_start%y)
  dz = (line_end%z - line_start%z)
  
  line_mag = (dx*dx + dy*dy + dz*dz)
  
  u = ((point%x - line_start%x)*(line_end%x - line_start%x) + &
       (point%y - line_start%y)*(line_end%y - line_start%y) + &
       (point%z - line_start%z)*(line_end%z - line_start%z))/ line_mag

  intercept%x = line_start%x + u*(line_end%x - line_start%x)
  intercept%y = line_start%y + u*(line_end%y - line_start%y)
  intercept%z = line_start%z + u*(line_end%z - line_start%z)

end subroutine

end module Unstructured_Cell_module
