module Material_module
 
  implicit none

  private
 
  type, public :: material_type
    real*8, pointer :: permeability(:,:,:)
    real*8, pointer :: porosity(:,:)
    real*8, pointer :: tortuosity(:,:)

    real*8 :: scale

    ! hydraulic properties
    integer :: iregperm, iran_por=0, iread_perm=0, iread_geom=1
    real** :: ran_fac = -1.d0
    integer*4, pointer :: i1reg(:),i2reg(:),j1reg(:),j2reg(:),k1reg(:),k2reg(:)
    real*8, pointer :: por_reg(:),tor_reg(:),perm_reg(:,:)
    
    
      
    ! thermal properties
    integer, pointer :: ithrm_reg(:)
    real*8, pointer :: rock_density(:)
    real*8, pointer :: cpr(:)
    real*8, pointer :: ckdry(:)
    real*8, pointer :: ckwet(:)
    real*8, pointer :: tau(:)
    real*8, pointer :: cdiff(:)
    real*8, pointer :: cexp(:)
    real*8, pointer :: dencpr(:)

    ! capillary function properties
    integer, pointer :: icap_reg(:)
    integer, pointer:: icaptype(:)
    real*8, pointer :: swir(:)
    real*8, pointer :: lambda(:)
    real*8, pointer :: alpha(:)
    real*8, pointer :: pckrm(:)
    real*8, pointer :: pcwmax(:)
    real*8, pointer :: pcbetac(:)
    real*8, pointer :: pwrprm(:)
    real*8, pointer :: sir(:)
    real*8, pointer :: fracture_aperture(:)
    real*8, pointer :: matrix_block(:)
    
    integer, pointer :: imat(:)
    
  end type material_type
  
end module Material_module
