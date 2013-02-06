! added by S. Karra 07/11/12

module Secondary_Continuum_Aux_module

  use Reactive_Transport_Aux_module

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
    PetscInt :: ncells                         ! number of secondary grid cells
    PetscReal :: aperture                      ! fracture aperture
    PetscReal :: epsilon                       ! vol. frac. of primary continuum
    type(sec_continuum_type) :: sec_continuum
    type(reactive_transport_auxvar_type), pointer :: sec_rt_auxvar(:)  ! for each secondary grid cell
    PetscReal, pointer :: sec_mnrl_volfrac(:,:)  ! array of mineral vol fraction at secondary grid cells for each species  (naqcomp x ncells)
    PetscInt, pointer :: sec_zeta(:,:)           ! array of zetas at secondary grid cells for each species  (naqcomp x ncells)
    PetscReal, pointer :: area(:)              ! surface area
    PetscReal, pointer :: vol(:)               ! volume     face      node       face
    PetscReal, pointer :: dm_plus(:)           ! see fig.    |----------o----------|
    PetscReal, pointer :: dm_minus(:)          ! see fig.      <dm_minus> <dm_plus>
    PetscReal :: interfacial_area              ! interfacial area between prim. and sec. per unit volume of prim.+sec.
    PetscReal, pointer :: sec_mnrl_area(:)     ! secondary mineral surface area for each species (size = naqcomp)
    PetscBool :: log_spacing                   ! flag to check if log spacing is set
    PetscReal :: outer_spacing                 ! value of the outer most grid cell spacing
    PetscReal, pointer :: updated_conc(:,:)    ! This stores the secondary concentration update values from secondary NR iteration  (naqcomp x ncells)
    PetscReal, pointer :: sec_jac(:,:)         ! stores the secondary continuum jacobian value (naqcomp x ncells)
    PetscBool :: sec_jac_update                ! flag to check if secondary jacobian is updated
    PetscReal, pointer :: cxm(:,:,:)           ! stores the coeff of left diag in block triag system (ncomp x ncomp x ncells-1)
    PetscReal, pointer :: cxp(:,:,:)           ! stores the coeff of right diag in block triag system (ncomp x ncomp x ncells-1)
    PetscReal, pointer :: cdl(:,:,:)           ! stores the coeff of central diag in block triag system (ncomp x ncomp x ncells)
    PetscReal, pointer :: r(:)                 ! stores the solution of the forward solve
  end type sec_transport_type  
        

  type, public :: sc_heat_type
    type(sec_heat_type), pointer :: sec_heat_vars(:)
  end type sc_heat_type

  type, public :: sc_rt_type
    type(sec_transport_type), pointer :: sec_transport_vars(:)
  end type sc_rt_type

  public :: SecondaryAuxHeatCreate, SecondaryAuxHeatDestroy, &
            SecondaryAuxRTCreate, SecondaryAuxRTDestroy, &
            SecondaryRTAuxVarComputeMulti, &
            THCSecHeatAuxVarCompute, &
            MphaseSecHeatAuxVarCompute
            
contains
  
  
! ************************************************************************** !
!
! SecondaryAuxHeatCreate: Allocate and initialize secondary continuum heat
! auxiliary object
! author: Satish Karra
! date: 01/10/13
!
! ************************************************************************** !
function SecondaryAuxHeatCreate(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(sc_heat_type), pointer :: SecondaryAuxHeatCreate
  
  type(sc_heat_type), pointer :: aux

  allocate(aux) 
  nullify(aux%sec_heat_vars)
  
  SecondaryAuxHeatCreate => aux
  
end function SecondaryAuxHeatCreate  
  
! ************************************************************************** !
!
! SecondaryAuxHeatDestroy: Deallocates a secondary continuum heat
! auxiliary object
! author: Satish Karra
! date: 01/10/13
!
! ************************************************************************** !
subroutine SecondaryAuxHeatDestroy(aux)

  implicit none

  type(sc_heat_type), pointer :: aux
   
  if (.not.associated(aux)) return
  
  deallocate(aux)
  nullify(aux)  

end subroutine SecondaryAuxHeatDestroy


! ************************************************************************** !
!
! SecondaryAuxRTCreate: Allocate and initialize secondary continuum
! reactive transport auxiliary object
! author: Satish Karra
! date: 01/10/13
!
! ************************************************************************** !
function SecondaryAuxRTCreate(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(sc_rt_type), pointer :: SecondaryAuxRTCreate
  
  type(sc_rt_type), pointer :: aux

  allocate(aux) 
  nullify(aux%sec_transport_vars)
  
  SecondaryAuxRTCreate => aux
  
end function SecondaryAuxRTCreate  
  
! ************************************************************************** !
!
! SecondaryAuxRTDestroy: Deallocates a secondary continuum reactive 
! transport auxiliary object
! author: Satish Karra
! date: 01/10/13
!
! ************************************************************************** !
subroutine SecondaryAuxRTDestroy(aux)

  implicit none

  type(sc_rt_type), pointer :: aux
   
  if (.not.associated(aux)) return
  
  deallocate(aux)
  nullify(aux)  

end subroutine SecondaryAuxRTDestroy

! ************************************************************************** !
!
! SecondaryRTAuxVarComputeMulti: Updates the secondary continuum 
! concentrations at end of each time step for multicomponent system
! author: Satish Karra
! date: 2/1/13
!
! ************************************************************************** !
subroutine SecondaryRTAuxVarComputeMulti(sec_transport_vars,aux_var, &
                                         global_aux_var,reaction, &
                                         diffusion_coefficient,porosity, &
                                         option)
                               
                            
  use Option_module 
  use Global_Aux_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use blksolv_module
  use Utility_module
  

  implicit none
  
  type(sec_transport_type) :: sec_transport_vars
  type(reactive_transport_auxvar_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var
  type(reaction_type) :: reaction
  type(option_type) :: option
  PetscReal :: coeff_left(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells-1)
  PetscReal :: coeff_diag(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells)
  PetscReal :: coeff_right(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells-1)
  PetscReal :: rhs(sec_transport_vars%ncells*reaction%naqcomp)
  
  PetscReal :: conc_upd(reaction%naqcomp,sec_transport_vars%ncells) 
  PetscReal :: sec_mnrl_volfrac(reaction%naqcomp,sec_transport_vars%ncells)
  PetscInt :: sec_zeta(reaction%naqcomp,sec_transport_vars%ncells)

  PetscReal :: res_transport(reaction%naqcomp)
  PetscReal :: kin_mnrl_rate(reaction%naqcomp)
  PetscReal :: mnrl_area(reaction%naqcomp)
  PetscReal :: mnrl_molar_vol(reaction%naqcomp)
  PetscReal :: equil_conc(reaction%naqcomp)
  PetscReal :: conc_primary_node(reaction%naqcomp)
  PetscReal :: Im
  
  PetscReal :: area(sec_transport_vars%ncells)
  PetscReal :: vol(sec_transport_vars%ncells)
  PetscReal :: dm_plus(sec_transport_vars%ncells)
  PetscReal :: dm_minus(sec_transport_vars%ncells)

  
  PetscInt :: i, j, k, n
  PetscInt :: ngcells, ncomp
  PetscReal :: area_fm
  PetscReal :: diffusion_coefficient
  PetscReal :: porosity
  PetscReal, parameter :: rgas = 8.3144621d-3
  PetscReal :: arrhenius_factor
  PetscReal :: pordt, pordiff
  
  PetscInt :: pivot(reaction%naqcomp,sec_transport_vars%ncells)
  PetscInt :: indx(reaction%naqcomp)
  PetscInt :: d
  
  ngcells = sec_transport_vars%ncells
  area = sec_transport_vars%area
  vol = sec_transport_vars%vol          
  dm_plus = sec_transport_vars%dm_plus
  dm_minus = sec_transport_vars%dm_minus
  area_fm = sec_transport_vars%interfacial_area
  sec_zeta = sec_transport_vars%sec_zeta
  conc_upd = sec_transport_vars%updated_conc
  ncomp = reaction%naqcomp
  ! Note that sec_transport_vars%sec_conc units are in mol/kg
  ! Need to convert to mol/L since the units of conc. in the Thomas 
  ! algorithm are in mol/L 
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  
  conc_primary_node = aux_var%total(:,1)                             ! in mol/L 
  pordt = porosity/option%tran_dt
  pordiff = porosity*diffusion_coefficient            
              
  
  if (reaction%mineral%nkinmnrl > 0) then
    kin_mnrl_rate = reaction%mineral%kinmnrl_rate                    ! in mol/m^2/s
    do i = 1, reaction%mineral%nkinmnrl
      ! Arrhenius factor
      arrhenius_factor = 1.d0
      if (reaction%mineral%kinmnrl_activation_energy(i) > 0.d0) then
        arrhenius_factor = exp(reaction%mineral%kinmnrl_activation_energy(i)/ &
        rgas*(1.d0/(25.d0+273.15d0)-1.d0/(global_aux_var%temp(1)+273.15d0)))
      endif    
      kin_mnrl_rate(i) = kin_mnrl_rate(i)*arrhenius_factor
      equil_conc(i) = (10.d0)**(reaction%mineral%mnrl_logK(i))          ! in mol/kg --> Note!
    enddo
    equil_conc = equil_conc*global_aux_var%den_kg(1)*1.d-3           ! in mol/L
    mnrl_molar_vol = reaction%mineral%kinmnrl_molar_vol              ! in m^3
    sec_mnrl_volfrac = sec_transport_vars%sec_mnrl_volfrac           ! dimensionless
    mnrl_area = sec_transport_vars%sec_mnrl_area                     ! in 1/m
  endif
   
  ! Use the stored coefficient matrices from LU decomposition of the
  ! block triagonal sytem
  coeff_left = sec_transport_vars%cxm
  coeff_right = sec_transport_vars%cxp
  coeff_diag = sec_transport_vars%cdl
  rhs = sec_transport_vars%r
        
  call bl3dsolb(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot,1,rhs)
  
  do j = 1, ncomp
    do i = 1, ngcells
      n = j + (i - 1)*ncomp
      conc_upd(j,i) = rhs(n) + conc_upd(j,i)
    enddo
  enddo

  ! Mineral linear kinetics  
  do j = 1, ncomp
    do i = 1, ngcells
      Im = kin_mnrl_rate(j)*mnrl_area(j)*(conc_upd(j,i)/equil_conc(j) - 1.d0) ! in mol/m^3/s
      if (Im > 0.d0) then 
        sec_mnrl_volfrac(j,i) = sec_mnrl_volfrac(j,i) + option%tran_dt* &
                              mnrl_molar_vol(j)*Im
        sec_zeta(j,i) = 1
      else
        if (sec_mnrl_volfrac(j,i) > 0.d0) then
          sec_mnrl_volfrac(j,i) = sec_mnrl_volfrac(j,i) + option%tran_dt* &
                                mnrl_molar_vol(j)*Im
          sec_zeta(j,i) = 1
        else
          Im = 0.d0
          sec_zeta(j,i) = 0
        endif
      endif
      if (sec_mnrl_volfrac(j,i) < 0.d0) then
        sec_mnrl_volfrac(j,i) = 0.d0
        sec_zeta(j,i) = 0
      endif
    enddo
  enddo

  ! Convert the units of sec_conc from mol/L to mol/kg before passing to
  ! sec_transport_vars
  do j = 1, ncomp
    do i = 1, ngcells
      sec_transport_vars%sec_rt_auxvar(i)%pri_molal(j) = conc_upd(j,i)/ &
        global_aux_var%den_kg(1)/1.d-3
    enddo
  enddo
  
  sec_transport_vars%sec_mnrl_volfrac = sec_mnrl_volfrac
  sec_transport_vars%sec_zeta = sec_zeta  

end subroutine SecondaryRTAuxVarComputeMulti

! ************************************************************************** !
! 
! THCSecHeatAuxVarCompute: Computes secondary auxillary variables for each
!                            grid cell for heat transfer only
! author: Satish Karra
! Date: 06/5/12
!
! ************************************************************************** !
subroutine THCSecHeatAuxVarCompute(sec_heat_vars,global_aux_var, &
                                   therm_conductivity,dencpr, &
                                   option)

  use Option_module 
  use Global_Aux_module
  
  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(global_auxvar_type) :: global_aux_var
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: sec_temp(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: temp_primary_node
  PetscReal :: m
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  dm_minus = sec_heat_vars%dm_minus
  area_fm = sec_heat_vars%interfacial_area
  temp_primary_node = global_aux_var%temp(1)
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  sec_temp = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

  
  ! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0

                        
  rhs = sec_heat_vars%sec_temp  ! secondary continuum values from previous time step
  rhs(ngcells) = rhs(ngcells) + & 
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 temp_primary_node
                
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! Calculate temperature in the secondary continuum
  sec_temp(ngcells) = rhs(ngcells)/coeff_diag(ngcells)
  do i = ngcells-1, 1, -1
    sec_temp(i) = (rhs(i) - coeff_right(i)*sec_temp(i+1))/coeff_diag(i)
  enddo
  
  sec_heat_vars%sec_temp = sec_temp
            
end subroutine THCSecHeatAuxVarCompute

! ************************************************************************** !
! 
! MphaseSecHeatAuxVarCompute: Computes secondary auxillary variables in each
!                             grid cell for heat transfer only
! author: Satish Karra
! Date: 06/28/12
!
! ************************************************************************** !
subroutine MphaseSecHeatAuxVarCompute(sec_heat_vars,aux_var,global_aux_var, &
                                   therm_conductivity,dencpr, &
                                   option)

  use Option_module 
  use Global_Aux_module
  use Mphase_Aux_module
  
  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(mphase_auxvar_elem_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: sec_temp(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: temp_primary_node
  PetscReal :: m
  
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  dm_minus = sec_heat_vars%dm_minus
  area_fm = sec_heat_vars%interfacial_area
  temp_primary_node = aux_var%temp

  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  sec_temp = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

  
  ! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0
                        
  rhs = sec_heat_vars%sec_temp  ! secondary continuum values from previous time step
  rhs(ngcells) = rhs(ngcells) + & 
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 temp_primary_node
                
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! Calculate temperature in the secondary continuum
  sec_temp(ngcells) = rhs(ngcells)/coeff_diag(ngcells)
  do i = ngcells-1, 1, -1
    sec_temp(i) = (rhs(i) - coeff_right(i)*sec_temp(i+1))/coeff_diag(i)
  enddo

! print *,'temp_dcdm= ',(sec_temp(i),i=1,ngcells)
  
  sec_heat_vars%sec_temp = sec_temp


end subroutine MphaseSecHeatAuxVarCompute


end module Secondary_Continuum_Aux_module
            