
! This subroutine works only for stuctured grid 
module nhydrostatic_module
  
private 
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

#include "definitions.h"

real*8, parameter ::  fmwnacl = 58.44277D0,  fmwh2o = 18.0153D0

real*8, pointer :: hys_pres(:), hys_temp(:)

#define WEST 1
#define EAST 2
#define SOUTH 3
#define NORTH 4
#define BOTTOM 5
#define TOP 6

public nhydrostatic
public nhydrostatic3

contains

subroutine Get_Hydrosta_Pres(nz, dz, pref0, tref, dtdz, gravity, m_nacl) 

  use water_eos_module

  implicit none
  
  integer nz
  real*8 pref0, tref, dz(1:nz+1), dtdz, gravity, m_nacl
  real*8,save :: pref 
  data pref/-10.D0/ 
     
  integer n, itrho
  real*8 depth, pres0, tmp, rho0, rho1, rho, pres, dw_kg
  real*8 dzz,zz, xm_nacl, p 
! first go with z
      
  if(dabs(pref-pref0)<1D-10) return   
      
    xm_nacl = m_nacl * fmwnacl
    xm_nacl = xm_nacl /(1.D3 + xm_nacl)

    pref = pref0
    depth = 0.D0
    pres0 = pref
    tmp = tref

    rho0=rho
    p=pres0

    call nacl_den(tmp, p*1D-6, xm_nacl, dw_kg) 
    rho0 = dw_kg * 1D3

    do n = 1, nz+1
      if (n.eq.1) then
        dzz = 0.5d0*dz(1)
        zz = dzz
      elseif(n.eq.nz+1)then 
        dzz = 0.5d0*dz(nz)
        zz= zz+ dzz 
      else
        dzz = 0.5d0*(dz(n)+ dz(n-1))
        zz = zz + dzz
      endif
      tmp = tref +  dTdz*zz
    
  !   call wateos(tmp, p, rho, dw_mol, dwp, &
  !    dum, dum, dum, dum, grid%scale, ierr)
      call nacl_den(tmp, p*1D-6, xm_nacl, dw_kg) 
      rho = dw_kg * 1D3
 
                  
      itrho= 0
          
              !betap = rho * option%gravity * option%beta
      if (n==1) then
        pres = p + rho0 * gravity * dzz !+  betap * horiz
         !        call wateos(tmp, pres, rho, dw_mol, dwp, &
         !               dum, dum, dum, dum, grid%scale, ierr)
        call nacl_den(tmp, pres*1D-6, xm_nacl, dw_kg) 
        rho = dw_kg * 1D3

      else
        do
          pres = p + (rho0*dz(n-1) + rho* dz(n))/(dz(n)+dz(n-1))&
                     * gravity * dzz 
              !  call wateos(tmp, pres, rho1, dw_mol, dwp, &
              !  dum, dum, dum, dum, grid%scale, ierr)
          call nacl_den(tmp, pres*1D-6, xm_nacl, dw_kg) 
          rho1 = dw_kg * 1D3
          if (abs(rho1-rho) < 1.d-10) exit
          rho = rho1
          itrho = itrho + 1
          if (itrho > 100) then
            print *,' no convergence in hydrostat-stop',itrho,rho1,rho
            stop
          endif
        enddo
    endif
    p = pres
    rho0=rho
    hys_pres(n)=p; hys_temp(n)= tmp
!    print *, 'nhydro: ', n, p, tmp
  enddo
end subroutine  Get_Hydrosta_Pres   
     
    


subroutine nhydrostatic(realization)

  use water_eos_module
  
  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Connection_module

  implicit none

  type(realization_type) :: realization
  PetscScalar, pointer :: xx_p(:) 
  
! integer :: itrho
  integer  ierr,i,j,k
  real*8 :: depth,horiz,dw_kg,p,dp
!           dum,dwp,rho,rho1,dw_mol,zz,dzz,tmp,pres,
            

! real*8 :: ibndtyp(grid%nblkbc),cbc(grid%nblkbc),sbc(grid%nblkbc)
! real*8 :: xxbc_rec(option%ndof,grid%nblkbc),iphasebc_rec(grid%nblkbc)
  real*8 :: rho_ref, xm_nacl, pref

! real*8 :: betap,dz1,dz2,dx1,dx2, 
  integer :: natural_id,nz, local_id, ghosted_id, iconn
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_type), pointer :: cur_connection_set
  
  option => realization%option
  grid => realization%grid
  field => realization%field
  
  !Vec :: temp1_nat_vec, temp2_nat_vec, temp3_nat_vec, temp4_nat_vec

!#include "definitions.h"

! set up hydrostatic initial and boundary conditions
! dp/dz = rho * g, dp/dx = -rho*g*beta, dT/dz = dTdz
  
! dTdz = geothermal gradient                    ! C/m
! beta = lateral hydrostatic gradient           ! m/m [dh/dx = -mu vx/(k rho g)]
! tref = reference temperature at top of domain ! C
! pref = reference pressure at top of domain    ! Pa

! depth = depth of domain top from surface   ! m
! length = width of system domain ! m

! this implementation is only approximate---need to take into account
! density and viscosity change with p and T: rho = rho(T,p)

!  vz = -k/mu (dp/dz - rho g z) = 0
!  vx = -k/mu dp/dx, h = p/(rho g)

 allocate(hys_pres(grid%structured_grid%nz+1))
 allocate(hys_temp(grid%structured_grid%nz+1))


  ! this 
  if (option%iread_init /= 1) then
    call VecGetArrayF90(field%xx, xx_p, ierr)

    !do nl=1, grid%nlmax
    xm_nacl = option%m_nacl * fmwnacl
    xm_nacl = xm_nacl /(1.D3 + xm_nacl)
    call nacl_den(option%tref, option%pref*1D-6, xm_nacl, dw_kg) 
    rho_ref = dw_kg * 1D3

    do i=grid%structured_grid%istart,grid%structured_grid%iend
      do j=grid%structured_grid%jstart,grid%structured_grid%jend
        do k=grid%structured_grid%kstart,grid%structured_grid%kend
            
          local_id = 1 + (i-grid%structured_grid%istart) + &
                         (j-grid%structured_grid%jstart)*grid%structured_grid%nlx + &
                         (k-grid%structured_grid%kstart)*grid%structured_grid%nlxy  
          natural_id=grid%nL2A(local_id)+1 ! the natural ordering start from 0
        ! Get nx,ny,nz
            
          nz= floor(((real(natural_id)-.5))/grid%structured_grid%nxy) + 1

          ghosted_id=grid%nL2G(local_id)
          depth = grid%z(ghosted_id)
          horiz = grid%x(ghosted_id)
          dp = rho_ref * option%gravity(3) * option%beta * &
               ((grid%x_max - &
                 0.5d0*grid%structured_grid%dx0(grid%structured_grid%nx)) - &
               grid%x(ghosted_id))
          pref =  option%pref + dp 
          call Get_Hydrosta_Pres(grid%structured_grid%nz,grid%structured_grid%dz0, &
                                 pref,option%tref,option%dtdz,option%gravity(3), &
                                 option%m_nacl) 
          xx_p(1+ (local_id-1)*option%ndof)=hys_pres(nz) 
          xx_p(2+ (local_id-1)*option%ndof)=hys_temp(nz)
        ! print *, i,j,k,nl,na,nz, hys_pres(nz), hys_temp(nz)
        enddo
      enddo  
    enddo
    call VecRestoreArrayF90(field%xx, xx_p, ierr)
  endif
  
  ! boundary conditions

  p = option%pref
  
  
  depth = depth + 0.5d0*grid%structured_grid%dz0(grid%structured_grid%nz)
  horiz = horiz + 0.5d0*grid%structured_grid%dx0(grid%structured_grid%nx)
  
  dp = rho_ref * option%gravity(3) * option%beta * horiz
  
  if (option%myrank == 0) then
    write(*,'(" --> hydrostatic: length= ",1pe11.4,"[m], depth= ",1pe11.4, &
 &          "[m], dp= ",1pe11.4,"[Pa]")') horiz,depth,dp
    write(IUNIT2,'(" --> hydrostatic: length= ",1pe11.4,"[m], depth= ", &
           & 1pe11.4, "[m], dp= ",1pe11.4,"[Pa]")') horiz,depth,dp
  endif
  
  !save input grid%ibndtyp values
  !note-bcon must be in the order: 3D-left, right, top, bottom, front, back
  !                                1D-top, bottom

  call VecGetArrayF90(field%xx, xx_p, ierr)
  
  boundary_condition => realization%boundary_conditions%first
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection
    
    do iconn = 1, cur_connection_set%num_connections    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      natural_id = grid%nL2A(local_id)+1 
        
      if (boundary_condition%iface /= BOTTOM .and. boundary_condition%iface /= TOP) then

        select case(boundary_condition%condition%itype(1))
          case(4)
          case(3)
            natural_id = grid%nL2A(local_id)+1 
            ghosted_id = grid%nL2G(local_id)
            nz = floor(((real(natural_id)-.5D0))/grid%structured_grid%nxy) + 1
            if (boundary_condition%condition%itype(1) == 1) then
              dp = rho_ref * option%gravity(3) * option%beta * &
                   (grid%x_max-0.5d0*grid%x(grid%structured_grid%nx))
            elseif (boundary_condition%condition%itype(1) == 2) then
              dp=0.D0
            else   
              dp = rho_ref * option%gravity(3) * option%beta * &
                   ((grid%x_max - &
                     0.5d0*grid%structured_grid%dx0(grid%structured_grid%nx)) - &
                   grid%x(ghosted_id))
            endif
            pref = option%pref + dp 
            call Get_Hydrosta_Pres(grid%structured_grid%nz,grid%structured_grid%dz0, &
                                   pref,option%tref,option%dtdz,option%gravity(3), &
                                   option%m_nacl)  
            boundary_condition%aux_real_var(1,iconn) = hys_pres(nz)
            boundary_condition%aux_real_var(2:option%ndof,iconn) = &
              xx_p(2+(local_id-1)*option%ndof:local_id*option%ndof)
          case(2)
            boundary_condition%aux_real_var(:,iconn) = &
              xx_p(1+(local_id-1)*option%ndof:local_id*option%ndof)
          case(1) 
            natural_id=grid%nL2A(local_id) +1
            ghosted_id = grid%nL2G(local_id)
            nz= floor(((real(natural_id)-.5))/grid%structured_grid%nxy) + 1
            if (boundary_condition%iface == 1)then 
              dp = rho_ref * option%gravity(3) * option%beta * &
                   (grid%x_max - &
                    0.5d0*grid%structured_grid%dx0(grid%structured_grid%nx))
            elseif(boundary_condition%iface == 2)then
              dp=0.D0
            else   
              dp = rho_ref * option%gravity(3) * option%beta * &
                   ((grid%x_max - &
                     0.5d0*grid%structured_grid%dx0(grid%structured_grid%nx)) - &
                   grid%x(ghosted_id))
            endif
            pref =  option%pref + dp
            call Get_Hydrosta_Pres(grid%structured_grid%nz,grid%structured_grid%dz0, &
                                   pref,option%tref,option%dtdz,option%gravity(3), &
                                   option%m_nacl)  
            boundary_condition%aux_real_var(1,iconn) = hys_pres(nz)
            boundary_condition%aux_real_var(2,iconn) = hys_temp(nz)
          end select 
        else if(boundary_condition%iface == 4) then
          select case(boundary_condition%condition%itype(1))
            case(4)
            case(3)
              nz = grid%structured_grid%nz
              ghosted_id = grid%nL2G(local_id)
              dp = rho_ref * option%gravity(3) * option%beta * &
                   ((grid%x_max - &
                     0.5d0*grid%structured_grid%dx0(grid%structured_grid%nx)) - &
                   grid%x(ghosted_id))
              pref =  option%pref + dp 
              call Get_Hydrosta_Pres(grid%structured_grid%nz,grid%structured_grid%dz0, &
                                     pref,option%tref,option%dtdz,option%gravity(3), &
                                     option%m_nacl)  
              boundary_condition%aux_real_var(1,iconn) = hys_pres(nz+1)
              boundary_condition%aux_real_var(2:option%ndof,iconn) = &
                xx_p(2+ (local_id-1)*option%ndof: local_id *option%ndof)
            case(2)
              boundary_condition%aux_real_var(:,iconn) = &
                xx_p(1+ (local_id-1)*option%ndof:local_id*option%ndof)
            case(1) 
              nz = grid%structured_grid%nz
              ghosted_id = grid%nL2G(local_id)
              dp = rho_ref * option%gravity(3) * option%beta * &
                   ((grid%x_max - &
                     0.5d0*grid%structured_grid%dx0(grid%structured_grid%nx)) - &
                   grid%x(natural_id))
              pref =  option%pref + dp 
              call Get_Hydrosta_Pres(grid%structured_grid%nz,grid%structured_grid%dz0, &
                                     pref,option%tref,option%dtdz,option%gravity(3), &
                                     option%m_nacl)  
              boundary_condition%aux_real_var(1,iconn) = hys_pres(nz+1)
              boundary_condition%aux_real_var(2,iconn) = hys_temp(nz+1)
          end select 
        elseif(boundary_condition%iface == 3) then
          select case(boundary_condition%condition%itype(1))
            case(4)
            case(3,1)
              natural_id=grid%nL2A(local_id)+1
              ghosted_id = grid%nL2G(local_id)
              dp = rho_ref * option%gravity(3) * option%beta * &
                   ((grid%x_max - &
                     0.5d0*grid%structured_grid%dx0(grid%structured_grid%nx)) - &
                   grid%x(ghosted_id))
              pref =  option%pref + dp 
              boundary_condition%aux_real_var(1,iconn) = pref
              boundary_condition%aux_real_var(2,iconn) = option%tref
            case(2)
              boundary_condition%aux_real_var(:,iconn) = &
                xx_p(1+ (local_id-1)*option%ndof:local_id*option%ndof)
        end select
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo
  call VecRestoreArrayF90(field%xx, xx_p, ierr)
 
  deallocate(hys_pres)
  deallocate(hys_temp)

end subroutine nhydrostatic
  
  
! clu: developed for 9M node SACROC case only.
!      it is a right hand coordinate system    
subroutine nhydrostatic3(grid)

  use pflow_gridtype_module
  use water_eos_module

  implicit none
  
! mostly copied from condition.F90 :ComputeHydrostaticPressure   
  type(pflowGrid) :: grid
#if 0  
  real*8 :: ref_coordinate(3), cell_coordinate(3)
  real*8 :: reference_pressure, pressure_gradient_X, pressure_gradient_Y,  &
            pressure_gradient_Z
  real*8 :: reference_temperature, temperature_gradient_X, &
            temperature_gradient_Y, temperature_gradient_Z

  integer :: i, num_increments, ierr
  real*8 :: dx, dy, dz, z, dum, rho, dw_mol, dwp
  real*8 :: temperature_at_xy, temperature_at_xyz
  real*8 :: pressure_at_xy, pressure_at_xyz
  real*8 :: increment, final_increment
  PetscScalar, pointer :: xx_p(:) 
  integer n, ng, nc, m, ibc
   
  ref_coordinate = grid%hydro_ref_xyz
  reference_pressure = grid% pref
  pressure_gradient_X = option%beta *1.01325D4
  pressure_gradient_Y =0.D0
  pressure_gradient_Z =0.D0
  reference_temperature = option%tref
  temperature_gradient_X = 0.D0
  temperature_gradient_Y = 0.D0
  temperature_gradient_Z = grid%dTdz
   
  
  call VecGetArrayF90(field%xx, xx_p, ierr)
  
  do n=1, grid%nlmax
       ng = grid%nL2G(n)
       cell_coordinate(1)= grid%x(ng)
       cell_coordinate(2)= grid%y(ng)
       cell_coordinate(3)= grid%z(ng)

        dx = cell_coordinate(1) - ref_coordinate(1)
        dy = cell_coordinate(2) - ref_coordinate(2)
        dz = cell_coordinate(3) - ref_coordinate(3)

        pressure_at_xy = reference_pressure + dx*pressure_gradient_X + &
                                        dy*pressure_gradient_Y
        temperature_at_xy = reference_temperature + dx*temperature_gradient_X + &
                                                dy*temperature_gradient_Y

          num_increments = int(abs(dz))  ! 1m increments
          final_increment = sign(abs(dz)-num_increments*1.d0,-dz)
          increment = 1.d0
         if (ref_coordinate(3) < cell_coordinate(3)) increment = -1.d0

          z = ref_coordinate(3)
          pressure_at_xyz = pressure_at_xy
          do i=1, num_increments
             z = z - increment
             temperature_at_xyz = temperature_at_xy + (z-ref_coordinate(3))* &
                                                temperature_gradient_Z
             call wateos(temperature_at_xyz,pressure_at_xyz,rho,dw_mol,dwp, &
                  dum,dum,dum,dum,grid%scale,ierr)
              pressure_at_xyz = pressure_at_xyz + increment*rho*option%gravity  
          enddo
         if (final_increment > 0.d0) then
          z = z - final_increment
          temperature_at_xyz = temperature_at_xy + (z-ref_coordinate(3))* &
                                               temperature_gradient_Z
         call wateos(temperature_at_xyz,pressure_at_xyz,rho,dw_mol,dwp, &
                  dum,dum,dum,dum,grid%scale,ierr)
          pressure_at_xyz = pressure_at_xyz + final_increment*rho*option%gravity  
         endif
          xx_p((n-1)*option%ndof +1) = pressure_at_xyz
          xx_p((n-1)*option%ndof +2) = temperature_at_xyz
     enddo    
       
     do nc = 1, grid%nconnbc
        m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
        ibc= grid%ibconn(nc)
        if(grid%ibndtyp(ibc) == 3)then
         field%xxbc(1,nc) = xx_p((m-1)*option%ndof +1)
         field%xxbc(2,nc) = xx_p((m-1)*option%ndof +2)     
     !    print *, nc,  field%xxbc(:,nc)
        endif
     enddo
  call VecRestoreArrayF90(field%xx, xx_p, ierr)
#endif  
end subroutine nhydrostatic3
  
end module nhydrostatic_module

