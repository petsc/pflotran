! ************************************************************************** !
!
! MiscibleResidualPatch: Computes the residual equation at patch level
!                      original version (not used)
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine MiscibleResidualPatch(snes,xx,r,realization,ierr)

  use Connection_module
  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: i, jn
  PetscInt :: ip1, ip2
  PetscInt :: local_id, ghosted_id, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer ::accum_p(:)

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:),&
               tortuosity_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
                          
               
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)

  PetscInt :: iphase
  PetscInt :: icap_up, icap_dn, ithrm_up, ithrm_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: dd, f_up, f_dn, ff
  PetscReal :: perm_up, perm_dn
  PetscReal :: D_up, D_dn  ! "Diffusion" constants at upstream, downstream faces.
  PetscReal :: dw_kg, dw_mol,dddt,dddp
  PetscReal :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1
  PetscReal :: rho, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt, dvdp, xphi
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof), v_darcy(realization%option%nphase)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscReal :: psrc(1:realization%option%nphase)
  PetscViewer :: viewer
  PetscInt :: nsrcpara
  PetscReal, pointer :: msrc(:)
  
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(Miscible_parameter_type), pointer :: Miscible_parameter
  type(Miscible_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscBool :: enthalpy_flag
  PetscInt :: ng
  PetscInt :: iconn, idof, istart, iend
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal, pointer :: Resold_AR(:), Resold_FL(:), delx(:)
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  Miscible_parameter => patch%aux%Miscible%Miscible_parameter
  aux_vars => patch%aux%Miscible%aux_vars
  aux_vars_bc => patch%aux%Miscible%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

 ! call MiscibleUpdateAuxVarsPatchNinc(realization)
  ! override flags since they will soon be out of date  
 ! patch%MiscibleAux%aux_vars_up_to_date = PETSC_FALSE 
 
  if (option%compute_mass_balance_new) then
    call MiscibleZeroMassBalDeltaPatch(realization)
  endif

! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr)
 
  call VecGetArrayF90(field%flow_yy,yy_p,ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
!  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  allocate(Resold_AR(option%nflowdof), Resold_FL(option%nflowdof), delx(option%nflowdof))
 
! Multiphase flash calculation is more expensive, so calculate once per iterration
#if 1
  ! Pertubations for aux terms --------------------------------
  do ng = 1, grid%ngmax
    if(grid%nG2L(ng)<0)cycle
    if (associated(patch%imat)) then
      if (patch%imat(ng) <= 0) cycle
    endif
    ghosted_id = ng   
    istart =  (ng-1) * option%nflowdof +1 ; iend = istart -1 + option%nflowdof
     ! iphase =int(iphase_loc_p(ng))
    call MiscibleAuxVarCompute_Ninc(xx_loc_p(istart:iend),aux_vars(ng)%aux_var_elem(0),&
          global_aux_vars(ng),&
          realization%fluid_properties,option)
!    print *,'flash ', xx_loc_p(istart:iend),aux_vars(ng)%aux_var_elem(0)%den
#if 1
    if (associated(global_aux_vars)) then
      global_aux_vars(ghosted_id)%pres(:) = aux_vars(ghosted_id)%aux_var_elem(0)%pres -&
               aux_vars(ghosted_id)%aux_var_elem(0)%pc(:)
!     global_aux_vars(ghosted_id)%temp(:) = aux_vars(ghosted_id)%aux_var_elem(0)%temp
      global_aux_vars(ghosted_id)%sat(:) = 1D0
!     global_aux_vars(ghosted_id)%sat_store =
      global_aux_vars(ghosted_id)%den(:) = aux_vars(ghosted_id)%aux_var_elem(0)%den(:)
      global_aux_vars(ghosted_id)%den_kg(:) = aux_vars(ghosted_id)%aux_var_elem(0)%den(:) &
                                          * aux_vars(ghosted_id)%aux_var_elem(0)%avgmw(:)
!     global_aux_vars(ghosted_id)%reaction_rate(:) = 0D0
!     global_aux_vars(ghosted_id)%pres(:)
    else
      print *,'Not associated global for Miscible'
    endif
#endif

    if (option%numerical_derivatives) then
      delx(1) = xx_loc_p((ng-1)*option%nflowdof+1)*dfac*1.D-3
 
      do idof = 2, option%nflowdof
        if(xx_loc_p((ng-1)*option%nflowdof+idof) <= 0.9) then
            delx(idof) = dfac*xx_loc_p((ng-1)*option%nflowdof+idof)*1D1 
        else
          delx(idof) = -dfac*xx_loc_p((ng-1)*option%nflowdof+idof)*1D1 
        endif
        if (delx(idof) < 1D-8 .and. delx(idof)>=0.D0) delx(idof) = 1D-8
        if (delx(idof) >-1D-8 .and. delx(idof)<0.D0) delx(idof) =-1D-8

        if ((delx(idof)+xx_loc_p((ng-1)*option%nflowdof+idof)) > 1.D0) then
          delx(idof) = (1.D0-xx_loc_p((ng-1)*option%nflowdof+idof))*1D-4
        endif
        if ((delx(idof)+xx_loc_p((ng-1)*option%nflowdof+idof)) < 0.D0) then
          delx(idof) = xx_loc_p((ng-1)*option%nflowdof+idof)*1D-4
        endif
      end do

      patch%aux%Miscible%delx(:,ng)=delx(:)
      call MiscibleAuxVarCompute_Winc(xx_loc_p(istart:iend),delx(:),&
            aux_vars(ng)%aux_var_elem(1:option%nflowdof),global_aux_vars(ng),&
            realization%fluid_properties,option)
!      if (aux_vars(ng)%aux_var_elem(option%nflowdof)%sat(2)>1D-8 .and. &
!           aux_vars(ng)%aux_var_elem(0)%sat(2)<1D-12)then
         
!      endif   
    endif
  enddo
#endif

  Resold_AR=0.D0; ResOld_FL=0.D0; r_p = 0.d0
  patch%aux%Miscible%Resold_AR=0.D0
  patch%aux%Miscible%Resold_BC=0.D0
  patch%aux%Miscible%ResOld_FL=0.D0
   
#if 1
  ! Accumulation terms ------------------------------------
  r_p = - accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    call MiscibleAccumulation(aux_vars(ghosted_id)%aux_var_elem(0),&
                            global_aux_vars(ghosted_id), &
                            porosity_loc_p(ghosted_id), &
                            volume_p(local_id), &
                            Miscible_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                            option,Res) 
    r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
    print *,'REs, acm: ', res
    patch%aux%Miscible%Resold_AR(local_id, :)= Res(1:option%nflowdof)
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  sum_connection = 0 
  do 
    if (.not.associated(source_sink)) exit
    !print *, 'RES s/s begin'
    ! check whether enthalpy dof is included
  !  if (source_sink%flow_condition%num_sub_conditions > 3) then
      enthalpy_flag = PETSC_TRUE
   ! else
   !   enthalpy_flag = PETSC_FALSE
   ! endif
    if (associated(source_sink%flow_condition%pressure)) then   
      psrc(:) = source_sink%flow_condition%pressure%flow_dataset%time_series%cur_value(:)
    endif 
!   qsrc1 = source_sink%flow_condition%pressure%flow_dataset%time_series%cur_value(1)
    tsrc1 = source_sink%flow_condition%temperature%flow_dataset%time_series%cur_value(1)
    csrc1 = source_sink%flow_condition%concentration%flow_dataset%time_series%cur_value(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%flow_dataset%time_series%cur_value(1)
!    hsrc1=0D0
!    qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
!    csrc1 = csrc1 / FMWCO2
!    msrc(1)=qsrc1; msrc(2) =csrc1
!    msrc(:)= psrc(:)

    select case(source_sink%flow_condition%itype(1))
      case(MASS_RATE_SS)
        msrc => source_sink%flow_condition%rate%flow_dataset%time_series%cur_value
        nsrcpara= 2
      case(WELL_SS) ! Well not implemented yet
        msrc => source_sink%flow_condition%well%flow_dataset%time_series%cur_value
        nsrcpara = 7 + option%nflowspec 
      case default
        print *, 'Flash mode does not support source/sink type: ', source_sink%flow_condition%itype(1)
        stop  
    end select

    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      call MiscibleSourceSink(msrc,nsrcpara,psrc,tsrc1,hsrc1,csrc1,aux_vars(ghosted_id)%aux_var_elem(0),&
                            source_sink%flow_condition%itype(1),Res, &
                            patch%ss_fluid_fluxes(:,sum_connection), &
                            enthalpy_flag, option)
 
      r_p((local_id-1)*option%nflowdof + jh2o) = r_p((local_id-1)*option%nflowdof + jh2o)-Res(jh2o)
      r_p((local_id-1)*option%nflowdof + jglyc) = r_p((local_id-1)*option%nflowdof + jglyc)-Res(jglyc)
      patch%aux%Miscible%Resold_AR(local_id,jh2o)= patch%aux%Miscible%Resold_AR(local_id,jh2o) - Res(jh2o)    
      patch%aux%Miscible%Resold_AR(local_id,jglyc)= patch%aux%Miscible%Resold_AR(local_id,jglyc) - Res(jglyc)    
#if 0
      if (enthalpy_flag)then
        r_p( local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - Res(option%nflowdof)
        patch%aux%Miscible%Resold_AR(local_id,option%nflowdof)=&
          patch%aux%Miscible%Resold_AR(local_id,option%nflowdof) - Res(option%nflowdof)
       endif 
#endif
  !  else if (qsrc1 < 0.d0) then ! withdrawal
  !  endif
    enddo
    source_sink => source_sink%next
  enddo
#endif
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
        
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn = int(ithrm_loc_p(ghosted_id))
      D_dn = Miscible_parameter%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))

      icap_dn = int(icap_loc_p(ghosted_id))
      
!     Then need fill up increments for BCs
      do idof = 1, option%nflowdof   
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(HYDROSTATIC_BC)
            xxbc(1) = boundary_condition%flow_aux_real_var(1,iconn)
            if (idof >= 2) then
              xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
            endif 
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
!           iphase = int(iphase_loc_p(ghosted_id))
        end select
      enddo

 
      call MiscibleAuxVarCompute_Ninc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0),&
           global_aux_vars_bc(sum_connection),&
           realization%fluid_properties, option)
#if 1
      if (associated(global_aux_vars_bc)) then
        global_aux_vars_bc(sum_connection)%pres(:)= aux_vars_bc(sum_connection)%aux_var_elem(0)%pres
!       global_aux_vars_bc(sum_connection)%temp(:)=aux_vars_bc(sum_connection)%aux_var_elem(0)%temp
        global_aux_vars_bc(sum_connection)%sat(:)= 1.0
      !    global_aux_vars(ghosted_id)%sat_store = 
        global_aux_vars_bc(sum_connection)%den(:)=aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:)
        global_aux_vars_bc(sum_connection)%den_kg = aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:) &
                                          * aux_vars_bc(sum_connection)%aux_var_elem(0)%avgmw(:)
  !     global_aux_vars(ghosted_id)%den_kg_store
      endif
#endif

      call MiscibleBCFlux(boundary_condition%flow_condition%itype, &
         boundary_condition%flow_aux_real_var(:,iconn), &
         aux_vars_bc(sum_connection)%aux_var_elem(0), &
         aux_vars(ghosted_id)%aux_var_elem(0), &
         porosity_loc_p(ghosted_id), &
         tortuosity_loc_p(ghosted_id), &
         Miscible_parameter%sir(:,icap_dn), &
         cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
         cur_connection_set%area(iconn), &
         distance_gravity,option, &
         v_darcy,Res)
      patch%boundary_velocities(:,sum_connection) = v_darcy(:)
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:option%nflowdof)
      patch%aux%Miscible%Resold_AR(local_id,1:option%nflowdof) = &
        patch%aux%Miscible%ResOld_AR(local_id,1:option%nflowdof) - Res(1:option%nflowdof)
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif

#if 1
  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id_up) <= 0 .or.  &
            patch%imat(ghosted_id_dn) <= 0) cycle
      endif

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)
        
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))

      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
   
      D_up = Miscible_parameter%ckwet(ithrm_up)
      D_dn = Miscible_parameter%ckwet(ithrm_dn)

      call MiscibleFlux(aux_vars(ghosted_id_up)%aux_var_elem(0),porosity_loc_p(ghosted_id_up), &
                          tortuosity_loc_p(ghosted_id_up),Miscible_parameter%sir(:,icap_up), &
                          dd_up,perm_up,D_up, &
                          aux_vars(ghosted_id_dn)%aux_var_elem(0),porosity_loc_p(ghosted_id_dn), &
                          tortuosity_loc_p(ghosted_id_dn),Miscible_parameter%sir(:,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                          cur_connection_set%area(iconn),distance_gravity, &
                          upweight,option,v_darcy,Res)

      patch%internal_velocities(:,sum_connection) = v_darcy(:)
      patch%aux%Miscible%Resold_FL(sum_connection,1:option%nflowdof)= Res(1:option%nflowdof)
 
     if (local_id_up>0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
      endif
   
      if (local_id_dn>0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
#endif

! adjust residual to R/dt
  select case (option%idt_switch) 
  case(1) 
     r_p(:) = r_p(:)/option%flow_dt
  case(-1)
     if(option%flow_dt>1.D0) r_p(:) = r_p(:)/option%flow_dt
  end select
  
  do local_id = 1, grid%nlmax
     if (associated(patch%imat)) then
        if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
     endif

     istart = 1 + (local_id-1)*option%nflowdof
     if(volume_p(local_id)>1.D0) r_p (istart:istart+2)=r_p(istart:istart+2)/volume_p(local_id)
     if(r_p(istart) >1E20 .or. r_p(istart) <-1E20) print *, r_p (istart:istart+2)
!     print *,'flash res', local_id, r_p (istart:istart+2)
  enddo

! print *,'finished rp vol scale'


  if (patch%aux%Miscible%inactive_cells_exist) then
    do i=1,patch%aux%Miscible%n_zero_rows
      r_p(patch%aux%Miscible%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(field%flow_yy, yy_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
!  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  deallocate(Resold_AR, Resold_FL, delx)
  
  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(option%mycomm,'Rresidual.out',viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(option%mycomm,'Rxx.out',viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
end subroutine MiscibleResidualPatch


! ************************************************************************** !
!
! MiscibleJacobianPatch: Computes the Jacobian
! author: Chuan Lu
! date: 10/13/08
!
! ************************************************************************** !
subroutine MiscibleJacobianPatch(snes,xx,A,B,flag,realization,ierr)

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_module
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  
  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr
  PetscInt :: nvar,neq,nr
  PetscInt :: ithrm_up, ithrm_dn, i
  PetscInt :: ip1, ip2 

  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), &
                          xx_loc_p(:), tortuosity_loc_p(:),&
                          perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)
  PetscInt :: icap,iphas,iphas_up,iphas_dn,icap_up,icap_dn
  PetscInt :: ii, jj
  PetscReal :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho
  PetscReal :: tsrc1,qsrc1,csrc1,hsrc1
  PetscReal :: dd_up, dd_dn, dd, f_up, f_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: dw_dp,dw_dt,hw_dp,hw_dt,dresT_dp,dresT_dt
  PetscReal :: D_up, D_dn  ! "Diffusion" constants upstream and downstream of a face.
  PetscReal :: zero, norm
  PetscReal :: upweight
  PetscReal :: max_dev  
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: natural_id_up,natural_id_dn
  
  PetscReal :: Jup(1:realization%option%nflowdof,1:realization%option%nflowdof), &
            Jdn(1:realization%option%nflowdof,1:realization%option%nflowdof)
  
  PetscInt :: istart, iend
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscBool :: enthalpy_flag
  PetscInt :: iconn, idof
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: Res(realization%option%nflowdof) 
  PetscReal :: xxbc(1:realization%option%nflowdof), delxbc(1:realization%option%nflowdof)
  PetscReal :: ResInc(realization%patch%grid%nlmax,realization%option%nflowdof, &
           realization%option%nflowdof)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(Miscible_parameter_type), pointer :: Miscible_parameter
  type(Miscible_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)

  PetscReal :: vv_darcy(realization%option%nphase), voltemp
  PetscReal :: ra(1:realization%option%nflowdof,1:realization%option%nflowdof*2)
  PetscInt nsrcpara 
  PetscReal, pointer :: msrc(:)
  PetscReal :: psrc(1:realization%option%nphase), ss_flow(1:realization%option%nphase)
  PetscReal :: dddt, dddp, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt,&
               dvdp, xphi
  PetscInt :: iphasebc                
  
  PetscViewer :: viewer
  Vec :: debug_vec
!-----------------------------------------------------------------------
! R stand for residual
!  ra       1              2              3              4          5              6            7      8
! 1: p     dR/dpi         dR/dTi          dR/dci        dR/dsi   dR/dpim        dR/dTim
! 2: T
! 3: c
! 4  s         
!-----------------------------------------------------------------------

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  Miscible_parameter => patch%aux%Miscible%Miscible_parameter
  aux_vars => patch%aux%Miscible%aux_vars
  aux_vars_bc => patch%aux%Miscible%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s
!  flag = SAME_NONZERO_PATTERN

#if 0
!  call MiscibleNumericalJacobianTest(xx,realization)
#endif

 ! print *,'*********** In Jacobian ********************** '
  call MatZeroEntries(A,ierr)

  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
!  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)

 ResInc = 0.D0
#if 1
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
     ghosted_id = grid%nL2G(local_id)
     !geh - Ignore inactive cells with inactive materials
     if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
     endif
     iend = local_id*option%nflowdof
     istart = iend-option%nflowdof+1
     icap = int(icap_loc_p(ghosted_id))
     
     do nvar =1, option%nflowdof
        call MiscibleAccumulation(aux_vars(ghosted_id)%aux_var_elem(nvar), &
             global_aux_vars(ghosted_id),& 
             porosity_loc_p(ghosted_id), &
             volume_p(local_id), &
             Miscible_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
             option,res) 
        ResInc( local_id,:,nvar) =  ResInc(local_id,:,nvar) + Res(:)
     enddo
     
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
  !  if (source_sink%flow_condition%num_sub_conditions > 3) then
      enthalpy_flag = PETSC_FALSE
   ! else
   !   enthalpy_flag = PETSC_FALSE
   ! endif
    if (associated(source_sink%flow_condition%pressure)) then
      psrc(:) = source_sink%flow_condition%pressure%flow_dataset%time_series%cur_value(:)
    endif
!   tsrc1 = source_sink%flow_condition%temperature%flow_dataset%time_series%cur_value(1)
    csrc1 = source_sink%flow_condition%concentration%flow_dataset%time_series%cur_value(1)
!   hsrc1=0.D0
!   if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%flow_dataset%time_series%cur_value(1)

   ! qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
   ! csrc1 = csrc1 / FMWCO2
    select case(source_sink%flow_condition%itype(1))
      case(MASS_RATE_SS)
        msrc => source_sink%flow_condition%rate%flow_dataset%time_series%cur_value
        nsrcpara= 2
      case(WELL_SS)
        msrc => source_sink%flow_condition%well%flow_dataset%time_series%cur_value
        nsrcpara = 7 + option%nflowspec 
      case default
        print *, 'Flash mode does not support source/sink type: ', source_sink%flow_condition%itype(1)
        stop  
    end select
 
      cur_connection_set => source_sink%connection_set
 
    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1 
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
!      if (enthalpy_flag) then
!        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - hsrc1 * option%flow_dt   
!      endif         
     do nvar =1, option%nflowdof
       call MiscibleSourceSink(msrc,nsrcpara,psrc,tsrc1,hsrc1,csrc1, aux_vars(ghosted_id)%aux_var_elem(nvar),&
                            source_sink%flow_condition%itype(1), Res,&
                            ss_flow, &
                            enthalpy_flag, option)
      
       ResInc(local_id,jh2o,nvar)=  ResInc(local_id,jh2o,nvar) - Res(jh2o)
       ResInc(local_id,jglyc,nvar)=  ResInc(local_id,jglyc,nvar) - Res(jglyc)
#if 0
       if (enthalpy_flag) & 
           ResInc(local_id,option%nflowdof,nvar)=&
           ResInc(local_id,option%nflowdof,nvar)- Res(option%nflowdof) 
#endif
     
     enddo 
    enddo
    source_sink => source_sink%next
  enddo
#endif
! Boundary conditions
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn = int(ithrm_loc_p(ghosted_id))
      D_dn = Miscible_parameter%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      icap_dn = int(icap_loc_p(ghosted_id))

! Then need fill up increments for BCs
    delxbc=0.D0;
    do idof =1, option%nflowdof   
       select case(boundary_condition%flow_condition%itype(idof))
       case(DIRICHLET_BC)
          xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          delxbc(idof)=0.D0
      case(HYDROSTATIC_BC)
          xxbc(1) = boundary_condition%flow_aux_real_var(1,iconn)
          if(idof>=2)then
             xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
             delxbc(idof)=patch%aux%Miscible%delx(idof,ghosted_id)
          endif 
       case(NEUMANN_BC, ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
          xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
          !iphasebc = int(iphase_loc_p(ghosted_id))
          delxbc(idof)=patch%aux%Miscible%delx(idof,ghosted_id)
       end select
    enddo
    !print *,'BC:',boundary_condition%flow_condition%itype, xxbc, delxbc

 
    call MiscibleAuxVarCompute_Ninc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0),&
         global_aux_vars_bc(sum_connection),&
         realization%fluid_properties, option)
    call MiscibleAuxVarCompute_Winc(xxbc,delxbc,&
         aux_vars_bc(sum_connection)%aux_var_elem(1:option%nflowdof),&
         global_aux_vars_bc(sum_connection),&
         realization%fluid_properties,option)
    
    do nvar=1,option%nflowdof
       call MiscibleBCFlux(boundary_condition%flow_condition%itype, &
         boundary_condition%flow_aux_real_var(:,iconn), &
         aux_vars_bc(sum_connection)%aux_var_elem(nvar), &
         aux_vars(ghosted_id)%aux_var_elem(nvar), &
         porosity_loc_p(ghosted_id), &
         tortuosity_loc_p(ghosted_id), &
         Miscible_parameter%sir(:,icap_dn), &
         cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
         cur_connection_set%area(iconn), &
         distance_gravity,option, &
         vv_darcy,Res)
       ResInc(local_id,1:option%nflowdof,nvar) = ResInc(local_id,1:option%nflowdof,nvar) - Res(1:option%nflowdof)
    enddo
 enddo
    boundary_condition => boundary_condition%next
 enddo
#endif
! Set matrix values related to single node terms: Accumulation, Source/Sink, BC
  do local_id = 1, grid%nlmax  ! For each local node do...
     ghosted_id = grid%nL2G(local_id)
     !geh - Ignore inactive cells with inactive materials
     if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
     endif

     ra=0.D0
     max_dev=0.D0
     do neq=1, option%nflowdof
        do nvar=1, option%nflowdof
           ra(neq,nvar)=(ResInc(local_id,neq,nvar)-patch%aux%Miscible%ResOld_AR(local_id,neq))&
              /patch%aux%Miscible%delx(nvar,ghosted_id)
           if(max_dev < dabs(ra(3,nvar))) max_dev = dabs(ra(3,nvar))
        enddo
     enddo
   
   select case(option%idt_switch)
      case(1) 
        ra(1:option%nflowdof,1:option%nflowdof) =ra(1:option%nflowdof,1:option%nflowdof) /option%flow_dt
      case(-1)
        if(option%flow_dt>1) ra(1:option%nflowdof,1:option%nflowdof) =ra(1:option%nflowdof,1:) /option%flow_dt
    end select

     Jup=ra(1:option%nflowdof,1:option%nflowdof)
     if(volume_p(local_id)>1.D0 ) Jup=Jup / volume_p(local_id)
   
!      if(local_id==1) print *, 'flash jac', volume_p(local_id), ra
     call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)
  end do

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_srcsink.out',viewer,ierr)
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 1
  ! Interior Flux Terms -----------------------------------  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0    
  ResInc = 0.D0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id_up) <= 0 .or. &
            patch%imat(ghosted_id_dn) <= 0) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
     ! natural_id_up = grid%nG2N(ghosted_id_up)
     ! natural_id_dn = grid%nG2N(ghosted_id_dn)
   
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)
    
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))
    
      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      D_up = Miscible_parameter%ckwet(ithrm_up)
      D_dn = Miscible_parameter%ckwet(ithrm_dn)
    
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
      
      do nvar = 1, option%nflowdof 
         call MiscibleFlux(aux_vars(ghosted_id_up)%aux_var_elem(nvar),porosity_loc_p(ghosted_id_up), &
                          tortuosity_loc_p(ghosted_id_up),Miscible_parameter%sir(:,icap_up), &
                          dd_up,perm_up,D_up, &
                          aux_vars(ghosted_id_dn)%aux_var_elem(0),porosity_loc_p(ghosted_id_dn), &
                          tortuosity_loc_p(ghosted_id_dn),Miscible_parameter%sir(:,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                          cur_connection_set%area(iconn),distance_gravity, &
                          upweight, option, vv_darcy, Res)
            ra(:,nvar)= (Res(:)-patch%aux%Miscible%ResOld_FL(iconn,:))&
              /patch%aux%Miscible%delx(nvar,ghosted_id_up)
         call MiscibleFlux(aux_vars(ghosted_id_up)%aux_var_elem(0),porosity_loc_p(ghosted_id_up), &
                          tortuosity_loc_p(ghosted_id_up),Miscible_parameter%sir(:,icap_up), &
                          dd_up,perm_up,D_up, &
                          aux_vars(ghosted_id_dn)%aux_var_elem(nvar),porosity_loc_p(ghosted_id_dn),&
                          tortuosity_loc_p(ghosted_id_dn),Miscible_parameter%sir(:,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                          cur_connection_set%area(iconn),distance_gravity, &
                          upweight, option, vv_darcy, Res)
         ra(:,nvar+option%nflowdof)= (Res(:)-patch%aux%Miscible%ResOld_FL(iconn,:))&
           /patch%aux%Miscible%delx(nvar,ghosted_id_dn)
    enddo

    select case(option%idt_switch)
    case(1)
       ra =ra / option%flow_dt
    case(-1)  
       if(option%flow_dt>1)  ra =ra / option%flow_dt
    end select
    
    if (local_id_up > 0) then
       voltemp=1.D0
       if(volume_p(local_id_up)>1.D0)then
         voltemp = 1.D0/volume_p(local_id_up)
       endif
       Jup(:,1:option%nflowdof)= ra(:,1:option%nflowdof)*voltemp !11
       jdn(:,1:option%nflowdof)= ra(:, 1 + option%nflowdof:2 * option%nflowdof)*voltemp !12

       call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
            Jup,ADD_VALUES,ierr)
       call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
            Jdn,ADD_VALUES,ierr)
    endif
    if (local_id_dn > 0) then
       voltemp=1.D0
       if(volume_p(local_id_dn)>1.D0)then
         voltemp=1.D0/volume_p(local_id_dn)
       endif
       Jup(:,1:option%nflowdof)= -ra(:,1:option%nflowdof)*voltemp !21
       jdn(:,1:option%nflowdof)= -ra(:, 1 + option%nflowdof:2 * option%nflowdof)*voltemp !22

 
       call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
            Jdn,ADD_VALUES,ierr)
       call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
            Jup,ADD_VALUES,ierr)
    endif
 enddo
    cur_connection_set => cur_connection_set%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
 ! print *,'end inter flux'
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_flux.out',viewer,ierr)
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 0
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_bcflux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#endif
  
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr)

   
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
! call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
! print *,'end jac'
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
 ! call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)

  if (patch%aux%Miscible%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,patch%aux%Miscible%n_zero_rows, &
                          patch%aux%Miscible%zero_rows_local_ghosted,f_up, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  endif

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(option%mycomm,'Rjacobian.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    call MatNorm(A,NORM_1,norm,ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(A,NORM_FROBENIUS,norm,ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(A,NORM_INFINITY,norm,ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option)
!    call GridCreateVector(grid,ONEDOF,debug_vec,GLOBAL)
!    call MatGetRowMaxAbs(A,debug_vec,PETSC_NULL_INTEGER,ierr)
!    call VecMax(debug_vec,i,norm,ierr)
!    call VecDestroy(debug_vec,ierr)
  endif
end subroutine MiscibleJacobianPatch

! ************************************************************************** !
!
! MiscibleBCFluxAdv: Computes the  boundary flux terms for the residual (not used)
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine MiscibleBCFluxAdv(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
     por_dn,tor_dn,dd_up,perm_dn,Dk_dn, &
     area,dist_gravity,option,vv_darcy,Res)
  use Option_module
  
  implicit none
  
  PetscInt :: ibndtype(:)
  type(Miscible_auxvar_elem_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: vv_darcy(:), area
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec, np
  PetscReal :: fluxm(option%nflowspec),fluxe,q,density_ave, v_darcy
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi
  
  fluxm = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

! Flow   
! diffdp = por_dn*tor_dn/dd_up*area
  do np = 1, option%nphase  
    select case(ibndtype(1))
        ! figure out the direction of flow
      case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
        Dq = perm_dn / dd_up
        ! Flow term
        ukvr=0.D0
        v_darcy=0.D0 
          upweight=1.D0
          if (aux_var_up%sat(np) < eps) then 
            upweight=0.d0
          else if (aux_var_dn%sat(np) < eps) then 
              upweight=1.d0
          endif
          density_ave = upweight*aux_var_up%den(np) + (1.D0-upweight)*aux_var_dn%den(np)
           
          gravity = (upweight*aux_var_up%den(np) * aux_var_up%avgmw(np) + &
                (1.D0-upweight)*aux_var_dn%den(np) * aux_var_dn%avgmw(np)) &
                * dist_gravity
       
          dphi = aux_var_up%pres - aux_var_dn%pres &
                - aux_var_up%pc(np) + aux_var_dn%pc(np) &
                + gravity
   
          if (dphi>=0.D0) then
            ukvr = aux_var_up%kvr(np)
          else
            ukvr = aux_var_dn%kvr(np)
          endif
     
          if (ukvr*Dq>floweps) then
            v_darcy = Dq * ukvr * dphi
          endif

      case(NEUMANN_BC)
        v_darcy = 0.D0
        if (dabs(aux_vars(1)) > floweps) then
          v_darcy = aux_vars(MIS_PRESSURE_DOF)
          if (v_darcy > 0.d0) then 
            density_ave = aux_var_up%den(np)
          else 
            density_ave = aux_var_dn%den(np)
          endif
        endif

    end select
     
    q = v_darcy*area
    vv_darcy(np) = v_darcy
    uxmol = 0.D0
     
    if (v_darcy >= 0.D0) then
      uxmol(:) = aux_var_up%xmol((np-1)*option%nflowspec+1:np*option%nflowspec)
    else
      uxmol(:) = aux_var_dn%xmol((np-1)*option%nflowspec+1:np*option%nflowspec)
    endif
    do ispec=1, option%nflowspec
      fluxm(ispec) = fluxm(ispec) + q*density_ave*uxmol(ispec)
    enddo 
  enddo

  Res(1:option%nflowspec) = fluxm(:)*option%flow_dt

end subroutine MiscibleBCFluxAdv

! ************************************************************************** !
!
! MiscibleBCFluxDiffusion: Computes the  boundary flux terms for the residual (not used)
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine MiscibleBCFluxDiffusion(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
     por_dn,tor_dn,dd_up,perm_dn,Dk_dn, &
     area,dist_gravity,option,vv_darcy,Res)
  use Option_module
  
  implicit none
  
  PetscInt :: ibndtype(:)
  type(Miscible_auxvar_elem_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: vv_darcy(:), area
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec, np
  PetscReal :: fluxm(option%nflowspec),q,density_ave,v_darcy
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi
  
  fluxm = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

! Diffusion term   
  diffdp = por_dn*tor_dn/dd_up*area
  select case(ibndtype(2))
    case(DIRICHLET_BC) 
!       diff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*(aux_var_up%den+aux_var_dn%den)
      do np = 1, option%nphase
!       if (aux_var_up%sat(np)>eps .and. aux_var_dn%sat(np)>eps) then
          diff = diffdp*0.25D0*(aux_var_up%sat(np)+aux_var_dn%sat(np))* &
                    (aux_var_up%den(np)+aux_var_up%den(np))
          do ispec = 1, option%nflowspec
            fluxm(ispec) = fluxm(ispec) + diff*aux_var_dn%diff((np-1)*option%nflowspec+ispec)* &
            (aux_var_up%xmol((np-1)* option%nflowspec+ispec) &
            -aux_var_dn%xmol((np-1)* option%nflowspec+ispec))
          enddo
      enddo
     
  end select

  Res(1:option%nflowspec) = fluxm(:)*option%flow_dt

end subroutine MiscibleBCFluxDiffusion


! ************************************************************************** !
!
! MiscibleFluxAdv: Computes the internal flux terms for the residual (not used)
! author: Chuan Lu
! date: 05/04/10
!
! ************************************************************************** ! 
subroutine MiscibleFluxAdv(aux_var_up,por_up,tor_up,dd_up,perm_up,Dk_up, &
                        aux_var_dn,por_dn,tor_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,vv_darcy,Res)
  use Option_module                              
  
  implicit none
  
  type(Miscible_auxvar_elem_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: vv_darcy(:),area
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
     
  PetscInt :: ispec, np, ind
  PetscReal :: fluxm(option%nflowspec),fluxe,q, v_darcy
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,difff,diffdp,DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
     
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
! diffdp = (por_up*tor_up*por_dn*tor_dn)/(dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  vv_darcy =0.D0 
  
! Flow term
  do np = 1, option%nphase
      upweight= dd_dn/(dd_up+dd_dn)
      if (aux_var_up%sat(np) <eps) then 
        upweight=0.d0
      else if (aux_var_dn%sat(np) <eps) then 
        upweight=1.d0
      endif
      density_ave = upweight*aux_var_up%den(np) + (1.D0-upweight)*aux_var_dn%den(np) 
        
      gravity = (upweight*aux_var_up%den(np) * aux_var_up%avgmw(np) + &
             (1.D0-upweight)*aux_var_dn%den(np) * aux_var_dn%avgmw(np)) &
             * dist_gravity

      dphi = aux_var_up%pres - aux_var_dn%pres &
             - aux_var_up%pc(np) + aux_var_dn%pc(np) &
             + gravity

      v_darcy = 0.D0
      ukvr = 0.D0
      uxmol = 0.D0

      ! note uxmol only contains one phase xmol
      if (dphi >= 0.D0) then
        ukvr = aux_var_up%kvr(np)
        uxmol(:) = aux_var_up%xmol((np-1)*option%nflowspec+1:np*option%nflowspec)
      else
        ukvr = aux_var_dn%kvr(np)
        uxmol(:) = aux_var_dn%xmol((np-1)*option%nflowspec+1:np*option%nflowspec)
      endif

      if (ukvr > floweps) then
        v_darcy = Dq*ukvr*dphi
        vv_darcy(np) = v_darcy
        q = v_darcy*area
        do ispec = 1, option%nflowspec
          fluxm(ispec) = fluxm(ispec) + q*density_ave*uxmol(ispec)
        enddo  
      endif
  end do
     
  Res(1:option%nflowspec) = fluxm(:)*option%flow_dt

end subroutine MiscibleFluxAdv

! ************************************************************************** !
!
! MiscibleFluxDiffusion: Computes the internal flux terms for the residual (not used)
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** ! 
subroutine MiscibleFluxDiffusion(aux_var_up,por_up,tor_up,dd_up,perm_up,Dk_up, &
                        aux_var_dn,por_dn,tor_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,vv_darcy,Res)
  use Option_module                              
  
  implicit none
  
  type(Miscible_auxvar_elem_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: vv_darcy(:),area
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
     
  PetscInt :: ispec, np, ind
  PetscReal :: fluxm(option%nflowspec),fluxe,q, v_darcy
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,difff,diffdp,DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
     

  diffdp = (por_up *tor_up * por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  vv_darcy = 0.D0 
  
! Flow term
  do np = 1, option%nphase
 
! Diffusion term   
! Note : average rule may not be correct  
        difff = diffdp * 0.25D0*(aux_var_up%sat(np) + aux_var_dn%sat(np))* &
             (aux_var_up%den(np) + aux_var_dn%den(np))
        do ispec=1, option%nflowspec
           ind = ispec + (np-1)*option%nflowspec
           fluxm(ispec) = fluxm(ispec) + difff * .5D0 * &
                (aux_var_up%diff(ind) + aux_var_dn%diff(ind))* &
                (aux_var_up%xmol(ind) - aux_var_dn%xmol(ind))
        enddo
  enddo

  Res(1:option%nflowspec) = fluxm(:)*option%flow_dt

end subroutine MiscibleFluxDiffusion
