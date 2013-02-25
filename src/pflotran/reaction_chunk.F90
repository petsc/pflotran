module Reaction_Chunk_module

  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"

  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: RReactChunk, &
            RPack, &
            RUnpack, &
            RTAuxVarChunkCreate


contains

! ************************************************************************** !
!
! RTAuxVarChunkCreate: Creates an RTAuxVar object with correct dimension
! author: Glenn Hammond
! date: 01/31/11
!
! ************************************************************************** !
function RTAuxVarChunkCreate(reaction,option)

  use Option_module

  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  type(react_tran_auxvar_chunk_type), pointer :: RTAuxVarChunkCreate
  
  type(react_tran_auxvar_chunk_type), pointer :: auxvar
  PetscInt :: chunk_size
  PetscInt :: num_threads

  allocate(auxvar)
  
  chunk_size = option%chunk_size
  num_threads = option%num_threads
  
  ! for global auxvar
  allocate(auxvar%den(chunk_size,num_threads,option%nphase))
  allocate(auxvar%temp(chunk_size,num_threads,option%nphase))
  allocate(auxvar%sat(chunk_size,num_threads,option%nphase))
  allocate(auxvar%vol(chunk_size,num_threads))
  allocate(auxvar%por(chunk_size,num_threads))
  
#ifdef CHUAN_CO2
  nullify(auxvar%pres)
  nullify(auxvar%xmass)
  nullify(auxvar%fugacoeff)
!  allocate(auxvar%pres(chunk_size,num_threads,option%nphase))
!  allocate(auxvar%xmass(chunk_size,num_threads,option%nphase))
!  allocate(auxvar%fugacoeff(chunk_size,num_threads,option%nphase))
#endif
  
  allocate(auxvar%pri_molal(chunk_size,num_threads,reaction%naqcomp))
  auxvar%pri_molal = 0.d0
  allocate(auxvar%ln_pri_molal(chunk_size,num_threads,reaction%naqcomp))
  auxvar%ln_pri_molal = 0.d0

  allocate(auxvar%total(chunk_size,num_threads,reaction%naqcomp,option%nphase))
  auxvar%total = 0.d0
  
  allocate(auxvar%dtotal(chunk_size,num_threads,reaction%naqcomp,reaction%naqcomp, &
                          option%nphase))
  auxvar%dtotal = 0.d0
  
  if (reaction%neqcplx > 0) then
    allocate(auxvar%sec_molal(chunk_size,num_threads,reaction%neqcplx))
    auxvar%sec_molal = 0.d0
  else
    nullify(auxvar%sec_molal)
  endif
  
  if (reaction%ngas > 0) then
    allocate(auxvar%gas_molal(chunk_size,num_threads,reaction%ngas))
    auxvar%gas_molal = 0.d0
  else
    nullify(auxvar%gas_molal)
  endif

  if (reaction%neqsorb > 0) then  
    allocate(auxvar%total_sorb_eq(chunk_size,num_threads,reaction%naqcomp))
    auxvar%total_sorb_eq = 0.d0
    if (reaction%kinmr_nrate <= 0) then
      allocate(auxvar%dtotal_sorb_eq(chunk_size,num_threads,reaction%naqcomp,reaction%naqcomp))
      auxvar%dtotal_sorb_eq = 0.d0
    else
      nullify(auxvar%dtotal_sorb_eq)
    endif
  else
    nullify(auxvar%total_sorb_eq)
    nullify(auxvar%dtotal_sorb_eq)
  endif    
  
  if (reaction%neqsrfcplxrxn > 0) then
    allocate(auxvar%eqsrfcplx_conc(chunk_size,num_threads,reaction%neqsrfcplx))
    auxvar%eqsrfcplx_conc = 0.d0
    
    allocate(auxvar%eqsrfcplx_free_site_conc(chunk_size,num_threads,reaction%neqsrfcplxrxn))
    auxvar%eqsrfcplx_free_site_conc = 1.d-9 ! initialize to guess
    
!   allocate(auxvar%eqsurf_site_density(reaction%neqsrfcplxrxn))
!   auxvar%eqsurf_site_density = 0.d0
  else
    nullify(auxvar%eqsrfcplx_conc)
    nullify(auxvar%eqsrfcplx_free_site_conc)
!   nullify(auxvar%eqsurf_site_density)
  endif
  
  if (reaction%nkinmnrl > 0) then
    allocate(auxvar%mnrl_volfrac(chunk_size,num_threads,reaction%nkinmnrl))
    auxvar%mnrl_volfrac = 0.d0
    allocate(auxvar%mnrl_area(chunk_size,num_threads,reaction%nkinmnrl))
    auxvar%mnrl_area = 0.d0
    allocate(auxvar%mnrl_rate(chunk_size,num_threads,reaction%nkinmnrl))
    auxvar%mnrl_rate = 0.d0
  else
    nullify(auxvar%mnrl_volfrac)
    nullify(auxvar%mnrl_area)
    nullify(auxvar%mnrl_rate)
  endif
  
  allocate(auxvar%pri_act_coef(chunk_size,num_threads,reaction%naqcomp))
  auxvar%pri_act_coef = 1.d0
  if (reaction%neqcplx > 0) then
    allocate(auxvar%sec_act_coef(chunk_size,num_threads,reaction%neqcplx))
    auxvar%sec_act_coef = 1.d0
  else
    nullify(auxvar%sec_act_coef)
  endif

! initialize ln activity H2O
  allocate(auxvar%ln_act_h2o(chunk_size,num_threads))
  auxvar%ln_act_h2o = 0.d0
  
  RTAuxVarChunkCreate => auxvar
  
end function RTAuxVarChunkCreate

! ************************************************************************** !
!
! RPack: Packs rt_auxvar_type data into dense chunked arrays
! author: Glenn Hammond
! date: 01/31/11
!
! ************************************************************************** !
subroutine RPack(rt_auxvar,global_auxvar,total,rt_auxvar_chunk,volume, &
                 porosity,ichunk,ithread,reaction)

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(react_tran_auxvar_chunk_type) :: rt_auxvar_chunk
  PetscReal :: volume
  PetscReal :: porosity
  PetscInt :: ichunk
  PetscInt :: ithread
  type(reaction_type) :: reaction
  PetscReal :: total(reaction%naqcomp)
  
  rt_auxvar_chunk%den(ichunk,ithread,:) = global_auxvar%den_kg(:)
  rt_auxvar_chunk%temp(ichunk,ithread,:) = global_auxvar%temp(:)
  rt_auxvar_chunk%sat(ichunk,ithread,:) = global_auxvar%sat(:)
  rt_auxvar_chunk%vol(ichunk,ithread) = volume
  rt_auxvar_chunk%por(ichunk,ithread) = porosity

#ifdef CHUAN_CO2
!  rt_auxvar_chunk%pres(ichunk,ithread,:) = global_auxvar%pres(:)
!  if (associated(global_auxvar%xmass)) &
!    rt_auxvar_chunk%xmass(ichunk,ithread,:) = global_auxvar%xmass(:)
!  if (associated(global_auxvar%fugacoeff)) &
!    rt_auxvar_chunk%fugacoeff(ichunk,ithread,:) = global_auxvar%fugacoeff(:)
#endif
  
  rt_auxvar_chunk%pri_molal(ichunk,ithread,:) = rt_auxvar%pri_molal(:)
  ! Nope
  !rt_auxvar_chunk%total(ichunk,ithread,:,:) = rt_auxvar%total(:,:)
  ! total(:) is the ichunk portion of tran_xx
  rt_auxvar_chunk%total(ichunk,ithread,:,ONE_INTEGER) = total(:)
  
  rt_auxvar_chunk%dtotal(ichunk,ithread,:,:,:) = rt_auxvar%aqueous%dtotal(:,:,:)
  rt_auxvar_chunk%total_sorb_eq(ichunk,ithread,:) = rt_auxvar%total_sorb_eq(:)
  rt_auxvar_chunk%dtotal_sorb_eq(ichunk,ithread,:,:) = rt_auxvar%dtotal_sorb_eq(:,:)
  rt_auxvar_chunk%sec_molal(ichunk,ithread,:) = rt_auxvar%sec_molal
  rt_auxvar_chunk%gas_molal(ichunk,ithread,:) = rt_auxvar%gas_molal
  rt_auxvar_chunk%eqsrfcplx_conc(ichunk,ithread,:) = rt_auxvar%eqsrfcplx_conc(:)
  rt_auxvar_chunk%eqsrfcplx_free_site_conc(ichunk,ithread,:) = &
    rt_auxvar%eqsrfcplx_free_site_conc(:)
  rt_auxvar_chunk%mnrl_volfrac(ichunk,ithread,:) = rt_auxvar%mnrl_volfrac(:)
  rt_auxvar_chunk%mnrl_area(ichunk,ithread,:) = rt_auxvar%mnrl_area(:)
  rt_auxvar_chunk%mnrl_rate(ichunk,ithread,:) = rt_auxvar%mnrl_rate(:)
  rt_auxvar_chunk%pri_act_coef(ichunk,ithread,:) = rt_auxvar%pri_act_coef(:)
  rt_auxvar_chunk%sec_act_coef(ichunk,ithread,:) = rt_auxvar%sec_act_coef(:)
  rt_auxvar_chunk%ln_act_h2o(ichunk,ithread) = rt_auxvar%ln_act_h2o

end subroutine RPack

! ************************************************************************** !
!
! RUnpack: Unpacks dense chunked arrays into rt_auxvar_type data 
! author: Glenn Hammond
! date: 01/31/11
!
! ************************************************************************** !
subroutine RUnpack(rt_auxvar,free_ion,rt_auxvar_chunk,ichunk,ithread,reaction)

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(react_tran_auxvar_chunk_type) :: rt_auxvar_chunk
  PetscInt :: ichunk
  PetscInt :: ithread
  type(reaction_type) :: reaction
  PetscReal :: free_ion(reaction%naqcomp)
  
  rt_auxvar%pri_molal(:) = rt_auxvar_chunk%pri_molal(ichunk,ithread,:)
  free_ion(:) = rt_auxvar%pri_molal(:)

  rt_auxvar%total(:,:) = rt_auxvar_chunk%total(ichunk,ithread,:,:)
  rt_auxvar%aqueous%dtotal(:,:,:) = rt_auxvar_chunk%dtotal(ichunk,ithread,:,:,:)
  rt_auxvar%total_sorb_eq(:) = rt_auxvar_chunk%total_sorb_eq(ichunk,ithread,:)
  rt_auxvar%dtotal_sorb_eq(:,:) = rt_auxvar_chunk%dtotal_sorb_eq(ichunk,ithread,:,:)
  rt_auxvar%sec_molal = rt_auxvar_chunk%sec_molal(ichunk,ithread,:)
  rt_auxvar%gas_molal = rt_auxvar_chunk%gas_molal(ichunk,ithread,:)
  rt_auxvar%eqsrfcplx_conc(:) = rt_auxvar_chunk%eqsrfcplx_conc(ichunk,ithread,:)
  rt_auxvar%eqsrfcplx_free_site_conc(:) = &
    rt_auxvar_chunk%eqsrfcplx_free_site_conc(ichunk,ithread,:)
!geh: since volfrac and area do not change, shouldn't need to update
!  rt_auxvar%mnrl_volfrac(:) = rt_auxvar_chunk%mnrl_volfrac(ichunk,ithread,:)
!  rt_auxvar%mnrl_area(:) = rt_auxvar_chunk%mnrl_area(ichunk,ithread,:)
  rt_auxvar%mnrl_rate(:) = rt_auxvar_chunk%mnrl_rate(ichunk,ithread,:)
  rt_auxvar%pri_act_coef(:) = rt_auxvar_chunk%pri_act_coef(ichunk,ithread,:)
  rt_auxvar%sec_act_coef(:) = rt_auxvar_chunk%sec_act_coef(ichunk,ithread,:)
  rt_auxvar%ln_act_h2o = rt_auxvar_chunk%ln_act_h2o(ichunk,ithread)

end subroutine RUnpack

! ************************************************************************** !
!
! RReactChunk: Solves reaction portion of operator splitting using 
!              Newton-Raphson
! author: Glenn Hammond
! date: 11/10/10
!
! ************************************************************************** !
subroutine RReactChunk(auxvar,num_iterations_,reaction,vector_length, &
                       ithread,option)

  use Option_module
  use Utility_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(option_type) :: option
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  PetscInt :: num_iterations_(option%chunk_size,option%num_threads)
  
  PetscReal :: residual(vector_length,option%num_threads,reaction%ncomp)
  PetscReal :: res(vector_length,option%num_threads,reaction%ncomp)
  PetscReal :: J(vector_length,option%num_threads,reaction%ncomp,reaction%ncomp)
  PetscReal :: one_over_dt
  PetscReal :: prev_molal(vector_length,option%num_threads,reaction%ncomp)
  PetscReal :: update(vector_length,option%num_threads,reaction%ncomp)
  PetscReal :: maximum_relative_change(vector_length,option%num_threads)
  PetscReal :: accumulation_coef
  PetscReal :: fixed_accum(vector_length,option%num_threads,reaction%ncomp)
  PetscInt :: num_iterations(vector_length,option%num_threads)
  PetscInt :: icomp

  ! for inlined RSolve() routine
  PetscReal :: norm
  PetscInt :: indices(vector_length,option%num_threads,reaction%ncomp)
  PetscReal :: rhs(vector_length,option%num_threads,reaction%ncomp)
  
  PetscInt :: ichunk
  PetscInt :: offset
  PetscInt :: d(vector_length,option%num_threads)
  
  PetscInt, parameter :: iphase = 1

  one_over_dt = 1.d0/option%tran_dt

  num_iterations(1:vector_length,ithread) = 0

  ! calculate fixed portion of accumulation term
  ! fixed_accum is overwritten in RTAccumulation
  ! Since RTAccumulation uses auxvar%total, we must overwrite the 
  ! auxvar total variables
  ! aqueous
  ! NOTE: This is now performed in RPack()
!  do ichunk = 1, vector_length
!    offset = (ichunk-1)*reaction%naqcomp
!    do icomp = 1, reaction%naqcomp
!      auxvar%total(ichunk,ithread,icomp,iphase) = total(offset+icomp)
!    enddo
!  enddo
  
  ! still need code to overwrite other phases
  call RTAccumulationChunk(auxvar,vector_length,ithread,reaction,option,fixed_accum)
  if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
    call RAccumulationSorbChunk(auxvar,vector_length,ithread,reaction,option,fixed_accum)  
  endif

  ! now update activity coefficients
  if (reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
    call RActivityCoefficientsChunk(auxvar,vector_length,ithread,reaction,option)
  endif
  
  do
  
    num_iterations(:,ithread) = num_iterations(:,ithread) + 1

    auxvar%ln_pri_molal(:,ithread,:) = log(auxvar%pri_molal(:,ithread,:))

    if (reaction%act_coef_update_frequency == ACT_COEF_FREQUENCY_NEWTON_ITER) then
      call RActivityCoefficientsChunk(auxvar,vector_length,ithread,reaction,option)
    endif
    call RTAuxVarComputeChunk(auxvar,vector_length,ithread,reaction,option)
    
    ! Accumulation
    ! residual is overwritten in RTAccumulation()
    call RTAccumulationChunk(auxvar,vector_length,ithread,reaction,option,residual)
                        
    do ichunk = 1, vector_length
      residual(ichunk,ithread,:) = &
        residual(ichunk,ithread,:) - &
        fixed_accum(ichunk,ithread,:)
    enddo

    ! J is overwritten in RTAccumulationDerivative()
    call RTAccumulationDerivativeChunk(auxvar,vector_length,ithread,reaction,option,J)

    if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
      call RAccumulationSorbChunk(auxvar,vector_length,ithread,reaction, &
                                  option,residual)
      call RAccumulationSorbDerivChunk(auxvar,vector_length,ithread,reaction,option,J)
    endif

                         ! derivative
    call RReactionChunk(residual,J,PETSC_TRUE,auxvar,vector_length,ithread,reaction,option)
    
    ! Manual inlining start
    ! call RSolveChunk(residual,J,auxvar%pri_molal,update,reaction%ncomp,vector_length)

    ! scale Jacobian
    do ichunk = 1, vector_length
    
      do icomp = 1, reaction%ncomp
        norm = max(1.d0,maxval(abs(J(ichunk,ithread,icomp,:))))
        norm = 1.d0/norm
        rhs(ichunk,ithread,icomp) = residual(ichunk,ithread,icomp)*norm
        J(ichunk,ithread,icomp,:) = J(ichunk,ithread,icomp,:)*norm
      enddo
      
      ! for derivatives with respect to ln conc
      do icomp = 1, reaction%ncomp
        J(ichunk,ithread,:,icomp) = J(ichunk,ithread,:,icomp) * &
                                  auxvar%pri_molal(ichunk,ithread,icomp)
      enddo

    enddo
   
    call ludcmp_chunk(J,reaction%ncomp,indices,d,vector_length,ithread, &
                      option%num_threads)
    call lubksb_chunk(J,reaction%ncomp,indices,rhs,vector_length,ithread, &
                      option%num_threads)
    
    update(:,ithread,:) = rhs(:,ithread,:)

    ! Manual inlining end

    do ichunk = 1, vector_length

      update(ichunk,ithread,:) = &
        dsign(1.d0,update(ichunk,ithread,:)) * &
        min(dabs(update(ichunk,ithread,:)),5.d0)

      prev_molal(ichunk,ithread,:) = auxvar%pri_molal(ichunk,ithread,:)
      auxvar%pri_molal(ichunk,ithread,:) = &
        auxvar%pri_molal(ichunk,ithread,:) * &
        exp(-update(ichunk,ithread,:))    
  
      maximum_relative_change(ichunk,ithread) = &
        maxval(abs((auxvar%pri_molal(ichunk,ithread,:) - &
                    prev_molal(ichunk,ithread,:))/ &
               prev_molal(ichunk,ithread,:)))
    enddo
  
!    if (maximum_relative_change < reaction%reaction_tolerance) exit
    if (maxval(maximum_relative_change(:,ithread)) < &
        reaction%max_relative_change_tolerance) exit
    
  enddo

  ! one last update
  call RTAuxVarComputeChunk(auxvar,vector_length,ithread,reaction,option)

  do ichunk = 1, vector_length
    num_iterations_(ichunk,ithread) = &
      num_iterations(ichunk,ithread)
  enddo
  
end subroutine RReactChunk
      
! ************************************************************************** !
!
! RReactionChunk: Computes reactions
! author: Glenn Hammond
! date: 09/30/08
!
! ************************************************************************** !
subroutine RReactionChunk(Res,Jac,derivative,auxvar,vector_length,ithread,reaction,option)

  use Option_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(option_type) :: option
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  PetscBool :: derivative
  PetscReal :: Res(vector_length,option%num_threads,reaction%ncomp)
  PetscReal :: Jac(vector_length,option%num_threads,reaction%ncomp,reaction%ncomp)

  if (reaction%nkinmnrl > 0) then
    call RKineticMineralChunk(Res,Jac,derivative,auxvar,vector_length,ithread,reaction,option)
  endif
  
  if (reaction%kinmr_nrate > 0) then
    option%io_buffer = "RMultiRateSorption not set up for chunking"
    call printErrMsg(option)
    !call RMultiRateSorption(Res,Jac,derivative,auxvar,reaction,option)
  endif
  
  if (reaction%nkinsrfcplxrxn > 0) then
    option%io_buffer = "RKineticSurfCplx not set up for chunking"
    call printErrMsg(option)
    !call RKineticSurfCplx(Res,Jac,derivative,auxvar,reaction,option)
  endif
  
  if (reaction%ngeneral_rxn > 0) then
    option%io_buffer = "RGeneral not set up for chunking"
    call printErrMsg(option)
    !call RGeneral(Res,Jac,derivative,auxvar,reaction,option)
  endif
  
  ! add new reactions here

end subroutine RReactionChunk

! ************************************************************************** !
!
! RActivityCoefficientsChunk: Computes the ionic strength and activity 
!                             coefficients
! author: Glenn Hammond
! date: 11/10/10
!
! ************************************************************************** !
subroutine RActivityCoefficientsChunk(auxvar,vector_length,ithread,reaction,option)

  use Option_module
  
  implicit none

  type(option_type) :: option
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  type(reaction_type) :: reaction
  
  PetscInt :: icplx, icomp, it, j, jcomp, ncomp
  PetscReal :: I, sqrt_I, II, sqrt_II, f, fpri, didi, dcdi, den, dgamdi, &
    lnQK, sum, sum_pri_molal, sum_sec_molal
  PetscReal :: sum_molality
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)

  PetscInt :: ichunk

  do ichunk = 1, vector_length

  if (reaction%use_activity_h2o) then
    sum_pri_molal = 0.d0
    do j = 1, reaction%naqcomp
      sum_pri_molal = sum_pri_molal + auxvar%pri_molal(ichunk,ithread,j)
    enddo
  endif

  if (reaction%act_coef_update_algorithm == ACT_COEF_ALGORITHM_NEWTON) then

    ln_conc = log(auxvar%pri_molal(ichunk,ithread,:))
    ln_act = ln_conc+log(auxvar%pri_act_coef(ichunk,ithread,:))
  
#if 0
    if (.not.option%use_isothermal) then
      call ReactionInterpolateLogK(reaction%eqcplx_logKcoef,reaction%eqcplx_logK, &
                                   auxvar%temp(ichunk,ithread,1),reaction%neqcplx)
    endif
#endif  
  
  ! compute primary species contribution to ionic strength
    fpri = 0.d0
    sum_molality = 0.d0
    do j = 1, reaction%naqcomp
      fpri = fpri + auxvar%pri_molal(ichunk,ithread,j)*reaction%primary_spec_Z(j)* &
                                               reaction%primary_spec_Z(j)
    enddo
  
    it = 0
    II = 0
    do
      it = it + 1
      
      if (it > 50) then
        print *,' too many iterations in computing activity coefficients-stop',it,f,I
        stop
      endif
    
  ! add secondary species contribution to ionic strength
      I = fpri
      do icplx = 1, reaction%neqcplx ! for each secondary species
        I = I + auxvar%sec_molal(ichunk,ithread,icplx)*reaction%eqcplx_Z(icplx)* &
                                                  reaction%eqcplx_Z(icplx)
      enddo
      I = 0.5d0*I
      f = I
    
      if (abs(I-II) < 1.d-6*I) exit
    
      if (reaction%neqcplx > 0) then
        didi = 0.d0
        sqrt_I = sqrt(I)
        do icplx = 1, reaction%neqcplx
          if (abs(reaction%eqcplx_Z(icplx)) > 0.d0) then
            sum = 0.5d0*reaction%debyeA*reaction%eqcplx_Z(icplx)* &
            reaction%eqcplx_Z(icplx) &
            /(sqrt_I*(1.d0+reaction%debyeB*reaction%eqcplx_a0(icplx)*sqrt_I)**2) &
            -reaction%debyeBdot
            ncomp = reaction%eqcplxspecid(0,icplx)
            do jcomp = 1, ncomp
              j = reaction%eqcplxspecid(jcomp,icplx)
              if (abs(reaction%primary_spec_Z(j)) > 0.d0) then
                dgamdi = -0.5d0*reaction%debyeA*reaction%primary_spec_Z(j)**2/(sqrt_I* &
                (1.d0+reaction%debyeB*reaction%primary_spec_a0(j)*sqrt_I)**2)+ &
                reaction%debyeBdot 
                sum = sum + reaction%eqcplxstoich(jcomp,icplx)*dgamdi
              endif
            enddo
            dcdi = auxvar%sec_molal(ichunk,ithread,icplx)*LOG_TO_LN*sum
            didi = didi+0.5d0*reaction%eqcplx_Z(icplx)*reaction%eqcplx_Z(icplx)*dcdi
          endif
        enddo
        den = 1.d0-didi
        if (abs(den) > 0.d0) then
          II = (f-I*didi)/den
        else
          II = f
        endif
      else
        II = f
      endif
    
      if (II < 0.d0) then
        write(option%io_buffer,*) 'ionic strength negative! it =',it, &
          ' I= ',I,II,den,didi,dcdi,sum
        call printErrMsg(option)        
      endif
    
  ! compute activity coefficients
  ! primary species
      I = II
      sqrt_I = sqrt(I)
      do icomp = 1, reaction%naqcomp
        if (abs(reaction%primary_spec_Z(icomp)) > 0.d0) then
          auxvar%pri_act_coef(ichunk,ithread,icomp) = exp((-reaction%primary_spec_Z(icomp)* &
                                        reaction%primary_spec_Z(icomp)* &
                                        sqrt_I*reaction%debyeA/ &
                                        (1.d0+reaction%primary_spec_a0(icomp)* &
                                        reaction%debyeB*sqrt_I)+ &
                                        reaction%debyeBdot*I)* &
                                        LOG_TO_LN)
        else
          auxvar%pri_act_coef(ichunk,ithread,icomp) = 1.d0
        endif
      enddo
                
  ! secondary species
      sum_sec_molal = 0.d0
      do icplx = 1, reaction%neqcplx
        if (abs(reaction%eqcplx_Z(icplx)) > 0.d0) then
          auxvar%sec_act_coef(ichunk,ithread,icplx) = exp((-reaction%eqcplx_Z(icplx)* &
                                        reaction%eqcplx_Z(icplx)* &
                                        sqrt_I*reaction%debyeA/ &
                                        (1.d0+reaction%eqcplx_a0(icplx)* &
                                        reaction%debyeB*sqrt_I)+ &
                                        reaction%debyeBdot*I)* &
                                        LOG_TO_LN)
        else
          auxvar%sec_act_coef(ichunk,ithread,icplx) = 1.d0
        endif
    
    ! compute secondary species concentration
        lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

    ! activity of water
        if (reaction%eqcplxh2oid(icplx) > 0) then
          lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*auxvar%ln_act_h2o(ichunk,ithread)
        endif

        ncomp = reaction%eqcplxspecid(0,icplx)
        do jcomp = 1, ncomp
          icomp = reaction%eqcplxspecid(jcomp,icplx)
          lnQK = lnQK + reaction%eqcplxstoich(jcomp,icplx)*ln_act(icomp)
        enddo
        auxvar%sec_molal(ichunk,ithread,icplx) = exp(lnQK)/auxvar%sec_act_coef(ichunk,ithread,icplx)
        sum_sec_molal = sum_sec_molal + auxvar%sec_molal(ichunk,ithread,icplx)
      
      enddo
      
      if (reaction%use_activity_h2o) then
        auxvar%ln_act_h2o(ichunk,ithread) = 1.d0-0.017d0*(sum_pri_molal+sum_sec_molal)
        if (auxvar%ln_act_h2o(ichunk,ithread) > 0.d0) then
          auxvar%ln_act_h2o(ichunk,ithread) = log(auxvar%ln_act_h2o(ichunk,ithread))
        else
          auxvar%ln_act_h2o(ichunk,ithread) = 0.d0
          write(option%io_buffer,*) 'activity of H2O negative! ln act H2O =', &
            auxvar%ln_act_h2o(ichunk,ithread)
          call printMsg(option)
        endif
      endif

    enddo
  
  else
  
  ! compute ionic strength
  ! primary species
    I = 0.d0
    do icomp = 1, reaction%naqcomp
      I = I + auxvar%pri_molal(ichunk,ithread,icomp)*reaction%primary_spec_Z(icomp)* &
                                                reaction%primary_spec_Z(icomp)
    enddo
  
  ! secondary species
    do icplx = 1, reaction%neqcplx ! for each secondary species
      I = I + auxvar%sec_molal(ichunk,ithread,icplx)*reaction%eqcplx_Z(icplx)* &
                                             reaction%eqcplx_Z(icplx)
    enddo
    I = 0.5d0*I
    sqrt_I = sqrt(I)
  
  ! compute activity coefficients
  ! primary species
    do icomp = 1, reaction%naqcomp
      if (abs(reaction%primary_spec_Z(icomp)) > 1.d-10) then
        auxvar%pri_act_coef(ichunk,ithread,icomp) = exp((-reaction%primary_spec_Z(icomp)* &
                                      reaction%primary_spec_Z(icomp)* &
                                      sqrt_I*reaction%debyeA/ &
                                      (1.d0+reaction%primary_spec_a0(icomp)* &
                                      reaction%debyeB*sqrt_I)+ &
                                      reaction%debyeBdot*I)* &
                                      LOG_TO_LN)
      else
        auxvar%pri_act_coef(ichunk,ithread,icomp) = 1.d0
      endif
    enddo
                
  ! secondary species
    sum_sec_molal = 0.d0
    do icplx = 1, reaction%neqcplx
      if (dabs(reaction%eqcplx_Z(icplx)) > 1.d-10) then
        auxvar%sec_act_coef(ichunk,ithread,icplx) = exp((-reaction%eqcplx_Z(icplx)* &
                                      reaction%eqcplx_Z(icplx)* &
                                      sqrt_I*reaction%debyeA/ &
                                      (1.d0+reaction%eqcplx_a0(icplx)* &
                                      reaction%debyeB*sqrt_I)+ &
                                      reaction%debyeBdot*I)* &
                                      LOG_TO_LN)
      else
        auxvar%sec_act_coef(ichunk,ithread,icplx) = 1.d0
      endif
      sum_sec_molal = sum_sec_molal + auxvar%sec_molal(ichunk,ithread,icplx)
    enddo
    
    if (reaction%use_activity_h2o) then
      auxvar%ln_act_h2o(ichunk,ithread) = 1.d0-0.017d0*(sum_pri_molal+sum_sec_molal)
      if (auxvar%ln_act_h2o(ichunk,ithread) > 0.d0) then
        auxvar%ln_act_h2o(ichunk,ithread) = log(auxvar%ln_act_h2o(ichunk,ithread))
      else
        auxvar%ln_act_h2o(ichunk,ithread) = 0.d0
      endif
    endif
  endif
  
  enddo ! over chunks
  
end subroutine RActivityCoefficientsChunk

! ************************************************************************** !
!
! RTotalChunk: Computes the total component concentrations and derivative with
!              respect to free-ion
! author: Glenn Hammond
! date: 11/10/10
!
! ************************************************************************** !
subroutine RTotalChunk(auxvar,vector_length,ithread,reaction,option)

  use Option_module
  use co2eos_module, only: Henry_duan_sun
  use Water_EOS_module
  
  type(option_type) :: option
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  type(reaction_type) :: reaction
  
  PetscInt :: i, j, icplx, icomp, jcomp, iphase, ncomp, ieqgas
  PetscErrorCode :: ierr
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: lnQK, tempreal
  PetscReal :: den_kg_per_L, xmass
  PetscReal :: pressure, temperature, xphico2, muco2, den, m_na, m_cl
  
  PetscInt :: ichunk

#ifdef CHUAN_CO2  
  PetscReal :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,&
               yco2,pco2,sat_pressure,lngamco2
#endif
  
  do ichunk = 1, vector_length
  
#ifdef CHUAN_CO2  
  auxvar%total(ichunk,ithread,:,:) = 0.d0 !debugging 
#endif
  
  iphase = 1           
!  den_kg_per_L = global_auxvar%den_kg(iphase)*1.d-3              
  xmass = 1.d0
#ifdef CHUAN_CO2  
  if (associated(auxvar%xmass)) xmass = auxvar%xmass(ichunk,ithread,iphase)
#endif  
  den_kg_per_L = auxvar%den(ichunk,ithread,iphase)*xmass*1.d-3

  ln_conc = log(auxvar%pri_molal(ichunk,ithread,:))
  ln_act = ln_conc+log(auxvar%pri_act_coef(ichunk,ithread,:))
  auxvar%total(ichunk,ithread,:,iphase) = auxvar%pri_molal(ichunk,ithread,:)
  ! initialize derivatives
  auxvar%dtotal(ichunk,ithread,:,:,:) = 0.d0
  do icomp = 1, reaction%naqcomp
    auxvar%dtotal(ichunk,ithread,icomp,icomp,iphase) = 1.d0
  enddo
  
#if 0
  if (.not.option%use_isothermal .and. reaction%neqcplx > 0) then
    call ReactionInterpolateLogK(reaction%eqcplx_logKcoef,reaction%eqcplx_logK, &
                                 auxvar%temp(ichunk,ithread,iphase),reaction%neqcplx)
  endif
#endif  
  
  do icplx = 1, reaction%neqcplx ! for each secondary species
    ! compute secondary species concentration
    lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

    ! activity of water
    if (reaction%eqcplxh2oid(icplx) > 0) then
      lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*auxvar%ln_act_h2o(ichunk,ithread)
    endif

    ncomp = reaction%eqcplxspecid(0,icplx)
    do i = 1, ncomp
      icomp = reaction%eqcplxspecid(i,icplx)
      lnQK = lnQK + reaction%eqcplxstoich(i,icplx)*ln_act(icomp)
    enddo
    auxvar%sec_molal(ichunk,ithread,icplx) = exp(lnQK)/auxvar%sec_act_coef(ichunk,ithread,icplx)
  
    ! add contribution to primary totals
    ! units of total = mol/L
    do i = 1, ncomp
      icomp = reaction%eqcplxspecid(i,icplx)
      auxvar%total(ichunk,ithread,icomp,iphase) = auxvar%total(ichunk,ithread,icomp,iphase) + &
                                      reaction%eqcplxstoich(i,icplx)* &
                                      auxvar%sec_molal(ichunk,ithread,icplx)
    enddo
    
    ! add contribution to derivatives of total with respect to free
    ! bear in mind that the water density portion is scaled below
    do j = 1, ncomp
      jcomp = reaction%eqcplxspecid(j,icplx)
      tempreal = reaction%eqcplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                 auxvar%sec_act_coef(ichunk,ithread,icplx)
      do i = 1, ncomp
        icomp = reaction%eqcplxspecid(i,icplx)
        auxvar%dtotal(ichunk,ithread,icomp,jcomp,iphase) = auxvar%dtotal(ichunk,ithread,icomp,jcomp,iphase) + &
                                                   reaction%eqcplxstoich(i,icplx)*tempreal
      enddo
    enddo
  enddo

  ! convert molality -> molarity
  ! unit of total = mol/L water
  auxvar%total(ichunk,ithread,:,iphase) = auxvar%total(ichunk,ithread,:,iphase)*den_kg_per_L
  
  ! units of dtotal = kg water/L water
  auxvar%dtotal(ichunk,ithread,:,:,:) = auxvar%dtotal(ichunk,ithread,:,:,:)*den_kg_per_L

  enddo ! loop over chunks

 !*********** Add SC phase contribution ***************************  

#ifdef CHUAN_CO2

  iphase = 2           
#if 0
  if (.not.option%use_isothermal .and. reaction%ngas > 0) then
    call ReactionInterpolateLogK(reaction%eqgas_logKcoef,reaction%eqgas_logK, &
                                 auxvar%temp(ichunk,ithread,1),reaction%ngas)
  endif
#endif  

  if (iphase > option%nphase) return 

  do ichunk = 1, vector_length

  auxvar%total(ichunk,ithread,:,iphase) = 0D0
  auxvar%dtotal(ichunk,ithread,:,:,iphase)=0D0
!  do icomp = 1, reaction%naqcomp
!    auxvar%dtotal(icomp,icomp,iphase) = 1.d0
!  enddo
    
!  den_kg_per_L = global_auxvar%den_kg(iphase)*1.d-3     
  if (auxvar%sat(ichunk,ithread,iphase)>1D-20) then
    do ieqgas = 1, reaction%ngas ! all gas phase species are secondary
   
      pressure = auxvar%pres(ichunk,ithread,2)
      temperature = auxvar%temp(ichunk,ithread,1)
      xphico2 = auxvar%fugacoeff(ichunk,ithread,1)
      den = auxvar%den(ichunk,ithread,2)
 
      call PSAT(temperature, sat_pressure, ierr)
      pco2 = pressure - sat_pressure
!     call co2_span_wagner(pressure*1.D-6,temperature+273.15D0,dg,dddt,dddp,fg, &
!              dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
!
!            fg = fg*1D6
!            xphico2 = fg / pco2
!            global_auxvar%fugacoeff(1) = xphico2


      if (abs(reaction%species_idx%co2_gas_id) == ieqgas ) then
!          call Henry_duan_sun_0NaCl(pco2*1D-5, temperature, henry)
        if (reaction%species_idx%na_ion_id /= 0 .and. reaction%species_idx%cl_ion_id /= 0) then
          m_na = auxvar%pri_molal(ichunk,ithread,reaction%species_idx%na_ion_id)
          m_cl = auxvar%pri_molal(ichunk,ithread,reaction%species_idx%cl_ion_id)
          call Henry_duan_sun(temperature,pressure*1D-5,muco2,xphico2, &
                lngamco2,m_na,m_cl,sat_pressure*1D-5)
        else
          call Henry_duan_sun(temperature,pressure*1D-5,muco2,xphico2, &
                lngamco2,option%m_nacl,option%m_nacl,sat_pressure*1D-5)
        endif
        !lnQk = - log(muco2) 
        lnQk = - log(muco2)-lngamco2
           
      else   
        lnQK = -reaction%eqgas_logK(ieqgas)*LOG_TO_LN
      endif 
          
      if (reaction%eqgash2oid(ieqgas) > 0) then
        lnQK = lnQK + reaction%eqgash2ostoich(ieqgas)*auxvar%ln_act_h2o(ichunk,ithread)
!       print *,'Ttotal', reaction%eqgash2ostoich(ieqgas), auxvar%ln_act_h2o
      endif
   
   ! contribute to %total          
   !     do i = 1, ncomp
   ! removed loop over species, suppose only one primary species is related
      icomp = reaction%eqgasspecid(1,ieqgas)
      pressure =pressure *1D-5
        
      auxvar%gas_molal(ichunk,ithread,ieqgas) = &
          exp(lnQK+lngamco2)*auxvar%pri_molal(ichunk,ithread,icomp)&
!          auxvar%pri_act_coef(ichunk,ithread,icomp)*exp(lnQK)*auxvar%pri_molal(ichunk,ithread,icomp)&
          /pressure /xphico2* den
      auxvar%total(ichunk,ithread,icomp,iphase) = auxvar%total(ichunk,ithread,icomp,iphase) + &
                                        reaction%eqgasstoich(1,ieqgas)* &
                                        auxvar%gas_molal(ichunk,ithread,ieqgas)
!       print *,'Ttotal',pressure, temperature, xphico2, den, lnQk,auxvar%pri_molal(ichunk,ithread,icomp),&
!        global_auxvar%sat(2),auxvar%gas_molal(ichunk,ithread,ieqgas)
   !     if (auxvar%total(ichunk,ithread,icomp,iphase) > den)auxvar%total(ichunk,ithread,icomp,iphase) = den* .99D0
   !     enddo

   ! contribute to %dtotal
   !      tempreal = exp(lnQK+lngamco2)/pressure /xphico2* den 
      tempreal = auxvar%pri_act_coef(ichunk,ithread,icomp)*exp(lnQK)/pressure /xphico2* den 
      auxvar%dtotal(ichunk,ithread,icomp,icomp,iphase) = &
        auxvar%dtotal(ichunk,ithread,icomp,icomp,iphase) + &
        reaction%eqgasstoich(1,ieqgas)*tempreal
    
    enddo
   ! auxvar%total(ichunk,ithread,:,iphase) = auxvar%total(ichunk,ithread,:,iphase)!*den_kg_per_L
    ! units of dtotal = kg water/L water
   ! auxvar%dtotal(ichunk,ithread,:, :,iphase) = auxvar%dtotal(ichunk,ithread,:,:,iphase)!*den_kg_per_L
  endif   
 
  enddo
  
#endif  
end subroutine RTotalChunk

! ************************************************************************** !
!
! RTotalSorbChunk: Computes the total sorbed component concentrations and 
!                  derivative with respect to free-ion
! author: Glenn Hammond
! date: 11/11/10
!
! ************************************************************************** !
subroutine RTotalSorbChunk(auxvar,vector_length,ithread,reaction,option)

  use Option_module
  
  type(option_type) :: option
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  type(reaction_type) :: reaction
  
  PetscInt :: ichunk
  
  ! initialize total sorbed concentrations and derivatives
  if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
    do ichunk = 1, vector_length
      auxvar%total_sorb_eq(ichunk,ithread,:) = 0.d0
      auxvar%dtotal_sorb_eq(ichunk,ithread,:,:) = 0.d0  
    enddo
  endif

  if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
    call RTotalSorbEqSurfCplxChunk(auxvar,vector_length,ithread,reaction,option)
  endif
  
  if (reaction%neqionxrxn > 0) then
    option%io_buffer = "RTotalSorbEqIonx not set up for chunking"
    call printErrMsg(option)
    !call RTotalSorbEqIonx(auxvar,vector_length,ithread,reaction,option)
  endif
  
  if (reaction%neqkdrxn > 0) then
    option%io_buffer = "RTotalSorbKD not set up for chunking"
    call printErrMsg(option)
    !call RTotalSorbKD(auxvar,vector_length,ithread,reaction,option)
  endif
  
end subroutine RTotalSorbChunk

! ************************************************************************** !
!
! RTotalSorbEqSurfCplxChunk: Computes the total sorbed component 
!                       concentrations and 
!                       derivative with respect to free-ion for equilibrium 
!                       surface complexation
! author: Glenn Hammond
! date: 11/11/10
!
! ************************************************************************** !
subroutine RTotalSorbEqSurfCplxChunk(auxvar,vector_length,ithread,reaction,option)

  use Option_module
  
  type(option_type) :: option
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  type(reaction_type) :: reaction
  
  PetscInt :: i, j, k, icplx, icomp, jcomp, ncomp, ncplx
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: srfcplx_conc(reaction%neqsrfcplx)
  PetscReal :: dSx_dmi(reaction%naqcomp)
  PetscReal :: nui_Si_over_Sx
  PetscReal :: free_site_conc
  PetscReal :: ln_free_site
  PetscReal :: lnQK, tempreal, tempreal1, tempreal2, total
  PetscInt :: irxn
  PetscInt, parameter :: iphase = 1
  PetscReal, parameter :: tol = 1.d-12
  PetscBool :: one_more
  PetscReal :: res, dres_dfree_site, dfree_site_conc
  PetscReal :: site_density(2)
  PetscReal :: mobile_fraction
  PetscInt :: num_types_of_sites
  PetscInt :: isite

  PetscInt :: ichunk

  do ichunk = 1, vector_length

  ln_conc = log(auxvar%pri_molal(ichunk,ithread,:))
  ln_act = ln_conc+log(auxvar%pri_act_coef(ichunk,ithread,:))

#if 0
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%eqsrfcplx_logKcoef,reaction%eqsrfcplx_logK, &
                               auxvar%temp(ichunk,ithread,iphase),reaction%neqsrfcplx)
  endif
#endif  

  ! Surface Complexation
  do irxn = 1, reaction%neqsrfcplxrxn
  
    ncplx = reaction%eqsrfcplx_rxn_to_complex(0,irxn)
    
    free_site_conc = auxvar%eqsrfcplx_free_site_conc(ichunk,ithread,irxn)

    select case(reaction%eqsrfcplx_rxn_surf_type(irxn))
      case(MINERAL_SURFACE)
        site_density(1) = reaction%eqsrfcplx_rxn_site_density(irxn)
!        site_density = reaction%eqsrfcplx_rxn_site_density(irxn)* &
!                       auxvar%mnrl_volfrac(ichunk,ithread,reaction%eqsrfcplx_rxn_to_surf(irxn))
        num_types_of_sites = 1
      case(COLLOID_SURFACE)
      case(NULL_SURFACE)
        site_density(1) = reaction%eqsrfcplx_rxn_site_density(irxn)
        num_types_of_sites = 1
    end select
    
    do isite=1, num_types_of_sites
      ! isite == 1 - immobile (colloids, minerals, etc.)
      ! isite == 2 - mobile (colloids)
    
      if (site_density(isite) < 1.d-40) cycle
    
      ! get a pointer to the first complex (there will always be at least 1)
      ! in order to grab free site conc
      one_more = PETSC_FALSE
      do

        total = free_site_conc
        ln_free_site = log(free_site_conc)
        do j = 1, ncplx
          icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
          ! compute secondary species concentration
          lnQK = -reaction%eqsrfcplx_logK(icplx)*LOG_TO_LN

          ! activity of water
          if (reaction%eqsrfcplxh2oid(icplx) > 0) then
            lnQK = lnQK + reaction%eqsrfcplxh2ostoich(icplx)*auxvar%ln_act_h2o(ichunk,ithread)
          endif

          lnQK = lnQK + reaction%eqsrfcplx_free_site_stoich(icplx)* &
                        ln_free_site
        
          ncomp = reaction%eqsrfcplxspecid(0,icplx)
          do i = 1, ncomp
            icomp = reaction%eqsrfcplxspecid(i,icplx)
            lnQK = lnQK + reaction%eqsrfcplxstoich(i,icplx)*ln_act(icomp)
          enddo
          srfcplx_conc(icplx) = exp(lnQK)
          total = total + reaction%eqsrfcplx_free_site_stoich(icplx)*srfcplx_conc(icplx) 
          
        enddo
        
        if (one_more) exit
        
        if (reaction%eqsrfcplx_rxn_stoich_flag(irxn)) then 
          ! stoichiometry for free sites in one of reactions is not 1, thus must
          ! use nonlinear iteration to solve
          res = site_density(isite)-total
          
          dres_dfree_site = 1.d0

          do j = 1, ncplx
            icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
            dres_dfree_site = dres_dfree_site + &
              reaction%eqsrfcplx_free_site_stoich(icplx)* &
              srfcplx_conc(icplx)/free_site_conc
          enddo

          dfree_site_conc = res / dres_dfree_site
          free_site_conc = free_site_conc + dfree_site_conc
        
          if (dabs(dfree_site_conc/free_site_conc) < tol) then
            one_more = PETSC_TRUE
          endif
        
        else
        
          total = total / free_site_conc
          free_site_conc = site_density(isite) / total  
          
          one_more = PETSC_TRUE 
        
        endif

      enddo ! generic do
      
      auxvar%eqsrfcplx_free_site_conc(ichunk,ithread,irxn) = free_site_conc
   
  !!!!!!!!!!!!
      ! 2.3-46

      ! Sx = free site
      ! mi = molality of component i
      dSx_dmi = 0.d0
      tempreal = 0.d0
      do j = 1, ncplx
        icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
        ncomp = reaction%eqsrfcplxspecid(0,icplx)
        do i = 1, ncomp
          icomp = reaction%eqsrfcplxspecid(i,icplx)
          ! sum of nu_li * nu_i * S_i
          dSx_dmi(icomp) = dSx_dmi(icomp) + reaction%eqsrfcplxstoich(i,icplx)* &
                                            reaction%eqsrfcplx_free_site_stoich(icplx)* &
                                            srfcplx_conc(icplx)
        enddo
        ! sum of nu_i^2 * S_i
        tempreal = tempreal + reaction%eqsrfcplx_free_site_stoich(icplx)* & 
                              reaction%eqsrfcplx_free_site_stoich(icplx)* &
                              srfcplx_conc(icplx)
      enddo 
      ! divide denominator by Sx
      tempreal = tempreal / free_site_conc
      ! add 1.d0 to denominator
      tempreal = tempreal + 1.d0
      ! divide numerator by denominator
      dSx_dmi = -dSx_dmi / tempreal
      ! convert from dlogm to dm
      dSx_dmi = dSx_dmi / auxvar%pri_molal(ichunk,ithread,:)
  !!!!!!!!!!!!
   
      do k = 1, ncplx
        icplx = reaction%eqsrfcplx_rxn_to_complex(k,irxn)

        auxvar%eqsrfcplx_conc(ichunk,ithread,icplx) = srfcplx_conc(icplx)

        ncomp = reaction%eqsrfcplxspecid(0,icplx)
        if (isite == 1) then ! immobile sites  
          do i = 1, ncomp
            icomp = reaction%eqsrfcplxspecid(i,icplx)
            auxvar%total_sorb_eq(ichunk,ithread,icomp) = auxvar%total_sorb_eq(ichunk,ithread,icomp) + &
              reaction%eqsrfcplxstoich(i,icplx)*srfcplx_conc(icplx)
          enddo
        else ! mobile sites
         
        endif
        
        ! for 2.3-47 which feeds into 2.3-50
        nui_Si_over_Sx = reaction%eqsrfcplx_free_site_stoich(icplx)* &
                         srfcplx_conc(icplx)/ &
                         free_site_conc

        do j = 1, ncomp
          jcomp = reaction%eqsrfcplxspecid(j,icplx)
          tempreal = reaction%eqsrfcplxstoich(j,icplx)*srfcplx_conc(icplx) / &
                     auxvar%pri_molal(ichunk,ithread,jcomp)+ &
                     nui_Si_over_Sx*dSx_dmi(jcomp)
          if (isite == 1) then ! immobile sites                  
            do i = 1, ncomp
              icomp = reaction%eqsrfcplxspecid(i,icplx)
              auxvar%dtotal_sorb_eq(ichunk,ithread,icomp,jcomp) = &
                auxvar%dtotal_sorb_eq(ichunk,ithread,icomp,jcomp) + &
                                                   reaction%eqsrfcplxstoich(i,icplx)* &
                                                   tempreal
            enddo ! i
          else ! mobile sites
   
          endif
        enddo ! j
      enddo ! k
    enddo ! isite
  enddo ! irxn
  
  ! units of total_sorb = mol/m^3
  ! units of dtotal_sorb = kg water/m^3 bulk
  
  enddo ! chunk loop
  
end subroutine RTotalSorbEqSurfCplxChunk

! ************************************************************************** !
!
! RKineticMineralChunk: Computes the kinetic mineral precipitation/dissolution
!                        rates
! author: Glenn Hammond
! date: 11/11/10
!
! ************************************************************************** !
subroutine RKineticMineralChunk(Res,Jac,compute_derivative,auxvar,vector_length,ithread,reaction, &
                                option)

  use Option_module
  
  type(option_type) :: option
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Res(vector_length,option%num_threads,reaction%ncomp)
  PetscReal :: Jac(vector_length,option%num_threads,reaction%ncomp,reaction%ncomp)
  
  PetscInt :: i, j, k, imnrl, icomp, jcomp, kcplx, iphase, ncomp, ipref
  PetscReal :: prefactor(10), sum_prefactor_rate
  PetscReal :: dIm_dsum_prefactor_rate, dIm_dprefactor_rate
  PetscReal :: dprefactor_dcomp_numerator, dprefactor_dcomp_denominator
  PetscReal :: tempreal, tempreal2
  PetscReal :: affinity_factor, sign_
  PetscReal :: Im, Im_const, dIm_dQK
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_sec(reaction%neqcplx)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: ln_sec_act(reaction%neqcplx)
  PetscReal :: QK, lnQK, dQK_dCj, dQK_dmj
  PetscBool :: prefactor_exists

  PetscInt, parameter :: needs_to_be_fixed = 1

  PetscInt :: ichunk

  do ichunk = 1, vector_length

  iphase = 1                         

  ln_conc = log(auxvar%pri_molal(ichunk,ithread,:))
  ln_act = ln_conc+log(auxvar%pri_act_coef(ichunk,ithread,:))

  if (reaction%neqcplx > 0) then
    ln_sec = log(auxvar%sec_molal(ichunk,ithread,:))
    ln_sec_act = ln_sec+log(auxvar%sec_act_coef(ichunk,ithread,:))
  endif

#if 0
  if (.not.option%use_isothermal) then
    call ReactionInterpolateLogK(reaction%kinmnrl_logKcoef,reaction%kinmnrl_logK, &
                               auxvar%temp(ichunk,ithread,iphase),reaction%nkinmnrl)
  endif
#endif  

  do imnrl = 1, reaction%nkinmnrl ! for each mineral
    ! compute ion activity product
    lnQK = -reaction%kinmnrl_logK(imnrl)*LOG_TO_LN

    ! activity of water
    if (reaction%kinmnrlh2oid(imnrl) > 0) then
      lnQK = lnQK + reaction%kinmnrlh2ostoich(imnrl)*auxvar%ln_act_h2o(ichunk,ithread)
    endif

    ncomp = reaction%kinmnrlspecid(0,imnrl)
    do i = 1, ncomp
      icomp = reaction%kinmnrlspecid(i,imnrl)
      lnQK = lnQK + reaction%kinmnrlstoich(i,imnrl)*ln_act(icomp)
    enddo
    QK = exp(lnQK)
    
    if (associated(reaction%kinmnrl_Tempkin_const)) then
      affinity_factor = 1.d0-QK**(1.d0/reaction%kinmnrl_Tempkin_const(imnrl))
    else
      affinity_factor = 1.d0-QK
    endif
    
    sign_ = sign(1.d0,affinity_factor)

    if (auxvar%mnrl_volfrac(ichunk,ithread,imnrl) > 0 .or. sign_ < 0.d0) then
    
!     check for supersaturation threshold for precipitation
      if (associated(reaction%kinmnrl_affinity_threshold)) then
        if (sign_ < 0.d0 .and. QK < reaction%kinmnrl_affinity_threshold(imnrl)) cycle
      endif

      ! compute prefactor
      if (reaction%kinmnrl_num_prefactors(imnrl) > 0) then
        print *, 'Kinetic mineral reaction prefactor calculations have not been verified.'
        stop
#if 0
        sum_prefactor_rate = 0
        do ipref = 1, reaction%kinmnrl_num_prefactors(imnrl)
          prefactor(ipref) = 1.d0
          do i = 1, reaction%kinmnrl_prefactor_id(0,ipref,imnrl) ! primary contribution
            icomp = reaction%kinmnrl_prefactor_id(i,ipref,imnrl)
            prefactor(ipref) = prefactor(ipref) * &
              exp(reaction%kinmnrl_pref_alpha(i,ipref,imnrl)* &
              ln_act(icomp))/ &
              (1.d0+reaction%kinmnrl_pref_atten_coef(i,ipref,imnrl)* &
              exp(reaction%kinmnrl_pref_beta(i,ipref,imnrl)* &
              ln_act(icomp)))
          enddo
          sum_prefactor_rate = sum_prefactor_rate + prefactor(ipref)*reaction%kinmnrl_rate(ipref,imnrl)
        enddo
#endif
      else
        sum_prefactor_rate = reaction%kinmnrl_rate(imnrl)
      endif

      ! compute rate
      ! rate = mol/m^2 mnrl/sec
      ! area = m^2 mnrl/m^3 bulk
      ! volume = m^3 bulk
      ! units = m^2 mnrl/m^3 bulk
      Im_const = -auxvar%mnrl_area(ichunk,ithread,imnrl)
      ! units = mol/sec/m^3 bulk
      if (associated(reaction%kinmnrl_affinity_power)) then
        Im = Im_const*sign_*abs(affinity_factor)**reaction%kinmnrl_affinity_power(imnrl)*sum_prefactor_rate
      else
        Im = Im_const*sign_*abs(affinity_factor)*sum_prefactor_rate
      endif
      auxvar%mnrl_rate(ichunk,ithread,imnrl) = Im ! mol/sec/m^3
    else
      auxvar%mnrl_rate(ichunk,ithread,imnrl) = 0.d0
      cycle
    endif
    
    ! units = m^2 mnrl
    Im_const = Im_const*auxvar%vol(ichunk,ithread)
    ! units = mol/sec
    Im = Im*auxvar%vol(ichunk,ithread)

    ncomp = reaction%kinmnrlspecid(0,imnrl)
    do i = 1, ncomp
      icomp = reaction%kinmnrlspecid(i,imnrl)
      Res(ichunk,ithread,icomp) = Res(ichunk,ithread,icomp) + reaction%kinmnrlstoich(i,imnrl)*Im
    enddo 
    
    if (.not. compute_derivative) cycle   

    ! calculate derivatives of rate with respect to free
    ! units = mol/sec
    if (associated(reaction%kinmnrl_affinity_power)) then
      dIm_dQK = -Im*reaction%kinmnrl_affinity_power(imnrl)/abs(affinity_factor)
    else
      dIm_dQK = -Im_const*sum_prefactor_rate
    endif
    if (associated(reaction%kinmnrl_Tempkin_const)) then
      dIm_dQK = dIm_dQK*(1.d0/reaction%kinmnrl_Tempkin_const(imnrl))/QK
    endif
    
    ! derivatives with respect to primary species in reaction quotient
    do j = 1, ncomp
      jcomp = reaction%kinmnrlspecid(j,imnrl)
      ! unit = L water/mol
      dQK_dCj = reaction%kinmnrlstoich(j,imnrl)*exp(lnQK-ln_conc(jcomp))
      ! units = (L water/mol)*(kg water/m^3 water)*(m^3 water/1000 L water) = kg water/mol
      dQK_dmj = dQK_dCj*auxvar%den(ichunk,ithread,iphase)*1.d-3 ! the multiplication by density could be moved
                                   ! outside the loop
      do i = 1, ncomp
        icomp = reaction%kinmnrlspecid(i,imnrl)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(ichunk,ithread,icomp,jcomp) = Jac(ichunk,ithread,icomp,jcomp) + &
                           reaction%kinmnrlstoich(i,imnrl)*dIm_dQK*dQK_dmj
      enddo
    enddo

    if (reaction%kinmnrl_num_prefactors(imnrl) > 0) then ! add contribution of derivative in prefactor - messy
      print *, 'Kinetic mineral reaction prefactor calculations have not been verified.'
      stop
      
    endif
  enddo  ! loop over minerals
  
  enddo ! loop over chunks
    
end subroutine RKineticMineralChunk

! ************************************************************************** !
!
! RAccumulationSorbChunk: Computes non-aqueous portion of the accumulation term in 
!                    residual function
! author: Glenn Hammond
! date: 11/11/10
!
! ************************************************************************** !
subroutine RAccumulationSorbChunk(auxvar,vector_length,ithread,reaction,option,Res)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  type(reaction_type) :: reaction
  PetscReal :: Res(vector_length,option%num_threads,reaction%ncomp)
  
  PetscReal :: v_t
  
  PetscInt :: ichunk
  
  ! units = (mol solute/m^3 bulk)*(m^3 bulk)/(sec) = mol/sec
  ! all residual entries should be in mol/sec
  do ichunk = 1, vector_length
  
  v_t = auxvar%vol(ichunk,ithread)/option%tran_dt
  Res(ichunk,ithread,1:reaction%naqcomp) = Res(ichunk,ithread,1:reaction%naqcomp) + &
    auxvar%total_sorb_eq(ichunk,ithread,:)*v_t

  enddo ! chunk loop

end subroutine RAccumulationSorbChunk

! ************************************************************************** !
!
! RAccumulationSorbDerivChunk: Computes derivative of non-aqueous portion
!                                   of the accumulation term in residual 
!                                   function 
! author: Glenn Hammond
! date: 11/11/10
!
! ************************************************************************** !
subroutine RAccumulationSorbDerivChunk(auxvar,vector_length,ithread,reaction,option,J)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  type(reaction_type) :: reaction
  PetscReal :: J(vector_length,option%num_threads,reaction%ncomp,reaction%ncomp)
  
  PetscInt :: icomp
  PetscReal :: v_t
  
  PetscInt :: ichunk
  
  do ichunk = 1, vector_length
  
  ! units = (kg water/m^3 bulk)*(m^3 bulk)/(sec) = kg water/sec
  ! all Jacobian entries should be in kg water/sec
  v_t = auxvar%vol(ichunk,ithread)/option%tran_dt
  J(ichunk,ithread,1:reaction%naqcomp,1:reaction%naqcomp) = &
    J(ichunk,ithread,1:reaction%naqcomp,1:reaction%naqcomp) + &
    auxvar%dtotal_sorb_eq(ichunk,ithread,:,:)*v_t
    
  enddo

end subroutine RAccumulationSorbDerivChunk

! ************************************************************************** !
!
! RTAuxVarComputeChunk: Computes secondary variables for each grid cell
! author: Glenn Hammond
! date: 11/10/10
!
! ************************************************************************** !
subroutine RTAuxVarComputeChunk(auxvar,vector_length,ithread,reaction,option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  
  ! any changes to the below must also be updated in 
  ! Reaction.F90:RReactionDerivative()
  
  call RTotalChunk(auxvar,vector_length,ithread,reaction,option)
#if 1  
  if (reaction%neqsorb > 0) then
    call RTotalSorbChunk(auxvar,vector_length,ithread,reaction,option)
  endif
#endif
  
end subroutine RTAuxVarComputeChunk

! ************************************************************************** !
!
! RTAccumulationChunk: Computes aqueous portion of the accumulation term in 
!                      residual function
! author: Glenn Hammond
! date: 11/10/10
!
! ************************************************************************** !
subroutine RTAccumulationChunk(auxvar,vector_length,ithread,reaction,option,Res)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  type(reaction_type) :: reaction
  PetscReal :: Res(vector_length,option%num_threads,reaction%ncomp)
  
  PetscInt :: iphase
  PetscInt :: istart, iend
  PetscInt :: idof
  PetscInt :: icoll
  PetscInt :: icollcomp
  PetscInt :: iaqcomp
  PetscReal :: psv_t(vector_length,option%num_threads)
  PetscReal :: v_t(vector_length,option%num_threads)
  
  PetscInt :: ichunk
  
  iphase = 1
  istart = 1
  iend = reaction%naqcomp
  
  ! units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  !         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  ! 1000.d0 converts vol from m^3 -> L
  ! all residual entries should be in mol/sec
  do ichunk = 1, vector_length
    psv_t(ichunk,ithread) = auxvar%por(ichunk,ithread)*auxvar%sat(ichunk,ithread,iphase)*1000.d0* &
                    auxvar%vol(ichunk,ithread)/option%tran_dt  
    Res(ichunk,ithread,istart:iend) = psv_t(ichunk,ithread)*auxvar%total(ichunk,ithread,:,iphase) 
  enddo
  
! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2
  do 
    iphase = iphase + 1
    if (iphase > option%nphase) exit

! super critical CO2 phase
    if (iphase == 2) then
      do ichunk = 1, vector_length
        psv_t(ichunk,ithread) = auxvar%por(ichunk,ithread)*auxvar%sat(ichunk,ithread,iphase)* &
                        1000.d0*auxvar%vol(ichunk,ithread)/option%tran_dt 
        Res(ichunk,ithread,istart:iend) = Res(ichunk,ithread,istart:iend) + &
                                psv_t(ichunk,ithread)*auxvar%total(ichunk,ithread,:,iphase) 
      enddo
      ! should sum over gas component only need more implementations
    endif 
! add code for other phases here
  enddo
#endif
  
end subroutine RTAccumulationChunk

! ************************************************************************** !
!
! RTAccumulationDerivativeChunk: Computes derivative of aqueous portion of the 
!                                accumulation term in residual function 
! author: Glenn Hammond
! date: 11/10/10
!
! ************************************************************************** !
subroutine RTAccumulationDerivativeChunk(auxvar,vector_length,ithread,reaction,option,J)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(react_tran_auxvar_chunk_type) :: auxvar
  PetscInt :: vector_length
  PetscInt :: ithread
  type(reaction_type) :: reaction
  PetscReal :: J(vector_length,option%num_threads,reaction%ncomp,reaction%ncomp)
  
  PetscInt :: icomp, iphase
  PetscInt :: istart, iendaq
  PetscInt :: idof
  PetscInt :: icoll
  PetscReal :: psvd_t(vector_length,option%num_threads)
  PetscReal :: v_t(vector_length,option%num_threads)

  PetscInt :: ichunk

  iphase = 1
  istart = 1
  iendaq = reaction%naqcomp 
  ! units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*(m^3 bulk)/(sec)
  !         *(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  ! all Jacobian entries should be in kg water/sec
  J(:,ithread,:,:) = 0.d0
  do ichunk = 1, vector_length
    ! this result of this conditional will be the same for all grid cells
    if (associated(auxvar%dtotal)) then ! units of dtotal = kg water/L water
      psvd_t(ichunk,ithread) = auxvar%por(ichunk,ithread)*auxvar%sat(ichunk,ithread,iphase)*1000.d0* &
                       auxvar%vol(ichunk,ithread)/option%tran_dt
      J(ichunk,ithread,istart:iendaq,istart:iendaq) = &
        auxvar%dtotal(ichunk,ithread,:,:,iphase)*psvd_t(ichunk,ithread)
    else
      psvd_t(ichunk,ithread) = auxvar%por(ichunk,ithread)*auxvar%sat(ichunk,ithread,iphase)* &
               auxvar%den(ichunk,ithread,iphase)*auxvar%vol(ichunk,ithread)/option%tran_dt ! units of den = kg water/m^3 water
      do icomp=istart,iendaq
        J(ichunk,ithread,icomp,icomp) = psvd_t(ichunk,ithread)
      enddo
    endif
  enddo

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2
  do
    iphase = iphase +1 
    if (iphase > option%nphase) exit
! super critical CO2 phase
    if (iphase == 2) then
      do ichunk = 1, vector_length
        ! this result of this conditional will be the same for all grid cells
        if (associated(auxvar%dtotal)) then
          psvd_t(ichunk,ithread) = auxvar%por(ichunk,ithread)*auxvar%sat(ichunk,ithread,iphase)* &
                           1000.d0*auxvar%vol(ichunk,ithread)/option%tran_dt  
          J(ichunk,ithread,istart:iendaq,istart:iendaq) = &
            J(ichunk,ithread,istart:iendaq,istart:iendaq) + &
            auxvar%dtotal(ichunk,ithread,:,:,iphase)*psvd_t(ichunk,ithread)
        else
          psvd_t(ichunk,ithread) = auxvar%por(ichunk,ithread)*auxvar%sat(ichunk,ithread,iphase)* &
            auxvar%den(ichunk,ithread,iphase)*auxvar%vol(ichunk,ithread)/option%tran_dt ! units of den = kg water/m^3 water
          do icomp=istart,iendaq
            J(ichunk,ithread,icomp,icomp) = J(ichunk,ithread,icomp,icomp) + psvd_t(ichunk,ithread)
          enddo
        endif   
      enddo
    endif
  enddo
#endif

end subroutine RTAccumulationDerivativeChunk

end module Reaction_Chunk_module
