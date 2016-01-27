module Init_Subsurface_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: InitSubsurfAssignMatIDsToRegns, &
            InitSubsurfAssignMatProperties, &
            InitSubsurfSetupRealization, &
            InitSubsurfaceReadRequiredCards, &
            InitSubsurfaceReadInput
  
contains

! ************************************************************************** !

subroutine InitSubsurfSetupRealization(realization)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Realization_class
  use Option_module
  use Logging_module
  use Waypoint_module
  use Init_Common_module
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF
  use Reaction_Database_module
  use EOS_Water_module
  use Dataset_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  PetscReal :: dum1
  PetscErrorCode :: ierr
  
  option => realization%option
  
  call PetscLogEventBegin(logging%event_setup,ierr);CHKERRQ(ierr)
  
  ! initialize reference density
  if (option%reference_water_density < 1.d-40) then
    call EOSWaterDensity(option%reference_temperature, &
                         option%reference_pressure, &
                         option%reference_water_density, &
                         dum1,ierr)    
  endif

  ! read reaction database
  if (associated(realization%reaction)) then
    if (realization%reaction%use_full_geochemistry) then
        call DatabaseRead(realization%reaction,option)
        call BasisInit(realization%reaction,option)
    else
      ! turn off activity coefficients since the database has not been read
      realization%reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
      allocate(realization%reaction%primary_species_print(option%ntrandof))
      realization%reaction%primary_species_print = PETSC_TRUE
    endif
  endif

  ! SK 09/30/13, Added to check if Mphase is called with OS
  if (option%transport%reactive_transport_coupling == OPERATOR_SPLIT .and. &
      option%iflowmode == MPH_MODE) then
    option%io_buffer = 'Operator split not implemented with MPHASE. ' // &
                       'Switching to Global Implicit.'
    !geh: We should force the user to switch without automatically switching
!    call printWrnMsg(option)
    call printErrMsg(option)
    option%transport%reactive_transport_coupling = GLOBAL_IMPLICIT
  endif
  
  ! create grid and allocate vectors
  call RealizationCreateDiscretization(realization)
  
  ! read any regions provided in external files
  call InitCommonReadRegionFiles(realization)
  ! clip regions and set up boundary connectivity, distance  
  call RealizationLocalizeRegions(realization)
  call RealizationPassPtrsToPatches(realization)
  call RealizationProcessDatasets(realization)
  ! link conditions with regions through couplers and generate connectivity
  call RealProcessMatPropAndSatFunc(realization)
  ! must process conditions before couplers in order to determine dataset types
  call RealizationProcessConditions(realization)
  call RealizationProcessCouplers(realization)
  call SubsurfSandboxesSetup(realization)
  call RealProcessFluidProperties(realization)
  call SubsurfInitMaterialProperties(realization)
  ! assignVolumesToMaterialAuxVars() must be called after 
  ! RealizInitMaterialProperties() where the Material object is created 
  call SubsurfAssignVolsToMatAuxVars(realization)
  call RealizationInitAllCouplerAuxVars(realization)
  if (option%ntrandof > 0) then
    call printMsg(option,"  Setting up TRAN Realization ")
    call RealizationInitConstraints(realization)
    call printMsg(option,"  Finished setting up TRAN Realization ")  
  endif
  call RealizationPrintCouplers(realization)
  if (.not.option%steady_state) then
    ! add waypoints associated with boundary conditions, source/sinks etc. to list
    call RealizationAddWaypointsToList(realization)
    ! fill in holes in waypoint data
!    call WaypointListFillIn(option,realization%waypoint_list)
!    call WaypointListRemoveExtraWaypnts(option,realization%waypoint_list)
  endif
  call PetscLogEventEnd(logging%event_setup,ierr);CHKERRQ(ierr)
  
#ifdef OS_STATISTICS
  call RealizationPrintGridStatistics(realization)
#endif

#if defined(PETSC_HAVE_HDF5)
#if !defined(HDF5_BROADCAST)
  call printMsg(option,"Default HDF5 method is used in Initialization")
#else
  call printMsg(option,"Glenn's HDF5 broadcast method is used in Initialization")
#endif
#endif

end subroutine InitSubsurfSetupRealization

! ************************************************************************** !

subroutine SubsurfInitMaterialProperties(realization)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/14
  ! 
  use Realization_class
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  call SubsurfAllocMatPropDataStructs(realization)
  call InitSubsurfAssignMatIDsToRegns(realization)
  call InitSubsurfAssignMatProperties(realization)
  
end subroutine SubsurfInitMaterialProperties

! ************************************************************************** !

subroutine SubsurfAllocMatPropDataStructs(realization)
  ! 
  ! Allocates data structures associated with storage of material properties
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/14
  ! 
  use Realization_class
  use Material_module
  use Option_module
  use Discretization_module
  use Grid_module
  use Patch_module
  use Material_Aux_class
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  PetscInt :: ghosted_id
  PetscInt :: istart, iend
  PetscInt :: i
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: cur_patch
  class(material_auxvar_type), pointer :: material_auxvars(:)

  option => realization%option

  ! initialize material auxiliary indices.  this does not have to be done
  ! for each patch, just once.
  call MaterialInitAuxIndices(realization%patch%material_property_array,option)

  ! create mappinging
  ! loop over all patches and allocation material id arrays
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    grid => cur_patch%grid
    if (.not.associated(cur_patch%imat)) then
      allocate(cur_patch%imat(grid%ngmax))
      ! initialize to "unset"
      cur_patch%imat = UNINITIALIZED_INTEGER
      ! also allocate saturation function id
      allocate(cur_patch%sat_func_id(grid%ngmax))
      cur_patch%sat_func_id = UNINITIALIZED_INTEGER
    endif
    
    cur_patch%aux%Material => MaterialAuxCreate()
    allocate(material_auxvars(grid%ngmax))
    do ghosted_id = 1, grid%ngmax
      call MaterialAuxVarInit(material_auxvars(ghosted_id),option)
    enddo
    cur_patch%aux%Material%num_aux = grid%ngmax
    cur_patch%aux%Material%auxvars => material_auxvars
    nullify(material_auxvars)
    
    cur_patch => cur_patch%next
  enddo

  ! Create Vec that holds compressibility
  if (soil_compressibility_index > 0) then
    call DiscretizationDuplicateVector(realization%discretization, &
                                       realization%field%work, &
                                       realization%field%compressibility0)
  endif

end subroutine SubsurfAllocMatPropDataStructs

! ************************************************************************** !

subroutine InitSubsurfAssignMatIDsToRegns(realization)
  ! 
  ! Assigns material properties to associated regions in the model
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 
  use Realization_class
  use Strata_module
  use Region_module
  use Option_module
  use Grid_module
  use Patch_module
  use Field_module
  use Material_module
  use Material_Aux_class
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  PetscInt :: icell, local_id, ghosted_id
  PetscInt :: istart, iend
  PetscInt :: local_min, global_min
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(strata_type), pointer :: strata
  type(patch_type), pointer :: cur_patch

  type(material_property_type), pointer :: material_property
  type(region_type), pointer :: region
  class(material_auxvar_type), pointer :: material_auxvars(:)
  
  option => realization%option

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    ! set material ids to uninitialized
    material_auxvars => cur_patch%aux%Material%auxvars
    cur_patch%imat = UNINITIALIZED_INTEGER
    grid => cur_patch%grid
    strata => cur_patch%strata_list%first
    do
      if (.not.associated(strata)) exit
      ! if not within time period specified, skip the strata.
      ! use a one second tolerance on the start time and end time
      if (StrataWithinTimePeriod(strata,option%time)) then
        ! Read in cell by cell material ids if they exist
        if (.not.associated(strata%region) .and. strata%active) then
          call SubsurfReadMaterialIDsFromFile(realization, &
                                              strata%realization_dependent, &
                                              strata%material_property_filename)
        ! Otherwise, set based on region
        else
          region => strata%region
          material_property => strata%material_property
          if (associated(region)) then
            istart = 1
            iend = region%num_cells
          else
            istart = 1
            iend = grid%nlmax
          endif
          do icell=istart, iend
            if (associated(region)) then
              local_id = region%cell_ids(icell)
            else
              local_id = icell
            endif
            ghosted_id = grid%nL2G(local_id)
            if (strata%active) then
              cur_patch%imat(ghosted_id) = material_property%internal_id
            else
              ! if not active, set material id to zero
              cur_patch%imat(ghosted_id) = 0
            endif
          enddo
        endif
      endif
      strata => strata%next
    enddo
    cur_patch => cur_patch%next
  enddo
  
  ! ensure that ghosted values for material ids are up to date
  call RealLocalToLocalWithArray(realization,MATERIAL_ID_ARRAY)
  
  ! set material ids in material auxvar.  this must come after the update of 
  ! ghost values.
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    ! set material ids to uninitialized
    material_auxvars => cur_patch%aux%Material%auxvars
    grid => cur_patch%grid
    do ghosted_id = 1, grid%ngmax
      material_auxvars(ghosted_id)%id = cur_patch%imat(ghosted_id)
    enddo
    cur_patch => cur_patch%next
  enddo

end subroutine InitSubsurfAssignMatIDsToRegns

 ! ************************************************************************** !

subroutine InitSubsurfAssignMatProperties(realization)
  ! 
  ! Assigns material properties based on material ids
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/14
  ! 
  use Realization_class
  use Grid_module
  use Discretization_module
  use Field_module
  use Patch_module
  use Material_Aux_class
  use Material_module
  use Option_module
  use Creep_Closure_module
  use Fracture_module
  use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, PERMEABILITY_XY, &
                               PERMEABILITY_YZ, PERMEABILITY_XZ, &
                               TORTUOSITY, POROSITY, SOIL_COMPRESSIBILITY
  use HDF5_module
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(realization_subsurface_type) :: realization
  
  PetscReal, pointer :: icap_loc_p(:)
  PetscReal, pointer :: ithrm_loc_p(:)
  PetscReal, pointer :: por0_p(:)
  PetscReal, pointer :: tor0_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: perm_xz_p(:)
  PetscReal, pointer :: perm_xy_p(:)
  PetscReal, pointer :: perm_yz_p(:)
  PetscReal, pointer :: perm_pow_p(:)
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: compress_p(:)
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(material_property_type), pointer :: material_property, null_material_property
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(material_type), pointer :: Material
  PetscInt :: local_id, ghosted_id, material_id
  PetscInt :: i
  PetscInt :: tempint
  PetscReal :: tempreal
  PetscErrorCode :: ierr
  
  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  
  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_material_property => MaterialPropertyCreate()
  if (option%nflowdof > 0) then
    call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%perm0_xx,perm_xx_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%perm0_yy,perm_yy_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%perm0_zz,perm_zz_p,ierr);CHKERRQ(ierr)
    if (soil_compressibility_index > 0) then
      call VecGetArrayF90(field%compressibility0,compress_p,ierr);CHKERRQ(ierr)
    endif
  endif
  call VecGetArrayF90(field%porosity0,por0_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%tortuosity0,tor0_p,ierr);CHKERRQ(ierr)
        
  ! have to use Material%auxvars() and not material_auxvars() due to memory
  ! errors in gfortran
  Material => patch%aux%Material

  !if material is associated with fracture, then allocate memory.  
  do ghosted_id = 1, grid%ngmax
    material_id = patch%imat(ghosted_id)
    if (material_id > 0) then
      material_property => &
        patch%material_property_array(material_id)%ptr
      call FractureAuxvarInit(material_property%fracture, &
        patch%aux%Material%auxvars(ghosted_id))
    endif
  enddo
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    material_id = patch%imat(ghosted_id)
    if (material_id == 0) then ! accommodate inactive cells
      material_property => null_material_property
    else if (material_id > 0 .and. &
              material_id <= &
              size(patch%material_property_array)) then
      material_property => &
        patch%material_property_array(material_id)%ptr
      if (.not.associated(material_property)) then
        write(string,*) &
          patch%imat_internal_to_external(material_id)
        option%io_buffer = 'No material property for material id ' // &
                            trim(adjustl(string)) &
                            //  ' defined in input file.'
        call printErrMsgByRank(option)
      endif
    else if (material_id < 0 .and. Initialized(material_id)) then 
      ! highjacking dataset_name and group_name for error processing
      write(string,*) grid%nG2A(ghosted_id)
      write(string2,*) -1*material_id
      option%io_buffer = 'Undefined material id ' // &
                          trim(adjustl(string2)) // &
                          ' at cell ' // &
                          trim(adjustl(string)) // '.'
      call printErrMsgByRank(option)
    else if (Uninitialized(material_id)) then 
      write(string,*) grid%nG2A(ghosted_id)
      option%io_buffer = 'Uninitialized material id in patch at cell ' // &
                          trim(adjustl(string))
      call printErrMsgByRank(option)
    else if (material_id > size(patch%material_property_array)) then
      write(option%io_buffer,*) patch%imat_internal_to_external(material_id)
      option%io_buffer = 'Unmatched material id in patch: ' // &
        adjustl(trim(option%io_buffer))
      call printErrMsgByRank(option)
    else
      option%io_buffer = 'Something messed up with material ids. ' // &
        ' Possibly material ids not assigned to all grid cells. ' // &
        ' Contact Glenn!'
      call printErrMsgByRank(option)
    endif
    if (option%nflowdof > 0) then
      patch%sat_func_id(ghosted_id) = &
        material_property%saturation_function_id
      icap_loc_p(ghosted_id) = material_property%saturation_function_id
      ithrm_loc_p(ghosted_id) = material_property%internal_id
      perm_xx_p(local_id) = material_property%permeability(1,1)
      perm_yy_p(local_id) = material_property%permeability(2,2)
      perm_zz_p(local_id) = material_property%permeability(3,3)
      if (soil_compressibility_index > 0) then
        compress_p(local_id) = material_property%soil_compressibility
      endif
    endif
    if (associated(Material%auxvars)) then
      call MaterialAssignPropertyToAux(Material%auxvars(ghosted_id), &
                                        material_property,option)
    endif
    por0_p(local_id) = material_property%porosity
    tor0_p(local_id) = material_property%tortuosity
  enddo
  call MaterialPropertyDestroy(null_material_property)

  if (option%nflowdof > 0) then
    call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_xx,perm_xx_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_yy,perm_yy_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_zz,perm_zz_p,ierr);CHKERRQ(ierr)
    if (soil_compressibility_index > 0) then
      call VecRestoreArrayF90(field%compressibility0,compress_p,ierr); &
                              CHKERRQ(ierr)
    endif
  endif
  call VecRestoreArrayF90(field%porosity0,por0_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%tortuosity0,tor0_p,ierr);CHKERRQ(ierr)
        
  ! read in any user-defined property fields
  do material_id = 1, size(patch%material_property_array)
    material_property => &
            patch%material_property_array(material_id)%ptr
    if (associated(material_property)) then
      if (associated(material_property%permeability_dataset)) then
        call SubsurfReadPermsFromFile(realization,material_property)
      endif
      if (associated(material_property%compressibility_dataset)) then
!        call SubsurfReadCompressFromFile(realization,material_property)
        call SubsurfReadDatasetToVecWithMask(realization, &
               material_property%compressibility_dataset, &
               material_property%internal_id,field%compressibility0)
      endif
      if (associated(material_property%porosity_dataset)) then
        call SubsurfReadDatasetToVecWithMask(realization, &
               material_property%porosity_dataset, &
               material_property%internal_id,field%porosity0)
      endif
    endif
  enddo
      
  ! update ghosted values
  if (option%nflowdof > 0) then
    call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                                     field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_X,ZERO_INTEGER)
    call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                                     field%work_loc,ONEDOF)  
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Y,ZERO_INTEGER)
    call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                                     field%work_loc,ONEDOF)   
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Z,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%icap_loc, &
                                    field%icap_loc,ONEDOF)   
    call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
                                    field%ithrm_loc,ONEDOF)
    call RealLocalToLocalWithArray(realization,SATURATION_FUNCTION_ID_ARRAY)
    
    if (soil_compressibility_index > 0) then
      call DiscretizationGlobalToLocal(discretization,field%compressibility0, &
                                       field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   SOIL_COMPRESSIBILITY,ZERO_INTEGER)
    endif
  endif
  
  call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                                   field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,POROSITY_MINERAL)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,POROSITY_CURRENT)
  call DiscretizationGlobalToLocal(discretization,field%tortuosity0, &
                                    field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               TORTUOSITY,ZERO_INTEGER)
  ! rock properties
  do i = 1, max_material_index
    call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    do local_id = 1, patch%grid%nlmax
      ghosted_id = patch%grid%nL2G(local_id)
      if (patch%imat(ghosted_id) > 0) then
        vec_p(local_id) = &
          Material%auxvars(patch%grid%nL2G(local_id))%soil_properties(i)
      else
        vec_p(local_id) = 1.d-40 ! some initialized value for inactive cells.
      endif
    enddo
    call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    call VecMin(field%work,tempint,tempreal,ierr)
    if (Uninitialized(tempreal)) then
      option%io_buffer = 'Incorrect assignment of soil properties. ' // &
        'Please send this error message and your input file to ' // &
        'pflotran-dev@googlegroups.com.'
        call printErrMsg(option)
    endif
    call DiscretizationGlobalToLocal(discretization,field%work, &
                                     field%work_loc,ONEDOF)
    call VecGetArrayF90(field%work_loc,vec_p,ierr);CHKERRQ(ierr)
    do ghosted_id = 1, patch%grid%ngmax
      Material%auxvars(ghosted_id)%soil_properties(i) = &
         vec_p(ghosted_id)
    enddo
    call VecRestoreArrayF90(field%work_loc,vec_p,ierr);CHKERRQ(ierr)
  enddo
  
  if (associated(creep_closure)) then
    material_property => &
      MaterialPropGetPtrFromArray(creep_closure%material_name, &
                                  patch%material_property_array)
    creep_closure%imat = material_property%internal_id
  endif
  
end subroutine InitSubsurfAssignMatProperties

! ************************************************************************** !

subroutine SubsurfReadMaterialIDsFromFile(realization,realization_dependent, &
                                          filename)
  ! 
  ! Reads in grid cell materials
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/03/08
  ! 

  use Realization_class
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
  use Discretization_module
  use Logging_module
  use Input_Aux_module
  use Material_module

  use HDF5_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
  
  class(realization_subsurface_type) :: realization
  PetscBool :: realization_dependent
  character(len=MAXSTRINGLENGTH) :: filename
  
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch   
  type(input_type), pointer :: input
  type(discretization_type), pointer :: discretization
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscBool :: append_realization_id
  PetscInt :: ghosted_id, natural_id, material_id
  PetscInt :: fid = 86
  PetscInt :: status
  PetscInt, pointer :: external_to_internal_mapping(:)
  Vec :: global_vec
  Vec :: local_vec
  PetscErrorCode :: ierr

  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  discretization => realization%discretization
  
  if (index(filename,'.h5') > 0) then
    group_name = 'Materials'
    dataset_name = 'Material Ids'
    call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                    option)
    call DiscretizationCreateVector(discretization,ONEDOF,local_vec,LOCAL, &
                                    option)
    call HDF5ReadCellIndexedIntegerArray(realization,global_vec, &
                                         filename,group_name, &
                                         dataset_name,realization_dependent)
    call DiscretizationGlobalToLocal(discretization,global_vec,local_vec,ONEDOF)
    call GridCopyVecToIntegerArray(grid,patch%imat,local_vec,grid%ngmax)
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
    call VecDestroy(local_vec,ierr);CHKERRQ(ierr)
  else
    call PetscLogEventBegin(logging%event_hash_map,ierr);CHKERRQ(ierr)
    call GridCreateNaturalToGhostedHash(grid,option)
    input => InputCreate(IUNIT_TEMP,filename,option)
    do
      call InputReadPflotranString(input,option)
      if (InputError(input)) exit
      call InputReadInt(input,option,natural_id)
      call InputErrorMsg(input,option,'natural id','STRATA')
      ! natural ids in hash are zero-based
      ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (ghosted_id > 0) then
        call InputReadInt(input,option,material_id)
        call InputErrorMsg(input,option,'material id','STRATA')
        patch%imat(ghosted_id) = material_id
      endif
    enddo
    call InputDestroy(input)
    call GridDestroyHashTable(grid)
    call PetscLogEventEnd(logging%event_hash_map,ierr);CHKERRQ(ierr)
  endif
  
  call MaterialCreateExtToIntMapping(patch%material_property_array, &
                                     external_to_internal_mapping)
  call MaterialApplyMapping(external_to_internal_mapping,patch%imat)
  deallocate(external_to_internal_mapping)
  nullify(external_to_internal_mapping)
  
end subroutine SubsurfReadMaterialIDsFromFile

! ************************************************************************** !

subroutine SubsurfReadPermsFromFile(realization,material_property)
  ! 
  ! Reads in grid cell permeabilities
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/19/09
  ! 

  use Realization_class
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
  use Discretization_module
  use Logging_module
  use Input_Aux_module
  use Material_module
  use HDF5_module
  use Dataset_Common_HDF5_class
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(realization_subsurface_type) :: realization
  type(material_property_type) :: material_property

  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(discretization_type), pointer :: discretization
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: local_id
  PetscInt :: idirection, temp_int
  PetscReal :: ratio, scale
  Vec :: global_vec
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)

  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  discretization => realization%discretization

  call VecGetArrayF90(field%perm0_xx,perm_xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm0_yy,perm_yy_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm0_zz,perm_zz_p,ierr);CHKERRQ(ierr)
  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  if (material_property%isotropic_permeability .or. &
      (.not.material_property%isotropic_permeability .and. &
        material_property%vertical_anisotropy_ratio > 0.d0)) then
!    material_property%permeability_dataset%name = 'Permeability'
    !geh: Pass in -1 so that entire dataset is read. The mask is applied below.
    call SubsurfReadDatasetToVecWithMask(realization, &
                                    material_property%permeability_dataset, &
                                    -1,global_vec)
    call VecGetArrayF90(global_vec,vec_p,ierr);CHKERRQ(ierr)
    ratio = 1.d0
    scale = 1.d0
    !TODO(geh): fix so that ratio and scale work for perms outside
    ! of dataset
    if (material_property%vertical_anisotropy_ratio > 0.d0) then
      ratio = material_property%vertical_anisotropy_ratio
    endif
    if (material_property%permeability_scaling_factor > 0.d0) then
      scale = material_property%permeability_scaling_factor
    endif
    do local_id = 1, grid%nlmax
      if (patch%imat(grid%nL2G(local_id)) == &
          material_property%internal_id) then
        perm_xx_p(local_id) = vec_p(local_id)*scale
        perm_yy_p(local_id) = vec_p(local_id)*scale
        perm_zz_p(local_id) = vec_p(local_id)*ratio*scale
      endif
    enddo
    call VecRestoreArrayF90(global_vec,vec_p,ierr);CHKERRQ(ierr)
  else
    temp_int = Z_DIRECTION
    do idirection = X_DIRECTION,temp_int
      select case(idirection)
        case(X_DIRECTION)
          word = 'X'
        case(Y_DIRECTION)
          word = 'Y'
        case(Z_DIRECTION)
          word = 'Z'
      end select
      select type(dataset => material_property%permeability_dataset)
        class is(dataset_common_hdf5_type)
          dataset%hdf5_dataset_name = trim(dataset%hdf5_dataset_name) // &
                                      trim(word)
      end select
      !geh: Pass in -1 so that entire dataset is read. The mask is applied 
      !     below.
      call SubsurfReadDatasetToVecWithMask(realization, &
                                      material_property%permeability_dataset, &
                                      -1,global_vec)
      call VecGetArrayF90(global_vec,vec_p,ierr);CHKERRQ(ierr)
      select case(idirection)
        case(X_DIRECTION)
          do local_id = 1, grid%nlmax
            if (patch%imat(grid%nL2G(local_id)) == &
                material_property%internal_id) then
              perm_xx_p(local_id) = vec_p(local_id)
            endif
          enddo
        case(Y_DIRECTION)
          do local_id = 1, grid%nlmax
            if (patch%imat(grid%nL2G(local_id)) == &
                material_property%internal_id) then
              perm_yy_p(local_id) = vec_p(local_id)
            endif
          enddo
        case(Z_DIRECTION)
          do local_id = 1, grid%nlmax
            if (patch%imat(grid%nL2G(local_id)) == &
                material_property%internal_id) then
              perm_zz_p(local_id) = vec_p(local_id)
            endif
          enddo
      end select
      call VecRestoreArrayF90(global_vec,vec_p,ierr);CHKERRQ(ierr)
    enddo
  endif
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  
  call VecRestoreArrayF90(field%perm0_xx,perm_xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm0_yy,perm_yy_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm0_zz,perm_zz_p,ierr);CHKERRQ(ierr)
  
end subroutine SubsurfReadPermsFromFile

! ************************************************************************** !

subroutine SubsurfReadDatasetToVecWithMask(realization,dataset,material_id, &
                                           vec)
  ! 
  ! Reads a dataset into a PETSc Vec using the material id as a mask
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/19/2016
  ! 

  use Realization_class
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
  use Logging_module
  use Input_Aux_module
  use Material_module
  use HDF5_module
  use Dataset_Base_class
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(realization_subsurface_type) :: realization
  class(dataset_base_type) :: dataset
  PetscInt :: material_id
  Vec :: vec

  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: local_id, ghosted_id, natural_id
  PetscReal :: tempreal
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: work_p(:)
  
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option

  call VecGetArrayF90(vec,vec_p,ierr);CHKERRQ(ierr)
  
  if (index(dataset%filename,'.h5') > 0) then
    group_name = ''
    dataset_name = dataset%name
    select type(dataset)
      class is(dataset_gridded_hdf5_type)
        call DatasetGriddedHDF5Load(dataset,option)
        do local_id = 1, grid%nlmax
          ghosted_id = grid%nL2G(local_id)
          if (material_id < 0 .or. &
              patch%imat(ghosted_id) == material_id) then
            call DatasetGriddedHDF5InterpolateReal(dataset, &
                   grid%x(ghosted_id),grid%y(ghosted_id),grid%z(ghosted_id), &
                   vec_p(local_id),option)
          endif
        enddo
        ! now we strip the dataset to save storage, saving only the name
        ! and filename
        filename = dataset%filename
        dataset_name = dataset%name
        call DatasetGriddedHDF5Strip(dataset)
        call DatasetGriddedHDF5Init(dataset)
      class is(dataset_common_hdf5_type)
        dataset_name = dataset%hdf5_dataset_name
        call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                          dataset%filename, &
                                          group_name,dataset_name, &
                                          dataset%realization_dependent)
        call VecGetArrayF90(field%work,work_p,ierr);CHKERRQ(ierr)
        if (material_id < 0) then
          vec_p(:) = work_p(:)
        else
          do local_id = 1, grid%nlmax
            if (patch%imat(grid%nL2G(local_id)) == material_id) then
              vec_p(local_id) = work_p(local_id)
            endif
          enddo
        endif
        call VecRestoreArrayF90(field%work,work_p,ierr);CHKERRQ(ierr)
      class default
        option%io_buffer = 'Dataset "' // trim(dataset%name) // '" is of ' // &
          'the wrong type for SubsurfReadDatasetToVecWithMask()'
        call printErrMsg(option)
    end select
  else
    call PetscLogEventBegin(logging%event_hash_map,ierr);CHKERRQ(ierr)
    call GridCreateNaturalToGhostedHash(grid,option)
    input => InputCreate(IUNIT_TEMP,dataset%filename,option)
    do
      call InputReadPflotranString(input,option)
      if (InputError(input)) exit
      call InputReadInt(input,option,natural_id)
      call InputErrorMsg(input,option,'ASCII natural id', &
                         'SubsurfReadDatasetToVecWithMask')
      ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (ghosted_id > 0) then
        if (material_id < 0 .or. &
            patch%imat(ghosted_id) == material_id) then
          local_id = grid%nG2L(ghosted_id)
          if (local_id > 0) then
            call InputReadDouble(input,option,tempreal)
            call InputErrorMsg(input,option,'dataset value', &
                               'SubsurfReadDatasetToVecWithMask')
            vec_p(local_id) = tempreal
          endif
        endif
      endif
    enddo
    call InputDestroy(input)
    call GridDestroyHashTable(grid)
    call PetscLogEventEnd(logging%event_hash_map,ierr);CHKERRQ(ierr)
  endif
  
  call VecRestoreArrayF90(vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine SubsurfReadDatasetToVecWithMask

! ************************************************************************** !

subroutine SubsurfAssignVolsToMatAuxVars(realization)
  ! 
  ! Assigns the cell volumes currently stored in field%volume0 to the 
  ! material auxiliary variable object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/14, 12/04/14
  ! 

  use Realization_class
  use Option_module
  use Material_module
  use Discretization_module
  use Field_module
  use Variables_module, only : VOLUME
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field

  option => realization%option
  field => realization%field

  call DiscretizationGlobalToLocal(realization%discretization,field%volume0, &
                                   field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                               field%work_loc,VOLUME,ZERO_INTEGER)

end subroutine SubsurfAssignVolsToMatAuxVars

! ************************************************************************** !

subroutine SubsurfSandboxesSetup(realization)
  ! 
  ! Initializes sandbox objects.
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/06/14, 12/04/14

  use Realization_class
  use SrcSink_Sandbox_module
  
  class(realization_subsurface_type) :: realization
  
   call SSSandboxSetup(realization%patch%region_list,realization%option)
  
end subroutine SubsurfSandboxesSetup

! ************************************************************************** !

subroutine InitSubsurfaceReadRequiredCards(simulation)
  ! 
  ! Reads required cards from input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07, refactored 08/20/14, refactored 12/10/14
  ! 

  use Option_module
  use Discretization_module
  use Grid_module
  use Input_Aux_module
  use String_module
  use Patch_module
  use Realization_class
  use HDF5_Aux_module

  use Simulation_Subsurface_class
  use General_module
  use Reaction_module  
  use Reaction_Aux_module  
  use Init_Common_module

  implicit none

  class(subsurface_simulation_type) :: simulation

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
  type(patch_type), pointer :: patch, patch2 
  type(grid_type), pointer :: grid
  class(realization_subsurface_type), pointer :: realization
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  
  realization => simulation%realization
  patch => realization%patch
  option => realization%option
  discretization => realization%discretization
  
  input => realization%input
  
! Read in select required cards
!.........................................................................
  
  ! GRID information - GRID is a required card for every simulation
  string = "GRID"
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)

  call DiscretizationReadRequiredCards(discretization,input,option)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      patch => PatchCreate()
      patch%grid => discretization%grid
      if (.not.associated(realization%patch_list)) then
        realization%patch_list => PatchCreateList()
      endif
      call PatchAddToList(patch,realization%patch_list)
      realization%patch => patch
  end select

  ! optional required cards - yes, an oxymoron, but we need to know if
  ! these exist before we can go any further.
  rewind(input%fid)  
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    select case(trim(card))

!....................
      case('DBASE_FILENAME')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'filename','DBASE_FILENAME')
        if (index(word,'.h5') > 0) then
#if defined(PETSC_HAVE_HDF5)
          call HDF5ReadDbase(word,option)
#endif
        else
          call InputReadASCIIDbase(word,option)
        endif

!....................
#if defined(SCORPIO)
      case('HDF5_WRITE_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_write_group_size)
        call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')
        call InputSkipToEnd(input,option,'HDF5_WRITE_GROUP_SIZE')

      case('HDF5_READ_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_read_group_size)
        call InputErrorMsg(input,option,'HDF5_READ_GROUP_SIZE','Group size')
#endif

!....................
      case('PROC')
        ! processor decomposition
        if (realization%discretization%itype == STRUCTURED_GRID) then
          grid => realization%patch%grid
          ! strip card from front of string
          call InputReadInt(input,option,grid%structured_grid%npx)
          call InputDefaultMsg(input,option,'npx')
          call InputReadInt(input,option,grid%structured_grid%npy)
          call InputDefaultMsg(input,option,'npy')
          call InputReadInt(input,option,grid%structured_grid%npz)
          call InputDefaultMsg(input,option,'npz')
 
          if (option%myrank == option%io_rank .and. &
              option%print_to_screen) then
            option%io_buffer = ' Processor Decomposition:'
            call printMsg(option)
            write(option%io_buffer,'("  npx   = ",3x,i4)') &
              grid%structured_grid%npx
            call printMsg(option)
            write(option%io_buffer,'("  npy   = ",3x,i4)') &
              grid%structured_grid%npy
            call printMsg(option)
            write(option%io_buffer,'("  npz   = ",3x,i4)') &
              grid%structured_grid%npz
            call printMsg(option)
          endif
  
          if (option%mycommsize /= grid%structured_grid%npx * &
                                 grid%structured_grid%npy * &
                                 grid%structured_grid%npz) then
            write(option%io_buffer,*) 'Incorrect number of processors specified: ', &
                           grid%structured_grid%npx*grid%structured_grid%npy* &
                           grid%structured_grid%npz,' commsize = ',option%mycommsize
            call printErrMsg(option)
          endif
        endif
  
!....................
      case('CHEMISTRY')
        if (.not.associated(simulation%rt_process_model_coupler)) then
          option%io_buffer = 'CHEMISTRY card included when no ' // &
            'SUBSURFACE_TRANSPORT process model included in SIMULATION block.'
          call printErrMsg(option)
        endif
        !geh: for some reason, we need this with CHEMISTRY read for 
        !     multicontinuum
 !       option%use_mc = PETSC_TRUE
        call ReactionInit(realization%reaction,input,option)
    end select
  enddo
  
#if defined(SCORPIO)
  call InitCommonCreateIOGroups(option)
#endif  

end subroutine InitSubsurfaceReadRequiredCards

! ************************************************************************** !

subroutine InitSubsurfaceReadInput(simulation)
  ! 
  ! Reads pflow input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  use Option_module
  use Field_module
  use Grid_module
  use Grid_Unstructured_Aux_module
  use Grid_Structured_module
  use Solver_module
  use Material_module
  use Saturation_Function_module  
  use Characteristic_Curves_module
  use Dataset_Base_class
  use Dataset_module
  use Dataset_Common_HDF5_class
  use Fluid_module
  use Realization_class
  use Realization_Base_class
  use Region_module
  use Condition_module
  use Transport_Constraint_module
  use Coupler_module
  use Strata_module
  use Observation_module
  use Integral_Flux_module
  use Waypoint_module
  use Debug_module
  use Patch_module
  use Reaction_module
  use Reaction_Aux_module
  use Discretization_module
  use Input_Aux_module
  use String_module
  use Units_module
  use Uniform_Velocity_module
  use Reaction_Mineral_module
  use Regression_module
  use Output_Aux_module
  use Output_Tecplot_module
  use Data_Mediator_Dataset_class
  use EOS_module
  use EOS_Water_module
  use SrcSink_Sandbox_module
  use Klinkenberg_module
  use WIPP_module
  
  use Simulation_Subsurface_class
  use PMC_Subsurface_class
  use Timestepper_BE_class
  use Timestepper_Steady_class
  
#ifdef SOLID_SOLUTION
  use Reaction_Solid_Solution_module, only : SolidSolutionReadFromInputFile
#endif
 
  implicit none
  
  class(subsurface_simulation_type) :: simulation

  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: units_category
    
  PetscBool :: continuation_flag
  
  character(len=1) :: backslash
  PetscReal :: temp_real, temp_real2
  PetscReal :: units_conversion
  PetscInt :: temp_int
  PetscInt :: count, id
  
  PetscBool :: vel_cent
  PetscBool :: vel_face
  PetscBool :: fluxes
  PetscBool :: mass_flowrate
  PetscBool :: energy_flowrate
  PetscBool :: aveg_mass_flowrate
  PetscBool :: aveg_energy_flowrate
  
  type(region_type), pointer :: region
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(tran_constraint_type), pointer :: tran_constraint
  type(tran_constraint_type), pointer :: sec_tran_constraint
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata
  type(observation_type), pointer :: observation
  type(integral_flux_type), pointer :: integral_flux
  
  type(waypoint_type), pointer :: waypoint
  
  type(material_property_type), pointer :: material_property
  type(fluid_property_type), pointer :: fluid_property
  type(saturation_function_type), pointer :: saturation_function
  class(characteristic_curves_type), pointer :: characteristic_curves

  class(realization_subsurface_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch   
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(uniform_velocity_dataset_type), pointer :: uniform_velocity_dataset
  class(dataset_base_type), pointer :: dataset
  class(data_mediator_dataset_type), pointer :: flow_data_mediator
  class(data_mediator_dataset_type), pointer :: rt_data_mediator
  type(input_type), pointer :: input, input_parent
  
  PetscReal :: dt_init
  PetscReal :: dt_min
  
  class(timestepper_BE_type), pointer :: flow_timestepper
  class(timestepper_BE_type), pointer :: tran_timestepper
  
  units_category = 'not_assigned'

  realization => simulation%realization
  patch => realization%patch
  
  if (associated(patch)) grid => patch%grid

  option => realization%option
  output_option => realization%output_option
  field => realization%field
  reaction => realization%reaction
  input => realization%input
  
  flow_timestepper => TimestepperBECreate()
  flow_timestepper%solver%itype = FLOW_CLASS
  tran_timestepper => TimestepperBECreate()
  tran_timestepper%solver%itype = TRANSPORT_CLASS

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
                              
  rewind(input%fid)  
  string = 'SUBSURFACE'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)  
      
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    option%io_buffer = 'pflotran card:: ' // trim(card)
    call printMsg(option)

    select case(trim(card))

!....................
      case ('GRID')
        call DiscretizationRead(realization%discretization,input,option)

!....................
      case ('CHEMISTRY')
        call ReactionReadPass2(reaction,input,option)

!....................
      case ('NONUNIFORM_VELOCITY')
        call InputReadNChars(input,option, &
                             realization%nonuniform_velocity_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','NONUNIFORM_VELOCITY')

      case ('UNIFORM_VELOCITY')
        uniform_velocity_dataset => UniformVelocityDatasetCreate()
        uniform_velocity_dataset%rank = 3
        uniform_velocity_dataset%interpolation_method = 1 ! 1 = STEP
        uniform_velocity_dataset%is_cyclic = PETSC_FALSE
        allocate(uniform_velocity_dataset%times(1))
        uniform_velocity_dataset%times = 0.d0
        allocate(uniform_velocity_dataset%values(3,1))
        uniform_velocity_dataset%values = 0.d0
        call InputReadDouble(input,option,uniform_velocity_dataset%values(1,1))
        call InputErrorMsg(input,option,'velx','UNIFORM_VELOCITY')
        call InputReadDouble(input,option,uniform_velocity_dataset%values(2,1))
        call InputErrorMsg(input,option,'vely','UNIFORM_VELOCITY')
        call InputReadDouble(input,option,uniform_velocity_dataset%values(3,1))
        call InputErrorMsg(input,option,'velz','UNIFORM_VELOCITY')
        ! read units, if present
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          units_category = 'length/time'
          units_conversion = UnitsConvertToInternal(word,units_category,option) 
          uniform_velocity_dataset%values(:,1) = &
            uniform_velocity_dataset%values(:,1) * units_conversion
        endif
        call UniformVelocityDatasetVerify(option,uniform_velocity_dataset)
        realization%uniform_velocity_dataset => uniform_velocity_dataset
      
      case ('VELOCITY_DATASET')
        uniform_velocity_dataset => UniformVelocityDatasetCreate()
        call UniformVelocityDatasetRead(uniform_velocity_dataset,input,option)
        realization%uniform_velocity_dataset => uniform_velocity_dataset

!....................
      case ('DEBUG')
        call DebugRead(realization%debug,input,option)
               
!....................
      case ('PRINT_PRIMAL_GRID')
        option%print_explicit_primal_grid = PETSC_TRUE
        
!....................
      case ('PRINT_DUAL_GRID')
        option%print_explicit_dual_grid = PETSC_TRUE

!....................
      case ('MAX_CHANGE')
        call InputReadDouble(input,option,option%dpmxe)
        call InputErrorMsg(input,option,'dpmxe','MAX_CHANGE')
        call InputReadDouble(input,option,option%dtmpmxe)
        call InputErrorMsg(input,option,'dtmpmxe','MAX_CHANGE')
        call InputReadDouble(input,option,option%dsmxe)
        call InputErrorMsg(input,option,'dsmxe','MAX_CHANGE')
        call InputReadDouble(input,option,option%dcmxe)
        call InputErrorMsg(input,option,'dcmxe','MAX_CHANGE')
        
!....................
      case ('PROC')
      
!....................
      case ('REGION')
        region => RegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','REGION') 
        call printMsg(option,region%name)
        call RegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call RegionAddToList(region,realization%region_list)   
        nullify(region)   

!....................
      case ('FLOW_CONDITION')
        flow_condition => FlowConditionCreate(option)
        call InputReadWord(input,option,flow_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'FLOW_CONDITION','name') 
        call printMsg(option,flow_condition%name)
        if (option%iflowmode == G_MODE) then
          call FlowConditionGeneralRead(flow_condition,input,option)
        else if(option%iflowmode == TOIL_IMS_MODE) then
          call FlowConditionTOilImsRead(flow_condition,input,option)
        else 
          call FlowConditionRead(flow_condition,input,option)
        endif
        call FlowConditionAddToList(flow_condition,realization%flow_conditions)
        nullify(flow_condition)
        
!....................
      case ('TRANSPORT_CONDITION')
        if (.not.associated(reaction)) then
          option%io_buffer = 'TRANSPORT_CONDITIONs not supported without ' // &
            'CHEMISTRY.'
          call printErrMsg(option)
        endif
        tran_condition => TranConditionCreate(option)
        call InputReadWord(input,option,tran_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'TRANSPORT_CONDITION','name') 
        call printMsg(option,tran_condition%name)
        call TranConditionRead(tran_condition,realization%transport_constraints, &
                               reaction,input,option)
        call TranConditionAddToList(tran_condition,realization%transport_conditions)
        nullify(tran_condition)

!....................
      case('CONSTRAINT')
        if (.not.associated(reaction)) then
          option%io_buffer = 'CONSTRAINTs not supported without CHEMISTRY.'
          call printErrMsg(option)
        endif
        tran_constraint => TranConstraintCreate(option)
        call InputReadWord(input,option,tran_constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'constraint','name') 
        call printMsg(option,tran_constraint%name)
        call TranConstraintRead(tran_constraint,reaction,input,option)
        call TranConstraintAddToList(tran_constraint,realization%transport_constraints)
        nullify(tran_constraint)


!....................
      case ('BOUNDARY_CONDITION')
        coupler => CouplerCreate(BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Boundary Condition name') 
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)
      
!....................
      case ('INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Initial Condition name') 
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)
      
!....................
      case ('SOURCE_SINK')
        coupler => CouplerCreate(SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name') 
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)        
      
!....................
      case ('SOURCE_SINK_SANDBOX')
        call SSSandboxInit(option)
        call SSSandboxRead(input,option)
      
!....................
      case ('FLOW_MASS_TRANSFER')
        flow_data_mediator => DataMediatorDatasetCreate()
        call InputReadWord(input,option,flow_data_mediator%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Flow Mass Transfer name') 
        call DataMediatorDatasetRead(flow_data_mediator,input,option)
        call flow_data_mediator%AddToList(realization%flow_data_mediator_list)
        nullify(flow_data_mediator)
      
!....................
      case ('RT_MASS_TRANSFER')
        rt_data_mediator => DataMediatorDatasetCreate()
        call InputReadWord(input,option,rt_data_mediator%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'RT Mass Transfer name')
        call DataMediatorDatasetRead(rt_data_mediator,input,option)
        call rt_data_mediator%AddToList(realization%tran_data_mediator_list)
        nullify(rt_data_mediator)
      
!....................
      case ('STRATIGRAPHY','STRATA')
        strata => StrataCreate()
        call StrataRead(strata,input,option)
        call RealizationAddStrata(realization,strata)
        nullify(strata)

!.....................
      case ('DATASET')
        nullify(dataset)
        call DatasetRead(input,dataset,option)
        call DatasetBaseAddToList(dataset,realization%datasets)
        nullify(dataset)
        
!....................

      case('REFERENCE_PRESSURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_pressure)
        call InputDefaultMsg(input,option,'Reference Pressure') 

!....................

      case('REFERENCE_DENSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_water_density)
        call InputDefaultMsg(input,option,'Reference Density') 

!....................

      case('MINIMUM_HYDROSTATIC_PRESSURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%minimum_hydrostatic_pressure)
        call InputDefaultMsg(input,option,'Minimum Hydrostatic Pressure') 

!......................

      case('REFERENCE_TEMPERATURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_temperature)
        call InputDefaultMsg(input,option,'Reference Temperature') 

!......................

      case('REFERENCE_POROSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_porosity)
        call InputDefaultMsg(input,option,'Reference Porosity') 

!......................

      case('REFERENCE_SATURATION')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_saturation)
        call InputDefaultMsg(input,option,'Reference Saturation') 

!......................

      case('NONISOTHERMAL')
        option%use_isothermal = PETSC_FALSE

!......................

      case('ISOTHERMAL')
        option%use_isothermal = PETSC_TRUE
        
!......................

      case('UPDATE_FLOW_PERMEABILITY')
        option%update_flow_perm = PETSC_TRUE
        
!......................

      case('DFN')
        grid%unstructured_grid%grid_type = TWO_DIM_GRID    
            
!......................

      case("MULTIPLE_CONTINUUM")
        option%use_mc = PETSC_TRUE
              
!......................

      case('SECONDARY_CONTINUUM_SOLVER')
        if (.not.option%use_mc) then
          option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can only be used ' // &
                             'with MULTIPLE_CONTINUUM keyword.'
          call printErrMsg(option)
        endif      
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('KEARST')
            option%secondary_continuum_solver = 1
          case('HINDMARSH')
            option%secondary_continuum_solver = 2
          case('THOMAS')
            option%secondary_continuum_solver = 3
          case default
            option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                               'HINDMARSH or KEARST. For single component'// &
                               'chemistry THOMAS can be used.'
          call printErrMsg(option)    
        end select        
!....................

      case('SECONDARY_CONSTRAINT')
        if (.not.option%use_mc) then
          option%io_buffer = 'SECONDARY_CONSTRAINT can only be used with ' // &
                             'MULTIPLE_CONTINUUM keyword.'
          call printErrMsg(option)
        endif
        if (.not.associated(reaction)) then
          option%io_buffer = 'SECONDARY_CONSTRAINT not supported without' // &
                             'CHEMISTRY.'
          call printErrMsg(option)
        endif
        sec_tran_constraint => TranConstraintCreate(option)
        call InputReadWord(input,option,sec_tran_constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'secondary constraint','name') 
        call printMsg(option,sec_tran_constraint%name)
        call TranConstraintRead(sec_tran_constraint,reaction,input,option)
        realization%sec_transport_constraint => sec_tran_constraint
        nullify(sec_tran_constraint)        

!......................

      case('BRIN','BRINE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%m_nacl)
        call InputDefaultMsg(input,option,'NaCl Concentration') 

        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word(1:len_trim(word)))
          case('MOLAL')
          case('MASS')
            option%m_nacl = option%m_nacl /FMWNACL/(1.D0-option%m_nacl)
          case('MOLE')    
            option%m_nacl = option%m_nacl /FMWH2O/(1.D0-option%m_nacl)
          case default
            print *, 'Wrong unit: ', word(1:len_trim(word))
            stop
         end select 
         if (OptionPrintToScreen(option)) print *, option%m_nacl
!......................

      case ('RESTART')
        option%restart_flag = PETSC_TRUE
        call InputReadNChars(input,option,option%restart_filename,MAXSTRINGLENGTH, &
                             PETSC_TRUE)
        call InputErrorMsg(input,option,'RESTART','Restart file name') 
        call InputReadDouble(input,option,option%restart_time)
        if (input%ierr == 0) then
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (input%ierr == 0) then
            units_category = 'time'
            option%restart_time = option%restart_time* &
                                  UnitsConvertToInternal(word,units_category,option)
          else
            call InputDefaultMsg(input,option,'RESTART, time units')
          endif
        endif

!......................

      case ('CHECKPOINT')
        option%checkpoint_flag = PETSC_TRUE
        call InputReadInt(input,option,option%checkpoint_frequency)

        if (input%ierr == 1) then
          option%checkpoint_frequency = 0
          do
            call InputReadPflotranString(input,option)
            call InputReadStringErrorMsg(input,option,card)
            if (InputCheckExit(input,option)) exit

            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','CHECKPOINT')
            call StringToUpper(word)

            select case(trim(word))
              case ('PERIODIC')
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'time increment', &
                                   'OUTPUT,PERIODIC')
                call StringToUpper(word)

                select case(trim(word))
                  case('TIME')
                    call InputReadDouble(input,option,temp_real)
                    call InputErrorMsg(input,option,'time increment', &
                                       'CHECKPOINT,PERIODIC,TIME')
                    call InputReadWord(input,option,word,PETSC_TRUE)
                    call InputErrorMsg(input,option,'time increment units', &
                                       'CHECKPOINT,PERIODIC,TIME')
                    units_category = 'time'
                    units_conversion = UnitsConvertToInternal(word,units_category,option)
                    output_option%periodic_checkpoint_time_incr = temp_real* &
                                                              units_conversion
                  case('TIMESTEP')
                    call InputReadInt(input,option,option%checkpoint_frequency)
                    call InputErrorMsg(input,option,'timestep increment', &
                                       'CHECKPOINT,PERIODIC,TIMESTEP')
                  case default
                    call InputKeywordUnrecognized(word,'CHECKPOINT,PERIODIC', &
                                                  option)
                end select

              case ('FORMAT')
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'format type', &
                                   'CHECKPOINT,FORMAT')
                call StringToUpper(word)
                select case(trim(word))
                  case('BINARY')
                    option%checkpoint_format_binary = PETSC_TRUE
                  case('HDF5')
                    option%checkpoint_format_hdf5 = PETSC_TRUE
                  case default
                    call InputKeywordUnrecognized(word,'CHECKPOINT,FORMAT', &
                                                  option)
                end select

              case default
                call InputKeywordUnrecognized(word,'CHECKPOINT',option)
            end select
          enddo
          if (output_option%periodic_checkpoint_time_incr /= 0.d0 .and. &
              option%checkpoint_frequency /= 0) then
            option%io_buffer = 'Both TIME and TIMESTEP cannot be specified ' // &
              'for CHECKPOINT,PERIODIC.'
            call printErrMsg(option)
          endif
          if (output_option%periodic_checkpoint_time_incr == 0.d0 .and. &
              option%checkpoint_frequency == 0) then
            option%io_buffer = 'Either, TIME and TIMESTEP need to be specified ' // &
              'for CHECKPOINT,PERIODIC.'
            call printErrMsg(option)
          endif
        endif

!......................

      case ('NUMERICAL_JACOBIAN_FLOW')
        option%numerical_derivatives_flow = PETSC_TRUE

!......................

      case ('NUMERICAL_JACOBIAN_RXN')
        option%numerical_derivatives_rxn = PETSC_TRUE

!......................

      case ('NUMERICAL_JACOBIAN_MULTI_COUPLE')
        option%numerical_derivatives_multi_coupling = PETSC_TRUE

!......................

      case ('COMPUTE_STATISTICS')
        option%compute_statistics = PETSC_TRUE

!....................

      case ('CO2_DATABASE')
        call InputReadNChars(input,option,option%co2_database_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'CO2_DATABASE','filename')
        
!....................

      case ('TIMESTEPPER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
             call flow_timestepper%ReadInput(input,option)
          case('TRAN','TRANSPORT')
            call tran_timestepper%ReadInput(input,option)
          case default
            option%io_buffer = 'TIMESTEPPER must specify FLOW or TRANSPORT.'
            call printErrMsg(option)
        end select

!....................

      case ('LINEAR_SOLVER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            call SolverReadLinear(flow_timestepper%solver,input,option)
          case('TRAN','TRANSPORT')
            call SolverReadLinear(tran_timestepper%solver,input,option)
          case default
            option%io_buffer = 'LINEAR_SOLVER must specify FLOW or TRANSPORT.'
            call printErrMsg(option)
        end select

!....................

      case ('NEWTON_SOLVER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            call SolverReadNewton(flow_timestepper%solver,input,option)
!TODO(geh): remove after 11/30/15 as inf_scaled_res_tol is no longer used
!            if (flow_timestepper%solver%check_post_convergence) then
!              option%flow%check_post_convergence = PETSC_TRUE
!            option%flow%inf_scaled_res_tol = &
!              flow_timestepper%solver%newton_inf_scaled_res_tol
!              option%flow%inf_rel_update_tol = &
!                flow_timestepper%solver%newton_inf_rel_update_tol
!            endif
          case('TRAN','TRANSPORT')
            call SolverReadNewton(tran_timestepper%solver,input,option)
            if (tran_timestepper%solver%check_post_convergence) then
              option%transport%check_post_convergence = PETSC_TRUE
              option%transport%inf_scaled_res_tol = &
                tran_timestepper%solver%newton_inf_scaled_res_tol
              option%transport%inf_rel_update_tol = &
                tran_timestepper%solver%newton_inf_rel_update_tol
            endif
          case default
            option%io_buffer = 'NEWTON_SOLVER must specify FLOW or TRANSPORT.'
            call printErrMsg(option)
        end select
!....................

      case ('FLUID_PROPERTY')

        fluid_property => FluidPropertyCreate()
        call FluidPropertyRead(fluid_property,input,option)
        call FluidPropertyAddToList(fluid_property,realization%fluid_properties)
        nullify(fluid_property)
        
!....................

      case ('EOS')
        call EOSRead(input,option)

!....................

      case ('SATURATION_FUNCTION')
        if (option%iflowmode == RICHARDS_MODE .or. &
            option%iflowmode == TOIL_IMS_MODE .or. &
            option%iflowmode == G_MODE) then
          option%io_buffer = &
            'Must compile with legacy_saturation_function=1 ' //&
            'to use the SATURATION_FUNCTION keyword.  Otherwise, use ' // &
            'CHARACTERISTIC_CURVES.'
          call printErrMsg(option)
        endif
        saturation_function => SaturationFunctionCreate(option)
        call InputReadWord(input,option,saturation_function%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','SATURATION_FUNCTION')
        call SaturationFunctionRead(saturation_function,input,option)
        call SatFunctionComputePolynomial(option,saturation_function)
        call PermFunctionComputePolynomial(option,saturation_function)
        call SaturationFunctionAddToList(saturation_function, &
                                         realization%saturation_functions)
        nullify(saturation_function)   

!....................

      case ('CHARACTERISTIC_CURVES')
      
        if (.not.(option%iflowmode == NULL_MODE .or. &
                  option%iflowmode == RICHARDS_MODE .or. &
                  option%iflowmode == TOIL_IMS_MODE .or. &
                  option%iflowmode == G_MODE)) then
          option%io_buffer = 'CHARACTERISTIC_CURVES not supported in flow ' // &
            'modes other than RICHARDS, TOIL_IMS,  or GENERAL.  Use ' // &
            'SATURATION_FUNCTION.'
          call printErrMsg(option)
        endif
        characteristic_curves => CharacteristicCurvesCreate()
        call InputReadWord(input,option,characteristic_curves%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','CHARACTERISTIC_CURVES')
        call CharacteristicCurvesRead(characteristic_curves,input,option)
!        call SatFunctionComputePolynomial(option,saturation_function)
!        call PermFunctionComputePolynomial(option,saturation_function)
        call CharacteristicCurvesAddToList(characteristic_curves, &
                                          realization%characteristic_curves)
        nullify(characteristic_curves)   

!....................
      
      case ('MATERIAL_PROPERTY')

        material_property => MaterialPropertyCreate()
        call InputReadWord(input,option,material_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','MATERIAL_PROPERTY')        
        call MaterialPropertyRead(material_property,input,option)
        call MaterialPropertyAddToList(material_property,realization%material_properties)
        nullify(material_property)

!....................

      case ('USE_TOUCH_OPTIONS')
        option%use_touch_options = PETSC_TRUE

      case ('MPI_IO')
!        call PetscOptionsInsertString('-viewer_binary_mpiio')

      case ('HANDSHAKE_IO')
        call InputReadInt(input,option,option%io_handshake_buffer_size)
        call InputErrorMsg(input,option,'io_handshake_buffer_size','HANDSHAKE_IO')

      case ('OVERWRITE_RESTART_TRANSPORT')
        option%overwrite_restart_transport = PETSC_TRUE

      case ('OVERWRITE_RESTART_FLOW_PARAMS')
        option%overwrite_restart_flow = PETSC_TRUE

      case ('INITIALIZE_FLOW_FROM_FILE')
        call InputReadNChars(input,option,option%initialize_flow_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','INITIALIZE_FLOW_FROM_FILE') 

      case ('INITIALIZE_TRANSPORT_FROM_FILE')
        call InputReadNChars(input,option,option%initialize_transport_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','INITIALIZE_TRANSPORT_FROM_FILE') 

      case ('CENTRAL_DIFFERENCE')
        option%use_upwinding = PETSC_FALSE

!....................
      case ('OBSERVATION')
        observation => ObservationCreate()
        call ObservationRead(observation,input,option)
        call ObservationAddToList(observation, &
                                  realization%patch%observation_list)
        nullify(observation)
      
!....................
      case ('INTEGRAL_FLUX')
        integral_flux => IntegralFluxCreate()
        call InputReadWord(input,option,integral_flux%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Integral Flux name') 
        call IntegralFluxRead(integral_flux,input,option)
        call IntegralFluxAddToList(integral_flux, &
                                   realization%patch%integral_flux_list)
        nullify(integral_flux)
      
!.....................
      case ('WALLCLOCK_STOP')
        option%wallclock_stop_flag = PETSC_TRUE
        call InputReadDouble(input,option,option%wallclock_stop_time)
        call InputErrorMsg(input,option,'stop time','WALLCLOCK_STOP') 

        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr /= 0) word = 'h'
        call InputDefaultMsg(input,option,'WALLCLOCK_STOP time units')
        units_category = 'time'
        units_conversion = UnitsConvertToInternal(word,units_category,option) 
        ! convert from hrs to seconds and add to start_time
        option%wallclock_stop_time = option%start_time + &
                                     option%wallclock_stop_time*units_conversion
      
!....................
      case ('OUTPUT')
        vel_cent = PETSC_FALSE
        vel_face = PETSC_FALSE
        fluxes = PETSC_FALSE
        mass_flowrate = PETSC_FALSE
        energy_flowrate = PETSC_FALSE
        aveg_mass_flowrate = PETSC_FALSE
        aveg_energy_flowrate = PETSC_FALSE
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','OUTPUT') 
          call StringToUpper(word)
          select case(trim(word))
            case('TIME_UNITS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Output Time Units','OUTPUT')
              realization%output_option%tunit = trim(word)
              units_category = 'time'
              realization%output_option%tconv = &
                UnitsConvertToInternal(word,units_category,option)
            case('NO_FINAL','NO_PRINT_FINAL')
              output_option%print_final = PETSC_FALSE
            case('NO_INITIAL','NO_PRINT_INITIAL')
              output_option%print_initial = PETSC_FALSE
            case('PROCESSOR_ID')
              option%io_buffer = 'PROCESSOR_ID output must now be entered under OUTPUT/VARIABLES card as PROCESS_ID.'
              call printErrMsg(option)
!              output_option%print_iproc = PETSC_TRUE
            case('PERMEABILITY')
              option%io_buffer = 'PERMEABILITY output must now be entered under OUTPUT/VARIABLES card.'
              call printErrMsg(option)
!              output_option%print_permeability = PETSC_TRUE
            case('POROSITY')
              option%io_buffer = 'POROSITY output must now be entered under OUTPUT/VARIABLES card.'
              call printErrMsg(option)
!              output_option%print_porosity = PETSC_TRUE
            case('TORTUOSITY')
              option%io_buffer = 'TORTUOSITY output must now be entered under OUTPUT/VARIABLES card.'
              call printErrMsg(option)
!              output_option%print_tortuosity = PETSC_TRUE
            case('VOLUME')
              option%io_buffer = 'VOLUME output must now be entered under OUTPUT/VARIABLES card.'
              call printErrMsg(option)
!              output_option%print_volume = PETSC_TRUE
            case('MASS_BALANCE')
              option%compute_mass_balance_new = PETSC_TRUE
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputDefaultMsg(input,option, &
                                   'OUTPUT,MASS_BALANCE,DETAILED')
              if (len_trim(word) > 0) then
                call StringToUpper(word)
                select case(trim(word))
                  case('DETAILED')
                    option%mass_bal_detailed = PETSC_TRUE
                  case default
                    call InputKeywordUnrecognized(word, &
                           'OUTPUT,MASS_BALANCE',option)
                end select
              endif
            case('PRINT_COLUMN_IDS')
              output_option%print_column_ids = PETSC_TRUE
            case('TIMES')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'units','OUTPUT')
              units_category = 'time'
              units_conversion = UnitsConvertToInternal(word,units_category,option) 
              continuation_flag = PETSC_TRUE
              do
                continuation_flag = PETSC_FALSE
                if (index(input%buf,backslash) > 0) &
                  continuation_flag = PETSC_TRUE
                input%ierr = 0
                do
                  if (InputError(input)) exit
                  call InputReadDouble(input,option,temp_real)
                  if (.not.InputError(input)) then
                    waypoint => WaypointCreate()
                    waypoint%time = temp_real*units_conversion
                    waypoint%print_output = PETSC_TRUE    
                    call WaypointInsertInList(waypoint,realization%waypoint_list)
                  endif
                enddo
                if (.not.continuation_flag) exit
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
              enddo
            case('OUTPUT_FILE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'OUTPUT,OUTPUT_FILE')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_file = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%output_file_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,OUTPUT_FILE,PERIODIC')
                case default
                  call InputKeywordUnrecognized(word, &
                         'OUTPUT,OUTPUT_FILE',option)
              end select
            case('SCREEN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment','OUTPUT,SCREEN')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_screen = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%screen_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC,SCREEN')
                case default
                  call InputKeywordUnrecognized(word, &
                         'OUTPUT,SCREEN',option)
              end select
            case('PERIODIC')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'OUTPUT,PERIODIC')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'OUTPUT,PERIODIC,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'OUTPUT,PERIODIC,TIME')
                  units_category = 'time'
                  units_conversion = UnitsConvertToInternal(word,units_category,option) 
                  output_option%periodic_output_time_incr = temp_real* &
                                                            units_conversion
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  if (input%ierr == 0) then
                    if (StringCompareIgnoreCase(word,'between')) then

                      call InputReadDouble(input,option,temp_real)
                      call InputErrorMsg(input,option,'start time', &
                                         'OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'start time units', &
                                         'OUTPUT,PERIODIC,TIME')
                      units_category = 'time'
                      units_conversion = UnitsConvertToInternal(word,units_category,option) 
                      temp_real = temp_real * units_conversion
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      if (.not.StringCompareIgnoreCase(word,'and')) then
                        input%ierr = 1
                      endif
                      call InputErrorMsg(input,option,'and', &
                                          'OUTPUT,PERIODIC,TIME"')
                      call InputReadDouble(input,option,temp_real2)
                      call InputErrorMsg(input,option,'end time', &
                                         'OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'end time units', &
                                         'OUTPUT,PERIODIC,TIME')
                      temp_real2 = temp_real2 * units_conversion
                      do
                        waypoint => WaypointCreate()
                        waypoint%time = temp_real
                        waypoint%print_output = PETSC_TRUE    
                        call WaypointInsertInList(waypoint,realization%waypoint_list)
                        temp_real = temp_real + output_option%periodic_output_time_incr
                        if (temp_real > temp_real2) exit
                      enddo
                      output_option%periodic_output_time_incr = 0.d0
                    else
                      input%ierr = 1
                      call InputErrorMsg(input,option,'between', &
                                          'OUTPUT,PERIODIC,TIME')
                    endif
                  endif                  
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC,TIMESTEP')
                case default
                  call InputKeywordUnrecognized(word, &
                         'OUTPUT,PERIODIC',option)
              end select
            case('PERIODIC_OBSERVATION')
              output_option%print_observation = PETSC_TRUE
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                'OUTPUT, PERIODIC_OBSERVATION')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIME')
                  units_category = 'time'
                  units_conversion = UnitsConvertToInternal(word,units_category,option) 
                  output_option%periodic_tr_output_time_incr = temp_real* &
                                                               units_conversion
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_tr_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIMESTEP')
                case default
                  call InputKeywordUnrecognized(word, &
                         'OUTPUT,PERIODIC_OBSERVATION',option)
              end select
            case('FORMAT')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'keyword','OUTPUT,FORMAT') 
              call StringToUpper(word)
              select case(trim(word))
                case ('HDF5')
                  output_option%print_hdf5 = PETSC_TRUE
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputDefaultMsg(input,option, &
                                       'OUTPUT,FORMAT,HDF5,# FILES')
                  if (len_trim(word) > 0) then
                    call StringToUpper(word)
                    select case(trim(word))
                      case('SINGLE_FILE')
                        output_option%print_single_h5_file = PETSC_TRUE
                      case('MULTIPLE_FILES')
                        output_option%print_single_h5_file = PETSC_FALSE
                        output_option%times_per_h5_file = 1
                        call InputReadWord(input,option,word,PETSC_TRUE)
                        if (len_trim(word)>0) then
                          select case(trim(word))
                            case('TIMES_PER_FILE')
                              call InputReadInt(input,option, &
                                              output_option%times_per_h5_file)
                              call InputErrorMsg(input,option, &
                                'timestep increment', &
                                'OUTPUT,FORMAT,HDF5,MULTIPLE_FILES,TIMES_PER_FILE')
                            case default
                              call InputKeywordUnrecognized(word, &
                                    'OUTPUT,FORMAT,HDF5,MULTIPLE_FILES',option)
                          end select
                        endif
                      case default
                        call InputKeywordUnrecognized(word, &
                               'OUTPUT,FORMAT,HDF5',option)
                    end select
                  endif
                case ('MAD')
                  output_option%print_mad = PETSC_TRUE
                case ('TECPLOT')
                  output_option%print_tecplot = PETSC_TRUE
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'TECPLOT','OUTPUT,FORMAT') 
                  call StringToUpper(word)
                  select case(trim(word))
                    case('POINT')
                      output_option%tecplot_format = TECPLOT_POINT_FORMAT
                    case('BLOCK')
                      output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
                    case('FEBRICK')
                      output_option%tecplot_format = TECPLOT_FEBRICK_FORMAT
                    case default
                      call InputKeywordUnrecognized(word, &
                               'OUTPUT,FORMAT,TECPLOT',option)
                  end select
                  if (output_option%tecplot_format == TECPLOT_POINT_FORMAT &
                      .and. option%mycommsize > 1) then
                    output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
                  endif
                  if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
                    output_option%tecplot_format = TECPLOT_FEBRICK_FORMAT
                  endif
                case ('VTK')
                  output_option%print_vtk = PETSC_TRUE
                case default
                  call InputKeywordUnrecognized(word,'OUTPUT,FORMAT',option)
              end select
            case('VELOCITY_AT_CENTER')
              vel_cent = PETSC_TRUE
            case('VELOCITY_AT_FACE')
              vel_face = PETSC_TRUE
            case('FLUXES')
              fluxes = PETSC_TRUE
            case('FLOWRATES','FLOWRATE')
              mass_flowrate = PETSC_TRUE
              energy_flowrate = PETSC_TRUE
            case('MASS_FLOWRATE')
              mass_flowrate = PETSC_TRUE
            case('ENERGY_FLOWRATE')
              energy_flowrate = PETSC_TRUE
            case('AVERAGE_FLOWRATES','AVERAGE_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
              aveg_energy_flowrate = PETSC_TRUE
            case('AVERAGE_MASS_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
            case('AVERAGE_ENERGY_FLOWRATE')
              aveg_energy_flowrate = PETSC_TRUE
            case ('HDF5_WRITE_GROUP_SIZE')
              call InputReadInt(input,option,option%hdf5_write_group_size)
              call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')
            case('VARIABLES')
              call OutputVariableRead(input,option,output_option%output_variable_list)
            case('AVERAGE_VARIABLES')
              call OutputVariableRead(input,option,output_option%aveg_output_variable_list)
            case('UNFILTER_NON_STATE_VARIABLES')
              output_option%filter_non_state_variables = PETSC_FALSE
            case default
              call InputKeywordUnrecognized(word,'OUTPUT',option)
          end select
        enddo
        if (vel_cent) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_vel_cent = PETSC_TRUE
          if (output_option%print_hdf5) &
            output_option%print_hdf5_vel_cent = PETSC_TRUE
          if (output_option%print_vtk) &
            output_option%print_vtk_vel_cent = PETSC_TRUE
        endif
        if (vel_face) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_vel_face = PETSC_TRUE
          if (output_option%print_hdf5) &
           output_option%print_hdf5_vel_face = PETSC_TRUE
        endif
        if (fluxes) then
          output_option%print_fluxes = PETSC_TRUE
        endif
        if(output_option%aveg_output_variable_list%nvars>0) then
          if(output_option%periodic_output_time_incr==0.d0) then
            option%io_buffer = 'Keyword: AVERAGE_VARIABLES defined without' // &
                               ' PERIODIC TIME being set.'
            call printErrMsg(option)
          endif
          if(.not.output_option%print_hdf5) then
            option%io_buffer = 'Keyword: AVERAGE_VARIABLES only defined for FORMAT HDF5'
            call printErrMsg(option)
          endif
        endif
        if (mass_flowrate.or.energy_flowrate.or.aveg_mass_flowrate.or.aveg_energy_flowrate) then
          if (output_option%print_hdf5) then
            output_option%print_hdf5_mass_flowrate = mass_flowrate
            output_option%print_hdf5_energy_flowrate = energy_flowrate
            output_option%print_hdf5_aveg_mass_flowrate = aveg_mass_flowrate
            output_option%print_hdf5_aveg_energy_flowrate = aveg_energy_flowrate
            if(aveg_mass_flowrate.or.aveg_energy_flowrate) then
              if(output_option%periodic_output_time_incr==0.d0) then
                option%io_buffer = 'Keyword: AVEGRAGE_FLOWRATES/ ' // &
                  'AVEGRAGE_MASS_FLOWRATE/ENERGY_FLOWRATE defined without' // &
                  ' PERIODIC TIME being set.'
                call printErrMsg(option)
              endif
            endif
           option%flow%store_fluxes = PETSC_TRUE
          endif
          if (associated(grid%unstructured_grid%explicit_grid)) then
           option%flow%store_fluxes = PETSC_TRUE
            output_option%print_explicit_flowrate = mass_flowrate
          endif
        
        endif

!.....................
      case ('REGRESSION')
        call RegressionRead(simulation%regression,input,option)

!.....................
      case ('TIME')
!        dt_init = UNINITIALIZED_DOUBLE
        dt_init = 1.d0
        dt_min = UNINITIALIZED_DOUBLE
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'word','TIME') 
          select case(trim(word))
            case('STEADY_STATE')
              option%steady_state = PETSC_TRUE
            case('FINAL_TIME')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Final Time','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Final Time Units','TIME')
              units_category = 'time'
              temp_real2 = UnitsConvertToInternal(word,units_category,option)
              if (len_trim(realization%output_option%tunit) == 0) then
                realization%output_option%tunit = trim(word)
                realization%output_option%tconv = temp_real2
              endif
              waypoint => WaypointCreate()
              waypoint%final = PETSC_TRUE
              waypoint%time = temp_real*temp_real2
              waypoint%print_output = PETSC_TRUE              
              call WaypointInsertInList(waypoint,realization%waypoint_list)
            case('INITIAL_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Initial Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Initial Timestep Size Time Units','TIME')
              ! convert to internal units
              units_category = 'time'
              dt_init = temp_real*UnitsConvertToInternal(word,units_category,option)
            case('MINIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Minimum Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Minimum Timestep Size Time Units','TIME')
              units_category = 'time'
              dt_min = temp_real*UnitsConvertToInternal(word,units_category,option)
            case('MAXIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Maximum Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Maximum Timestep Size Time Units','TIME')
              waypoint => WaypointCreate()
              units_category = 'time'
              waypoint%dt_max = temp_real*UnitsConvertToInternal(word,units_category,option)
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call StringToUpper(word)
                if (StringCompare(word,'AT',TWO_INTEGER)) then
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'Maximum Timestep Size Update Time','TIME') 
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'Maximum Timestep Size Update Time Units','TIME')
                  units_category = 'time'
                  waypoint%time = temp_real*UnitsConvertToInternal(word,units_category,option)
                else
                  option%io_buffer = 'Keyword under "MAXIMUM_TIMESTEP_SIZE" after ' // &
                                     'maximum timestep size should be "at".'
                  call printErrMsg(option)
                endif
              else
                waypoint%time = 0.d0
              endif     
              call WaypointInsertInList(waypoint,realization%waypoint_list)
            case default
              call InputKeywordUnrecognized(word,'TIME',option)
          end select
        enddo
        if (Initialized(dt_init)) then
          if (associated(flow_timestepper)) then
            flow_timestepper%dt_init = dt_init
            option%flow_dt = dt_init
          endif
          if (associated(tran_timestepper)) then
            tran_timestepper%dt_init = dt_init
            option%tran_dt = dt_init
          endif
        endif
        if (Initialized(dt_min)) then
          if (associated(flow_timestepper)) then
            flow_timestepper%dt_min = dt_min
          endif
          if (associated(tran_timestepper)) then
            tran_timestepper%dt_min = dt_min
          endif
        endif

!......................
      case ('HDF5_READ_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_read_group_size)
        call InputErrorMsg(input,option,'HDF5_READ_GROUP_SIZE','Group size')

!......................
      case ('HDF5_WRITE_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_write_group_size)
        call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')

!....................
      case('WIPP')
        wipp => WIPPGetPtr()
        call WIPPRead(input,option)
        
!....................
      case('KLINKENBERG_EFFECT')
        wipp => WIPPGetPtr()
        call KlinkenbergInit()
        klinkenberg => KlinkenbergCreate()
        call Klinkenberg%Read(input,option)
        
!....................
      case ('ONLY_VERTICAL_FLOW')
        option%flow%only_vertical_flow = PETSC_TRUE
        if (option%iflowmode /= TH_MODE .and. &
            option%iflowmode /= RICHARDS_MODE) then
          option%io_buffer = 'ONLY_VERTICAL_FLOW implemented in RICHARDS and TH mode.'
          call printErrMsg(option)
        endif

!....................
      case ('RELATIVE_PERMEABILITY_AVERAGE')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (trim(word))
          case ('UPWIND')
            option%rel_perm_aveg = UPWIND
          case ('HARMONIC')
            option%rel_perm_aveg = HARMONIC
          case ('DYNAMIC_HARMONIC')
            option%rel_perm_aveg = DYNAMIC_HARMONIC
          case default
            option%io_buffer = 'Cannot identify the specificed ' // &
              'RELATIVE_PERMEABILITY_AVERAGE.'
            call printErrMsg(option)
          end select
          
!....................
      case ('DBASE_FILENAME')

!....................
      case ('END_SUBSURFACE')
        exit

!....................
      case ('MIN_ALLOWABLE_SCALE')
        call InputReadDouble(input,option,option%min_allowable_scale)
        call InputErrorMsg(input,option,'minimium allowable scaling factor', &
                           'InitSubsurface')

!....................
      case default
        call InputKeywordUnrecognized(word,'InitSubsurfaceReadInput()',option)
    end select

  enddo

  if (associated(simulation%flow_process_model_coupler)) then
    flow_timestepper%name = 'FLOW'
    if (option%steady_state) call TimestepperSteadyCreateFromBE(flow_timestepper)
    simulation%flow_process_model_coupler%timestepper => flow_timestepper
  else
    call flow_timestepper%Destroy()
    deallocate(flow_timestepper)
    nullify(flow_timestepper)
  endif
  if (associated(simulation%rt_process_model_coupler)) then
    tran_timestepper%name = 'TRAN'
    if (option%steady_state) call TimestepperSteadyCreateFromBE(tran_timestepper)    
    simulation%rt_process_model_coupler%timestepper => tran_timestepper
  else
    call tran_timestepper%Destroy()
    deallocate(tran_timestepper)
    nullify(tran_timestepper)
  endif
                                      
end subroutine InitSubsurfaceReadInput

end module Init_Subsurface_module
