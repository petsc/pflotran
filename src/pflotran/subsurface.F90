module Subsurface_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  public :: SubsurfInitMaterialProperties, &
            SubsurfAssignMatIDsToRegions, &
            SubsurfAssignMaterialProperties
  
contains

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
  
  type(realization_type) :: realization
  
  call SubsurfAllocMatPropDataStructs(realization)
  call SubsurfAssignMatIDsToRegions(realization)
  call SubsurfAssignMaterialProperties(realization)
  
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
  use Grid_module
  use Patch_module
  use Material_Aux_class
  
  implicit none
  
  type(realization_type) :: realization
  
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

end subroutine SubsurfAllocMatPropDataStructs

! ************************************************************************** !

subroutine SubsurfAssignMatIDsToRegions(realization)
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

  implicit none
  
  type(realization_type) :: realization
  
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
  
  option => realization%option

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    ! set material ids to uninitialized
    cur_patch%imat = UNINITIALIZED_INTEGER
    grid => cur_patch%grid
    strata => cur_patch%strata%first
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

end subroutine SubsurfAssignMatIDsToRegions

 ! ************************************************************************** !

subroutine SubsurfAssignMaterialProperties(realization)
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
  use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, PERMEABILITY_XY, &
                               PERMEABILITY_YZ, PERMEABILITY_XZ, &
                               TORTUOSITY, POROSITY
  use HDF5_module
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization
  
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
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(material_property_type), pointer :: material_property, null_material_property
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  class(material_auxvar_type), pointer :: material_auxvars(:)
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
    if (option%mimetic) then
      call VecGetArrayF90(field%perm0_xz,perm_xz_p,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(field%perm0_xy,perm_xy_p,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(field%perm0_yz,perm_yz_p,ierr);CHKERRQ(ierr)
    endif
  endif
  call VecGetArrayF90(field%porosity0,por0_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%tortuosity0,tor0_p,ierr);CHKERRQ(ierr)
        
  material_auxvars => patch%aux%Material%auxvars

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
      if (option%mimetic) then
        perm_xz_p(local_id) = material_property%permeability(1,3)
        perm_xy_p(local_id) = material_property%permeability(1,2)
        perm_yz_p(local_id) = material_property%permeability(2,3)
      endif
    endif
    if (associated(material_auxvars)) then
      call MaterialAssignPropertyToAux(material_auxvars(ghosted_id), &
                                        material_property,option)
    endif
    por0_p(local_id) = material_property%porosity
    tor0_p(local_id) = material_property%tortuosity
  enddo

  if (option%nflowdof > 0) then
    call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_xx,perm_xx_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_yy,perm_yy_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_zz,perm_zz_p,ierr);CHKERRQ(ierr)
    if (option%mimetic) then
      call VecRestoreArrayF90(field%perm0_xz,perm_xz_p,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(field%perm0_xy,perm_xy_p,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(field%perm0_yz,perm_yz_p,ierr);CHKERRQ(ierr)
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
      if (associated(material_property%porosity_dataset)) then
        string = ''
        string2 = 'Porosity'
        call HDF5ReadCellIndexedRealArray(realization,field%work, &
                    material_property%porosity_dataset%filename, &
                    string, &
                    string2, &
                    material_property%porosity_dataset%realization_dependent)
        call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
        call VecGetArrayF90(field%porosity0,por0_p,ierr);CHKERRQ(ierr)
        do local_id = 1, grid%nlmax
          if (patch%imat(grid%nL2G(local_id)) == &
              material_property%internal_id) then
            por0_p(local_id) = vec_p(local_id)
          endif
        enddo
        call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
        call VecRestoreArrayF90(field%porosity0,por0_p,ierr);CHKERRQ(ierr)
      endif
    endif
  enddo
      
  ! update ghosted values
  if (option%nflowdof > 0) then
    call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                                     field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_X,0)
    call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                                     field%work_loc,ONEDOF)  
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Y,0)
    call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                                     field%work_loc,ONEDOF)   
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Z,0)
    if (option%mimetic) then
      call DiscretizationGlobalToLocal(discretization,field%perm0_xz, &
                                       field%work_loc,ONEDOF)  
      call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_XZ,0)
      call DiscretizationGlobalToLocal(discretization,field%perm0_xy, &
                                       field%work_loc,ONEDOF)  
      call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_YZ,0)
      call DiscretizationGlobalToLocal(discretization,field%perm0_yz, &
                                       field%work_loc,ONEDOF)   
      call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_YZ,0)
    endif
    call DiscretizationLocalToLocal(discretization,field%icap_loc, &
                                    field%icap_loc,ONEDOF)   
    call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
                                    field%ithrm_loc,ONEDOF)
    call RealLocalToLocalWithArray(realization,SATURATION_FUNCTION_ID_ARRAY)
  endif
  
  call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                                   field%porosity_mnrl_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%porosity_mnrl_loc, &
                               POROSITY,0)
  call DiscretizationGlobalToLocal(discretization,field%tortuosity0, &
                                    field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               TORTUOSITY,0)
  ! rock properties
  do i = 1, max_material_index
    call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    do local_id = 1, patch%grid%nlmax
      ghosted_id = patch%grid%nL2G(local_id)
      vec_p(local_id) = &
        patch%aux%Material%auxvars(patch%grid%nL2G(local_id))% &
        soil_properties(i)
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
      patch%aux%Material%auxvars(ghosted_id)%soil_properties(i) = &
         vec_p(ghosted_id)
    enddo
    call VecRestoreArrayF90(field%work_loc,vec_p,ierr);CHKERRQ(ierr)
  enddo
  
end subroutine SubsurfAssignMaterialProperties

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

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(realization_type) :: realization
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
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization
  type(material_property_type) :: material_property

  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(discretization_type), pointer :: discretization
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: local_id, ghosted_id, natural_id
  PetscReal :: permeability
  PetscBool :: append_realization_id
  PetscInt :: fid = 86
  PetscInt :: status
  PetscInt :: idirection
  PetscInt :: temp_int
  PetscReal :: ratio, scale
  Vec :: global_vec
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: perm_xyz_p(:)

  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  discretization => realization%discretization

  call VecGetArrayF90(field%perm0_xx,perm_xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm0_yy,perm_yy_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%perm0_zz,perm_zz_p,ierr);CHKERRQ(ierr)
  
  if (index(material_property%permeability_dataset%filename,'.h5') > 0) then
    group_name = ''
    if (material_property%permeability_dataset%realization_dependent) then
      append_realization_id = PETSC_TRUE
    else
      append_realization_id = PETSC_FALSE
    endif

    call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                    option)
    if (material_property%isotropic_permeability .or. &
        (.not.material_property%isotropic_permeability .and. &
         material_property%vertical_anisotropy_ratio > 0.d0)) then
      dataset_name = 'Permeability'
      call HDF5ReadCellIndexedRealArray(realization,global_vec, &
                          material_property%permeability_dataset%filename, &
                          group_name,dataset_name,append_realization_id)
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
      if (grid%itype == STRUCTURED_GRID_MIMETIC) temp_int = YZ_DIRECTION
      do idirection = X_DIRECTION,temp_int
        select case(idirection)
          case(X_DIRECTION)
            dataset_name = 'PermeabilityX'
          case(Y_DIRECTION)
            dataset_name = 'PermeabilityY'
          case(Z_DIRECTION)
            dataset_name = 'PermeabilityZ'
          case(XY_DIRECTION)
            dataset_name = 'PermeabilityXY'
            call VecGetArrayF90(field%perm0_xy,perm_xyz_p,ierr);CHKERRQ(ierr)
          case(XZ_DIRECTION)
            dataset_name = 'PermeabilityXZ'
            call VecGetArrayF90(field%perm0_xz,perm_xyz_p,ierr);CHKERRQ(ierr)
          case(YZ_DIRECTION)
            dataset_name = 'PermeabilityYZ'
            call VecGetArrayF90(field%perm0_yz,perm_xyz_p,ierr);CHKERRQ(ierr)
        end select          
        call HDF5ReadCellIndexedRealArray(realization,global_vec, &
                                          material_property%permeability_dataset%filename, &
                                          group_name, &
                                          dataset_name,append_realization_id)
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
          case(XY_DIRECTION,XZ_DIRECTION,YZ_DIRECTION)
            do local_id = 1, grid%nlmax
              if (patch%imat(grid%nL2G(local_id)) == &
                  material_property%internal_id) then
                perm_xyz_p(local_id) = vec_p(local_id)
              endif
            enddo
            select case(idirection)
              case(XY_DIRECTION)
                call VecRestoreArrayF90(field%perm0_xy,perm_xyz_p, &
                                            ierr);CHKERRQ(ierr)
              case(XZ_DIRECTION)
                call VecRestoreArrayF90(field%perm0_xz,perm_xyz_p, &
                                            ierr);CHKERRQ(ierr)
              case(YZ_DIRECTION)
                call VecRestoreArrayF90(field%perm0_yz,perm_xyz_p, &
                                            ierr);CHKERRQ(ierr)
            end select
        end select
        call VecRestoreArrayF90(global_vec,vec_p,ierr);CHKERRQ(ierr)
      enddo
    endif
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  else

    call PetscLogEventBegin(logging%event_hash_map,ierr);CHKERRQ(ierr)
    call GridCreateNaturalToGhostedHash(grid,option)
    input => InputCreate(IUNIT_TEMP, &
                material_property%permeability_dataset%filename,option)
    do
      call InputReadPflotranString(input,option)
      if (InputError(input)) exit
      call InputReadInt(input,option,natural_id)
      call InputErrorMsg(input,option,'natural id','STRATA')
      ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (ghosted_id > 0) then
        if (patch%imat(ghosted_id) /= &
            material_property%internal_id) cycle
        local_id = grid%nG2L(ghosted_id)
        if (local_id > 0) then
          call InputReadDouble(input,option,permeability)
          call InputErrorMsg(input,option,'permeability','STRATA')
          perm_xx_p(local_id) = permeability
          perm_yy_p(local_id) = permeability
          perm_zz_p(local_id) = permeability
        endif
      endif
    enddo

    call InputDestroy(input)
    call GridDestroyHashTable(grid)
    call PetscLogEventEnd(logging%event_hash_map,ierr);CHKERRQ(ierr)
  endif
  
  call VecRestoreArrayF90(field%perm0_xx,perm_xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm0_yy,perm_yy_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%perm0_zz,perm_zz_p,ierr);CHKERRQ(ierr)
  
end subroutine SubsurfReadPermsFromFile

end module Subsurface_module
