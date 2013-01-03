!=======================================================================
! PFLOTRAN v2.0 LA-CC-09-047
!=======================================================================

!Copyright 2009. Los Alamos National Security, LLC. This material was produced under U.S. 
!Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated 
!by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has 
!rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS 
!NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE 
!USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software 
!should be clearly marked, so as not to confuse it with the version available from LANL.
!Additionally, this library is free software; you can redistribute it and/or modify it under the 
!terms of the GNU Lesser General Public License as published by the Free Software Foundation; 
!either version 2.1 of the License, or (at your option) any later version. Accordingly, this 
!library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
!the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser 
!General Public License for more details.

! Send all bug reports/questions/comments to:
!
! Peter C. Lichtner
! Los Alamos National Laboratory
! Earth and Environmental Sciences
! EES-16, MS: D469
! (505) 667-3420
! lichtner@lanl.gov
! Los Alamos, NM

! or

! Glenn E. Hammond
! Pacific Northwest National Laboratory
! Energy and Environment Directorate
! MSIN K9-36
! (509) 375-3875
! glenn.hammond@pnnl.gov
! Richland, WA

!=======================================================================
program pflotran_rxn
  
  use Reaction_module
  use Reaction_Aux_module
  use Database_module
  use Option_module
  use Input_module
  
  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  PetscBool :: option_found  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: filename_out
  type(reaction_type), pointer :: reaction
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  PetscErrorCode :: ierr
  
  option => OptionCreate()
  option%fid_out = OUT_UNIT

  call MPI_Init(ierr)
  option%global_comm = MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD,option%global_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,option%global_commsize,ierr)
  call MPI_Comm_group(MPI_COMM_WORLD,option%global_group,ierr)
  option%mycomm = option%global_comm
  option%myrank = option%global_rank
  option%mycommsize = option%global_commsize
  option%mygroup = option%global_group

  ! check for non-default input filename
  option%input_filename = "pflotran.in"
  string = '-pflotranin'
  call InputGetCommandLineString(string, option%input_filename, option_found,option)

  string = '-output_prefix'
  call InputGetCommandLineString(string,option%global_prefix,option_found,option)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  input => InputCreate(IN_UNIT, option%input_filename, option)

  filename_out = trim(option%global_prefix) // trim(option%group_prefix) // &
                 '.out'

  if (option%myrank == option%io_rank .and. option%print_to_file) then
    open(option%fid_out, file=filename_out, action="write", status="unknown")
  endif

  ! first pass through the input file to check for a CHEMISTRY block
  string = "CHEMISTRY"
  call InputFindStringInFile(input,option,string)
  if (.not.InputError(input)) then
    ! found a chemistry block, initialize the chemistry
    reaction => ReactionCreate()
    call ReactionRead(reaction, input, option)
    reaction%primary_species_names => GetPrimarySpeciesNames(reaction)
    ! PCL add in colloid dofs
    option%ntrandof = GetPrimarySpeciesCount(reaction)
    option%ntrandof = option%ntrandof + GetColloidCount(reaction)
    reaction%ncomp = option%ntrandof
  endif
    
  if (associated(reaction)) then
    if (reaction%use_full_geochemistry) then
       call DatabaseRead(reaction, option)
       call BasisInit(reaction, option)    
    else
      ! NOTE(bja): do we need this for the batch chemistry driver?

      ! turn off activity coefficients since the database has not been read
      reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
      allocate(reaction%primary_species_print(option%ntrandof))
      reaction%primary_species_print = PETSC_TRUE
    endif
  endif

  ! the second pass through the input file to read the remaining blocks
  ! NOTE(bja) : how do we read chemistry info w/o requiring a full simualtion object?
  !call InitReadInput(simulation)


  ! cleanup
  call ReactionDestroy(reaction)
  call InputDestroy(input)
  call OptionDestroy(option)
  call PetscFinalize(ierr)
  call MPI_Finalize(ierr)

end program pflotran_rxn
