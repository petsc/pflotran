
program lookup_table_test
  

  use PFLOTRAN_Constants_module
  use Lookup_Table_module

  implicit none

#include "finclude/petscsys.h"

  class(lookup_table_uniform_type), pointer :: lookup_table
  class(lookup_table_general_type), pointer :: lookup_table_gen
  PetscReal :: pert

  ! 1D evenly gridded
  lookup_table => LookupTableCreateUniform(1)
  lookup_table%dims(1) = 5
  allocate(lookup_table%data(lookup_table%dims(1)))
  lookup_table%data = [500.d0, 400.d0, 450.d0, 550.d0, 605.d0]
  allocate(lookup_table%axis1%values(lookup_table%dims(1)))
  lookup_table%axis1%values = [-3.d0, -1.d0, 1.d0, 3.d0, 5.d0]

  print *
  print *, 'uniform'
  call LookupTableTest(lookup_table,-3.d0,500.d0)
  call LookupTableTest(lookup_table,-1.d0,400.d0)
  call LookupTableTest(lookup_table,1.d0,450.d0)
  call LookupTableTest(lookup_table,3.d0,550.d0)
  call LookupTableTest(lookup_table,5.d0,605.d0)
  
  ! 2.d-2 is 1 % of 2m distance
  pert = 2.d0/100.d0
  call LookupTableTest(lookup_table,-3.d0-pert,500.d0)
  call LookupTableTest(lookup_table,-3.d0+pert,499.d0)
  call LookupTableTest(lookup_table,-1.d0-pert,401.d0)
  call LookupTableTest(lookup_table,-1.d0+pert,400.5d0)
  call LookupTableTest(lookup_table,1.d0-pert,449.5d0)
  call LookupTableTest(lookup_table,1.d0+pert,451.d0)
  call LookupTableTest(lookup_table,3.d0-pert,549.d0)
  call LookupTableTest(lookup_table,3.d0+pert,550.55d0)
  call LookupTableTest(lookup_table,5.d0-pert,604.45d0)
  call LookupTableTest(lookup_table,5.d0+pert,605.d0)

  call LookupTableTest(lookup_table,0.75d0,443.75d0)
  call LookupTableTest(lookup_table,3.3d0,558.25d0)
  call LookupTableTest(lookup_table,-0.2d0,420.d0)
  call LookupTableTest(lookup_table,-3.2d0,500.d0)
  call LookupTableTest(lookup_table,1000.d0,605.d0)
 
  call LookupTableUniformDestroy(lookup_table)

  ! 1D general gridded
  lookup_table_gen => LookupTableCreateGeneral(1)
  lookup_table_gen%dims(1) = 5
  allocate(lookup_table_gen%data(lookup_table_gen%dims(1)))
  lookup_table_gen%data = [500.d0, 400.d0, 450.d0, 550.d0, 605.d0]
  allocate(lookup_table_gen%axis1%values(lookup_table_gen%dims(1)))
  lookup_table_gen%axis1%values = [-3.d0, -1.d0, 1.d0, 3.d0, 5.d0]  
  
  print *
  print *, 'general'
  call LookupTableTest(lookup_table_gen,-3.d0,500.d0)
  call LookupTableTest(lookup_table_gen,-1.d0,400.d0)
  call LookupTableTest(lookup_table_gen,1.d0,450.d0)
  call LookupTableTest(lookup_table_gen,3.d0,550.d0)
  call LookupTableTest(lookup_table_gen,5.d0,605.d0)
  
  ! 2.d-2 is 1 % of 2m distance
  pert = 2.d0/100.d0
  call LookupTableTest(lookup_table_gen,-3.d0-pert,500.d0)
  call LookupTableTest(lookup_table_gen,-3.d0+pert,499.d0)
  call LookupTableTest(lookup_table_gen,-1.d0-pert,401.d0)
  call LookupTableTest(lookup_table_gen,-1.d0+pert,400.5d0)
  call LookupTableTest(lookup_table_gen,1.d0-pert,449.5d0)
  call LookupTableTest(lookup_table_gen,1.d0+pert,451.d0)
  call LookupTableTest(lookup_table_gen,3.d0-pert,549.d0)
  call LookupTableTest(lookup_table_gen,3.d0+pert,550.55d0)
  call LookupTableTest(lookup_table_gen,5.d0-pert,604.45d0)
  call LookupTableTest(lookup_table_gen,5.d0+pert,605.d0)

  call LookupTableTest(lookup_table_gen,0.75d0,443.75d0)
  call LookupTableTest(lookup_table_gen,3.3d0,558.25d0)
  call LookupTableTest(lookup_table_gen,-0.2d0,420.d0)
  call LookupTableTest(lookup_table_gen,-3.2d0,500.d0)
  call LookupTableTest(lookup_table_gen,1000.d0,605.d0)
 
  call LookupTableGeneralDestroy(lookup_table_gen) 
  
  ! 2D evenly gridded
  lookup_table => LookupTableCreateUniform(2)
  lookup_table%dims(1) = 5
  lookup_table%dims(2) = 3
  allocate(lookup_table%data(lookup_table%dims(1)*lookup_table%dims(2)))
  lookup_table%data = [500.d0, 400.d0, 450.d0, 550.d0, 605.d0, &
                       -100.d0, 2.d0, 37.d0, 2.d0, 3.d0, &
                       -200.d0, 30.d0, 78.d0, 3.d0, 99.d0]
  allocate(lookup_table%axis1%values(lookup_table%dims(1)))
  lookup_table%axis1%values = [-3.d0, -1.d0, 1.d0, 3.d0, 5.d0]  
  allocate(lookup_table%axis2%values(lookup_table%dims(2)))
  lookup_table%axis2%values = [7.5d0, 10.d0, 12.5d0]  
  
  print *
  print *, 'uniform 2D'
  call LookupTableTest(lookup_table,-3.d0,7.5d0,500.d0)
  call LookupTableTest(lookup_table,-3.d0,8.5d0,260.d0)
  call LookupTableTest(lookup_table,-4.d0,8.5d0,260.d0)
  call LookupTableTest(lookup_table,-1.65d0,8.95d0,163.583d0)
  call LookupTableTest(lookup_table,-1.2d0,8.6d0,225.992d0)
  call LookupTableTest(lookup_table,1.25d0,11.d0,47.025d0)
  call LookupTableTest(lookup_table,2.5d0,12.05d0,19.77d0)
  call LookupTableTest(lookup_table,1.75d0,12.5d0,49.875d0)
  call LookupTableUniformDestroy(lookup_table) 
  
    ! 2D evenly gridded
  lookup_table_gen => LookupTableCreateGeneral(2)
  lookup_table_gen%dims(1) = 5
  lookup_table_gen%dims(2) = 3
  allocate(lookup_table_gen%data(lookup_table_gen%dims(1)*lookup_table_gen%dims(2)))
  lookup_table_gen%data = [500.d0, -100.d0, -200.d0, &
                           400.d0, 2.d0, 30.d0, &
                           450.d0, 37.d0, 78.d0, &
                           550.d0, 2.d0, 3.d0, &
                           605.d0, 3.d0, 99.d0]
  allocate(lookup_table_gen%axis1%values(lookup_table_gen%dims(1)))
  lookup_table_gen%axis1%values = [-3.d0, -1.d0, 1.d0, 3.d0, 5.d0]  
  allocate(lookup_table_gen%axis2%values(lookup_table_gen%dims(1)*lookup_table_gen%dims(2)))
  lookup_table_gen%axis2%values = [7.5d0, 10.d0, 12.5d0, &
                                   7.5d0, 10.d0, 14.5d0, &  ! note the 14.5
                                   7.5d0, 10.d0, 12.5d0, &  
                                   7.5d0, 10.d0, 12.5d0, &  
                                   -0.25d0, 10.d0, 12.5d0]  ! note the -0.25
                               
  print *
  print *, 'general 2D'
  call LookupTableTest(lookup_table_gen,-3.d0,7.5d0,500.d0)
  call LookupTableTest(lookup_table_gen,-3.d0,8.5d0,260.d0)
  call LookupTableTest(lookup_table_gen,-4.d0,8.5d0,260.d0)
  call LookupTableTest(lookup_table_gen,-1.65d0,8.95d0,163.583d0)

  call LookupTableTest(lookup_table_gen,-2.1d0,12.7d0,-101.54d0)
  call LookupTableTest(lookup_table_gen,4.9d0,-0.05d0,591.090975609756d0)

  call LookupTableTest(lookup_table_gen,-1.2d0,8.6d0,225.992d0)
  call LookupTableTest(lookup_table_gen,1.25d0,11.d0,47.025d0)
  call LookupTableTest(lookup_table_gen,2.5d0,12.05d0,19.77d0)
  call LookupTableTest(lookup_table_gen,1.75d0,12.5d0,49.875d0)
  call LookupTableGeneralDestroy(lookup_table_gen) 
  
end program lookup_table_test
