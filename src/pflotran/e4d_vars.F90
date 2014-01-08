module vars

  implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscviewer.h"
#include "finclude/petscviewer.h90"  
#include "finclude/petscksp.h"
#include "finclude/petscksp.h90"

  integer, dimension(:), allocatable :: e4d_ranks,pf_e4d_ranks
  integer :: mpi_comm_grp,mpi_e4d_grp,mpi_pfe4d_grp,i
  integer :: my_wrank,my_pfe4d_rank,n_pfe4drank             !!my mpi rank
  integer :: tn_rank                                        !!number of processes 


  character*40 :: mshfile                                  !!file containing the mesh options
  character*40 :: efile                                    !!survey configuration file
  character*40 :: mapfile

  integer :: my_rank                                       !!my mpi rank
  integer :: ierr                                          !!generall error
  integer :: n_rank                                        !!number of processes 
  integer :: E4D_COMM                                      !!E4D_COMMUNICATOR
  integer :: PFE4D_MASTER_COMM                             !!PF to E4D MASTER COMMUNICATOR
  integer :: ne,tne                                        !!num my electrodes, total electrodes
  integer :: nsig                                          !!number of elements, sigma values
  integer :: nm                                            !!number of measurements
  integer :: nnodes                                        !!number of nodes
  integer :: nelem                                         !!number of elements (same as nsig)
  integer :: nfaces                                        !!number of faces
  integer :: nvals                                         !!number of non-zeros in coupling mat
  integer :: my_ne                                         !!number of electrodes I'm assigned
  integer :: nmy_drows                                     !!number of data in my assembly vector
  integer :: nmap

  integer, dimension(:,:), allocatable :: map_inds
  integer, dimension(:,:), allocatable :: s_conf           !!abmn survey configuration
  integer, dimension(:,:), allocatable :: eind             !!electrode assignments
  integer, dimension(:), allocatable :: nbounds,zones      !!node boundaries and element zones
  integer, dimension(:,:), allocatable :: elements         !!elements connections
  integer, dimension(:,:), allocatable :: faces            !!face connections
  integer, dimension(:), allocatable :: e_nods             !!indices of electrode nodes
  integer, dimension(:), allocatable :: rows,cols
  integer, dimension(:), allocatable :: trows,tcols
  integer, dimension(:), allocatable :: A_map              !!coupling matrix mapping vector
  integer, dimension(:), allocatable :: S_map              !!Sigma mapping vector
  integer, dimension(:), allocatable :: my_drows           !!rows of my data assemble vector

  real, dimension(:,:), allocatable :: e_pos
  real, dimension(:,:), allocatable :: nodes               !!node positions
  real, dimension(:,:), allocatable :: poles               !!pole solutions
  real, dimension(:), allocatable :: pf_sol                !!pflotran solution
  real, dimension(:), allocatable :: sigma                 !!element conductivities
  real, dimension(:), allocatable :: dpred                 !!simulated data vector
  real, dimension(:), allocatable :: my_dvals              !!values in my data assembly vector
  real, dimension(:), allocatable :: map

  PetscInt, dimension(:), allocatable :: d_nnz             !!petsc prealloc vec (diag blocks)
  PetscReal, dimension(:), allocatable :: delA
  Mat :: A,Ai
  PetscErrorCode :: perr
  MatType :: tp
  PetscInt :: prn(1),pcn(1)
  PetscReal :: val(1)
  PetscInt :: d_nz,o_nz
  Vec :: psol
  Vec :: X
  Vec :: B
  KSP :: KS,KSi
  PC :: P
  logical :: nzero_flag=.true.
  PetscScalar, pointer :: vloc(:)
  Vec :: pflotran_solution_vec_mpi
  Vec :: pflotran_solution_vec_seq
  VecScatter :: pflotran_scatter
  PetscInt :: pflotran_solution_vec_size
 
end module vars
