#if (PETSC_VERSION_RELEASE == 1)

! Name differences between petsc-dev and petsc-release.
#define SNESMonitorSet SNESSetMonitor
#define MatCreateMFFD MatCreateSNESMF
#define MatMFFDSetType MatSNESMFSetType
#define MATMFFD_WP MATSNESMF_WP
#define MatMFFDGetH MatSNESMFGetH
#define MatMFFDSetHHistory MatSNESMFSetHHistory
#define IS_COLORING_GLOBAL IS_COLORING_LOCAL

! Wrappers to handle prototype differences between petsc-dev and petsc-release.
#define VecScatterBegin VecScatterBegin_wrap
#define VecScatterEnd VecScatterEnd_wrap

#endif
