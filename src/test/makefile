include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

ifdef tau
#FC = /ccs/proj/geo002/tau/tau-2.19.2/craycnl/bin/tau_f90.sh -optCompInst -tau_makefile=$(TAU_MAKEFILE)
#FLINKER = /ccs/proj/geo002/tau/tau-2.19.2/craycnl/bin/tau_f90.sh -optCompInst -tau_makefile=$(TAU_MAKEFILE)
FC = /ccs/proj/geo002/tau/tau-2.19.2/craycnl/bin/tau_f90.sh -tau_makefile=$(TAU_MAKEFILE)
FLINKER = /ccs/proj/geo002/tau/tau-2.19.2/craycnl/bin/tau_f90.sh -tau_makefile=$(TAU_MAKEFILE)
endif

testc: testc.o
	${CLINKER} -o testc testc.o $(LIBS) $(PETSC_LIB)

testf90: testf90.o
	${FLINKER} -o testf90 testf90.o $(LIBS) $(PETSC_LIB)

# Should add this to default PETSc targets as well?
%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<
