include ../../../src/pflotran/makefile
include ${TAO_DIR}/bmake/tao_common
VPATH=../../../src/pflotran
common_src=${VPATH}/
MYFLAGS = -I${objdir} $(HDF5_FLAGS) -I${common_src} -g


pflotran_tao: $(util_obj) $(pflotran_obj) pflotran_tao.o
	${FLINKER}   -o pflotran_tao $(util_obj) $(pflotran_obj) pflotran_tao.o  ${TAO_FORTRAN_LIB} ${TAO_LIB} -ltaopetsc ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${LIBS} 

generate_reference: $(util_obj) $(pflotran_obj) generate_reference.o
	${FLINKER}   -o generate_reference $(util_obj) $(pflotran_obj) generate_reference.o  ${TAO_FORTRAN_LIB} ${TAO_LIB} -ltaopetsc ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${LIBS} 


