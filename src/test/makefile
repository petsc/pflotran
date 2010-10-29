# $Id: makefile,v 1.3 2004/07/31 01:16:44 lichtner Exp $

objdir = ./${PETSC_ARCH}/obj
srcdir = ./
pflotran_src = ./
common_src   = ./
bindir = ./${PETSC_ARCH}/bin

MYFLAGS          = 

CFLAGS		 = 
FFLAGS		 =
CPPFLAGS         = ${MYFLAGS}
FPPFLAGS         = ${MYFLAGS}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

petsc_test_f : petsc_test_f.o
	${FLINKER}   -o petsc_test_f petsc_test_f.o ${PETSC_LIB} ${LIBS} 

petsc_test_c: petsc_test_c.o
	${CLINKER}   -o petsc_test_c petsc_test_c.o ${PETSC_LIB} ${LIBS} 

# Should add this to default PETSc targets as well?
%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<
