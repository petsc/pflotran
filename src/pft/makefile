objdir = ./${PETSC_ARCH}/obj
srcdir = ./
pft_src = ./
common_src   = ./
bindir = ./${PETSC_ARCH}/bin

#MYFLAGS = -I/usr/software/intel/compiler91/x86_64/include/c++ -I${objdir} 
MYFLAGS = -I${objdir} 
LIBS   = 

CFLAGS		 = 
FFLAGS		 = 
CPPFLAGS         = ${MYFLAGS}
FPPFLAGS         = ${MYFLAGS}
LOCDIR		 = /home/hammond/workspace/pft++
CLEANFILES       = pft

include ${PETSC_DIR}/bmake/common/base

pft_obj = ${common_src}Grid.o \
	${common_src}StructuredGrid.o \
	${common_src}GridConnection.o \
	${common_src}BoundaryCondition.o \
	${common_src}Source.o \
	${common_src}GridCell.o \
	${common_src}Output.o \
	${common_src}Flow.o \

pft: $(pft_obj) \
	${pft_src}Pft.o
	${CLINKER} -o pft ${pft_src}Pft.o ${pft_obj} \
	${PETSC_LIB}

