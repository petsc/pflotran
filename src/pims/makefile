# $Id: makefile,v 1.3 2004/07/31 01:16:44 lichtner Exp $
F90=f90
FLINKER=f90

objdir = ./${PETSC_ARCH}/obj
srcdir = ./
pflotran_src = ./
pflow_src    = ./
ptran_src    = ./
common_src   = ./
bindir = ./${PETSC_ARCH}/bin

CFLAGS		 = 
FFLAGS		 = 
CPPFLAGS         =
FPPFLAGS         =
LOCDIR		 = /home/clu
#MYFLAGS = ${PETSCFLAGS} -I${objdir} -DDEBUG -DUSE_COMPRESSIBILITY
#MYFLAGS = ${PETSCFLAGS} -I${objdir} -DUSE_COMPRESSIBILITY -DUSE_PETSC216
MYFLAGS = ${PETSCFLAGS} -O3 -I${objdir} -DUSE_COMPRESSIBILITY -DUSE_PETSC221
MYFLAGS = ${PETSCFLAGS} -g -I${objdir} -DUSE_COMPRESSIBILITY -DUSE_PETSC221
#MYFLAGS = ${PETSCFLAGS} -I${objdir} -DUSE_COMPRESSIBILITY 
#                                   -DHAVE_MPITOMPIZERO
LIBS   = 

include ${PETSC_DIR}/bmake/common/base

util_obj  = ${common_src}fileio.o ${common_src}utilities.o \
            ${common_src}water_eos.o ${common_src}gaseos_mod.o 

pflow_obj = ${common_src}pflow_gridtype.o \
            ${common_src}pflow_read_input.o \
            ${common_src}Readfield.o \
	    ${common_src}translator_fluid_ims.o \
	    ${common_src}pflow_ims_ResJac_patch.o\
	    ${common_src}pflow_ims_ResJac.o\
	    ${common_src}pflow_solve.o \
            ${common_src}pflow_output.o \
	    ${common_src}pflow_step.o \
	    ${common_src}pflowgrid_mod.o


pims : $(util_obj) $(pflow_obj) ${pflow_src}pims.o stubs.o
	${FLINKER}   -o pims $(util_obj) $(pflow_obj) ${pflow_src}pims.o stubs.o\
	${MYFLAGS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

amrpims : $(util_obj) $(pflow_obj) ${pflow_src}pims.o

%.o : %.F90
	${FC} ${FOPTFLAGS} $< ${MYFLAGS} -c ${PETSC_INCLUDE}

%.mod : %.F90
	${FLINKER} $< -c ${PETSC_INCLUDE} 

# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
# ps.o : pflowgrid_mod.o
# pflowgrid_mod.o : pflow_gridtype.o water_eos.o fileio.o
# pflow_output.o : pflow_gridtype.o
