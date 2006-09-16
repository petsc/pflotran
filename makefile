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
#MYFLAGS = ${PETSCFLAGS} -I${objdir} -DUSE_COMPRESSIBILITY 
#                                   -DHAVE_MPITOMPIZERO
LIBS   = 

include ${PETSC_DIR}/bmake/common/base

util_obj  = ${common_src}fileio.o ${common_src}water_eos.o \
            ${common_src}co2eos.o  ${common_src}gaseos_mod.o \
            ${common_src}co2_sw_rtsafe.o ${common_src}co2_span_wagner.o \
            ${common_src}oil_eos.o ${common_src}oil_pckr.o \
			${common_src}utilities.o \
	        ${common_src}ptran_global.o ${common_src}trdynmem.o

pflow_obj = ${common_src}pflow_vector_ops.o \
            ${common_src}pflow_read_gridsize.o \
            ${common_src}pflow_gridtype.o \
			${common_src}Readfield.o \
	        ${common_src}pflow_pckr_mod.o \
	        ${common_src}mixed_fluid_eos.o \
			${common_src}translator_mixed_fluid_mph.o \
			${common_src}translator_mixed_fluid_OWG.o \
		    ${common_src}hydrostat.o \
            ${common_src}pflow_LIQUID.o \
            ${common_src}pflow_COND.o \
            ${common_src}pflow_TH.o \
            ${common_src}pflow_THC.o \
            ${common_src}pflow_TWOPH.o \
            ${common_src}pflow_2pha_ResJac.o\
			${common_src}pflow_owg_ResJac.o\
			${common_src}pflow_mphase_ResJac.o\
			${common_src}pflow_solve.o \
            ${common_src}pflow_output.o \
			${common_src}pflowgrid_mod.o

ptran_obj = ${common_src}ptran_psi.o     ${common_src}ptran_dbase.o \
            ${common_src}trgamdh.o       ${common_src}ptran_speciation.o \
            ${common_src}ptran_setbnd.o  ${common_src}ptran_conn.o \
            ${common_src}ptran_destroy.o ${common_src}ptran_dt.o \
            ${common_src}ptran_init.o    ${common_src}ptran_out.o \
            ${common_src}ptran_read.o \
            ${common_src}trkinmin.o      ${common_src}trionexc.o \
            ${common_src}ptran_multi.o   ${common_src}ptran_update.o \
            ${common_src}trstdyst.o \
            ${common_src}ptran_solv.o
			
pflotran_obj = ${common_src}rock_react.o    ${common_src}pflotran_couple.o

pflow : $(util_obj) $(pflow_obj) ${pflow_src}pflow.o
	${FLINKER}   -o pflow $(util_obj) $(pflow_obj) ${pflow_src}pflow.o \
	${MYFLAGS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

ptran : $(util_obj) $(ptran_obj) ${ptran_src}ptran.o
	${FLINKER}   -o ptran $(util_obj) $(ptran_obj) ${ptran_src}ptran.o \
	${MYFLAGS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

pflotran: $(util_obj) $(ptran_obj) $(pflow_obj) $(pflotran_obj) \
	${pflotran_src}pflotran.o
	${FLINKER}   -o pflotran $(util_obj) $(ptran_obj) $(pflow_obj) \
	$(pflotran_obj) ${pflotran_src}pflotran.o \
	${MYFLAGS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

pflotran_fc: $(util_obj) $(ptran_obj) $(pflow_obj) $(pflotran_obj) \
	${pflotran_src}pflotran_fc.o
	${FLINKER}   -o pflotran_fc $(util_obj) $(ptran_obj) $(pflow_obj) \
	$(pflotran_obj) ${pflotran_src}pflotran_fc.o \
	${MYFLAGS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

%.o : %.F90
	${FC} ${FOPTFLAGS} $< ${MYFLAGS} -c ${PETSC_INCLUDE}
%.mod : %.F90
	${FLINKER} $< -c ${PETSC_INCLUDE} 

# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
pflow.o : pflowgrid_mod.o
pflowgrid_mod.o : pflow_gridtype.o water_eos.o pflow_vector_ops.o \
                  pflow_COND.o pflow_LIQUID.o pflow_TH.o pflow_THC.o \
                  pflow_TWOPH.o pflow_output.o ptran_global.o \
                  trdynmem.o fileio.o
pflow_vector_ops.o : water_eos.o
pflow_output.o : pflow_gridtype.o
pflow_COND.o : pflow_gridtype.o 
pflow_LIQUID.o : pflow_gridtype.o water_eos.o
pflow_TH.o : pflow_gridtype.o water_eos.o
pflow_THC.o : pflow_gridtype.o water_eos.o
pflow_TWOPH.o : pflow_gridtype.o water_eos.o
