# $Id: makefile,v 1.3 2004/07/31 01:16:44 lichtner Exp $

objdir = ./${PETSC_ARCH}/obj
srcdir = ./
pflotran_src = ./
pflow_src    = ./
ptran_src    = ./
common_src   = ./
bindir = ./${PETSC_ARCH}/bin

# To use the definitions below, set them equal to something on the command
# line, i.e.
#
# make pflow custom_convergence=1 hdf5=1
#
# Ensure that environment variables are defined (e.g. HDF5_INCLUDE, HDF5_LIB)
# It is probably best to define them in your .cshrc or .bashrc file.

ifdef custom_convergence
  CUSTOM_CONVERGENCE_FLAGS = -DCHUAN
endif
ifdef debug_convergence
  DEBUG_CONVERGENCE_FLAGS = -DGLENN
endif
ifdef hdf5
  HDF5_FLAGS = -DUSE_HDF5 -I$(HDF5_INCLUDE)
  HDF5_LIBS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5 -lz  
endif
MYFLAGS = -I${objdir} -DUSE_COMPRESSIBILITY $(CUSTOM_CONVERGENCE_FLAGS) \
          $(DEBUG_CONVERGENCE_FLAGS) $(HDF5_FLAGS)
LIBS = $(HDF5_LIBS)

CFLAGS		 = 
FFLAGS		 = #-Rb
CPPFLAGS         = ${MYFLAGS}
FPPFLAGS         = ${MYFLAGS}
CLEANFILES       = pflow ptran pflotran pflotran_fc

include ${PETSC_DIR}/bmake/common/base

util_obj  = ${common_src}fileio.o \
	${common_src}water_eos.o \
	${common_src}co2eos.o \
	${common_src}gaseos_mod.o \
	${common_src}co2_sw_rtsafe.o \
	${common_src}co2_span_wagner.o \
	${common_src}oil_eos.o \
	${common_src}oil_pckr.o \
	${common_src}utilities.o \
	${common_src}ptran_global.o \
	${common_src}trdynmem.o \
	${common_src}petscrelwrappers.o

pflow_obj = ${common_src}pflow_vector_ops.o \
	${common_src}pflow_read_gridsize.o \
	${common_src}pflow_gridtype.o \
	${common_src}hdf5_output.o \
	${common_src}Readfield.o \
	${common_src}pflow_pckr_mod.o \
	${common_src}mixed_fluid_eos.o \
	${common_src}translator_mixed_fluid_mph.o \
	${common_src}translator_mixed_fluid_OWG.o \
	${common_src}translator_mixed_fluid_vad.o \
	${common_src}translator_mixed_fluid_richards.o \
	${common_src}translator_mixed_fluid_flash.o \
	${common_src}condition.o \
	${common_src}pflow_vadose_ResJac.o\
	${common_src}unstructured_grid.o \
	${common_src}hydrostat.o \
	${common_src}nhydrostatic.o \
	${common_src}pflow_LIQUID.o \
	${common_src}pflow_COND.o \
	${common_src}pflow_TH.o \
	${common_src}pflow_THC.o \
	${common_src}pflow_2pha_ResJac.o\
	${common_src}pflow_owg_ResJac.o\
	${common_src}pflow_mphase_ResJac.o\
	${common_src}pflow_flash_ResJac.o\
	${common_src}pflow_richards_ResJac.o\
	${common_src}pflow_solve.o \
	${common_src}vectompizero.o \
	${common_src}pflow_output.o \
	${common_src}pflow_checkpoint.o \
	${common_src}pflow_output_new.o \
	${common_src}pflow_convergence.o \
	${common_src}pflowgrid_mod.o

ptran_obj = ${common_src}ptran_psi.o \
	${common_src}ptran_dbase.o \
	${common_src}trgamdh.o \
	${common_src}ptran_speciation.o \
	${common_src}ptran_setbnd.o \
	${common_src}ptran_conn.o \
	${common_src}ptran_destroy.o \
	${common_src}ptran_dt.o \
	${common_src}ptran_init.o \
	${common_src}ptran_out.o \
	${common_src}ptran_read.o \
	${common_src}trkinmin.o \
	${common_src}trionexc.o \
	${common_src}ptran_multi.o \
	${common_src}ptran_update.o \
	${common_src}trstdyst.o \
	${common_src}ptran_solv.o

pflotran_obj = ${common_src}rock_react.o \
	${common_src}pflotran_couple.o

pflow : $(util_obj) $(pflow_obj) ${pflow_src}pflow.o
	${FLINKER}   -o pflow $(util_obj) $(pflow_obj) ${pflow_src}pflow.o ${PETSC_LIB} $(LIBS)

ptran : $(util_obj) $(ptran_obj) ${ptran_src}ptran.o
	${FLINKER}   -o ptran $(util_obj) $(ptran_obj) ${ptran_src}ptran.o ${PETSC_LIB}

pflotran: $(util_obj) $(ptran_obj) $(pflow_obj) $(pflotran_obj) \
	${pflotran_src}pflotran.o
	${FLINKER}   -o pflotran $(util_obj) $(ptran_obj) $(pflow_obj) \
	$(pflotran_obj) ${pflotran_src}pflotran.o ${PETSC_LIB}

pflotran_fc: $(util_obj) $(ptran_obj) $(pflow_obj) $(pflotran_obj) \
	${pflotran_src}pflotran_fc.o
	${FLINKER}   -o pflotran_fc $(util_obj) $(ptran_obj) $(pflow_obj) \
	$(pflotran_obj) ${pflotran_src}pflotran_fc.o ${PETSC_LIB}

# Should add this to default PETSc targets as well?
%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<

# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
pflow.o : pflowgrid_mod.o
pflowgrid_mod.o : pflow_gridtype.o water_eos.o pflow_vector_ops.o \
                  pflow_COND.o pflow_LIQUID.o pflow_TH.o pflow_THC.o \
                  pflow_output.o ptran_global.o \
                  trdynmem.o fileio.o
pflow_vector_ops.o : water_eos.o
pflow_output.o : pflow_gridtype.o petscrelwrappers.o
pflow_output_new.o : pflow_gridtype.o petscrelwrappers.o hdf5_output.o
hdf5_output.o : pflow_gridtype.o 
pflow_COND.o : pflow_gridtype.o 
pflow_LIQUID.o : pflow_gridtype.o water_eos.o
pflow_TH.o : pflow_gridtype.o water_eos.o
pflow_THC.o : pflow_gridtype.o water_eos.o
