# $Id: makefile,v 1.3 2004/07/31 01:16:44 lichtner Exp $

objdir = ./${PETSC_ARCH}/obj
srcdir = ./
pflotran_src = ./
pflow_src    = ./
ptran_src    = ./
common_src   = ./
bindir = ./${PETSC_ARCH}/bin

ifdef hdf5
  ifndef jaguar
    HDF5_FLAGS = -I$(HDF5_INCLUDE) -I$(HDF5_LIB) -DUSE_HDF5
    HDF5_LIBS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5 -lz  
  else
# the below is for jaguar
    HDF5_FLAGS = $(HDF5_FLIB) -DUSE_HDF5
    HDF5_LIBS = $(HDF5_FLIB) -lhdf5_fortran -lhdf5 -lz
  endif
endif
ifdef compile_broken
  COMPILE_BROKEN_FLAGS = -DCOMPILE_BROKEN
endif

MYFLAGS = -I${objdir} $(HDF5_FLAGS) \
          $(COMPILE_BROKEN_FLAGS) \
          $(TRANSPORT_FLAGS)
LIBS = $(HDF5_LIBS)

CFLAGS		 = 
FFLAGS		 = #-Rb
CPPFLAGS         = ${MYFLAGS}
FPPFLAGS         = ${MYFLAGS}
ifdef ABSOFT
FPPFLAGS         += ${MYFLAGS} -p${HDF5_LIB}
endif

CLEANFILES       = pflow ptran pflotran pflotran_fc

ifdef petscrel
 include ${PETSC_DIR}/bmake/common/base
else
 include ${PETSC_DIR}/conf/base
endif

util_obj  = ${common_src}fileio.o \
	${common_src}water_eos.o \
	${common_src}co2eos.o \
	${common_src}gaseos_mod.o \
	${common_src}co2_sw_rtsafe.o \
	${common_src}co2_span_wagner.o \
	${common_src}oil_eos.o \
	${common_src}oil_pckr.o \
	${common_src}utility.o \
	${common_src}ptran_global.o \
	${common_src}trdynmem.o \
	${common_src}petscrelwrappers.o

pflow_obj = ${common_src}pflow_pckr_mod.o \
	${common_src}option.o \
	${common_src}reactive_transport_aux.o\
	${common_src}field.o \
	${common_src}debug.o \
	${common_src}waypoint.o \
	${common_src}region.o \
	${common_src}material.o \
	${common_src}strata.o \
	${common_src}connection.o \
	${common_src}breakthrough.o \
	${common_src}solver.o \
	${common_src}unstructured_grid.o \
	${common_src}structured_grid.o \
	${common_src}condition.o \
	${common_src}coupler.o \
	${common_src}grid.o \
	${common_src}hydrostatic.o \
	${common_src}realization.o \
	${common_src}hdf5.o \
	${common_src}mixed_fluid_eos.o \
	${common_src}translator_mixed_fluid_mph.o \
	${common_src}pflow_THC.o \
	${common_src}pflow_mphase_ResJac.o\
	${common_src}transport.o\
	${common_src}reaction.o\
	${common_src}reactive_transport.o\
	${common_src}richards.o\
	${common_src}richards_lite.o\
	${common_src}general_grid.o \
	${common_src}pflow_checkpoint.o \
	${common_src}convergence.o \
	${common_src}timestepper.o \
	${common_src}simulation.o \
	${common_src}pflow_output2.o \
	${common_src}pflow_init.o

ptran_obj = ${common_src}option.o \
	${common_src}debug.o \
	${common_src}waypoint.o \
	${common_src}region.o \
	${common_src}material.o \
	${common_src}strata.o \
	${common_src}connection.o \
	${common_src}breakthrough.o \
	${common_src}transport_solver.o \
	${common_src}unstructured_grid.o \
	${common_src}structured_grid.o \
	${common_src}condition.o \
	${common_src}coupler.o \
	${common_src}grid.o \
	${common_src}transport_realization.o \
	${common_src}hdf5.o \
	${common_src}transport_timestepper.o \
	${common_src}simulation.o \
	${common_src}ptran_output.o \
	${common_src}ptran_init.o

pflotran_obj = ${common_src}rock_react.o \
	${common_src}pflotran_couple.o

pflowamr: $(util_obj) $(pflow_obj)

pflow : $(util_obj) $(pflow_obj) ${pflow_src}pflow.o
	${FLINKER}   -o pflow $(util_obj) $(pflow_obj) ${pflow_src}pflow.o ${PETSC_LIB} ${LIBS}

ptran : TRANSPORT_FLAGS = -DTRANSPORT
ptran : $(util_obj) $(ptran_obj) ${ptran_src}ptran.o
	${FLINKER}   -o ptran $(util_obj) $(ptran_obj) ${ptran_src}ptran.o ${PETSC_LIB} ${LIBS}

pflotran: $(util_obj) $(ptran_obj) $(pflow_obj) $(pflotran_obj) \
	${pflotran_src}pflotran.o
	${FLINKER}   -o pflotran $(util_obj) $(ptran_obj) $(pflow_obj) \
	$(pflotran_obj) ${pflotran_src}pflotran.o ${PETSC_LIB}

pflotran_fc: $(util_obj) $(ptran_obj) $(pflow_obj) $(pflotran_obj) \
	${pflotran_src}pflotran_fc.o
	${FLINKER}   -o pflotran_fc $(util_obj) $(ptran_obj) $(pflow_obj) \
	$(pflotran_obj) ${pflotran_src}pflotran_fc.o ${PETSC_LIB}

ascii_to_hdf5: fileio.o hdf5_dataset.o ascii_to_hdf5.o 
	${FLINKER} -o ascii_to_hdf5 ascii_to_hdf5.o hdf5_dataset.o fileio.o \
	$(LIBS)

# Should add this to default PETSc targets as well?
%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<

# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
breakthrough.o : region.o option.o fileio.o
pflow.o : pflow_init.o timestepper.o simulation.o
ptran.o : ptran_init.o transport_timestepper.o simulation.o
#simulation.o : timestepper.o realization.o transport_realization.o \
#               transport_timestepper.o
pflow_init.o : option.o grid.o debug.o realization.o solver.o timestepper.o \
	simulation.o fileio.o material.o region.o condition.o convergence.o \
	waypoint.o field.o general_grid.o co2_span_wagner.o utility.o \
	coupler.o strata.o pflow_pckr_mod.o field.o pflow_mphase_ResJac.o \
	richards_lite.o richards.o
unstructured_grid.o : option.o connection.o solver.o
structured_grid.o : connection.o option.o fileio.o region.o solver.o
general_grid.o : grid.o realization.o hdf5.o
grid.o : structured_grid.o unstructured_grid.o option.o region.o connection.o \
	 solver.o
solver.o : option.o field.o
region.o : fileio.o utility.o
material.o : region.o option.o fileio.o field.o
coupler.o : condition.o connection.o region.o fileio.o grid.o
convergence.o : option.o solver.o
condition.o : fileio.o utility.o option.o field.o
realization.o : grid.o option.o region.o condition.o coupler.o material.o \
	     strata.o waypoint.o field.o hydrostatic.o debug.o breakthrough.o
transport_realization.o : grid.o option.o region.o condition.o coupler.o \
                          material.o strata.o waypoint.o debug.o breakthrough.o
strata.o : fileio.o region.o material.o
hydrostatic.o : option.o grid.o coupler.o condition.o connection.o region.o \
                structured_grid.o
waypoint.o : option.o 
pflow_vector_ops.o : water_eos.o
pflow_output2.o : grid.o structured_grid.o option.o realization.o field.o \
                  coupler.o connection.o hdf5.o breakthrough.o
pflow_mphase_ResJac.o : realization.o option.o field.o grid.o coupler.o \
                          translator_mixed_fluid_mph.o
reaction.o : reactive_transport_aux.o
transport.o : reactive_transport_aux.o
reactive_transport.o : reaction.o transport.o reactive_transport_aux.o
timestepper.o : \
	        mixed_fluid_eos.o \
                translator_mixed_fluid_mph.o \
                richards.o \
                richards_lite.o \
                pflow_mphase_ResJac.o \
                realization.o \
                grid.o \
                option.o \
                field.o \
                solver.o \
                pflow_output2.o \
                convergence.o \
                pflow_checkpoint.o \
		waypoint.o
transport_timestepper.o : \
                transport_realization.o \
                grid.o \
                option.o \
                field.o \
                transport_solver.o \
                ptran_output.o \
		waypoint.o
hdf5.o : realization.o option.o grid.o
pflow_THC.o : water_eos.o
pflow.o : simulation.o realization.o timestepper.o solver.o grid.o option.o \
	  pflow_init.o pflow_output2.o field.o
