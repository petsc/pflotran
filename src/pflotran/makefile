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
    HDF5_LIBS = $(HDF5_FLIB) -lhdf5_fortran -lhdf5 -lz -lsz
  endif
endif
ifdef compile_broken
  COMPILE_BROKEN_FLAGS = -DCOMPILE_BROKEN
endif

MYFLAGS = -I${objdir} $(HDF5_FLAGS) \
          $(COMPILE_BROKEN_FLAGS)
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
	${common_src}gaseos_mod.o \
	${common_src}co2eos.o \
	${common_src}co2_sw_rtsafe.o \
	${common_src}co2_span_wagner.o \
	${common_src}spline.o \
	${common_src}co2_span_wagner_spline.o \
	${common_src}co2_sw.o \
	${common_src}utility.o

pflotran_obj = ${common_src}option.o \
	${common_src}logging.o \
	${common_src}reactive_transport_aux.o\
	${common_src}reaction_aux.o\
	${common_src}thc_aux.o\
	${common_src}richards_aux.o\
	${common_src}mphase_aux.o\
	${common_src}immis_aux.o\
	${common_src}auxilliary.o\
	${common_src}field.o \
	${common_src}debug.o \
	${common_src}waypoint.o \
	${common_src}region.o \
	${common_src}material.o \
	${common_src}pckr_mod.o \
	${common_src}strata.o \
	${common_src}connection.o \
	${common_src}breakthrough.o \
	${common_src}solver.o \
	${common_src}cf90bridge.o\
	${common_src}unstructured_grid.o \
	${common_src}structured_grid.o \
	${common_src}condition.o \
	${common_src}coupler.o \
	${common_src}grid.o \
	${common_src}hydrostatic.o \
	${common_src}patch.o \
	${common_src}level.o \
	${common_src}amrgrid.o\
	${common_src}discretization.o \
	${common_src}database.o \
	${common_src}realization.o \
	${common_src}hdf5.o \
	${common_src}mphase.o\
	${common_src}immis.o\
	${common_src}transport.o\
	${common_src}reaction.o\
	${common_src}reactive_transport.o\
	${common_src}thc.o\
	${common_src}richards.o\
	${common_src}mass_balance.o\
	${common_src}general_grid.o \
	${common_src}checkpoint.o \
	${common_src}convergence.o \
	${common_src}timestepper.o \
	${common_src}simulation.o \
	${common_src}output.o \
	${common_src}init.o\
	${common_src}amrstubs.o


pflotranamr_obj = ${common_src}option.o \
	${common_src}logging.o \
	${common_src}reactive_transport_aux.o\
	${common_src}reaction_aux.o\
	${common_src}thc_aux.o\
	${common_src}richards_aux.o\
	${common_src}mphase_aux.o\
	${common_src}immis_aux.o\
	${common_src}auxilliary.o\
	${common_src}field.o \
	${common_src}debug.o \
	${common_src}waypoint.o \
	${common_src}region.o \
	${common_src}material.o \
        ${common_src}pckr_mod.o \
	${common_src}strata.o \
	${common_src}connection.o \
	${common_src}breakthrough.o \
	${common_src}cf90bridge.o \
	${common_src}solver.o \
	${common_src}unstructured_grid.o \
	${common_src}structured_grid.o \
	${common_src}condition.o \
	${common_src}coupler.o \
	${common_src}grid.o \
	${common_src}hydrostatic.o \
	${common_src}patch.o \
	${common_src}level.o \
	${common_src}amrgrid.o\
	${common_src}discretization.o \
	${common_src}database.o \
	${common_src}realization.o \
	${common_src}hdf5.o \
	${common_src}mphase.o\
	${common_src}immis.o\
	${common_src}transport.o\
	${common_src}reaction.o\
	${common_src}reactive_transport.o\
	${common_src}thc.o\
	${common_src}richards.o\
	${common_src}mass_balance.o\
	${common_src}general_grid.o \
	${common_src}checkpoint.o \
	${common_src}convergence.o \
	${common_src}timestepper.o \
	${common_src}simulation.o \
	${common_src}output.o \
	${common_src}init.o

pflotranamr: $(util_obj) $(pflotranamr_obj)
	ar ru libpflotranamr.a $(util_obj) $(pflotranamr_obj)

pflotran : $(util_obj) $(pflotran_obj) ${pflotran_src}pflotran.o
	${FLINKER}   -o pflotran $(util_obj) $(pflotran_obj) ${pflotran_src}pflotran.o ${PETSC_LIB} ${LIBS}

ascii_to_hdf5: fileio.o hdf5_dataset.o ascii_to_hdf5.o 
	${FLINKER} -o ascii_to_hdf5 ascii_to_hdf5.o hdf5_dataset.o fileio.o \
	$(LIBS)

test: functions.o test.o
	${FLINKER} -o test test.o functions.o $(LIBS) $(PETSC_LIB)

# Should add this to default PETSc targets as well?
%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<

# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
auxilliary.o : thc_aux.o richards_aux.o reactive_transport_aux.o \
	mphase_aux.o
pflotran.o : simulation.o output.o init.o timestepper.o logging.o
simulation.o : timestepper.o realization.o 
init.o : option.o grid.o debug.o realization.o solver.o timestepper.o \
	simulation.o fileio.o material.o region.o condition.o convergence.o \
	waypoint.o field.o general_grid.o co2_span_wagner.o utility.o \
	co2_span_wagner_spline.o co2_sw.o \
	coupler.o strata.o pckr_mod.o field.o mphase.o \
	richards.o thc.o discretization.o logging.o mass_balance.o \
	database.o reaction_aux.o
unstructured_grid.o : option.o connection.o 
structured_grid.o : connection.o option.o fileio.o 
general_grid.o : grid.o realization.o hdf5.o logging.o
discretization.o : grid.o structured_grid.o unstructured_grid.o option.o \
                   fileio.o
grid.o : structured_grid.o unstructured_grid.o option.o region.o connection.o
solver.o : option.o field.o
region.o : fileio.o utility.o logging.o
material.o : region.o option.o fileio.o field.o
coupler.o : condition.o connection.o region.o fileio.o grid.o
breakthrough.o : region.o option.o fileio.o
convergence.o : option.o solver.o
condition.o : fileio.o utility.o option.o field.o logging.o reaction_aux.o
realization.o : grid.o option.o region.o condition.o coupler.o material.o \
	     strata.o waypoint.o field.o hydrostatic.o debug.o breakthrough.o \
             patch.o level.o discretization.o reaction.o reaction_aux.o
level.o : patch.o option.o
patch.o : option.o grid.o coupler.o breakthrough.o strata.o region.o \
          auxilliary.o
strata.o : fileio.o region.o material.o
hydrostatic.o : option.o grid.o coupler.o condition.o connection.o region.o \
                structured_grid.o
waypoint.o : option.o 
output.o : grid.o structured_grid.o option.o realization.o field.o \
           coupler.o connection.o hdf5.o breakthrough.o \
					 reactive_transport.o \
           thc.o richards.o mphase.o \
					 discretization.o logging.o
reaction.o : reactive_transport_aux.o reaction_aux.o option.o
reaction_aux.o : option.o
transport.o : reactive_transport_aux.o option.o
reactive_transport.o : realization.o reaction.o transport.o \
                       reactive_transport_aux.o
reactive_transport_aux.o : option.o
thc.o : thc_aux.o realization.o
richards.o : richards_aux.o realization.o
thc_aux.o : option.o material.o water_eos.o
richards_aux.o : option.o material.o water_eos.o
mphase.o : mphase_aux.o realization.o 
mphase.o : immis_aux.o realization.o 
mphase_aux.o :  option.o material.o water_eos.o co2_span_wagner.o \
		co2eos.o co2_span_wagner_spline.o co2_sw.o pckr_mod.o
immis_aux.o :  option.o material.o water_eos.o co2_span_wagner.o \
		co2eos.o co2_span_wagner_spline.o co2_sw.o pckr_mod.o
timestepper.o : thc.o richards.o mphase.o \
                realization.o grid.o option.o field.o solver.o \
                output.o convergence.o checkpoint.o waypoint.o logging.o \
                mass_balance.o
hdf5.o : realization.o option.o logging.o
checkpoint.o : realization.o logging.o
mass_balance.o : realization.o option.o field.o solver.o logging.o thc.o
database.o : reaction.o fileio.o option.o utility.o reaction_aux.o
