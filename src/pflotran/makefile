# $Id: makefile,v 1.3 2004/07/31 01:16:44 lichtner Exp $

objdir = ./${PETSC_ARCH}/obj
srcdir = ./
pflotran_src = ./
common_src   = ./
bindir = ./${PETSC_ARCH}/bin
pflotran_test_dir = ../../test

MYFLAGS = -I${objdir} 

ifdef hdf5-vamsi
    MYFLAGS += -DVAMSI_HDF5 -DVAMSI_HDF5_READ -DVAMSI_HDF5_WRITE -DVAMSI_STAGE_BARRIER
endif   

ifdef compile_broken
  MYFLAGS += -DCOMPILE_BROKEN
endif

ifdef uo2
  MYFLAGS += -DCOMPUTE_INTERNAL_MASS_FLUX
endif

ifdef glenn
  #MYFLAGS += -DGLENN -DUGRID_DEBUG
  MYFLAGS += -DUGRID_DEBUG -UGRID_NEW
  #MYFLAGS += -DGLENN -DTVD_DEBUG
endif

ifdef solid_solution
  MYFLAGS += -DSOLID_SOLUTION
endif

ifdef radioactive_decay
  MYFLAGS += -DRADIOACTIVE_DECAY
endif

ifdef matcreate_old
  MYFLAGS += -DMATCREATE_OLD
endif

ifdef mualem_spline
  MYFLAGS += -DMUALEM_SPLINE
endif

ifdef snes_old
  MYFLAGS += -DHAVE_SNES_API_3_2
endif

ifdef snes_api_3_2
  MYFLAGS += -DHAVE_SNES_API_3_2
endif

ifdef petsc_api_3_3
  MYFLAGS += -DHAVE_PETSC_API_3_3
endif

ifdef dbl
  MYFLAGS += -DDOUBLE_LAYER
endif

ifdef pcl
  MYFLAGS += -DPCL -DMUALEM_SPLINE
endif

ifdef tab
  MYFLAGS += -DTABLE
endif

ifdef temp
  MYFLAGS += -DTEMP_DEPENDENT_LOGK
endif

ifdef ice
  MYFLAGS += -DICE
endif

ifdef remove_sat
  MYFLAGS += -DREMOVE_SATURATION
endif

ifdef compressibility
  MYFLAGS += -DUSE_COMPRESSIBILITY
endif

ifdef no_vapor_diffusion
  MYFLAGS += -DNO_VAPOR_DIFFUSION
endif

ifdef scco2
  MYFLAGS += -DCHUAN_CO2
endif

ifdef debug_gen
  MYFLAGS += -DDEBUG_GENERAL
endif 

ifdef coll
  MYFLAGS += -DCOLL
endif

ifdef dont_use_wateos
  MYFLAGS += -DDONT_USE_WATEOS
endif

ifdef vamsi_hdf5
  MYFLAGS += -DVAMSI_HDF5 -DVAMSI_DEFAULT -DVAMSI_HDF5_READ -DVAMSI_HDF5_WRITE -DVAMSI_STAGE_BARRIER
endif

ifdef have_hdf5
  MYFLAGS += -I$(HDF5_INCLUDE) -I$(HDF5_LIB) -DPETSC_HAVE_HDF5
endif

ifdef gpu
  MYFLAGS += -DCHUNK
endif

ifdef mfd
  MYFLAGS += -DDASVYAT
endif

ifdef dasvyat_debug
  MYFLAGS += -DDASVYAT_DEBUG
endif

ifdef dasvyat_test_cut
  MYFLAGS += -DDASVYAT_TEST_CUT
endif

ifdef gautam
  MYFLAGS += -w -g #-DSURFACE_FLOW #-DUGRID_DEBUG #-DPARALLELIO_LIB -DOS_STATISTICS 
endif

ifdef surface_flow
  MYFLAGS += -DSURFACE_FLOW
endif

ifdef amanzi_bgd
  MYFLAGS += -DAMANZI_BGD
endif

CFLAGS		 = 
FFLAGS		 =
CPPFLAGS         = ${MYFLAGS}
FPPFLAGS         = ${MYFLAGS}

CLEANFILES       = pflotran 

ifdef petscrel
 include ${PETSC_DIR}/bmake/common/base
else
 include ${PETSC_DIR}/conf/variables
 include ${PETSC_DIR}/conf/rules
endif

ifdef have_hdf5
LIBS +=  -L${HDF5_LIB} -lhdf5_fortran -lhdf5 -lz 
endif

ifdef have_parallelio_lib
LIBS += -L${PARALLELIO_LIB} -lparallelio
endif

ifdef tau
FC = /ccs/proj/geo002/tau/tau-2.19.1/craycnl/bin/tau_f90.sh -optCompInst -tau_makefile=$(TAU_MAKEFILE)
FLINKER = /ccs/proj/geo002/tau/tau-2.19.1/craycnl/bin/tau_f90.sh -optCompInst -tau_makefile=$(TAU_MAKEFILE)
endif

util_obj  = ${common_src}input.o \
	${common_src}string.o \
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
	${common_src}geometry.o\
	${common_src}units.o\
	${common_src}logging.o \
	${common_src}output_aux.o \
	${common_src}hdf5_aux.o \
	${common_src}global_aux.o\
	${common_src}material_aux.o\
	${common_src}matrix_block_aux.o\
	${common_src}database_aux.o\
	${common_src}secondary_continuum.o\
	${common_src}surface_complexation_aux.o\
	${common_src}mineral_aux.o\
	${common_src}microbial_aux.o\
	${common_src}solid_solution_aux.o\
	${common_src}reaction_aux.o\
	${common_src}reactive_transport_aux.o\
	${common_src}thc_aux.o\
	${common_src}thmc_aux.o\
	${common_src}richards_aux.o\
	${common_src}general_aux.o\
	${common_src}mphase_aux.o\
	${common_src}miscible_aux.o\
	${common_src}immis_aux.o\
	${common_src}mfd_aux.o\
	${common_src}flash2_aux.o\
	${common_src}auxilliary.o\
	${common_src}field.o \
	${common_src}debug.o \
	${common_src}waypoint.o \
	${common_src}region.o \
	${common_src}fluid.o \
	${common_src}time_series.o \
	${common_src}dataset_aux.o \
	${common_src}dataset.o \
	${common_src}saturation_function.o \
	${common_src}material.o \
	${common_src}pckr_mod.o \
	${common_src}strata.o \
	${common_src}connection.o \
	${common_src}observation.o \
	${common_src}solver.o \
	${common_src}cf90bridge.o\
	${common_src}unstructured_cell.o \
	${common_src}unstructured_explicit.o \
	${common_src}unstructured_grid_aux.o \
	${common_src}unstructured_grid.o \
	${common_src}structured_grid.o \
	${common_src}condition.o \
	${common_src}coupler.o \
	${common_src}grid.o \
	${common_src}matrix_buffer.o \
	${common_src}hydrostatic.o \
	${common_src}saturation.o \
	${common_src}patch.o \
	${common_src}mfd.o\
	${common_src}level.o \
	${common_src}discretization.o \
	${common_src}database.o \
	${common_src}database_hpt.o \
	${common_src}velocity.o \
	${common_src}realization.o \
	${common_src}hdf5.o \
	${common_src}global.o \
	${common_src}mphase.o\
	${common_src}immis.o\
	${common_src}miscible.o\
	${common_src}flash2.o\
	${common_src}transport.o\
	${common_src}surface_complexation.o\
	${common_src}mineral.o\
	${common_src}microbial.o\
	${common_src}solid_solution.o\
	${common_src}reaction.o\
	${common_src}reactive_transport.o\
	${common_src}thc.o\
	${common_src}thmc.o\
	${common_src}richards.o\
	${common_src}general.o\
	${common_src}condition_control.o \
	${common_src}general_grid.o \
	${common_src}checkpoint.o \
	${common_src}convergence.o \
	${common_src}timestepper.o \
	${common_src}simulation.o \
	${common_src}stochastic_aux.o \
	${common_src}stochastic.o \
	${common_src}surface_field.o \
	${common_src}surface_flow_aux.o \
	${common_src}surface_flow.o \
	${common_src}surface_material.o \
	${common_src}surface_realization.o \
	${common_src}output.o \
	${common_src}regression.o \
	${common_src}init.o
#	${common_src}reaction_chunk.o\

pflotran : $(util_obj) ${pflotran_obj} ${pflotran_src}pflotran.o
	${FLINKER} -o pflotran $(util_obj) ${pflotran_obj} ${pflotran_src}pflotran.o ${PETSC_LIB} ${LIBS} 

ascii_to_hdf5: input.o hdf5_dataset.o ascii_to_hdf5.o 
	${FLINKER} -o ascii_to_hdf5 ascii_to_hdf5.o hdf5_dataset.o input.o \
	$(LIBS)

# do these files exist anymore?
#test: functions.o test.o
#	${FLINKER} -o test test.o functions.o $(LIBS) $(PETSC_LIB)

# Should add this to default PETSc targets as well?
%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<

# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
auxilliary.o : thc_aux.o thmc_aux.o richards_aux.o reactive_transport_aux.o \
               mphase_aux.o flash2_aux.o global_aux.o immis_aux.o \
               miscible_aux.o
checkpoint.o : realization.o logging.o output_aux.o
condition.o : input.o utility.o option.o field.o logging.o reaction_aux.o \
              string.o hdf5_aux.F90 dataset.o time_series.o \
              reactive_transport_aux.o
condition_control.o : realization.o discretization.o region.o option.o grid.o \
                      field.o coupler.o condition.o connection.o reaction.o \
                      hdf5.o mfd.o dataset_aux.o surface_realization.o
convergence.o : option.o solver.o grid.o
coupler.o : condition.o connection.o region.o input.o grid.o mfd_aux.o \
            auxilliary.o
database.o : reaction.o input.o option.o utility.o reaction_aux.o \
             database_aux.o surface_complexation_aux.o string.o mineral_aux.o \
             microbial_aux.o
database_aux.o : option.o string.o
database_hpt.o : database.o reaction.o database_aux.o mineral_aux.o
dataset.o : hdf5_aux.o utility.o option.o dataset_aux.o
dataset_aux.o : string.o utility.o option.o input.o
discretization.o : grid.o structured_grid.o unstructured_grid.o option.o \
                   input.o string.o mfd_aux.o mfd.o unstructured_grid_aux.o
flash2.o : flash2_aux.o realization.o fluid.o global_aux.o
flash2_aux.o : option.o material.o water_eos.o co2_span_wagner.o \
               co2eos.o co2_span_wagner_spline.o co2_sw.o pckr_mod.o fluid.o \
               saturation_function.o
fluid.o : option.o input.o string.o
general.o : general_aux.o global_aux.o realization.o option.o level.o patch.o \
            field.o coupler.o grid.o connection.o material.o material_aux.o \
            discretization.o water_eos.o debug.o saturation_function.o
general_aux.o : option.o water_eos.o gaseos_mod.o saturation_function.o \
                global_aux.o
general_grid.o : grid.o realization.o hdf5.o logging.o
geometry.o : input.o option.o
global.o : realization.o level.o patch.o global_aux.o connection.o grid.o \
           option.o coupler.o
global_aux.o : option.o
grid.o : structured_grid.o unstructured_grid.o option.o region.o connection.o \
         input.o mfd_aux.o unstructured_grid_aux.o unstructured_explicit.o
hdf5.o : realization.o discretization.o grid.o region.o patch.o field.o
hdf5_aux.o : option.o logging.o dataset_aux.o
hydrostatic.o : option.o grid.o coupler.o condition.o connection.o region.o \
                structured_grid.o dataset.o time_series.o
immis.o : immis_aux.o realization.o  global_aux.o
immis_aux.o : option.o material.o water_eos.o co2_span_wagner.o \
              co2eos.o co2_span_wagner_spline.o co2_sw.o pckr_mod.o fluid.o \
              saturation_function.o
init.o : option.o grid.o debug.o realization.o solver.o timestepper.o \
         simulation.o input.o material.o region.o condition.o convergence.o \
         waypoint.o field.o general_grid.o co2_span_wagner.o utility.o \
         co2_span_wagner_spline.o co2_sw.o fluid.o saturation_function.o \
         coupler.o strata.o pckr_mod.o field.o mphase.o flash2.o \
         richards.o thc.o thmc.o discretization.o logging.o \
         database.o reaction_aux.o string.o units.o stochastic_aux.o \
         velocity.o dataset.o condition_control.o surface_flow.o \
         database_hpt.o output_aux.o
input.o : string.o option.o
level.o : patch.o option.o
material.o : option.o input.o string.o dataset_aux.o material_aux.o \
             surface_material.o
material_aux.o : option.o
matrix_buffer.o : grid.o structured_grid.o
mfd.o     : option.o mfd_aux.o grid.o
mfd_aux.o : connection.o
microbial.o : microbial_aux.o reaction_aux.o reactive_transport_aux.o \
              global_aux.o option.o string.o input.o utility.o
microbial_aux.o : database_aux.o utility.o
mineral.o : mineral_aux.o reaction_aux.o reactive_transport_aux.o global_aux.o \
            option.o string.o input.o utility.o solid_solution_aux.o database_aux.o
mineral_aux.o : database_aux.o option.o string.o
miscible.o : miscible_aux.o realization.o  global_aux.o
miscible_aux.o : option.o material.o
mphase.o : mphase_aux.o realization.o fluid.o global_aux.o \
           secondary_continuum.o
mphase_aux.o : option.o material.o water_eos.o co2_span_wagner.o \
               co2eos.o co2_span_wagner_spline.o co2_sw.o pckr_mod.o \
               fluid.o saturation_function.o secondary_continuum.o
observation.o : region.o option.o input.o
output.o : grid.o structured_grid.o option.o realization.o field.o \
           coupler.o connection.o hdf5.o observation.o \
           reactive_transport.o immis.o miscible.o \
           thc.o thmc.o richards.o mphase.o flash2.o general.o \
           discretization.o logging.o unstructured_grid_aux.o \
           surface_realization.o surface_flow.o output_aux.o
patch.o : option.o grid.o coupler.o observation.o strata.o region.o \
          auxilliary.o reaction.o condition.o connection.o hydrostatic.o \
          field.o saturation.o surface_field.o mineral.o output_aux.o
pflotran.o : simulation.o output.o init.o timestepper.o logging.o \
             stochastic_aux.o stochastic.o surface_flow_aux.o regression.o
reaction.o : reactive_transport_aux.o reaction_aux.o global_aux.o \
             surface_complexation_aux.o mineral_aux.o solid_solution_aux.o \
             surface_complexation.o mineral.o solid_solution.o option.o \
             string.o input.o condition.o co2eos.o water_eos.o microbial.o \
             microbial_aux.o
reaction_aux.o : option.o string.o database_aux.o mineral_aux.o \
                 surface_complexation_aux.o solid_solution_aux.o utility.o \
                 microbial_aux.o
#reaction_chunk.o : reactive_transport_aux.o reaction_aux.o option.o string.o
reactive_transport.o : realization.o reaction.o transport.o matrix_block_aux.o \
                       reactive_transport_aux.o fluid.o logging.o global.o \
                       output_aux.o
reactive_transport_aux.o : option.o reaction_aux.o surface_complexation_aux.o \
                           matrix_block_aux.o
realization.o : grid.o option.o region.o condition.o coupler.o material.o \
                strata.o waypoint.o field.o hydrostatic.o debug.o \
                observation.o velocity.o patch.o level.o discretization.o \
                reaction.o reaction_aux.o fluid.o saturation_function.o \
                dataset.o surface_field.o surface_material.o output_aux.o
region.o : input.o utility.o logging.o string.o geometry.o
regression.o : realization.o discretization.o timestepper.o option.o string.o \
               utility.o grid.o output.o output_aux.o
richards.o : richards_aux.o realization.o saturation_function.o \
             matrix_buffer.o logging.o output_aux.o
richards_aux.o : option.o material.o water_eos.o saturation_function.o \
                matrix_buffer.o global_aux.o
saturation.o : option.o grid.o coupler.o condition.o connection.o region.o
saturation_function.o : option.o input.o string.o utility.o
simulation.o : timestepper.o realization.o regression.o
solid_solution.o : solid_solution_aux.o option.o string.o reaction_aux.o mineral_aux.o \
                   utility.o mineral.o input.o
solid_solution_aux.o : mineral_aux.o
solver.o : option.o field.o input.o string.o
stochastic.o : simulation.o option.o realization.o timestepper.o \
               init.o logging.o regression.o
stochastic_aux.o : simulation.o input.o option.o
strata.o : input.o region.o material.o
structured_grid.o : connection.o option.o input.o cf90bridge.o
surface_complexation.o : surface_complexation_aux.o reaction_aux.o \
                         global_aux.o  reactive_transport_aux.o option.o \
                         string.o input.o utility.o
surface_complexation_aux.o : option.o string.o database_aux.o 
surface_flow.o : surface_flow_aux.o surface_realization.o surface_field.o \
                 surface_material.o
surface_material.o : string.o option.o input.o
surface_realization.o : condition.o debug.o discretization.o input.o level.o \
                        option.o patch.o region.o surface_field.o \
                        surface_material.o waypoint.o dataset_aux.o coupler.o \
                        string.o grid.o strata.o unstructured_grid.o dataset.o \
                        unstructured_grid_aux.o
thc.o : secondary_continuum.o thc_aux.o realization.o saturation_function.o
thc_aux.o : secondary_continuum.o option.o material.o water_eos.o \
            saturation_function.o grid.o
thmc.o : thmc_aux.o realization.o saturation_function.o
thmc_aux.o : option.o material.o water_eos.o saturation_function.o grid.o
timestepper.o : thc.o thmc.o richards.o mphase.o flash2.o general.o \
                realization.o grid.o option.o field.o solver.o \
                output.o convergence.o checkpoint.o waypoint.o logging.o \
                input.o string.o reactive_transport.o condition_control.o \
                surface_flow.o output_aux.o
time_series.o : string.o option.o input.o
transport.o : reactive_transport_aux.o option.o matrix_block_aux.o
units.o : option.o
unstructured_cell.o : option.o utility.o 
unstructured_explicit.o : geometry.o unstructured_grid_aux.o
unstructured_grid.o : option.o connection.o unstructured_cell.o \
                      unstructured_grid_aux.o geometry.o
unstructured_grid_aux.o : unstructured_cell.o
utility.o : string.o input.o option.o
velocity.o : option.o input.o string.o logging.o units.o
waypoint.o : option.o output_aux.o


# developer level regression testing
test : FORCE
	@if [ -d $(pflotran_test_dir) ]; then \
		$(MAKE) --directory=$(pflotran_test_dir) $(MAKECMDGOALS); \
	fi

# user lever tests to verify pflotran is built correctly
check : FORCE
	@if [ -d $(pflotran_test_dir) ]; then \
		$(MAKE) --directory=$(pflotran_test_dir) $(MAKECMDGOALS); \
	fi

# null rule to force things to happen
FORCE :

