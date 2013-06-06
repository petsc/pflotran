# $Id: makefile,v 1.3 2004/07/31 01:16:44 lichtner Exp $

pflotran_src = ./
common_src   = ./
pflotran_test_dir = ../../regression_tests


MYFLAGS = -I.

ifdef geomech
  MYFLAGS += -DGEOMECH
endif

ifdef satish
  MYFLAGS += -DGEOMECH_DEBUG
endif

ifdef uo2
  MYFLAGS += -DCOMPUTE_INTERNAL_MASS_FLUX
endif

ifdef glenn
  MYFLAGS += -DXINGYUAN_BC
endif

ifdef solid_solution
  MYFLAGS += -DSOLID_SOLUTION
endif

ifdef radioactive_decay
  MYFLAGS += -DRADIOACTIVE_DECAY
endif

ifdef mualem_spline
  MYFLAGS += -DMUALEM_SPLINE
endif

ifdef petsc_api_3_3
  MYFLAGS += -DHAVE_PETSC_API_3_3
endif

ifdef have_petscpcmg_h
	MYFLAGS += -DHAVE_PETSCPCMG_H
endif

ifdef dbl
  MYFLAGS += -DDOUBLE_LAYER
endif

ifdef pcl
  MYFLAGS += -DPCL -DMUALEM_SPLINE
endif

ifdef condnr
  MYFLAGS += -DCONDNR
endif

ifdef perm
  MYFLAGS += -DPERM
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

ifdef ugrid_debug
  MYFLAGS += -DUGRID_DEBUG
endif

ifdef surface_flow
  MYFLAGS += -DSURFACE_FLOW -w
endif

ifdef amanzi_bgd
  MYFLAGS += -DAMANZI_BGD
endif

ifdef store_flowrates
  MYFLAGS += -DSTORE_FLOWRATES
endif

ifdef mfd_ugrid
  MYFLAGS += -DMFD_UGRID
endif

CFLAGS   = 
FFLAGS   =
CPPFLAGS = ${MYFLAGS}
FPPFLAGS = ${MYFLAGS}

CLEANFILES       = pflotran pflotran_rxn libpflotranchem.a make.log

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

ifdef have_hdf5
LIBS +=  -L${HDF5_LIB} -lhdf5_fortran -lhdf5 -lz 
endif

# Set this accordingly on your platform
# SCORPIO_DIR=${HOME}/soft/scorpio
ifdef scorpio
  LIBS += -L${SCORPIO_DIR}/lib -lscorpio
  MYFLAGS += -DSCORPIO
  MYFLAGS += -DSCORPIO_WRITE 
  MYFLAGS += -I${SCORPIO_DIR}/include
endif

ifdef tau
FC = /ccs/proj/geo002/tau/tau-2.19.1/craycnl/bin/tau_f90.sh -optCompInst -tau_makefile=$(TAU_MAKEFILE)
FLINKER = /ccs/proj/geo002/tau/tau-2.19.1/craycnl/bin/tau_f90.sh -optCompInst -tau_makefile=$(TAU_MAKEFILE)
endif

util_obj  = \
	${common_src}input.o \
	${common_src}logging.o \
	${common_src}matrix_block_aux.o \
	${common_src}option.o \
	${common_src}output_aux.o \
	${common_src}string.o \
	${common_src}units.o \
	${common_src}utility.o \
	${common_src}variables.o

eos_obj = \
	${common_src}co2eos.o \
	${common_src}co2_sw.o \
	${common_src}co2_sw_rtsafe.o \
	${common_src}co2_span_wagner.o \
	${common_src}co2_span_wagner_spline.o \
	${common_src}gas_eos_mod.o \
	${common_src}water_eos.o 

mode_aux_obj = \
	${common_src}auxiliary.o \
	${common_src}flash2_aux.o \
	${common_src}general_aux.o \
	${common_src}immis_aux.o \
	${common_src}miscible_aux.o \
	${common_src}mphase_aux.o \
	${common_src}richards_aux.o \
	${common_src}th_aux.o \
	${common_src}thc_aux.o \
	${common_src}thmc_aux.o

mode_obj = \
	${common_src}flash2.o \
	${common_src}general.o \
	${common_src}global.o \
	${common_src}immis.o \
	${common_src}miscible.o \
	${common_src}mphase.o \
	${common_src}reactive_transport.o \
	${common_src}richards.o \
	${common_src}richards_common.o \
	${common_src}richards_lsm.o \
	${common_src}richards_mfd.o \
	${common_src}th.o \
	${common_src}thc.o \
	${common_src}thmc.o 

shared_mode_aux_obj = \
	${common_src}global_aux.o \
	${common_src}reactive_transport_aux.o 

chem_obj = \
	${common_src}constraint.o \
	${common_src}database.o \
	${common_src}database_aux.o \
	${common_src}database_hpt.o \
	${common_src}immobile.o \
	${common_src}immobile_aux.o \
	${common_src}microbial.o \
	${common_src}microbial_aux.o \
	${common_src}mineral.o \
	${common_src}mineral_aux.o \
	${common_src}reaction.o \
	${common_src}reaction_aux.o \
	${common_src}reaction_sandbox.o \
	${common_src}reaction_sandbox_base.o \
	${common_src}reaction_sandbox_clm_cn.o \
	${common_src}reaction_sandbox_example.o \
	${common_src}solid_solution.o \
	${common_src}solid_solution_aux.o \
	${common_src}surface_complexation.o \
	${common_src}surface_complexation_aux.o 

grid_obj = \
	${common_src}connection.o \
	${common_src}discretization.o \
	${common_src}general_grid.o \
	${common_src}grid.o \
	${common_src}mfd.o \
	${common_src}mfd_aux.o \
	${common_src}unstructured_cell.o \
	${common_src}unstructured_explicit.o \
	${common_src}unstructured_grid_aux.o \
	${common_src}unstructured_grid.o \
	${common_src}structured_grid.o

relations_obj = \
	${common_src}saturation.o \
	${common_src}saturation_function.o \
	${common_src}mphase_pckr_mod.o

properties_obj = \
	${common_src}fluid.o \
	${common_src}material_aux.o \
	${common_src}material.o 

dataset_obj = \
	${common_src}dataset.o \
	${common_src}dataset_aux.o \
	${common_src}dataset_base.o \
	${common_src}dataset_global.o \
	${common_src}time_series.o \
	${common_src}time_storage.o 

io_obj = \
	${common_src}hdf5.o \
	${common_src}hdf5_aux.o \
	${common_src}output.o \
	${common_src}output_common.o \
	${common_src}output_hdf5.o \
	${common_src}output_observation.o \
	${common_src}output_surface.o \
	${common_src}output_tecplot.o \
	${common_src}output_vtk.o

misc_obj = \
	${common_src}block_solve.o \
	${common_src}block_tridiag.o \
	${common_src}checkpoint.o \
	${common_src}condition.o \
	${common_src}condition_control.o \
	${common_src}convergence.o \
	${common_src}coupler.o \
	${common_src}debug.o \
	${common_src}field.o \
	${common_src}geometry.o \
	${common_src}hydrostatic.o \
	${common_src}init.o \
	${common_src}level.o \
	${common_src}mass_transfer.o \
	${common_src}matrix_buffer.o \
	${common_src}observation.o \
	${common_src}patch.o \
	${common_src}realization.o \
	${common_src}realization_base.o \
	${common_src}region.o \
	${common_src}regression.o \
	${common_src}secondary_continuum_aux.o\
	${common_src}secondary_continuum.o \
	${common_src}simulation.o \
	${common_src}solver.o \
	${common_src}spline.o \
	${common_src}stochastic.o \
	${common_src}stochastic_aux.o \
	${common_src}strata.o \
	${common_src}timestepper.o \
	${common_src}transport.o \
	${common_src}uniform_velocity.o \
	${common_src}waypoint.o 
#	${common_src}reaction_chunk.o\

surface_obj = \
	${common_src}surface_auxiliary.o \
	${common_src}surface_field.o \
	${common_src}surface_flow.o \
	${common_src}surface_flow_aux.o \
	${common_src}surface_init.o \
	${common_src}surface_global_aux.o \
	${common_src}surface_global.o \
	${common_src}surface_material.o \
	${common_src}surface_realization.o \
	${common_src}surface_th_aux.o \
	${common_src}surface_th.o
  
geomech_obj = \
  ${common_src}geomechanics_init.o \
  ${common_src}geomechanics_realization.o \
  ${common_src}geomechanics_discretization.o \
  ${common_src}geomech_grid.o \
  ${common_src}geomech_grid_aux.o \
  ${common_src}geomechanics_material.o \
  ${common_src}geomechanics_patch.o\
  ${common_src}geomechanics_field.o


pflotran_obj = $(util_obj) $(eos_obj) $(mode_aux_obj) $(mode_obj) \
	$(shared_mode_aux_obj) $(chem_obj) $(grid_obj) \
	$(relations_obj) $(properties_obj) $(dataset_obj) \
	$(io_obj) $(misc_obj) $(surface_obj) ${geomech_obj} \
  ${pflotran_src}pflotran.o

pflotran_rxn_obj = ${pflotran_src}pflotran_rxn.o

pflotran : $(pflotran_obj)
	${FLINKER} -o pflotran $(pflotran_obj) ${PETSC_LIB} ${LIBS} 

libpflotranchem.a : $(util_obj) $(shared_mode_aux_obj) $(chem_obj)
	$(AR) $(ARFLAGS) $@ $(util_obj) $(shared_mode_aux_obj) $(chem_obj)

pflotran_rxn : libpflotranchem.a $(pflotran_rxn_obj)
	${FLINKER} -o $@ $(pflotran_rxn_obj) -L. -lpflotranchem ${PETSC_LIB} ${LIBS} 

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
auxiliary.o : th_aux.o thc_aux.o thmc_aux.o richards_aux.o \
              reactive_transport_aux.o \
              mphase_aux.o flash2_aux.o global_aux.o immis_aux.o \
              miscible_aux.o general_aux.o material_aux.o \
              secondary_continuum_aux.o surface_global_aux.o
block_solve.o:
block_tridiag.o:
checkpoint.o : realization.o logging.o output_aux.o variables.o \
               realization_base.o
condition.o : input.o utility.o option.o field.o logging.o reaction_aux.o \
              string.o hdf5_aux.o dataset.o time_series.o \
              reactive_transport_aux.o constraint.o general_aux.o
condition_control.o : realization.o discretization.o region.o option.o grid.o \
                      field.o coupler.o condition.o connection.o reaction.o \
                      hdf5.o mfd.o dataset_aux.o surface_realization.o \
                      constraint.o general_aux.o
convergence.o : option.o solver.o grid.o
ifdef PFLOTRAN_RXN_FLAG
constraint.o : reaction_aux.o reactive_transport_aux.o global_aux.o \
               surface_complexation_aux.o mineral_aux.o option.o input.o \
               immobile.o \
               units.o string.o
else
constraint.o : reaction_aux.o reactive_transport_aux.o global_aux.o \
               surface_complexation_aux.o mineral_aux.o option.o input.o \
               immobile.o \
               units.o string.o logging.o
endif
coupler.o : condition.o connection.o region.o input.o grid.o mfd_aux.o
co2_span_wagner.o : co2_sw_rtsafe.o
co2_span_wagner_spline.o : spline.o
co2_sw.o : co2_span_wagner.o
database.o : reaction.o input.o option.o utility.o reaction_aux.o \
             database_aux.o surface_complexation_aux.o string.o mineral_aux.o \
             microbial.o microbial_aux.o immobile.o immobile_aux.o
database_aux.o : option.o string.o input.o
database_hpt.o : database.o reaction.o database_aux.o mineral_aux.o
dataset.o : hdf5_aux.o utility.o option.o dataset_aux.o
dataset_aux.o : string.o utility.o option.o input.o
dataset_base.o : time_storage.o option.o utility.o
dataset_global.o : dataset_base.o discretization.o grid.o option.o \
                   time_storage.o logging.o 
discretization.o : grid.o structured_grid.o unstructured_grid.o option.o \
                   input.o string.o mfd_aux.o mfd.o unstructured_grid_aux.o
flash2.o : flash2_aux.o realization.o fluid.o global_aux.o variables.o \
           output_aux.o
flash2_aux.o : option.o material.o water_eos.o co2_span_wagner.o \
               co2eos.o co2_span_wagner_spline.o co2_sw.o mphase_pckr_mod.o fluid.o \
               saturation_function.o
fluid.o : option.o input.o string.o
general.o : general_aux.o global_aux.o realization.o option.o level.o patch.o \
            field.o coupler.o grid.o connection.o material.o material_aux.o \
            discretization.o water_eos.o debug.o saturation_function.o
general_aux.o : option.o water_eos.o gas_eos_mod.o saturation_function.o \
                global_aux.o
general_grid.o : grid.o realization.o hdf5.o logging.o hdf5_aux.o connection.o \
                 option.o logging.o hdf5_aux.o discretization.o field.o \
                 patch.o
geometry.o : input.o option.o
geomechanics_init.o : geomechanics_discretization.o geomech_grid.o \
                      geomechanics_realization.o geomech_grid_aux.o \
                      input.o string.o geomechanics_patch.o level.o \
                      option.o geomechanics_material.o 
geomech_grid_aux.o : unstructured_cell.o utility.o
geomech_grid.o : geomech_grid_aux.o unstructured_cell.o utility.o \
                 input.o option.o string.o 
geomechanics_discretization: geomech_grid.o geomech_grid_aux.o
geomech_patch.o : option.o geomech_grid.o geomech_grid_aux.o \
                  geomechanics_material.o region.o strata.o
geomechanics_realization.o : geomech_grid.o geomechanics_discretization.o \
                             geomech_grid_aux.o input.o \
	                     geomechanics_material.o geomechanics_patch.o\
	                     waypoint.o geomechanics_field.o option.o
geomechanics_material.o : option.o string.o input.o
global.o : realization.o level.o patch.o global_aux.o connection.o grid.o \
           option.o coupler.o variables.o realization_base.o
global_aux.o : option.o
grid.o : structured_grid.o unstructured_grid.o option.o region.o connection.o \
         input.o mfd_aux.o unstructured_grid_aux.o unstructured_explicit.o
hdf5.o : realization.o discretization.o grid.o region.o patch.o field.o \
         realization_base.o
hdf5_aux.o : option.o logging.o dataset_aux.o units.o
hydrostatic.o : option.o grid.o coupler.o condition.o connection.o region.o \
                structured_grid.o dataset.o time_series.o water_eos.o \
                utility.o dataset_aux.o general_aux.o
immis.o : immis_aux.o realization.o global_aux.o variables.o output_aux.o
immis_aux.o : option.o material.o water_eos.o co2_span_wagner.o \
              co2eos.o co2_span_wagner_spline.o co2_sw.o mphase_pckr_mod.o fluid.o \
              saturation_function.o
immobile.o : immobile_aux.o option.o string.o input.o utility.o
immobile_aux.o : option.o utility.o
init.o : option.o grid.o debug.o realization.o solver.o timestepper.o \
         simulation.o input.o material.o region.o condition.o convergence.o \
         waypoint.o field.o general_grid.o co2_span_wagner.o utility.o \
         co2_span_wagner_spline.o co2_sw.o fluid.o saturation_function.o \
         coupler.o strata.o mphase_pckr_mod.o mphase.o flash2.o \
         richards.o th.o thc.o thmc.o discretization.o logging.o \
         database.o reaction_aux.o string.o units.o stochastic_aux.o \
         uniform_velocity.o dataset.o condition_control.o surface_flow.o \
         database_hpt.o output_aux.o variables.o constraint.o output_tecplot.o \
         output.o output_surface.o reaction_aux.o richards_mfd.o \
         surface_init.o surface_th.o surface_global.o \
		 geomechanics_init.o
input.o : string.o option.o
level.o : patch.o option.o
mass_transfer.o : dataset_global.o option.o input.o string.o discretization.o \
                  grid.o
material.o : option.o input.o string.o dataset_aux.o material_aux.o \
             surface_material.o
material_aux.o : option.o
matrix_block.o : option.o
matrix_buffer.o : grid.o structured_grid.o
mfd.o     : option.o mfd_aux.o grid.o
mfd_aux.o : connection.o
microbial.o : microbial_aux.o reaction_aux.o reactive_transport_aux.o \
              global_aux.o option.o string.o input.o utility.o
microbial_aux.o : database_aux.o utility.o
mineral.o : mineral_aux.o reaction_aux.o reactive_transport_aux.o global_aux.o \
            option.o string.o input.o utility.o solid_solution_aux.o \
            database_aux.o units.o
mineral_aux.o : database_aux.o option.o string.o utility.o
miscible.o : miscible_aux.o realization.o global_aux.o variables.o \
             output_aux.o
miscible_aux.o : option.o material.o
mphase.o : mphase_aux.o realization.o fluid.o global_aux.o \
           secondary_continuum.o secondary_continuum_aux.o variables.o \
           output_aux.o
mphase_aux.o : option.o material.o water_eos.o co2_span_wagner.o \
               co2eos.o co2_span_wagner_spline.o co2_sw.o mphase_pckr_mod.o \
               fluid.o saturation_function.o 
observation.o : region.o option.o input.o
output_aux.o : variables.o
option.o :
output.o : realization_base.o output_aux.o output_hdf5.o output_observation.o \
           output_tecplot.o output_vtk.o option.o grid.o discretization.o \
           field.o patch.o reaction_aux.o variables.o hdf5.o connection.o \
           coupler.o 
output_aux.o :
output_common.o : output_aux.o option.o grid.o field.o connection.o coupler.o \
                  patch.o realization_base.o variables.o logging.o \
                  unstructured_grid_aux.o unstructured_cell.o 
output_hdf5.o : realization_base.o output_aux.o output_common.o logging.o \
                option.o discretization.o grid.o field.o patch.o \
                hdf5.o hdf5_aux.o variables.o connection.o \
                unstructured_grid_aux.o
output_observation.o : realization_base.o logging.o output_aux.o option.o \
                       discretization.o grid.o field.o patch.o observation.o \
                       utility.o region.o reactive_transport.o reaction_aux.o \
                       structured_grid.o connection.o coupler.o variables.o \
                       global_aux.o mphase.o immis.o th.o thc.o thmc.o richards.o \
                       flash2.o miscible.o general.o reactive_transport.o 
output_surface.o : realization_base.o logging.o output_common.o output_hdf5.o \
                   output_tecplot.o option.o surface_realization.o \
                   realization.o discretization.o grid.o variables.o \
                   unstructured_grid_aux.o surface_field.o patch.o coupler.o \
                   connection.o utility.o hdf5.o hdf5_aux.o output_aux.o
output_tecplot.o : realization_base.o logging.o output_aux.o output_common.o \
                   grid.o option.o patch.o discretization.o structured_grid.o \
                   unstructured_grid_aux.o field.o variables.o \
                   unstructured_cell.o
output_vtk.o : realization_base.o output_aux.o output_common.o grid.o \
               discretization.o option.o structured_grid.o field.o patch.o \
               reaction_aux.o variables.o 
patch.o : option.o grid.o coupler.o observation.o strata.o region.o \
          auxiliary.o reaction.o condition.o connection.o hydrostatic.o \
          field.o saturation.o surface_field.o mineral.o output_aux.o \
          variables.o constraint.o general_aux.o reaction_aux.o richards_aux.o \
          reactive_transport_aux.o dataset_aux.o thc_aux.o thmc_aux.o \
          global_aux.o mphase_aux.o immis_aux.o miscible_aux.o mineral_aux.o \
          surface_complexation_aux.o surface_auxiliary.o \
		  geomech_grid_aux.o
pflotran.o : simulation.o output.o init.o timestepper.o logging.o \
             stochastic_aux.o stochastic.o surface_flow_aux.o regression.o
pflotran_rxn.o : reaction.o reaction_aux.o option.o input.o
ifdef PFLOTRAN_RXN_FLAG
reaction.o : reactive_transport_aux.o reaction_aux.o global_aux.o \
             surface_complexation_aux.o mineral_aux.o solid_solution_aux.o \
             surface_complexation.o mineral.o solid_solution.o option.o \
             string.o input.o microbial.o microbial_aux.o variables.o \
             reaction_sandbox.o constraint.o
else
reaction.o : reactive_transport_aux.o reaction_aux.o global_aux.o \
             surface_complexation_aux.o mineral_aux.o solid_solution_aux.o \
             surface_complexation.o mineral.o solid_solution.o option.o \
             string.o input.o microbial.o microbial_aux.o variables.o \
             reaction_sandbox.o constraint.o immobile.o immobile_aux.o \
             co2eos.o water_eos.o 
endif
reaction_aux.o : option.o string.o database_aux.o mineral_aux.o \
                 surface_complexation_aux.o solid_solution_aux.o utility.o \
                 microbial_aux.o immobile_aux.o
reaction_sandbox.o : option.o reaction_aux.o reactive_transport_aux.o string.o \
                     input.o utility.o global_aux.o reaction_sandbox_base.o \
                     reaction_sandbox_clm_cn.o reaction_sandbox_example.o
reaction_sandbox_base.o : reaction_aux.o reactive_transport_aux.o \
                          global_aux.o
reaction_sandbox_clm_cn.o : option.o reaction_aux.o reactive_transport_aux.o \
                            string.o input.o utility.o global_aux.o \
                            reaction_sandbox_base.o database_aux.o \
                            immobile_aux.o
reaction_sandbox_example.o : reaction_aux.o option.o reactive_transport_aux.o \
                             global_aux.o units.o string.o input.o
#reaction_chunk.o : reactive_transport_aux.o reaction_aux.o option.o string.o
reactive_transport.o : realization.o reaction.o transport.o matrix_block_aux.o \
                       reactive_transport_aux.o fluid.o logging.o global.o \
                       output_aux.o variables.o secondary_continuum_aux.o\
                       secondary_continuum.o mass_transfer.o 
reactive_transport_aux.o : option.o reaction_aux.o surface_complexation_aux.o \
                           matrix_block_aux.o
realization.o : grid.o option.o region.o condition.o coupler.o material.o \
                strata.o waypoint.o field.o hydrostatic.o debug.o \
                observation.o uniform_velocity.o patch.o level.o \
                discretization.o \
                reaction.o reaction_aux.o fluid.o saturation_function.o \
                dataset.o surface_field.o surface_material.o output_aux.o \
                constraint.o realization_base.o mass_transfer.o
realization_base.o : level.o patch.o discretization.o option.o input.o \
                     debug.o output_aux.o field.o reaction_aux.o \
                     mass_transfer.o 
region.o : input.o utility.o logging.o string.o geometry.o option.o \
           unstructured_cell.o structured_grid.o
regression.o : realization.o discretization.o timestepper.o option.o string.o \
               utility.o grid.o output.o output_aux.o output_common.o
richards.o : richards_aux.o realization.o saturation_function.o \
             matrix_buffer.o logging.o output_aux.o variables.o richards_lsm.o \
             richards_mfd.o richards_common.o
richards_aux.o : option.o material.o water_eos.o saturation_function.o \
                matrix_buffer.o global_aux.o
richards_common.o : richards_aux.o global_aux.o saturation_function.o option.o \
                    water_eos.o
richards_lsm.o : richards_aux.o global_aux.o realization.o patch.o option.o \
                 field.o grid.o coupler.o connection.o material.o logging.o \
                 saturation_function.o richards_common.o
richards_mfd.o : richards_aux.o realization.o saturation_function.o \
                 matrix_buffer.o logging.o mfd.o mfd_aux.o richards_common.o
saturation.o : option.o grid.o coupler.o condition.o connection.o region.o
saturation_function.o : option.o input.o string.o utility.o
secondary_continuum.o : string.o secondary_continuum_aux.o reaction_aux.o\
                        material.o option.o reactive_transport_aux.o\
                        block_solve.o global_aux.o mphase_aux.o reaction.o\
						block_tridiag.o
secondary_continuum_aux.o : option.o reactive_transport_aux.o utility.o
simulation.o : timestepper.o realization.o regression.o \
               geomechanics_realization.o
solid_solution.o : solid_solution_aux.o option.o string.o reaction_aux.o \
                   mineral_aux.o utility.o mineral.o input.o
solid_solution_aux.o : mineral_aux.o
solver.o : option.o field.o input.o string.o
stochastic.o : simulation.o option.o realization.o timestepper.o \
               init.o logging.o regression.o
stochastic_aux.o : simulation.o input.o option.o
strata.o : input.o region.o material.o
string.o :
structured_grid.o : connection.o option.o input.o
surface_auxiliary.o : surface_global_aux.o surface_flow_aux.o \
                      surface_th_aux.o
surface_complexation.o : surface_complexation_aux.o reaction_aux.o \
                         global_aux.o  reactive_transport_aux.o option.o \
                         string.o input.o utility.o
surface_complexation_aux.o : option.o string.o database_aux.o 
surface_flow.o : surface_flow_aux.o surface_realization.o surface_field.o \
                 surface_material.o solver.o
surface_init.o : option.o input.o string.o surface_material.o \
                 surface_realization.o grid.o structured_grid.o \
                 unstructured_grid.o unstructured_grid_aux.o region.o \
                 condition.o dataset_aux.o coupler.o strata.o units.o \
                 waypoint.o debug.o patch.o output_aux.o
surface_global_aux.o : option.o
surface_global.o : surface_global_aux.o surface_realization.o level.o \
                   patch.o option.o coupler.o connection.o grid.o \
                   variables.o patch.o field.o 
surface_material.o : string.o option.o input.o
surface_realization.o : condition.o debug.o discretization.o input.o level.o \
                        option.o patch.o region.o surface_field.o \
                        surface_material.o waypoint.o dataset_aux.o coupler.o \
                        string.o grid.o strata.o unstructured_grid.o dataset.o \
                        unstructured_grid_aux.o realization_base.o
surface_th_aux.o: option.o material.o water_eos.o \
            saturation_function.o grid.o surface_global_aux.o
surface_th.o : surface_global_aux.o surface_realization.o output_aux.o variables.o \
               surface_material.o surface_th_aux.o
th.o : secondary_continuum.o th_aux.o realization.o saturation_function.o \
       variables.o secondary_continuum_aux.o 
th_aux.o : option.o material.o water_eos.o \
            saturation_function.o grid.o global_aux.o
thc.o : secondary_continuum.o thc_aux.o realization.o saturation_function.o \
        variables.o secondary_continuum_aux.o
thc_aux.o : option.o material.o water_eos.o \
            saturation_function.o grid.o global_aux.o
thmc.o : thmc_aux.o realization.o saturation_function.o
thmc_aux.o : option.o material.o water_eos.o saturation_function.o grid.o \
             global_aux.o
time_storage.o : option.o
timestepper.o : th.o thc.o thmc.o richards.o mphase.o flash2.o general.o \
                realization.o grid.o option.o field.o solver.o \
                output.o convergence.o checkpoint.o waypoint.o logging.o \
                input.o string.o reactive_transport.o condition_control.o \
                surface_flow.o output_aux.o output_surface.o reaction_aux.o \
                surface_th.o
time_series.o : string.o option.o input.o
transport.o : reactive_transport_aux.o option.o matrix_block_aux.o
units.o : option.o
unstructured_cell.o : option.o utility.o 
unstructured_explicit.o : geometry.o unstructured_grid_aux.o
unstructured_grid.o : option.o connection.o unstructured_cell.o \
                      unstructured_grid_aux.o geometry.o reaction_aux.o hdf5_aux.o \
                      mfd_aux.o
unstructured_grid_aux.o : unstructured_cell.o geometry.o option.o utility.o
utility.o : string.o input.o option.o
variables.o :
uniform_velocity.o : option.o input.o string.o logging.o units.o
waypoint.o : option.o output_aux.o


# developer level regression testing
test : FORCE
	@if [ -d $(pflotran_test_dir) ]; then \
		$(MAKE) --directory=$(pflotran_test_dir) $(MAKECMDGOALS); \
	fi

clean-tests : FORCE
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

