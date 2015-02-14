# PFLOTRAN makefile

pflotran_src = ./
common_src   = ./
regression_test_dir = ../../regression_tests
unit_test_dir = ./unittests


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

MYFLAGS = -I.

# if compiling pflotran_rxn or libpflotranchem.a, need to toggle 
# PFLOTRAN_RXN_FLAG 
ifneq (,$(filter libpflotranchem.a pflotran_rxn,$(MAKECMDGOALS)))
  PFLOTRAN_RXN_FLAG := 1
endif

ifdef uo2
  MYFLAGS += ${FC_DEFINE_FLAG}COMPUTE_INTERNAL_MASS_FLUX
endif

ifdef glenn
#  MYFLAGS += ${FC_DEFINE_FLAG}XINGYUAN_BC
#  MYFLAGS += ${FC_DEFINE_FLAG}BROADCAST_DATASET
  MYFLAGS += ${FC_DEFINE_FLAG}REFACTOR_CHARACTERISTIC_CURVES
endif

ifdef ye_flux
  MYFLAGS += ${FC_DEFINE_FLAG}YE_FLUX
endif

ifdef solid_solution
  MYFLAGS += ${FC_DEFINE_FLAG}SOLID_SOLUTION
endif

ifdef dbl
  MYFLAGS += ${FC_DEFINE_FLAG}DOUBLE_LAYER
endif

ifdef pcl
  MYFLAGS += ${FC_DEFINE_FLAG}PCL ${FC_DEFINE_FLAG}MUALEM_SPLINE
endif

ifdef condnr
  MYFLAGS += ${FC_DEFINE_FLAG}CONDNR
endif

ifdef temp
  MYFLAGS += ${FC_DEFINE_FLAG}TEMP_DEPENDENT_LOGK
endif

ifdef remove_sat
  MYFLAGS += ${FC_DEFINE_FLAG}REMOVE_SATURATION
endif

ifdef no_vapor_diffusion
  MYFLAGS += ${FC_DEFINE_FLAG}NO_VAPOR_DIFFUSION
endif

ifdef debug_gen
  MYFLAGS += ${FC_DEFINE_FLAG}DEBUG_GENERAL
  ifeq ($(debug_gen),2)
    MYFLAGS += ${FC_DEFINE_FLAG}DEBUG_GENERAL_INFO
  endif
endif 

ifdef coll
  MYFLAGS += ${FC_DEFINE_FLAG}COLL
endif

ifdef have_hdf5
  MYFLAGS += -I$(HDF5_INCLUDE) -I$(HDF5_LIB) ${FC_DEFINE_FLAG}PETSC_HAVE_HDF5
endif

ifdef ugrid_debug
  MYFLAGS += ${FC_DEFINE_FLAG}UGRID_DEBUG
endif

ifdef amanzi_bgd
  MYFLAGS += ${FC_DEFINE_FLAG}AMANZI_BGD
endif

ifdef surface_th_debug
  MYFLAGS += ${FC_DEFINE_FLAG}SURFACE_TH_DEBUG
endif

ifdef th_bc_hack
  MYFLAGS += ${FC_DEFINE_FLAG}THDIRICHLET_TEMP_BC_HACK
endif

UPDATE_PROVENANCE=0
ifdef provenance
  UPDATE_PROVENANCE=1
endif

CFLAGS   = 
FFLAGS   =
CPPFLAGS = ${MYFLAGS}
FPPFLAGS = ${MYFLAGS}

CLEANFILES       = pflotran pflotran_rxn libpflotranchem.a make.log pflotran_provenance.F90

ifdef have_hdf5
LIBS +=  -L${HDF5_LIB} -lhdf5_fortran -lhdf5 -lz 
endif

# Set this accordingly on your platform
# SCORPIO_DIR=${HOME}/soft/scorpio
ifdef scorpio
  LIBS += -L${SCORPIO_DIR}/lib -lscorpio
  MYFLAGS += ${FC_DEFINE_FLAG}SCORPIO
  MYFLAGS += ${FC_DEFINE_FLAG}SCORPIO_WRITE
  MYFLAGS += -I${SCORPIO_DIR}/include
endif

ifdef tau
FC = /ccs/proj/geo002/tau/tau-2.19.1/craycnl/bin/tau_f90.sh -optCompInst -tau_makefile=$(TAU_MAKEFILE)
FLINKER = /ccs/proj/geo002/tau/tau-2.19.1/craycnl/bin/tau_f90.sh -optCompInst -tau_makefile=$(TAU_MAKEFILE)
endif

# Begin Source Block

util_obj  = \
	${common_src}input_aux.o \
	${common_src}logging.o \
	${common_src}matrix_block_aux.o \
	${common_src}option.o \
	${common_src}option_flow.o \
	${common_src}option_transport.o \
	${common_src}output_aux.o \
	${common_src}pflotran_constants.o \
	${common_src}pflotran_provenance.o \
	${common_src}string.o \
	${common_src}units.o \
	${common_src}utility.o \
	${common_src}variables.o

eos_obj = \
	${common_src}co2_sw.o \
	${common_src}co2_span_wagner_spline.o \
	${common_src}eos.o \
	${common_src}eos_gas.o \
	${common_src}gas_eos_mod.o 

mode_aux_obj = \
	${common_src}auxiliary.o \
	${common_src}flash2_aux.o \
	${common_src}general_aux.o \
	${common_src}immis_aux.o \
	${common_src}miscible_aux.o \
	${common_src}mphase_aux.o \
	${common_src}pm_base.o \
	${common_src}pm_base_pointer.o \
	${common_src}pm_general.o \
	${common_src}pm_flash2.o \
	${common_src}pm_immis.o \
	${common_src}pm_mphase.o \
	${common_src}pm_miscible.o \
	${common_src}pm_richards.o \
	${common_src}pm_rt.o \
	${common_src}pm_subsurface.o \
	${common_src}pm_surface_flow.o \
	${common_src}pm_surface_th.o \
	${common_src}pm_surface.o \
	${common_src}pm_th.o \
	${common_src}pm_waste_form.o \
	${common_src}richards_aux.o \
	${common_src}th_aux.o

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
	${common_src}th.o

shared_mode_aux_obj = \
	${common_src}co2eos.o \
	${common_src}co2_span_wagner.o \
	${common_src}co2_sw_rtsafe.o \
	${common_src}eos_water.o \
	${common_src}global_aux.o \
	${common_src}material_aux.o \
	${common_src}reactive_transport_aux.o 

chem_obj = \
	${common_src}transport_constraint.o \
	${common_src}reaction_database.o \
	${common_src}reaction_database_aux.o \
	${common_src}reaction_database_hpt.o \
	${common_src}reaction_immobile.o \
	${common_src}reaction_immobile_aux.o \
	${common_src}reaction_microbial.o \
	${common_src}reaction_microbial_aux.o \
	${common_src}reaction_mineral.o \
	${common_src}reaction_mineral_aux.o \
	${common_src}reaction.o \
	${common_src}reaction_aux.o \
	${common_src}reaction_clm.o \
	${common_src}reaction_sandbox.o \
	${common_src}reaction_sandbox_base.o \
	${common_src}reaction_sandbox_clm_cn.o \
	${common_src}reaction_sandbox_example.o \
	${common_src}reaction_sandbox_ufd_wp.o \
	${common_src}reaction_solid_solution.o \
	${common_src}reaction_solid_solution_aux.o \
	${common_src}reaction_surf_complex.o \
	${common_src}reaction_surf_complex_aux.o 

grid_obj = \
	${common_src}connection.o \
	${common_src}communicator_base.o \
	${common_src}communicator_structured.o \
	${common_src}communicator_unstructured.o \
	${common_src}discretization.o \
	${common_src}dm_kludge.o \
	${common_src}grid.o \
	${common_src}grid_structured.o \
	${common_src}grid_unstructured_cell.o \
	${common_src}grid_unstructured_explicit.o \
	${common_src}grid_unstructured_aux.o \
	${common_src}grid_unstructured.o \
	${common_src}grid_unstructured_polyhedra.o

relations_obj = \
	${common_src}saturation.o \
	${common_src}saturation_function.o \
	${common_src}characteristic_curves.o \
	${common_src}mphase_pckr_mod.o

properties_obj = \
	${common_src}fluid.o \
	${common_src}material.o 

dataset_obj = \
	${common_src}dataset_ascii.o \
	${common_src}dataset_base.o \
	${common_src}dataset_common_hdf5.o \
	${common_src}dataset_global_hdf5.o \
	${common_src}dataset_map_hdf5.o \
	${common_src}dataset.o \
	${common_src}dataset_gridded_hdf5.o \
	${common_src}time_storage.o 

srcsink_obj = \
	${common_src}srcsink_sandbox.o \
	${common_src}srcsink_sandbox_base.o \
	${common_src}srcsink_sandbox_mass_rate.o \
	${common_src}srcsink_sandbox_downreg.o \
	${common_src}srcsink_sandbox_wipp_gas.o \
	${common_src}srcsink_sandbox_wipp_well.o

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
	${common_src}e4d_mat_inv.o \
	${common_src}e4d_run.o \
	${common_src}e4d_setup.o \
	${common_src}e4d_vars.o \
	${common_src}factory_hydrogeophysics.o \
	${common_src}factory_pflotran.o \
	${common_src}factory_subsurface.o \
	${common_src}field.o \
	${common_src}geometry.o \
	${common_src}hydrostatic.o \
	${common_src}init_common.o \
	${common_src}init_geomechanics.o \
	${common_src}init_subsurface.o \
	${common_src}init_subsurface_flow.o \
	${common_src}init_subsurface_transport.o \
	${common_src}init_surface.o \
	${common_src}integral_flux.o \
	${common_src}lookup_table.o \
	${common_src}mass_transfer.o \
	${common_src}matrix_buffer.o \
	${common_src}multisimulation.o \
	${common_src}observation.o \
	${common_src}patch.o \
	${common_src}pmc_base.o \
	${common_src}pmc_hydrogeophysics.o \
	${common_src}pmc_material.o \
	${common_src}pmc_subsurface.o \
	${common_src}pmc_third_party.o \
	${common_src}realization.o \
	${common_src}realization_base.o \
	${common_src}region.o \
	${common_src}regression.o \
	${common_src}secondary_continuum_aux.o\
	${common_src}secondary_continuum.o \
	${common_src}simulation_aux.o \
	${common_src}simulation_base.o \
	${common_src}simulation_hydrogeophysics.o \
	${common_src}simulation_subsurface.o \
	${common_src}solver.o \
	${common_src}spline.o \
	${common_src}strata.o \
	${common_src}timestepper_base.o \
	${common_src}timestepper_BE.o \
	${common_src}transport.o \
	${common_src}uniform_velocity.o \
	${common_src}waypoint.o \
	${common_src}wipp.o \
	${common_src}wrapper_hydrogeophysics.o

surface_obj = \
	${common_src}checkpoint_surface.o \
	${common_src}factory_surface.o \
	${common_src}factory_surfsubsurface.o \
	${common_src}pmc_surface.o \
	${common_src}simulation_surface.o \
	${common_src}simulation_surfsubsurface.o \
	${common_src}surface_auxiliary.o \
	${common_src}surface_field.o \
	${common_src}surface_flow.o \
	${common_src}surface_init.o \
	${common_src}surface_global_aux.o \
	${common_src}surface_global.o \
	${common_src}surface_material.o \
	${common_src}surface_realization.o \
	${common_src}surface_th_aux.o \
	${common_src}surface_th.o \
	${common_src}timestepper_surface.o

geomech_obj = \
	${common_src}factory_geomechanics.o \
	${common_src}gauss.o \
	${common_src}geomechanics_auxiliary.o \
	${common_src}geomechanics_init.o \
	${common_src}geomechanics_realization.o \
	${common_src}geomechanics_discretization.o \
	${common_src}geomech_grid.o \
	${common_src}geomech_grid_aux.o \
	${common_src}geomechanics_condition.o \
	${common_src}geomechanics_coupler.o \
	${common_src}geomechanics_debug.o \
	${common_src}geomechanics_field.o \
	${common_src}geomechanics_force.o \
	${common_src}geomechanics_global.o \
	${common_src}geomechanics_global_aux.o\
	${common_src}geomechanics_logging.o \
	${common_src}geomechanics_material.o \
	${common_src}geomechanics_patch.o \
	${common_src}geomechanics_region.o \
	${common_src}geomechanics_strata.o \
	${common_src}output_geomechanics.o \
	${common_src}pm_geomechanics_force.o \
	${common_src}pmc_geomechanics.o \
	${common_src}simulation_geomechanics.o \
	${common_src}shape_function.o \
	${common_src}timestepper_geomechanics.o

# End Source Block

pflotran_obj = $(util_obj) $(eos_obj) $(mode_aux_obj) $(mode_obj) \
	$(shared_mode_aux_obj) $(chem_obj) $(geomech_obj) $(grid_obj) \
	$(relations_obj) $(properties_obj) $(dataset_obj) $(srcsink_obj) \
	$(io_obj) $(misc_obj) $(surface_obj) ${pflotran_src}pflotran.o

pflotran_rxn_obj = ${pflotran_src}pflotran_rxn.o

pflotran : $(pflotran_obj)
	${FLINKER} -o pflotran $(pflotran_obj) ${PETSC_LIB} ${LIBS} 

libpflotranchem.a : $(util_obj) $(shared_mode_aux_obj) $(chem_obj)
	$(AR) $(ARFLAGS) $@ $(util_obj) $(shared_mode_aux_obj) $(chem_obj)

pflotran_rxn : libpflotranchem.a $(pflotran_rxn_obj)
	${FLINKER} -o $@ $(pflotran_rxn_obj) -L. -lpflotranchem ${PETSC_LIB} ${LIBS} 

pflotran_provenance.F90 : FORCE
ifeq ($(UPDATE_PROVENANCE),1)
	@-rm pflotran_provenance.o pflotran_provenance_module.mod pflotran_provenance.F90
	$(PYTHON) ./pflotran-provenance.py \
		--pflotran-fflags "${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS}" \
		--petsc-config "${CONFIGURE_OPTIONS}"
else
	@if test ! -f pflotran_provenance.F90; then cp pflotran_no_provenance.F90 pflotran_provenance.F90; fi
endif

# do these files exist anymore?
#test: functions.o test.o
#	${FLINKER} -o test test.o functions.o $(LIBS) $(PETSC_LIB)

# Should add this to default PETSc targets as well?
%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<

# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
auxiliary.o : flash2_aux.o general_aux.o global_aux.o immis_aux.o \
              material_aux.o miscible_aux.o mphase_aux.o pflotran_constants.o \
              reactive_transport_aux.o richards_aux.o \
              secondary_continuum_aux.o th_aux.o
block_solve.o :
block_tridiag.o : pflotran_constants.o
characteristic_curves.o : input_aux.o option.o pflotran_constants.o string.o \
                          utility.o
checkpoint.o : discretization.o field.o global.o grid.o material.o option.o \
               pflotran_constants.o realization.o string.o variables.o
checkpoint_surface.o : discretization.o grid.o option.o output_aux.o \
                       pflotran_constants.o surface_field.o \
                       surface_realization.o
co2_span_wagner.o : co2_sw_rtsafe.o option.o pflotran_constants.o
co2_span_wagner_spline.o : pflotran_constants.o spline.o
co2_sw.o : co2_span_wagner.o co2_span_wagner_spline.o option.o \
           pflotran_constants.o
co2_sw_rtsafe.o : pflotran_constants.o
co2eos.o : pflotran_constants.o
communicator_base.o : dm_kludge.o option.o pflotran_constants.o
communicator_structured.o : communicator_base.o dm_kludge.o grid_structured.o \
                            pflotran_constants.o
communicator_unstructured.o : communicator_base.o dm_kludge.o \
                              grid_unstructured.o grid_unstructured_aux.o \
                              grid_unstructured_explicit.o \
                              pflotran_constants.o
condition.o : dataset.o dataset_ascii.o dataset_base.o general_aux.o \
              global_aux.o hdf5_aux.o input_aux.o logging.o option.o \
              pflotran_constants.o reaction_aux.o string.o time_storage.o \
              transport_constraint.o units.o utility.o
condition_control.o : condition.o connection.o coupler.o dataset_base.o \
                      dataset_common_hdf5.o dataset_gridded_hdf5.o \
                      discretization.o eos_water.o field.o general_aux.o \
                      global.o global_aux.o grid.o hdf5.o material_aux.o \
                      option.o patch.o pflotran_constants.o reaction.o \
                      reaction_aux.o reactive_transport.o \
                      reactive_transport_aux.o realization.o region.o \
                      surface_field.o surface_global_aux.o \
                      surface_realization.o surface_th_aux.o \
                      transport_constraint.o variables.o
connection.o : pflotran_constants.o utility.o
convergence.o : grid.o option.o pflotran_constants.o solver.o
coupler.o : condition.o connection.o dataset_base.o dataset_gridded_hdf5.o \
            grid.o grid_unstructured_aux.o grid_unstructured_explicit.o \
            input_aux.o option.o pflotran_constants.o region.o string.o \
            utility.o
dataset.o : dataset_ascii.o dataset_base.o dataset_common_hdf5.o \
            dataset_global_hdf5.o dataset_gridded_hdf5.o dataset_map_hdf5.o \
            dm_kludge.o input_aux.o option.o pflotran_constants.o string.o \
            time_storage.o
dataset_ascii.o : dataset_base.o input_aux.o option.o pflotran_constants.o \
                  string.o time_storage.o units.o utility.o
dataset_base.o : option.o pflotran_constants.o string.o time_storage.o \
                 utility.o
dataset_common_hdf5.o : dataset_base.o hdf5_aux.o input_aux.o logging.o \
                        option.o pflotran_constants.o string.o time_storage.o \
                        units.o
dataset_global_hdf5.o : dataset_base.o dataset_common_hdf5.o dm_kludge.o \
                        logging.o option.o pflotran_constants.o time_storage.o \
                        utility.o
dataset_gridded_hdf5.o : dataset_base.o dataset_common_hdf5.o logging.o \
                         option.o pflotran_constants.o string.o time_storage.o \
                         units.o utility.o
dataset_map_hdf5.o : dataset_base.o dataset_common_hdf5.o input_aux.o \
                     logging.o option.o pflotran_constants.o string.o \
                     time_storage.o units.o utility.o
debug.o : input_aux.o option.o pflotran_constants.o
discretization.o : dm_kludge.o grid.o grid_structured.o grid_unstructured.o \
                   grid_unstructured_aux.o grid_unstructured_explicit.o \
                   grid_unstructured_polyhedra.o input_aux.o option.o \
                   pflotran_constants.o string.o
dm_kludge.o : grid_unstructured_aux.o pflotran_constants.o
e4d_mat_inv.o :
e4d_run.o : e4d_vars.o
e4d_setup.o : e4d_mat_inv.o e4d_vars.o
e4d_vars.o :
eos.o : eos_gas.o eos_water.o input_aux.o option.o pflotran_constants.o \
        string.o
eos_gas.o : pflotran_constants.o
eos_water.o : pflotran_constants.o
factory_geomechanics.o : condition_control.o factory_subsurface.o \
                         geomechanics_discretization.o geomechanics_force.o \
                         geomechanics_realization.o init_common.o \
                         init_geomechanics.o input_aux.o logging.o option.o \
                         output.o output_aux.o pflotran_constants.o pm_base.o \
                         pm_base_pointer.o pm_geomechanics_force.o pmc_base.o \
                         pmc_geomechanics.o simulation_aux.o simulation_base.o \
                         simulation_geomechanics.o simulation_subsurface.o \
                         timestepper_base.o timestepper_geomechanics.o
factory_hydrogeophysics.o : discretization.o factory_subsurface.o input_aux.o \
                            logging.o option.o pflotran_constants.o pm_base.o \
                            pmc_hydrogeophysics.o simulation_base.o \
                            simulation_hydrogeophysics.o string.o \
                            wrapper_hydrogeophysics.o
factory_pflotran.o : eos.o factory_hydrogeophysics.o factory_subsurface.o \
                     factory_surfsubsurface.o input_aux.o logging.o \
                     multisimulation.o option.o pflotran_constants.o pm_base.o \
                     pm_surface_th.o pmc_base.o simulation_base.o string.o
factory_subsurface.o : condition_control.o general.o global.o init_common.o \
                       init_subsurface.o init_subsurface_flow.o \
                       init_subsurface_transport.o input_aux.o logging.o \
                       option.o output.o output_aux.o pflotran_constants.o \
                       pm_base.o pm_base_pointer.o pm_flash2.o pm_general.o \
                       pm_immis.o pm_miscible.o pm_mphase.o pm_richards.o \
                       pm_rt.o pm_subsurface.o pm_th.o pm_waste_form.o \
                       pmc_base.o pmc_material.o pmc_subsurface.o \
                       pmc_third_party.o reactive_transport.o realization.o \
                       regression.o simulation_base.o simulation_subsurface.o \
                       solver.o strata.o string.o timestepper_BE.o \
                       timestepper_base.o waypoint.o
factory_surface.o : checkpoint_surface.o condition_control.o init_common.o \
                    input_aux.o logging.o option.o output.o output_aux.o \
                    output_surface.o pflotran_constants.o simulation_base.o \
                    simulation_surface.o surface_realization.o \
                    timestepper_base.o timestepper_surface.o
factory_surfsubsurface.o : coupler.o discretization.o dm_kludge.o \
                           factory_subsurface.o factory_surface.o field.o \
                           grid.o grid_unstructured.o grid_unstructured_aux.o \
                           grid_unstructured_cell.o init_common.o \
                           init_surface.o input_aux.o logging.o option.o \
                           patch.o pflotran_constants.o pm_base.o \
                           pm_base_pointer.o pm_surface.o pm_surface_th.o \
                           pmc_base.o pmc_surface.o realization.o region.o \
                           simulation_aux.o simulation_base.o \
                           simulation_subsurface.o simulation_surface.o \
                           simulation_surfsubsurface.o string.o \
                           surface_field.o surface_flow.o surface_init.o \
                           surface_realization.o surface_th.o \
                           timestepper_surface.o waypoint.o
field.o : pflotran_constants.o
flash2.o : co2_span_wagner.o co2_span_wagner_spline.o co2_sw.o co2eos.o \
           connection.o coupler.o debug.o discretization.o eos_water.o field.o \
           flash2_aux.o global_aux.o grid.o logging.o material.o \
           material_aux.o option.o output_aux.o patch.o pflotran_constants.o \
           realization.o variables.o
flash2_aux.o : co2_span_wagner.o co2_span_wagner_spline.o co2_sw.o co2eos.o \
               eos_water.o fluid.o gas_eos_mod.o global_aux.o \
               mphase_pckr_mod.o option.o pflotran_constants.o \
               saturation_function.o
fluid.o : input_aux.o option.o pflotran_constants.o string.o
gas_eos_mod.o : pflotran_constants.o
gauss.o : grid_unstructured_cell.o pflotran_constants.o utility.o
general.o : connection.o coupler.o debug.o discretization.o eos_gas.o \
            eos_water.o field.o fluid.o general_aux.o global_aux.o grid.o \
            input_aux.o material.o material_aux.o option.o output_aux.o \
            patch.o pflotran_constants.o realization.o saturation_function.o \
            srcsink_sandbox.o srcsink_sandbox_base.o \
            srcsink_sandbox_wipp_gas.o srcsink_sandbox_wipp_well.o string.o \
            variables.o
general_aux.o : characteristic_curves.o eos_gas.o eos_water.o global_aux.o \
                material_aux.o option.o pflotran_constants.o string.o \
                utility.o wipp.o
geomech_grid.o : gauss.o geomech_grid_aux.o geomechanics_region.o \
                 grid_unstructured_aux.o grid_unstructured_cell.o input_aux.o \
                 logging.o option.o pflotran_constants.o utility.o
geomech_grid_aux.o : gauss.o grid_unstructured_cell.o option.o \
                     pflotran_constants.o utility.o
geomechanics_auxiliary.o : geomechanics_global_aux.o pflotran_constants.o
geomechanics_condition.o : condition.o dataset.o dataset_ascii.o \
                           dataset_base.o geomechanics_logging.o input_aux.o \
                           option.o pflotran_constants.o string.o \
                           time_storage.o
geomechanics_coupler.o : geomechanics_condition.o geomechanics_region.o \
                         input_aux.o option.o pflotran_constants.o string.o
geomechanics_debug.o : input_aux.o option.o pflotran_constants.o
geomechanics_discretization.o : geomech_grid.o geomech_grid_aux.o option.o \
                                pflotran_constants.o
geomechanics_field.o : pflotran_constants.o
geomechanics_force.o : discretization.o field.o geomech_grid.o \
                       geomech_grid_aux.o geomechanics_auxiliary.o \
                       geomechanics_coupler.o geomechanics_debug.o \
                       geomechanics_discretization.o geomechanics_field.o \
                       geomechanics_global_aux.o geomechanics_logging.o \
                       geomechanics_material.o geomechanics_patch.o \
                       geomechanics_realization.o geomechanics_region.o grid.o \
                       grid_unstructured_cell.o material.o material_aux.o \
                       option.o output_aux.o patch.o pflotran_constants.o \
                       realization.o shape_function.o string.o utility.o \
                       variables.o
geomechanics_global.o : geomech_grid.o geomech_grid_aux.o \
                        geomechanics_coupler.o geomechanics_discretization.o \
                        geomechanics_field.o geomechanics_global_aux.o \
                        geomechanics_patch.o geomechanics_realization.o \
                        option.o pflotran_constants.o variables.o
geomechanics_global_aux.o : option.o pflotran_constants.o utility.o
geomechanics_init.o : convergence.o discretization.o geomech_grid.o \
                      geomech_grid_aux.o geomechanics_condition.o \
                      geomechanics_coupler.o geomechanics_debug.o \
                      geomechanics_discretization.o geomechanics_field.o \
                      geomechanics_force.o geomechanics_global.o \
                      geomechanics_material.o geomechanics_patch.o \
                      geomechanics_realization.o geomechanics_region.o \
                      geomechanics_strata.o grid_unstructured.o \
                      grid_unstructured_aux.o input_aux.o option.o \
                      output_aux.o output_tecplot.o patch.o \
                      pflotran_constants.o realization.o \
                      simulation_geomechanics.o solver.o string.o units.o \
                      waypoint.o
geomechanics_logging.o : pflotran_constants.o
geomechanics_material.o : input_aux.o option.o pflotran_constants.o string.o
geomechanics_patch.o : dataset_base.o geomech_grid.o geomech_grid_aux.o \
                       geomechanics_auxiliary.o geomechanics_condition.o \
                       geomechanics_coupler.o geomechanics_field.o \
                       geomechanics_global_aux.o geomechanics_material.o \
                       geomechanics_region.o geomechanics_strata.o option.o \
                       output_aux.o pflotran_constants.o variables.o
geomechanics_realization.o : dataset.o dataset_base.o geomech_grid.o \
                             geomech_grid_aux.o geomechanics_condition.o \
                             geomechanics_coupler.o geomechanics_debug.o \
                             geomechanics_discretization.o \
                             geomechanics_field.o geomechanics_material.o \
                             geomechanics_patch.o geomechanics_region.o \
                             geomechanics_strata.o grid.o input_aux.o option.o \
                             output_aux.o pflotran_constants.o realization.o \
                             realization_base.o string.o time_storage.o \
                             waypoint.o
geomechanics_region.o : geometry.o grid_unstructured_cell.o input_aux.o \
                        logging.o option.o pflotran_constants.o string.o \
                        utility.o
geomechanics_strata.o : geomechanics_material.o geomechanics_region.o \
                        input_aux.o option.o pflotran_constants.o string.o
geometry.o : input_aux.o option.o pflotran_constants.o
global.o : communicator_base.o connection.o coupler.o field.o global_aux.o \
           grid.o material.o option.o patch.o pflotran_constants.o \
           realization.o realization_base.o variables.o
global_aux.o : option.o pflotran_constants.o utility.o
grid.o : connection.o dm_kludge.o grid_structured.o grid_unstructured.o \
         grid_unstructured_aux.o grid_unstructured_explicit.o \
         grid_unstructured_polyhedra.o logging.o option.o pflotran_constants.o \
         region.o utility.o
grid_structured.o : connection.o input_aux.o option.o pflotran_constants.o \
                    string.o utility.o
grid_unstructured.o : connection.o geometry.o grid_unstructured_aux.o \
                      grid_unstructured_cell.o hdf5_aux.o input_aux.o option.o \
                      pflotran_constants.o string.o utility.o
grid_unstructured_aux.o : geometry.o grid_unstructured_cell.o option.o \
                          pflotran_constants.o utility.o
grid_unstructured_cell.o : option.o pflotran_constants.o utility.o
grid_unstructured_explicit.o : connection.o geometry.o grid_unstructured_aux.o \
                               input_aux.o option.o pflotran_constants.o \
                               string.o utility.o
grid_unstructured_polyhedra.o : connection.o geometry.o grid_unstructured.o \
                                grid_unstructured_aux.o \
                                grid_unstructured_cell.o input_aux.o option.o \
                                pflotran_constants.o string.o utility.o
hdf5.o : discretization.o field.o grid.o grid_unstructured_cell.o hdf5_aux.o \
         logging.o option.o patch.o pflotran_constants.o realization.o \
         realization_base.o region.o utility.o
hdf5_aux.o : input_aux.o logging.o option.o pflotran_constants.o region.o \
             string.o
hydrostatic.o : condition.o connection.o coupler.o dataset_ascii.o \
                dataset_gridded_hdf5.o eos_water.o general_aux.o grid.o \
                grid_structured.o option.o pflotran_constants.o region.o \
                utility.o
immis.o : co2_span_wagner.o co2_span_wagner_spline.o co2_sw.o co2eos.o \
          connection.o coupler.o debug.o discretization.o eos_water.o field.o \
          global_aux.o grid.o immis_aux.o material.o material_aux.o option.o \
          output_aux.o patch.o pflotran_constants.o realization.o variables.o
immis_aux.o : co2_span_wagner.o co2_span_wagner_spline.o co2_sw.o co2eos.o \
              eos_water.o fluid.o gas_eos_mod.o mphase_pckr_mod.o option.o \
              pflotran_constants.o saturation_function.o
init_common.o : condition.o connection.o coupler.o discretization.o field.o \
                general.o grid.o hdf5.o hdf5_aux.o input_aux.o logging.o \
                option.o output.o output_tecplot.o patch.o \
                pflotran_constants.o reaction.o reaction_aux.o realization.o \
                region.o string.o
init_geomechanics.o : convergence.o discretization.o geomech_grid.o \
                      geomech_grid_aux.o geomechanics_discretization.o \
                      geomechanics_field.o geomechanics_force.o \
                      geomechanics_global.o geomechanics_material.o \
                      geomechanics_patch.o geomechanics_realization.o \
                      geomechanics_region.o geomechanics_strata.o option.o \
                      pflotran_constants.o realization.o \
                      simulation_geomechanics.o solver.o waypoint.o
init_subsurface.o : characteristic_curves.o condition.o coupler.o dataset.o \
                    dataset_base.o dataset_common_hdf5.o debug.o \
                    discretization.o eos.o eos_water.o field.o fluid.o \
                    general.o grid.o grid_structured.o grid_unstructured_aux.o \
                    hdf5.o hdf5_aux.o init_common.o input_aux.o \
                    integral_flux.o logging.o mass_transfer.o material.o \
                    material_aux.o observation.o option.o output_aux.o \
                    output_tecplot.o patch.o pflotran_constants.o \
                    pmc_subsurface.o reaction.o reaction_aux.o \
                    reaction_database.o reaction_mineral.o \
                    reaction_solid_solution.o realization.o region.o \
                    regression.o saturation_function.o simulation_subsurface.o \
                    solver.o srcsink_sandbox.o strata.o string.o \
                    timestepper_BE.o transport_constraint.o uniform_velocity.o \
                    units.o variables.o waypoint.o wipp.o
init_subsurface_flow.o : co2_sw.o condition_control.o convergence.o \
                         discretization.o field.o flash2.o general.o grid.o \
                         hdf5.o immis.o init_common.o material.o miscible.o \
                         mphase.o option.o patch.o pflotran_constants.o \
                         realization.o richards.o solver.o th.o
init_subsurface_transport.o : condition_control.o convergence.o \
                              discretization.o global.o init_common.o option.o \
                              pflotran_constants.o reactive_transport.o \
                              realization.o secondary_continuum.o solver.o \
                              variables.o
init_surface.o : condition_control.o convergence.o discretization.o \
                 eos_water.o field.o grid.o hdf5.o material.o option.o patch.o \
                 pflotran_constants.o realization.o region.o solver.o strata.o \
                 surface_field.o surface_flow.o surface_global.o \
                 surface_material.o surface_realization.o surface_th.o \
                 timestepper_base.o waypoint.o
input_aux.o : option.o pflotran_constants.o string.o
integral_flux.o : geometry.o input_aux.o option.o pflotran_constants.o \
                  string.o utility.o
logging.o : pflotran_constants.o
lookup_table.o : pflotran_constants.o utility.o
mass_transfer.o : dataset_base.o dataset_common_hdf5.o dataset_global_hdf5.o \
                  discretization.o grid.o input_aux.o option.o \
                  pflotran_constants.o string.o
material.o : characteristic_curves.o communicator_base.o dataset_common_hdf5.o \
             field.o global_aux.o input_aux.o material_aux.o option.o \
             pflotran_constants.o string.o variables.o
material_aux.o : option.o pflotran_constants.o utility.o variables.o
matrix_block_aux.o : option.o pflotran_constants.o utility.o
matrix_buffer.o : grid.o grid_structured.o pflotran_constants.o
miscible.o : co2_span_wagner.o co2_span_wagner_spline.o co2_sw.o co2eos.o \
             connection.o coupler.o debug.o discretization.o eos_water.o \
             field.o global_aux.o grid.o logging.o material.o miscible_aux.o \
             option.o output_aux.o patch.o pflotran_constants.o realization.o \
             variables.o
miscible_aux.o : fluid.o global_aux.o option.o pflotran_constants.o
mphase.o : co2_span_wagner.o co2_span_wagner_spline.o co2_sw.o co2eos.o \
           connection.o coupler.o debug.o discretization.o eos_water.o field.o \
           gas_eos_mod.o global_aux.o grid.o material.o material_aux.o \
           mphase_aux.o option.o output_aux.o patch.o pflotran_constants.o \
           realization.o secondary_continuum.o secondary_continuum_aux.o \
           variables.o
mphase_aux.o : co2_span_wagner.o co2_span_wagner_spline.o co2_sw.o co2eos.o \
               eos_water.o fluid.o gas_eos_mod.o global_aux.o \
               mphase_pckr_mod.o option.o pflotran_constants.o \
               saturation_function.o
mphase_pckr_mod.o : material.o option.o pflotran_constants.o \
                    saturation_function.o
multisimulation.o : input_aux.o option.o pflotran_constants.o
observation.o : connection.o input_aux.o option.o pflotran_constants.o \
                region.o string.o
option.o : logging.o option_flow.o option_transport.o pflotran_constants.o
option_flow.o : pflotran_constants.o
option_transport.o : pflotran_constants.o
output.o : connection.o coupler.o debug.o discretization.o field.o \
           general_aux.o grid.o hdf5.o input_aux.o logging.o option.o \
           output_aux.o output_common.o output_hdf5.o output_observation.o \
           output_tecplot.o output_vtk.o patch.o pflotran_constants.o \
           reaction_aux.o realization_base.o variables.o
output_aux.o : input_aux.o option.o pflotran_constants.o string.o variables.o
output_common.o : connection.o coupler.o field.o global_aux.o grid.o \
                  grid_unstructured_aux.o grid_unstructured_cell.o hdf5.o \
                  hdf5_aux.o logging.o option.o output_aux.o patch.o \
                  pflotran_constants.o realization_base.o richards_aux.o \
                  variables.o
output_geomechanics.o : geomech_grid.o geomech_grid_aux.o \
                        geomechanics_discretization.o geomechanics_field.o \
                        geomechanics_logging.o geomechanics_patch.o \
                        geomechanics_realization.o grid_unstructured_cell.o \
                        hdf5.o hdf5_aux.o option.o output_aux.o \
                        output_common.o output_hdf5.o output_tecplot.o \
                        pflotran_constants.o string.o variables.o
output_hdf5.o : connection.o coupler.o discretization.o field.o grid.o \
                grid_unstructured_aux.o grid_unstructured_cell.o hdf5.o \
                hdf5_aux.o input_aux.o logging.o option.o output_aux.o \
                output_common.o patch.o pflotran_constants.o \
                pflotran_provenance.o reaction_aux.o realization_base.o \
                string.o variables.o
output_observation.o : connection.o coupler.o discretization.o field.o \
                       flash2.o general.o global_aux.o grid.o \
                       grid_structured.o immis.o integral_flux.o logging.o \
                       material_aux.o miscible.o mphase.o observation.o \
                       option.o output_aux.o output_common.o patch.o \
                       pflotran_constants.o reaction_aux.o \
                       reactive_transport.o reactive_transport_aux.o \
                       realization.o realization_base.o region.o richards.o \
                       th.o utility.o variables.o
output_surface.o : connection.o coupler.o discretization.o field.o grid.o \
                   grid_unstructured_aux.o grid_unstructured_cell.o hdf5.o \
                   hdf5_aux.o input_aux.o logging.o option.o output_aux.o \
                   output_common.o output_hdf5.o output_tecplot.o patch.o \
                   pflotran_constants.o reaction_aux.o realization.o string.o \
                   surface_field.o surface_realization.o utility.o variables.o
output_tecplot.o : connection.o discretization.o field.o grid.o \
                   grid_structured.o grid_unstructured_aux.o \
                   grid_unstructured_cell.o logging.o observation.o option.o \
                   output_aux.o output_common.o patch.o pflotran_constants.o \
                   reaction_aux.o realization_base.o region.o \
                   secondary_continuum_aux.o string.o utility.o variables.o
output_vtk.o : discretization.o field.o grid.o grid_structured.o logging.o \
               option.o output_aux.o output_common.o patch.o \
               pflotran_constants.o reaction_aux.o realization_base.o string.o \
               variables.o
patch.o : auxiliary.o characteristic_curves.o condition.o connection.o \
          coupler.o dataset.o dataset_ascii.o dataset_base.o \
          dataset_common_hdf5.o dataset_gridded_hdf5.o dataset_map_hdf5.o \
          eos_water.o field.o general_aux.o geometry.o global_aux.o grid.o \
          hydrostatic.o immis_aux.o integral_flux.o material.o material_aux.o \
          miscible_aux.o mphase_aux.o observation.o option.o output_aux.o \
          pflotran_constants.o reaction.o reaction_aux.o reaction_mineral.o \
          reaction_mineral_aux.o reaction_surf_complex_aux.o \
          reactive_transport_aux.o region.o richards_aux.o saturation.o \
          saturation_function.o strata.o surface_auxiliary.o surface_field.o \
          surface_material.o th_aux.o transport_constraint.o utility.o \
          variables.o
pflotran.o : factory_geomechanics.o factory_hydrogeophysics.o \
             factory_pflotran.o factory_subsurface.o factory_surface.o \
             factory_surfsubsurface.o multisimulation.o option.o output_aux.o \
             pflotran_constants.o pflotran_provenance.o pm_base.o \
             simulation_base.o
pflotran_constants.o :
pflotran_provenance.o :
pm_base.o : option.o output_aux.o pflotran_constants.o realization_base.o \
            solver.o
pm_base_pointer.o : option.o pflotran_constants.o pm_base.o realization.o
pm_flash2.o : flash2.o pflotran_constants.o pm_base.o pm_subsurface.o
pm_general.o : field.o general.o general_aux.o global.o global_aux.o grid.o \
               material.o option.o pflotran_constants.o pm_base.o \
               pm_subsurface.o realization.o realization_base.o solver.o \
               variables.o
pm_geomechanics_force.o : checkpoint.o communicator_base.o \
                          communicator_structured.o \
                          communicator_unstructured.o condition.o \
                          geomechanics_discretization.o geomechanics_force.o \
                          geomechanics_realization.o global.o grid.o option.o \
                          pflotran_constants.o pm_base.o
pm_immis.o : immis.o pflotran_constants.o pm_base.o pm_subsurface.o
pm_miscible.o : miscible.o pflotran_constants.o pm_base.o pm_subsurface.o
pm_mphase.o : mphase.o pflotran_constants.o pm_base.o pm_subsurface.o
pm_richards.o : pflotran_constants.o pm_base.o pm_subsurface.o richards.o \
                solver.o
pm_rt.o : communicator_base.o communicator_structured.o \
          communicator_unstructured.o condition.o condition_control.o \
          discretization.o field.o global.o grid.o integral_flux.o \
          mass_transfer.o material.o option.o pflotran_constants.o pm_base.o \
          reaction_aux.o reactive_transport.o realization.o realization_base.o \
          variables.o
pm_subsurface.o : checkpoint.o communicator_base.o communicator_structured.o \
                  communicator_unstructured.o condition.o condition_control.o \
                  discretization.o global.o grid.o integral_flux.o material.o \
                  material_aux.o option.o pflotran_constants.o pm_base.o \
                  realization.o solver.o variables.o
pm_surface.o : checkpoint_surface.o communicator_base.o \
               communicator_unstructured.o condition.o discretization.o \
               global.o grid.o option.o pflotran_constants.o pm_base.o \
               surface_realization.o
pm_surface_flow.o : condition.o discretization.o pflotran_constants.o \
                    pm_base.o pm_surface.o surface_field.o surface_flow.o
pm_surface_th.o : communicator_base.o condition.o discretization.o grid.o \
                  option.o pflotran_constants.o pm_base.o pm_surface.o \
                  realization.o surface_field.o surface_realization.o \
                  surface_th.o
pm_th.o : communicator_base.o communicator_structured.o \
          communicator_unstructured.o discretization.o eos_water.o global.o \
          grid.o input_aux.o option.o pflotran_constants.o pm_base.o \
          pm_subsurface.o realization.o solver.o string.o th.o utility.o
pm_waste_form.o : geometry.o global.o global_aux.o grid.o grid_structured.o \
                  grid_unstructured.o input_aux.o option.o \
                  pflotran_constants.o pm_base.o reaction_aux.o \
                  reactive_transport_aux.o realization.o string.o utility.o
pmc_base.o : checkpoint.o logging.o option.o output.o output_aux.o \
             pflotran_constants.o pm_base.o pm_base_pointer.o simulation_aux.o \
             timestepper_base.o utility.o waypoint.o
pmc_geomechanics.o : geomechanics_discretization.o geomechanics_force.o \
                     geomechanics_realization.o grid.o option.o \
                     pflotran_constants.o pm_base.o pmc_base.o realization.o \
                     timestepper_base.o
pmc_hydrogeophysics.o : option.o pflotran_constants.o pmc_base.o realization.o \
                        realization_base.o string.o timestepper_base.o \
                        utility.o variables.o wrapper_hydrogeophysics.o
pmc_material.o : init_subsurface.o option.o pflotran_constants.o pmc_base.o \
                 realization.o timestepper_base.o
pmc_subsurface.o : connection.o coupler.o discretization.o eos_water.o field.o \
                   grid.o material.o material_aux.o option.o patch.o \
                   pflotran_constants.o pm_base.o pm_base_pointer.o pmc_base.o \
                   realization.o realization_base.o string.o timestepper_BE.o \
                   timestepper_base.o variables.o
pmc_surface.o : connection.o coupler.o eos_water.o grid.o option.o output.o \
                output_surface.o patch.o pflotran_constants.o pm_base.o \
                pm_surface_flow.o pmc_base.o realization.o string.o \
                surface_flow.o surface_global_aux.o surface_realization.o \
                surface_th.o surface_th_aux.o timestepper_base.o \
                timestepper_surface.o
pmc_third_party.o : option.o pflotran_constants.o pmc_base.o realization.o \
                    timestepper_base.o utility.o
ifdef PFLOTRAN_RXN_FLAG
reaction.o : global_aux.o input_aux.o material_aux.o option.o \
             pflotran_constants.o reaction_aux.o reaction_clm.o \
             reaction_immobile.o reaction_immobile_aux.o reaction_microbial.o \
             reaction_microbial_aux.o reaction_mineral.o \
             reaction_mineral_aux.o reaction_sandbox.o \
             reaction_solid_solution.o reaction_solid_solution_aux.o \
             reaction_surf_complex.o reaction_surf_complex_aux.o \
             reactive_transport_aux.o string.o transport_constraint.o units.o \
             utility.o variables.o
else
reaction.o : co2_span_wagner.o co2eos.o eos_water.o global_aux.o input_aux.o \
             material_aux.o option.o pflotran_constants.o reaction_aux.o \
             reaction_clm.o reaction_immobile.o reaction_immobile_aux.o \
             reaction_microbial.o reaction_microbial_aux.o reaction_mineral.o \
             reaction_mineral_aux.o reaction_sandbox.o \
             reaction_solid_solution.o reaction_solid_solution_aux.o \
             reaction_surf_complex.o reaction_surf_complex_aux.o \
             reactive_transport_aux.o string.o transport_constraint.o units.o \
             utility.o variables.o
endif
reaction_aux.o : option.o pflotran_constants.o reaction_database_aux.o \
                 reaction_immobile_aux.o reaction_microbial_aux.o \
                 reaction_mineral_aux.o reaction_solid_solution_aux.o \
                 reaction_surf_complex_aux.o string.o utility.o
reaction_clm.o : global_aux.o input_aux.o material_aux.o option.o \
                 pflotran_constants.o reaction_aux.o reaction_immobile_aux.o \
                 reactive_transport_aux.o string.o units.o utility.o
reaction_database.o : input_aux.o option.o pflotran_constants.o reaction.o \
                      reaction_aux.o reaction_clm.o reaction_database_aux.o \
                      reaction_immobile.o reaction_immobile_aux.o \
                      reaction_microbial.o reaction_microbial_aux.o \
                      reaction_mineral.o reaction_mineral_aux.o \
                      reaction_sandbox.o reaction_solid_solution.o \
                      reaction_surf_complex_aux.o string.o utility.o
reaction_database_aux.o : input_aux.o option.o pflotran_constants.o string.o \
                          utility.o
reaction_database_hpt.o : input_aux.o option.o pflotran_constants.o reaction.o \
                          reaction_aux.o reaction_database.o \
                          reaction_database_aux.o reaction_mineral_aux.o \
                          reaction_surf_complex_aux.o string.o utility.o
reaction_immobile.o : input_aux.o option.o pflotran_constants.o \
                      reaction_immobile_aux.o string.o utility.o
reaction_immobile_aux.o : option.o pflotran_constants.o \
                          reaction_database_aux.o string.o utility.o
reaction_microbial.o : global_aux.o input_aux.o material_aux.o option.o \
                       pflotran_constants.o reaction_aux.o \
                       reaction_immobile_aux.o reaction_microbial_aux.o \
                       reactive_transport_aux.o string.o utility.o
reaction_microbial_aux.o : pflotran_constants.o reaction_database_aux.o \
                           utility.o
reaction_mineral.o : global_aux.o input_aux.o material_aux.o option.o \
                     pflotran_constants.o reaction_aux.o \
                     reaction_database_aux.o reaction_mineral_aux.o \
                     reaction_solid_solution_aux.o reactive_transport_aux.o \
                     string.o units.o utility.o
reaction_mineral_aux.o : option.o pflotran_constants.o reaction_database_aux.o \
                         string.o utility.o
reaction_sandbox.o : global_aux.o input_aux.o material_aux.o option.o \
                     pflotran_constants.o reaction_aux.o \
                     reaction_sandbox_base.o reaction_sandbox_clm_cn.o \
                     reaction_sandbox_example.o reaction_sandbox_ufd_wp.o \
                     reactive_transport_aux.o string.o utility.o
reaction_sandbox_base.o : global_aux.o input_aux.o material_aux.o option.o \
                          pflotran_constants.o reaction_aux.o \
                          reactive_transport_aux.o
reaction_sandbox_clm_cn.o : global_aux.o input_aux.o material_aux.o option.o \
                            pflotran_constants.o reaction_aux.o \
                            reaction_immobile_aux.o reaction_sandbox_base.o \
                            reactive_transport_aux.o string.o units.o \
                            utility.o
reaction_sandbox_example.o : global_aux.o input_aux.o material_aux.o option.o \
                             pflotran_constants.o reaction_aux.o \
                             reaction_sandbox_base.o reactive_transport_aux.o \
                             string.o units.o
reaction_sandbox_ufd_wp.o : global_aux.o input_aux.o material_aux.o option.o \
                            pflotran_constants.o reaction_aux.o \
                            reaction_immobile_aux.o reaction_sandbox_base.o \
                            reactive_transport_aux.o string.o units.o
reaction_solid_solution.o : input_aux.o option.o pflotran_constants.o \
                            reaction_aux.o reaction_mineral.o \
                            reaction_mineral_aux.o \
                            reaction_solid_solution_aux.o string.o utility.o
reaction_solid_solution_aux.o : pflotran_constants.o reaction_mineral_aux.o
reaction_surf_complex.o : global_aux.o input_aux.o material_aux.o \
                          matrix_block_aux.o option.o pflotran_constants.o \
                          reaction_aux.o reaction_surf_complex_aux.o \
                          reactive_transport_aux.o string.o utility.o
reaction_surf_complex_aux.o : option.o pflotran_constants.o \
                              reaction_database_aux.o utility.o
reactive_transport.o : co2_span_wagner.o co2eos.o condition.o connection.o \
                       coupler.o debug.o discretization.o eos_water.o field.o \
                       fluid.o global.o global_aux.o grid.o logging.o \
                       mass_transfer.o material.o material_aux.o option.o \
                       output_aux.o output_tecplot.o patch.o \
                       pflotran_constants.o reaction.o reaction_aux.o \
                       reaction_surf_complex_aux.o reactive_transport_aux.o \
                       realization.o region.o secondary_continuum.o \
                       secondary_continuum_aux.o transport.o \
                       transport_constraint.o variables.o
reactive_transport_aux.o : matrix_block_aux.o option.o pflotran_constants.o \
                           reaction_aux.o reaction_surf_complex_aux.o \
                           utility.o
realization.o : characteristic_curves.o communicator_structured.o \
                communicator_unstructured.o condition.o coupler.o dataset.o \
                dataset_base.o dataset_common_hdf5.o debug.o discretization.o \
                dm_kludge.o field.o fluid.o grid.o grid_structured.o \
                grid_unstructured.o grid_unstructured_aux.o \
                grid_unstructured_cell.o input_aux.o mass_transfer.o \
                material.o material_aux.o option.o output_aux.o patch.o \
                pflotran_constants.o reaction.o reaction_aux.o \
                reactive_transport_aux.o realization_base.o region.o \
                saturation_function.o strata.o string.o time_storage.o \
                transport_constraint.o uniform_velocity.o variables.o \
                waypoint.o
realization_base.o : communicator_base.o debug.o discretization.o field.o \
                     input_aux.o mass_transfer.o option.o output_aux.o patch.o \
                     pflotran_constants.o reaction_aux.o
region.o : geometry.o grid_structured.o grid_unstructured_cell.o input_aux.o \
           logging.o option.o pflotran_constants.o string.o utility.o
regression.o : discretization.o grid.o input_aux.o option.o output.o \
               output_aux.o output_common.o pflotran_constants.o realization.o \
               string.o timestepper_BE.o utility.o
richards.o : connection.o coupler.o debug.o discretization.o eos_water.o \
             field.o global_aux.o grid.o logging.o mass_transfer.o material.o \
             material_aux.o matrix_buffer.o option.o output_aux.o patch.o \
             pflotran_constants.o realization.o realization_base.o \
             richards_aux.o richards_common.o saturation_function.o \
             srcsink_sandbox.o srcsink_sandbox_base.o \
             srcsink_sandbox_downreg.o string.o utility.o variables.o
richards_aux.o : eos_water.o global_aux.o material_aux.o matrix_buffer.o \
                 option.o pflotran_constants.o saturation_function.o utility.o
richards_common.o : connection.o eos_water.o global_aux.o material_aux.o \
                    option.o pflotran_constants.o richards_aux.o \
                    saturation_function.o utility.o
saturation.o : condition.o connection.o coupler.o grid.o option.o \
               pflotran_constants.o region.o saturation_function.o
saturation_function.o : input_aux.o option.o pflotran_constants.o string.o \
                        utility.o
secondary_continuum.o : block_solve.o block_tridiag.o condition.o coupler.o \
                        eos_water.o global_aux.o grid.o material.o \
                        material_aux.o mphase_aux.o option.o \
                        pflotran_constants.o reaction.o reaction_aux.o \
                        reactive_transport_aux.o realization.o \
                        secondary_continuum_aux.o string.o \
                        transport_constraint.o utility.o
secondary_continuum_aux.o : option.o pflotran_constants.o \
                            reactive_transport_aux.o utility.o
shape_function.o : gauss.o grid_unstructured_cell.o pflotran_constants.o
simulation_aux.o : option.o pflotran_constants.o
simulation_base.o : input_aux.o logging.o option.o output.o output_aux.o \
                    pflotran_constants.o pm_base.o pmc_base.o simulation_aux.o \
                    timestepper_base.o waypoint.o
simulation_geomechanics.o : geomechanics_realization.o option.o output.o \
                            pflotran_constants.o pmc_base.o pmc_geomechanics.o \
                            pmc_subsurface.o realization.o regression.o \
                            simulation_base.o simulation_subsurface.o \
                            timestepper_base.o
simulation_hydrogeophysics.o : option.o pflotran_constants.o \
                               pmc_hydrogeophysics.o simulation_base.o \
                               simulation_subsurface.o \
                               wrapper_hydrogeophysics.o
simulation_subsurface.o : logging.o option.o output.o output_aux.o \
                          pflotran_constants.o pmc_base.o pmc_subsurface.o \
                          pmc_third_party.o reaction_clm.o reaction_sandbox.o \
                          realization.o regression.o simulation_base.o \
                          srcsink_sandbox.o timestepper_BE.o \
                          timestepper_base.o wipp.o
simulation_surface.o : logging.o option.o output.o pflotran_constants.o \
                       pmc_base.o pmc_surface.o regression.o simulation_base.o \
                       surface_realization.o timestepper_base.o
simulation_surfsubsurface.o : logging.o option.o output.o pflotran_constants.o \
                              pmc_base.o pmc_subsurface.o pmc_surface.o \
                              realization.o regression.o simulation_aux.o \
                              simulation_base.o simulation_subsurface.o \
                              simulation_surface.o surface_realization.o \
                              timestepper_base.o
solver.o : input_aux.o option.o pflotran_constants.o string.o
spline.o : pflotran_constants.o
srcsink_sandbox.o : grid.o input_aux.o material_aux.o option.o \
                    pflotran_constants.o region.o srcsink_sandbox_base.o \
                    srcsink_sandbox_downreg.o srcsink_sandbox_mass_rate.o \
                    srcsink_sandbox_wipp_gas.o srcsink_sandbox_wipp_well.o \
                    string.o utility.o
srcsink_sandbox_base.o : input_aux.o material_aux.o option.o \
                         pflotran_constants.o region.o
srcsink_sandbox_downreg.o : condition.o dataset.o dataset_ascii.o \
                            dataset_base.o general_aux.o input_aux.o \
                            material_aux.o option.o pflotran_constants.o \
                            reaction_aux.o region.o srcsink_sandbox_base.o \
                            string.o units.o
srcsink_sandbox_mass_rate.o : general_aux.o input_aux.o material_aux.o \
                              option.o pflotran_constants.o reaction_aux.o \
                              region.o srcsink_sandbox_base.o string.o units.o
srcsink_sandbox_wipp_gas.o : eos_gas.o input_aux.o material_aux.o option.o \
                             pflotran_constants.o reaction_aux.o region.o \
                             srcsink_sandbox_base.o string.o units.o
srcsink_sandbox_wipp_well.o : input_aux.o material_aux.o option.o \
                              pflotran_constants.o reaction_aux.o region.o \
                              srcsink_sandbox_base.o string.o units.o
strata.o : input_aux.o material.o option.o pflotran_constants.o region.o \
           string.o surface_material.o units.o
string.o : pflotran_constants.o
surface_auxiliary.o : pflotran_constants.o surface_global_aux.o \
                      surface_th_aux.o
surface_field.o : pflotran_constants.o
surface_flow.o : connection.o coupler.o debug.o discretization.o dm_kludge.o \
                 eos_water.o global_aux.o grid.o logging.o option.o \
                 output_aux.o patch.o pflotran_constants.o realization.o \
                 realization_base.o string.o surface_field.o \
                 surface_global_aux.o surface_material.o surface_realization.o \
                 variables.o
surface_global.o : connection.o coupler.o discretization.o grid.o option.o \
                   patch.o pflotran_constants.o surface_field.o \
                   surface_global_aux.o surface_realization.o variables.o
surface_global_aux.o : option.o pflotran_constants.o utility.o
surface_init.o : condition.o coupler.o dataset.o dataset_base.o \
                 dataset_common_hdf5.o debug.o discretization.o grid.o \
                 grid_structured.o grid_unstructured.o grid_unstructured_aux.o \
                 input_aux.o option.o output_aux.o output_surface.o \
                 output_tecplot.o patch.o pflotran_constants.o region.o \
                 solver.o strata.o string.o surface_auxiliary.o \
                 surface_material.o surface_realization.o units.o waypoint.o
surface_material.o : input_aux.o option.o pflotran_constants.o string.o
surface_realization.o : condition.o coupler.o dataset.o dataset_base.o debug.o \
                        discretization.o dm_kludge.o field.o grid.o \
                        grid_unstructured.o grid_unstructured_aux.o \
                        grid_unstructured_cell.o input_aux.o option.o \
                        output_aux.o patch.o pflotran_constants.o \
                        reaction_aux.o realization.o realization_base.o \
                        region.o strata.o string.o surface_field.o \
                        surface_material.o time_storage.o waypoint.o
surface_th.o : connection.o coupler.o debug.o discretization.o dm_kludge.o \
               eos_water.o fluid.o grid.o option.o output_aux.o patch.o \
               pflotran_constants.o realization.o realization_base.o region.o \
               string.o surface_field.o surface_global_aux.o \
               surface_material.o surface_realization.o surface_th_aux.o \
               variables.o
surface_th_aux.o : eos_water.o option.o pflotran_constants.o \
                   surface_global_aux.o
th.o : condition.o connection.o coupler.o debug.o discretization.o eos_water.o \
       field.o fluid.o global_aux.o grid.o logging.o material.o material_aux.o \
       option.o output_aux.o patch.o pflotran_constants.o realization.o \
       region.o saturation_function.o secondary_continuum.o \
       secondary_continuum_aux.o string.o th_aux.o utility.o variables.o
th_aux.o : eos_water.o global_aux.o material_aux.o option.o \
           pflotran_constants.o saturation_function.o
time_storage.o : option.o pflotran_constants.o
timestepper_BE.o : convergence.o input_aux.o option.o output.o \
                   pflotran_constants.o pm_base.o solver.o string.o \
                   timestepper_base.o utility.o
timestepper_base.o : input_aux.o option.o output.o pflotran_constants.o \
                     pm_base.o string.o waypoint.o
timestepper_geomechanics.o : convergence.o option.o output.o \
                             pflotran_constants.o pm_base.o \
                             pm_geomechanics_force.o solver.o \
                             timestepper_base.o waypoint.o
timestepper_surface.o : option.o output.o pflotran_constants.o pm_base.o \
                        pm_surface_flow.o solver.o surface_flow.o \
                        timestepper_base.o waypoint.o
transport.o : connection.o global_aux.o material_aux.o matrix_block_aux.o \
              option.o pflotran_constants.o reactive_transport_aux.o
transport_constraint.o : global_aux.o input_aux.o logging.o option.o \
                         pflotran_constants.o reaction_aux.o \
                         reaction_immobile_aux.o reaction_mineral_aux.o \
                         reaction_surf_complex_aux.o reactive_transport_aux.o \
                         string.o units.o
uniform_velocity.o : input_aux.o logging.o option.o pflotran_constants.o \
                     string.o units.o utility.o
units.o : option.o pflotran_constants.o
utility.o : input_aux.o option.o pflotran_constants.o string.o
variables.o :
waypoint.o : option.o output_aux.o pflotran_constants.o utility.o
wipp.o : input_aux.o lookup_table.o option.o pflotran_constants.o string.o \
         units.o utility.o
wrapper_hydrogeophysics.o : e4d_run.o e4d_setup.o e4d_vars.o option.o \
                            pflotran_constants.o

# developer level regression testing
test : utest rtest

utest : pflotran
	@if [ -d $(unit_test_dir) ]; then \
		$(MAKE) --directory=$(unit_test_dir) test; \
	fi

rtest : pflotran
	@if [ -d $(regression_test_dir) ]; then \
		$(MAKE) --directory=$(regression_test_dir) test; \
	fi

clean-tests : FORCE
	@if [ -d $(unit_test_dir) ]; then \
		$(MAKE) --directory=$(unit_test_dir) $(MAKECMDGOALS); \
	fi
	@if [ -d $(regression_test_dir) ]; then \
		$(MAKE) --directory=$(regression_test_dir) $(MAKECMDGOALS); \
	fi

# user level tests to verify pflotran is built correctly
check : FORCE
	@if [ -d $(regression_test_dir) ]; then \
		$(MAKE) --directory=$(regression_test_dir) $(MAKECMDGOALS); \
	fi

# null rule to force things to happen
FORCE :

