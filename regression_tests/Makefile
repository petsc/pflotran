# Makefile for running standard regression tests
#
# include the make variables from PETSc so we can use PYTHON and do
# some conditional testing, e.g. test if unstructured mesh is available.
#
# If PYTHON defined from petsc is not correct, override it on the
# command line with: make PYTHON=python3.3 test

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


TEST_MANAGER = regression_tests.py
PFLOTRAN = ../src/pflotran/pflotran

TEST_OPTIONS =

# make PERFORMANCE=true check
ifdef PERFORMANCE
	TEST_OPTIONS += --check-performance
endif

ifdef BACKTRACE
	TEST_OPTIONS += --backtrace
endif

ifdef TIMEOUT
	TEST_OPTIONS += --timeout $(TIMEOUT)
endif

ifneq ($(strip $(MPIEXEC)),)
	TEST_OPTIONS += --mpiexec $(MPIEXEC)
endif

#
# standard tests that are run to verify pflotran is built correctly
#
STANDARD_CFG = \
	ascem/batch/batch.cfg \
	ascem/1d/1d-calcite/1d-calcite.cfg \
	default/input_format/input_format.cfg \
	default/543/543.cfg \
	default/anisothermal/anisothermal.cfg \
	default/column/column.cfg \
	default/batch/batch.cfg \
	default/eos/eos.cfg \
	default/dispersion/dispersion.cfg \
	default/materials/materials.cfg \
	default/materials/soil_compressibility/soil_compressibility.cfg \
	default/materials/mineral_volume/mineral_volume.cfg \
	default/infiltrometer/infiltrometer.cfg \
	default/numerical/numerical.cfg \
	default/output/output.cfg \
	default/restart/restart.cfg \
	default/steady/steady.cfg \
	default/condition/condition.cfg \
	default/condition/dataset/dataset.cfg \
	default/multicontinuum/multicontinuum.cfg \
	default/scco2/mphase/mphase.cfg \
	default/scco2/flash2/flash2.cfg \
	default/regions/regions.cfg \
	default/regions/cell_mapping/cell_mapping.cfg \
	default/salinity/salinity.cfg \
	default/timestepping/timestepping.cfg \
	geothermal_hpt/1D_Calcite/calcite_hpt.cfg\
	general/general.cfg \
	hydrogeophysics/hydrogeophysics.cfg \
	ngee/ngee.cfg \
	shortcourse/1D_Calcite/calcite.cfg \
	shortcourse/1D_variably_saturated_flow/vsat_flow.cfg \
	shortcourse/regional_doublet/regional_doublet.cfg \
	shortcourse/copper_leaching/cu_leaching.cfg \
	toil_ims/toil_ims.cfg \
	towg/towg.cfg \
	ufd/ufd.cfg \
	wipp/nuts/nuts.cfg \
	wipp/bragflo/bragflo.cfg

ifneq ($(strip $(PARMETIS_LIB)),)
	STANDARD_CFG += \
		default/discretization/discretization.cfg \
		geomech/uniaxial_strain/uniaxial_strain.cfg \
	        geomech/uniaxial_loading/uniaxial_loading.cfg \
        	geomech/uniaxial_loading_biot/uniaxial_loading_biot.cfg \
	        surfaceflow/surface_flow_only/surface_only.cfg \
		surfaceflow/surface_th/surface_th.cfg \
		surfaceflow/inlinesurface/inlinesurface.cfg
endif

#
# domain specific problems
#
GEOCHEMISTRY_CFG = \
	ascem/1d/1d-calcite/1d-calcite.cfg \
	ascem/batch/batch.cfg \
	default/543/543.cfg \
	default/anisothermal/anisothermal.cfg \
	default/column/column.cfg \
	default/batch/batch.cfg \
	default/multicontinuum/multicontinuum.cfg \
	geothermal_hpt/1D_Calcite/calcite_hpt.cfg\
        ngee/ngee.cfg \
	shortcourse/1D_Calcite/calcite.cfg \
	shortcourse/copper_leaching/cu_leaching.cfg \
	ufd/ufd.cfg \
	wipp/nuts/nuts.cfg

GEOMECHANICS_CFG = \
	geomech/uniaxial_strain/uniaxial_strain.cfg \
	geomech/uniaxial_loading/uniaxial_loading.cfg \
	geomech/uniaxial_loading_biot/uniaxial_loading_biot.cfg

FLOW_CFG = \
	default/543/543.cfg \
	default/condition/condition.cfg \
	default/eos/eos.cfg \
	default/infiltrometer/infiltrometer.cfg \
	default/steady/steady.cfg \
	default/scco2/mphase/mphase.cfg \
	default/scco2/flash2/flash2.cfg \
	general/general.cfg \
        toil_ims/toil_ims.cfg \
	towg/towg.cfg \
	shortcourse/1D_variably_saturated_flow/vsat_flow.cfg \
	shortcourse/copper_leaching/cu_leaching.cfg \
	wipp/bragflo/bragflo.cfg

TRANSPORT_CFG = \
	default/543/543.cfg \
	default/column/column.cfg \
	default/dispersion/dispersion.cfg \
	default/multicontinuum/multicontinuum.cfg \

MESH_CFG = \
	default/discretization/discretization.cfg

MFD_CFG = \
	mfd/mfd.cfg

CHECK_CFG = \
	ascem/1d/1d-calcite/1d-calcite.cfg

NOHDF5_CFG = \
	default/column/column.cfg \
	default/scco2/flash2/flash2.cfg \
	geothermal_hpt/1D_Calcite/calcite_hpt.cfg \
	default/materials/soil_compressibility/soil_compressibility.cfg \
	wipp/nuts/nuts.cfg \
	default/multicontinuum/multicontinuum.cfg \
	wipp/bragflo/bragflo.cfg \
	ufd/ufd.cfg \
	hydrogeophysics/hydrogeophysics.cfg \
	ascem/1d/1d-calcite/1d-calcite.cfg \
	default/materials/mineral_volume/mineral_volume.cfg \
	default/batch/batch.cfg \
	shortcourse/1D_Calcite/calcite.cfg \
	ascem/batch/batch.cfg \
	default/scco2/mphase/mphase.cfg \
	default/input_format/input_format.cfg \
	shortcourse/copper_leaching/cu_leaching.cfg \
	default/eos/eos.cfg \
	default/regions/regions.cfg \
	shortcourse/1D_variably_saturated_flow/vsat_flow.cfg

test : standard

standard :
ifeq ($(strip $(PARMETIS_LIB)),)
	@echo "********************************************************"
	@echo "  PFLOTRAN does not appear to be compiled with Parmetis."
	@echo "  Skipping unstructured mesh tests."
	@echo "********************************************************"
endif
	$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard standard_parallel \
		--config-files $(STANDARD_CFG)


standard_parallel :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard_parallel \
		--config-files $(STANDARD_CFG)

check : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard standard_parallel \
		--config-files $(CHECK_CFG)

geochemistry :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(GEOCHEMISTRY_CFG) --suite geochemistry 

geomechanics :
ifneq ($(strip $(PARMETIS_LIB)),)
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(GEOMECHANICS_CFG) --suite geomechanics 
else
	@echo "********************************************************"
	@echo "  PFLOTRAN does not appear to be compiled with Parmetis."
	@echo "  Skipping geomechanics tests."
	@echo "********************************************************"
endif

flow :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(FLOW_CFG) --suite flow

transport :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(TRANSPORT_CFG) --suite transport

mesh :
ifneq ($(strip $(PARMETIS_LIB)),)
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(MESH_CFG) --suite mesh
else
	@echo "********************************************************"
	@echo "  PFLOTRAN does not appear to be compiled with Parmetis."
	@echo "  Skipping unstructured mesh tests."
	@echo "********************************************************"
endif

ngee-biogeochemistry : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		-r ngee --suite biogeochemistry 

test_mimetic :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(MFD_CFG)

nonhdf5_test :
	$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard standard_parallel \
		--config-files $(NOHDF5_CFG)

#
# Try to cleanup the regression test output.
#
# NOTE: do NOT add a "*.h5" glob. The .h5 extension is used for both
# input and output files, and a simple glob can not distinguish
# between them. Manually add rm commands for each file that needs to be
# removed.
#
clean-tests :
	-find . -type f -name '*.testlog' -print0 | xargs -0 rm
	-find . -type f -name '*.out' -print0 | xargs -0 rm
	-find . -type f -name '*.tec' -print0 | xargs -0 rm
	-find . -type f -name '*.regression' -print0 | xargs -0 rm
	-find . -type f -name '*.stdout' -print0 | xargs -0 rm
	-find . -type f -name '*.old' -print0 | xargs -0 rm
	-find . -type f -name '*~' -print0 | xargs -0 rm
	-find . -type f -name 'tmp-restart-*' -print0 | xargs -0 rm
	-find . -type f -name '*.chk' -print0 | xargs -0 rm
	-find . -type f -name '*.pyc' -print0 | xargs -0 rm
	-find . -type f -name 'sigma*.txt' -print0 | xargs -0 rm 
	-find . -type f -name '*.dpd' -print0 | xargs -0 rm
	-find . -type f -name 'e4d*.log' -print0 | xargs -0 rm
	-find . -type d -name '__pycache__' -print0 | xargs -0 rm -r
	-rm -f default/anisothermal/thc_1d.h5
	-rm -f default/anisothermal/th_1d.h5
	-rm -f mfd/mfd-mas.dat
	-rm -f mfd/mfd.h5
	-rm -f default/543/543_flow_and_tracer-np8.h5
	-rm -f default/543/543_flow_and_tracer.h5
	-rm -f default/discretization/mixed_implicit-np4.h5
	-rm -f default/discretization/mixed_implicit.h5
	-rm -f general/calcite_flow_and_tran_general.h5
	-rm -f ngee/sink-mas.dat 
	-rm -f ngee/sink.h5 
	-rm -f ngee/source-mas.dat
	-rm -f ngee/source.h5
	-rm -f default/543/543_flow_and_tracer_dbase.h5
	-rm -f wipp/bragflo/gas_generation_2d.h5
	-rm -f wipp/bragflo/gas_generation_one_cell.h5
	-rm -f ascem/batch/calcite-kinetics-int.dat
	-rm -f ascem/batch/calcite-kinetics-volume-fractions-int.dat
	-rm -f ascem/batch/general-reaction-int.dat
	-rm -f default/multicontinuum/0D_heat_MPH-int.dat
	-rm -f default/multicontinuum/0D_heat_TH-int.dat
	-rm -f ngee/CLM-CN-int.dat
	-rm -f shortcourse/1D_variably_saturated_flow/vsat_flow_pulse_2layer-int.dat
	-rm -f shortcourse/copper_leaching/cu_leaching-int.dat
	-rm -f wipp/nuts/case3-int.dat
	-rm -f wipp/nuts/case4r1_modified-int.dat
	-rm -f wipp/bragflo/gas_generation_2d-np4.h5
	-rm -f default/543/543_hanford_srfcplx_base*.h5
	-rm -f default/condition/dataset/east_surface_cell_centered.h5
	-rm -f default/condition/dataset/east_surface_node_centered.h5
	-rm -f default/condition/dataset/north_surface_cell_centered.h5
	-rm -f default/condition/dataset/north_surface_node_centered.h5
	-rm -f default/condition/dataset/top_surface_cell_centered.h5
	-rm -f default/condition/dataset/top_surface_node_centered.h5
	-rm -f default/condition/dataset/x_line_cell_centered.h5
	-rm -f default/condition/dataset/x_line_node_centered.h5
	-rm -f default/condition/dataset/y_line_cell_centered.h5
	-rm -f default/condition/dataset/y_line_node_centered.h5
	-rm -f default/condition/dataset/z_line_cell_centered.h5
	-rm -f default/condition/dataset/z_line_node_centered.h5
	-rm -f default/salinity/20x20_xz.h5
	-rm -f default/timestepping/20x20_xz.h5
	-rm -f hydrogeophysics/e4d__000000000003600.dpd
	-rm -f hydrogeophysics/e4d__000000000010800.dpd
	-rm -f hydrogeophysics/e4d__000000000025200.dpd
	-rm -f hydrogeophysics/e4d__000000000054000.dpd
	-rm -f hydrogeophysics/e4d__000000000111600.dpd
	-rm -f hydrogeophysics/sigma__000000000003600.txt
	-rm -f hydrogeophysics/sigma__000000000010800.txt
	-rm -f hydrogeophysics/sigma__000000000025200.txt
	-rm -f hydrogeophysics/sigma__000000000054000.txt
	-rm -f hydrogeophysics/sigma__000000000111600.txt
	-rm -f ufd/glass-wf_mass-0.dat
	-rm -f ufd/fmdm_dummy-wf_mass-0.dat
	-rm -f ufd/wf_doubles-wf_mass-0.dat
	-rm -f ufd/wf_region-wf_mass-*.dat
	-rm -f wipp/bragflo/well_production2-ss_mass-0.dat
	-rm -f ufd/glass-*.vtk
	-rm -f ufd/wf_doubles-*.vtk
	-rm -f ufd/wf_region-*.vtk
	-rm -f ufd/glass_full_geochem_np4-*.vtk

