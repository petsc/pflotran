# Makefile for running standard regression tests
#
# include the make variables from PETSc so we can use PYTHON and do
# some conditional testing, e.g. test if unstructured mesh is available.
#
# If PYTHON definied from petsc is not correct, override it on the
# command line with: make PYTHON=python3.3 test

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


TEST_MANAGER = regression-tests.py
PFLOTRAN = ../src/pflotran/pflotran

TEST_OPTIONS =

# make VERBOSE=true check
ifdef VERBOSE
	TEST_OPTIONS += --verbose
endif

# make PERFORMANCE=true check
ifdef PERFORMANCE
	TEST_OPTIONS += --check-performance
endif

ifdef BACKTRACE
	TEST_OPTIONS += --backtrace
endif

#
# standard tests that are run to verify pflotran is built correctly
#
STANDARD_CFG = \
	ascem/batch/batch.cfg \
	ascem/1d/1d-calcite/1d-calcite.cfg \
	default/543/543.cfg \
	default/anisothermal/anisothermal.cfg \
	default/column/column.cfg \
	default/infiltrometer/infiltrometer.cfg \
	default/condition/condition.cfg \
	default/multicontinuum/multicontinuum.cfg \
	ngee/ngee.cfg \
	shortcourse/copper_leaching/cu_leaching.cfg

STANDARD_PARALLEL_CFG = \
	ascem/1d/1d-calcite/1d-calcite.cfg \
	default/543/543.cfg

ifneq ($(strip $(PARMETIS_LIB)),)
	STANDARD_CFG += \
		default/discretization/discretization.cfg

	STANDARD_PARALLEL_CFG = \
		default/discretization/discretization.cfg
else
	@echo "********************************************************"
	@echo "  PFloTran does not appear to be compiled with Parmetis."
	@echo "  Skipping unstructured mesh tests."
	@echo "********************************************************"
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
	default/multicontinuum/multicontinuum.cfg \
	ngee/ngee.cfg \
	shortcourse/1D_Calcite/calcite.cfg \
	shortcourse/copper_leaching/cu_leaching.cfg

FLOW_CFG = \
	default/543/543.cfg \
	default/condition/condition.cfg \
	default/infiltrometer/infiltrometer.cfg \
	shortcourse/1D_variably_saturated_flow/vsat_flow.cfg \
	shortcourse/copper_leaching/cu_leaching.cfg

TRANSPORT_CFG = \
	default/543/543.cfg \
	default/column/column.cfg \
	default/multicontinuum/multicontinuum.cfg \

MESH_CFG = \
	default/discretization/discretization.cfg


check : standard standard_parallel

test : standard standard_parallel

standard :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		 --suite standard --config-files $(STANDARD_CFG)

standard_parallel :
ifneq ($(strip $(MPIEXEC)),)
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--mpiexec $(MPIEXEC)  --suite standard_parallel \
	--config-files $(STANDARD_PARALLEL_CFG)
else
	@echo "********************************************************"
	@echo "  Could not find mpiexec."
	@echo "  Skipping parallel tests."
	@echo "********************************************************"
endif

geochemistry :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(GEOCHEMISTRY_CFG) --suite geochemistry 

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
	@echo "  PFloTran does not appear to be compiled with Parmetis."
	@echo "  Skipping unstructured mesh tests."
	@echo "********************************************************"
endif

ngee-biogeochemistry : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		-r ngee --suite biogeochemistry 

