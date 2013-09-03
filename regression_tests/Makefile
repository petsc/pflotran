# Makefile for running standard regression tests
#
# include the make variables from PETSc so we can use PYTHON and do
# some conditional testing, e.g. test if unstructured mesh is available.
#
# If PYTHON defined from petsc is not correct, override it on the
# command line with: make PYTHON=python3.3 test

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


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

ifneq ($(strip $(MPIEXEC)),)
	TEST_OPTIONS += --mpiexec $(MPIEXEC)
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
	default/mphase/mphase.cfg \
	ngee/ngee.cfg \
	shortcourse/1D_Calcite/calcite.cfg \
	shortcourse/1D_variably_saturated_flow/vsat_flow.cfg \
	shortcourse/regional_doublet/regional_doublet.cfg \
	shortcourse/copper_leaching/cu_leaching.cfg 

ifneq ($(strip $(PARMETIS_LIB)),)
	STANDARD_CFG += \
		default/discretization/discretization.cfg

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

MFD_CFG = \
	mfd/mfd.cfg

SCCO2_CFG = \
	scco2/scco2.cfg

CHECK_CFG = \
	ascem/1d/1d-calcite/1d-calcite.cfg

test : standard

standard :
ifeq ($(strip $(PARMETIS_LIB)),)
	@echo "********************************************************"
	@echo "  PFloTran does not appear to be compiled with Parmetis."
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

test_mimetic :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(MFD_CFG)

test_scco2 :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(SCCO2_CFG)

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
	-rm -f default/anisothermal/thc_1d.h5
	-rm -f default/anisothermal/th_1d.h5
	-rm -f mfd/mfd-mas.dat
	-rm -f mfd/mfd.h5
	-rm -f scco2/scco2.h5
	-rm -f scco2/scco2-mas.dat
