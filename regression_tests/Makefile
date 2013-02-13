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

# tests that are run to verify pflotran is built correctly
STANDARD_TESTS = \
	ascem-geochemistry \
	ascem-geochemistry-parallel \
	default-flow \
	default-transport \
	default-geochemistry \
	default-parallel \
	default-discretization-umesh \
	ngee-biogeochemistry

# regression tests for developers changing pflotran
DEV_TESTS = \
	cu-leaching

check : $(STANDARD_TESTS)

test : $(STANDARD_TESTS) $(DEV_TESTS)

ascem-geochemistry : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		-r ascem --suite geochemistry 

ascem-geochemistry-parallel :
ifneq ($(strip $(MPIEXEC)),)
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--mpiexec $(MPIEXEC) -r ascem/1d --suite parallel
else
	@echo "********************************************************"
	@echo "  Could not find mpiexec."
	@echo "  Skipping parallel tests."
	@echo "********************************************************"
endif

default-discretization-umesh :
ifneq ($(strip $(PARMETIS_LIB)),)
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		-r default/discretization --suite mesh
else
	@echo "********************************************************"
	@echo "  PFloTran does not appear to be compiled with Parmetis."
	@echo "  Skipping unstructured mesh tests."
	@echo "********************************************************"
endif

default-flow : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		-r default --suite flow 

default-transport : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		-r default --suite transport 

default-geochemistry : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		-r default --suite geochemistry 

default-parallel :
ifneq ($(strip $(MPIEXEC)),)
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--mpiexec $(MPIEXEC) -r default --suite parallel
else
	@echo "********************************************************"
	@echo "  Could not find mpiexec."
	@echo "  Skipping parallel tests."
	@echo "********************************************************"
endif

cu-leaching :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		-c shortcourse/copper_leaching/cu_leaching.cfg --test cu_leaching

ngee-biogeochemistry : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		-r ngee --suite biogeochemistry 

