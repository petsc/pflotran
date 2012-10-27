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

# make VERBOSE=true check
ifdef VERBOSE
	VERBOSE_FLAG=--verbose
else
	VERBOSE_FLAG=
endif

# tests that are run to verify pflotran is built correctly
STANDARD_TESTS = \
	ascem-geochemistry \
	ascem-geochemistry-parallel \
	default-discretization-umesh

# regression tests for developers changing pflotran
DEV_TESTS = \
	cu-leaching

check : $(STANDARD_TESTS)

test : $(STANDARD_TESTS) $(DEV_TESTS)

ascem-geochemistry : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(VERBOSE_FLAG) \
		-r ascem --suite geochemistry 

ascem-geochemistry-parallel :
ifneq ($(strip $(MPIEXEC)),)
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(VERBOSE_FLAG) \
		--mpiexec $(MPIEXEC) -r ascem/1d --suite parallel
else
	@echo "Could not find mpiexec."
	@echo "Skipping parallel tests."
endif

default-discretization-umesh :
ifneq ($(strip $(PARMETIS_LIB)),)
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(VERBOSE_FLAG) \
		-r default/discretization --suite mixed
else
	@echo "PFloTran does not appear to be compiled with Parmetis."
	@echo "Skipping unstructured mesh tests."
endif


cu-leaching :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(VERBOSE_FLAG) \
		-c shortcourse/copper_leaching/cu_leaching.cfg --test cu_leaching
