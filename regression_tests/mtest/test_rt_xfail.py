#!/usr/bin/env python
"""Expected fail functionality tests for the pflotran regression test
manager

"""

from __future__ import print_function

import os
import re
import subprocess
import unittest

import regression_tests


class ExpectedFailTests_IndividualMain(unittest.TestCase):
    """tests: higher level functionality for the complete test
    manager (stand alone functions, and the RegressionTestManager and
    RegressionTest classes)

    Run the expected fail suite and verify that all tests fail for the
    correct reason.

    Notes:

    - This test class works by creating a dummy commandline options
      object and then calling regression_tests.main() for each
      expected fail test.

    - This method is prefered over using a subprocess to call the test
      manager, because it allows us to generate more meaningful code
      coverage reports.

    - Running each test separately (instead of running the entire test
      suite at once) allows us to search each log file individually
      for the correct fail rather than just verifing that it is in the
      file someplace.

    - 2013-06-26: some tests truncate the search string (full_string
      and short_string) because python2.7 and python3.3 have different
      roundoff...

    """

    class Options(object):
        def __init__(self):
            self.advanced = False
            self.backtrace = False
            self.check_only = False
            self.check_performance = False
            self.config_files = ['mtest/expected-fail.cfg']
            self.debug = False
            self.dry_run = False
            self.executable = ['../src/pflotran/pflotran']
            self.list_suites = False
            self.list_tests = False
            self.mpiexec = None
            self.new_tests = False
            self.recursive_search = None
            self.suites = []
            self.tests = []
            self.timeout = None
            self.update = False

    def setUp(self):
        self.test_stdout = open("xfail-stdout.testlog", 'w')
        # run make clean-tests just incase it hasn't been done
        self.clean_cmd = ["make", "clean-tests"]
        subprocess.call(self.clean_cmd, shell=False,
                        stdout=self.test_stdout, stderr=subprocess.STDOUT)

        self.cmdl_options = self.Options()

        # construct the basic commanad to run tests
        self.run_cmd = ["/usr/bin/env",
                        "python",
                        "regression_tests.py",
                        "--executable",
                        "--config-files",
                        "mtest/expected-fail.cfg",
                        "--tests"]

    def tearDown(self):
        pass
        #subprocess.call(self.clean_cmd, shell=False, stdout=self.test_stdout)
        #self.test_stdout.close()

    def find_testlog(self):
        testlog_re = re.compile(r"^pflotran-tests-[\w_-]+\.testlog")
        testlog = None
        for entry in os.listdir(os.getcwd()):
            match = testlog_re.search(entry)
            if match is not None:
                testlog = entry
                break
        return testlog

    def search_testlog(self, testlog_name, failure_string):
        logdata = None
        with open(testlog_name, 'r') as testlog:
            logdata = testlog.readlines()
        failure_found = False
        for line in logdata:
            if line.find(failure_string) != -1:
                failure_found = True
                break
        return failure_found

    def check_testlog(self, testlog, failure_string):
        failure_found = self.search_testlog(testlog, failure_string)
        self.assertEqual(failure_found, True, failure_string)

    def run_regression_test_manager(self):
        try:
            status = regression_tests.main(self.cmdl_options)
        except:
            pass
        testlog = self.find_testlog()
        self.assertNotEqual(testlog, None)
        return testlog

    def test_xfail_discrete_mean(self):
        self.cmdl_options.tests.append("fail-discrete-mean")
        testlog = self.run_regression_test_manager()
        self.check_testlog(testlog,
                           "FAIL: Material ID:Mean : 0.5 > 1e-12 [absolute]")
        

    def test_xfail_discrete(self):
        self.cmdl_options.tests.append("fail-discrete")
        testlog = self.run_regression_test_manager()
        self.check_testlog(testlog, "FAIL: Material ID:Max : 1 > 0 [absolute]")
        self.check_testlog(testlog, "FAIL: Material ID:1 : 1 > 0 [absolute]")

    def test_xfail_volume_fraction(self):
        self.cmdl_options.tests.append("fail-volume-fraction")
        testlog = self.run_regression_test_manager()
        full_string = "FAIL: Calcite VF:Min : 0.197999999817 > 1e-12 [absolute]"
        short_string = "FAIL: Calcite VF:Min : 0.197999999817"
        self.check_testlog(testlog, short_string)
        full_string = "FAIL: Calcite VF:1 : 0.179999999834 > 1e-12 [absolute]"
        short_string = "FAIL: Calcite VF:1 : 0.179999999834"
        self.check_testlog(testlog, short_string)

    def test_xfail_skipped_solution_s2norm(self):
        self.cmdl_options.tests.append("fail-solution-s2norm")
        testlog = self.run_regression_test_manager()
        self.check_testlog(testlog, "Skipping SOLUTION : Transport")

    def test_xfail_solution_s2norm(self):
        self.cmdl_options.tests.append("fail-solution-s2norm")
        self.cmdl_options.check_performance = True
        testlog = self.run_regression_test_manager()
        full_string = "FAIL: Transport:Solution 2-Norm : 0.00137252117371 > 1e-12 [absolute]"
        short_string = "FAIL: Transport:Solution 2-Norm : 0.0013725211"
        self.check_testlog(testlog, short_string)

    def test_xfail_solution_linear(self):
        self.cmdl_options.tests.append("fail-solution-linear")
        self.cmdl_options.check_performance = True
        testlog = self.run_regression_test_manager()
        self.check_testlog(
            testlog, "FAIL: Transport:Solver Iterations : 1800 > 0 [absolute]")

    def test_xfail_vector_length(self):
        self.cmdl_options.tests.append("fail-vector-length")
        testlog = self.run_regression_test_manager()
        self.check_testlog(
            testlog,
            "FAIL: Gamma HCO3- : 1 : vector lengths not equal. gold 2, current 1")

    def test_xfail_solution_steps(self):
        self.cmdl_options.tests.append("fail-solution-steps")
        self.cmdl_options.check_performance = True
        testlog = self.run_regression_test_manager()
        self.check_testlog(testlog,
                           "FAIL: Transport:Time Steps : 90 > 0.0 [absolute]")

    def test_xfail_section_missing(self):
        self.cmdl_options.tests.append("fail-section-missing")
        testlog = self.run_regression_test_manager()
        self.check_testlog(
            testlog,
            "FAIL: section 'Total Ca++' is in the gold output, but not the current output.")

    def test_xfail_timeout(self):
        self.cmdl_options.tests.append("fail-timeout")
        testlog = self.run_regression_test_manager()
        self.check_testlog(testlog,
                           "FAIL : fail-timeout : pflotran return an error code")

    def test_xfail_solution_r2norm(self):
        self.cmdl_options.tests.append("fail-solution-r2norm")
        self.cmdl_options.check_performance = True
        testlog = self.run_regression_test_manager()
        full_string = "FAIL: Transport:Residual 2-Norm : 1.45584930368e+18 > 1e-12 [absolute]"
        short_string = "FAIL: Transport:Residual 2-Norm : 1.45584930368"
        self.check_testlog(testlog, short_string)

    def test_xfail_generic(self):
        self.cmdl_options.tests.append("fail-generic")
        testlog = self.run_regression_test_manager()
        full_string = "FAIL: pH:Min : 3.67633734668 > 1e-12 [absolute]"
        short_string = "FAIL: pH:Min : 3.676337346"
        self.check_testlog(testlog, short_string)
        self.check_testlog(testlog, "FAIL: pH:1 : 1.0 > 1e-12 [absolute]")

    def test_xfail_pressure(self):
        self.cmdl_options.tests.append("fail-pressure")
        testlog = self.run_regression_test_manager()
        full_string = "FAIL: Liquid Pressure:61 : 55133.6565938 > 1e-05 [absolute]"
        short_string = "FAIL: Liquid Pressure:61 : 55133.6565938"
        self.check_testlog(testlog, short_string)

    def test_xfail_rate(self):
        self.cmdl_options.tests.append("fail-rate")
        testlog = self.run_regression_test_manager()
        self.check_testlog(testlog,
                           "FAIL: Calcite Rate:Max : 2.0 > 1e-07 [relative]")
        self.check_testlog(testlog,
                           "FAIL: Calcite Rate:1 : 2.0 > 1e-07 [relative]")

    def test_xfail_section_extra(self):
        self.cmdl_options.tests.append("fail-section-extra")
        testlog = self.run_regression_test_manager()
        self.check_testlog(
            testlog,
            "FAIL: section 'Gamma H+' is in the current output, but not the gold output.")

    def test_xfail_velocity(self):
        self.cmdl_options.tests.append("fail-velocity")
        testlog = self.run_regression_test_manager()
        self.check_testlog(
            testlog,
            "FAIL: LIQUID VELOCITY [m/y]:81 : 1.0 > 1e-12 [absolute]")

    def test_xfail_solution_cuts(self):
        self.cmdl_options.tests.append("fail-solution-cuts")
        self.cmdl_options.check_performance = True
        testlog = self.run_regression_test_manager()
        self.check_testlog(testlog,
                           "FAIL: Transport:Time Step Cuts : 2 > 0 [absolute]")

    def test_xfail_solution_time(self):
        self.cmdl_options.tests.append("fail-solution-time")
        self.cmdl_options.check_performance = True
        testlog = self.run_regression_test_manager()
        full_string = "FAIL: Transport:Time (seconds) : 99.9999999995 > 1.0 [percent]"
        short_string = "FAIL: Transport:Time (seconds) : 99.999999"
        self.check_testlog(testlog, short_string)

    def test_xfail_concentration(self):
        self.cmdl_options.tests.append("fail-concentration")
        testlog = self.run_regression_test_manager()
        full_string = "FAIL: Total HCO3-:1 : 0.0001 > 1e-12 [absolute]"
        short_string = "FAIL: Total HCO3-:1 : 0.0001"
        self.check_testlog(testlog, short_string)
        full_string = "FAIL: Total HCO3-:Mean : 0.0001 > 1e-12 [absolute]"
        short_string = "FAIL: Total HCO3-:Mean : 0.0001"
        self.check_testlog(testlog, short_string)

    def test_xfail_solution_newton(self):
        self.cmdl_options.tests.append("fail-solution-newton")
        self.cmdl_options.check_performance = True
        testlog = self.run_regression_test_manager()
        self.check_testlog(
            testlog, "FAIL: Transport:Newton Iterations : 180 > 0 [absolute]")

    def test_xfail_saturation(self):
        self.cmdl_options.tests.append("fail-saturation")
        testlog = self.run_regression_test_manager()
        self.check_testlog(
            testlog, "FAIL: Liquid Saturation:81 : 0.1 > 1e-12 [absolute]")

    if False:
        def test_xfail_(self):
            self.cmdl_options.tests.append("fail-")
            testlog = self.run_regression_test_manager()
            self.check_testlog(testlog, "???")


@unittest.skip("Skipping expected fails all in one.")
class ExpectedFailTests_SubprocessAllInOne(unittest.TestCase):
    """tests: runtime behavior of regression_tests.py

    Run the expected fail suite and verify that all tests fail.

    Notes:

    - This runs the entire suite with a single subprocess call to the
      test manager

    - This is does NOT verify that each xfail test fails for the
      correct reason, just that all the expected fail messages are in
      the test log.

    """

    def setUp(self):
        self.test_stdout = open("xfail-stdout.testlog", 'w')
        # run make clean-tests just incase it hasn't been done
        self.clean_cmd = ["make", "clean-tests"]
        subprocess.call(self.clean_cmd, shell=False,
                        stdout=self.test_stdout, stderr=subprocess.STDOUT)

        # construct the basic commanad to run tests
        self.run_cmd = ["/usr/bin/env",
                        "python",
                        "regression_tests.py",
                        "--executable",
                        "../src/pflotran/pflotran",
                        "--config-files",
                        "mtest/expected-fail.cfg"]

        self.failure_strings_general = [
            "FAIL: Material ID:Mean : 0.5 > 1e-12 [absolute]",
            "FAIL: Material ID:Max : 1 > 0 [absolute]",
            "FAIL: Material ID:1 : 1 > 0 [absolute]",
            "FAIL: Calcite VF:Min : 0.197999999817 > 1e-12 [absolute]",
            "FAIL: Calcite VF:1 : 0.179999999834 > 1e-12 [absolute]",
            "FAIL: Gamma HCO3- : 1 : vector lengths not equal. gold 2, current 1",
            "FAIL: section 'Total Ca++' is in the gold output, but not the current output.",
            "FAIL : fail-timeout : pflotran return an error code",
            "FAIL: pH:Min : 3.67633734668 > 1e-12 [absolute]",
            "FAIL: pH:1 : 1.0 > 1e-12 [absolute]",
            "FAIL: Liquid Pressure:61 : 55133.6565938 > 1e-05 [absolute]",
            "FAIL: Calcite Rate:Max : 2.0 > 1e-07 [relative]",
            "FAIL: Calcite Rate:1 : 2.0 > 1e-07 [relative]",
            "FAIL: section 'Gamma H+' is in the current output, but not the gold output.",
            "FAIL: LIQUID VELOCITY [m/y]:81 : 1.0 > 1e-12 [absolute]",
            "FAIL: Total HCO3-:1 : 0.0001 > 1e-12 [absolute]",
            "FAIL: Total HCO3-:Mean : 0.0001 > 1e-12 [absolute]",
            "FAIL: Liquid Saturation:81 : 0.1 > 1e-12 [absolute]",
        ]

        self.failure_strings_performance = [
            "FAIL: Transport:Solution 2-Norm : 0.00137252117371 > 1e-12 [absolute]",
            "FAIL: Transport:Solver Iterations : 1800 > 0 [absolute]",
            "FAIL: Transport:Time Steps : 90 > 0.0 [absolute]",
            "FAIL: Transport:Residual 2-Norm : 1.45584930368e+18 > 1e-12 [absolute]",
            "FAIL: Transport:Time Step Cuts : 2 > 0 [absolute]",
            "FAIL: Transport:Time (seconds) : 99.9999999995 > 1.0 [percent]",
            "FAIL: Transport:Newton Iterations : 180 > 0 [absolute]",
        ]

    def tearDown(self):
        pass
        #subprocess.call(self.clean_cmd, shell=False, stdout=self.test_stdout)
        #self.test_stdout.close()

    def find_testlog(self):
        testlog_re = re.compile(r"^pflotran-tests-[\w_-]+\.testlog")
        testlog = None
        for entry in os.listdir(os.getcwd()):
            match = testlog_re.search(entry)
            if match is not None:
                testlog = entry
                break
        return testlog

    def search_testlog(self, testlog_name, failure_string):
        logdata = None
        with open(testlog_name, 'r') as testlog:
            logdata = testlog.readlines()
        failure_found = False
        for line in logdata:
            if line.find(failure_string) != -1:
                failure_found = True
                break
        return failure_found

    def check_testlog(self, testlog, failure_string):
        failure_found = self.search_testlog(testlog, failure_string)
        self.assertEqual(failure_found, True, failure_string)

    def run_regression_test_manager(self):
        subprocess.call(self.run_cmd, shell=False, stdout=self.test_stdout)
        testlog = self.find_testlog()
        self.assertNotEqual(testlog, None)
        return testlog

    def test_xfail_no_performance(self):
        """
        Run the test manager without performance checks on "solution" blocks
        """
        testlog = self.run_regression_test_manager()
        for failure in self.failure_strings_general:
            self.check_testlog(testlog, failure)
        self.check_testlog(testlog, "Skipping SOLUTION : Transport")

    def test_xfail_performance(self):
        """
        Run the test manager with performance check turned on
        """
        self.run_cmd.append("--check-performance")
        testlog = self.run_regression_test_manager()
        for failure in self.failure_strings_general:
            self.check_testlog(testlog, failure)
        for failure in self.failure_strings_performance:
            self.check_testlog(testlog, failure)

    # NOTE: tell nose to ignore these non-test functions
    find_testlog.__test__ = False
    search_testlog.__test__ = False
    check_testlog.__test__ = False
    run_regression_test_manager.__test__ = False


if __name__ == '__main__':
    unittest.main(buffer=True)
    #unittest.main()
