#!/usr/bin/env python
"""
Unit test suite for the pflotran regression test manager

"""

from __future__ import print_function

import os
import re
import subprocess
import unittest

from regression_tests import RegressionTest


class RegressionTest_SetTestData(unittest.TestCase):
    """Tests to verify that the default parameters are correctly set and
    overwritten by user supplied values when appropriate.

    """

    def setUp(self):
        self.rt = RegressionTest()
        # test_data is the test section from a config file
        self.test_data = {"name": "dummy_test_name", }
        # timeout as specified by the command line arg
        self.timeout = None
        # check performance as specified by the command line arg
        self.check_performance = False
        # criteria is the default-test-criteria section from the config file
        self.criteria = {}
        self.testlog = open("dummy.testlog", 'w')

    def tearDown(self):
        self.testlog.close()

    #--- test name ------------------------------------------------------------
    def test_no_test_name(self):
        self.test_data = {}
        self.assertRaises(KeyError, self.rt.setup,
                          self.criteria, self.test_data, self.timeout,
                          self.check_performance, self.testlog)

    def test_set_test_name(self):
        self.rt.setup(self.criteria, self.test_data, self.timeout,
                      self.check_performance, self.testlog)
        self.assertEqual(self.rt.name(), "dummy_test_name")

    #--- timeout --------------------------------------------------------------
    def test_timeout_default(self):
        self.assertEqual(self.rt._timeout, 60.0)

    def test_timeout_cfg(self):
        # verify that test specifc timeout is assigned
        self.test_data["timeout"] = 123.0
        self.rt.setup(self.criteria, self.test_data, self.timeout,
                      self.check_performance, self.testlog)
        self.assertEqual(self.rt._timeout, 123.0)

    def test_timeout_cmdline(self):
        # verify that commandline timeout is assigned
        self.timeout = [234.0]
        self.rt.setup(self.criteria, self.test_data, self.timeout,
                      self.check_performance, self.testlog)
        self.assertEqual(self.rt._timeout, 234.0)

    def test_timeout_cmdline_override(self):
        # verify that commandline timeout over rides the test specific
        # timeout
        self.test_data["timeout"] = 123.0
        self.timeout = [234.0]
        self.rt.setup(self.criteria, self.test_data, self.timeout,
                      self.check_performance, self.testlog)
        self.assertEqual(self.rt._timeout, 234.0)

    #--- check-performance flag -----------------------------------------------
    def test_check_performance_default(self):
        self.assertEqual(self.rt._check_performance, False)

    def test_check_performance_cmdline(self):
        # verify that command line flag overrides the default
        self.check_performance = True
        self.rt.setup(self.criteria, self.test_data, self.timeout,
                      self.check_performance, self.testlog)
        self.assertEqual(self.rt._check_performance, True)

    #--- np -------------------------------------------------------------------
    def test_np_default(self):
        self.assertEqual(self.rt._np, None)

    def test_np_cfg(self):
        # verify that np from config file get set correctly
        self.test_data["np"] = '2'
        self.rt.setup(self.criteria, self.test_data, self.timeout,
                      self.check_performance, self.testlog)
        self.assertEqual(self.rt._np, '2')

    #--- pflotran commandline args --------------------------------------------
    def test_pflotran_args_default(self):
        self.assertEqual(self.rt._pflotran_args, None)

    def test_pflotran_args_flag(self):
        self.test_data['input_arguments'] = "-foo"
        self.rt.setup(self.criteria, self.test_data, self.timeout,
                      self.check_performance, self.testlog)
        self.assertEqual(self.rt._pflotran_args, ['-foo'])

    def test_pflotran_args_flag_param(self):
        self.test_data['input_arguments'] = "-bar baz"
        self.rt.setup(self.criteria, self.test_data, self.timeout,
                      self.check_performance, self.testlog)
        self.assertEqual(self.rt._pflotran_args, ['-bar', 'baz'])

    def test_pflotran_args_stochastic_no_realizations(self):
        # stochastic w/o num_realizations is an error
        self.test_data['input_arguments'] = "-stochastic"
        self.assertRaises(Exception, self.rt.setup,
                          self.criteria, self.test_data, self.timeout,
                          self.check_performance, self.testlog)

    def test_pflotran_args_stochastic_empty_realizations(self):
        # stochastic requires an int num_realizations
        self.test_data['input_arguments'] = "-stochastic -num_realizations"
        self.assertRaises(Exception, self.rt.setup,
                          self.criteria, self.test_data, self.timeout,
                          self.check_performance, self.testlog)

    def test_pflotran_args_stochastic_str_realizations(self):
        # stochastic requires an int num_realizations
        self.test_data['input_arguments'] = "-stochastic -num_realizations foo"
        self.assertRaises(ValueError, self.rt.setup,
                          self.criteria, self.test_data, self.timeout,
                          self.check_performance, self.testlog)

    def test_pflotran_args_stochastic_realizations(self):
        # stochastic requires an int num_realizations
        self.test_data['input_arguments'] = "-stochastic -num_realizations 3"
        self.rt.setup(self.criteria, self.test_data, self.timeout,
                      self.check_performance, self.testlog)
        self.assertEqual(self.rt._stochastic_realizations, 3)

    #--- set_criteria ---------------------------------------------------------
    # verify that RegressionTest._set_criteria() is assigning test
    # criteria with the expected priority
    #
    # for these tests "criteria" is the default criteria from the
    # config file
    def test_set_criteria_invalid(self):
        key = "bicycles"
        self.assertRaises(Exception, self.rt._set_criteria, key,
                          self.criteria, self.test_data)

    def test_set_criteria_default(self):
        key = "concentration"
        self.rt._set_criteria(key, self.criteria, self.test_data)
        self.assertEqual(self.rt._tolerance[key], [1.0e-12, "absolute"])

    def test_set_criteria_cfg_default(self):
        # cfg criteria overrides default
        key = "concentration"
        self.criteria[key] = "1.0e-5 relative"
        self.rt._set_criteria(key, self.criteria, self.test_data)
        self.assertEqual(self.rt._tolerance[key], [1.0e-5, "relative"])

    def test_set_criteria_test_default(self):
        # individual test criteria overrides default
        key = "concentration"
        self.test_data[key] = "1.0e-7 percent"
        self.rt._set_criteria(key, self.criteria, self.test_data)
        self.assertEqual(self.rt._tolerance[key], [1.0e-7, "percent"])

    def test_set_criteria_test_cfg_default(self):
        # individual test criteria overrides cfg default overrides
        # global default
        key = "concentration"
        self.test_data[key] = "1.0e-7 percent"
        self.criteria[key] = "1.0e-5 relative"
        self.rt._set_criteria(key, self.criteria, self.test_data)
        self.assertEqual(self.rt._tolerance[key], [1.0e-7, "percent"])

    #--- validate_criteria ----------------------------------------------------
    def test_validate_criteria_absolute(self):
        key = "velocity"
        criteria_str = "1.0e-4 absolute"
        value, criteria = self.rt._validate_criteria(key, criteria_str)
        self.assertAlmostEqual(value, 1.0e-4, delta=1.0e-16)
        self.assertEqual(criteria, "absolute")

    def test_validate_criteria_relative(self):
        key = "velocity"
        criteria_str = "1.0e-4 relative"
        value, criteria = self.rt._validate_criteria(key, criteria_str)
        self.assertAlmostEqual(value, 1.0e-4, delta=1.0e-16)
        self.assertEqual(criteria, "relative")

    def test_validate_criteria_percent(self):
        key = "velocity"
        criteria_str = "20 percent"
        value, criteria = self.rt._validate_criteria(key, criteria_str)
        self.assertAlmostEqual(value, 20.0, delta=1.0e-16)
        self.assertEqual(criteria, "percent")

    def test_validate_criteria_incorrect(self):
        key = "velocity"
        criteria_str = "1.0e-4 junk"
        self.assertRaises(Exception,
                          self.rt._validate_criteria, key, criteria_str)

    def test_validate_criteria_value_str(self):
        key = "velocity"
        criteria_str = "absolute 1.0e-4"
        self.assertRaises(Exception,
                          self.rt._validate_criteria, key, criteria_str)


class RegressionTest_GetSections(unittest.TestCase):
    """Tests to verify that get_sections is correctly splitting up
    regression files into sections.

    get_sections takes a list of lines (file.readlines()) as input and
    returns a dictionary of sections.

    """

    def setUp(self):
        self.rt = RegressionTest()
        # test_data is the test section from a config file
        self.test_data = {"name": "dummy_test_name", }
        self.testlog = open("dummy.testlog", 'w')

    def tearDown(self):
        self.testlog.close()

    def test_get_section(self):
        """Verify that the compare_sections raises an exception when
        different sections are compared.

        """
        filelines = ["-- GENERIC: pH --",
                     "Max:  4.6763373466787E+00",
                     "Min:  4.6763373466787E+00",
                     "Mean:  4.6763373466787E+00",
                     "1:  4.6763373466787E+00",
                     "-- CONCENTRATION: Total H+ --",
                     "Max:  1.0000000000000E-03",
                     "Min:  1.0000000000000E-03",
                     "Mean:  1.0000000000000E-03",
                     "1:  1.0000000000000E-03",
                     "-- SOLUTION: Transport --",
                     "Time (seconds):  0.0000000000000E+00",
                     "Time Steps:            0",
                     "Newton Iterations:            0",
                     "Solver Iterations:            0",
                     "Time Step Cuts:            0",
                     "Solution 2-Norm:  2.9796944879051E-05",
                     "Residual 2-Norm:  0.0000000000000E+00",
                     "-- GENERIC: LIQUID VELOCITY [m/y] --",
                     "1:   0.0000000000000E+00  0.0000000000000E+00 -9.9999573084086E-02",
                     "21:   0.0000000000000E+00  0.0000000000000E+00 -9.9999990931782E-02",
                     "41:   0.0000000000000E+00  0.0000000000000E+00 -9.9999991630680E-02",
                     "61:   0.0000000000000E+00  0.0000000000000E+00 -9.9999994712625E-02",
                     "81:   0.0000000000000E+00  1.0000000000000E+00 -9.9999998335505E-02",
                 ]
        pH = {"name": "pH",
              "type": "GENERIC",
              "Max": "4.6763373466787E+00",
              "Min": "4.6763373466787E+00",
              "Mean": "4.6763373466787E+00",
              "1": "4.6763373466787E+00"
        }
        total_h = {"type": "CONCENTRATION",
                   "name": "Total H+",
                   "Max": "1.0000000000000E-03",
                   "Min": "1.0000000000000E-03",
                   "Mean": "1.0000000000000E-03",
                   "1": "1.0000000000000E-03",
        }
        transport = {"type": "SOLUTION",
                     "name": "Transport",
                     "Time (seconds)": "0.0000000000000E+00",
                     "Time Steps": "0",
                     "Newton Iterations": "0",
                     "Solver Iterations": "0",
                     "Time Step Cuts": "0",
                     "Solution 2-Norm": "2.9796944879051E-05",
                     "Residual 2-Norm": "0.0000000000000E+00",
                 }
        velocity = {"type": "GENERIC",
                    "name": "LIQUID VELOCITY [m/y]",
                    "1": "0.0000000000000E+00  0.0000000000000E+00 -9.9999573084086E-02",
                    "21": "0.0000000000000E+00  0.0000000000000E+00 -9.9999990931782E-02",
                    "41": "0.0000000000000E+00  0.0000000000000E+00 -9.9999991630680E-02",
                    "61": "0.0000000000000E+00  0.0000000000000E+00 -9.9999994712625E-02",
                    "81": "0.0000000000000E+00  1.0000000000000E+00 -9.9999998335505E-02",
                }
        sections = self.rt._get_sections(filelines)
        self.assertTrue("pH" in sections)
        self.assertDictEqual(sections["pH"], pH)
        self.assertTrue("Total H+" in sections)
        self.assertDictEqual(sections["Total H+"], total_h)
        self.assertTrue("Transport" in sections)
        self.assertDictEqual(sections["Transport"], transport)
        self.assertTrue("LIQUID VELOCITY [m/y]" in sections)
        self.assertDictEqual(sections["LIQUID VELOCITY [m/y]"], velocity)



class RegressionTest_CompareSections(unittest.TestCase):
    """Tests to verify that section comparisons between regression and gold
    data are correct.

    """

    def setUp(self):
        self.rt = RegressionTest()
        # test_data is the test section from a config file
        self.test_data = {"name": "dummy_test_name", }
        self.testlog = open("dummy.testlog", 'w')

    def tearDown(self):
        self.testlog.close()

    def test_compare_different_sections(self):
        """Verify that the compare_sections raises an exception when
        different sections are compared.

        """
        gold_section = {'name': 'Cats',
                        'type': 'generic',
                        'Mean': '0.025'}
        current_section = {'name': 'Dogs',
                           'type': 'generic',
                           'Mean': '0.025'}
        self.assertRaises(Exception, self.rt._compare_sections,
                          gold_section, current_section, self.testlog)


    def test_same_section(self):
        """Verify that the comparison passes if the same
        section info is passed for both current and gold.

        Note: return value "status" == number of failed comparisons
        """
        gold_section = {'name': 'Cats',
                        'type': 'generic',
                        'Min': '1.0e-2',
                        'Max': '1.0e-1',
                        '1': '0.05',
                        'Mean': '0.025'}
        num_failures = self.rt._compare_sections(
            gold_section, gold_section, self.testlog)
        self.assertEqual(num_failures, 0)

    def test_current_section_missing_field(self):
        """Failure if a section is missing a field
        """
        gold_section = {'name': 'Cats',
                        'type': 'generic',
                        'Min': '1.0e-2',
                        'Max': '1.0e-1'}
        current_section = {'name': 'Cats',
                           'type': 'generic',
                           'Min': '1.0e-2'}
        num_failures = self.rt._compare_sections(gold_section, current_section,
                                                 self.testlog)
        self.assertEqual(num_failures, 1)

    def test_current_section_extra_field(self):
        """Failure if a section has an extra field
        """
        gold_section = {'name': 'Cats',
                        'type': 'generic',
                        'Max': '1.0e-1'}
        current_section = {'name': 'Cats',
                           'type': 'generic',
                           'Min': '1.0e-2',
                           'Max': '1.0e-1'}
        num_failures = self.rt._compare_sections(gold_section, current_section,
                                                 self.testlog)
        self.assertEqual(num_failures, 1)


    def test_compare_skip_performance_pass(self):
        """Passes when performance sections are skipped, even if the results
        are different.

        """
        gold_section = {'name': 'Transport',
                        'type': 'SOLUTION',
                        'Time (seconds)': '1.3227415084839E-01',
                        'Newton Iterations': '200'}
        current_section = {'name': 'Transport',
                           'type': 'SOLUTION',
                           'Time (seconds)': '1.3227415084839E-01',
                           'Newton Iterations': '100'}
        self.rt._check_performance = False
        num_failures = self.rt._compare_sections(gold_section, current_section,
                                                 self.testlog)
        self.assertEqual(num_failures, 0)

    def test_compare_performance_fail(self):
        """Fails when performance sections are different.

        """
        gold_section = {'name': 'Transport',
                        'type': 'SOLUTION',
                        'Time (seconds)': '1.3227415084839E-01',
                        'Newton Iterations': '200'}
        current_section = {'name': 'Transport',
                           'type': 'SOLUTION',
                           'Time (seconds)': '1.3227415084839E-01',
                           'Newton Iterations': '100'}
        self.rt._check_performance = True
        num_failures = self.rt._compare_sections(gold_section, current_section,
                                                 self.testlog)
        self.assertEqual(num_failures, 1)

    def test_compare_vectors(self):
        """Failure if vector lengths are different
        """
        gold_section = {
            'name': 'LIQUID VELOCITY [m/y]',
            'type': 'GENERIC',
            '41': '0.0000000000000E+00  0.0000000000000E+00 -9.9999991630680E-02',
            '1': '0.0000000000000E+00  0.0000000000000E+00 -9.9999573084086E-02'}
        current_section = {
            'name': 'LIQUID VELOCITY [m/y]',
            'type': 'GENERIC',
            '41': '0.0000000000000E+00  0.0000000000000E+00 -9.9999991630680E-02',
            '1': '0.0000000000000E+00  0.0000000000000E+00'}
        num_failures = self.rt._compare_sections(gold_section, current_section,
                                                 self.testlog)
        self.assertEqual(num_failures, 1)

    def test_compare_values_error(self):
        """Failure when compare_values throws an exception. 

        type: pets will cause compare_values to throw an exception.

        """
        gold_section = {'name': 'Cats',
                        'type': 'pets',
                        'Max': '1.0e-1'}
        self.rt._check_performance = False
        num_failures = self.rt._compare_sections(gold_section, gold_section,
                                                 self.testlog)
        self.assertEqual(num_failures, 1)


class RegressionTest_CompareValues(unittest.TestCase):
    """Tests to verify that value comparisons between regression and gold
    data are correct.

    """

    def setUp(self):
        self.rt = RegressionTest()
        # test_data is the test section from a config file
        self.test_data = {"name": "dummy_test_name", }
        self.testlog = open("dummy.testlog", 'w')

    def tearDown(self):
        self.testlog.close()

    def test_compare_invalid_type(self):
        """Failure if unknown data type.
        """
        name_str = "Cats"
        data_type = "pets"
        gold_value = "2.6"
        self.assertRaises(Exception, self.rt._compare_values, name_str,
                          data_type, gold_value, gold_value, self.testlog)

    def test_compare_invalid_tolerance_type(self):
        """Failure if unknown data type.
        """
        name_str = "Cats"
        data_type = "generic"
        gold_value = "2.6"
        self.rt._tolerance["generic"] = [1.0e-5, "pets"]
        self.assertRaises(Exception, self.rt._compare_values, name_str,
                          data_type, gold_value, gold_value, self.testlog)

    def test_compare_generic_different(self):
        """Failure if generic values differ.
        """
        name_str = "Cats"
        data_type = "generic"
        gold_value = "2.6"
        current_value = "6.2"
        status = self.rt._compare_values(name_str, data_type, gold_value,
                                         current_value, self.testlog)
        self.assertEqual(status, 1)

    def test_compare_generic_same(self):
        """Pass if generic values differ.
        """
        name_str = "Cats"
        data_type = "generic"
        gold_value = "4.5"
        current_value = "4.5"
        status = self.rt._compare_values(name_str, data_type, gold_value,
                                         current_value, self.testlog)
        self.assertEqual(status, 0)

    def test_compare_relative_gold_zero(self):
        """Pass if gold value is zero.
        """
        name_str = "Cats"
        data_type = "generic"
        gold_value = "0.0"
        current_value = "1.0e-9"
        self.rt._tolerance["generic"] = [1.0e-5, "relative"]
        status = self.rt._compare_values(name_str, data_type, gold_value,
                                         current_value, self.testlog)
        self.assertEqual(status, 1)

    def test_compare_relative_current_zero(self):
        """Pass if current value is zero
        """
        name_str = "Cats"
        data_type = "generic"
        gold_value = "1.0e-9"
        current_value = "0.0"
        self.rt._tolerance["generic"] = [1.0e-5, "relative"]
        status = self.rt._compare_values(name_str, data_type, gold_value,
                                         current_value, self.testlog)
        self.assertEqual(status, 1)

    def test_compare_relative_both_zero(self):
        """Pass if both values zero
        """
        name_str = "Cats"
        data_type = "generic"
        gold_value = "0.0"
        current_value = "0.0"
        self.rt._tolerance["generic"] = [1.0e-5, "relative"]
        status = self.rt._compare_values(name_str, data_type, gold_value,
                                         current_value, self.testlog)
        self.assertEqual(status, 0)

    def test_compare_values_discrete_int(self):
        """Correctly compare discrete values
        """
        name_str = "Cats"
        data_type = "discrete"
        gold_value = "0"
        current_value = "0"
        status = self.rt._compare_values(name_str, data_type, gold_value,
                                         current_value, self.testlog)
        self.assertEqual(status, 0)

    def test_compare_values_discrete_int_different(self):
        """Correctly compare discrete values
        """
        name_str = "Cats"
        data_type = "discrete"
        gold_value = "0"
        current_value = "1"
        status = self.rt._compare_values(name_str, data_type, gold_value,
                                         current_value, self.testlog)
        self.assertEqual(status, 1)

    def test_compare_values_discrete_float(self):
        """Correctly compare discrete values
        """
        name_str = "Mean"
        data_type = "discrete"
        gold_value = "0.1"
        current_value = "0.1"
        status = self.rt._compare_values(name_str, data_type, gold_value,
                                         current_value, self.testlog)
        self.assertEqual(status, 0)

    def test_compare_values_discrete_float_different(self):
        """Correctly compare discrete values
        """
        name_str = "Mean"
        data_type = "discrete"
        gold_value = "0.1"
        current_value = "0.11"
        status = self.rt._compare_values(name_str, data_type, gold_value,
                                         current_value, self.testlog)
        self.assertEqual(status, 1)

class RegressionTest_CompareDiscrete(unittest.TestCase):
    """Tests to verify correct identification of discrete and discrete
    mean values.

    """

    def setUp(self):
        self.rt = RegressionTest()
        # test_data is the test section from a config file
        self.test_data = {"name": "dummy_test_name", }
        self.testlog = open("dummy.testlog", 'w')

    def tearDown(self):
        self.testlog.close()

    def test_compare_discrete_int(self):
        """Correctly identify discrete int values, set tolerance appropriately.
        """
        name_str = "Cats"
        gold_value = "1"
        previous, current, tolerance_type, tolerance = (
            self.rt._compare_discrete(name_str, gold_value, gold_value))
        self.assertEqual(tolerance_type, "absolute")
        self.assertEqual(tolerance, 0)
        self.assertTrue(isinstance(previous, int))
        self.assertTrue(isinstance(current, int))

    def test_compare_discrete_mean(self):
        """Correcty identify discrete mean values, set tolerance appropriately.
        """
        name_str = "Mean"
        gold_value = "1.23"
        previous, current, tolerance_type, tolerance = (
            self.rt._compare_discrete(name_str, gold_value, gold_value))
        self.assertEqual(tolerance_type, "absolute")
        self.assertEqual(tolerance, 1.0e-12)
        self.assertTrue(isinstance(previous, float))
        self.assertTrue(isinstance(current, float))

    def test_compare_discrete_float_gold(self):
        """Fail if discrete gold is not an int and not flagged as Mean
        """
        name_str = "Cats"
        gold_value = "2.34"
        current_value = "1"
        self.assertRaises(Exception, self.rt._compare_discrete, name_str,
                          gold_value, current_value)

    def test_compare_discrete_float_current(self):
        """Fail if discrete current is not an int and not flagged as Mean
        """
        name_str = "Cats"
        gold_value = "2"
        current_value = "4.5"
        self.assertRaises(Exception, self.rt._compare_discrete, name_str,
                          gold_value, current_value)


class RegressionTest_CompareSolution(unittest.TestCase):
    """Tests to verify correct identification of solution block types
    are correct.

    """

    def setUp(self):
        self.rt = RegressionTest()
        # test_data is the test section from a config file
        self.test_data = {"name": "dummy_test_name", }
        self.testlog = open("dummy.testlog", 'w')

    def tearDown(self):
        self.testlog.close()

    def test_compare_solution_time(self):
        """Correctly identify solution time and set tolerance appropriately.
        """
        name_str = "Transport : Time (seconds)"
        gold_value = "1.0"
        previous, current, tolerance_type, tolerance = (
            self.rt._compare_solution(name_str, gold_value, gold_value))
        self.assertEqual(tolerance_type, "percent")
        self.assertEqual(tolerance, 5.0)
        self.assertTrue(isinstance(previous, float))
        self.assertTrue(isinstance(current, float))

    def test_compare_solution_time_steps(self):
        """Correctly identify time steps, set tolerance appropriately.
        """
        name_str = "Transport : Time Steps"
        gold_value = "10"
        previous, current, tolerance_type, tolerance = (
            self.rt._compare_solution(name_str, gold_value, gold_value))
        self.assertEqual(tolerance_type, "absolute")
        self.assertEqual(tolerance, 0)
        self.assertTrue(isinstance(previous, int))
        self.assertTrue(isinstance(current, int))

    def test_compare_solution_solver_iterations(self):
        """Correctly identify solver iterations and set tolerance appropriately.
        """
        name_str = "Transport : Solver Iterations"
        gold_value = "10"
        previous, current, tolerance_type, tolerance = (
            self.rt._compare_solution(name_str, gold_value, gold_value))
        self.assertEqual(tolerance_type, "absolute")
        self.assertEqual(tolerance, 0)
        self.assertTrue(isinstance(previous, int))
        self.assertTrue(isinstance(current, int))

    def test_compare_solution_time_step_cuts(self):
        """Correctly identify time step cuts, set tolerance appropriately.
        """
        name_str = "Transport : Time Step Cuts"
        gold_value = "10"
        previous, current, tolerance_type, tolerance = (
            self.rt._compare_solution(name_str, gold_value, gold_value))
        self.assertEqual(tolerance_type, "absolute")
        self.assertEqual(tolerance, 0)
        self.assertTrue(isinstance(previous, int))
        self.assertTrue(isinstance(current, int))

    def test_compare_solution_solution_2norm(self):
        """Correctly identify solution 2norm and set tolerance appropriately.
        """
        name_str = "Transport : Solution 2-Norm"
        gold_value = "1.0"
        previous, current, tolerance_type, tolerance = (
            self.rt._compare_solution(name_str, gold_value, gold_value))
        self.assertEqual(tolerance_type, "absolute")
        self.assertEqual(tolerance, 1.0e-12)
        self.assertTrue(isinstance(previous, float))
        self.assertTrue(isinstance(current, float))

    def test_compare_solution_residual_2norm(self):
        """Correctly identify residual 2norm and set tolerance appropriately.
        """
        name_str = "Transport : Residual 2-Norm"
        gold_value = "1.0"
        previous, current, tolerance_type, tolerance = (
            self.rt._compare_solution(name_str, gold_value, gold_value))
        self.assertEqual(tolerance_type, "absolute")
        self.assertEqual(tolerance, 1.0e-12)
        self.assertTrue(isinstance(previous, float))
        self.assertTrue(isinstance(current, float))

    def test_compare_solution_invalid(self):
        """Correctly identify residual 2norm and set tolerance appropriately.
        """
        name_str = "Transport : Bicycles"
        gold_value = "1.0"
        self.assertRaises(Exception, self.rt._compare_solution, name_str,
                          gold_value, gold_value)


if __name__ == '__main__':
    #unittest.main(buffer=True)
    unittest.main()
