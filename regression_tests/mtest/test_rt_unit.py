#!/usr/bin/env python
"""
Unit test suite for the pflotran regression test manager

"""

from __future__ import print_function

import os
import re
import struct
import subprocess
import sys
import unittest

import h5py

from regression_tests import RegressionTest, TestStatus


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

    #--- restart timestep -----------------------------------------------------
    def test_restart_timestep_default(self):
        self.assertEqual(self.rt._pflotran_args, None)

    def test_restart_timestep_default2(self):
        self.rt.setup(self.criteria, self.test_data, self.timeout,
                      self.check_performance, self.testlog)
        self.assertEqual(self.rt._pflotran_args, None)

    def test_restart_timestep_int(self):
        self.test_data['restart_timestep'] = 10
        self.rt.setup(self.criteria, self.test_data, self.timeout,
                      self.check_performance, self.testlog)
        self.assertEqual(self.rt._restart_timestep, 10)

    def test_restart_timestep_junk(self):
        self.test_data['restart_timestep'] = "foo bar baz"
        self.assertRaises(ValueError, self.rt.setup,
                          self.criteria, self.test_data, self.timeout,
                          self.check_performance, self.testlog)

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
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_VALUE], 1.0e-12)
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_TYPE], "absolute")
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_MIN_THRESHOLD], 0.0)
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_MAX_THRESHOLD], sys.float_info.max)

    def test_set_criteria_cfg_default(self):
        # cfg criteria overrides default
        key = "concentration"
        self.criteria[key] = "1.0e-5 relative , min_threshold 1.0e-2 , max_threshold 34.5"
        self.rt._set_criteria(key, self.criteria, self.test_data)
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_VALUE], 1.0e-5)
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_TYPE], "relative")
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_MIN_THRESHOLD], 1.0e-2)
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_MAX_THRESHOLD], 34.5)

    def test_set_criteria_test_default(self):
        # individual test criteria overrides default
        key = "concentration"
        self.test_data[key] = "1.0e-7 percent , min_threshold 1.0e-11"
        self.rt._set_criteria(key, self.criteria, self.test_data)
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_VALUE], 1.0e-7)
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_TYPE], "percent")
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_MIN_THRESHOLD], 1.0e-11)
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_MAX_THRESHOLD], sys.float_info.max)

    def test_set_criteria_test_cfg_default(self):
        # individual test criteria overrides cfg default overrides
        # global default
        key = "concentration"
        self.test_data[key] = "1.0e-7 percent"
        self.criteria[key] = "1.0e-5 relative , max_threshold 1234.0"
        self.rt._set_criteria(key, self.criteria, self.test_data)
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_VALUE], 1.0e-7)
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_TYPE], "percent")
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_MIN_THRESHOLD], 0.0)
        self.assertEqual(self.rt._tolerance[key][self.rt._TOL_MAX_THRESHOLD], sys.float_info.max)

    def test_set_criteria_cfg_empty_default(self):
        # an exception is thrown when cfg default criteria is not iterable (a dict)
        key = "concentration"
        default_cfg_criteria = None
        self.assertRaises(Exception,
                          self.rt._set_criteria, key, default_cfg_criteria, "cfg default criteria not a dict")

    #--- validate_criteria ----------------------------------------------------
    def test_validate_criteria_absolute(self):
        key = "velocity"
        criteria_str = "1.0e-4 absolute"
        criteria = self.rt._validate_criteria(key, criteria_str)
        self.assertAlmostEqual(criteria[self.rt._TOL_VALUE], 1.0e-4, delta=1.0e-16)
        self.assertEqual(criteria[self.rt._TOL_TYPE], "absolute")
        self.assertEqual(criteria[self.rt._TOL_MIN_THRESHOLD], None)
        self.assertEqual(criteria[self.rt._TOL_MAX_THRESHOLD], None)

    def test_validate_criteria_relative(self):
        key = "velocity"
        criteria_str = "1.0e-4 relative"
        criteria = self.rt._validate_criteria(key, criteria_str)
        self.assertAlmostEqual(criteria[self.rt._TOL_VALUE], 1.0e-4, delta=1.0e-16)
        self.assertEqual(criteria[self.rt._TOL_TYPE], "relative")
        self.assertEqual(criteria[self.rt._TOL_MIN_THRESHOLD], None)
        self.assertEqual(criteria[self.rt._TOL_MAX_THRESHOLD], None)

    def test_validate_criteria_percent(self):
        key = "velocity"
        criteria_str = "20 percent"
        criteria = self.rt._validate_criteria(key, criteria_str)
        self.assertAlmostEqual(criteria[self.rt._TOL_VALUE], 20.0, delta=1.0e-16)
        self.assertEqual(criteria[self.rt._TOL_TYPE], "percent")
        self.assertEqual(criteria[self.rt._TOL_MIN_THRESHOLD], None)
        self.assertEqual(criteria[self.rt._TOL_MAX_THRESHOLD], None)

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

    def test_validate_criteria_value_error(self):
        key = "velocity"
        criteria_str = "1.23,4 absolute"
        self.assertRaises(Exception,
                          self.rt._validate_criteria, key, criteria_str)

    def test_validate_criteria_min_threshold(self):
        key = "velocity"
        criteria_str = "20 percent , min_threshold 1.0e-11"
        criteria = self.rt._validate_criteria(key, criteria_str)
        self.assertAlmostEqual(criteria[self.rt._TOL_VALUE], 20.0, delta=1.0e-16)
        self.assertEqual(criteria[self.rt._TOL_TYPE], "percent")
        self.assertEqual(criteria[self.rt._TOL_MIN_THRESHOLD], 1.0e-11)
        self.assertEqual(criteria[self.rt._TOL_MAX_THRESHOLD], None)

    def test_validate_criteria_max_threshold(self):
        key = "velocity"
        criteria_str = "20 percent , max_threshold 1.23e4"
        criteria = self.rt._validate_criteria(key, criteria_str)
        self.assertAlmostEqual(criteria[self.rt._TOL_VALUE], 20.0, delta=1.0e-16)
        self.assertEqual(criteria[self.rt._TOL_TYPE], "percent")
        self.assertEqual(criteria[self.rt._TOL_MIN_THRESHOLD], None)
        self.assertEqual(criteria[self.rt._TOL_MAX_THRESHOLD], 1.23e4)

    def test_validate_criteria_min_max_threshold(self):
        key = "velocity"
        criteria_str = "20 percent , max_threshold 1.23e4 , min_threshold 2.34e-5"
        criteria = self.rt._validate_criteria(key, criteria_str)
        self.assertAlmostEqual(criteria[self.rt._TOL_VALUE], 20.0, delta=1.0e-16)
        self.assertEqual(criteria[self.rt._TOL_TYPE], "percent")
        self.assertEqual(criteria[self.rt._TOL_MIN_THRESHOLD], 2.34e-5)
        self.assertEqual(criteria[self.rt._TOL_MAX_THRESHOLD], 1.23e4)

    def test_validate_criteria_unknown_threshold_error(self):
        key = "velocity"
        criteria_str = "20 percent , cat 1234"
        self.assertRaises(Exception, self.rt._validate_criteria,
                          key ,criteria_str)

    def test_validate_criteria_threshold_no_number_error(self):
        key = "velocity"
        criteria_str = "20 percent , min_threshold whiskers"
        self.assertRaises(Exception, self.rt._validate_criteria,
                          key ,criteria_str)


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
        self.rt._tolerance["generic"] = [1.0e-5, "relative", 0.0, sys.float_info.max]
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
        self.rt._tolerance["generic"] = [1.0e-5, "relative", 0.0, sys.float_info.max]
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
        self.rt._tolerance["generic"] = [1.0e-5, "relative", 0.0, sys.float_info.max]
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

    def test_compare_values_skip_below_min_threshold(self):
        """Correctly pass (skip) when gold value is less than min_threshold,
        even if it would otherwise fail.

        """
        name_str = "Cats"
        data_type = "generic"
        gold_value = "1.0e-12"
        current_value = "1.0e-9"
        self.rt._tolerance["generic"] = [1.0e-5, "relative", 1.0e-10, sys.float_info.max]
        status = self.rt._compare_values(name_str, data_type, gold_value,
                                         current_value, self.testlog)
        self.assertEqual(status, 0)

    def test_compare_values_skip_above_max_threshold(self):
        """Correctly pass (skip) when gold value is larger than max_threshold,
        even if it would otherwise fail.

        """
        name_str = "Cats"
        data_type = "generic"
        gold_value = "10.0"
        current_value = "0.1"
        self.rt._tolerance["generic"] = [1.0e-5, "absolute", 0.0, 1.0]
        status = self.rt._compare_values(name_str, data_type, gold_value,
                                         current_value, self.testlog)
        self.assertEqual(status, 0)

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
        previous, current, tolerance = \
            self.rt._compare_discrete(name_str, gold_value, gold_value)
        self.assertEqual(tolerance[self.rt._TOL_TYPE], "absolute")
        self.assertEqual(tolerance[self.rt._TOL_VALUE], 0)
        self.assertTrue(isinstance(previous, int))
        self.assertTrue(isinstance(current, int))

    def test_compare_discrete_mean(self):
        """Correcty identify discrete mean values, set tolerance appropriately.
        """
        name_str = "Mean"
        gold_value = "1.23"
        previous, current, tolerance = \
            self.rt._compare_discrete(name_str, gold_value, gold_value)
        self.assertEqual(tolerance[self.rt._TOL_TYPE], "absolute")
        self.assertEqual(tolerance[self.rt._TOL_VALUE], 1.0e-12)
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
        previous, current, tolerance = \
            self.rt._compare_solution(name_str, gold_value, gold_value)
        self.assertEqual(tolerance[self.rt._TOL_TYPE], "percent")
        self.assertEqual(tolerance[self.rt._TOL_VALUE], 5.0)
        self.assertTrue(isinstance(previous, float))
        self.assertTrue(isinstance(current, float))

    def test_compare_solution_time_steps(self):
        """Correctly identify time steps, set tolerance appropriately.
        """
        name_str = "Transport : Time Steps"
        gold_value = "10"
        previous, current, tolerance = \
            self.rt._compare_solution(name_str, gold_value, gold_value)
        self.assertEqual(tolerance[self.rt._TOL_TYPE], "absolute")
        self.assertEqual(tolerance[self.rt._TOL_VALUE], 0)
        self.assertTrue(isinstance(previous, int))
        self.assertTrue(isinstance(current, int))

    def test_compare_solution_solver_iterations(self):
        """Correctly identify solver iterations and set tolerance appropriately.
        """
        name_str = "Transport : Solver Iterations"
        gold_value = "10"
        previous, current, tolerance = \
            self.rt._compare_solution(name_str, gold_value, gold_value)
        self.assertEqual(tolerance[self.rt._TOL_TYPE], "absolute")
        self.assertEqual(tolerance[self.rt._TOL_VALUE], 0)
        self.assertTrue(isinstance(previous, int))
        self.assertTrue(isinstance(current, int))

    def test_compare_solution_time_step_cuts(self):
        """Correctly identify time step cuts, set tolerance appropriately.
        """
        name_str = "Transport : Time Step Cuts"
        gold_value = "10"
        previous, current, tolerance = \
            self.rt._compare_solution(name_str, gold_value, gold_value)
        self.assertEqual(tolerance[self.rt._TOL_TYPE], "absolute")
        self.assertEqual(tolerance[self.rt._TOL_VALUE], 0)
        self.assertTrue(isinstance(previous, int))
        self.assertTrue(isinstance(current, int))

    def test_compare_solution_solution_2norm(self):
        """Correctly identify solution 2norm and set tolerance appropriately.
        """
        name_str = "Transport : Solution 2-Norm"
        gold_value = "1.0"
        previous, current, tolerance = \
            self.rt._compare_solution(name_str, gold_value, gold_value)
        self.assertEqual(tolerance[self.rt._TOL_TYPE], "absolute")
        self.assertEqual(tolerance[self.rt._TOL_VALUE], 1.0e-12)
        self.assertTrue(isinstance(previous, float))
        self.assertTrue(isinstance(current, float))

    def test_compare_solution_residual_2norm(self):
        """Correctly identify residual 2norm and set tolerance appropriately.
        """
        name_str = "Transport : Residual 2-Norm"
        gold_value = "1.0"
        previous, current, tolerance = \
            self.rt._compare_solution(name_str, gold_value, gold_value)
        self.assertEqual(tolerance[self.rt._TOL_TYPE], "absolute")
        self.assertEqual(tolerance[self.rt._TOL_VALUE], 1.0e-12)
        self.assertTrue(isinstance(previous, float))
        self.assertTrue(isinstance(current, float))

    def test_compare_solution_invalid(self):
        """Correctly identify residual 2norm and set tolerance appropriately.
        """
        name_str = "Transport : Bicycles"
        gold_value = "1.0"
        self.assertRaises(Exception, self.rt._compare_solution, name_str,
                          gold_value, gold_value)


class RegressionTest_CleanupGeneratedFiles(unittest.TestCase):
    """Test logic for pre-run file cleanup. cleanup should only touch
    files generated by the current test.

    """

    def setUp(self):
        self.testlog = open("dummy.testlog", 'w')
        self.status = TestStatus()
        self.rt = RegressionTest()
        # test_data is the test section from a config file
        test_data = {"name": "dummy_test_name"}
        # timeout as specified by the command line arg
        timeout = None
        # check performance as specified by the command line arg
        check_performance = False
        # criteria is the default-test-criteria section from the config file
        criteria = {}
        self.rt.setup(criteria, test_data, timeout,
                      check_performance, self.testlog)

        # list of files that should be preserved (not caputured by generic rules)
        self._save_suffixes = (".in", ".h5", ".txt", "-np8.in")
        self._generate_files(self._save_suffixes)

        # list of files that should always be moved
        self._rename_suffixes = [".regression", ".out", ".stdout", ]

    def tearDown(self):
        self.testlog.close()
        self._remove_files(self._save_suffixes)

    def _generate_files(self, suffixes):
        """
        Generate a bunch of files.
        """
        data = "12345"
        for suffix in suffixes:
            name = "{0}{1}".format(self.rt.name(), suffix)
            with open(name, "w") as tmpfile:
                tmpfile.write(data)
            name = "{0}-{1}{2}".format(self.rt._RESTART_PREFIX, self.rt.name(), suffix)
            with open(name, "w") as tmpfile:
                tmpfile.write(data)

    def _remove_files(self, suffixes):
        """Remove any remaining files
        """
        for suffix in suffixes:
            name = "{0}{1}".format(self.rt.name(), suffix)
            if os.path.isfile(name):
                os.remove(name)
            name = "{0}{1}.old".format(self.rt.name(), suffix)
            if os.path.isfile(name):
                os.remove(name)
            name = "{0}-{1}{2}".format(self.rt._RESTART_PREFIX, self.rt.name(), suffix)
            if os.path.isfile(name):
                os.remove(name)
            name = "{0}-{1}{2}.old".format(self.rt._RESTART_PREFIX, self.rt.name(), suffix)
            if os.path.isfile(name):
                os.remove(name)

    def test_cleanup_generated_files_standard(self):
        """test cleanup files for a standard run
        """
        self._generate_files(self._rename_suffixes)
        self.rt._cleanup_generated_files()

        for suffix in self._save_suffixes:
            # verify the "save" files have been preserved
            name = "{0}{1}".format(self.rt.name(), suffix)
            self.assertTrue(os.path.isfile(name), name)

        for suffix in self._rename_suffixes:
            # verify the generated files have been moved to ".old"
            name = "{0}{1}".format(self.rt.name(), suffix)
            self.assertFalse(os.path.isfile(name), name)
            self.assertTrue(os.path.isfile(name + ".old"), name)

        self._remove_files(self._rename_suffixes)

    def test_cleanup_generated_files_stochastic(self):
        """test cleanup files for a stochastic run
        """
        self.rt._stochastic_realizations = 2
        self._rename_suffixes.extend(["R1.regression", "R2.regression",
                                      "R1.out", "R2.out"])
        self._generate_files(self._rename_suffixes)
        self.rt._cleanup_generated_files()

        for suffix in self._save_suffixes:
            # verify the "save" files have been preserved
            name = "{0}{1}".format(self.rt.name(), suffix)
            self.assertTrue(os.path.isfile(name), name)

        for suffix in self._rename_suffixes:
            # verify the generated files have been moved to ".old"
            name = "{0}{1}".format(self.rt.name(), suffix)
            self.assertFalse(os.path.isfile(name), name)
            self.assertTrue(os.path.isfile(name + ".old"), name)
            
        self._remove_files(self._rename_suffixes)

    def test_cleanup_generated_files_restart(self):
        """test file cleanup for a restart run
        """
        self.rt._restart_timestep = 10
        self._rename_suffixes.extend(["-restart.chk", "-5.chk", "-10.chk"])
        self._generate_files(self._rename_suffixes)
        self.rt._cleanup_generated_files()

        for suffix in self._save_suffixes:
            # verify the "save" files have been preserved
            name = "{0}{1}".format(self.rt.name(), suffix)
            self.assertTrue(os.path.isfile(name), name)

        # temp input file has been moved
        name = "{0}-{1}.in".format(self.rt._RESTART_PREFIX, self.rt.name())
        self.assertFalse(os.path.isfile(name), name)
        self.assertTrue(os.path.isfile(name + ".old"), name)
        for suffix in self._rename_suffixes:
            # verify the generated files (orig run) have been moved to ".old"
            name = "{0}{1}".format(self.rt.name(), suffix)
            self.assertFalse(os.path.isfile(name), name)
            self.assertTrue(os.path.isfile(name + ".old"), name)
            # verify the generated files (restart run) have been moved to ".old"
            name = "{0}-{1}{2}".format(self.rt._RESTART_PREFIX, self.rt.name(), suffix)
            self.assertFalse(os.path.isfile(name), name)
            self.assertTrue(os.path.isfile(name + ".old"), name)

        self._remove_files(self._rename_suffixes)


class RegressionTest_CompareRestart(unittest.TestCase):
    """Tests to verify comparison of restart files.

    """

    def setUp(self):
        self.testlog = open("dummy.testlog", 'w')
        self.status = TestStatus()
        self.rt = RegressionTest()
        # test_data is the test section from a config file
        test_data = {"name": "dummy_test_name", }
        # timeout as specified by the command line arg
        timeout = None
        # check performance as specified by the command line arg
        check_performance = False
        # criteria is the default-test-criteria section from the config file
        criteria = {}
        self.rt.setup(criteria, test_data, timeout,
                      check_performance, self.testlog)

    def tearDown(self):
        self.testlog.close()

    def test_get_hash_missing_file(self):
        """Mark a test fail if restart file is missing.
        """
        tmp_filename = "tmp.file"
        file_hash = self.rt._get_binary_restart_hash(tmp_filename, self.status, self.testlog)
        self.assertEqual(self.status.fail, 1)

    def test_get_hash(self):
        """Get the hash of a file.
        """
        tmp_filename = "tmp.file"
        with open(tmp_filename, 'w') as tmpfile:
            tmpfile.write(str("This is a file."))
        file_hash = self.rt._get_binary_restart_hash(tmp_filename, self.status, self.testlog)
        self.assertEqual(file_hash, "679a60d9c581d48aa428a323879edc40a01e2a9d")
        os.remove(tmp_filename)

    def test_compare_hashes_same(self):
        """status.fail not set when hashing indentical files.
        """
        data = struct.pack('<HHHH',
                           0x0000, 0xFFFF, 0x0000, 0xFFFF)
        tmp_filename_1 = "{0}-restart.chk".format(self.rt.name())
        with open(tmp_filename_1, 'wb') as tmpfile_1:
            tmpfile_1.write(data)

        tmp_filename_2 = "tmp-restart-{0}-restart.chk".format(self.rt.name())
        with open(tmp_filename_2, 'wb') as tmpfile_2:
            tmpfile_2.write(data)

        self.rt._check_restart(self.status, self.testlog)
        self.assertEqual(self.status.fail, 0)

        os.remove(tmp_filename_1)
        os.remove(tmp_filename_2)

    def test_compare_hashes_different(self):
        """Fail when hashing files that differ by a single bit.
        """
        data = struct.pack('<HHHH',
                            0x0000, 0xFFFF, 0x0000, 0xFFFF)
        tmp_filename_1 = "{0}-restart.chk".format(self.rt.name())
        with open(tmp_filename_1, 'wb') as tmpfile_1:
            tmpfile_1.write(data)

        data = struct.pack('<HHHH',
                            0x0000, 0xFFFF, 0x0010, 0xFFFF)
        tmp_filename_2 = "tmp-restart-{0}-restart.chk".format(self.rt.name())
        with open(tmp_filename_2, 'wb') as tmpfile_2:
            tmpfile_2.write(data)

        self.rt._check_restart(self.status, self.testlog)
        self.assertEqual(self.status.fail, 1)

        os.remove(tmp_filename_1)
        os.remove(tmp_filename_2)

class RegressionTest_CompareHDF5(unittest.TestCase):
    """Tests to verify comparison of meta-data in hdf5 files.

    """

    def setUp(self):
        self.testlog = open("dummy.testlog", 'w')
        self.status = TestStatus()
        self.rt = RegressionTest()
        # test_data is the test section from a config file
        test_data = {"name": "dummy_test_name", }
        # timeout as specified by the command line arg
        timeout = None
        # check performance as specified by the command line arg
        check_performance = False
        # criteria is the default-test-criteria section from the config file
        criteria = {}
        self.rt.setup(criteria, test_data, timeout,
                      check_performance, self.testlog)
        self.h5_gold = h5py.File("dummy.h5.gold", 'w')
        self.h5_gold.create_group("Provenance")
        self.h5_gold.create_dataset("/Provenance/Bikes", shape=(2,))
        self.h5_gold.create_group("Cat")
        self.h5_gold.create_dataset("/Cat/Paws", shape=(2, 3))
        self.h5_gold.create_dataset("/Cat/Whiskers", shape=(2, 4, 6))
        self.h5_gold.create_group("Dog")
        self.h5_gold.create_dataset("/Dog/Tail", shape=(3,4))
        self.h5_gold.close()
        self.h5_gold = h5py.File("dummy.h5.gold", 'r')

    def tearDown(self):
        self.testlog.close()
        self.h5_gold.close()

    def test_different_num_groups(self):
        """Mark a test fail if the number of groups is different
        """
        h5_current = h5py.File("dummy.h5", 'w')
        h5_current.create_group("Provenance")
        h5_current.create_group("Cat")
        h5_current.create_dataset("/Cat/Paws", shape=(2, 3))
        h5_current.create_dataset("/Cat/Whiskers", shape=(2, 4, 6))
        h5_current.close()
        h5_current = h5py.File("dummy.h5", 'r')

        self.rt._compare_hdf5_data(h5_current, self.h5_gold, self.status, self.testlog)
        self.assertEqual(self.status.fail, 1)

    def test_different_num_datasets(self):
        """Mark a test fail if the number of datasets in a group is different
        """
        h5_current = h5py.File("dummy.h5", 'w')
        h5_current.create_group("Provenance")
        h5_current.create_group("Cat")
        h5_current.create_dataset("/Cat/Paws", shape=(2, 3))
        h5_current.create_dataset("/Cat/Whiskers", shape=(2, 4, 6))
        h5_current.create_dataset("/Cat/Ears", shape=(2, 4, 3))
        h5_current.create_group("Dog")
        h5_current.create_dataset("/Dog/Tail", shape=(3,4))
        h5_current.close()
        h5_current = h5py.File("dummy.h5", 'r')

        self.rt._compare_hdf5_data(h5_current, self.h5_gold, self.status, self.testlog)
        self.assertEqual(self.status.fail, 1)
        h5_current.close()

    def test_different_dataset_names(self):
        """Mark a test fail if the datasets in a group are different
        """
        h5_current = h5py.File("dummy.h5", 'w')
        h5_current.create_group("Provenance")
        h5_current.create_group("Cat")
        h5_current.create_dataset("/Cat/Paws", shape=(2, 3))
        h5_current.create_dataset("/Cat/Ears", shape=(2, 4, 6))
        h5_current.create_group("Dog")
        h5_current.create_dataset("/Dog/Tail", shape=(3,4))
        h5_current.close()
        h5_current = h5py.File("dummy.h5", 'r')

        self.rt._compare_hdf5_data(h5_current, self.h5_gold, self.status, self.testlog)
        self.assertEqual(self.status.fail, 1)
        h5_current.close()

    def test_different_dataset_shapes(self):
        """Mark a test fail if the dataset shapes are different
        """
        h5_current = h5py.File("dummy.h5", 'w')
        h5_current.create_group("Provenance")
        h5_current.create_group("Cat")
        h5_current.create_dataset("/Cat/Paws", shape=(2, 3))
        h5_current.create_dataset("/Cat/Whiskers", shape=(2, 4, ))
        h5_current.create_group("Dog")
        h5_current.create_dataset("/Dog/Tail", shape=(3,4))
        h5_current.close()
        h5_current = h5py.File("dummy.h5", 'r')

        self.rt._compare_hdf5_data(h5_current, self.h5_gold, self.status, self.testlog)
        self.assertEqual(self.status.fail, 1)
        h5_current.close()

    def test_different_dataset_dtype(self):
        """Mark a test fail if the dataset dtypes are different
        """
        h5_current = h5py.File("dummy.h5", 'w')
        h5_current.create_group("Provenance")
        h5_current.create_group("Cat")
        h5_current.create_dataset("/Cat/Paws", shape=(2, 3))
        h5_current.create_dataset("/Cat/Whiskers", shape=(2, 4, 6), dtype='i8')
        h5_current.create_group("Dog")
        h5_current.create_dataset("/Dog/Tail", shape=(3,4))
        h5_current.close()
        h5_current = h5py.File("dummy.h5", 'r')

        self.rt._compare_hdf5_data(h5_current, self.h5_gold, self.status, self.testlog)
        self.assertEqual(self.status.fail, 1)
        h5_current.close()

    def test_same_hdf5(self):
        """Mark a test success if the hdf5 files have the same groups,
        datasets and metadata.

        """
        h5_current = h5py.File("dummy.h5.gold", 'r')

        self.rt._compare_hdf5_data(h5_current, self.h5_gold, self.status, self.testlog)
        self.assertEqual(self.status.fail, 0)
        h5_current.close()


if __name__ == '__main__':
    #unittest.main(buffer=True)
    unittest.main()
