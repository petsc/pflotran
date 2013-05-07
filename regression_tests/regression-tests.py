#!/bin/env python
#
# Program to manage and run PFloTran regression tests
#
#
#

from __future__ import print_function
from __future__ import division

import argparse
from collections import deque
import datetime
import math
import os
import pprint
import re
import shutil
import subprocess
import sys
import textwrap
import time
import traceback

if sys.version_info[0] == 2:
    import ConfigParser as config_parser
else:
    import configparser as config_parser


class TestStatus(object):
    def __init__(self):
        self.fail = 0
        self.warning = 0
        self.error = 0
        self.skipped = 0
        self.test_count = 0

    def __str__(self):
        message = "fail = {0}\n".format(self.fail)
        message += "warning = {0}\n".format(self.warning)
        message += "error = {0}\n".format(self.error)
        message += "skipped = {0}\n".format(self.skipped)
        message += "test_count = {0}\n".format(self.test_count)
        return message


class RegressionTest(object):
    """
    Class to collect data about a test problem, run the problem, and
    compare the results to a known result.
    """
    _wall_time_re = re.compile("Time \(seconds\)")

    def __init__(self):
        # define some constants
        self._ABSOLUTE = "absolute"
        self._RELATIVE = "relative"
        self._PERCENT = "percent"
        self._TIME = "time"
        self._CONCENTRATION = "concentration"
        self._GENERIC = "generic"
        self._DISCRETE = "discrete"
        self._RATE = "rate"
        self._VOLUME_FRACTION = "volume_fraction"
        self._PRESSURE = "pressure"
        self._SATURATION = "saturation"
        self._SOLUTION = "solution"
        self._RESIDUAL = "residual"
        self._TOL_VALUE = 0
        self._TOL_TYPE = 1
        self._PFLOTRAN_SUCCESS = 86
        # misc test parameters
        self._pprint = pprint.PrettyPrinter(indent=2)
        self._txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")
        self._debug = False
        self._executable = None
        self._input_arg = "-input_prefix"
        self._input_suffix = "in"
        self._np = None
        self._timeout = 60.0
        self._check_performance = False
        self._num_failed = 0
        self._test_name = None
        # assign default tolerances for different classes of variables
        self._tolerance = {}
        self._tolerance[self._TIME] = [5.0, self._PERCENT]
        self._tolerance[self._CONCENTRATION] = [1.0e-12, self._ABSOLUTE]
        self._tolerance[self._GENERIC] = [1.0e-12, self._ABSOLUTE]
        self._tolerance[self._DISCRETE] = [0, self._ABSOLUTE]
        self._tolerance[self._RATE] = [1.0e-12, self._ABSOLUTE]
        self._tolerance[self._VOLUME_FRACTION] = [1.0e-12, self._ABSOLUTE]
        self._tolerance[self._PRESSURE] = [1.0e-12, self._ABSOLUTE]
        self._tolerance[self._SATURATION] = [1.0e-12, self._ABSOLUTE]
        self._tolerance[self._RESIDUAL] = [1.0e-12, self._ABSOLUTE]

    def __str__(self):
        message = "  {0} :\n".format(self.name())
        message += "    timeout = {0}\n".format(self._timeout)
        message += "    np = {0}\n".format(self._np)
        message += "    executable args :\n"
        message += "        input arg : {0}\n".format(self._input_arg)
        message += "        input suffix : {0}\n".format(self._input_suffix)
        message += "    test criteria :\n"
        for k in self._tolerance:
            message += "        {0} : {1} [{2}]\n".format(
                k,
                self._tolerance[k][self._TOL_VALUE],
                self._tolerance[k][self._TOL_TYPE])

        return message

    def setup(self, executable_args, default_criteria, test_data,
              timeout, check_performance):
        self._test_name = test_data["name"]

        if executable_args is not None:
            self._set_executable_args(executable_args)

        self._set_test_data(default_criteria, test_data,
                            timeout, check_performance)

    def name(self):
        return self._test_name

    def run(self, mpiexec, executable, dry_run, status, testlog):
        """
        Build up the run command, including mpiexec, np, pflotran,
        input file, output file. Then run the job as a subprocess.

        * NOTE(bja) - starting in python 3.3, we can use:

          subprocess.Popen(...).wait(timeout)

          to catch hanging jobs, but for python < 3.3 we have to
          manually manage the timeout...?
        """
        # TODO(bja) : need to think more about the desired behavior if
        # mpiexec is passed for a serial test or not passed for a
        # parallel test.
        command = []
        if self._np is not None:
            if mpiexec:
                command.append(mpiexec)
                command.append("-np")
                command.append(self._np)
            else:
                # parallel test, but don't have mpiexec, we mark the
                # test as skipped and bail....
                message = self._txtwrap.fill(
                    "WARNING : mpiexec was not provided for a parallel test '{0}'.\n"
                    "This test was skipped!".format(self.name()))
                print(message, file=testlog)
                status.skipped = 1
                return None

        command.append(executable)
        #geh: kludge for -malloc 0
        command.append("-malloc")
        command.append("0")
        if self._input_arg is not None:
            command.append(self._input_arg)
            command.append(self.name())

        if os.path.isfile(self.name() + ".regression"):
            os.rename(self.name() + ".regression",
                      self.name() + ".regression.old")

        if os.path.isfile(self.name() + ".out"):
            os.rename(self.name() + ".out",
                      self.name() + ".out.old")

        if os.path.isfile(self.name() + ".stdout"):
            os.rename(self.name() + ".stdout",
                      self.name() + ".stdout.old")

        if not dry_run:
            print("    cd {0}".format(os.getcwd()), file=testlog)
            print("    {0}".format(" ".join(command)), file=testlog)
            run_stdout = open(self.name() + ".stdout", 'w')
            start = time.time()
            proc = subprocess.Popen(command,
                                    shell=False,
                                    stdout=run_stdout,
                                    stderr=subprocess.STDOUT)
            while proc.poll() is None:
                time.sleep(0.1)
                if time.time() - start > self._timeout:
                    proc.kill()
                    time.sleep(0.1)
                    message = self._txtwrap.fill(
                        "ERROR: job '{0}' has exceeded timeout limit of "
                        "{1} seconds.".format(self.name(), self._timeout))
                    print(''.join(['\n', message, '\n']), file=testlog)
            pflotran_status = abs(proc.returncode)
            run_stdout.close()
        # pflotran returns 0 on an error (e.g. can't find an input
        # file), 86 on success. 59 for timeout errors?
        if pflotran_status != self._PFLOTRAN_SUCCESS:
            message = self._txtwrap.fill(
                "FAIL : {name} : pflotran return an error "
                "code ({status}) indicating the simulation may have "
                "failed. Please check '{name}.out' and '{name}.stdout' "
                "for error messages.".format(
                    name=self.name(), status=pflotran_status))
            print("".join(['\n', message, '\n']), file=testlog)
            status.fail = 1

    def check(self, status, testlog):
        """
        Test the output from the run against the known "gold standard"
        output and determine if the test succeeded or failed.

        We return zero on success, one on failure so that the test
        manager can track how many tests succeeded and failed.
        """
        gold_filename = self.name() + ".regression.gold"
        if not os.path.isfile(gold_filename):
            message = self._txtwrap.fill(
                "FAIL: could not find regression test gold file "
                "'{0}'. If this is a new test, please create "
                "it with '--new-test'.".format(gold_filename))
            print("".join(['\n', message, '\n']), file=testlog)
            status.fail = 1
            return
        else:
            with open(gold_filename, 'rU') as gold_file:
                gold_output = gold_file.readlines()

        current_filename = self.name() + ".regression"
        if not os.path.isfile(current_filename):
            message = self._txtwrap.fill(
                "FAIL: could not find regression test file '{0}'."
                " Please check the standard output file for "
                "errors.".format(current_filename))
            print("".join(['\n', message, '\n']), file=testlog)
            status.fail = 1
            return
        else:
            with open(current_filename, 'rU') as current_file:
                current_output = current_file.readlines()

        print("    diff {0} {1}".format(gold_filename, current_filename), file=testlog)

        gold_sections = self._get_sections(gold_output)
        current_sections = self._get_sections(current_output)
        if self._debug:
            print("--- Gold sections:")
            self._pprint.pprint(gold_sections)
            print("--- Current sections:")
            self._pprint.pprint(current_sections)

        # look for sections that are in gold but not current
        for s in gold_sections:
            if s not in current_sections:
                self._num_failed += 1
                print("    FAIL: section '{0}' is in the gold output, but "
                      "not the current output.".format(s), file=testlog)

        # look for sections that are in current but not gold
        for s in current_sections:
            if s not in gold_sections:
                self._num_failed += 1
                print("    FAIL: section '{0}' is in the current output, "
                      "but not the gold output.".format(s), file=testlog)

        # compare common sections
        for s in gold_sections:
            if s in current_sections:
                self._num_failed += self._compare_sections(gold_sections[s],
                                                           current_sections[s], testlog)

        if self._num_failed > 0:
            status.fail = 1

    def update(self, status, testlog):
        """
        Update the gold standard test results to the current
        output. Both the current regression output and a gold file
        must exist.
        """
        gold_name = self.name() + ".regression.gold"
        current_name = self.name() + ".regression"

        # verify that the gold file exists
        if not os.path.isfile(gold_name):
            print("ERROR: test '{0}' results can not be updated "
                  "because a gold file does not "
                  "exist!".format(self.name()), file=testlog)
            status.error = 1

        # verify that the regression file exists
        if not os.path.isfile(current_name):
            print("ERROR: test '{0}' results can not be updated "
                  "because no regression file "
                  "exists!".format(self.name()), file=testlog)
            status.error = 1
        try:
            print("  updating test '{0}'... ".format(self.name()),
                  end='', file=testlog)
            os.rename(current_name, gold_name)
            print("done", file=testlog)
        except Exception as e:
            status = 1
            message = str(e)
            message += "\nFAIL : Could not rename '{0}' to '{1}'. "
            message += "Please rename the file manually!".format(current_name,
                                                                 gold_name)
            message += "    mv {0} {1}".format(current_name, gold_name)
            print(message, file=testlog)
            status.fail = 1

    def new_test(self, status, testlog):
        """
        A new test does not have a gold standard regression test. We
        will check to see if a gold standard file exists (an error),
        then create the gold file by copying the current regression
        file to gold.
        """
        gold_name = self.name() + ".regression.gold"
        current_name = self.name() + ".regression"

        # check if the gold file exists already
#        if os.path.isfile(gold_name):
#            raise Exception("ERROR: test '{0}' was classified as new, "
#                            "but a gold file already "
#                            "exists!".format(self.name()))

        # check that the regression file was created.
        if not os.path.isfile(current_name):
            print("ERROR: could not create new gold file for "
                  "test '{0}' because no regression file "
                  "exists!".format(self.name()), file=testlog)
            status.error = 1

        try:
            print("  creating gold file '{0}'... ".format(self.name()),
                  end='', file=testlog)

            os.rename(current_name, gold_name)
            print("done", file=testlog)
        except Exception as e:
            status = 1
            message = str(e)
            message += "\nFAIL : Could not rename '{0}' to '{1}'. "
            message += "Please rename the file manually!".format(current_name,
                                                                 gold_name)
            message += "    mv {0} {1}".format(current_name, gold_name)
            print(message, file=testlog)
            status.fail = 1

    def _get_sections(self, output):
        """
        Each section in the regression test file looks like:
        -- TYPE: NAME --
           key: value
           key: value
        """
        name_re = re.compile("^--[\s]+([\w]+):(.*)[\s]+--$")
        sections = {}
        s = {}
        for line in output:
            match = name_re.match(line)
            if match:
                # save the old section, if any
                if 'name' in s:
                    sections[s['name']] = s
                name = match.group(2).strip()
                data_type = match.group(1)
                s = {}
                s['name'] = name
                s['type'] = data_type
            else:
                temp = line.split(':')
                name = temp[0].strip()
                value = temp[1].strip()
                s[name] = value
        # add the final section
        if 'name' in s:
            sections[s['name']] = s

        return sections

    def _compare_sections(self, gold_section, current_section, testlog):
        name = gold_section['name']
        data_type = gold_section['type']
        section_status = 0
        if self._check_performance is False and data_type.lower() == self._SOLUTION:
            # solution blocks contain platform dependent performance
            # metrics. We skip them unless they are explicitly
            # requested.
            print("    Skipping {0} : {1}".format(data_type, name), file=testlog)
        else:
            # if key in gold but not in current --> failed test
            for k in gold_section:
                if k not in current_section:
                    section_status += 1
                    print("    FAIL: key '{0}' in section '{1}' found in gold "
                          "output but not current".format(
                              k, gold_section['name']), file=testlog)

            # if key in current but not gold --> failed test
            for k in current_section:
                if k not in gold_section:
                    section_status += 1
                    print("    FAIL: key '{0}' in section '{1}' found in current "
                          "output but not gold".format(k, current_section['name']),
                          file=testlog)

            # now compare the keys that are in both...
            for k in gold_section:
                if k == "name" or k == 'type':
                    pass
                elif k in current_section:
                    name_str = name + ":" + k
                    # the data may be vector
                    gold = gold_section[k].split()
                    current = current_section[k].split()
                    if len(gold) != len(current):
                        section_status += 1
                        print("    FAIL: {0} : {1} : vector lengths not "
                              "equal. gold {2}, current {3}".format(
                                  name, k, len(gold), len(current)), file=testlog)
                    else:
                        for i in range(len(gold)):
                            try:
                                status = self._compare_values(
                                    name_str, data_type, gold[i], current[i], testlog)
                                section_status += status
                            except Exception as e:
                                section_status += 1
                                print("ERROR: {0} : {1}.\n  {2}".format(
                                    self.name(), k, str(e)), file=testlog)

        if False:
            print("    {0} : status : {1}".format(
                name, section_status), file=testlog)
        return section_status

    def _compare_values(self, name, key, previous, current, testlog):
        """
        NOTE(bja): previous and current come into this function as
        strings. We don't know if they should be floats or ints (or
        possibly strings) until we know what the data type is. For
        'discrete' or 'solution' variables, we have to do further work
        to figure it out!
        """
        status = 0
        tolerance_type = None
        tolerance = None
        key = key.lower()
        if (key == self._CONCENTRATION or
            key == self._GENERIC or
            key == self._RATE or
            key == self._VOLUME_FRACTION or
            key == self._PRESSURE or
                key == self._SATURATION):
            previous = float(previous)
            current = float(current)
            tolerance_type = self._tolerance[key][self._TOL_TYPE]
            tolerance = self._tolerance[key][self._TOL_VALUE]
        elif key.lower() == self._SOLUTION:
            previous, current, tolerance_type, tolerance = \
                self._compare_solution(name, previous, current)
        elif key.lower() == self._DISCRETE:
            previous, current, tolerance_type, tolerance = \
                self._compare_discrete(name, previous, current)
        else:
            # should be an error? We'll fail anyway in the checks
            # because nothing is set.
            print("WARNING: the data caterogy '{0}' for '{1}' is not a known "
                  "data category.".format(key, name), file=testlog)

        if tolerance_type == self._ABSOLUTE:
            delta = abs(previous - current)
        elif (tolerance_type == self._RELATIVE or
              tolerance_type == self._PERCENT):
            if previous != 0:
                delta = abs(previous - current) / previous
            elif current != 0:
                delta = abs(previous - current) / current
            else:
                # both are zero
                delta = 0.0
            if tolerance_type == self._PERCENT:
                delta *= 100.0
        else:
            raise Exception("ERROR: unknown test tolerance_type '{0}' for "
                            "variable '{1}, {2}.'".format(tolerance_type,
                                                          name, key))
        if delta > tolerance:
            status = 1
            print("    FAIL: {0} : {1} > {2} [{3}]".format(
                name, delta, tolerance,
                tolerance_type), file=testlog)
        elif self._debug:
            print("    PASS: {0} : {1} <= {2} [{3}]".format(
                name, delta, tolerance,
                tolerance_type))

        return status

    def _compare_solution(self, name, previous, current):
        # NOTE(bja): hard coding this for now until we decide how do
        # deal with all the different requirements.
        section = name.split(':')[0]
        param = name.split(':')[1]
        param = param.strip()
        if param == "Time (seconds)":
            previous = float(previous)
            current = float(current)
            tolerance = self._tolerance[self._TIME][self._TOL_VALUE]
            tolerance_type = self._tolerance[self._TIME][self._TOL_TYPE]
        elif param == "Time Steps":
            previous = int(previous)
            current = int(current)
            tolerance = self._tolerance[self._DISCRETE][self._TOL_VALUE]
            tolerance_type = self._tolerance[self._DISCRETE][self._TOL_TYPE]
        elif param == "Newton Iterations":
            previous = int(previous)
            current = int(current)
            tolerance = self._tolerance[self._DISCRETE][self._TOL_VALUE]
            tolerance_type = self._tolerance[self._DISCRETE][self._TOL_TYPE]
        elif param == "Solver Iterations":
            previous = int(previous)
            current = int(current)
            tolerance = self._tolerance[self._DISCRETE][self._TOL_VALUE]
            tolerance_type = self._tolerance[self._DISCRETE][self._TOL_TYPE]
        elif param == "Time Step Cuts":
            previous = int(previous)
            current = int(current)
            tolerance = self._tolerance[self._DISCRETE][self._TOL_VALUE]
            tolerance_type = self._tolerance[self._DISCRETE][self._TOL_TYPE]
        elif param == "Solution 2-Norm":
            previous = float(previous)
            current = float(current)
            tolerance = self._tolerance[self._GENERIC][self._TOL_VALUE]
            tolerance_type = self._tolerance[self._GENERIC][self._TOL_TYPE]
        elif param == "Residual 2-Norm":
            previous = float(previous)
            current = float(current)
            tolerance = self._tolerance[self._RESIDUAL][self._TOL_VALUE]
            tolerance_type = self._tolerance[self._RESIDUAL][self._TOL_TYPE]
        else:
            raise Exception("ERROR: unknown variable '{0}' in solution "
                            "section '{1}'".format(param, section))

        return previous, current, tolerance_type, tolerance

    def _compare_discrete(self, name, previous, current):
        # NOTE(bja): discrete values are integers, except when we are
        # looking at the mean of a discrete variable. Then we
        # may(probably) have a floating point value!

        mean_re = re.compile("Mean")
        have_mean = mean_re.search(name)
        if not have_mean:
            previous = int(previous)
            current = int(current)
            tolerance = self._tolerance[self._DISCRETE][self._TOL_VALUE]
            tolerance_type = self._tolerance[self._DISCRETE][self._TOL_TYPE]
        else:
            previous = float(previous)
            current = float(current)
            tolerance = self._tolerance[self._GENERIC][self._TOL_VALUE]
            tolerance_type = self._tolerance[self._GENERIC][self._TOL_TYPE]

        return previous, current, tolerance_type, tolerance

    def _set_executable_args(self, executable_args):
        if "input arg" in executable_args:
            self._input_arg = executable_args["input arg"]

        if "input suffix" in executable_args:
            self._input_suffix = executable_args["input suffix"]

    def _set_test_data(self, default_criteria, test_data, timeout, check_performance):
        """
        Set the test criteria for different categories of variables.
        """
        self._np = test_data.pop('np', None)

        self._check_performance = check_performance

        # timeout : preference (1) command line (2) test data (3) class default
        self._timeout = float(test_data.pop('timeout', self._timeout))
        if timeout:
            self._timeout = float(timeout[0])

        self._set_criteria(self._TIME, default_criteria, test_data)

        self._set_criteria(self._CONCENTRATION, default_criteria, test_data)

        self._set_criteria(self._GENERIC, default_criteria, test_data)

        self._set_criteria(self._DISCRETE, default_criteria, test_data)

        self._set_criteria(self._RATE, default_criteria, test_data)

        self._set_criteria(self._VOLUME_FRACTION, default_criteria, test_data)

        self._set_criteria(self._PRESSURE, default_criteria, test_data)

        self._set_criteria(self._SATURATION, default_criteria, test_data)

    def _set_criteria(self, key, default_criteria, test_data):
        """
        Our prefered order for selecting test criteria is:
        (1) test data section of the config file
        (2) config-file wide defaults
        (3) hard coded class default
        """
        if key in test_data:
            self._tolerance[key] = \
                self._validate_criteria(key, test_data[key])
        elif key in default_criteria:
            self._tolerance[key] = \
                self._validate_criteria(key, default_criteria[key])
        elif key in self._tolerance:
            # already a correctly formatted list and stored
            pass
        else:
            raise Exception("ERROR : tolerance for data type '{0}' could "
                            "not be determined from the config file or "
                            "default values!".format(key))

    def _validate_criteria(self, key, criteria):
        value = criteria.split()[0]
        try:
            value = float(value)
        except Exception as e:
            raise Exception("ERROR : Could not convert '{0}' test criteria "
                            "value '{1}' into a float!".format(key, value))

        criteria_type = criteria.split()[1]
        if (criteria_type.lower() != self._PERCENT and
            criteria_type.lower() != self._ABSOLUTE and
                criteria_type.lower() != self._RELATIVE):
            raise Exception("ERROR : invalid test criteria string '{0}' "
                            "for '{1}'".format(criteria_type, key))
        return [value, criteria_type]


class RegressionTestManager(object):
    """
    Class to open a configuration file, process it into a group of
    tests, and manange running the tests.

    Notes:

        * The ConfigParser class converts all section names and key
          names into lower case. This means we need to preprocess user
          input names to lower case.

        * Foo
    """

    NO_TESTS_RUN = -1000

    def __init__(self):
        self._debug = False
        self._file_status = TestStatus()
        self._config_filename = None
        self._executable_args = None
        self._default_test_criteria = None
        self._available_tests = {}
        self._available_suites = {}
        self._tests = []
        self._txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")
        self._pprint = pprint.PrettyPrinter(indent=2)

    def __str__(self):
        data = "Regression Test Manager :\n"
        data += "    configuration file : {0}\n".format(self._config_filename)
        data += "    executable_args :\n"
        data += self._dict_to_string(self._executable_args)
        data += "    default test criteria :\n"
        data += self._dict_to_string(self._default_test_criteria)
        data += "    suites :\n"
        data += self._dict_to_string(self._available_suites)
        data += "    available tests :\n"
        data += self._dict_to_string(self._available_tests)

        data += "Tests :\n"
        for t in self._tests:
            data += t.__str__()

        return data

    def num_tests(self):
        return len(self._tests)

    def generate_tests(self, config_file, user_suites, user_tests,
                       timeout, check_performance, testlog):
        self._read_config_file(config_file)
        self._validate_suites()
        user_suites, user_tests = self._validate_user_lists(user_suites,
                                                            user_tests, testlog)
        self._create_tests(user_suites, user_tests, timeout, check_performance)

    def run_tests(self, mpiexec, executable,
                  dry_run, update, new_test, check_only, testlog):
        """
        Run the tests specified in the config file.

        * dry_run - flag indicates that the test is setup then print
          the command that would be used, but don't actually run
          anything or compare results.

        * new_test - flag indicates that the test is a new test, and
          there should not be a gold standard regression file
          present. Run the executable and create the gold file.

        * update - flag indicates that the output from pflotran has
          changed, and we want to update the gold standard regression
          file to reflect this. Run the executable and replace the
          gold file.

        * check_only - flag to indicate just diffing the existing
          regression files without rerunning pflotran.
        """
        if self.num_tests() > 0:
            if new_test:
                self._run_new(mpiexec, executable, dry_run, testlog)
            elif update:
                self._run_update(mpiexec, executable, dry_run, testlog)
            elif check_only:
                self._check_only(dry_run, testlog)
            else:
                self._run_check(mpiexec, executable, dry_run, testlog)
        else:
            self._file_status.test_count = 0

    def _run_check(self, mpiexec, executable, dry_run, testlog):
        if dry_run:
            print("Dry run:")
        print("Running tests from '{0}':".format(self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for t in self._tests:
            status = TestStatus()
            self._test_header(t.name(), testlog)

            t.run(mpiexec, executable, dry_run, status, testlog)

            if not dry_run and status.skipped == 0:
                t.check(status, testlog)

            self._add_to_file_status(status)

            self._test_summary(t.name(), status, dry_run,
                               "passed", "failed", testlog)

        self._print_file_summary(dry_run, "passed", "failed", testlog)

    def _check_only(self, dry_run, testlog):
        if dry_run:
            print("Dry run:")
        print("Checking existing test results from '{0}':".format(
            self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for t in self._tests:
            status = TestStatus()
            self._test_header(t.name(), testlog)

            if not dry_run and status.skipped == 0:
                t.check(status, testlog)

            self._add_to_file_status(status)

            self._test_summary(t.name(), status, dry_run,
                               "passed", "failed", testlog)

        self._print_file_summary(dry_run, "passed", "failed", testlog)

    def _run_new(self, mpiexec, executable, dry_run, testlog):
        if dry_run:
            print("Dry run:")

        print("New tests from '{0}':".format(self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for t in self._tests:
            status = TestStatus()
            self._test_header(t.name(), testlog)

            t.run(mpiexec, executable, dry_run, status, testlog)

            if not dry_run and status.skipped == 0:
                t.new_test(status, testlog)
            self._add_to_file_status(status)
            self._test_summary(t.name(), status, dry_run,
                               "created", "error creating new test files.", testlog)

        self._print_file_summary(dry_run, "created", "could not be created", testlog)

    def _run_update(self, mpiexec, executable, dry_run, testlog):
        if dry_run:
            print("Dry run:")
        print("Updating tests from '{0}':".format(self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for t in self._tests:
            status = TestStatus()
            self._test_header(t.name(), testlog)
            t.run(mpiexec, executable, dry_run, status, testlog)

            if not dry_run and status.skipped == 0:
                t.update(status, testlog)
            self._add_to_file_status(status)
            self._test_summary(t.name(), status, dry_run,
                               "updated", "error updating test.", testlog)

        self._print_file_summary(dry_run, "updated", "could not be updated", testlog)

    def _test_header(self, name, testlog):
        print(40 * '-', file=testlog)
        print("{0}... ".format(name), file=testlog)

    def _test_summary(self, name, status, dry_run,
                      success_message, fail_message, testlog):
        if dry_run:
            print("D", end='', file=sys.stdout)
            print(" dry run.", file=testlog)
        else:
            if (status.fail == 0 and
                    status.warning == 0 and
                    status.error == 0 and
                    status.skipped == 0):
                print(".", end='', file=sys.stdout)
                print("{0}... {1}.".format(name, success_message), file=testlog)
            elif status.fail != 0:
                print("F", end='', file=sys.stdout)
                print("{0}... {1}.".format(name, fail_message), file=testlog)
            elif status.warning != 0:
                print("W", end='', file=sys.stdout)
            elif status.error != 0:
                print("E", end='', file=sys.stdout)
            elif status.skipped != 0:
                print("S", end='', file=sys.stdout)
                print("{0}... skipped.".format(name), file=testlog)
            else:
                print("?", end='', file=sys.stdout)

        sys.stdout.flush()

    def _print_file_summary(self, dry_run, success_message, fail_message, testlog):
        # print a summary of the results for this config file
        print(50 * '-', file=testlog)
        if dry_run:
            print("{0} : dry run.".format(self._config_filename), file=testlog)
        elif self._file_status.test_count == 0:
            print("{0} : no tests run.".format(self._config_filename), file=testlog)
        else:
            line = "{0} : {1} tests : ".format(self._config_filename,
                                               self._file_status.test_count)
            if self._file_status.fail > 0:
                line = "{0} {1} tests {2}, ".format(
                    line, self._file_status.fail, fail_message)
            if self._file_status.skipped > 0:
                line = "{0} {1} tests {2}, ".format(
                    line, self._file_status.skipped, "skipped")
            num_passed = (self._file_status.test_count -
                          self._file_status.fail - self._file_status.skipped)
            line = "{0} {1} tests {2}".format(line, num_passed, success_message)
            print(line, file=testlog)

    def _add_to_file_status(self, status):
        self._file_status.fail += status.fail
        self._file_status.warning += status.warning
        self._file_status.error += status.error
        self._file_status.skipped += status.skipped
        self._file_status.test_count += 1

    def run_status(self):
        return self._file_status

    def display_available_tests(self):
        print("Available tests: ")
        for t in sorted(self._available_tests.keys()):
            print("    {0}".format(t))

    def display_available_suites(self):
        print("Available test suites: ")
        for s in self._available_suites:
            print("    {0} :".format(s))
            for t in self._available_suites[s].split():
                print("        {0}".format(t))

    def _read_config_file(self, config_file):
        """
        Read the configuration file.

        Sections : The config file will have know sections:
        "executable", "suites", "default-test-criteria".

        All other sections are assumed to be test names.
        """
        if config_file is None:
            raise Exception("Error, must provide a config filename")
        self._config_filename = config_file
        config = config_parser.SafeConfigParser()
        config.read(self._config_filename)

        if config.has_section("executable"):
            self._executable_args = \
                self._list_to_dict(config.items("executable"))

        if config.has_section("default-test-criteria"):
            self._default_test_criteria = \
                self._list_to_dict(config.items("default-test-criteria"))

        if config.has_section("suites"):
            self._available_suites = \
                self._list_to_dict(config.items("suites"))

        self._identify_tests(config)

    def _identify_tests(self, config):
        # section names are test names
        test_names = config.sections()

        # remove the fixed section names
        if config.has_section("executable"):
            test_names.remove("executable")
        if config.has_section("default-test-criteria"):
            test_names.remove("default-test-criteria")
        if config.has_section("suites"):
            test_names.remove("suites")

        # all remaining sections should be individual tests
        for t in test_names:
            self._available_tests[t] = self._list_to_dict(config.items(t))
            self._available_tests[t]['name'] = t

    def _dict_to_string(self, data):
        temp = ""
        for k, v in data.items():
            temp += "        {0} : {1}\n".format(k, v)
        return temp

    def _list_to_dict(self, input_list):
        output_dict = {}
        for item in input_list:
            output_dict[item[0]] = item[1]
        return output_dict

    def _validate_suites(self):
        """
        Validates the suites defined in configuration file by
        checking that each test in a suite is one of the available
        tests.

        If the config file has an empty suite, we report that to the
        user then remove it from the list.
        """
        invalid_tests = []
        empty_suites = []
        for s in self._available_suites:
            suite_tests = self._available_suites[s].split()
            if len(suite_tests) == 0:
                empty_suites.append(s)
            else:
                # validate the list
                for t in suite_tests:
                    if t not in self._available_tests:
                        name = "suite : '{0}' --> test : '{1}'".format(s, t)
                        invalid_tests.append(name)

        for s in empty_suites:
            # empty suite, warn the user and remove it from the list
            del self._available_suites[s]
            print("DEV WARNING : {0} : cfg validation : empty suite "
                  ": '{1}'".format(self._config_filename, s))

        if len(invalid_tests) != 0:
            raise Exception("ERROR : suites contain unknown tests in "
                            "configuration file '{0}' : {1}".format(
                                self._config_filename, invalid_tests))

    def _validate_user_lists(self, user_suites, user_tests, testlog):
        """
        Check that the list of suites or tests passed from the command
        line are valid.
        """
        # if no suites or tests is specified, use all available tests
        if len(user_suites) == 0 and len(user_tests) == 0:
            u_suites = []
            u_tests = self._available_tests
        else:
            # check that the processed user supplied names are valid
            # convert user supplied names to lower case
            u_suites = []
            for s in user_suites:
                if s.lower() in self._available_suites:
                    u_suites.append(s.lower())
                else:
                    message = self._txtwrap.fill(
                        "WARNING : {0} : Skipping requested suite '{1}' (not "
                        "present, misspelled or empty).".format(
                            self._config_filename, s))
                    print(message, file=testlog)

            u_tests = []
            for t in user_tests:
                if t in self._available_tests:
                    u_tests.append(t.lower())
                else:
                    message = self._txtwrap.fill(
                        "WARNING : {0} : Skipping test '{1}' (not present or "
                        "misspelled).".format(self._config_filename, t))
                    print(message, file=testlog)

        return u_suites, u_tests

    def _create_tests(self, user_suites, user_tests, timeout, check_performance):
        all_tests = user_tests
        for s in user_suites:
            for t in self._available_suites[s].split():
                all_tests.append(t)

        for t in all_tests:
            try:
                test = RegressionTest()
                test.setup(self._executable_args, self._default_test_criteria,
                           self._available_tests[t], timeout, check_performance)
                self._tests.append(test)
            except Exception as e:
                raise Exception("ERROR : could not create test '{0}' from "
                                "config file '{1}'. {2}".format(
                                    t, self._config_filename, str(e)))


def commandline_options():
    parser = argparse.ArgumentParser(description='Run a pflotran regression '
                                     'tests or suite of tests.')

    parser.add_argument('--backtrace', action='store_true',
                        help='show exception backtraces as extra debugging '
                        'output')

    parser.add_argument('--advanced', action='store_true',
                        help="enable advanced options for developers")

    parser.add_argument('-c', '--config-files', nargs="+", default=None,
                        help='test configuration file to use')

    parser.add_argument('--check-only', action='store_true', default=False,
                        help="diff the existing regression files without "
                        "running pflotran again.")

    parser.add_argument('--check-performance', action='store_true', default=False,
                        help="include the performance metrics ('SOLUTION' blocks) "
                        "in regression checks.")

    parser.add_argument('--debug', action='store_true',
                        help='extra debugging output')

    parser.add_argument('-d', '--dry-run',
                        default=False, action='store_true',
                        help='perform a dry run, setup the test commands but '
                        'don\'t run them')

    parser.add_argument('-e', '--executable', nargs=1, default=['executable'],
                        help='path to executable to use for testing')

    parser.add_argument('--list-suites', default=False, action='store_true',
                        help='print the list of test suites from the config '
                        'file and exit')

    parser.add_argument('--list-tests', default=False, action='store_true',
                        help='print the list of tests from the config file '
                        'and exit')

    parser.add_argument('-m', '--mpiexec', nargs=1, default=None,
                        help="path to the executable for mpiexec (mpirun, etc)"
                        "on the current machine.")

    parser.add_argument('-n', '--new-tests',
                        action="store_true", default=False,
                        help="indicate that there are new tests being run. "
                        "Skips the output check and creates a new gold file.")

    parser.add_argument('-r', '--recursive-search', nargs='*', default=None,
                        help='recursively search the current directory and '
                        'all sub-directories, using any configuration files '
                        'in those directories.')

    parser.add_argument('-s', '--suites', nargs="+", default=[],
                        help='space separated list of test suite names')

    parser.add_argument('-t', '--tests', nargs="+", default=[],
                        help='space separated list of test names')

    parser.add_argument('--timeout', nargs=1, default=None,
                        help="test timeout (for assuming a job has hung and "
                        "needs to be killed)")

    parser.add_argument('-u', '--update',
                        action="store_true", default=False,
                        help='update the tests listed by the "--tests" '
                        'option, with the current output becoming the new '
                        'gold standard')

    options = parser.parse_args()
    return options


def generate_config_file_list(options):
    """
    Try to generate a list of configuration files from the commandline
    options.
    """
    config_file_list = []
    # search for config files
    if options.recursive_search is not None:
        if options.recursive_search == []:
            # if we have an empty list, use the cwd as the starting point
            options.recursive_search.append(os.getcwd())
        for base_dir in options.recursive_search:
            if not os.path.isabs(base_dir):
                # if we received a relative path, make it absolute
                base_dir = os.path.abspath(base_dir)

            if os.path.isdir(base_dir):
                search_for_config_files(base_dir, config_file_list)
            else:
                raise Exception("ERROR: can not search for config files "
                                "in '{0}' because it is not a "
                                "directory.".format(base_dir))

    # add the explicitly listed config files
    if options.config_files is not None:
        for f in options.config_files:
            if not os.path.isabs(f):
                f = os.path.abspath(f)
            if os.path.isfile(f):
                config_file_list.append(f)
            else:
                raise Exception("ERROR: specified config file '{0}' is not a "
                                "file!".format(f))

    if options.debug:
        print("\nKnown config files:")
        for c in config_file_list:
            print("    {0}".format(c))

    if len(config_file_list) == 0:
        raise Exception("ERROR: no config files were found. Please specify a "
                        "config file with '--config' or search for files "
                        "with '--recursive-search'.")

    return config_file_list


def search_for_config_files(base_dir, config_file_list):
    """
    recursively search the directory tree, creating a list of all config files
    """
    subdirlist = []
    for entry in os.listdir(base_dir):
        # append entry to path
        file_path = os.path.join(base_dir, entry)
         # determine whether a file or a directory
        if os.path.isfile(file_path):
            # is the file a config file?
            if file_path.endswith('.cfg'):
                # if it's a config file, add it to the list
                config_file_list.append(file_path)
        else:
            # add path to list of subdirectories
            subdirlist.append(file_path)

    # recursive search the subdirectories
    for sub_dir in subdirlist:
        search_for_config_files(sub_dir, config_file_list)


def check_options(options):
    """
    Run some sanity checks on the commandline options.
    """
    # prevent the user from updating regression output during a
    # recursive search for config files
    if options.update and options.recursive_search is not None:
        if not options.advanced:
            raise Exception("ERROR: can not update gold regression files "
                            "during a recursive search for config files.")

    if options.update and options.new_tests:
        raise Exception("ERROR: can not create new tests and update gold "
                        "regression files at the same time.")


def check_for_executable(options):
    """
    Try to verify that we have something reasonable for the executable
    """
    # check the executable
    if options.executable == ['executable']:
        options.dry_run = True
        executable = "/usr/bin/false"
    else:
        # absolute path to the executable
        executable = os.path.abspath(options.executable[0])
        # is it a valid file?
        if not os.path.isfile(executable):
            raise Exception("ERROR: executable is not a valid file: "
                            "'{0}'".format(executable))
    return executable


def check_for_mpiexec(options, testlog):
    """
    Try to verify that we have something reasonable for the mpiexec executable

    Notes:

    geh: need to add code to determine full path of mpiexec if not specified

    bja: the problem is that we don't know how to get the correct
    mpiexec. On the mac, there is a system mpiexe that shows up in the
    path, but this is provided by apple and it is not the correct one
    to use because it doesn't include a fortran compiler. We need the
    exact mpiexec/mpirun that was used to compile petsc, which may
    come from a system installed package in /usr/bin or /opt/local/bin
    or maybe petsc compiled it. This is best handled outside the test
    manager, e.g. use make to identify mpiexec from the petsc
    variables
    """

    # check for mpiexec
    mpiexec = None
    if options.mpiexec is not None:
        # mpiexec = os.path.abspath(options.mpiexec[0])
        mpiexec = options.mpiexec[0]
        # try to log some info about mpiexec
        print("MPI information :", file=testlog)
        print("-----------------", file=testlog)
        tempfile = "{0}/tmp-pflotran-regression-test-info.txt".format(os.getcwd())
        command = [mpiexec, "--version"]
        append_command_to_log(command, testlog, tempfile)
        print("\n\n", file=testlog)
        os.remove(tempfile)
        # is it a valid file?
###       if not os.path.isfile(mpiexec):
###           raise Exception("ERROR: mpiexec is not a valid file: "
###                           "'{0}'".format(mpiexec))
    else:
        message = ("\n** WARNING ** : mpiexec was not provided on the command line.\n"
                   "                All parallel tests will be skipped!\n")
        print(message, file=sys.stdout)
        print(message, file=testlog)

    return mpiexec


def summary_report_by_file(report, outfile):
    print(70 * '-', file=outfile)
    print("Regression test file summary:", file=outfile)
    for t in report:
        line = "    {0}... {1} tests : ".format(t, report[t].test_count)
        if report[t].warning > 0:
            line = "{0} {1} test warnings, ".format(line, report[t].warning)
        if report[t].error > 0:
            line = "{0} {1} test errors, ".format(line, report[t].error)

        if report[t].test_count == 0:
            line = "{0}... no tests were run.".format(line)
        else:
            if report[t].fail > 0:
                line = "{0} {1} tests failed, ".format(line, report[t].fail)
            if report[t].skipped > 0:
                line = "{0} {1} tests skipped, ".format(line, report[t].skipped)
            if report[t].fail == 0 and report[t].skipped == 0:
                line = "{0} all tests passed".format(line)
            else:
                num_passed = (report[t].test_count - report[t].fail -
                              report[t].skipped)
                line = "{0} {1} passed.".format(line, num_passed)

        print("{0}".format(line), file=outfile)

    print("\n", file=outfile)


def summary_report(run_time, report, outfile):
    print(70 * '-', file=outfile)
    print("Regression test summary:", file=outfile)
    print("    Total run time: {0:4g} [s]".format(run_time), file=outfile)
    test_count = 0
    num_failures = 0
    num_errors = 0
    num_warnings = 0
    num_skipped = 0
    for t in report:
        test_count += report[t].test_count
        num_failures += report[t].fail
        num_errors += report[t].error
        num_warnings += report[t].warning
        num_skipped += report[t].skipped

    print("    Total tests : {0}".format(test_count), file=outfile)

    if num_skipped > 0:
        print("    Skipped : {0}".format(num_skipped), file=outfile)
        success = False

    print("    Tests run : {0}".format(test_count - num_skipped), file=outfile)

    success = True
    if num_failures > 0:
        print("    Failed : {0}".format(num_failures), file=outfile)
        success = False

    if num_errors > 0:
        print("    Errors : {0}".format(num_errors), file=outfile)
        success = False

    if num_warnings > 0:
        print("    Warnings : {0}".format(num_warnings), file=outfile)
        success = False

    if success:
        print("    All tests passed.", file=outfile)

    print("\n", file=outfile)
    return num_failures


def append_command_to_log(command, testlog, tempfile):
    print("$ {0}".format(" ".join(command)), file=testlog)
    testlog.flush()
    with open(tempfile, "w") as tempinfo:
        proc = subprocess.Popen(command, shell=False,
                                stdout=tempinfo,
                                stderr=subprocess.STDOUT)
        time.sleep(0.5)
    shutil.copyfileobj(open(tempfile,'r'), testlog)


def setup_testlog(txtwrap):
    now = datetime.datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
    filename = "pflotran-tests-{0}.testlog".format(now)
    testlog = open(filename, 'w')
    print("  Test log file : {0}".format(filename))

    # try to report some useful information about the environment,
    # petsc, pflotran...
    print("PFloTran Regression Test Log", file=testlog)
    print("Date : {0}".format(now), file=testlog)
    print("System Info :", file=testlog)
    print("    platform : {0}".format(sys.platform), file=testlog)
    test_dir = os.getcwd()
    print("Test directory : ", file=testlog)
    print("    {0}".format(test_dir), file=testlog)
    
    tempfile = "{0}/tmp-pflotran-regression-test-info.txt".format(test_dir)

    print("\nPFLOTRAN repository status :", file=testlog)
    print("----------------------------", file=testlog)
    if os.path.isdir("{0}/../.hg".format(test_dir)):
        cmd = ["hg", "parent"]
        append_command_to_log(cmd, testlog, tempfile)
        cmd = ["hg", "status", "-q"]
        append_command_to_log(cmd, testlog, tempfile)
        print("\n\n", file=testlog)
    else:
        print("    unknown", file=testlog)

    print("PETSc information :", file=testlog)
    print("-------------------", file=testlog)
    petsc_dir = os.getenv("PETSC_DIR", None)
    if petsc_dir:
        message = txtwrap.fill(
            "* WARNING * This information may be incorrect if you have more "
            "than one version of petsc installed.\n")
        print(message, file=testlog)
        print("    PETSC_DIR : {0}".format(petsc_dir), file=testlog)
        petsc_arch = os.getenv("PETSC_ARCH", None)
        if petsc_arch:
            print("    PETSC_ARCH : {0}".format(petsc_arch), file=testlog)
        
        os.chdir(petsc_dir)
        print("    petsc repository status :", file=testlog)
        if os.path.isdir("{0}/.git".format(petsc_dir)):
            cmd = ["git", "log", "-1", "HEAD"]
            append_command_to_log(cmd, testlog, tempfile)
            cmd = ["git", "status", "-u", "no"]
            append_command_to_log(cmd, testlog, tempfile)
        elif os.path.isdir("{0}/.hg".format(petsc_dir)):
            cmd = ["hg", "parent"]
            append_command_to_log(cmd, testlog, tempfile)
            cmd = ["hg", "status", "-q"]
            append_command_to_log(cmd, testlog, tempfile)
        else:
            print("    No git or hg directory was found in your PETSC_DIR",
                  file=testlog)
        os.chdir(test_dir)
        print("\n\n", file=testlog)
    else:
        print("    PETSC_DIR was not defined.", file=testlog)

    os.remove(tempfile)
    return testlog


def main(options):
    txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")
    testlog = setup_testlog(txtwrap)

    check_options(options)
    executable = check_for_executable(options)
    mpiexec = check_for_mpiexec(options, testlog)
    config_file_list = generate_config_file_list(options)

    print("Running pflotran regression tests :")

    # loop through config files, cd into the appropriate directory,
    # read the appropriate config file and run the various tests.
    start = time.time()
    report = {}
    for f in config_file_list:
        try:
            # NOTE(bja): the try block is inside this loop so that if
            # a single test throws an exception in a large batch of
            # tests, we can recover and at least try running the other
            # config files.
            print(80 * '=', file=testlog)

            # get the absolute path of the directory
            test_dir = os.path.dirname(f)
            # cd into the test directory so that the relative paths in
            # test files are correct
            os.chdir(test_dir)
            if options.debug:
                print("Changed to working directory: {0}".format(test_dir))

            test_manager = RegressionTestManager()

            if options.debug:
                test_manager._debug = True

            # get the relative file name
            filename = os.path.basename(f)

            test_manager.generate_tests(filename,
                                        options.suites,
                                        options.tests,
                                        options.timeout,
                                        options.check_performance,
                                        testlog)

            if options.debug:
                print(70 * '-')
                print(test_manager)

            if options.list_suites:
                test_manager.display_available_suites()

            if options.list_tests:
                test_manager.display_available_tests()

            test_manager.run_tests(mpiexec,
                                   executable,
                                   options.dry_run,
                                   options.update,
                                   options.new_tests,
                                   options.check_only,
                                   testlog)

            report[filename] = test_manager.run_status()
        except Exception as e:
            message = txtwrap.fill(
                "ERROR: a problem occured in file '{0}'.  This is "
                "probably an error with commandline options, the "
                "configuration file, or an internal error.  The "
                "error is:\n{1}".format(f, str(e)))
            print(''.join(['\n', message, '\n']), file=testlog)
            if options.backtrace:
                traceback.print_exc()
            print('F', end='', file=sys.stdout)
            report[filename] = TestStatus()
            report[filename].fail = 1

    stop = time.time()
    status = 0
    if not options.dry_run and not options.update:
        print("")
        run_time = stop - start
        summary_report_by_file(report, testlog)
        summary_report(run_time, report, testlog)
        status = summary_report(run_time, report, sys.stdout)

    if options.update:
        message = txtwrap.fill(
            "Test results were updated! Please document why you modified the "
            "gold standard test results in your revision control commit message!\n")
        print(''.join(['\n', message, '\n']))

    testlog.close()

    return status

if __name__ == "__main__":
    options = commandline_options()
    try:
        status = main(options)
        sys.exit(status)
    except Exception as e:
        print(str(e))
        if options.backtrace:
            traceback.print_exc()
        sys.exit(1)
