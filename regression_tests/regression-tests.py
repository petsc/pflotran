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
import math
import os
import pprint
import re
import subprocess
import sys
import time
import traceback

if sys.version_info[0] == 2:
    import ConfigParser as config_parser
else:
    import configparser as config_parser


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
        self._debug = False
        self._verbose = False
        self._executable = None
        self._input_arg = "-pflotranin"
        self._input_suffix = "in"
        self._output_arg = "-output_prefix"
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
        message += "        output arg : {0}\n".format(self._output_arg)
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

    def run(self, mpiexec, executable, dry_run, verbose):
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
        if mpiexec:
            command.append(mpiexec)
            command.append("-np")
            if self._np is None:
                self._np = '1'
                if verbose:
                    print("WARNING : mpiexec specified for test '{0}', "
                          "but the test section does not specify the number "
                          "of parallel jobs! Running test as "
                          "serial.".format(self.name()))
            command.append(self._np)
        else:
            if self._np is not None:
                raise Exception("ERROR : test '{0}' : np was specified in "
                                "the test data, but mpiexec was not "
                                "provided.".format(self.name()))

        command.append(executable)
        input_file_name = self.name() + '.' + self._input_suffix
        if self._input_arg != None:
            command.append(self._input_arg)
            command.append(input_file_name)
        if self._output_arg != None:
            command.append(self._output_arg)
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

        status = -1
        if dry_run:
            print("\n    {0}".format(" ".join(command)))
        else:
            if verbose:
                print("    {0}".format(" ".join(command)))
            run_stdout = open(self.name() + ".stdout", 'w')
            start = time.time()
            proc = subprocess.Popen(command,
                                    shell=False,
                                    stdout=run_stdout,
                                    stderr=subprocess.STDOUT)
            while proc.poll() is None:
                time.sleep(0.1)
                if time.time() - start > self._timeout:
                    proc.terminate()
                    time.sleep(0.1)
                    print("ERROR: job '{0}' has exceeded timeout limit of "
                          "{1} seconds.".format(self.name(), self._timeout))
            status = abs(proc.returncode)
            run_stdout.close()
        # pflotran returns 0 on an error (e.g. can't find an input
        # file), 86 on success. 59 for timeout errors?
        if status != self._PFLOTRAN_SUCCESS:
            print("\nWARNING : {name} : pflotran return an error "
                  "code ({status}) indicating the simulation may have "
                  "failed. Please check '{name}.out' and '{name}.stdout' "
                  "for error messages.\n".format(
                    name=self.name(), status=status))
        return status

    def check(self, verbose):
        """
        Test the output from the run against the known "gold standard"
        output and determine if the test succeeded or failed.

        We return zero on success, one on failure so that the test
        manager can track how many tests succeeded and failed.
        """
        self._verbose = verbose
        gold_filename = self.name() + ".regression.gold"
        if not os.path.isfile(gold_filename):
            print("ERROR: could not find regression test gold file "
                  "'{0}'. If this is a new test, please create "
                  "it with '--new-test'.".format(gold_filename))
            return 1
        else:
            with open(gold_filename, 'rU') as gold_file:
                gold_output = gold_file.readlines()

        current_filename = self.name() + ".regression"
        if not os.path.isfile(current_filename):
            print("ERROR: could not find regression test file '{0}'."
                  " Please check the standard output file for "
                  "errors.".format(current_filename))
            return 1
        else:
            with open(current_filename, 'rU') as current_file:
                current_output = current_file.readlines()

        if verbose:
            print("    diff {0} {1}".format(gold_filename, current_filename))

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
                if self._verbose:
                    print("    FAIL: section '{0}' is in the gold output, but "
                          "not the current output.".format(s))

        # look for sections that are in current but not gold
        for s in current_sections:
            if s not in gold_sections:
                self._num_failed += 1
                if self._verbose:
                    print("    FAIL: section '{0}' is in the current output, "
                          "but not the gold output.".format(s))

        # compare common sections
        for s in gold_sections:
            if s in current_sections:
                self._num_failed += self._compare_sections(gold_sections[s],
                                                           current_sections[s])

        status = 0
        if self._num_failed > 0:
            status = 1
        return status

    def update(self, verbose):
        """
        Update the gold standard test results to the current
        output. Both the current regression output and a gold file
        must exist.
        """
        status = 0
        gold_name = self.name() + ".regression.gold"
        current_name = self.name() + ".regression"

        # verify that the gold file exists
        if not os.path.isfile(gold_name):
            raise Exception("ERROR: test '{0}' results can not be updated "
                            "because a gold file does not "
                            "exist!".format(self.name()))

        # verify that the regression file exists
        if not os.path.isfile(current_name):
            raise Exception("ERROR: test '{0}' results can not be updated "
                            "because no regression file "
                            "exists!".format(self.name()))

        try:
            if verbose:
                print("  updating test '{0}'... ".format(self.name()),
                      end='')
            os.rename(current_name, gold_name)
            if verbose:
                print("done")
        except Exception as e:
            status = 1
            message = str(e)
            message += "\nERROR : Could not rename '{0}' to '{1}'. "
            message += "Please rename the file manually!".format(current_name,
                                                                 gold_name)
            message += "    mv {0} {1}".format(current_name, gold_name)
            print(message)
            # should we rethrow this exception, or continue?
            #raise Exception(message)
        return status

    def new_test(self, verbose):
        """
        A new test does not have a gold standard regression test. We
        will check to see if a gold standard file exists (an error),
        then create the gold file by copying the current regression
        file to gold.
        """
        status = 0
        gold_name = self.name() + ".regression.gold"
        current_name = self.name() + ".regression"

        # check if the gold file exists already
#        if os.path.isfile(gold_name):
#            raise Exception("ERROR: test '{0}' was classified as new, "
#                            "but a gold file already "
#                            "exists!".format(self.name()))

        # check that the regression file was created.
        if not os.path.isfile(current_name):
            raise Exception("ERROR: could not create new gold file for "
                            "test '{0}' because no regression file "
                            "exists!".format(self.name()))

        try:
            if verbose:
                print("  creating gold file '{0}'... ".format(self.name()),
                      end='')

            os.rename(current_name, gold_name)
            if verbose:
                print("done")
        except Exception as e:
            status = 1
            message = str(e)
            message += "\nERROR : Could not rename '{0}' to '{1}'. "
            message += "Please rename the file manually!".format(current_name,
                                                                 gold_name)
            message += "    mv {0} {1}".format(current_name, gold_name)
            print(message)
            # should we rethrow this exception, or continue?
            #raise Exception(message)
        return status

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

    def _compare_sections(self, gold_section, current_section):
        name = gold_section['name']
        data_type = gold_section['type']
        section_status = 0
        if self._check_performance == False and data_type.lower() == self._SOLUTION:
            # solution blocks contain platform dependent performance
            # metrics. We skip them unless they are explicitly
            # requested.
            if self._verbose:
                print("    Skipping {0} : {1}".format(data_type, name))
        else:
            # if key in gold but not in current --> failed test
            for k in gold_section:
                if k not in current_section:
                    section_status += 1
                    if self._verbose:
                        print("    FAIL: key '{0}' in section '{1}' found in gold "
                              "output but not current".format(
                                k, gold_section['name']))

            # if key in current but not gold --> failed test
            for k in current_section:
                if k not in gold_section:
                    section_status += 1
                    print("    FAIL: key '{0}' in section '{1}' found in current "
                          "output but not gold".format(k, current_section['name']))

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
                        if self._verbose:
                            print("    FAIL: {0} : {1} : vector lengths not "
                                  "equal. gold {2}, current {3}".format(
                                    name, k, len(gold), len(current)))
                    else:
                        for i in range(len(gold)):
                            try:
                                status = self._compare_values(name_str, data_type,
                                                              gold[i], current[i])
                                section_status += status
                            except Exception as e:
                                section_status += 1
                                if self._verbose:
                                    print("ERROR: {0} : {1}.\n  {2}".format(
                                            self.name(), k, str(e)))


        if self._verbose and False:
            print("    {0} : status : {1}".format(name, section_status))
        return section_status

    def _compare_values(self, name, key, previous, current):
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
                  "data category.".format(key, name))

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
            if self._verbose:
                print("    FAIL: {0} : {1} > {2} [{3}]".format(
                        name, delta, tolerance,
                        tolerance_type))
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

        if "output arg" in executable_args:
            self._output_arg = executable_args["output arg"]

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

    def __init__(self):
        self._debug = False
        self._verbose = False
        self._num_failed = 0
        self._config_filename = None
        self._executable_args = None
        self._default_test_criteria = None
        self._available_tests = {}
        self._available_suites = None
        self._tests = []

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

    def generate_tests(self, config_file, user_suites, user_tests,
                       timeout, check_performance):
        self._read_config_file(config_file)
        self._validate_suites()
        user_suites, user_tests = self._validate_user_lists(user_suites,
                                                            user_tests)
        self._create_tests(user_suites, user_tests, timeout, check_performance)

    def run_tests(self, mpiexec, executable, verbose,
                  dry_run, update, new_test, check_only):
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

        if new_test:
            self._run_new(mpiexec, executable, dry_run, verbose)
        elif update:
            self._run_update(mpiexec, executable, dry_run, verbose)
        elif check_only:
            self._check_only(dry_run, verbose)
        else:
            self._run_check(mpiexec, executable, dry_run, verbose)

    def _run_check(self, mpiexec, executable, dry_run, verbose):
        if dry_run:
            print("Dry run:")
        else:
            print("Running tests from '{0}':".format(self._config_filename))
        print(50 * '-')

        for t in self._tests:
            self._test_header(t.name(), verbose)

            t.run(mpiexec, executable, dry_run, verbose)

            status = 0
            if not dry_run:
                status = t.check(verbose)

            self._num_failed += status

            self._test_summary(t.name(), status, verbose, dry_run,
                               "passed", "failed")

        self._print_file_summary(dry_run, "passed", "failed")

    def _check_only(self, dry_run, verbose):
        if dry_run:
            print("Dry run:")
        else:
            print("Diffing tests from '{0}':".format(self._config_filename))
        print(50 * '-')

        for t in self._tests:
            self._test_header(t.name(), verbose)

            status = 0
            if not dry_run:
                status = t.check(verbose)

            self._num_failed += status

            self._test_summary(t.name(), status, verbose, dry_run,
                               "passed", "failed")

        self._print_file_summary(dry_run, "passed", "failed")

    def _run_new(self, mpiexec, executable, dry_run, verbose):
        if dry_run:
            print("Dry run:")
        else:
            print("New tests from '{0}':".format(self._config_filename))
        print(50 * '-')

        for t in self._tests:
            self._test_header(t.name(), verbose)

            t.run(mpiexec, executable, dry_run, verbose)

            status = 0
            if not dry_run:
                status = t.new_test(verbose)
            self._num_failed += status
            self._test_summary(t.name(), status, verbose, dry_run,
                               "created", "error creating new test files.")

        self._print_file_summary(dry_run, "created", "could not be created")

    def _run_update(self, mpiexec, executable, dry_run, verbose):
        if dry_run:
            print("Dry run:")
        else:
            print("Updating tests from '{0}':".format(self._config_filename))
        print(50 * '-')

        for t in self._tests:
            self._test_header(t.name(), verbose)
            t.run(mpiexec, executable, dry_run, verbose)

            status = 0
            if not dry_run:
                status = t.update(verbose)
            self._num_failed += status
            self._test_summary(t.name(), status, verbose, dry_run,
                               "updated", "error updating test.")

        self._print_file_summary(dry_run, "updated", "could not be updated")

    def _test_header(self, name, verbose):
        if verbose:
            print(40 * '-')
        print("{0}... ".format(name), end='')
        if verbose:
            print()

    def _test_summary(self, name, status, verbose, dry_run,
                      success_message, fail_message):
        if status == 0:
            if not dry_run:
                if verbose:
                    print("{0}... {1}.".format(name, success_message))
                else:
                    print(" {0}.".format(success_message))
            else:
                print(" skipped.")
        else:
            if verbose:
                print("{0}... {1}.".format(name, fail_message))
            else:
                print(" {0}.".format(fail_message))

    def _print_file_summary(self, dry_run, success_message, fail_message):
        # print a summary of the results for this config file
        print(50 * '-')
        if self._num_failed > 0:
            print("{0} : {1} of {2} tests {3}".format(
                    self._config_filename, self._num_failed, len(self._tests),
                    fail_message))
        else:
            if not dry_run:
                print("{0} : {1} tests {2}".format(self._config_filename,
                                                   len(self._tests),
                                                   success_message))
            else:
                print("{0} : no tests run.".format(self._config_filename))

    def status(self):
        return self._num_failed

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
        if config_file == None:
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
        invalid_tests = []
        for s in self._available_suites:
            suite_tests = self._available_suites[s].split()
            for t in suite_tests:
                if t not in self._available_tests:
                    name = "suite : '{0}' --> test : '{1}'".format(s, t)
                    invalid_tests.append(name)

        if len(invalid_tests) != 0:
            raise Exception("ERROR : suites contain unknown tests in "
                            "configuration file '{0}' : {1}".format(
                    self._config_filename, invalid_tests))

    def _validate_user_lists(self, user_suites, user_tests):

        # if no suites or tests is specified, use all tests
        if len(user_suites) == 0 and len(user_tests) == 0:
            u_suites = []
            u_tests = self._available_tests
        else:
            # convert user supplied names to lower case
            u_suites = []
            for s in user_suites:
                u_suites.append(s.lower())

            u_tests = []
            for t in user_tests:
                u_tests.append(t.lower())

        # now check that the processed user supplied names are valid
        invalid_user_names = []
        for s in u_suites:
            if s.lower() not in self._available_suites:
                invalid_user_names.append("suite : '{0}'".format(s))

        for t in u_tests:
            if t not in self._available_tests:
                invalid_user_names.append("test : '{0}'".format(t))

        if len(invalid_user_names) != 0:
            # exception or print a warning and continue...?
            raise Exception("ERROR : {0} : unknown suite or test provided "
                            "on command line : {1}".format(
                    self._config_filename, invalid_user_names))

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

    parser.add_argument('-c', '--config-file', nargs=1, default=None,
                        help='test configuration file to use')

    parser.add_argument('--check-only', action='store_true', default=False,
                        help="diff the existing regression files without "
                        "running pflotran again.")

    parser.add_argument('--check-performance', action='store_true', default=False,
                        help="include the performance metrics ('SOLUTION' blocks)"
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

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose output')

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
    if options.config_file is not None:
        for f in options.config_file:
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


def check_for_mpiexec(options):
    """
    Try to verify that we have something reasonable for the mpiexec executable
    """
    # check for mpiexec
    mpiexec = None
    if options.mpiexec is not None:
        # absolute path to the executable
#geh: need to add code to determine full path of mpiexec if not specified
#        mpiexec = os.path.abspath(options.mpiexec[0])
        mpiexec = options.mpiexec[0]
        # is it a valid file?
#        if not os.path.isfile(mpiexec):
#            raise Exception("ERROR: mpiexec is not a valid file: "
#                            "'{0}'".format(mpiexec))
    return mpiexec


def main(options):
    check_options(options)
    executable = check_for_executable(options)
    mpiexec = check_for_mpiexec(options)
    config_file_list = generate_config_file_list(options)

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
            print(70 * '-')

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
                                        options.check_performance)

            if options.debug:
                print(70 * '-')
                print(test_manager)

            if options.list_suites:
                test_manager.display_available_suites()

            if options.list_tests:
                test_manager.display_available_tests()

            test_manager.run_tests(mpiexec,
                                   executable,
                                   options.verbose,
                                   options.dry_run,
                                   options.update,
                                   options.new_tests,
                                   options.check_only)

            report[filename] = test_manager.status()
        except Exception as e:
            message = ("\nERROR: a problem occured in file '{0}'.\n  This is "
                       "probably an error with commandline options, the "
                       "configuration file, or an enternal error.\n  The "
                       "error is:\n{1}".format(f, str(e)))
            print(message)
            if options.backtrace:
                traceback.print_exc()
            report[filename] = -1

    stop = time.time()
    status = 0
    if not options.dry_run and not options.update:
        print(70 * '-')
        print("Regression test summary:")
        print("    Total run time: {0:4g} [s]".format(stop - start))
        for t in report:
            status += report[t]
            if report[t] > 0:
                print("    {0}... {1} tests failed".format(t, report[t]))
            elif report[t] == 0:
                print("    {0}... all tests passed".format(t))
            else:
                print("    {0}... could not be run.".format(t, report[t]))
        print("\n\n")

    if options.update:
        print("\nTest results were updated!\n"
              "Please document why you modified the gold standard test "
              "results in your revision control commit message!\n")

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
