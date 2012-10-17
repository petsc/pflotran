#!/bin/env python
#
# Program to manage and run PFloTran regression tests
#
#
#

import argparse
from collections import deque

import difflib
import math
import os
import pprint
import re
import subprocess
import sys
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

    def __init__(self):
        self._pprint = pprint.PrettyPrinter(indent=2)
        self._debug = None
        self._verbose = None
        self._executable = None
        self._input_arg = None
        self._input_suffix = None
        self._output_arg = None
        self._test_name = None
        self._time_tolerance = None
        self._time_type = None
        self._concentration_tolerance = None
        self._concentration_type = None
        self._generic_tolerance = None
        self._generic_type = None
        self._discrete_tolerance = None
        self._discrete_type = None
        self._num_failed = 0


    def __str__(self):
        message = "    {0} :\n".format(self._test_name)
        message += "        test criteria :\n"
        message += "            time : {0} [{1}]\n".format(self._time_tolerance, self._time_type)
        message += "            concentration : {0} [{1}]\n".format(self._concentration_tolerance, self._concentration_type)
        message += "            generic : {0} [{1}]\n".format(self._generic_tolerance, self._generic_type)
        message += "            discrete : {0} [{1}]\n".format(self._discrete_tolerance, self._discrete_type)
        message += "        executable args :\n"
        message += "            input arg : {0}\n".format(self._input_arg)
        message += "            input suffix : {0}\n".format(self._input_suffix)
        message += "            output arg : {0}\n".format(self._output_arg)


        return message

    def setup(self, executable_args, default_criteria, test_data):
        self._test_name = test_data["name"]

        if executable_args is not None:
            self._set_executable_args(executable_args)

        self._set_test_criteria(default_criteria, test_data)

    def run(self, executable, dry_run, verbose):
        # need some os specific magic here...?
        input_file_name = self._test_name + '.' + self._input_suffix
        command = "{0} ".format(executable)
        if self._input_arg != None:
            command += "{0} {1} ".format(self._input_arg, input_file_name)
        if self._output_arg != None:
            command += "{0} {1} ".format(self._output_arg, self._test_name)
        if dry_run:
            print("    dry run: \"{0}\"".format(command))
        else:
            if verbose:
                print("    {0}".format(command))
            run_stdout = open(self._test_name + ".stdout", 'w')
            subprocess.Popen(command.split(), 
                             stdout=run_stdout, stderr=subprocess.STDOUT).wait()
            run_stdout.close()

    def diff(self, verbose):
        """
        Test the output from the run against the known "gold standard"
        output and determine if the test succeeded or failed.

        We return one on success, zero on failure so that the test
        manager can track how many tests succeeded and failed.
        """
        self._verbose = verbose
        gold_filename = self._test_name + ".regression.gold"
        with open(gold_filename, 'rU') as gold_file:
            gold_output = gold_file.readlines()

        current_filename = self._test_name + ".regression"
        with open(current_filename, 'rU') as current_file:
            current_output = current_file.readlines()

        if verbose:
            print("    diff {0} {1}".format(current_filename, gold_filename))

        gold_sections = self._get_sections(gold_output)
        current_sections = self._get_sections(current_output)
        for s in gold_sections.keys():
            # compare sections
            if s in current_sections.keys():
                self._num_failed += self._compare_sections(gold_sections[s], current_sections[s])
            else:
                # what is the correct way to handle this...?
                print("ERROR: section '{0}' is in the gold output, but not the current output!".format(s))
        status = 0    
        if self._num_failed > 0:
            status = 1
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
        if self._debug:
            self._pprint.pprint(sections)
        return sections

    def _compare_sections(self, gold_section, current_section):
        name = gold_section['name']
        data_type = gold_section['type']
        section_status = 0
        for k in gold_section.keys():
            if k == "name" or k == 'type':
                break
            if k in current_section.keys():
                gold = gold_section[k]
                current = current_section[k]
                name_str = name + " --> " + k 
                status = self._compare_values(name_str, data_type, gold, current)
                section_status += status
            else:
                # correct way to handle this?
                print("ERROR: key '{0}' in gold section '{1}' but not in current section".format(k, gold_section['name']))
        if self._verbose and False:
            print("    {0} : status : {1}".format(name, section_status))
        return section_status

    def _compare_values(self, name, data_type, previous, current):
        status = 0
        if data_type.lower() == "concentration":
            status = self._compare_concentration(name, float(previous), float(current))
        elif data_type.lower() == "solution":
            status = self._compare_solution(name, previous, current)
        elif data_type.lower() == "generic":
            status = self._compare_generic(name, float(previous), float(current))
        elif data_type.lower() == "discrete":
            status = self._compare_discrete(name, int(previous), int(previous))
        return status

    def _compare_concentration(self, name, previous, current):
        if self._concentration_type == "absolute":
            delta = math.fabs(previous - current)
        elif self._concentration_type == "relative":
            delta = math.fabs(previous - current) / previous
        elif self._concentration_type == "percent":
            delta = 100.0 * math.abs(previous - current) / previous

        status = 0
        if delta > self._concentration_tolerance:
            status = 1
            if self._verbose:
                print("    FAILED: {0} : {1} > {2}".format(name, delta, self._concentration_tolerance))
        elif self._debug:
            print("    SUCCESS: {0} : {1} < {2}".format(name, delta, self._concentration_tolerance))
            
        return status

    def _compare_solution(self, name, previous, current):
        return 0

    def _compare_generic(self, name, previous, current):
        return 0

    def _compare_discrete(self, name, previous, current):
        return 0

    def _set_executable_args(self, executable_args):
        if "input arg" in executable_args.keys():
            self._input_arg = executable_args["input arg"]

        if "input suffix" in executable_args.keys():
            self._input_suffix = executable_args["input suffix"]

        if "output arg" in executable_args.keys():
            self._output_arg = executable_args["output arg"]
        

    def _set_test_criteria(self, default_criteria, test_data):
        criteria = self._get_criteria("time", 
                                      default_criteria, test_data)
        self._time_type = criteria[1]
        self._time_tolerance = criteria[0]

        criteria = self._get_criteria("concentration",
                                      default_criteria, test_data)
        self._concentration_type = criteria[1]
        self._concentration_tolerance = criteria[0]

        criteria = self._get_criteria("generic",
                                      default_criteria, test_data)
        self._generic_type = criteria[1]
        self._generic_tolerance = criteria[0]

        criteria = self._get_criteria("discrete",
                                      default_criteria, test_data)
        self._discrete_type = criteria[1]
        self._discrete_tolerance = criteria[0]

    def _get_criteria(self, key, default_criteria, test_data):
        criteria = None
        if key in test_data.keys():
            criteria = test_data[key]
        elif key in default_criteria.keys():
            criteria = default_criteria[key]
        else:
            raise Exception("ERROR : tolerance for '{0}' must be specified in either the default-test-criteria or test section!".format(key))

        return self._validate_criteria(key, criteria)

    def _validate_criteria(self, key, criteria):
        value = criteria.split()[0]
        try:
            value = float(value)
        except Exception as e:
            message = "ERROR : Could not convert '{0}' test criteria value '{1}' into a float!".format(key, value)
            raise Exception(message)

        criteria_type = criteria.split()[1]
        if (criteria_type.lower() != "percent" and
            criteria_type.lower() != "absolute" and
            criteria_type.lower() != "relative"):
            raise Exception("ERROR : invalid test criteria string '{0}' for '{1}'".format(criteria_type, key))
        return (value, criteria_type)

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

    def generate_tests(self, config_file, user_suites, user_tests):
        self._read_config_file(config_file)
        self._validate_suites()
        user_suites, user_names = self._validate_user_lists(user_suites, user_tests)
        self._create_tests(user_suites, user_tests)

    def run_tests(self, executable, dry_run, verbose, update):
        print(70*"-")
        print("Running tests:")
        for t in self._tests:
            if verbose:
                print(40*'-')
            print("{0}...".format(t._test_name),)
            if verbose:
                print()
            t.run(executable, dry_run, verbose)
            if update:
              # geh: if we are updating, skip the diff
              continue
            status = t.diff(verbose)
            self._num_failed += status
            if status == 0:
                if verbose:
                    print("{0}... passed.".format(t._test_name))
                else:
                    print(" passed.")
            else:
                if verbose:
                    print("{0}... failed.".format(t._test_name))
                else:
                    print(" failed.")

        print(70*"-")
        if not update:
            # geh: if we are updating, skip reporting the number of failed tests
            print("{0} of {1} tests failed".format(self._num_failed, len(self._tests)))

    def update_test_results(self, user_tests):
        pass

    def status(self):
        return self._num_failed

    def display_available_tests(self):
        print("Available tests: ")
        for t in sorted(self._available_tests.keys()):
            print("    {0}".format(t))

    def display_available_suites(self):
        print("Available test suites: ")
        for s in self._available_suites.keys():
            print("    {0} :".format(s))
            for t in self._available_suites[s].split():
                print("        {0}".format(t))

    def _read_config_file(self, config_file):
        """
        Read the configuration file.

        Sections : The config file will have know sections:
        "executable", "suites", "test-criteria".

        All other sections are assumed to be test names.
        """
        if config_file == None:
            raise Exception("Error, must provide a config filename")
        self._config_filename = config_file
        config = config_parser.SafeConfigParser()
        config.read(self._config_filename)

        if config.has_section("executable"):
            self._executable_args = self._list_to_dict(config.items("executable"))

        if config.has_section("default-test-criteria"):
            self._default_test_criteria = self._list_to_dict(config.items("default-test-criteria"))

        if config.has_section("suites"):
            self._available_suites = self._list_to_dict(config.items("suites"))

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
        for k, v in data.iteritems():
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
            raise Exception("ERROR : suites contain unknown tests in configuration file '{0}' : {1}".format(self._config_filename, invalid_tests))

    def _validate_user_lists(self, user_suites, user_tests):

        #geh: if no suites or tests is specified, use all suites
        if len(user_suites) == 0 and len(user_tests) == 0:
          u_suites = self._available_suites
          u_tests = []
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
            if s.lower() not in self._available_suites.keys():
                invalid_user_names.append("suite : '{0}'".format(s))

        for t in u_tests:
            if t not in self._available_tests.keys():
                invalid_user_names.append("test : '{0}'".format(t))

        if len(invalid_user_names) != 0:
            raise Exception("ERROR : unknown suite or test provided on command line : {0}".format(invalid_user_names))

        return u_suites, u_tests

    def _create_tests(self, user_suites, user_tests):
        all_tests = user_tests
        for s in user_suites:
            for t in self._available_suites[s].split():
                all_tests.append(t)

        for t in all_tests:
            test = RegressionTest()
            test.setup(self._executable_args, self._default_test_criteria, self._available_tests[t])
            self._tests.append(test)


def commandline_options():
    parser = argparse.ArgumentParser(description='Run a pflotran regression tests or suite of tests.')
    parser.add_argument('--backtrace', action='store_true',
                        help='show exception backtraces as extra debugging output')
    parser.add_argument('-c', '--config-file', nargs=1, required=True,
                        help='test configuration file to use')
    parser.add_argument('--debug', action='store_true',
                        help='extra debugging output')
    parser.add_argument('-d', '--dry-run',
                        default=False, action='store_true',
                        help='perform a dry run, setup the test commands but don\'t run them')
    parser.add_argument('-e', '--executable', nargs=1, default=None,
                        help='path to executable to use for testing')
    parser.add_argument('--list-suites', action='store_true',
                        help='print the list of test suites from the config file and exit')
    parser.add_argument('--list-tests', action='store_true',
                        help='print the list of tests from the config file and exit')
    parser.add_argument('-s', '--suites', nargs="+", default=[],
                        help='space separated list of test suite names')
    parser.add_argument('-t', '--tests', nargs="+", default=[],
                        help='space separated list of test names')
    parser.add_argument('-u', '--update',
                        action="store_true", default="store_false",
                        help='update the tests listed by the "--tests" option, with the current output becoming the new gold standard')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose output')

    options = parser.parse_args()
    return options


def main(options):
    test_manager = RegressionTestManager()

    if options.debug:
        test_manager._debug = True

    test_manager.generate_tests(options.config_file[0],
                                options.suites,
                                options.tests)

    if options.debug:
        print(70*'-')
        print(test_manager)

    if options.list_suites == True:
        test_manager.display_available_suites()

    if options.list_tests == True:
        test_manager.display_available_tests()

    if options.executable == None:
        options.dry_run = True
        # geh: kludge to permit call to test_manager.run_tests() below
        options.executable = []
        options.executable.append('executable')

    test_manager.run_tests(options.executable[0],
                           options.dry_run,
                           options.verbose,
                           options.update)

    if options.update:
        test_manager.update_test_results(options.tests)

    return test_manager.status()

if __name__ == "__main__":
    options = commandline_options()
    try:
        status = main(options)
        sys.exit(status)
    except Exception as e:
        print(str(e))
        if options.backtrace:
            traceback.print_exc()
        os._exit(1)
