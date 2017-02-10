#!/bin/env python
"""
Program to manage and run PFLOTRAN regression tests.

"""

from __future__ import print_function
from __future__ import division

import sys

if sys.hexversion < 0x02070000:
    print(70*"*")
    print("ERROR: PFLOTRAN's regression test manager requires python >= 2.7.x. ")
    print("It appears that you are running python {0}.{1}.{2}".format(
        sys.version_info[0], sys.version_info[1], sys.version_info[2]))
    print(70*"*")
    sys.exit(1)

import argparse
import datetime
import hashlib
import os
import pprint
import re
import shutil
import subprocess
import textwrap
import time
import traceback

if sys.version_info[0] == 2:
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser

# optional libraries
try:
    import h5py
except Exception as e:
    h5py = None



class TestStatus(object):
    """
    Simple class to hold status info.
    """
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
    _wall_time_re = re.compile(r"Time \(seconds\)")

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
        self._DISPLACEMENT = "displacement"
        self._STRAIN = "strain"
        self._STRESS = "stress"
        self._SOLUTION = "solution"
        self._RESIDUAL = "residual"
        self._TOL_VALUE = 0
        self._TOL_TYPE = 1
        self._TOL_MIN_THRESHOLD = 2
        self._TOL_MAX_THRESHOLD = 3
        self._PFLOTRAN_SUCCESS = 86
        self._RESTART_PREFIX = "tmp-restart"
        # misc test parameters
        self._pprint = pprint.PrettyPrinter(indent=2)
        self._txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")
        self._debug = False
        self._executable = None
        self._input_arg = "-input_prefix"
        self._input_suffix = "in"
        self._np = None
        self._pflotran_args = None
        self._stochastic_realizations = None
        self._restart_timestep = None
        self._compare_hdf5 = False
        self._timeout = 60.0
        self._skip_check_gold = False
        self._check_performance = False
        self._num_failed = 0
        self._test_name = None
        # assign default tolerances for different classes of variables
        # absolute min and max thresholds for determining whether to
        # compare to baseline, i.e. if (min_threshold <= abs(value) <=
        # max_threshold) then compare values. By default we use the
        # python definitions for this platform
        self._tolerance = {}
        self._tolerance[self._TIME] = [5.0, self._PERCENT, \
                                       0.0, sys.float_info.max]
        self._tolerance[self._DISCRETE] = [0, self._ABSOLUTE, 0, sys.maxsize]
        common = [self._CONCENTRATION, self._GENERIC, self._RATE, self._VOLUME_FRACTION, \
                  self._PRESSURE, self._SATURATION, self._RESIDUAL, \
                  self._DISPLACEMENT, self._STRESS, self._STRAIN]
        for t in common:
            self._tolerance[t] = [1.0e-12, self._ABSOLUTE, \
                                  0.0, sys.float_info.max]


    def __str__(self):
        message = "  {0} :\n".format(self.name())
        message += "    timeout = {0}\n".format(self._timeout)
        message += "    np = {0}\n".format(self._np)
        message += "    pflotran args :\n"
        message += "        input : {0}\n".format(self._input_arg)
        message += "        optional : {0}\n".format(self._pflotran_args)
        message += "    test criteria :\n"
        for k in self._tolerance:
            message += "        {0} : {1} [{2}] : {3} <= abs(value) <= {4}\n".format(
                k,
                self._tolerance[k][self._TOL_VALUE],
                self._tolerance[k][self._TOL_TYPE],
                self._tolerance[k][self._TOL_MIN_THRESHOLD],
                self._tolerance[k][self._TOL_MAX_THRESHOLD])

        return message

    def setup(self, cfg_criteria, test_data,
              timeout, check_performance, testlog):
        """
        Setup the test object

        cfg_criteria - dict from cfg file, all tests in file
        test_data - dict from cfg file, test specific
        timeout - list(?) from command line option
        check_performance - bool from command line option
        """
        self._test_name = test_data["name"]

        self._set_test_data(cfg_criteria, test_data,
                            timeout, check_performance, testlog)

    def name(self):
        return self._test_name

    def run(self, mpiexec, executable, dry_run, status, testlog):
        """Run the test.

        If this is a restart test, then we will copy the input file,
        append the restart flag and run the test a second time.

        """
        self._cleanup_generated_files()

        self._run_test(mpiexec, executable, self.name(), dry_run, status,
                       testlog)

        if self._restart_timestep is not None:
            restart_file = "{0}-ts{1}.chk".format(self.name(),
                                                self._restart_timestep)
            if os.path.isfile(restart_file):
                restart_name = "{0}-{1}".format(self._RESTART_PREFIX,
                                                self.name())
                shutil.copy("{0}.in".format(self.name()),
                            "{0}.in".format(restart_name))
                with open("{0}.in".format(restart_name), 'a') as tempfile:
                    tempfile.write("RESTART {0}\n".format(restart_file))
                self._run_test(mpiexec, executable, restart_name, dry_run,
                               status, testlog)
            elif not status.skipped:
                status.fail = 1
                message = self._txtwrap.fill(
                    "ERROR: restart test '{0}' did not generate a checkpoint "
                    "file at the specified step: '{1}' ('{2}'). This can "
                    "occur if the checkpoint interval and restart step are "
                    "not consistent or the run failed.".format(
                        self.name(), self._restart_timestep, restart_file))
                print("".join(['\n', message, '\n']), file=testlog)

    def _run_test(self, mpiexec, executable, test_name, dry_run, status, testlog):
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
        #geh: we now set the successful exit code through the command line
        #     so that users are not confused by any error codes reported by
        #     wrapper libraries (e.g. MPICH2) due to the non-zero
        #     PFLOTRAN_SUCCESS.
        command.append("-successful_exit_code")
        command.append("%d" % self._PFLOTRAN_SUCCESS)
        if self._input_arg is not None:
            command.append(self._input_arg)
            command.append(test_name)

        if self._pflotran_args is not None:
            # assume that we already called split() on the
            # pflotran_args string. That will always return a list (it
            # may be empty)
            command = command + self._pflotran_args

        if not dry_run:
            print("    cd {0}".format(os.getcwd()), file=testlog)
            print("    {0}".format(" ".join(command)), file=testlog)
            run_stdout = open(test_name + ".stdout", 'w')
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
                        "{1} seconds.".format(test_name, self._timeout))
                    print(''.join(['\n', message, '\n']), file=testlog)
            finish = time.time()
            print("    # {0} : run time : {1:.2f} seconds".format(test_name, finish - start), file=testlog)
            pflotran_status = abs(proc.returncode)
            run_stdout.close()
        # pflotran returns 0 on an error (e.g. can't find an input
        # file), 86 on success. 59 for timeout errors?
        if pflotran_status != self._PFLOTRAN_SUCCESS:
            status.fail = 1
            message = self._txtwrap.fill(
                "FAIL : {name} : pflotran return an error "
                "code ({status}) indicating the simulation may have "
                "failed. Please check '{name}.out' and '{name}.stdout' "
                "for error messages (included below).".format(
                    name=test_name, status=pflotran_status))
            print("".join(['\n', message, '\n']), file=testlog)
            print("~~~~~ {0}.stdout ~~~~~".format(test_name), file=testlog)
            try:
                with open("{0}.stdout".format(test_name), 'r') as tempfile:
                    shutil.copyfileobj(tempfile, testlog)
            except Exception as e:
                print("   Error opening file: {0}.stdout\n    {1}".format(test_name, e))
            print("~~~~~ {0}.out ~~~~~".format(test_name), file=testlog)
            try:
                with open("{0}.out".format(test_name), 'r') as tempfile:
                    shutil.copyfileobj(tempfile, testlog)
            except Exception as e:
                print("   Error opening file: {0}.out\n    {1}".format(test_name, e))
            print("~~~~~~~~~~", file=testlog)

    def _cleanup_generated_files(self):
        """Cleanup old generated files that may be hanging around from a
        previous run by renaming them with .old appended to the name.

        NOTE:

          - Do NOT match files with something like "if self.name() in
            filename" or "if filenamename.startswith(self.name())"
            because this will capture files from another test with a
            similar name, e.g.  "test_flow" and "test_flow_np4" would
            both be captured.

          - This assumes that all files with the listed suffixes are
            old. That means checkpoint files must be generated, not
            saved.

        """
        # files from a normal run
        suffixes = ["regression", "out", "stdout"]
        for suffix in suffixes:
            name = "{0}.{1}".format(self.name(), suffix)
            if os.path.isfile(name):
                os.rename(name, name + ".old")

        # files from a stochastic run
        if self._stochastic_realizations is not None:
            for i in range(1, self._stochastic_realizations + 1):
                run_id = "R{0}".format(i)
                for suffix in suffixes:
                    name = "{0}{1}.{2}".format(self.name(), run_id, suffix)
                    if os.path.isfile(name):
                        os.rename(name, name + ".old")

        # temp files from a restart run
        if self._restart_timestep is not None:
            suffixes.append("in")
            for suffix in suffixes:
                name = "{0}-{1}.{2}".format(self._RESTART_PREFIX, self.name(), suffix)
                if os.path.isfile(name):
                    os.rename(name, name + ".old")

            # checkpoint/restart files, both from the original and restart run
            cwd = os.getcwd()
            for entry in os.listdir(cwd):
                if os.path.isfile(entry):
                    search_checkpoint = "^({0}-)?{1}-(ts[\d]+|restart).chk$".format(
                        self._RESTART_PREFIX, self.name())
                    if re.search(search_checkpoint, entry):
                        os.rename(entry, entry + ".old")

    def check(self, status, testlog):
        """
        Check the test results against the gold standard

        Some tests generate multiple regression files, i.e. stochastic
        realizations. In that case, we need to loop over a group of
        regression files.
        """
        run_id = ''
        if self._stochastic_realizations is not None:
            print("\n    {0} : Test has {1} stochastic realizations\n".format(
                self.name(), self._stochastic_realizations), file=testlog)
            for i in range(1, self._stochastic_realizations + 1):
                run_id = "R{0}".format(i)
                self._check_gold(status, run_id, testlog)
        else:
            self._check_gold(status, run_id, testlog)

        if self._restart_timestep is not None:
            self._check_restart(status, testlog)

        if self._compare_hdf5:
            if h5py is not None:
                self._check_hdf5(status, testlog)
            else:
                print("    h5py not in python path. Skipping hdf5 check.", file=testlog)

    def _check_gold(self, status, run_id, testlog):
        """
        Test the output from the run against the known "gold standard"
        output and determine if the test succeeded or failed.

        We return zero on success, one on failure so that the test
        manager can track how many tests succeeded and failed.
        """
        if self._skip_check_gold:
            message = "    Skipping comparison to regression gold file (only test if model runs to completion)."
            print("".join(['\n', message, '\n']), file=testlog)
            return

        gold_filename = self.name() + run_id + ".regression.gold"
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

        current_filename = self.name() + run_id + ".regression"
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
        for section in gold_sections:
            if section not in current_sections:
                self._num_failed += 1
                print("    FAIL: section '{0}' is in the gold output, but "
                      "not the current output.".format(section), file=testlog)

        # look for sections that are in current but not gold
        for section in current_sections:
            if section not in gold_sections:
                self._num_failed += 1
                print("    FAIL: section '{0}' is in the current output, "
                      "but not the gold output.".format(section), file=testlog)

        # compare common sections
        for section in gold_sections:
            if section in current_sections:
                try:
                    self._num_failed += self._compare_sections(
                        gold_sections[section], current_sections[section],
                        testlog)
                except Exception as error:
                    self._num_failed += 1
                    print(error, file=testlog)

        if self._num_failed > 0:
            status.fail = 1

    def _check_restart(self, status, testlog):
        """Check that binary restart files are bit for bit after a restart.

        The entire restart file should be bit for bit, both the "big"
        data and the metadata. This means we can take a hash of the
        file and report even a single bit difference between the files
        as a failure. If the meta data is allowd to be different, we
        will need a more sophisticated way of diffing the files.

        """

        orig_filename="{0}-restart.chk".format(self.name())
        restart_filename="{0}-{1}".format(self._RESTART_PREFIX, orig_filename)

        print("\n    comparing restart files:\n        {0}\n        {1}".format(
            orig_filename, restart_filename), file=testlog)

        orig_hash = self._get_binary_restart_hash(orig_filename,
                                                  status, testlog)
        restart_hash = self._get_binary_restart_hash(restart_filename,
                                                     status, testlog)

        if orig_hash is not False and restart_hash is not False:
            if orig_hash != restart_hash:
                print("    FAIL: final restart files are not bit for bit "
                      "identical.", file=testlog)
                status.fail = 1
            else:
                print("    bit for bit restart test passed.\n", file=testlog)

    def _get_binary_restart_hash(self, filename, status, testlog):
        """Get the sha1 hash of a restart file. The hash should be different
        if even one bit is different in the file.

        """
        hash_digest = False
        if os.path.isfile(filename):
            restart_hash = hashlib.sha1()
            with open(filename, mode='rb') as restart_file:
                restart_hash.update(restart_file.read())
            hash_digest = restart_hash.hexdigest()
        else:
            message = self._txtwrap.fill(
                "FAIL: could not find restart file '{0}'".format(filename))
            print("".join(['\n', message, '\n']), file=testlog)
            status.fail = 1
        return hash_digest

    def _check_hdf5(self, status, testlog):
        """Check that output hdf5 file has not changed from the baseline.

        Note: we open the files here and do the comparison in another
        function so it can be unit tested.

        """

        filename="{0}.h5".format(self.name())
        try:
            h5_current = h5py.File(filename, 'r')
        except Exception as e:
            print("    FAIL: Could not open file: '{0}'".format(filename), file=testlog)
            status.fail = 1
            h5_current = None

        filename = "{0}.gold".format(filename)
        try:
            h5_gold = h5py.File(filename, 'r')
        except Exception as e:
            print("    FAIL: Could not open file: '{0}'".format(filename), file=testlog)
            status.fail = 1
            h5_gold = None

        if h5_gold is not None and h5_current is not None:
            self._compare_hdf5_data(h5_current, h5_gold, status, testlog)

    def _compare_hdf5_data(self, h5_current, h5_gold, status, testlog):
        """Check that output hdf5 file has not changed from the baseline.

        The focus is on the meta-data, e.g. groups, datasets, vector
        sizes, etc. We assume that the *data* in the hdf5 file
        is already being checked the regression files.

        Note: We can't just h5dump the file and diff or campare
        hashes because the provenance information may have changed!

        """

        if len(h5_current.keys()) != len(h5_gold.keys()):
            status.fail = 1
            print("    FAIL: current and gold hdf5 files do not have the same number of groups!", file=testlog)
            print("    h5_gold : {0}".format(h5_gold.keys()), file=testlog)
            print("    h5_current : {0}".format(h5_current.keys()), file=testlog)
            return

        for group in h5_gold:
            if group == "Provenance":
                continue
            if len(h5_current[group].keys()) != len(h5_gold[group].keys()):
                status.fail = 1
                print("    FAIL: group '{0}' does not have the same number of datasets!".format(group), file=testlog)
                print("    h5_gold : {0} : {1}".format(
                    group, h5_gold[group].keys()), file=testlog)
                print("    h5_current : {0} : {1}".format(
                    group, h5_current[group].keys()), file=testlog)
            else:
                for dataset in h5_gold[group].keys():
                    if not h5_current[group].get(dataset):
                        status.fail = 1
                        print("    FAIL: current group '{0}' does not have dataset '{1}'!".format(group, dataset), file=testlog)
                    else:
                        if h5_gold[group][dataset].shape != h5_current[group][dataset].shape:
                            status.fail = 1
                            print("    FAIL: current dataset '/{0}/{1}' does not have the correct shape!".format(group, dataset), file=testlog)
                            print("        gold : {0}".format(
                                h5_gold[group][dataset].shape), file=testlog)
                            print("        current : {0}".format(
                                h5_current[group][dataset].shape), file=testlog)
                        if h5_gold[group][dataset].dtype != h5_current[group][dataset].dtype:
                            status.fail = 1
                            print("    FAIL: current dataset '/{0}/{1}' does not have the correct data type!".format(group, dataset), file=testlog)
                            print("        gold : {0}".format(
                                h5_gold[group][dataset].dtype), file=testlog)
                            print("        current : {0}".format(
                                h5_current[group][dataset].dtype), file=testlog)


        if status.fail == 0:
            print("    Passed hdf5 check.", file=testlog)


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
        except Exception as error:
            status = 1
            message = str(error)
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
#            raise RuntimeError("ERROR: test '{0}' was classified as new, "
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
        except Exception as error:
            status = 1
            message = str(error)
            message += "\nFAIL : Could not rename '{0}' to '{1}'. "
            message += "Please rename the file manually!".format(current_name,
                                                                 gold_name)
            message += "    mv {0} {1}".format(current_name, gold_name)
            print(message, file=testlog)
            status.fail = 1

    def _get_sections(self, output):
        """
        Get the sections from the regression file.

        Each section in the regression test file looks like:
        -- TYPE: NAME --
           key: value
           key: value
        """
        name_re = re.compile(r"^--[\s]+([\w]+):(.*)[\s]+--$")
        sections = {}
        sect = {}
        for line in output:
            match = name_re.match(line)
            if match:
                # save the old section, if any
                if 'name' in sect:
                    sections[sect['name']] = sect
                name = match.group(2).strip()
                data_type = match.group(1)
                sect = {}
                sect['name'] = name
                sect['type'] = data_type
            else:
                temp = line.split(':')
                name = temp[0].strip()
                value = temp[1].strip()
                sect[name] = value
        # add the final section
        if 'name' in sect:
            sections[sect['name']] = sect

        return sections

    def _compare_sections(self, gold_section, current_section, testlog):
        """
        Compare the fields of the current section.
        """
        if gold_section["name"] != current_section["name"]:
            raise RuntimeError("ERROR: an internal error occured. "
                               "compare_sections receive different sections. "
                               "gold = '{0}' current = '{1}'".format(
                                   gold_section["name"], current_section["name"]))
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
                            except Exception as error:
                                section_status += 1
                                print("ERROR: {0} : {1}.\n  {2}".format(
                                    self.name(), k, str(error)), file=testlog)

        if False:
            print("    {0} : status : {1}".format(
                name, section_status), file=testlog)
        return section_status

    def _compare_values(self, name, key, previous, current, testlog):
        """
        Compare the values using the appropriate tolerance and criteria.

        NOTE(bja): previous and current come into this function as
        strings. We don't know if they should be floats or ints (or
        possibly strings) until we know what the data type is. For
        'discrete' or 'solution' variables, we have to do further work
        to figure it out!
        """
        status = 0
        tol = None
        key = key.lower()
        if (key == self._CONCENTRATION or
            key == self._GENERIC or
            key == self._RATE or
            key == self._VOLUME_FRACTION or
            key == self._PRESSURE or
            key == self._SATURATION or
            key == self._DISPLACEMENT or
            key == self._STRAIN or
            key == self._STRESS):
            previous = float(previous)
            current = float(current)
            tol = self._tolerance[key]
        elif key.lower() == self._SOLUTION:
            previous, current, tol = \
                self._compare_solution(name, previous, current)
        elif key.lower() == self._DISCRETE:
            previous, current, tol = \
                self._compare_discrete(name, previous, current)
        else:
            raise RuntimeError(
                "WARNING: the data caterogy '{0}' for '{1}' is not a known "
                "data category.".format(key, name))

        tolerance_type = tol[self._TOL_TYPE]
        tolerance = tol[self._TOL_VALUE]
        min_threshold = tol[self._TOL_MIN_THRESHOLD]
        max_threshold = tol[self._TOL_MAX_THRESHOLD]

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
            # should never get here....
            raise RuntimeError("ERROR: unknown test tolerance_type '{0}' for "
                               "variable '{1}, {2}.'".format(tolerance_type,
                                                          name, key))

        # base comparison to threshold on previous (gold) because it
        # is the known correct value!
        if min_threshold <= abs(previous) and abs(previous) <= max_threshold:
            if delta > tolerance:
                status = 1
                print("    FAIL: {0} : {1} > {2} [{3}]".format(
                    name, delta, tolerance, tolerance_type), file=testlog)
            elif self._debug:
                print("    PASS: {0} : {1} <= {2} [{3}]".format(
                    name, delta, tolerance, tolerance_type), file=testlog)
        else:
            print("    SKIP: {0} : gold value ({1}) outside threshold range: "
                  "{2} <= abs(value) <= {3}".format(
                      name, previous, min_threshold, max_threshold),
                  file=testlog)

        return status

    def _compare_solution(self, name, previous, current):
        """
        Tolerances for the solution variables depend on what field we are
        testing.

        """
        section = name.split(':')[0]
        param = name.split(':')[1]
        param = param.strip()
        if param == "Time (seconds)":
            previous = float(previous)
            current = float(current)
            tolerance = self._tolerance[self._TIME]
        elif param == "Time Steps":
            previous = int(previous)
            current = int(current)
            tolerance = self._tolerance[self._DISCRETE]
        elif param == "Newton Iterations":
            previous = int(previous)
            current = int(current)
            tolerance = self._tolerance[self._DISCRETE]
        elif param == "Solver Iterations":
            previous = int(previous)
            current = int(current)
            tolerance = self._tolerance[self._DISCRETE]
        elif param == "Time Step Cuts":
            previous = int(previous)
            current = int(current)
            tolerance = self._tolerance[self._DISCRETE]
        elif param == "Solution 2-Norm":
            previous = float(previous)
            current = float(current)
            tolerance = self._tolerance[self._GENERIC]
        elif param == "Residual 2-Norm":
            previous = float(previous)
            current = float(current)
            tolerance = self._tolerance[self._RESIDUAL]
        else:
            raise RuntimeError("ERROR: unknown variable '{0}' in solution "
                               "section '{1}'".format(param, section))

        return previous, current, tolerance

    def _compare_discrete(self, name, previous, current):
        """
        Get tolerances for a "discrete" section, e.g. Material IDs.

        Discrete values are integers, except when we are
        looking at the mean of a discrete variable. Then we
        may (probably) have a floating point value!
        """
        mean_re = re.compile("Mean")
        have_mean = mean_re.search(name)
        if not have_mean:
            # previous and current must both be ints...
            try:
                previous = int(previous)
            except ValueError:
                raise RuntimeError(
                    "ERROR: discrete gold value must be an int: '{0}' = {1}.".format(
                        name, previous))
            try:
                current = int(current)
            except ValueError:
                raise RuntimeError(
                    "ERROR: discrete current value must be an int: '{0}' = {1}.".format(
                        name, current))

            tolerance = self._tolerance[self._DISCRETE]
        else:
            previous = float(previous)
            current = float(current)
            tolerance = self._tolerance[self._GENERIC]

        return previous, current, tolerance

    def _set_test_data(self, cfg_criteria, test_data, timeout,
                       check_performance, testlog):
        """
        Set the test criteria for different categories of variables.
        """
        self._np = test_data.pop('np', None)

        self._pflotran_args = test_data.pop('input_arguments', None)
        if self._pflotran_args is not None:
            # save the arg list so we can append it to the run command
            self._pflotran_args = self._pflotran_args.split()
            # additional processing that may change the test manager behavior
            if "-stochastic" in self._pflotran_args:
                if "-num_realizations" in self._pflotran_args:
                    index = self._pflotran_args.index("-num_realizations")
                    if len(self._pflotran_args) > index + 1:
                        self._stochastic_realizations = int(self._pflotran_args[index + 1])
                    else:
                        raise RuntimeError(
                            "ERROR: num_realizations requires an integer parameter N")
                else:
                    raise RuntimeError(
                        "ERROR : stochastic simulations require a "
                        "num_realizations flag as well. "
                        "test : {0}".format(self.name()))

        self._restart_timestep = test_data.pop('restart_timestep', None)
        if self._restart_timestep is not None:
            try:
                self._restart_timestep = int(self._restart_timestep)
            except ValueError as error:
                raise ValueError("ERROR: restart_timestep must be an integer value. "
                                "test : {0}".format(self.name()))

        self._compare_hdf5 = test_data.pop('compare_hdf5', None)
        if self._compare_hdf5 is not None:
            if self._compare_hdf5[0] in ["T", "t", "Y", "y"]:
                self._compare_hdf5 = True
            else:
                self._compare_hdf5 = False

        self._skip_check_gold = test_data.pop('skip_check_gold', None)
        if self._skip_check_gold is not None:
            if self._skip_check_gold[0] in ["T", "t", "Y", "y"]:
                self._skip_check_gold = True

        self._check_performance = check_performance

        # timeout : preference (1) command line (2) test data (3) class default
        self._timeout = float(test_data.pop('timeout', self._timeout))
        if timeout:
            self._timeout = float(timeout[0])

        self._set_criteria(self._TIME, cfg_criteria, test_data)

        self._set_criteria(self._CONCENTRATION, cfg_criteria, test_data)

        self._set_criteria(self._GENERIC, cfg_criteria, test_data)

        self._set_criteria(self._DISCRETE, cfg_criteria, test_data)

        self._set_criteria(self._RATE, cfg_criteria, test_data)

        self._set_criteria(self._VOLUME_FRACTION, cfg_criteria, test_data)

        self._set_criteria(self._PRESSURE, cfg_criteria, test_data)

        self._set_criteria(self._SATURATION, cfg_criteria, test_data)

        self._set_criteria(self._DISPLACEMENT, cfg_criteria, test_data)

        self._set_criteria(self._STRAIN, cfg_criteria, test_data)

        self._set_criteria(self._STRESS, cfg_criteria, test_data)

    def _set_criteria(self, key, cfg_criteria, test_data):
        """
        Our prefered order for selecting test criteria is:
        (1) test data section of the config file
        (2) config-file wide defaults
        (3) hard coded class default
        """
        if key in test_data:
            criteria = self._validate_criteria(key, test_data[key])
        elif key in cfg_criteria:
            criteria = self._validate_criteria(key, cfg_criteria[key])
        elif key in self._tolerance:
            # already a correctly formatted list and stored
            criteria = [None]
        else:
            raise RuntimeError("ERROR : tolerance for data type '{0}' could "
                               "not be determined from the config file or "
                               "default values!".format(key))
        for i, c in enumerate(criteria):
            if c is not None:
                self._tolerance[key][i] = c

    def _validate_criteria(self, key, test_data):
        """
        Validate the criteria string from a config file.

        Valid input configurations are:

        * key = tolerance type

        * key = tolerance type [, min_threshold value] [, max_threshold value]

        where min_threshold and max_threshold are optional

        """
        criteria = 4*[None]
        test_data = test_data.split(",")
        test_criteria = test_data[0]
        value = test_criteria.split()[0]
        try:
            value = float(value)
        except Exception:
            raise RuntimeError("ERROR : Could not convert '{0}' test criteria "
                               "value '{1}' into a float!".format(key, value))
        criteria[self._TOL_VALUE] = value

        criteria_type = test_criteria.split()[1]
        if (criteria_type.lower() != self._PERCENT and
            criteria_type.lower() != self._ABSOLUTE and
                criteria_type.lower() != self._RELATIVE):
            raise RuntimeError("ERROR : invalid test criteria string '{0}' "
                               "for '{1}'".format(criteria_type, key))
        criteria[self._TOL_TYPE] = criteria_type

        thresholds = {}
        for t in range(1, len(test_data)):
            name = test_data[t].split()[0].strip()
            value = test_data[t].split()[1].strip()
            try:
                value = float(value)
            except Exception:
                raise RuntimeError(
                    "ERROR : Could not convert '{0}' test threshold '{1}'"
                    "value '{2}' into a float!".format(key, name, value))
            thresholds[name] = value
        value = thresholds.pop("min_threshold", None)
        criteria[self._TOL_MIN_THRESHOLD] = value
        value = thresholds.pop("max_threshold", None)
        criteria[self._TOL_MAX_THRESHOLD] = value
        if len(thresholds) > 0:
            raise RuntimeError("ERROR: test {0} : unknown criteria threshold: {1}",
                               key, thresholds)

        return criteria




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
        self._default_test_criteria = {}
        self._available_tests = {}
        self._available_suites = {}
        self._tests = []
        self._txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")
        self._pprint = pprint.PrettyPrinter(indent=2)

    def __str__(self):
        data = "Regression Test Manager :\n"
        data += "    configuration file : {0}\n".format(self._config_filename)
        data += "    default test criteria :\n"
        data += self._dict_to_string(self._default_test_criteria)
        data += "    suites :\n"
        data += self._dict_to_string(self._available_suites)
        data += "    available tests :\n"
        data += self._dict_to_string(self._available_tests)

        data += "Tests :\n"
        for test in self._tests:
            data += test.__str__()

        return data

    def debug(self, debug):
        self._debug = debug

    def num_tests(self):
        return len(self._tests)

    def generate_tests(self, config_file, user_suites, user_tests,
                       timeout, check_performance, testlog):
        """
        Read the config file, validate the input and generate the test objects.
        """
        self._read_config_file(config_file)
        self._validate_suites()
        user_suites, user_tests = self._validate_user_lists(user_suites,
                                                            user_tests, testlog)
        self._create_tests(user_suites, user_tests, timeout, check_performance,
                           testlog)

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
        """
        Run the test and check the results.
        """
        if dry_run:
            print("Dry run:")
        print("Running tests from '{0}':".format(self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for test in self._tests:
            status = TestStatus()
            self._test_header(test.name(), testlog)

            test.run(mpiexec, executable, dry_run, status, testlog)

            if not dry_run and status.skipped == 0:
                test.check(status, testlog)

            self._add_to_file_status(status)

            self._test_summary(test.name(), status, dry_run,
                               "passed", "failed", testlog)

        self._print_file_summary(dry_run, "passed", "failed", testlog)

    def _check_only(self, dry_run, testlog):
        """
        Recheck the regression files from a previous run.
        """
        if dry_run:
            print("Dry run:")
        print("Checking existing test results from '{0}':".format(
            self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for test in self._tests:
            status = TestStatus()
            self._test_header(test.name(), testlog)

            if not dry_run and status.skipped == 0:
                test.check(status, testlog)

            self._add_to_file_status(status)

            self._test_summary(test.name(), status, dry_run,
                               "passed", "failed", testlog)

        self._print_file_summary(dry_run, "passed", "failed", testlog)

    def _run_new(self, mpiexec, executable, dry_run, testlog):
        """
        Run the tests and create new gold files.
        """
        if dry_run:
            print("Dry run:")

        print("New tests from '{0}':".format(self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for test in self._tests:
            status = TestStatus()
            self._test_header(test.name(), testlog)

            test.run(mpiexec, executable, dry_run, status, testlog)

            if not dry_run and status.skipped == 0:
                test.new_test(status, testlog)
            self._add_to_file_status(status)
            self._test_summary(test.name(), status, dry_run,
                               "created", "error creating new test files.", testlog)

        self._print_file_summary(dry_run, "created", "could not be created", testlog)

    def _run_update(self, mpiexec, executable, dry_run, testlog):
        """
        Run the tests and update the gold file with the current output
        """
        if dry_run:
            print("Dry run:")
        print("Updating tests from '{0}':".format(self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for test in self._tests:
            status = TestStatus()
            self._test_header(test.name(), testlog)
            test.run(mpiexec, executable, dry_run, status, testlog)

            if not dry_run and status.skipped == 0:
                test.update(status, testlog)
            self._add_to_file_status(status)
            self._test_summary(test.name(), status, dry_run,
                               "updated", "error updating test.", testlog)

        self._print_file_summary(dry_run, "updated", "could not be updated", testlog)

    def _test_header(self, name, testlog):
        """
        Write a header to the log file to separate tests.
        """
        print(40 * '-', file=testlog)
        print("{0}... ".format(name), file=testlog)

    def _test_summary(self, name, status, dry_run,
                      success_message, fail_message, testlog):
        """
        Write the test status information to stdout and the test log.
        """
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
        """
        Print a summary of the results for this config file
        """
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
        """
        Add the current test status to the overall status for the file.
        """
        self._file_status.fail += status.fail
        self._file_status.warning += status.warning
        self._file_status.error += status.error
        self._file_status.skipped += status.skipped
        self._file_status.test_count += 1

    def run_status(self):
        return self._file_status

    def display_available_tests(self):
        print("Available tests: ")
        for test in sorted(self._available_tests.keys()):
            print("    {0}".format(test))

    def display_available_suites(self):
        print("Available test suites: ")
        for suite in self._available_suites:
            print("    {0} :".format(suite))
            for test in self._available_suites[suite].split():
                print("        {0}".format(test))

    def _read_config_file(self, config_file):
        """
        Read the configuration file.

        Sections : The config file will have known sections:
        "suites", "default-test-criteria".

        All other sections are assumed to be test names.
        """
        if config_file is None:
            raise RuntimeError("Error, must provide a config filename")
        self._config_filename = config_file
        config = config_parser()
        config.read(self._config_filename)

        if config.has_section("default-test-criteria"):
            self._default_test_criteria = \
                self._list_to_dict(config.items("default-test-criteria"))

        if config.has_section("suites"):
            self._available_suites = \
                self._list_to_dict(config.items("suites"))

        self._identify_tests(config)

    def _identify_tests(self, config):
        """
        Create a list of all tests in a config file.

        Assumes every section is a test except for some fixed section
        names

        """
        # section names are test names
        test_names = config.sections()

        # remove the fixed section names
        if config.has_section("default-test-criteria"):
            test_names.remove("default-test-criteria")
        if config.has_section("suites"):
            test_names.remove("suites")

        # all remaining sections should be individual tests
        for test in test_names:
            self._available_tests[test] = self._list_to_dict(config.items(test))
            self._available_tests[test]['name'] = test

    def _dict_to_string(self, data):
        """
        Format dictionary key-value pairs in a string
        """
        temp = ""
        for key, value in data.items():
            temp += "        {0} : {1}\n".format(key, value)
        return temp

    def _list_to_dict(self, input_list):
        """
        Convert a list of key-value pairs into a dictionary.
        """
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
        for suite in self._available_suites:
            suite_tests = self._available_suites[suite].split()
            if len(suite_tests) == 0:
                empty_suites.append(suite)
            else:
                # validate the list
                for test in suite_tests:
                    if test not in self._available_tests:
                        name = "suite : '{0}' --> test : '{1}'".format(
                            suite, test)
                        invalid_tests.append(name)

        for suite in empty_suites:
            # empty suite, warn the user and remove it from the list
            del self._available_suites[suite]
            print("DEV WARNING : {0} : cfg validation : empty suite "
                  ": '{1}'".format(self._config_filename, suite))

        if len(invalid_tests) != 0:
            raise RuntimeError("ERROR : suites contain unknown tests in "
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
            for suite in user_suites:
                if suite.lower() in self._available_suites:
                    u_suites.append(suite.lower())
                else:
                    message = self._txtwrap.fill(
                        "WARNING : {0} : Skipping requested suite '{1}' (not "
                        "present, misspelled or empty).".format(
                            self._config_filename, suite))
                    print(message, file=testlog)

            u_tests = []
            for test in user_tests:
                if test in self._available_tests:
                    u_tests.append(test.lower())
                else:
                    message = self._txtwrap.fill(
                        "WARNING : {0} : Skipping test '{1}' (not present or "
                        "misspelled).".format(self._config_filename, test))
                    print(message, file=testlog)

        return u_suites, u_tests

    def _create_tests(self, user_suites, user_tests, timeout, check_performance,
                      testlog):
        """
        Create the test objects for all user specified suites and tests.
        """
        all_tests = user_tests
        for suite in user_suites:
            for test in self._available_suites[suite].split():
                all_tests.append(test)

        for test in all_tests:
            try:
                new_test = RegressionTest()
                new_test.setup(self._default_test_criteria,
                               self._available_tests[test], timeout,
                               check_performance, testlog)
                self._tests.append(new_test)
            except Exception as error:
                raise RuntimeError("ERROR : could not create test '{0}' from "
                                   "config file '{1}'. {2}".format(
                                       test, self._config_filename, str(error)))


def commandline_options():
    """
    Process the command line arguments and return them as a dict.
    """
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
                raise RuntimeError("ERROR: can not search for config files "
                                   "in '{0}' because it is not a "
                                   "directory.".format(base_dir))

    # add the explicitly listed config files
    if options.config_files is not None:
        for config_file in options.config_files:
            if not os.path.isabs(config_file):
                config_file = os.path.abspath(config_file)
            if os.path.isfile(config_file):
                config_file_list.append(config_file)
            else:
                raise RuntimeError("ERROR: specified config file '{0}' is not a "
                                   "file!".format(config_file))

    if options.debug:
        print("\nKnown config files:")
        for config_file in config_file_list:
            print("    {0}".format(config_file))

    if len(config_file_list) == 0:
        raise RuntimeError("ERROR: no config files were found. Please specify a "
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
            raise RuntimeError("ERROR: can not update gold regression files "
                               "during a recursive search for config files.")

    if options.update and options.new_tests:
        raise RuntimeError("ERROR: can not create new tests and update gold "
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
            raise RuntimeError("ERROR: executable is not a valid file: "
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
###           raise RuntimeError("ERROR: mpiexec is not a valid file: "
###                           "'{0}'".format(mpiexec))
    else:
        message = ("\n** WARNING ** : mpiexec was not provided on the command line.\n"
                   "                All parallel tests will be skipped!\n")
        print(message, file=sys.stdout)
        print(message, file=testlog)

    return mpiexec


def summary_report_by_file(report, outfile):
    """
    Summarize the results for each config file.
    """
    print(70 * '-', file=outfile)
    print("Regression test file summary:", file=outfile)
    for filename in report:
        line = "    {0}... {1} tests : ".format(filename, report[filename].test_count)
        if report[filename].warning > 0:
            line = "{0} {1} test warnings, ".format(line, report[filename].warning)
        if report[filename].error > 0:
            line = "{0} {1} test errors, ".format(line, report[filename].error)

        if report[filename].test_count == 0:
            line = "{0}... no tests were run.".format(line)
        else:
            if report[filename].fail > 0:
                line = "{0} {1} tests failed, ".format(line, report[filename].fail)
            if report[filename].skipped > 0:
                line = "{0} {1} tests skipped, ".format(line, report[filename].skipped)
            if report[filename].fail == 0 and report[filename].skipped == 0:
                line = "{0} all tests passed".format(line)
            else:
                num_passed = (report[filename].test_count - report[filename].fail -
                              report[filename].skipped)
                line = "{0} {1} passed.".format(line, num_passed)

        print("{0}".format(line), file=outfile)

    print("\n", file=outfile)


def summary_report(run_time, report, outfile):
    """
    Overall summary of test results
    """
    print(70 * '-', file=outfile)
    print("Regression test summary:", file=outfile)
    print("    Total run time: {0:4g} [s]".format(run_time), file=outfile)
    test_count = 0
    num_failures = 0
    num_errors = 0
    num_warnings = 0
    num_skipped = 0
    for filename in report:
        test_count += report[filename].test_count
        num_failures += report[filename].fail
        num_errors += report[filename].error
        num_warnings += report[filename].warning
        num_skipped += report[filename].skipped

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
    """
    Append the results of a shell command to the test log
    """
    print("$ {0}".format(" ".join(command)), file=testlog)
    testlog.flush()
    with open(tempfile, "w") as tempinfo:
        subprocess.call(command, shell=False,
                        stdout=tempinfo,
                        stderr=subprocess.STDOUT)
        # NOTE(bja) 2013-06 : need a short sleep to ensure the
        # contents get written...?
        time.sleep(0.01)
    with open(tempfile, 'r') as tempinfo:
        shutil.copyfileobj(tempinfo, testlog)


def setup_testlog(txtwrap):
    """
    Create the test log and try to add some useful information about
    the environment, petsc and pflotran.
    """
    now = datetime.datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
    filename = "pflotran-tests-{0}.testlog".format(now)
    testlog = open(filename, 'w')
    print("  Test log file : {0}".format(filename))

    print("PFLOTRAN Regression Test Log", file=testlog)
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

    root_dir = os.getcwd()

    check_options(options)
    executable = check_for_executable(options)
    mpiexec = check_for_mpiexec(options, testlog)
    config_file_list = generate_config_file_list(options)

    print("Running pflotran regression tests :")

    # loop through config files, cd into the appropriate directory,
    # read the appropriate config file and run the various tests.
    start = time.time()
    report = {}
    for config_file in config_file_list:
        try:
            # NOTE(bja): the try block is inside this loop so that if
            # a single test throws an exception in a large batch of
            # tests, we can recover and at least try running the other
            # config files.
            print(80 * '=', file=testlog)

            # get the absolute path of the directory
            test_dir = os.path.dirname(config_file)
            # cd into the test directory so that the relative paths in
            # test files are correct
            os.chdir(test_dir)
            if options.debug:
                print("Changed to working directory: {0}".format(test_dir))

            test_manager = RegressionTestManager()

            if options.debug:
                test_manager.debug(True)

            # get the relative file name
            filename = os.path.basename(config_file)

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

            report_entry = config_file.split('regression_tests/')[1]
            report[report_entry] = test_manager.run_status()
            os.chdir(root_dir)
        except Exception as error:
            message = txtwrap.fill(
                "ERROR: a problem occured in file '{0}'.  This is "
                "probably an error with commandline options, the "
                "configuration file, or an internal error.  Please send "
                "this log file to pflotran-dev mailing list. The "
                "error is:\n{1}".format(config_file, str(error)))
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
    cmdl_options = commandline_options()
    try:
        suite_status = main(cmdl_options)
        sys.exit(suite_status)
    except Exception as error:
        print(str(error))
        if cmdl_options.backtrace:
            traceback.print_exc()
        sys.exit(1)
