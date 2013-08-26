#!/usr/bin/env python

"""python script to find all pflotran source files that #include
definitions.h and refactor them to "use PFLOTRAN_Constants_module"

"""

from __future__ import print_function

import fileinput
import os
import sys
import traceback

def get_source_list(cwd):
    """
    Get a list of all F90 source files in the directory
    """
    source = []
    for f in os.listdir(cwd):
        if os.path.isfile(f):
            if f.endswith(".F90"):
                source.append(f)
    return source

def get_include_definitions(source):
    """
    Make a list of the F90 source files that include definitions.h
    """
    include_definitions = []
    for f in source:
        with open(f, 'r') as src:
            for line in src.readlines():
                if line.find("definitions.h") != -1:
                    include_definitions.append(f)
                    break
    return include_definitions

def replace_include_with_use(source):
    """Open each file and place a use pflotran constants statement at the
    begining of a module, remove include definitions.h

    Special logic required for checkpoint because there are two module
    is the same file.

    """
    for s in source:
        b = "{0}.bck".format(s)
        with open(s, 'r') as src, open(b, 'w') as bck:
            if (s.find("checkpoint.F90") > -1 or
                s.find("surface_checkpoint.F90") > -1):
                finished = True
            else:
                finished = False
            for line in src:
                if (line.find("end module Checkpoint_Header_module") > -1 or
                    line.find("end module Surface_Checkpoint_Header_module") > -1):
                    finished = False
                if finished is False:
                    if line.find("implicit none") != -1:
                        print("  use PFLOTRAN_Constants_module\n", file=bck)
                        finished = True
                if line.find('#include "definitions.h"') is -1:
                    print(line, end='', file=bck)
                else:
                    print('#include "finclude/petscsys.h"', file=bck)
        os.rename(b, s)

    
def main():
    if len(sys.argv) > 0:
        cwd = sys.argv[1]
        if not os.path.isdir(cwd):
            raise Exception("ERROR: commandline arguments must be directory names.")
    else:
        cwd = os.getcwd()
    cwd = os.path.abspath(cwd)
    source = get_source_list(cwd)
    source = get_include_definitions(source)
    replace_include_with_use(source)

    return 0

if __name__ == "__main__":
    try:
        status = main()
        sys.exit(status)
    except Exception as error:
        print(str(error))
        traceback.print_exc()
        sys.exit(1)
