#-----------------------------------
# Driver for all zone files
# to vset files
# Author: Satish Karra
# Date: July 11, 2013
#-----------------------------------

import sys

import os
print ('Top')
cmd = '/scratch/er/dharp/source/epd-7.3-1-rh3-x86_64/bin/python zone2vset.py top.zone top.vset'
failure = os.system(cmd)
if failure:
    print 'Failed to run zone2pflotran.py for top.zone'; sys.exit(1)

print ('Bottom')
cmd = '/scratch/er/dharp/source/epd-7.3-1-rh3-x86_64/bin/python zone2vset.py bottom.zone bottom.vset'
failure = os.system(cmd)
if failure:
    print 'Failed to run zone2pflotran.py for top.zone'; sys.exit(1)

print ('East')
cmd = '/scratch/er/dharp/source/epd-7.3-1-rh3-x86_64/bin/python zone2vset.py east.zone east.vset'
failure = os.system(cmd)
if failure:
    print 'Failed to run zone2pflotran.py for top.zone'; sys.exit(1)

print ('West')
cmd = '/scratch/er/dharp/source/epd-7.3-1-rh3-x86_64/bin/python zone2vset.py west.zone west.vset'
failure = os.system(cmd)
if failure:
    print 'Failed to run zone2pflotran.py for top.zone'; sys.exit(1)

print ('North')
cmd = '/scratch/er/dharp/source/epd-7.3-1-rh3-x86_64/bin/python zone2vset.py north.zone north.vset'
failure = os.system(cmd)
if failure:
    print 'Failed to run zone2pflotran.py for top.zone'; sys.exit(1)

print ('South')
cmd = '/scratch/er/dharp/source/epd-7.3-1-rh3-x86_64/bin/python zone2vset.py south.zone south.vset'
failure = os.system(cmd)
if failure:
    print 'Failed to run zone2pflotran.py for top.zone'; sys.exit(1)


