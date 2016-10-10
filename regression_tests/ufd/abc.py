import sys
import os
try:
  pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
  print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
  sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import pflotran as pft

path = []
path.append('.')

#files = pft.get_default_tec_filenames([5])
files = []
files.append('ufd_abc-obs-0.tec')
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(10,6))
plt.subplot(1,1,1)
#f.suptitle("ABC at 25 Years",fontsize=16)
plt.xlabel('Time [y]')
plt.ylabel('Concentration [M]')

#plt.grid(True)
plt.yscale('log')

minval = 1.e20
maxval = -1.e20
for ifile in range(len(filenames)):
  columns = [2,3,4,5,6,7,8]
  for icol in range(len(columns)):
    data = pft.Dataset(filenames[ifile],1,columns[icol])
    y = data.get_array('y')
    minval = min(minval,y.min())
    maxval = max(maxval,y.max())
    plt.plot(data.get_array('x'),y,label=data.get_name('yname'))

for ifile in range(len(filenames)):
  columns = [9,10,11,12,13,14]
  for icol in range(len(columns)):
    data = pft.Dataset(filenames[ifile],1,columns[icol])
    y = data.get_array('y')
    minval = min(minval,y.min())
    maxval = max(maxval,y.max())
    plt.plot(data.get_array('x'),y,ls='--',label=data.get_name('yname'))

#plt.xlim(0.,1.)
plt.ylim(0.5*minval,1.5*maxval)

#'best'         : 0, (only implemented for axis legends)
#'upper right'  : 1,
#'upper left'   : 2,
#'lower left'   : 3,
#'lower right'  : 4,
#'right'        : 5,
#'center left'  : 6,
#'center right' : 7,
#'lower center' : 8,
#'upper center' : 9,
#'center'       : 10,
plt.legend(loc=1)
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#      plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))
#        plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

f.subplots_adjust(hspace=0.2,wspace=0.2,
                  bottom=.1,top=.96,
                  left=.1,right=.96)

plt.show()
