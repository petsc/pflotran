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

files = pft.get_default_tec_filenames(range(1,6))
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(6,6))
plt.subplot(1,1,1)
f.suptitle("1D Calcite Variably-Saturated Flow",fontsize=16)
plt.xlabel('pH')
plt.ylabel('Z [m]')

#plt.xlim(0.,1.)
plt.xlim(4.8,8.2)
#plt.grid(True)

for ifile in range(len(filenames)):
  data = pft.Dataset(filenames[ifile],6,3)
  plt.plot(data.get_array('x'),data.get_array('y'),label=data.title)

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
plt.legend(loc=3,title='Time [y]')
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

plt.twiny()
plt.xlabel('Saturation [-]')
plt.xlim(0.,1.05)
data = pft.Dataset(filenames[ifile],5,3)
plt.plot(data.get_array('x'),data.get_array('y'), \
         label='Saturation',color='black',ls='--')
plt.legend(loc=1)
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))

f.subplots_adjust(hspace=0.2,wspace=0.2,
                  bottom=.12,top=.85,
                  left=.12,right=.9)

plt.show()
