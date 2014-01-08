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

files = []
files.append('case3-obs-0.tec')
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(6,6))
plt.subplot(1,1,1)
f.suptitle("Case 3_0",fontsize=14)
plt.xlabel('Time [y]')
plt.ylabel('Concentration [-]')

#plt.xlim(0.,1.)
#plt.ylim(0.,1.)
#plt.grid(True)

data = pft.Dataset(filenames[0],1,2)
times = data.get_array('x')
times_y = [0.]*len(times)
for i in range(len(times)):
  times_y[i] = times[i]/3600./24./365.
conc_M = data.get_array('y')

observed_concentration = [0.]*len(times)
for i in range(len(times)):
  observed_concentration[i] = conc_M[i]*1000.
plt.plot(times_y,observed_concentration,label='observed')

analytical_concentration = [0.]*len(times)
for i in range(len(times)):
#  analytical_concentration[i] = 4.1631973e-9*times[i]
  analytical_concentration[i] = 5.1894e-3*(1.-math.exp(-8.0225e-7*times[i]))
plt.plot(times_y,analytical_concentration,ls='--',color='red',label='analytical')

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
plt.legend(loc=2,title='Solution')
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))
#plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

f.subplots_adjust(hspace=0.2,wspace=0.2,
                  bottom=.12,top=.9,
                  left=.16,right=.9)

plt.show()
