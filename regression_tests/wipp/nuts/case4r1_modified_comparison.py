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
import analytical_solutions as analytical

path = []
path.append('.')

files = []
files.append('case4r1_modified-002.tec')
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(6,6))
plt.subplot(1,1,1)
f.suptitle("NUTS Case 4 Run 1 (modified)",fontsize=14)
plt.xlabel('X [m]')
plt.ylabel('C [-]')

plt.xlim(0.,30.)
#plt.ylim(0.,1.)
#plt.grid(True)

data = pft.Dataset(filenames[0],1,4)

x = data.get_array('x')
y = data.get_array('y')

plt.plot(x,y,label='PFLOTRAN')

n = len(x)
n = 30
x = [0.]*n
y_ogata = [0.]*n
y_demarsily_no_rxn = [0.]*n
y_demarsily = [0.]*n
y_bear = [0.]*n

time = 8.64e6
porosity = 0.2
tortuosity = .1
vdarcy = 1.e-7
saturation = 1.
dispersivity = 10.
diffusion = 7.5e-6
half_life = 1.e6
anal = analytical.AnalyticalSolution(0,1.,vdarcy,saturation,diffusion,
                                     dispersivity,tortuosity,porosity,
                                     1.,half_life)
for i in range(n):
  x[i] = (i+0.5)*1.
  y_bear[i] = anal.bear(x[i],time)
  y_demarsily_no_rxn[i] = anal.de_marsily_no_reaction(x[i],time)
  y_demarsily[i] = anal.de_marsily(x[i],time)
  y_ogata[i] = anal.ogata_banks(x[i],time)

plt.plot(x,y_bear,label='Bear',ls='-.')
plt.plot(x,y_demarsily,label='de Marsily',ls='--')

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
plt.legend(loc=1,title='Solution')
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
