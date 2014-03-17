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

files = pft.get_tec_filenames('tracer_1D_MC-sec-rank0-obs0',range(5))
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(6,6))
plt.subplot(1,1,1)
f.suptitle("1D Tracer - Matrix Concentration Profiles",fontsize=16)
plt.xlabel('Time [s]')
plt.ylabel('Concentration [mol/L]')


for ifile in range(len(filenames)):
  data = pft.Dataset(filenames[ifile],1,2)
  plt.plot(data.get_array('x'),data.get_array('y'),label=data.title)

plt.legend(loc=3,title='Time [s]')

plt.show()
