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

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator

path = []
path.append('.')

files = []
files.append('pflotran-005.tec')
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(10,6))
f.suptitle("Copper Leaching",fontsize=14)
ax = f.gca(projection='3d')

ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Volume Fraction')

#plt.xlim(0.,1.)
#plt.ylim(0.,1.2)
#plt.grid(True)

#data = pft.Dataset(filenames[0],6,0)
data = pft.Dataset(filenames[0],'Jurbanite_vf',0)
X,Y = np.meshgrid(data.get_array('x'),data.get_array('y'))
Z = data.get_array('z')

nx = len(data.get_array('X'))
ny = len(data.get_array('Y'))
ZZ = np.zeros((nx,ny),'=f8')
for j in range(ny):
  for i in range(nx):
    ZZ[i][j] = Z[i+j*nx]

#surf = ax.plot_surface(X,Y,ZZ,rstride=1,cstride=1,cmap=cm.jet)
surf = ax.plot_surface(X,Y,ZZ,rstride=1,cstride=1,cmap=cm.jet,linewidth=0,antialiased=False)

#f.subplots_adjust(hspace=0.2,wspace=0.2,
#                  bottom=.12,top=.9,
#                  left=.14,right=.9)

plt.show()
