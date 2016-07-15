import sys
from h5py import *
import numpy as np
import math
import matplotlib.pyplot as plt

filename = 'dataset.h5'
h5file = File(filename,mode='w')

T0 = 1.        # [C]
L = 100.       # [m]
dx = 1.0       # [m]
dy = 1.0       # [m]
dz = 1.0       # [m]
K = 0.5787037  # [W/m-C]
rho = 2000.    # [kg/m^3]
Cp = 0.01      # [J/kg-C]
chi = K/(rho*Cp)

# 2D Surface:
# -----------------------------------------------------------------------------
# Temperature initial condition z=0 face
h5grp = h5file.create_group('initial')
h5grp.attrs['Dimension'] = np.string_('XY')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [dx,dy]
# Location of origin
h5grp.attrs['Origin'] = [0.,0.]
# Load the dataset values
nx = L*dx + 1
ny = L*dy + 1
rarray = np.zeros((nx,ny),'=f8')

# Use the analytical solution at t=0 to get the initial condition:
print 'Creating dataset . . .'
t = 0.

y = np.linspace(0,L,ny)
T2y = np.zeros(int(ny))
sum_term_y = np.zeros(int(ny))
# truncate infinite sum to 5000
for n in range(5000):
  n = n + 1
  sum_term_y = sum_term_y + (np.cos(n*math.pi*y/L)*np.exp(-chi*pow(n,2)*pow(math.pi,2)*t/pow(L,2))*(80./(3.*pow((n*math.pi),2)))*np.cos(n*math.pi/2.)*np.sin(n*math.pi/4.)*np.sin(3.*n*math.pi/20.))
T2y = 0.5 + sum_term_y

x = np.linspace(0,L,nx)
T1x = np.zeros(int(nx))
sum_term_x = np.zeros(int(nx))
# truncate infinite sum to 5000
for n in range(5000):
  n = n + 1
  sum_term_x = sum_term_x + (np.sin(n*math.pi*x/L)*np.exp(-chi*pow(n,2)*pow(math.pi,2)*t/pow(L,2))*(80./(3.*pow((n*math.pi),2)))*np.sin(n*math.pi/2.)*np.sin(n*math.pi/4.)*np.sin(3.*n*math.pi/20.))
T1x = sum_term_x

for i in range(int(nx)):
  for j in range(int(ny)):
    rarray[i][j] = T0*T1x[i]*T2y[j] + 0.10
    if rarray[i][j] < 1.0e-5:
      rarray[i][j] = 0.
h5dset = h5grp.create_dataset('Data', data=rarray)



