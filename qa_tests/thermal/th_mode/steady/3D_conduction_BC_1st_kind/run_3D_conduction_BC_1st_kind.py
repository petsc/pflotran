#!/bin/env python

# *****************************************************************************
# Test Description: 
# 3D heat conduction with constant thermal conductivity in time and space, with 
# steady temperature boundary conditions at each end of plate. Steady state
# solution.
#
# Kolditz, et al. (2015) Thermo-Hydro-Mechanical-Chemical Processes in 
# Fractured Porous Media: Modelling and Benchmarking, Closed Form Solutions,
# Springer International Publishing, Switzerland.
# Section 2.1.5, pg.18
# "A 3D Steady-State Temperature Distribution"
#
# Author: Jenn Frederick
# Date: 06/28/2016
# *****************************************************************************
#
# Domain
# ------
# 10x10x10 cubic elements
# 0.1x0.1x0.1 meter elements (1x1x1 m cube domain)
#
# Parameters
# ----------
# 1 W/m-C thermal conductivity everywhere, constant in time
#
# Initial Conditions
# ------------------
# At t=0, T=1.0C everywhere 
#
# Boundary Conditions
# -------------------
# At x=0 face, T=T(t=0)*(0+(y/L)+(z/L))
# At y=0 face, T=T(t=0)*((x/L)+0+(z/L))
# At z=0 face, T=T(t=0)*((x/L)+(y/L)+0)
# At x=L face, T=T(t=0)*(1+(y/L)+(z/L))
# At y=L face, T=T(t=0)*((x/L)+1+(z/L))
# At z=L face, T=T(t=0)*((x/L)+(y/L)+1)
#
# Simulation Time
# ---------------
# Solution at t = 400 yr assumed steady state solution
# 
# Solution (steady state)
# -----------------------
# T(x,y) = T=T(t=0)*((x/L)+(y/L)+(z/L))
#
# *****************************************************************************

import numpy as np
import matplotlib.pyplot as plt

# Option parameters:
plot_flag = 0  # [1 = make plots, 0 = do not make plots]
passing_crit = 2.  # [% error]

# Create the analytical solution
L = 1.   # [m]
T0 = 1.  # [C]
T_soln = np.zeros((10,10,10))  # [C]
x_soln = np.linspace(0.05,0.95,10)  # [m]
y_soln = np.linspace(0.05,0.95,10)  # [m]
z_soln = np.linspace(0.05,0.95,10)  # [m]
x = -1
for i in x_soln:
  x = x + 1
  y = -1
  for j in y_soln:
    y = y + 1
    z = -1
    for k in z_soln:
      z = z + 1
      T_soln[x,y,z] = T0*((i/L)+(j/L)+(k/L))

# Read PFLOTRAN output file containing the temperature solution
f = open('3D_conduction_BC_1st_kind-001.vtk', 'r')

# Temperature solution is contained starting at the 3344th line 
# and continues in chunks of 10 values.
T_pflotran = np.zeros((10,10,10))  # [C]
lines = f.readlines()
i = 3343
k = 0
while i < 3443:
  lines[i] = lines[i].strip()
  temperature = np.array(lines[i].split())
  T_pflotran[:,0,k] = temperature
  lines[i+1] = lines[i+1].strip()
  temperature = np.array(lines[i+1].split())
  T_pflotran[:,1,k] = temperature
  lines[i+2] = lines[i+2].strip()
  temperature = np.array(lines[i+2].split())
  T_pflotran[:,2,k] = temperature
  lines[i+3] = lines[i+3].strip()
  temperature = np.array(lines[i+3].split())
  T_pflotran[:,3,k] = temperature
  lines[i+4] = lines[i+4].strip()
  temperature = np.array(lines[i+4].split())
  T_pflotran[:,4,k] = temperature
  lines[i+5] = lines[i+5].strip()
  temperature = np.array(lines[i+5].split())
  T_pflotran[:,5,k] = temperature
  lines[i+6] = lines[i+6].strip()
  temperature = np.array(lines[i+6].split())
  T_pflotran[:,6,k] = temperature
  lines[i+7] = lines[i+7].strip()
  temperature = np.array(lines[i+7].split())
  T_pflotran[:,7,k] = temperature
  lines[i+8] = lines[i+8].strip()
  temperature = np.array(lines[i+8].split())
  T_pflotran[:,8,k] = temperature
  lines[i+9] = lines[i+9].strip()
  temperature = np.array(lines[i+9].split())
  T_pflotran[:,9,k] = temperature
  i = i + 10
  k = k + 1

x_pflotran = np.linspace(0.05,0.95,10)  # [m]
y_pflotran = np.linspace(0.05,0.95,10)  # [m]
z_pflotran = np.linspace(0.05,0.95,10)  # [m]
T_pflotran = T_pflotran.astype(np.float)
  
# Close PFLOTRAN output file because we are done with it
f.close()

# Plot the PFLOTRAN and analytical solutions
if plot_flag == 1:
  # temperature values to contour against and compare visually
  levels = np.linspace(0.,3.,31)
  slice_ind = 6
  X,Y = np.meshgrid(x_soln,y_soln)
  plt.contourf(X,Y,T_soln[:,:,slice_ind],levels)
  C = plt.contour(X,Y,T_pflotran[:,:,slice_ind],levels,colors='black',linewidth=0.5)
  plt.clabel(C,inline=True,fontsize=10)
  plt.xlabel('Distance-x (m)')
  plt.ylabel('Distance-y (m)')
  plt.title('Analytical (fill) vs. PFLOTRAN (contour) Solution, Temperature Slice')
  plt.show()
  slice_ind = 4
  plt.contourf(X,Y,T_soln[:,slice_ind,:],levels)
  C = plt.contour(X,Y,T_pflotran[:,slice_ind,:],levels,colors='black',linewidth=0.5)
  plt.clabel(C,inline=True,fontsize=10)
  plt.xlabel('Distance-x (m)')
  plt.ylabel('Distance-z (m)')
  plt.title('Analytical (fill) vs. PFLOTRAN (contour) Solution, Temperature Slice')
  plt.show()

# Calculate error between analytical and PFLOTRAN solutions
# Shift solution so we avoid dividing by zero
T_pflotran = T_pflotran + 0.50
T_soln = T_soln + 0.50
percent_error = 100.0*(T_pflotran-T_soln)/T_soln
max_percent_error = np.nanmax(abs(percent_error))
print 'Percent Error (temperature):'
print percent_error
print 'Maximum Error:'
print max_percent_error

# Decide if test passes
if abs(max_percent_error) > passing_crit:
  print '-- Test FAIL --\n'
else: 
  print '-- Test PASS --\n'
  



