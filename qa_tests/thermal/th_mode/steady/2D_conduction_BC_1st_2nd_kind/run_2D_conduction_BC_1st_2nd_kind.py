#!/bin/env python

# *****************************************************************************
# Test Description: 
# 2D heat conduction with constant thermal conductivity in time and space, with 
# steady temperature boundary conditions at three ends of plate, and a steady
# heat flux boundary condition at the fourth end. Steady state solution.
#
# Kolditz, et al. (2015) Thermo-Hydro-Mechanical-Chemical Processes in 
# Fractured Porous Media: Modelling and Benchmarking, Closed Form Solutions,
# Springer International Publishing, Switzerland.
# Section 2.1.4, pg.17
# "A 2D Steady-State Temperature Distribution, Boundary Conditions of 1st and
# 2nd Kind"
#
# Author: Jenn Frederick
# Date: 06/29/2016
# *****************************************************************************
#
# Domain
# ------
# 20x10x1 hexahedral elements
# 0.1x0.1x1 meter elements (2x1x1 m domain, L=1m)
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
# At x=0,y=0:L,  dT/dx=T(t=0)/L, or q=-1 W/m2
# At x=0:2L,y=0, T=(T(t=0)/L)*x
# At x=2L,y=0:L, T=(T(t=0)/L)*(2L+2y)
# At x=0:2L,y=L, T=(T(t=0)/L)*(2L+x)
#
# Simulation Time
# ---------------
# Solution at t = 100 yr assumed steady state solution
# 
# Solution (steady state)
# -----------------------
# T(x,y) = (T(t=0)/L)*(x + 2y)
#
# *****************************************************************************

import numpy as np
import matplotlib.pyplot as plt
import sys

# Default options:
plot_flag = False
print_error = True
passing_crit = 2.  # [% error]

# Option parameters:
options = sys.argv[1:]
num_options = len(options)
for n in options:
  if n == 'print_error=false':
    print_error = False
  if n == 'print_error=true':
    print_error = True
  if n == 'plot_flag=false':
    plot_flag = False
  if n == 'plot_flag=true':
    plot_flag = True

# Create the analytical solution
L = 1.   # [m]
T0 = 1.  # [C]
T_soln = np.zeros((10,20))  # [C]
x_soln = np.linspace(0.05,1.95,20)  # [m]
y_soln = np.linspace(0.05,0.95,10)  # [m]
x = -1
for i in x_soln:
  x = x + 1
  y = -1
  for j in y_soln:
    y = y + 1
    T_soln[y,x] = (T0/L)*(i + (2*j))

# Read PFLOTRAN output file containing the temperature solution
f = open('2D_conduction_BC_1st_2nd_kind-001.vtk', 'r')

# Temperature solution is contained starting at the 875th line 
# and continues in chunks of 10 values, so we need to read 2 lines
# at a time (10*2 = 20)
T_pflotran = np.zeros((10,20))  # [C]
lines = f.readlines()
i = 874
j = -1
while i < 893:
  lines[i] = lines[i].strip()
  temperature = np.array(lines[i].split())
  lines[i+1] = lines[i+1].strip()
  temperature = np.concatenate((temperature,np.array(lines[i+1].split())),axis=0)
  i = i + 2
  j = j + 1
  T_pflotran[j,:] = temperature
x_pflotran = np.linspace(0.05,1.95,20)
y_pflotran = np.linspace(0.05,0.95,10)
T_pflotran = T_pflotran.astype(np.float)
  
# Close PFLOTRAN output file because we are done with it
f.close()

# Plot the PFLOTRAN and analytical solutions
if plot_flag:  
  X,Y = np.meshgrid(x_soln,y_soln)
  # temperature values to contour against and compare visually
  levels = np.linspace(0.,4.,11)
  plt.contourf(X,Y,T_soln,levels,alpha=0.75)
  C = plt.contour(X,Y,T_pflotran,levels,colors='black',linewidth=0.5)
  plt.clabel(C,inline=True,fontsize=10)
  plt.xlabel('Distance (m)')
  plt.ylabel('Distance (m)')
  plt.title('Analytical (fill) vs. PFLOTRAN (contours) Solution, Temperature')
  plt.show()

# Calculate error between analytical and PFLOTRAN solutions
# Shift solution so we avoid dividing by zero
T_pflotran = T_pflotran + 0.50
T_soln = T_soln + 0.50
percent_error = 100.0*(T_pflotran-T_soln)/T_soln
max_percent_error = np.nanmax(abs(percent_error))
if print_error:
  print 'Percent Error (temperature):'
  print percent_error
print 'Maximum Error:'
print max_percent_error

# Decide if test passes
if abs(max_percent_error) > passing_crit:
  print '-- Test FAIL --\n'
else: 
  print '-- Test PASS --\n'
  



