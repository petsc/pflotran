#!/bin/env python

# *****************************************************************************
# Test Description: 
# 1D heat conduction with constant thermal conductivity in time and space, with 
# constant temperature boundary conditions at each end of beam. Steady state
# solution.
#
# Kolditz, et al. (2015) Thermo-Hydro-Mechanical-Chemical Processes in 
# Fractured Porous Media: Modelling and Benchmarking, Closed Form Solutions,
# Springer International Publishing, Switzerland.
# Section 2.1.1, pg.14
# "A 1D Steady-State Temperature Distribution, Boundary Conditions of 1st Kind"
#
# Author: Jenn Frederick
# Date: 06/27/2016
# *****************************************************************************
#
# Domain
# ------
# 10x1x1 cubic elements
# 10 meter elements
#
# Parameters
# ----------
# 1 W/m-C thermal conductivity everywhere, constant in time
#
# Initial Conditions
# ------------------
# At t=0, T=1.5C everywhere 
#
# Boundary Conditions
# -------------------
# At x=0, T=1C
# At x=L=100m, T=2C
#
# Simulation Time
# ---------------
# Solution at t = 500 yr assumed steady state solution
# 
# Solution (steady state)
# -----------------------
# T(x) = a*x + b
# a = 1/100
# b = 1
# T(x) = x/100 + 1
#
# *****************************************************************************

import numpy as np
import matplotlib.pyplot as plt

# Option parameters:
plot_flag = 0  # [1 = make plots, 0 = do not make plots]
passing_crit = 2.  # [% error]

# Create the analytical solution
x_soln = np.array([0.,5.,15.,25.,35.,45.,55.,65.,75.,85.,95.,100.])
T_soln = np.array(x_soln/100. + 1.)

# Read PFLOTRAN output file containing the temperature solution
f = open('1D_conduction_BC_1st_kind-001.vtk', 'r')

# Temperature solution is contained in the 77th line 
for i, line in enumerate(f):
  if i == 76:
    temperature = line.strip() 
  elif i > 76:
    break

# Close PFLOTRAN output file because we are done with it
f.close()

# Extract the temperature from the output string into an array
x_pflotran = np.array([5.,15.,25.,35.,45.,55.,65.,75.,85.,95.])
T_pflotran = np.array(temperature.split())
T_pflotran = T_pflotran.astype(np.float)

# Add boundary temperature values
x_pflotran = np.concatenate(([0.],x_pflotran),axis=0)
T_pflotran = np.concatenate(([1.],T_pflotran),axis=0)
x_pflotran = np.concatenate((x_pflotran,[100.]),axis=0)
T_pflotran = np.concatenate((T_pflotran,[2.]),axis=0)

# Plot the PFLOTRAN and analytical solutions
if plot_flag == 1:
  plt.plot(x_pflotran,T_pflotran,'o',x_soln,T_soln)
  plt.xlabel('Distance (m)')
  plt.ylabel('Temperature (C)')
  plt.title('Analytical vs. PFLOTRAN Solution')
  plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
  plt.show()

# Calculate error between analytical and PFLOTRAN solutions
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
  



