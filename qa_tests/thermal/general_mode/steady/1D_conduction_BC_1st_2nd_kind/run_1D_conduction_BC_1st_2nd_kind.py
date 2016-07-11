#!/bin/env python

# *****************************************************************************
# Test Description: 
# 1D heat conduction with two thermal conductivity values in space that are
# constant in time, with a constant temperature boundary condition at one end
# of the beam, and a constant heat flux at the other end. Steady state
# solution.
#
# Kolditz, et al. (2015) Thermo-Hydro-Mechanical-Chemical Processes in 
# Fractured Porous Media: Modelling and Benchmarking, Closed Form Solutions,
# Springer International Publishing, Switzerland.
# Section 2.1.2, pg.15
# "A 1D Steady-State Temperature Distribution, Boundary Conditions of 1st and
# 2nd Kind"
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
# 100 W/m-C thermal conductivity x < 2L/5, constant in time, K1
# 300 W/m-C thermal conductivity x > 2L/5, constant in time, K2
#
# Initial Conditions
# ------------------
# At t=0, T=1.0C everywhere 
#
# Boundary Conditions
# -------------------
# At x=0, T=1C
# At x=L=100m, q=-1.5 W/m2
#
# Simulation Time
# ---------------
# Solution at t = 100 yr assumed steady state solution
# 
# Solution (steady state)
# -----------------------
# T(x) | x <= 2L/5 = -(q/K1)*x + T(x=0)
# T(x) | x > 2L/5 = -(q/K2)*x + q*(2L/5)*(1/K2 - 1/K1)
#
# *****************************************************************************

import numpy as np
import matplotlib.pyplot as plt

# Option parameters:
plot_flag = 0
passing_crit = 2.

# Create the analytical solution
x_soln = np.array([0.,5.,15.,25.,35.,45.,55.,65.,75.,85.,95.])  # [m]
K1 = 100.  # [W/m-C]
K2 = 300.  # [W/m-C]
q = -1.5   # [W/m^2]
L = 100.   # [m]
T_soln = np.zeros(11)
k = -1
for j in x_soln:
  k = k + 1
  if j <= (2.*L/5.):
    T_soln[k] = np.array( -((q/K1)*j) + 1. )
  else:
    T_soln[k] = np.array( -((q/K2)*j) + 1. + (q*(2.*L/5.)*((1/K2)-(1/K1))) )

# Read PFLOTRAN output file containing the temperature solution
f = open('1D_conduction_BC_1st_2nd_kind-001.vtk', 'r')

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
  



