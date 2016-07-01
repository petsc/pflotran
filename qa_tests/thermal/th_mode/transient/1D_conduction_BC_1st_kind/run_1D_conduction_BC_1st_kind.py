#!/bin/env python

# *****************************************************************************
# Test Description: 
# 1D heat conduction with constant thermal conductivity in time and space, with 
# transient temperature boundary conditions at each end of a beam. Solution 
# evaluated at t=0.25 day and t=0.50 day.
#
# Kolditz, et al. (2015) Thermo-Hydro-Mechanical-Chemical Processes in 
# Fractured Porous Media: Modelling and Benchmarking, Closed Form Solutions,
# Springer International Publishing, Switzerland.
# Section 2.1.6, pg.19
# "A Transient 1D Temperature Distribution, Time-Dependent Boundary Conditions 
# of 1st Kind"
# With some modification of the parameter values given in Kolditz (2015)
#
# Author: Jenn Frederick
# Date: 06/30/2016
# *****************************************************************************
#
# Domain
# ------
# 200x1x1 hexahedral elements
# 0.5x1x1 meter elements
# L = 50 m, domain x = -L:L
#
# Parameters
# ----------
# 5.0 W/m-C thermal conductivity everywhere, constant in time
# 0.001 J/kg-C heat capacity everywhere, constant in time
# 2500 kg/m^3 density everywhere, constant in time
#
# Initial Conditions
# ------------------
# At t=0, T=0C everywhere 
#
# Boundary Conditions
# -------------------
# At x=-L=-50m, T=Tb*t, Tb=2C/day, t-units=days
# At x=L=50m, T=Tb*t, Tb=2C/day, t-units=days
#
# Time Stepping
# -------------
# Maximum timestep enforced: 0.01 day
#
# Simulation Time
# ---------------
# Simulations runs for 1 day
# Solution at t=0.25 day and t=0.50 day evaluated
# 
# Solution (time-dependent)
# -------------------------
# T(x,t) = Tb*t + ((Tb*(x^2-L^2))/(2chi)) + ...
#          ((16Tb*L^2)/(chi*pi^3))*sum_n=0^inf[((-1^n)/(2n+1)^3)*cos(((2n+1)*pi*x)/(2L))*exp(-chi*(2n+1)^2*pi^2*(t/(4L^2)))]
# chi = K/(Cp*rho)
# 
# *****************************************************************************

import numpy as np
import math
import matplotlib.pyplot as plt

# Option parameters:
plot_flag = 1  # [1 = make plots, 0 = do not make plots]
passing_crit = 2.  # [% error]

# Create the analytical solution
K = 5.0  #250.0    # [W/m-C]
Cp = 0.001     # [J/kg-C]
rho = 2500.    # [kg/m^3]
Tb = 2.0       # [C/day]
L = 50.0       # [m]
chi = K/(Cp*rho)
x_soln1 = np.linspace(0.5,49.5,100)    # [m]
x_soln2 = np.linspace(-49.5,-0.5,100)  # [m]
x_soln = np.concatenate((x_soln2,x_soln1),axis=0)
# Add boundary values to analytical solution
x_soln = np.concatenate(([-50.0],x_soln),axis=0)
x_soln = np.concatenate((x_soln,[50.0]),axis=0)
t_soln = np.array([0.0,0.25,0.50,1.0]) # [day]
T_soln = np.zeros((4,202))
for time in range(4):
  t = t_soln[time]
  i = 0
  for x in x_soln:
    T_soln[time,i] = Tb*t + ((Tb*(pow(x,2)-pow(L,2)))/(2.*chi)) 
    sum_term = 0
    # infinite sum truncated to 1000:
    for n in range(1000):
      sum_term = sum_term + (((pow(-1.,n))/(pow(((2*n)+1),3)))*math.cos((math.pi*x*((2*n)+1))/(2*L))*math.exp(-chi*pow((2*n)+1,2)*pow(math.pi,2)*(t/(4*pow(L,2))))) 
    T_soln[time,i] = T_soln[time,i] + ((16.*Tb*pow(L,2))/(chi*pow(math.pi,3)))*sum_term
    i = i + 1

# Read PFLOTRAN output file containing the temperature solution
# There are four output files
T_pflotran = np.zeros((4,202))  # [C]

# 1:=========================================================================
f = open('1D_conduction_BC_1st_kind-000.vtk', 'r')
# Temperature solution is contained starting at the 1217th line 
# and continues in chunks of 10 values, so we need to read all lines
# and concatenate them.
lines = f.readlines()
i = 1216
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 1235:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
# Add boundary values
temperature = np.concatenate(([Tb*t_soln[0]],temperature),axis=0)
temperature = np.concatenate((temperature,[Tb*t_soln[0]]),axis=0)
T_pflotran[0,:] = temperature
T_pflotran = T_pflotran.astype(np.float)
# Close PFLOTRAN output file because we are done with it
f.close()

# 2:=========================================================================
f = open('1D_conduction_BC_1st_kind-001.vtk', 'r')
# Temperature solution is contained starting at the 1217th line 
# and continues in chunks of 10 values, so we need to read all lines
# and concatenate them.
lines = f.readlines()
i = 1216
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 1235:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
# Add boundary values
temperature = np.concatenate(([Tb*t_soln[1]],temperature),axis=0)
temperature = np.concatenate((temperature,[Tb*t_soln[1]]),axis=0)
T_pflotran[1,:] = temperature
T_pflotran = T_pflotran.astype(np.float)
# Close PFLOTRAN output file because we are done with it
f.close()

# 3:=========================================================================
f = open('1D_conduction_BC_1st_kind-002.vtk', 'r')
# Temperature solution is contained starting at the 1217th line 
# and continues in chunks of 10 values, so we need to read all lines
# and concatenate them.
lines = f.readlines()
i = 1216
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 1235:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
# Add boundary values
temperature = np.concatenate(([Tb*t_soln[2]],temperature),axis=0)
temperature = np.concatenate((temperature,[Tb*t_soln[2]]),axis=0)
T_pflotran[2,:] = temperature
T_pflotran = T_pflotran.astype(np.float)
# Close PFLOTRAN output file because we are done with it
f.close()

# 4:=========================================================================
f = open('1D_conduction_BC_1st_kind-003.vtk', 'r')
# Temperature solution is contained starting at the 1217th line 
# and continues in chunks of 10 values, so we need to read all lines
# and concatenate them.
lines = f.readlines()
i = 1216
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 1235:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
# Add boundary values
temperature = np.concatenate(([Tb*t_soln[3]],temperature),axis=0)
temperature = np.concatenate((temperature,[Tb*t_soln[3]]),axis=0)
T_pflotran[3,:] = temperature
T_pflotran = T_pflotran.astype(np.float)
# Close PFLOTRAN output file because we are done with it
f.close()

# Plot the PFLOTRAN and analytical solutions
if plot_flag == 1:
  t_max = 2.0
  plt.subplot(221)
  plt.plot(x_soln,T_pflotran[0,:],'o',x_soln,T_soln[0,:])
  plt.xlabel('Distance (m)')
  plt.ylabel('Temperature (C)')
  plt.ylim([0,t_max])
  plt.title('Analytical vs. PFLOTRAN Solution, t=0day')
  plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
  plt.subplot(222)
  plt.plot(x_soln,T_pflotran[1,:],'o',x_soln,T_soln[1,:])
  plt.xlabel('Distance (m)')
  plt.ylabel('Temperature (C)')
  plt.ylim([0,t_max])
  plt.title('Analytical vs. PFLOTRAN Solution, t=0.25day')
  plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
  plt.subplot(223)
  plt.plot(x_soln,T_pflotran[2,:],'o',x_soln,T_soln[2,:])
  plt.xlabel('Distance (m)')
  plt.ylabel('Temperature (C)')
  plt.ylim([0,t_max])
  plt.title('Analytical vs. PFLOTRAN Solution, t=0.5day')
  plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
  plt.subplot(224)
  plt.plot(x_soln,T_pflotran[3,:],'o',x_soln,T_soln[3,:])
  plt.xlabel('Distance (m)')
  plt.ylabel('Temperature (C)')
  plt.ylim([0,t_max])
  plt.title('Analytical vs. PFLOTRAN Solution, t=1.0day')
  plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
  plt.show()

# Calculate error between analytical and PFLOTRAN solutions
# Shift both solutions to avoid dividing by zero
T_pflotran = T_pflotran + 0.5
T_soln = T_soln + 0.5
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
  



