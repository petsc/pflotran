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
# 0.01 J/kg-C heat capacity everywhere, constant in time
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
K = 5.0              # [W/m-C]
Cp = 0.01            # [J/kg-C]
rho = 2500.          # [kg/m^3]
Tb_day = 2.0         # [C/day]
Tb_sec = 2.0/(24.0*3600.0)   # [C/sec]
L = 50.0             # [m]
chi = K/(Cp*rho)
x_soln1 = np.linspace(0.5,49.5,2*L)    # [m]
x_soln2 = np.linspace(-49.5,-0.5,2*L)  # [m]
x_soln = np.concatenate((x_soln2,x_soln1),axis=0)
# Add boundary values to analytical solution
x_soln = np.concatenate(([-50.0],x_soln),axis=0)
x_soln = np.concatenate((x_soln,[50.0]),axis=0)
t_soln = np.array([0.0,0.25,0.50,0.75]) # [day]
T_soln = np.zeros((4,(4*L+2)))

for time in range(4):
  t = t_soln[time]*(24.0*3600.0)   # [sec]
  T_soln[time,:] = Tb_sec*t + ((Tb_sec*(pow(x_soln,2)-pow(L,2)))/(2.*chi))
  sum_term = np.zeros(4*L+2)
  sum_term_old = np.zeros(4*L+2)
  n = 0
  epsilon = 1.0
  # infinite sum truncated once max(epsilon) < 1e-10:
  while epsilon > 1e-10:
    sum_term_old = sum_term
    sum_term = sum_term_old + (((pow(-1.,n))/(pow(((2*n)+1),3)))*np.cos((math.pi*x_soln*((2*n)+1))/(2*L))*np.exp(-chi*pow((2*n)+1,2)*pow(math.pi,2)*(t/(4*pow(L,2))))) 
    epsilon = np.max(np.abs(sum_term_old-sum_term))
    n = n + 1
  T_soln[time,:] = T_soln[time,:] + ((16.*Tb_sec*pow(L,2))/(chi*pow(math.pi,3)))*sum_term

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
temperature = np.concatenate(([Tb_day*t_soln[0]],temperature),axis=0)
temperature = np.concatenate((temperature,[Tb_day*t_soln[0]]),axis=0)
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
temperature = np.concatenate(([Tb_day*t_soln[1]],temperature),axis=0)
temperature = np.concatenate((temperature,[Tb_day*t_soln[1]]),axis=0)
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
temperature = np.concatenate(([Tb_day*t_soln[2]],temperature),axis=0)
temperature = np.concatenate((temperature,[Tb_day*t_soln[2]]),axis=0)
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
temperature = np.concatenate(([Tb_day*t_soln[3]],temperature),axis=0)
temperature = np.concatenate((temperature,[Tb_day*t_soln[3]]),axis=0)
T_pflotran[3,:] = temperature
T_pflotran = T_pflotran.astype(np.float)
# Close PFLOTRAN output file because we are done with it
f.close()

# Plot the PFLOTRAN and analytical solutions
t_max = 2.0
plt.figure(figsize=(10,10))
plt.subplot(221)
plt.plot(x_soln,T_pflotran[0,:],'o',x_soln,T_soln[0,:])
plt.xlabel('Distance (m)')
plt.ylabel('Temperature (C)')
plt.ylim([0,t_max])
plt.title('Solution, t=0day')
plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
plt.subplot(222)
plt.plot(x_soln,T_pflotran[1,:],'o',x_soln,T_soln[1,:])
plt.xlabel('Distance (m)')
plt.ylabel('Temperature (C)')
plt.ylim([0,t_max])
plt.title('Solution, t=0.25day')
plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
plt.subplot(223)
plt.plot(x_soln,T_pflotran[2,:],'o',x_soln,T_soln[2,:])
plt.xlabel('Distance (m)')
plt.ylabel('Temperature (C)')
plt.ylim([0,t_max])
plt.title('Solution, t=0.5day')
plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
plt.subplot(224)
plt.plot(x_soln,T_pflotran[3,:],'o',x_soln,T_soln[3,:])
plt.xlabel('Distance (m)')
plt.ylabel('Temperature (C)')
plt.ylim([0,t_max])
plt.title('Solution, t=0.75day')
plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
plt.savefig('comparison_plot.png')
  
if plot_flag:
  plt.show()

# Calculate error between analytical and PFLOTRAN solutions
# Shift both solutions to avoid dividing by zero
T_pflotran = T_pflotran + 0.5
T_soln = T_soln + 0.5
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
  



