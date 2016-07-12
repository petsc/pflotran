#!/bin/env python

# *****************************************************************************
# Test Description: 
# 1D heat conduction with constant thermal conductivity in time and space, with 
# transient heat flux boundary condition at right end of them beam, and steady
# temperature boundary condition at left end of the beam. Solution 
# evaluated at t=0.045 day and t=0.090 day.
#
# Kolditz, et al. (2015) Thermo-Hydro-Mechanical-Chemical Processes in 
# Fractured Porous Media: Modelling and Benchmarking, Closed Form Solutions,
# Springer International Publishing, Switzerland.
# Section 2.1.7, pg.20
# "A Transient 1D Temperature Distribution, Time-Dependent Boundary Conditions 
# of 1st and 2nd Kind"
# With some modification of the parameter values given in Kolditz (2015)
#
# Author: Jenn Frederick
# Date: 07/11/2016
# *****************************************************************************
#
# Domain
# ------
# 25x1x1 hexahedral elements
# 1x1x1 meter elements
# L = 25 m, domain x = 0:L
#
# Parameters
# ----------
# 1.16 W/m-C thermal conductivity everywhere, constant in time
# 0.01 J/kg-C heat capacity everywhere, constant in time
# 2000 kg/m^3 density everywhere, constant in time
#
# Initial Conditions
# ------------------
# At t=0, T=0C everywhere 
#
# Boundary Conditions
# -------------------
# At x=0m, T=0C
# At x=L=25m, q=0.385802*t, t-units=days
#
# Time Stepping
# -------------
# Maximum timestep enforced: 0.01 day
#
# Simulation Time
# ---------------
# Simulations runs for 0.5 days
# Solution at t=0.04 day and t=0.09 day evaluated
# 
# Solution (time-dependent)
# -------------------------
# T(x,t) = ((8*q*sqrt(chi*t^3))/K)*sum_n=0^n-inf[(i^3erfc{G1})+(i^3erfc{G2})]
# chi = K/(Cp*rho)
# G1 = ((2*n+1)*L-x)/(2*sqrt(chi*t))
# G2 = ((2*n-1)*L-x)/(2*sqrt(chi*t))
# where i^3ercf{G} =
# (2/sqrt(pi))*(int_G^inf{(s-G)^3/3! * e^(-s^2)}ds)
# 
# *****************************************************************************

import numpy as np
import math
import sympy as sym
import matplotlib.pyplot as plt
import sys
from scipy.special import erf

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
K = 1.16       # [W/m-C]
Cp = 0.01      # [J/kg-C]
rho = 2000.    # [kg/m^3]
q = 0.385802   # [C/day]
L = 25.0       # [m]
chi = K/(Cp*rho)
x_soln = np.linspace(0.5,24.5,25)               # [m]
x_soln = np.concatenate((x_soln,[25.]),axis=0)  # [m]
t_soln = np.array([0.01,0.04,0.09,0.12])        # [day]
T_soln = np.zeros((4,26))
p1 = np.zeros(26)
p2 = np.zeros(26)
i3erfc1 = np.zeros(26)
i3erfc2 = np.zeros(26)
for time in range(4):
  t = t_soln[time]*24.0*3600.0  # [sec]
  sum_term = np.zeros(26)
  # truncate infinite sum to 500
  for n in range(500):
    p1 = ((((2*n)+1)*L)-x_soln)/(2.*math.sqrt(chi*t))
    # the definite integral is below
    i3erfc1 = (1./12.)*np.sqrt(math.pi)*p1**3*erf(1.0*p1) - (1./12.)*np.sqrt(math.pi)*p1**3 + (1./3.)*p1**2*np.exp(-1.0*p1**2) + 0.5*p1*(-0.5*p1*np.exp(-1.0*p1**2) + 0.25*np.sqrt(math.pi)*erf(1.0*p1)) - 0.125*np.sqrt(math.pi)*p1 + (1./12.)*np.exp(-1.0*p1**2)
    p2 = ((((2*n)+1)*L)+x_soln)/(2.*math.sqrt(chi*t))
    i3erfc2 = (1./12.)*np.sqrt(math.pi)*p2**3*erf(1.0*p2) - (1./12.)*np.sqrt(math.pi)*p2**3 + (1./3.)*p2**2*np.exp(-1.0*p2**2) + 0.5*p2*(-0.5*p2*np.exp(-1.0*p2**2) + 0.25*np.sqrt(math.pi)*erf(1.0*p2)) - 0.125*np.sqrt(math.pi)*p2 + (1./12.)*np.exp(-1.0*p2**2)
    sum_term = sum_term + (2./math.pi)*i3erfc1 + (2./math.pi)*i3erfc2
  T_soln[time,:] = ((8*q*(1./(24.*3600.))*math.sqrt(chi*pow(t,3)))/K)*sum_term
    
# To calculate the definite integral, use sympy:
# p1, s = sym.symbols('p1 s')
# i3erfc1 = sym.integrate(((pow((s-p1),3))/(3.*2.*1.))*sym.exp(-1.*pow(s,2)),(s,p1,sym.oo))
# p2, s = sym.symbols('p2 s')
# i3erfc2 = sym.integrate(((pow((s-p2),3))/(3.*2.*1.))*sym.exp(-1.*pow(s,2)),(s,p2,sym.oo))
  
    
# Read PFLOTRAN output file containing the temperature solution
# There are four output files
T_pflotran = np.zeros((4,25))  # [C]
x_pflotran = x_soln[0:25]

# 1:=========================================================================
f = open('1D_conduction_BC_1st_2nd_kind-001.vtk', 'r')
# Temperature solution is contained starting at the 167th line 
# and continues in chunks of 10 values, so we need to read all lines
# and concatenate them.
lines = f.readlines()
i = 166
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 169:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
temperature = temperature.astype(np.float)
T_pflotran[0,:] = temperature
# Close PFLOTRAN output file because we are done with it
f.close()

# 2:=========================================================================
f = open('1D_conduction_BC_1st_2nd_kind-002.vtk', 'r')
# Temperature solution is contained starting at the 167th line 
# and continues in chunks of 10 values, so we need to read all lines
# and concatenate them.
lines = f.readlines()
i = 166
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 169:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
temperature = temperature.astype(np.float)
T_pflotran[1,:] = temperature
# Close PFLOTRAN output file because we are done with it
f.close()

# 3:=========================================================================
f = open('1D_conduction_BC_1st_2nd_kind-003.vtk', 'r')
# Temperature solution is contained starting at the 167th line 
# and continues in chunks of 10 values, so we need to read all lines
# and concatenate them.
lines = f.readlines()
i = 166
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 169:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
temperature = temperature.astype(np.float)
T_pflotran[2,:] = temperature
# Close PFLOTRAN output file because we are done with it
f.close()

# 4:=========================================================================
f = open('1D_conduction_BC_1st_2nd_kind-004.vtk', 'r')
# Temperature solution is contained starting at the 167th line 
# and continues in chunks of 10 values, so we need to read all lines
# and concatenate them.
lines = f.readlines()
i = 166
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 169:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
temperature = temperature.astype(np.float)
T_pflotran[3,:] = temperature
# Close PFLOTRAN output file because we are done with it
f.close()

# Plot the PFLOTRAN and analytical solutions
if plot_flag:
  t_max = 0.6
  plt.subplot(221)
  plt.plot(x_pflotran,T_pflotran[0,:],'o',x_soln,T_soln[0,:])
  plt.xlabel('Distance (m)')
  plt.ylabel('Temperature (C)')
  plt.ylim([0,t_max])
  plt.xlim([0,25.5])
  plt.title('Analytical vs. PFLOTRAN Solution, t=0.01day')
  plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
  plt.subplot(222)
  plt.plot(x_pflotran,T_pflotran[1,:],'o',x_soln,T_soln[1,:])
  plt.xlabel('Distance (m)')
  plt.ylabel('Temperature (C)')
  plt.ylim([0,t_max])
  plt.xlim([0,25.5])
  plt.title('Analytical vs. PFLOTRAN Solution, t=0.04day')
  plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
  plt.subplot(223)
  plt.plot(x_pflotran,T_pflotran[2,:],'o',x_soln,T_soln[2,:])
  plt.xlabel('Distance (m)')
  plt.ylabel('Temperature (C)')
  plt.ylim([0,t_max])
  plt.xlim([0,25.5])
  plt.title('Analytical vs. PFLOTRAN Solution, t=0.09day')
  plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
  plt.subplot(224)
  plt.plot(x_pflotran,T_pflotran[3,:],'o',x_soln,T_soln[3,:])
  plt.xlabel('Distance (m)')
  plt.ylabel('Temperature (C)')
  plt.ylim([0,t_max])
  plt.xlim([0,25.5])
  plt.title('Analytical vs. PFLOTRAN Solution, t=0.12day')
  plt.legend(('PFLOTRAN','analytical'),'best',numpoints=1)
  plt.show()

# Calculate error between analytical and PFLOTRAN solutions
# Shift both solutions to avoid dividing by zero
T_pflotran = T_pflotran + 0.5
T_soln = T_soln + 0.5
percent_error = 100.0*(T_pflotran-T_soln[:,0:25])/T_soln[:,0:25]
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
  



