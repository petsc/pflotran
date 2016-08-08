#!/bin/env python

# *****************************************************************************
# Test Description: 
# 2D heat conduction with constant thermal conductivity in time and space, with 
# steady temperature boundary conditions at two ends of plate, and a steady
# heat flux boundary condition at the other two ends. Transient solution 
# sampled and compared at t=0.02 day and t=0.04 day.
#
# Kolditz, et al. (2015) Thermo-Hydro-Mechanical-Chemical Processes in 
# Fractured Porous Media: Modelling and Benchmarking, Closed Form Solutions,
# Springer International Publishing, Switzerland.
# Section 2.1.9, pg.23
# "A Transient 2D Temperature Distribution, Non-Zero Initial Temperature,
# Boundary Conditions of 1st and 2nd Kind"
#
# Author: Jenn Frederick
# Date: 07/13/2016
# *****************************************************************************
#
# Domain
# ------
# 100x100x1 hexahedral elements
# 1x1x1 meter elements (100x100x1 m domain, L=100m)
#
# Parameters
# ----------
# 0.5787037 W/m-C thermal conductivity everywhere, constant in time
# 0.01 J/kg-C heat capacity everywhere, constant in time
# 2000 kg/m^3 density everywhere, constant in time
#
# Initial Conditions
# ------------------
# At t=0, T=1C*f(x)*f(y) everywhere
# f(x) =
# 0C for [0 <= x <= L/10]
# ((10/3L)x - (1/3))C for [L/10 <= x <= 4L/10]
# 1C for [4L/10 <= x <= 6L/10]
# (3 - (10/3L)x)C for [6L/10 <= x <= 9L/10]
# 0C for [9L/10 <= x <= L]
# f(y) = same function
#
# Boundary Conditions
# -------------------
# At x=0,y=0:L,  T=0C
# At x=0:L,y=0,  q=0 W/m2
# At x=L,y=0:L,  T=0C
# At x=0:L,y=L,  q=0 W/m2
#
# Time Stepping
# -------------
# Maximum timestep enforced: 0.005 day
#
# Simulation Time
# ---------------
# Simulations runs for 0.10 days
# Solution at t=0.05 day and t=0.10 day evaluated
# 
# Solution (time-dependent)
# -------------------------
# T1(x,t) = sum_n=1^n=inf{sin(n*pi*x/L)*exp(-chi*n^2*pi^2*t/L^2)*...
#           (80/(3*(n*pi)^2))*sin(n*pi/2)*sin(n*pi/4)*sin(3*n*pi/20)}
# T2(y,t) = 0.5 + sum_n=1^n=inf{cos(n*pi*y/L)*exp(-chi*n^2*pi^2*t/L^2)*...
#           (80/(3*(n*pi)^2))*cos(n*pi/2)*sin(n*pi/4)*sin(3*n*pi/20)}
# chi = K/(rho*Cp) 
# T(x,y,t) = T0*T1(x,t)*T2(y,t)
# T0 = 1C
#
# *****************************************************************************

import numpy as np
import matplotlib.pyplot as plt
import sys
import math

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
L = 100.   # [m]
T0 = 1.  # [C]
nx = 50.
ny = 50.
K = 0.5787037  # [W/m-C]
rho = 2000.    # [kg/m^3]
Cp = 0.01      # [J/kg-C]
chi = K/(rho*Cp)
T_soln = np.zeros((4,ny,nx))          # [C]
t_soln = np.array([0.0,0.02,0.04,0.10]) # [day]
mid = L/nx/2.
x_soln = np.linspace(0.+mid,L-mid,nx)   # [m]
y_soln = np.linspace(0.+mid,L-mid,ny)   # [m]
T1x = np.zeros(int(nx))
T2y = np.zeros(int(ny))
k = 0
for t in t_soln:
  t = t*(24.*3600.)
  sum_term_y = np.zeros(int(ny))
  # truncate infinite sum to 500
  for n in range(500):
    n = n + 1
    sum_term_y = sum_term_y + (np.cos(n*math.pi*y_soln/L)*np.exp(-chi*pow(n,2)*pow(math.pi,2)*t/pow(L,2))*(80./(3.*pow((n*math.pi),2)))*np.cos(n*math.pi/2.)*np.sin(n*math.pi/4.)*np.sin(3.*n*math.pi/20.))
  T2y = 0.5 + sum_term_y
  sum_term_x = np.zeros(int(nx))
  # truncate infinite sum to 500
  for n in range(500):
    n = n + 1
    sum_term_x = sum_term_x + (np.sin(n*math.pi*x_soln/L)*np.exp(-chi*pow(n,2)*pow(math.pi,2)*t/pow(L,2))*(80./(3.*pow((n*math.pi),2)))*np.sin(n*math.pi/2.)*np.sin(n*math.pi/4.)*np.sin(3.*n*math.pi/20.))
  T1x = sum_term_x
  for i in range(int(nx)):
    for j in range(int(ny)):
      T_soln[k,j,i] = T0*T1x[i]*T2y[j]
  k = k + 1

# Read PFLOTRAN output file containing the temperature solution
# There are four output files
T_pflotran = np.zeros((4,ny,nx))   # [C]

# 1:=========================================================================
f = open('2D_conduction_BC_1st_2nd_kind-000.vtk', 'r')
# Temperature solution is contained starting at the 10215th line 
# and continues in chunks of 10 values, and we need to read all lines
lines = f.readlines()
i = 10214
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 10464:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
temperature = temperature.astype(np.float)
for i in range(int(nx)):
  for j in range(int(ny)):
    pt = ny*(i)+(j)
    T_pflotran[0,i,j] = temperature[pt]
#Close PFLOTRAN output file because we are done with it
f.close()
    
# 2:=========================================================================
f = open('2D_conduction_BC_1st_2nd_kind-001.vtk', 'r')
# Temperature solution is contained starting at the 10215th line 
# and continues in chunks of 10 values, and we need to read all lines
lines = f.readlines()
i = 10214
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 10464:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
temperature = temperature.astype(np.float)
for i in range(int(nx)):
  for j in range(int(ny)):
    pt = ny*(i)+(j)
    T_pflotran[1,i,j] = temperature[pt]
#Close PFLOTRAN output file because we are done with it
f.close()

# 3:=========================================================================
f = open('2D_conduction_BC_1st_2nd_kind-002.vtk', 'r')
# Temperature solution is contained starting at the 10215th line 
# and continues in chunks of 10 values, and we need to read all lines
lines = f.readlines()
i = 10214
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 10464:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
temperature = temperature.astype(np.float)
for i in range(int(nx)):
  for j in range(int(ny)):
    pt = ny*(i)+(j)
    T_pflotran[2,i,j] = temperature[pt]
#Close PFLOTRAN output file because we are done with it
f.close()

# 4:=========================================================================
f = open('2D_conduction_BC_1st_2nd_kind-003.vtk', 'r')
# Temperature solution is contained starting at the 10215th line 
# and continues in chunks of 10 values, and we need to read all lines
lines = f.readlines()
i = 10214
# Get first line
lines[i] = lines[i].strip()
temperature = np.array(lines[i].split())
while i < 10464:
  i = i + 1
  lines[i] = lines[i].strip()
  temperature = np.concatenate((temperature,np.array(lines[i].split())),axis=0)
temperature = temperature.astype(np.float)
for i in range(int(nx)):
  for j in range(int(ny)):
    pt = ny*(i)+(j)
    T_pflotran[3,i,j] = temperature[pt]
#Close PFLOTRAN output file because we are done with it
f.close()

# Plot the PFLOTRAN and analytical solutions
if plot_flag:  
  X,Y = np.meshgrid(x_soln,y_soln)
  plt.subplot(221)
  # temperature values to contour against and compare visually
  levels = np.linspace(0.,1.,11)
  plt.contourf(X,Y,T_soln[0,:,:],levels,alpha=0.75)
  C = plt.contour(X,Y,T_pflotran[0,:,:],levels,colors='black',linewidth=0.5)
  plt.clabel(C,inline=True,fontsize=10)
  plt.xlabel('Distance (m)')
  plt.ylabel('Distance (m)')
  plt.title('Analytical (fill) vs. PFLOTRAN (contours) Solution, Temperature, t=0d')
  plt.subplot(222)
  # temperature values to contour against and compare visually
  levels = np.linspace(0.,1.,11)
  plt.contourf(X,Y,T_soln[1,:,:],levels,alpha=0.75)
  C = plt.contour(X,Y,T_pflotran[1,:,:],levels,colors='black',linewidth=0.5)
  plt.clabel(C,inline=True,fontsize=10)
  plt.xlabel('Distance (m)')
  plt.ylabel('Distance (m)')
  plt.title('Analytical (fill) vs. PFLOTRAN (contours) Solution, Temperature, t=0.02d')
  plt.subplot(223)
  # temperature values to contour against and compare visually
  levels = np.linspace(0.,1.,11)
  plt.contourf(X,Y,T_soln[2,:,:],levels,alpha=0.75)
  C = plt.contour(X,Y,T_pflotran[2,:,:],levels,colors='black',linewidth=0.5)
  plt.clabel(C,inline=True,fontsize=10)
  plt.xlabel('Distance (m)')
  plt.ylabel('Distance (m)')
  plt.title('Analytical (fill) vs. PFLOTRAN (contours) Solution, Temperature, t=0.04d')
  plt.subplot(224)
  # temperature values to contour against and compare visually
  levels = np.linspace(0.,1.,11)
  plt.contourf(X,Y,T_soln[3,:,:],levels,alpha=0.75)
  C = plt.contour(X,Y,T_pflotran[3,:,:],levels,colors='black',linewidth=0.5)
  plt.clabel(C,inline=True,fontsize=10)
  plt.xlabel('Distance (m)')
  plt.ylabel('Distance (m)')
  plt.title('Analytical (fill) vs. PFLOTRAN (contours) Solution, Temperature, t=0.10d')
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
  


