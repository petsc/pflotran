import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math

filename = 'pflotran.out'
f = open(filename,'r')

time = []
dt = []
newton_iterations = []
linear_iterations = []
max_pressure_change = []
max_temperature_change = []
for line in f:
  if line.find('Dt=') > -1:
    dt.append(float(line[37:48]))
    time.append(float(line[20:32]))
  if line.find('newton =') > -1:
    newton_iterations.append(int(line[10:14]))
    linear_iterations.append(int(line[34:40]))
  if line.find('dpl=') > -1:
    max_pressure_change.append(float(line[20:33]))
  if line.find('dt=') > -1:
    max_temperature_change.append(float(line[38:51]))

f = plt.figure(figsize=(30,10))
f.suptitle("Convergence over Time",fontsize=16)

plt.subplot(2,3,1)
plt.xlabel('Time')
plt.ylabel('DT')
plt.ylim(0.,1.1*max(dt))
plt.plot(time,dt)

plt.subplot(2,3,2)
plt.xlabel('Time')
plt.ylabel('# Newton It')
plt.ylim(0.,1.1*max(newton_iterations))
plt.plot(time,newton_iterations)

plt.subplot(2,3,3)
plt.xlabel('Time')
plt.ylabel('# Linear It')
plt.ylim(0.,1.1*max(linear_iterations))
plt.plot(time,linear_iterations)

plt.subplot(2,3,4)
plt.xlabel('Time')
plt.ylabel('Maximum Pressure Change')
plt.ylim(0.,1.1*max(max_pressure_change))
plt.plot(time,max_pressure_change)

plt.subplot(2,3,5)
plt.xlabel('Time')
plt.ylabel('Maximum Temperature Change')
plt.ylim(0.,1.1*max(max_temperature_change))
plt.plot(time,max_temperature_change)

f.subplots_adjust(hspace=0.2,wspace=0.2,
                  bottom=.12,top=.9,
                  left=.12,right=.9)

plt.show()
