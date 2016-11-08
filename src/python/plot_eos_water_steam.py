from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

f = open('eos_water_steam_test.txt')

line = f.readline()
column_headings = line.split(',')
line = f.readline()
w = line.split()
ntemp = int(w[0])
npres = int(w[1])
nvar = 2

X = np.zeros((ntemp,npres),np.float)
Y = np.zeros((ntemp,npres),np.float)
Z = np.zeros((nvar,ntemp,npres),np.float)
for itemp in range(ntemp):
  if ntemp*npres > 99999:
    if itemp % ntemp/10 == 0:
      print('%.0f%% done' % (float(itemp)/(ntemp/1000.)))
  for ipres in range(npres):
    line = f.readline()
    w = line.split()
    X[itemp][ipres] = float(w[0])
    Y[itemp][ipres] = float(w[1])
    for ivar in range(nvar):
      word = w[ivar+2]
      if word.startswith('NaN'):
        Z[ivar][itemp][ipres] = 999.e-10
      else:
        Z[ivar][itemp][ipres] = float(word)
f.close()

f = open('eos_water_steam_psat_test.txt')

line = f.readline()
column_headings_psat = line.split(',')
line = f.readline()
w = line.split()
ntemp_psat = int(w[0])
nvar_psat = 3

X_psat = np.zeros((ntemp_psat),np.float)
Z_psat = np.zeros((nvar_psat,ntemp_psat),np.float)
for itemp in range(ntemp_psat):
  if ntemp_psat > 99999:
    if itemp % ntemp_psat/10 == 0:
      print('%.0f%% done' % (float(itemp)/(ntemp_psat/1000.)))
  line = f.readline()
  w = line.split()
  X_psat[itemp] = float(w[0])
  for ivar in range(nvar_psat):
    word = w[ivar+1]
    if word.startswith('NaN'):
      Z_psat[ivar][itemp] = 999.e-10
    else:
      Z_psat[ivar][itemp] = float(word)
f.close()

fig = plt.figure(figsize=(20,16))
fig.suptitle('Water EOS Steam Constitutive Relations vs. Temperature/Pressure',fontsize=24)
for ivar in range(nvar):
  ax = fig.add_subplot(2,3,ivar+1,projection='3d')
  ax.set_title(column_headings[ivar+2].split(' [')[0])
  ax.set_xlabel('Temperature [C]')
  ax.set_ylabel('Pressure [Pa]')
  w0 = column_headings[ivar+2].split('(')[0]
  w1 = column_headings[ivar+2].split(')')[1]
  heading = w0+w1
  ax.set_zlabel(heading)
  surf = ax.plot_surface(X, Y, Z[ivar][:][:], rstride=2, cstride=2, 
                         cmap=cm.coolwarm,
                         linewidth=0, antialiased=False)
  fig.colorbar(surf, shrink=0.5, aspect=5)

# Saturaton pressure
offset = nvar+2
for ivar in range(nvar_psat):
  ax = fig.add_subplot(2,3,ivar+offset)
  ax.set_title(column_headings_psat[ivar+1].split(' [')[0])
  ax.set_xlabel('Temperature [C]')
  w0 = column_headings_psat[ivar+1].split('(')[0]
  w1 = column_headings_psat[ivar+1].split(')')[1]
  heading = w0+w1
  ax.set_ylabel(heading)
  # for some reason, cannot index numpy arrays passed to plot
  #surf = ax.plot(X[:][0],Z[ivar][:][0])
  # therefore, have to repack.
  XX = np.zeros((ntemp_psat),np.float)
  YY = np.zeros((ntemp_psat),np.float)
  for i in range(ntemp_psat):
    XX[i] = X_psat[i]
    YY[i] = Z_psat[ivar][i]
  surf = ax.plot(XX,YY)

fig.subplots_adjust(hspace=0.12,wspace=0.12,
                    bottom=0.05,top = 0.93,
                    left=0.02,right=0.98)
plt.show()
  
