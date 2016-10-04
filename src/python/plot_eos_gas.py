from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

f = open('eos_gas_test.txt')

line = f.readline()
column_headings = line.split(',')
line = f.readline()
w = line.split()
ntemp = int(w[0])
npres = int(w[1])
nvar = 6

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

fig = plt.figure(figsize=(16,13))
fig.suptitle('Gas EOS Constitutive Relations vs. Temperature/Pressure',fontsize=24)
for ivar in range(nvar-2):
  ax = fig.add_subplot(3,2,ivar+1,projection='3d')
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

for ivar in range(nvar-2,nvar):
  ax = fig.add_subplot(3,2,ivar+1)
  ax.set_title(column_headings[ivar+2].split(' [')[0])
  ax.set_xlabel('Temperature [C]')
  w0 = column_headings[ivar+2].split('(')[0]
  w1 = column_headings[ivar+2].split(')')[1]
  heading = w0+w1
  ax.set_ylabel(heading)
# for some reason, cannot index numpy arrays passed to plot
#surf = ax.plot(X[:][0],Z[ivar][:][0])
# therefore, have to repack.
  XX = np.zeros((ntemp),np.float)
  YY = np.zeros((ntemp),np.float)
  for i in range(ntemp):
    XX[i] = X[i][0]
    YY[i] = Z[ivar][i][0]
  surf = ax.plot(XX,YY)

fig.subplots_adjust(hspace=0.17,wspace=0.15,
                    bottom=0.05,top = 0.93,
                    left=0.05,right=0.98)
plt.show()
  
