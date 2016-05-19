from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

f = open('eos_water_density_test.txt')

line = f.readline()
column_headings = line.split()
line = f.readline()
w = line.split()
ntemp = int(w[0])
npres = int(w[1])
neos = int(w[2])

X = np.zeros((ntemp,npres),np.float)
Y = np.zeros((ntemp,npres),np.float)
Z = np.zeros((neos,ntemp,npres),np.float)
for itemp in range(ntemp):
  if itemp % 100 == 0:
    print(itemp)
  for ipres in range(npres):
    line = f.readline()
    w = line.split()
    X[itemp][ipres] = float(w[0])
    Y[itemp][ipres] = float(w[1])
    for ieos in range(neos):
      word = w[ieos+2]
      if word.startswith('NaN'):
        Z[ieos][itemp][ipres] = -999.
      else:
        Z[ieos][itemp][ipres] = float(word)

print('done reading file')

fig = plt.figure(figsize=(20,16))
fig.suptitle('Water Density vs. Temperature/Pressure vs. EOS',fontsize=24)
for ieos in range(4):
  ax = fig.add_subplot(2,2,ieos+1,projection='3d')
  ax.set_title(column_headings[ieos+2])
  ax.set_xlabel('Temperature [C]')
  ax.set_ylabel('Pressure [Pa]')
  ax.set_zlabel('Density [kg/m^3]')
  surf = ax.plot_surface(X, Y, Z[ieos][:][:], rstride=2, cstride=2, 
                         cmap=cm.coolwarm,
#                         shade=True, 
                         linewidth=0, antialiased=False)
  fig.colorbar(surf, shrink=0.5, aspect=5)
fig.subplots_adjust(hspace=0.02,wspace=0.02,
                    bottom=0.02,top = 0.94,
                    left=0.02,right=0.98)
plt.show()
  
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))


