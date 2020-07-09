from scipy import interpolate
import numpy as np

x = np.linspace(0,10)
y = np.linspace(0,10)

xx, yy = np.meshgrid(x,y)

z = 0*xx + 0*yy

for i in range(0,len(z[:][0])):
    z[i][0] = 1

xnew = np.linspace(1,10,100)
ynew = np.linspace(1,10,100)

f = interpolate.interp2d(x,y,z)

znew = f(xnew,ynew)
print(z)
print(znew)
