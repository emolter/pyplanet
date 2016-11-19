import img
import numpy as np
import matplotlib.pyplot as plt

jp = np.loadtxt('jupprfull.dat')
b = jp[:,0]
TB = jp[:,1]
kern = np.exp(-b**2/0.0002)

coslat = []
for bb in b:
    if abs(bb)>1.0:
        cl = 0.0
    else:
        cl = np.sqrt(1.0 - bb**2)
    coslat.append(cl)
    
lat = np.sign(b)*np.arccos(coslat)*180.0/np.pi
emp2 = max(TB)*np.power(coslat,0.2)
emp08= max(TB)*np.power(coslat,0.08)

x = lat
plt.plot(x,TB,label='hi')
plt.plot(x,kern*280.,label='kernel')
plt.plot(x,emp08,label='0.08')
plt.plot(x,emp2,label='0.20')

cvl = img.convolvend(TB,kern,normalize_kernel=True)

plt.plot(x,cvl,label='cnvld')
plt.legend(loc='upper left')
