#reads in file to look at limb

import numpy as np
import matplotlib.pyplot as plt

imke = np.loadtxt('jup_outputSolar_limbtest.dat.txt')

mu = np.delete(imke[0],0,0)
b = np.sqrt(1.0 - mu**2)
binp = [0.0]
for bbb in b:
    binp.append(bbb)

freq = list(np.delete(imke[:,0],0,0))


data = np.delete(np.delete(imke,0,0),0,1)


plt.figure(1)
for i in range(len(freq)):
    plt.plot(b,data[i])
plt.xlabel('b')

plt.figure(2)
for i in range(len(freq)):
    plt.plot(freq,data[:,i])
plt.xlabel('freq')
