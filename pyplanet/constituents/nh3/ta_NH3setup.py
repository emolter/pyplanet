import nh3_kd
import nh3_bg
import nh3_sjs
import nh3_sjsd
import nh3_bell
import nh3_dbs
import nh3_dbs_sjs
import matplotlib.pyplot as plt
import math
import numpy as np

P_dict = {'H2':0,'HE':1,'NH3':2} # standard columns
otherPar = []                    # just need it defined

print 'Loading in Jupiter "solar":  ',
print '   use P_Jup, T_Jup, X_Jup'
jupdat = np.loadtxt('jupiter.paulSolar')
P_Jup = jupdat[:,2]
T_Jup = jupdat[:,1]
XH2_Jup  = jupdat[:,3]
XHe_Jup  = jupdat[:,4]
XNH3_Jup = jupdat[:,6]
X_Jup = []
for i in range(len(P_Jup)):
    X_Jup.append( [XH2_Jup[i],XHe_Jup[i],XNH3_Jup[i]] )






