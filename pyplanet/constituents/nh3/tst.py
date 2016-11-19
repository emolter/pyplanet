from ta_NH3setup import *

###Wrapper to test eliminating lines in NH3 formalism

f = [300.0]
il=2
P_1 = P_Jup[il]
T_1 = T_Jup[il]
X_partial=X_Jup[il]
print 'ATM LEVEL'
print P_1, T_1, X_partial
print

iv = np.arange(1300)
plt.figure()
alpha = []
for ii in iv:
    nh3_dbs.readInputFiles('./',im=0,rm=ii,vm=0,verbose=True)
    alpha.append(nh3_dbs.alpha(f,T_1,P_1,X_partial,P_dict,otherPar))
alpha = np.array(alpha)
plt.plot(iv,alpha/alpha[0])

