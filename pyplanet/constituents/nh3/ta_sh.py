from ta_NH3setup import *

f = [1.25]
P = 1.
Tr = np.arange(200.,1000.,1.0)
T = 200.
fr = np.arange(10, 60, .1)
X_nh3=0.0015
X_he = 0.132002
X_h2 = 1.0 - (X_nh3+X_he)
X_partial=[X_h2,X_he,X_nh3]
print 100.0*np.array(X_partial)
P_dict = {'H2':0,'HE':1,'NH3':2}

akd = []
abell = []
asjsd = []
#for f in fr:
    # = nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar)
    #akd.append(a)
    #a = nh3_bell.alpha(f,T,P,X_partial,P_dict,otherPar)
    #abell.append(a)
    #a = nh3_sjsd.alpha(f,T,P,X_partial,P_dict,otherPar)
    #asjsd.append(a)
asjsd = nh3_sjsd.alpha(fr,T,P,X_partial,P_dict,otherPar)
plt.semilogy(fr,asjsd)
#plt.semilogy(fr,abell)


