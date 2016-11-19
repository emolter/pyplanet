from ta_NH3setup import *

f         = [   15.,    3.,   2.5,   1.5,    75.,   75.,    4.,   0.6] #GHz
P         = [    1.,   10.,   50.,   80.,    22.,   12.,   12.,  250.] #bars
T         = [  160.,  325.,  550.,  600.,   330.,  330.,  330.,  900.] #K
XNH3      = [  2E-4,  2E-4,  2E-4,  2E-4,   2E-4,  2E-4,  2E-4,  2E-4]
XH2       = [  0.87,  0.87,  0.87,  0.87,   0.87,  0.87,  0.87,  0.87]
XHE       = [  0.13,  0.13,  0.13,  0.13,   0.13,  0.13,  0.13,  0.13]
aBellotti = [0.6873,0.1150,0.1922,0.0802,16.2537,5.1798,0.2562,0.0111] #dB/km
aDevaraj  = [0.7399,0.1136,0.1902,0.0813,16.7760,3.7556,0.2645,0.0125] #dB/km

akd = []
abell = []
adbs = []
db = []
dd = []
for i,ff in enumerate(f):
    X_partial = [XH2[i],XHE[i],XNH3[i]]
    a = nh3_kd.alpha(ff,T[i],P[i],X_partial,P_dict,otherPar)
    akd.append(a[0])
    dd.append(abs(aDevaraj[i] - a[0]))
    a = nh3_bell.alpha(ff,T[i],P[i],X_partial,P_dict,otherPar)
    abell.append(a[0])
    a = nh3_dbs.alpha(ff,T[i],P[i],X_partial,P_dict,otherPar)
    adbs.append(a[0])
    db.append(abs(aBellotti[i] - a[0]))
print 'VALUES'
print adbs
print akd

dd = np.array(dd)
db = np.array(db)

o = T
plt.figure(1)
plt.semilogy(o,akd,'o',label='devaraj')
plt.semilogy(o,aDevaraj,'x',label='devaraj_from_amadeo')
plt.semilogy(o,adbs,'o',label='new_bellotti')
plt.semilogy(o,aBellotti,'x',label='bellotti_from_amadeo')
plt.semilogy(o,abell,'o',label='old_bellotti')

pt = 'diff'
if pt == 'frac':
    ab = np.array(aBellotti)
    ad = np.array(aDevaraj)
else:
    ab = np.ones(np.shape(aBellotti))
    ad = np.ones(np.shape(aDevaraj))
plt.figure(2)
plt.loglog(P,db/ab,'gd',label='P_dbell')
plt.loglog(P,dd/ad,'bd',label='P_ddev')
plt.loglog(T,db/ab,'gs',label='T_dbell')
plt.loglog(T,dd/ad,'bs',label='T_ddev')
plt.loglog(f,db/ab,'go',label='f_dbell')
plt.loglog(f,dd/ad,'bo',label='f_ddev')
print 'DIFFERENTIALS ',ab
print db
print dd

