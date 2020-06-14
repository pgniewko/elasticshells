#! /usr/bin/env python


import matplotlib.pyplot as plt
import numpy as np


fv = '/Users/pawel/Desktop/single-shell/V/1_shells_seed_1111_1.0.out'
fp = '/Users/pawel/Desktop/single-shell/P/1_shells_seed_1111_1.0.out'
fn = '/Users/pawel/Desktop/single-shell/N/1_shells_seed_1111_1.0.out'


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

dv = np.loadtxt(fv)
dp = np.loadtxt(fp)
dn = np.loadtxt(fn)

ax.plot(dv[:,5]/dv[:,14],dv[:,11]/dv[:,2],label='V-const',lw=3,color='green')
ax.plot(dp[:,5]/dp[:,14],dp[:,11]/dp[:,2],'--',label='P-const',lw=2)
ax.plot(dn[:,5]/dn[:,14],dn[:,11]/dn[:,2],'-.',label='N-const',lw=2)

ax.set_ylabel('$\phi_P$',fontsize=15)
ax.set_xlabel('$P_P/<\Pi>$',fontsize=15)
ax.set_xscale('log')

plt.legend()
plt.show()
