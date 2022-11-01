#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

soccer_diff_advect_weno3_err = np.array([ \
[ 16, 2.39541065112017890e+02,5.90267456054687500e+02,4.91427570581436157e-01,1.54292511940002441e+00,1.13716201782226562e+02, ], \
[ 32, 2.14247784775595761e+02,5.55230346679687500e+02,4.60284084081649780e-01,8.78463327884674072e-01,1.42950958251953125e+02, ], \
[ 64, 8.97210167120488506e+01,2.67552246093750000e+02,1.97596877813339233e-01,3.59150648117065430e-01,9.53643112182617188e+01, ], \
[ 128, 2.62928036435352404e+01,8.97547302246093750e+01,5.86575828492641449e-02,1.91223993897438049e-01,1.51416076660156250e+02, ], \
[ 256, 6.02156071300540496e+00,2.77148437500000000e+01,1.35223967954516411e-02,3.03839687258005142e-02,1.91376739501953125e+02, ], \
[ 512, 9.20525354921992500e-01,6.27313232421875000e+00,2.07407074049115181e-03,3.21510410867631435e-03,7.77862167358398438e+01, ], \
[ 1024, 1.07500602310178897e-01,1.00473022460937500e+00,2.42618567426688969e-04,3.48547531757503748e-04,2.44288940429687500e+01, ], \
[ 2048, 1.14801766200861820e-02,1.14349365234375000e-01,2.59313674177974463e-05,7.59834001655690372e-05,3.11767692565917969e+01, ], \
])

#  ADVECTION-DIFFUSION:
idx=3

xvec_space = np.pi/np.squeeze(soccer_diff_advect_weno3_err[idx:,0])
xvec_space_max = np.max(xvec_space)
xvec_space_min = np.min(xvec_space)
xvec_space_log = np.log10(xvec_space)
yvec_space = np.squeeze(soccer_diff_advect_weno3_err[idx:,4])
yvec_space_max = np.max(yvec_space)
yvec_space_min = np.min(yvec_space)
yvec_space_log = np.log10(yvec_space)
yvec=yvec_space_log
xvec=xvec_space_log
n=len(xvec)
coef_test_vec=(yvec[1:n]-yvec[0:n-1])/(xvec[1:n]-xvec[0:n-1])
mean_coef=np.mean(coef_test_vec)
print('coefs: ',coef_test_vec)
print('mean coef: ',mean_coef)


# PLOTTING

fsize=30
ms=500
lw=4.0

fig, ax = plt.subplots(num=None, figsize=[18,14], dpi=120, facecolor='w')
h_space = plt.scatter(xvec_space,yvec_space,s=ms,c='Blue',edgecolors='Blue',zorder=3,marker="o")
h_space2_fit = ax.plot([xvec_space_max,xvec_space_min], \
[10.0**(np.log10(yvec_space_min)+3.0*(np.log10(xvec_space_max)-np.log10(xvec_space_min))),\
yvec_space_min],'k--',linewidth=lw)
xlabel='$\Delta \\theta$, $\Delta \phi$'
ylabel='(CV)RMSD'
ax.set_ylabel(ylabel,fontsize=fsize)
ax.set_xlabel(xlabel,fontsize=fsize)
ax.tick_params(axis='y',labelsize=fsize)
ax.tick_params(axis='x',labelsize=fsize)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_aspect('equal')
ax.grid(zorder=0)
plt.xlim(left=0.001,right=0.05)
plt.ylim(bottom=0.00006,top=0.35)
plt.xticks(np.flipud(xvec_space),('$\\frac{\pi}{2048}$','$\\frac{\pi}{1024}$','$\\frac{\pi}{512}$','$\\frac{\pi}{256}$','$\\frac{\pi}{128}$'))
#plt.yticks((5e-2,1e-1,2e-1,4e-1),('$0.05$','$0.1$','$0.2$','$0.4$'))


plt.legend([h_space2_fit[0]],['$O(\Delta x^3)$'],\
            fontsize=fsize,loc='upper left')
fig.tight_layout()
fig.savefig("validation_advect_diff.png",bbox_inches='tight',dpi=120)
fig.savefig("validation_advect_diff.eps",bbox_inches='tight',dpi=120)
plt.close('all')














