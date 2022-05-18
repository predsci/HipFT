#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

diff_only_scaling_space_err = np.array( \
[ 4096, 512, 5.95045776173929852e-03,1.43986368977948587e-02,1.34070324295705770e-05,2.29270079514263343e-05,2.18120942140924728e-01, ], \
)
diff_only_scaling_time_err = np.array( \
[ 32, 4096, 1.13954910979923333e-04,2.68808869350323221e-04,2.57508623242862414e-07,6.00821461087308830e-07,6.56498962920585072e-01, ], \
)
diff_only_scaling_space_data = np.array([ \
[ 32, 16, 6.18100059904120425e+00,1.08050866606477030e+01,1.25270957924716892e-02,2.31811837684234856e-02,2.11415238609474709e+00, ], \
[ 32, 32, 1.57911856270545758e+00,3.46126156065838586e+00,3.38148886266376484e-03,3.69360142028437536e-03,1.96006908973816507e-01, ], \
[ 32, 64, 3.94390868495534952e-01,9.05052434204890233e-01,8.67861196301150076e-04,1.09940212822741781e-03,1.57161455968773733e-01, ], \
[ 32, 128, 9.71711428748575912e-02,2.29287617961404067e-01,2.16737602438677498e-04,3.71879815788404204e-04,1.17031982785045396e-01, ], \
[ 32, 256, 2.40025544254899546e-02,5.76093319190249531e-02,5.38988651869625103e-05,9.72830514521579200e-05,4.39657864836592527e-01, ], \
[ 32, 512, 5.97045287024687305e-03,1.44416515670400258e-02,1.34520829647893923e-05,2.29670379514495549e-05,2.18120665295730876e-01, ], \
[ 32, 1024, 1.50158584188168432e-03,3.64090843925168883e-03,3.38892496995667301e-06,6.91456905866778980e-06,6.10528860755100489e-01, ], \
[ 32, 2048, 3.89707304293749203e-04,9.42176078297052300e-04,8.80268136042378642e-07,2.09113954269389770e-06,4.03575189499616849e-01, ], \
])
diff_only_scaling_time_data = np.array([ \
[ 1, 512, 2.81109326506552466e-02,6.05258268191164461e-02,6.33338170348833645e-05,7.07933042405072427e-05,2.17810423335038000e-01, ], \
[ 2, 512, 1.12905078309708320e-02,2.56586635532585206e-02,2.54384366848029414e-05,3.41062930344266139e-05,2.18046679592697440e-01, ], \
[ 4, 512, 7.25342263522701924e-03,1.71803500793430430e-02,1.63427052009892926e-05,2.55855157316813597e-05,2.18102788397424086e-01, ], \
[ 8, 512, 6.27248324994617978e-03,1.50899474374455167e-02,1.41325804177784432e-05,2.35752959597580001e-05,2.18116456466640490e-01, ], \
[ 16, 512, 6.03059458587616639e-03,1.45709592786715803e-02,1.35875869051296718e-05,2.30876332039732105e-05,2.18119827051921594e-01, ], \
[ 32, 512, 5.97045287024687305e-03,1.44416515670400258e-02,1.34520829647893923e-05,2.29670379514495549e-05,2.18120665295730876e-01, ], \
[ 64, 512, 5.95545167183604736e-03,1.44093810646381826e-02,1.34182840979570983e-05,2.29370023666564525e-05,2.18120872849552977e-01, ], \
[ 128, 512, 5.95170512960375726e-03,1.44013213099469795e-02,1.34098428465411530e-05,2.29295042223140619e-05,2.18120924680792250e-01, ], \
])

adv_only_scaling_space_err = np.array(\
[ 0.00234375, 512, 1.03802019544919631e+02,2.64777381800534386e+02,2.33364512064043461e-01,2.05924863779493927e+00,9.54020364255495515e+04, ], \
)
adv_only_scaling_time_err = np.array(\
[ 0.035, 4096, 9.36366694231118935e+00,2.37500889657619609e+01,1.92481648350241383e-02,1.10332867920165459e-01,1.09497810580730467e+05, ], \
)
adv_only_scaling_space_data = np.array([ \
[ 0.035, 16, 2.65832586509107728e+02,6.50308294962716673e+02,5.63718527916724388e-01,1.92706022758112239e+00,1.58759625198685743e+01, ], \
[ 0.035, 32, 2.68168852369625370e+02,6.86065449841636678e+02,6.06804970534637667e-01,1.63038639880939495e+00,2.18943421255221651e+01, ], \
[ 0.035, 64, 2.64379119822003020e+02,6.78843247643279710e+02,6.20045773367936426e-01,4.45464578596424854e+00,2.99131644863746442e+03, ], \
[ 0.035, 128, 2.31069074119812427e+02,5.90658009369955153e+02,5.49838284918585996e-01,3.17649366522541809e+00,2.56105433187001245e+03, ], \
[ 0.035, 256, 1.65664409458004201e+02,4.22764521770646866e+02,3.87518484316619893e-01,4.02706433048700241e+00,3.75092936258885966e+04, ], \
[ 0.035, 512, 9.97916489705713303e+01,2.54476075437499503e+02,2.23575098646921044e-01,1.66281369689299252e+00,1.20188972153800554e+05, ], \
[ 0.035, 1024, 5.28057110793683719e+01,1.34571836098023255e+02,1.13469706253124650e-01,6.01373844037326633e-01,1.21414445453928623e+05, ], \
[ 0.035, 2048, 2.47129241288685577e+01,6.29091121310374319e+01,5.16296097015795583e-02,2.67101528098655228e-01,3.31821854492790226e+04, ], \
])
adv_only_scaling_time_data = np.array([ \
[ 0.3, 512, 6.36869298114669178e+01,1.61490089790942335e+02,1.37942328815178983e-01,1.11573047520846735e+00,1.76767840479520179e+05, ], \
[ 0.15, 512, 8.49143170693997149e+01,2.16218162532605135e+02,1.87721846206592741e-01,9.82505585208350740e-01,9.48429310883712424e+03, ], \
[ 0.075, 512, 9.47517810438186388e+01,2.41524070156279095e+02,2.11345921408214760e-01,1.59732052002630542e+00,1.30435845177633601e+05, ], \
[ 0.0375, 512, 9.94807997100546828e+01,2.53677475218182167e+02,2.22818438447038358e-01,1.74253576275738054e+00,5.82902832128454684e+04, ], \
[ 0.01875, 512, 1.01798834755947595e+02,2.59632413126873416e+02,2.28468351856660118e-01,1.92100256918447654e+00,2.06479336630882579e+05, ], \
[ 0.009375, 512, 1.02946365302656133e+02,2.62579622058606333e+02,2.31271598360553166e-01,1.66310457292409009e+00,2.74311404027590252e+04, ], \
[ 0.0046875, 512, 1.03517275821900370e+02,2.64046069127667238e+02,2.32667781718723793e-01,1.81554036486517090e+00,8.56508555529466539e+04, ], \
[ 0.00234375, 512, 1.03802019544919631e+02,2.64777381800534386e+02,2.33364512064043461e-01,2.05924863779493927e+00,9.54020364255495515e+04, ], \
])

adv_diff_scaling_space_data_sc1 = np.array([ \
[ 1, 16, 2.39442988556413809e+02,5.90072610575660178e+02,5.56331862486284745e-01,1.90035274290316170e+00,1.56446237896795779e+01, ], \
[ 1, 32, 2.39147837775512443e+02,6.12427537321157615e+02,5.98702447935466897e-01,1.67451618881548070e+00,3.15414004600994673e+01, ], \
[ 1, 64, 2.14652579457495989e+02,5.48213381004065013e+02,5.56753846793051399e-01,5.66087492339979459e+00,1.86412620974176971e+04, ], \
[ 1, 128, 1.58418268824059595e+02,4.02551537186259850e+02,4.09470442448182315e-01,1.99558046634180553e+00,5.99816089263607296e+02, ], \
[ 1, 256, 9.89698561335192153e+01,2.51125227892461567e+02,2.46952466327332532e-01,1.75441721417343310e+00,2.16888038347284128e+04, ], \
[ 1, 512, 5.57168335196979143e+01,1.41268601090014982e+02,1.33781031203049405e-01,7.40790799604525763e-01,2.25011852471673665e+04, ], \
[ 1, 1024, 2.96144039286726155e+01,7.50600731312115954e+01,6.92126588772629392e-02,4.91202749227063862e-01,1.87036528445763746e+05, ], \
[ 1, 2048, 1.52733437126358567e+01,3.87057466548935167e+01,3.51363299861303935e-02,1.82550640371447415e-01,1.97396186278373279e+04, ], \
])

adv_diff_scaling_space_data_sc30 = np.array([ \
[ 30, 16, 2.39442988350971035e+02,5.90072605707771913e+02,5.56331873483073669e-01,1.90035276600588432e+00,1.56446240214351970e+01, ], \
[ 30, 32, 2.39147837770872059e+02,6.12427535693798291e+02,5.98702451830708360e-01,1.67451619672073759e+00,3.15414006085000409e+01, ], \
[ 30, 64, 2.14652579482171461e+02,5.48213380806650207e+02,5.56753847411173730e-01,5.66087435608420719e+00,1.86412573547277243e+04, ], \
[ 30, 128, 1.58418268838626091e+02,4.02551537173928466e+02,4.09470442559987213e-01,1.99558046256688448e+00,5.99816088497475221e+02, ], \
[ 30, 256, 9.89698561400292078e+01,2.51125227896538604e+02,2.46952466355034927e-01,1.75441721635081072e+00,2.16888037624182361e+04, ], \
[ 30, 512, 5.57168335235948717e+01,1.41268601094202211e+02,1.33781031215191221e-01,7.40790846992025132e-01,2.25011856222101087e+04, ], \
[ 30, 1024, 2.96144039338841623e+01,7.50600731385135873e+01,6.92126588911213120e-02,4.91202764152502525e-01,1.87036550863870681e+05, ], \
])


#  DIFFUSION-ONLY:

xvec_space = np.pi/np.squeeze(diff_only_scaling_space_data[:,1])
xvec_space_max = np.max(xvec_space)
xvec_space_min = np.min(xvec_space)
xvec_space_log = np.log10(xvec_space)

yvec_space = np.squeeze(diff_only_scaling_space_data[:,4])
yvec_space_max = np.max(yvec_space)
yvec_space_min = np.min(yvec_space)
yvec_space_log = np.log10(yvec_space)

yvec_space2 = np.squeeze(diff_only_scaling_space_data[:,4]) - diff_only_scaling_time_err[4]
yvec_space2_max = np.max(yvec_space2)
yvec_space2_min = np.min(yvec_space2)
yvec_space2_log = np.log10(yvec_space2)


yvec=yvec_space_log
xvec=xvec_space_log
n=len(xvec)
coef_test_vec=(yvec[1:n]-yvec[0:n-1])/(xvec[1:n]-xvec[0:n-1])
mean_coef=np.mean(coef_test_vec)
print('SPACE coefs: ',coef_test_vec)
print('SPACE mean coef: ',mean_coef)

yvec=yvec_space2_log
xvec=xvec_space_log
n=len(xvec)
coef_test_vec=(yvec[1:n]-yvec[0:n-1])/(xvec[1:n]-xvec[0:n-1])
mean_coef=np.mean(coef_test_vec)
print('SPACE-CORR coefs: ',coef_test_vec)
print('SPACE-CORR mean coef: ',mean_coef)


xvec_time = 672.0/np.squeeze(diff_only_scaling_time_data[:,0])
xvec_time_max = np.max(xvec_time)
xvec_time_min = np.min(xvec_time)
xvec_time_log = np.log10(xvec_time)

yvec_time = np.squeeze(diff_only_scaling_time_data[:,4])
yvec_time_max = np.max(yvec_time)
yvec_time_min = np.min(yvec_time)
yvec_time_log = np.log10(yvec_time)

yvec_time2 = np.squeeze(diff_only_scaling_time_data[:,4]) \
                      - diff_only_scaling_space_err[4]
yvec_time2_max = np.max(yvec_time2)
yvec_time2_min = np.min(yvec_time2)
yvec_time2_log = np.log10(yvec_time2)

yvec=yvec_time_log
xvec=xvec_time_log
n=len(xvec)
coef_test_vec=(yvec[1:n]-yvec[0:n-1])/(xvec[1:n]-xvec[0:n-1])
mean_coef=np.mean(coef_test_vec)
print('TIME coefs: ',coef_test_vec)
print('TIME mean coef: ',mean_coef)

yvec=yvec_time2_log
xvec=xvec_time_log
n=len(xvec)
coef_test_vec=(yvec[1:n]-yvec[0:n-1])/(xvec[1:n]-xvec[0:n-1])
mean_coef=np.mean(coef_test_vec)
print('TIME-CORR coefs: ',coef_test_vec)
print('TIME-CORR mean coef: ',mean_coef)


#PLOTTING

fsize=30
ms=500
lw=4.0




fig, ax = plt.subplots(num=None, figsize=[18,14], dpi=120, facecolor='w')
h_space = plt.scatter(xvec_space,yvec_space,s=ms,c='cornflowerblue',edgecolors='cornflowerblue',zorder=3,marker="o")
h_space_corr = plt.scatter(xvec_space,yvec_space2,s=ms-0.1*ms,c='b',edgecolors='b',zorder=3,marker="o")
h_space_fit = ax.plot([xvec_space_max,xvec_space_min], \
[10.0**(np.log10(yvec_space2_min)+2.0*(np.log10(xvec_space_max)-np.log10(xvec_space_min))),\
yvec_space2_min],'k--',linewidth=lw,zorder=5)
plt.scatter(xvec_space[5],yvec_space2[5],s=ms,c='Yellow',edgecolors='b',zorder=4,marker="o",linewidth=lw)
xlabel='$\Delta \\theta$, $\Delta \phi$ '
ylabel='(CV)RMSD'
ax.set_ylabel(ylabel,fontsize=fsize)
ax.set_xlabel(xlabel,fontsize=fsize)
ax.tick_params(axis='y',labelsize=fsize)
ax.tick_params(axis='x',labelsize=fsize)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_aspect('equal')
ax.grid(zorder=0)
plt.xlim(left=7e-4,right=0.3)
plt.ylim(bottom=3e-7,top=3e-1)
plt.xticks((np.pi/2048,np.pi/512,np.pi/128,np.pi/32),('$\\frac{\pi}{2048}$','$\\frac{\pi}{512}$','$\\frac{\pi}{128}$','$\\frac{\pi}{32}$'))
plt.legend([h_space,h_space_corr,h_space_fit[0]],[\
                                   'E($\Delta x$,$\Delta t_{0}$)', \
                                   'E($\Delta x$,$\Delta t_{0}$)-$\\widetilde{E(\Delta t_{0})}$', \
                                   '$O(\Delta x^2)$', \
                                   ],fontsize=fsize,loc='upper left')
fig.tight_layout()
fig.savefig("validation_diff_space.png",bbox_inches='tight',dpi=120)
fig.savefig("validation_diff_space.eps",bbox_inches='tight',dpi=120)
plt.close('all')



fig, ax = plt.subplots(num=None, figsize=[18,14], dpi=120, facecolor='w')
h_time = plt.scatter(xvec_time,yvec_time,s=ms,c='Pink',edgecolors='Pink',zorder=3,marker="d")
h_time_corr = plt.scatter(xvec_time,yvec_time2,s=ms,c='r',edgecolors='r',zorder=3,marker="d")
h_time_fit = ax.plot([xvec_time_max,xvec_time_min], \
[10.0**(np.log10(yvec_time2_min)+2.0*(np.log10(xvec_time_max)-np.log10(xvec_time_min))),\
yvec_time2_min],'k--',linewidth=lw,zorder=5)
#plt.scatter(xvec_time[4],yvec_time2[4],s=ms,c='Yellow',edgecolors='r',zorder=4,marker="d",linewidth=lw)
xlabel='$\Delta t$ (hours)'
ylabel='(CV)RMSD'
ax.set_ylabel(ylabel,fontsize=fsize)
ax.set_xlabel(xlabel,fontsize=fsize)
ax.tick_params(axis='y',labelsize=fsize)
ax.tick_params(axis='x',labelsize=fsize)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_aspect('equal')
ax.grid(zorder=0)
plt.xlim(left=2.5,right=1000)
plt.ylim(bottom=1e-9,top=2e-3)
plt.legend([h_time,h_time_corr,h_time_fit[0]],[\
                                 'E($\Delta x_{0}$,$\Delta t$)', \
                                 'E($\Delta x_{0}$,$\Delta t$)-$\\widetilde{E(\Delta x_{0})}$', \
                                 '$O(\Delta t^2)$', \
                                ],fontsize=fsize)
fig.tight_layout()
fig.savefig("validation_diff_time.png",bbox_inches='tight',dpi=120)
fig.savefig("validation_diff_time.eps",bbox_inches='tight',dpi=120)
plt.close('all')





#  ADVECTION-ONLY:
idx=4

xvec_space = np.pi/np.squeeze(adv_only_scaling_space_data[idx:,1])
xvec_space_max = np.max(xvec_space)
xvec_space_min = np.min(xvec_space)
xvec_space_log = np.log10(xvec_space)

yvec_space = np.squeeze(adv_only_scaling_space_data[idx:,4])
yvec_space_max = np.max(yvec_space)
yvec_space_min = np.min(yvec_space)
yvec_space_log = np.log10(yvec_space)

yvec_space2 = np.squeeze(adv_only_scaling_space_data[idx:,4]) - adv_only_scaling_time_err[4]
yvec_space2_max = np.max(yvec_space2)
yvec_space2_min = np.min(yvec_space2)
yvec_space2_log = np.log10(yvec_space2)


yvec=yvec_space_log
xvec=xvec_space_log
n=len(xvec)
coef_test_vec=(yvec[1:n]-yvec[0:n-1])/(xvec[1:n]-xvec[0:n-1])
mean_coef=np.mean(coef_test_vec)
print('SPACE coefs: ',coef_test_vec)
print('SPACE mean coef: ',mean_coef)

yvec=yvec_space2_log
xvec=xvec_space_log
n=len(xvec)
coef_test_vec=(yvec[1:n]-yvec[0:n-1])/(xvec[1:n]-xvec[0:n-1])
mean_coef=np.mean(coef_test_vec)
print('SPACE-CORR coefs: ',coef_test_vec)
print('SPACE-CORR mean coef: ',mean_coef)



xvec_time = np.squeeze(adv_only_scaling_time_data[:,0])
xvec_time_max = np.max(xvec_time)
xvec_time_min = np.min(xvec_time)
xvec_time_log = np.log10(xvec_time)

yvec_time = np.squeeze(adv_only_scaling_time_data[:,4])
yvec_time_max = np.max(yvec_time)
yvec_time_min = np.min(yvec_time)
yvec_time_log = np.log10(yvec_time)

yvec_time2 = np.squeeze(adv_only_scaling_time_data[:,4]) - adv_only_scaling_space_err[4]
yvec_time2_max = np.max(yvec_time2)
yvec_time2_min = np.min(yvec_time2)
yvec_time2_log = np.log10(yvec_time2)

yvec=yvec_time_log
xvec=xvec_time_log
n=len(xvec)
coef_test_vec=(yvec[1:n]-yvec[0:n-1])/(xvec[1:n]-xvec[0:n-1])
mean_coef=np.mean(coef_test_vec)
print('TIME coefs: ',coef_test_vec)
print('TIME mean coef: ',mean_coef)

yvec=yvec_time2_log
xvec=xvec_time_log
n=len(xvec)
coef_test_vec=(yvec[1:n]-yvec[0:n-1])/(xvec[1:n]-xvec[0:n-1])
mean_coef=np.mean(coef_test_vec)
print('TIME-CORR coefs: ',coef_test_vec)
print('TIME-CORR mean coef: ',mean_coef)

fsize=40
ms=800
lw=6
fig, ax = plt.subplots(num=None, figsize=[18,14], dpi=120, facecolor='w')
h_space = plt.scatter(xvec_space,yvec_space,s=ms,c='b',edgecolors='b',zorder=3,marker="o")
h_space_corr = plt.scatter(xvec_space,yvec_space2,s=ms,c='c',edgecolors='c',zorder=3,marker="o")
h_space_fit = ax.plot([xvec_space_max,xvec_space_min], \
[10.0**(np.log10(yvec_space_min)+1.0*(np.log10(xvec_space_max)-np.log10(xvec_space_min))),\
yvec_space_min],'k--',linewidth=lw)
plt.scatter(xvec_space[1],yvec_space[1],s=ms,c='Yellow',edgecolors='b',zorder=4,marker="o",linewidth=lw)
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
plt.xlim(left=0.00125,right=0.014)
plt.ylim(bottom=0.035,top=0.5)
plt.xticks(np.flipud(xvec_space),('$\\frac{\pi}{2048}$','$\\frac{\pi}{1024}$','$\\frac{\pi}{512}$','$\\frac{\pi}{256}$'))
plt.yticks((5e-2,1e-1,2e-1,4e-1),('$0.05$','$0.1$','$0.2$','$0.4$'))
plt.legend([h_space,h_space_fit[0]],[\
                                   'E($\Delta x$,$\Delta t_{0}$)', \
                                   '$O(\Delta x)$', \
                                    ],fontsize=fsize,loc='upper left')
fig.tight_layout()
fig.savefig("validation_advect_space.png",bbox_inches='tight',dpi=120)
fig.savefig("validation_advect_space.eps",bbox_inches='tight',dpi=120)
plt.close('all')



fig, ax = plt.subplots(num=None, figsize=[18,14], dpi=120, facecolor='w')
h_time = plt.scatter(xvec_time,yvec_time,s=ms,c='Pink',edgecolors='Pink',zorder=3,marker="d")
h_time_corr = plt.scatter(xvec_time,yvec_time2,s=ms,c='r',edgecolors='r',zorder=3,marker="d")
h_time_fit = ax.plot([xvec_time_max,xvec_time_min], \
[10.0**(np.log10(yvec_time2_min)+1.0*(np.log10(xvec_time_max)-np.log10(xvec_time_min))),\
yvec_time2_min],'k--',linewidth=lw)
xlabel='$\Delta t$'
ylabel='(CV)RMSD'
ax.set_ylabel(ylabel,fontsize=fsize)
ax.set_xlabel(xlabel,fontsize=fsize)
ax.tick_params(axis='y',labelsize=fsize)
ax.tick_params(axis='x',labelsize=fsize)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_aspect('equal')
ax.grid(zorder=0)
#plt.xlim(left=0.00234375,right=0.3)
#plt.ylim(bottom=1e-9,top=2e-3)
plt.legend([h_time,h_time_corr,h_time_fit[0]],[\
                                 'E($\Delta x_{0}$,$\Delta t$)', \
                                 'E($\Delta x_{0}$,$\Delta t$)-$\\widetilde{E(\Delta x_{0})}$', \
                                 '$O(\Delta t^2)$', \
                                ],fontsize=fsize)
fig.tight_layout()
fig.savefig("validation_advect_time.png",bbox_inches='tight',dpi=120)
fig.savefig("validation_advect_time.eps",bbox_inches='tight',dpi=120)
plt.close('all')




#  ADVECTION-DIFFUSION:
idx=4

xvec_space = np.pi/np.squeeze(adv_diff_scaling_space_data_sc1[idx:,1])
xvec_space_max = np.max(xvec_space)
xvec_space_min = np.min(xvec_space)
xvec_space_log = np.log10(xvec_space)
yvec_space = np.squeeze(adv_diff_scaling_space_data_sc1[idx:,4])
yvec_space_max = np.max(yvec_space)
yvec_space_min = np.min(yvec_space)
yvec_space_log = np.log10(yvec_space)
yvec=yvec_space_log
xvec=xvec_space_log
n=len(xvec)
coef_test_vec=(yvec[1:n]-yvec[0:n-1])/(xvec[1:n]-xvec[0:n-1])
mean_coef=np.mean(coef_test_vec)
print('SPACE_SC1 coefs: ',coef_test_vec)
print('SPACE_SC1 mean coef: ',mean_coef)
xvec_space2 = np.pi/np.squeeze(adv_diff_scaling_space_data_sc30[:,1])
xvec_space2_max = np.max(xvec_space)
xvec_space2_min = np.min(xvec_space)
xvec_space2_log = np.log10(xvec_space)
yvec_space2 = np.squeeze(adv_diff_scaling_space_data_sc30[:,4])
yvec_space2_max = np.max(yvec_space)
yvec_space2_min = np.min(yvec_space)
yvec_space2_log = np.log10(yvec_space)
yvec=yvec_space2_log
xvec=xvec_space2_log
n=len(xvec)
coef_test_vec=(yvec[1:n]-yvec[0:n-1])/(xvec[1:n]-xvec[0:n-1])
mean_coef=np.mean(coef_test_vec)
print('SPACE_SC30 coefs: ',coef_test_vec)
print('SPACE_SC30 mean coef: ',mean_coef)

# PLOTTING

fig, ax = plt.subplots(num=None, figsize=[18,14], dpi=120, facecolor='w')
h_space = plt.scatter(xvec_space,yvec_space,s=ms,c='Purple',edgecolors='Purple',zorder=3,marker="o")
#h_space_fit = ax.plot([xvec_space_min,xvec_space_max], \
#[10.0**(np.log10(yvec_space_max)-1.0*(np.log10(xvec_space_max)-np.log10(xvec_space_min))),yvec_space_max],'k--',linewidth=2.0,markersize=1.0)
#h_space2 = plt.scatter(xvec_space2,yvec_space2,s=300,c='b',edgecolors='b',zorder=3,marker="o")
h_space2_fit = ax.plot([xvec_space_max,xvec_space_min], \
[10.0**(np.log10(yvec_space_min)+1.0*(np.log10(xvec_space_max)-np.log10(xvec_space_min))),\
yvec_space_min],'k--',linewidth=lw)
plt.scatter(xvec_space[1],yvec_space[1],s=ms,c='Yellow',edgecolors='Purple',zorder=4,marker="o",linewidth=lw)
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

plt.xlim(left=0.00125,right=0.014)
plt.ylim(bottom=0.03,top=0.4)

plt.xticks(np.flipud(xvec_space),('$\\frac{\pi}{2048}$','$\\frac{\pi}{1024}$','$\\frac{\pi}{512}$','$\\frac{\pi}{256}$'))
plt.yticks((5e-2,1e-1,2e-1,4e-1),('$0.05$','$0.1$','$0.2$','$0.4$'))


plt.legend([h_space,h_space2_fit[0]],['E($\Delta x$,$\Delta t_{stb}$)', \
                     '$O(\Delta x)$', \
                     ],fontsize=fsize,loc='upper left')
fig.tight_layout()
fig.savefig("validation_advect_diff.png",bbox_inches='tight',dpi=120)
fig.savefig("validation_advect_diff.eps",bbox_inches='tight',dpi=120)
plt.close('all')














