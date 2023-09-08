#!/usr/bin/env python3
from scipy.io.idl import readsav
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
import os

width = 0.6
fsize = 22
mksz = 15
lsz = 3.0
fc = 'w'
tc = 'k'
dpi=120  
tfac = 1/24.0
flux_fac = 1e-21

time_start = Time('2011-01-01T00:00:00.00',scale='utc')
time_end = Time('2022-01-01T00:00:00.00',scale='utc')

# Get years jd and string array:
x_labels = []
x_ticks = np.asarray(range(2011,2023))
ii=0
for i in range(2011,2023):
  date_str = str(i)+'-01-01T00:00:00.00'
  curr_time = Time(date_str,scale='utc')
  x_ticks[ii] = curr_time.jd
  x_labels = x_labels + [str(i)]
  ii = ii+1

#
# READ IN AFT DATA:
#
data_aft = readsav('AFT_Analysis_Data_v2.sav',python_dict=True)
jd_aft = np.asarray(data_aft['jddate'])
flux_tot_unsigned_aft = np.asarray(data_aft['tuf'])
flux_tot_signed_aft = np.asarray(data_aft['tsf'])
eq_dipole_aft = np.asarray(data_aft['edipole'])
axial_dipole_aft = np.asarray(data_aft['adipole'])
pole_n_aft = np.asarray(data_aft['northpole60'])
pole_s_aft = np.asarray(data_aft['southpole60'])

print(np.min(flux_tot_signed_aft))
print(np.max(flux_tot_signed_aft))


#
# READ IN OFT DATA:
#
data_oft = pd.read_table('/home/sumseq/Dropbox/PSI/SWQU/OFT/codes/hipft/tests/conflow_dev/test_11yrs_nodiff_weno3_conflow1cr_v4/history_sol.dat',header=0,sep='\s+')
hr_oft  = np.asarray(data_oft['TIME'])
init_frame_time = Time('2010-12-31T23:59:52.00',scale='utc')
jd_oft = init_frame_time.jd + hr_oft/24.0
fluxp = np.array(data_oft['FLUX_POSITIVE'])
fluxm = np.array(data_oft['FLUX_NEGATIVE'])
fluxp_pn = np.array(data_oft['NPOLE_FLUX_POSITIVE'])
fluxm_pn = np.array(data_oft['NPOLE_FLUX_NEGATIVE'])
area_pn = np.array(data_oft['NPOLE_AREA'])
fluxp_ps = np.array(data_oft['SPOLE_FLUX_POSITIVE'])
fluxm_ps = np.array(data_oft['SPOLE_FLUX_NEGATIVE'])
area_ps = np.array(data_oft['SPOLE_AREA'])
axial_dipole_oft = np.array(data_oft['AX_DIPOLE'])
eq_dipole_oft = np.array(data_oft['EQ_DIPOLE'])
flux_tot_unsigned_oft = np.abs(fluxp) + np.abs(fluxm)
flux_tot_signed_oft = fluxp + fluxm
pole_n_oft = (fluxp_pn+fluxm_pn)/area_pn
pole_s_oft = (fluxp_ps+fluxm_ps)/area_ps


print(np.min(flux_tot_signed_oft))
print(np.max(flux_tot_signed_oft))

#
# ############## PLOT ##############
#


#
# ****** Total unsigned flux.
#
fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, \
                 facecolor=fc, frameon=True)
ax = plt.gca()
h = plt.plot(jd_aft, flux_fac*flux_tot_unsigned_aft, \
             'k-',linewidth=lsz, markersize=mksz)
hoft = plt.plot(jd_oft, flux_fac*flux_tot_unsigned_oft, \
             'b-',linewidth=0.8*lsz, markersize=mksz)
plt.xlim(xmin=time_start.jd,xmax=time_end.jd)
plt.title('Total Unsigned Flux', {'fontsize': fsize, 'color': tc})
plt.xlabel('Year', {'fontsize': fsize, 'color': tc})
plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
ax.tick_params(axis='y',labelsize=fsize)
ax.tick_params(axis='x',labelsize=fsize)
plt.xticks(x_ticks,x_labels)
ax.grid(zorder=0)
plt.legend([h[0],hoft[0]],["AFT","OFT"],loc='upper right',fontsize=fsize)
fig.tight_layout()
fig.savefig('oft_vs_aft_flux_total_unsigned.png', bbox_inches="tight",   \
            pad_inches=0.1, dpi=dpi, facecolor=fig.get_facecolor(), \
            edgecolor=None)
plt.close('all')


#
# ****** Total signed flux.
#
fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, \
                 facecolor=fc, frameon=True)
ax = plt.gca()
h = plt.plot(jd_aft, flux_fac*flux_tot_signed_aft, \
             'k-',linewidth=lsz, markersize=mksz)
hoft = plt.plot(jd_oft, flux_fac*flux_tot_signed_oft, \
             'b-',linewidth=0.8*lsz, markersize=mksz)
plt.xlim(xmin=time_start.jd,xmax=time_end.jd)
plt.title('Total Signed Flux', {'fontsize': fsize, 'color': tc})
plt.xlabel('Year', {'fontsize': fsize, 'color': tc})
plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
ax.tick_params(axis='y',labelsize=fsize)
ax.tick_params(axis='x',labelsize=fsize)
plt.xticks(x_ticks,x_labels)
ax.grid(zorder=0)
plt.legend([h[0],hoft[0]],["AFT","OFT"],loc='upper right',fontsize=fsize)
fig.tight_layout()
fig.savefig('oft_vs_aft_flux_total_signed.png', bbox_inches="tight",   \
            pad_inches=0.1, dpi=dpi, facecolor=fig.get_facecolor(), \
            edgecolor=None)
plt.close('all')


#
# ****** Axial dipole.
#
fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, \
                 facecolor=fc, frameon=True)
ax = plt.gca()
h = plt.plot(jd_aft, axial_dipole_aft, \
             'k-',linewidth=lsz, markersize=mksz)
hoft = plt.plot(jd_oft, axial_dipole_oft, \
             'b-',linewidth=0.8*lsz, markersize=mksz)             
plt.xlim(xmin=time_start.jd,xmax=time_end.jd)
plt.title('Axial Dipole Strength', {'fontsize': fsize, 'color': tc})
plt.xlabel('Year', {'fontsize': fsize, 'color': tc})
plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
ax.tick_params(axis='y',labelsize=fsize)
ax.tick_params(axis='x',labelsize=fsize)
plt.xticks(x_ticks,x_labels)
ax.grid(zorder=0)
plt.legend([h[0],hoft[0]],["AFT","OFT"],loc='upper right',fontsize=fsize)
fig.tight_layout()
fig.savefig('oft_vs_aft_dipole_axial.png', bbox_inches="tight",   \
            pad_inches=0.1, dpi=dpi, facecolor=fig.get_facecolor(), \
            edgecolor=None)
plt.close('all')


#
# ****** Equatorial dipole.
#
fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, \
                 facecolor=fc, frameon=True)
ax = plt.gca()
h = plt.plot(jd_aft, eq_dipole_aft, \
             'k-',linewidth=lsz, markersize=mksz)
hoft = plt.plot(jd_oft, eq_dipole_oft, \
             'b-',linewidth=0.8*lsz, markersize=mksz)
plt.xlim(xmin=time_start.jd,xmax=time_end.jd)
plt.title('Equatorial Dipole Strength', {'fontsize': fsize, 'color': tc})
plt.xlabel('Year', {'fontsize': fsize, 'color': tc})
plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
ax.tick_params(axis='y',labelsize=fsize)
ax.tick_params(axis='x',labelsize=fsize)
plt.xticks(x_ticks,x_labels)
ax.grid(zorder=0)
plt.legend([h[0],hoft[0]],["AFT","OFT"],loc='upper right',fontsize=fsize)
fig.tight_layout()
fig.savefig('oft_vs_aft_dipole_eq.png', bbox_inches="tight",   \
            pad_inches=0.1, dpi=dpi, facecolor=fig.get_facecolor(), \
            edgecolor=None)
plt.close('all')


#
# ****** Polar average field strengths.
#
fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, \
                 facecolor=fc, frameon=True)
ax = plt.gca()
h = plt.plot(jd_aft,pole_n_aft,'k-',linewidth=lsz,markersize=mksz)
h2 = plt.plot(jd_aft,pole_s_aft,'r-',linewidth=lsz,markersize=mksz) 
hoft = plt.plot(jd_oft,pole_n_oft,'b-',linewidth=0.8*lsz,markersize=mksz)
h2oft = plt.plot(jd_oft,pole_s_oft,'m-',linewidth=0.8*lsz,markersize=mksz) 
plt.xlim(xmin=time_start.jd,xmax=time_end.jd)
plt.title('Polar Average Field (within 30 degrees of poles)', {'fontsize': fsize, 'color': tc})
plt.xlabel('Year', {'fontsize': fsize, 'color': tc})
plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
ax.tick_params(axis='y',labelsize=fsize)
ax.tick_params(axis='x',labelsize=fsize)
plt.xticks(x_ticks,x_labels)
plt.legend([h[0],h2[0],hoft[0],h2oft[0]],["AFT North","AFT South","OFT North","OFT South"],loc='upper right',fontsize=fsize)
ax.grid(zorder=0)
fig.tight_layout()
fig.savefig('oft_vs_aft_polar_fields.png', bbox_inches="tight",   \
            pad_inches=0.1, dpi=dpi, facecolor=fig.get_facecolor(), \
            edgecolor=None)
plt.close('all')















