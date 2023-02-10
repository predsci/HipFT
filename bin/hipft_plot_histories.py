#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
import os

#history_sol.dat
#history_num.dat

# - Add skip option
# - Add initial UT time and plot by year, MM/YY, etc.
# - Add ability to compare two runs...


def argParsing():
  parser = argparse.ArgumentParser(description='HipFt History Plots.')

  parser.add_argument('-valrun',
    help='Plot validation run errors',
    dest='valrun',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-utstart',
    help='Start Date in UT: YYYYMMDDTHH:MM:SS',
    dest='utstart',
    default='00000000T00:00:00',
    required=False)

  parser.add_argument('-xaxunitstr',
    help='Plot axis units',
    dest='xaxunitstr',
    default='days',
    required=False)

  return parser.parse_args()


###  Need to make Time object, convert to datetime, use %s in order ot get "unix time"
###


def run(args):

  flux_fac=1e-21
  area_fac=1e-21

  # Read the data:
  hist_sol = pd.read_table('history_sol.dat',header=0,sep='\s+')
  hist_num = pd.read_table('history_num.dat',header=0,sep='\s+')
  
  step  = np.array(hist_sol['STEP'])
  time  = np.array(hist_sol['TIME'])
  fluxp = np.array(hist_sol['FLUX_POSITIVE'])
  fluxm = np.array(hist_sol['FLUX_NEGATIVE'])
  
  fluxp_pn = np.array(hist_sol['NPOLE_FLUX_POSITIVE'])
  fluxm_pn = np.array(hist_sol['NPOLE_FLUX_NEGATIVE'])
  area_pn = np.array(hist_sol['NPOLE_AREA'])
  
  fluxp_ps = np.array(hist_sol['SPOLE_FLUX_POSITIVE'])
  fluxm_ps = np.array(hist_sol['SPOLE_FLUX_NEGATIVE'])
  area_ps = np.array(hist_sol['SPOLE_AREA'])
  
  ax_dipole = np.array(hist_sol['AX_DIPOLE'])
  eq_dipole = np.array(hist_sol['EQ_DIPOLE'])

  brmin = np.array(hist_sol['BR_MIN'])
  brmax = np.array(hist_sol['BR_MAX'])
  brabsmin = np.array(hist_sol['BR_ABS_MIN'])
  valerr = np.array(hist_sol['VALIDATION_ERR_CVRMSD'])

  dtime = np.array(hist_num['DTIME'])
  dtime_advection_stable = np.array(hist_num['DTIME_ADV_STB'])
  dtime_advection_used = np.array(hist_num['DTIME_ADV_USED'])
  dtime_diffusion_stable = np.array(hist_num['DTIME_DIFF_STB'])
  dtime_diffusion_used = np.array(hist_num['DTIME_DIFF_USED'])
  n_sts_iter_per_step = np.array(hist_num['N_DIFF_PER_STEP'])

  #Compute derived quantities:
  flux_tot_un = np.abs(fluxp) + np.abs(fluxm)
  flux_tot_s = fluxp + fluxm
  flux_imb = flux_tot_un*0.0
  flux_imb[flux_tot_un < 1e-15] = 0.0
  flux_imb[flux_tot_un >= 1e-15] = 100.0*flux_tot_s[flux_tot_un >= 1e-15]/flux_tot_un[flux_tot_un >= 1e-15]
  pole_n_avg_field = (fluxp_pn+fluxm_pn)/area_pn
  pole_s_avg_field = (fluxp_ps+fluxm_ps)/area_ps
  
  ###### PLOTTING ######
  
  width = 0.6
  fsize = 25
  mksz = 15
  lsz = 3.0
  fc = 'w'
  tc = 'k'
  dpi=120
  
  tfac = 1/24.0
  
#
# ****** Total flux imbalance.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  h = plt.plot(time*tfac,flux_imb,'k-',linewidth=lsz,markersize=mksz)
  xmn = np.min(time*tfac)
  xmx = np.max(time*tfac)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.title('Relative Flux Imbalance', {'fontsize': fsize, 'color': tc})
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('%', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  fig.tight_layout()
  fig.savefig('history_flux_imb_pm.png', bbox_inches="tight", pad_inches=0, dpi=dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Total (+) and (-) flux.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  h = ax.plot(time*tfac,flux_fac*np.abs(fluxm),'-',linewidth=lsz,color='Blue',markersize=mksz)
  h2 = ax.plot(time*tfac,flux_fac*fluxp,'-',linewidth=lsz,color='Red',markersize=mksz)
  xmn = np.min(time*tfac)
  xmx = np.max(time*tfac)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.title('Total Positive and Negative Flux', {'fontsize': fsize, 'color': tc})
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  plt.legend([h[0],h2[0]],["|Flux (-)|","Flux (+)"],loc='upper right',fontsize=fsize)
  fig.tight_layout()
  fig.savefig('history_flux_pm.png', bbox_inches="tight", pad_inches=0, dpi=dpi, \
              facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Total unsigned flux.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  h = plt.plot(time*tfac,flux_fac*flux_tot_un,'k-',linewidth=lsz,markersize=mksz)
  xmn = np.min(time*tfac)
  xmx = np.max(time*tfac)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.title('Total Unsigned Flux', {'fontsize': fsize, 'color': tc})
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  fig.tight_layout()
  fig.savefig('history_flux_total_unsigned.png', bbox_inches="tight", pad_inches=0, dpi=dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Total signed flux.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  h = plt.plot(time*tfac,flux_fac*flux_tot_s,'b-',linewidth=lsz,markersize=mksz)
  xmn = np.min(time*tfac)
  xmx = np.max(time*tfac)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.title('Total Signed Flux', {'fontsize': fsize, 'color': tc})
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  fig.tight_layout()
  fig.savefig('history_flux_total_signed.png', bbox_inches="tight", pad_inches=0, dpi=dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')  
#
# ****** Polar +/- fluxes.
#
  ymax = 1.1*np.amax([np.amax(np.abs(flux_fac*fluxp_pn)),np.amax(np.abs(flux_fac*fluxm_pn)),\
                      np.amax(np.abs(flux_fac*fluxp_ps)),np.amax(np.abs(flux_fac*fluxm_ps))])
  ymin = -ymax 
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  xmn = np.min(time*tfac)
  xmx = np.max(time*tfac)
  h = plt.plot(time*tfac,flux_fac*fluxp_pn,'r-',linewidth=lsz,markersize=mksz)
  h2 = plt.plot(time*tfac,flux_fac*fluxm_pn,'b-',linewidth=lsz,markersize=mksz) 
  h3 = plt.plot(time*tfac,flux_fac*fluxp_ps,'-',color='firebrick',linewidth=lsz,markersize=mksz)
  h4 = plt.plot(time*tfac,flux_fac*fluxm_ps,'-',color='navy',linewidth=lsz,markersize=mksz)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.ylim(ymin=ymin,ymax=ymax)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
  plt.title('Polar Flux (within 30 degrees of poles)', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  plt.legend([h[0],h2[0],h3[0],h4[0]],["N (+)","N (-)","S (+)","S (-)"],loc='upper right',fontsize=fsize)
  fig.tight_layout()
  fig.savefig('history_flux_poles_30.png', bbox_inches="tight", pad_inches=0, dpi=dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Polar average field strengths.
#
  ymax = 1.1*np.amax([np.amax(np.abs(pole_n_avg_field)),np.amax(np.abs(pole_s_avg_field))])
  ymin = -ymax 
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  xmn = np.min(time*tfac)
  xmx = np.max(time*tfac)
  h = plt.plot(time*tfac,pole_n_avg_field,'k-',linewidth=lsz,markersize=mksz)
  h2 = plt.plot(time*tfac,pole_s_avg_field,'b-',linewidth=lsz,markersize=mksz) 
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.ylim(ymin=ymin,ymax=ymax)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  plt.title('Polar Average Field (within 30 degrees of poles)', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  plt.legend([h[0],h2[0]],["North","South"],loc='upper right',fontsize=fsize)
  fig.tight_layout()
  fig.savefig('history_field_poles_30.png', bbox_inches="tight", pad_inches=0, dpi=dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Min and max field.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  h = plt.plot(time*tfac,brmax,'b-',linewidth=lsz,markersize=mksz)
  xmn = np.min(time*tfac)
  xmx = np.max(time*tfac)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  plt.title('Min and Max Br', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  h1 = plt.plot(time*tfac,np.abs(brmin),'r-',linewidth=lsz,markersize=mksz)
  ax.grid(zorder=0)
  plt.legend([h[0],h1[0]],["max(Br)","|min(Br)|"],loc='upper right',fontsize=fsize)
  fig.tight_layout()
  fig.savefig('history_br.png', bbox_inches="tight", dpi=dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Axial Dipole strength
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  xmn = np.min(time*tfac)
  xmx = np.max(time*tfac)
  h = plt.plot(time*tfac,ax_dipole,'k-',linewidth=lsz,markersize=mksz)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.title('Axial Dipole Strength', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  fig.tight_layout()
  fig.savefig('history_dipole_axial.png', bbox_inches="tight", pad_inches=0, dpi=dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Equatorial Dipole strength
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  xmn = np.min(time*tfac)
  xmx = np.max(time*tfac)
  h = plt.plot(time*tfac,eq_dipole,'k-',linewidth=lsz,markersize=mksz)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.title('Equatorial Dipole Strength', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  fig.tight_layout()
  fig.savefig('history_dipole_eq.png', bbox_inches="tight", pad_inches=0, dpi=dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
  
#
# ****** Validation
#
  if (args.valrun):
    fig = plt.figure(num=None, figsize=(14, 7), dpi=120, facecolor=fc,frameon=True)
    ax = plt.gca()
    h = plt.plot(time*tfac,(1e5)*valerr,'k-',linewidth=lsz,markersize=mksz)
    plt.title('Validation Error', {'fontsize': fsize, 'color': tc})
    plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
    plt.ylabel('(CV)RMSD ($10^{-5}$)', {'fontsize': fsize, 'color': tc})
    ax.tick_params(axis='y',labelsize=fsize)
    ax.tick_params(axis='x',labelsize=fsize)
    ax.grid(zorder=0)
    fig.tight_layout()  
    fig.savefig('history_val.png', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  
def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()

