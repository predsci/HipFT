#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

#history_sol.dat
#history_num.dat

# Add option to indicate a validation run, otherwise ignore that.
# Plot each quantity 
# Add skip option
# Plot versus time or iters?

def argParsing():
  parser = argparse.ArgumentParser(description='HipFt History Plots.')

 # parser.add_argument('file', help='Name of file')

  return parser.parse_args()


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
  flux_imb=(1e2)*flux_tot_s/(flux_tot_un + 1e-200)
  pole_n_avg_field = (fluxp_pn+fluxm_pn)/area_pn
  pole_s_avg_field = (fluxp_ps+fluxm_ps)/area_ps
  
  ###### PLOTTING ######
  
  
  width = 0.6
  fsize = 25
  mksz = 15
  lsz = 2.0
  fc = 'w'
  tc = 'k'
#
# ****** Total flux imbalance, and +/- flux.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=120, facecolor=fc,frameon=True)
  ax = plt.gca()
  h = plt.plot(time/24.0,flux_imb,'k-',linewidth=lsz,markersize=mksz)
  xmx = np.max(time/24.0)
  plt.xlim(xmin=0,xmax=xmx)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Relative Flux Imbalance ($10^{-2}$)', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax2 = ax.twinx()  
  h1 = ax2.plot(time/24.0,flux_fac*np.abs(fluxm),'-',linewidth=1.75*lsz,color='Blue',markersize=mksz)
  h2 = ax2.plot(time/24.0,flux_fac*fluxp,'-',linewidth=0.75*lsz,color='Red',markersize=mksz)
  ax2.set_ylabel('Flux ($10^{21}$ Mx)', color='Red',fontsize=fsize)
  ax2.tick_params(axis='y', labelcolor='Red',labelsize=fsize)
  ax.grid(zorder=0)
  plt.legend([h2[0],h1[0],h[0]],["Flux (+)","|Flux (-)|","Flux Imbalance"],loc='upper center',fontsize=fsize)
  fig.tight_layout()
  fig.savefig('history_flux_imb_pm.png', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  #fig.savefig('history_flux.eps', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Total signed and unsigned flux.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=120, facecolor=fc,frameon=True)
  ax = plt.gca()
  h = plt.plot(time/24.0,flux_fac*flux_tot_un,'k-',linewidth=lsz,markersize=mksz)
  xmx = np.max(time/24.0)
  plt.xlim(xmin=0,xmax=xmx)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax2 = ax.twinx()  
  h1 = ax2.plot(time/24.0,flux_fac*flux_tot_s,'--',linewidth=1.75*lsz,color='Blue',markersize=mksz)
  ax2.set_ylabel('$10^{21}$ Mx', color='Blue',fontsize=fsize)
  ax2.tick_params(axis='y', labelcolor='Blue',labelsize=fsize)
  ax.grid(zorder=0)
  plt.legend([h[0],h1[0]],["Total Unsigned Flux","Total Signed Flux"],loc='upper center',fontsize=fsize)
  fig.tight_layout()
  fig.savefig('history_flux_total_u_s.png', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  #fig.savefig('history_flux.eps', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Polar +/- fluxes.
#
  ymax = 1.1*np.amax([np.amax(np.abs(flux_fac*fluxp_pn)),np.amax(np.abs(flux_fac*fluxm_pn)),\
                      np.amax(np.abs(flux_fac*fluxp_ps)),np.amax(np.abs(flux_fac*fluxm_ps))])
  ymin = -ymax 
  fig = plt.figure(num=None, figsize=(14, 7), dpi=120, facecolor=fc,frameon=True)
  ax = plt.gca()
  xmx = np.max(time/24.0)
  h = plt.plot(time/24.0,flux_fac*fluxp_pn,'r-',linewidth=lsz,markersize=mksz)
  h2 = plt.plot(time/24.0,flux_fac*fluxm_pn,'b-',linewidth=lsz,markersize=mksz) 
  h3 = plt.plot(time/24.0,flux_fac*fluxp_ps,'r--',linewidth=lsz,markersize=mksz)
  h4 = plt.plot(time/24.0,flux_fac*fluxm_ps,'b--',linewidth=lsz,markersize=mksz)
  plt.xlim(xmin=0,xmax=xmx)
  plt.ylim(ymin=ymin,ymax=ymax)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Polar Flux ($10^{21}$ Mx)', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  plt.legend([h[0],h2[0],h3[0],h4[0]],["N (+)","N (-)","S (+)","S (-)"],loc='upper left',fontsize=fsize)
  fig.tight_layout()
  fig.savefig('history_flux_poles.png', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Polar average field strengths.
#
  ymax = 1.1*np.amax([np.amax(np.abs(pole_n_avg_field)),np.amax(np.abs(pole_s_avg_field))])
  ymin = -ymax 
  fig = plt.figure(num=None, figsize=(14, 7), dpi=120, facecolor=fc,frameon=True)
  ax = plt.gca()
  xmx = np.max(time/24.0)
  h = plt.plot(time/24.0,pole_n_avg_field,'m-',linewidth=lsz,markersize=mksz)
  h2 = plt.plot(time/24.0,pole_s_avg_field,'m--',linewidth=lsz,markersize=mksz) 
  plt.xlim(xmin=0,xmax=xmx)
  plt.ylim(ymin=ymin,ymax=ymax)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Avg polar field (Gauss)', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  plt.legend([h[0],h2[0]],["North","South"],loc='upper left',fontsize=fsize)
  fig.tight_layout()
  fig.savefig('history_field_poles.png', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Min and max field.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=120, facecolor=fc,frameon=True)
  ax = plt.gca()
  h = plt.plot(time/24.0,brmax,'b-',linewidth=lsz,markersize=mksz)
  xmx = np.max(time/24.0)
  plt.xlim(xmin=0,xmax=xmx)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  h1 = plt.plot(time/24.0,np.abs(brmin),'r-',linewidth=lsz,markersize=mksz)
  ax.grid(zorder=0)
  plt.legend([h[0],h1[0]],["max(Br)","|min(Br)|"],loc='upper center',fontsize=fsize)
  fig.tight_layout()
  fig.savefig('history_br.png', bbox_inches="tight", dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Axial Dipole strength
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=120, facecolor=fc,frameon=True)
  ax = plt.gca()
  xmx = np.max(time/24.0)
  h = plt.plot(time/24.0,ax_dipole,'k-',linewidth=lsz,markersize=mksz)
  plt.xlim(xmin=0,xmax=xmx)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Axial Dipole Strength', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  fig.tight_layout()
  fig.savefig('history_dipole_axial.png', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Equatorial Dipole strength
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=120, facecolor=fc,frameon=True)
  ax = plt.gca()
  xmx = np.max(time/24.0)
  h = plt.plot(time/24.0,eq_dipole,'k-',linewidth=lsz,markersize=mksz)
  plt.xlim(xmin=0,xmax=xmx)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Equatorial Dipole Strength', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  fig.tight_layout()
  fig.savefig('history_dipole_eq.png', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Validation
#
  # Check sum of val - if exactly 0, assume no validation run performed...
  
  fig = plt.figure(num=None, figsize=(14, 7), dpi=120, facecolor=fc,frameon=True)
  ax = plt.gca()
  h = plt.plot(time/24.0,(1e5)*valerr,'k-',linewidth=lsz,markersize=mksz)
  #plt.xlim(xmin=0,xmax=28)
  #plt.xticks((0, 4, 8, 12, 16, 20, 24, 28))
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('(CV)RMSD ($10^{-5}$)', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  fig.tight_layout()
  
  fig.savefig('history_val.png', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
#  fig.savefig('history_val.eps', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  
def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()

