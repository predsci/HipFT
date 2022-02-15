#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


def argParsing():
  parser = argparse.ArgumentParser(description='HipFt History Plots.')

  parser.add_argument('file', help='Name of file')

  return parser.parse_args()


def run(args):

  #Get base filename for automatic output names:
  filename = os.path.basename(args.file)
  file_root_str = filename[0:len(filename) - 4]
  
  # Read the data:
  hist_sol = pd.read_table(args.file,header=0,sep='\s+')

  step  = np.array(hist_sol['STEP'])
  time  = np.array(hist_sol['TIME'])
  fluxp = np.array(hist_sol['FLUX_POSITIVE'])
  fluxm = np.array(hist_sol['FLUX_NEGATIVE'])
  brmin = np.array(hist_sol['BR_MIN'])
  brmax = np.array(hist_sol['BR_MAX'])
  brabsmin = np.array(hist_sol['BR_ABS_MIN'])
  valerr = np.array(hist_sol['VALIDATION_ERR_CVRMSD'])
  
  #Compute derived quantities:
  flux_tot_un = np.abs(fluxp) + np.abs(fluxm)
  flux_tot_s = fluxp + fluxm
  
  
  ###### PLOTTING ######
  
  
  width = 0.6
  fsize = 25
  mksz = 15
  lsz = 5.0
  fc = 'w'
  tc = 'k'

  fig = plt.figure(num=None, figsize=(14, 7), dpi=120, facecolor=fc,frameon=True)
  ax = plt.gca()
  h = plt.plot(time/24.0,(1e6)*flux_tot_s/(flux_tot_un),'k-',linewidth=lsz,markersize=mksz)
  plt.xlim(xmin=0,xmax=28)
  plt.xticks((0, 4, 8, 12, 16, 20, 24, 28))
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Relative Flux Imbalance ($10^{-6}$)', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax2 = ax.twinx()  
  h1 = ax2.plot(time/24.0,abs(fluxm),'-',linewidth=1.75*lsz,color='Red',markersize=mksz)
  ax2.set_ylabel('Flux ($10^{21}$ Mx)', color='Blue',fontsize=fsize)
  h2 = ax2.plot(time/24.0,fluxp,'-',linewidth=0.75*lsz,color='Blue',markersize=mksz)
  ax2.tick_params(axis='y', labelcolor='Blue',labelsize=fsize)
  ax.grid(zorder=0)
  plt.legend([h2[0],h1[0],h[0]],["Flux (+)","Flux (-)","Flux Imbalance"],loc='upper center',fontsize=fsize)
  fig.tight_layout()

  fig.savefig(file_root_str+'_flux.png', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  fig.savefig(file_root_str+'_flux.eps', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  
  plt.close('all')
  
  # VALIDATION
  
  fig = plt.figure(num=None, figsize=(14, 7), dpi=120, facecolor=fc,frameon=True)
  ax = plt.gca()

  h = plt.plot(time/24.0,(1e5)*valerr,'k-',linewidth=lsz,markersize=mksz)
  plt.xlim(xmin=0,xmax=28)
  plt.xticks((0, 4, 8, 12, 16, 20, 24, 28))
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('(CV)RMSD ($10^{-5}$)', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.grid(zorder=0)
  fig.tight_layout()
  
  fig.savefig(file_root_str+'_valerr.png', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  fig.savefig(file_root_str+'_valerr.eps', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)
  

def main():
  ## Get input agruments:
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()












