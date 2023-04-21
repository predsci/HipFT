#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
import os

def argParsing():
  parser = argparse.ArgumentParser(description='HipFt History Plots.')
  parser.add_argument('-name',help='Name of the set of plots',type=str, default='')
  parser.add_argument('-valrun', help='Plot validation run errors', dest='valrun', action='store_true', default=False, required=False)
  parser.add_argument('-samples',help='Number of points to plot, this helps with larger files (default:all)',default=-1)
  parser.add_argument('-dirs',help='A comma separated list of paths to the directory each run is located in', type=str, required=True)
  parser.add_argument('-labels',help='A comma separated list of labels to display in the plots for each run (default:"Run 1","Run 2",...)',type=str,default=' ')

  return parser.parse_args()


###  Need to make Time object, convert to datetime, use %s in order ot get "unix time"
###


def run(args):

  arg_dict = vars(args)
  dir_list = arg_dict['dirs'].split(',')
  label_list = arg_dict['labels'].split(',')

  # Build the default list of labels based on the number runs the user entered
  # if no label list is entered
  if arg_dict['labels'] == ' ':
      label_list = []
      def_label = 'Run '
      for i,dire in enumerate(dir_list):
          label_list.append(def_label+str(i+1))

  LABEL_LEN = len(label_list)
  # Validate the list arguments:
  if not len(dir_list) == LABEL_LEN:
      print('ERROR: Number of runs, dirs, and labels must match. Use -h for more information.')
      sys.exit()

  if LABEL_LEN > 5:
      print('ERROR: Can only compare a maximum of 5 runs.')

  print("==> Reading history files...")
  time_list=[]
  fluxp_list=[]
  fluxm_list=[]
  fluxp_pn_list=[]
  fluxm_pn_list=[]
  fluxp_ps_list=[]
  fluxm_ps_list=[]
  ax_dipole_list=[]
  eq_dipole_list=[]
  brmin_list=[]
  brmax_list=[]
  valerr_list=[]
  flux_tot_un_list=[]
  flux_tot_s_list=[]
  flux_imb_list=[]
  pole_n_avg_field_list=[]
  pole_s_avg_field_list=[]

  for i,dire in enumerate(dir_list):
    h_file_name = "{}/{}".format(dire,'history_sol.dat')
    hist_sol_full = pd.read_table(h_file_name,header=0,sep='\s+')

    samples = int(args.samples)

    number_of_data_points = len(hist_sol_full.iloc[:,0])

    if samples == -1:
        skip = 1
    else:
        skip = int(np.floor(np.amax([1,number_of_data_points/samples])))

    #thin out data in hist_sol
    hist_sol = hist_sol_full[::skip]

    time_list.append(np.array(hist_sol['TIME']))
    fluxp_list.append(np.array(hist_sol['FLUX_POSITIVE']))
    fluxm_list.append(np.array(hist_sol['FLUX_NEGATIVE']))
  
    fluxp_pn_list.append(np.array(hist_sol['NPOLE_FLUX_POSITIVE']))
    fluxm_pn_list.append(np.array(hist_sol['NPOLE_FLUX_NEGATIVE']))
    area_pn = np.array(hist_sol['NPOLE_AREA'])
  
    fluxp_ps_list.append(np.array(hist_sol['SPOLE_FLUX_POSITIVE']))
    fluxm_ps_list.append(np.array(hist_sol['SPOLE_FLUX_NEGATIVE']))
    area_ps = np.array(hist_sol['SPOLE_AREA'])
  
    ax_dipole_list.append(np.array(hist_sol['AX_DIPOLE']))
    eq_dipole_list.append(np.array(hist_sol['EQ_DIPOLE']))

    brmin_list.append(np.array(hist_sol['BR_MIN']))
    brmax_list.append(np.array(hist_sol['BR_MAX']))
    valerr_list.append(np.array(hist_sol['VALIDATION_ERR_CVRMSD']))

    #Compute derived quantities:
    flux_tot_un = np.abs(fluxp_list[i]) + np.abs(fluxm_list[i])
    flux_tot_un_list.append(flux_tot_un)
    flux_tot_s = fluxp_list[i] + fluxm_list[i]
    flux_tot_s_list.append(flux_tot_s)
    flux_imb = flux_tot_un*0.0
    flux_imb[flux_tot_un < 1e-15] = 0.0
    flux_imb[flux_tot_un >= 1e-15] = 100.0*flux_tot_s[flux_tot_un >= 1e-15]/flux_tot_un[flux_tot_un >= 1e-15]
    flux_imb_list.append(flux_imb)
    pole_n_avg_field_list.append((fluxp_pn_list[i]+fluxm_pn_list[i])/area_pn)
    pole_s_avg_field_list.append((fluxp_ps_list[i]+fluxm_ps_list[i])/area_ps)


  COLORS = ['#0000ff','#ff0000','#008000','#800080','#b79f00']
  COLOR_NAMES = ['Blue', 'Red', 'Green', 'Purple', 'Gold']
  MARKERS = ['s','o','d','^','>']
  MARKERS_NAMES = ['Sq','Cir','Dia','^','>']


  title = arg_dict['name']
  title2 = arg_dict['name']
  if LABEL_LEN > 1:
      for i, label in enumerate(label_list[:-1]):
          title2 = title2 + ' ' + '{} ({})'.format(label, MARKERS_NAMES[i])
          title2 = title2 + ' ' + ' vs. '
          title = title + ' ' + '{} ({})'.format(label, COLOR_NAMES[i])
          title = title + ' ' + ' vs. '

      title = title + '{} ({})'.format(label_list[-1], COLOR_NAMES[LABEL_LEN - 1])
      title2 = title2 + '{} ({})'.format(label_list[-1], MARKERS_NAMES[LABEL_LEN - 1])

  flux_fac=1e-21

  ###### PLOTTING ######

  width = 0.3
  fsize = 25
  MS = 6
  LW = 1.0
  LWM = LW/LABEL_LEN
  fc = 'w'
  tc = 'k'
  dpi=120

  tfac = 1/24.0

#
# ****** Total flux imbalance.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  fig.suptitle(title, fontsize=20)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    plt.plot(time_list[run]*tfac, flux_imb_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)
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
  fig.suptitle(title2, fontsize=20)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    h=plt.plot(time_list[run]*tfac, flux_fac*np.abs(fluxm_list[run]),color='Blue',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor='Blue')
    h2=plt.plot(time_list[run]*tfac, flux_fac*fluxp_list[run],color='Red',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor='Red')

  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)
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
  fig.suptitle(title, fontsize=20)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    plt.plot(time_list[run]*tfac, flux_fac*flux_tot_un_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)
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
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  fig.suptitle(title, fontsize=20)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    plt.plot(time_list[run]*tfac, flux_fac*flux_tot_s_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)
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
  ymax = 1.1*np.amax([np.amax(np.abs(flux_fac*np.array(fluxp_pn_list))),np.amax(np.abs(flux_fac*np.array(fluxm_pn_list))),\
                      np.amax(np.abs(flux_fac*np.array(fluxp_ps_list))),np.amax(np.abs(flux_fac*np.array(fluxm_ps_list)))])
  ymin = -ymax 
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  fig.suptitle(title2, fontsize=20)
  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    h=plt.plot(time_list[run]*tfac, flux_fac*fluxp_pn_list[run],color='Red',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor='Red')
    h2=plt.plot(time_list[run]*tfac, flux_fac*fluxm_pn_list[run],color='Blue',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor='Blue')
    h3=plt.plot(time_list[run]*tfac, flux_fac*fluxp_ps_list[run],color='firebrick',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor='firebrick')
    h4=plt.plot(time_list[run]*tfac, flux_fac*fluxm_ps_list[run],color='navy',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor='navy')

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
  ymax = 1.1*np.amax([np.amax(np.abs(np.array(pole_n_avg_field_list))),np.amax(np.abs(np.array(pole_s_avg_field_list)))])
  ymin = -ymax 
  fig = plt.figure(num=None, figsize=(14, 7), dpi=dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  fig.suptitle(title2, fontsize=20)
  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    h=plt.plot(time_list[run]*tfac, pole_n_avg_field_list[run],color='Black',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor='Black')
    h2=plt.plot(time_list[run]*tfac, pole_s_avg_field_list[run],color='Blue',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor='Blue')

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
  fig.suptitle(title2, fontsize=20)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    h=plt.plot(time_list[run]*tfac, brmax_list[run],color='blue',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor='blue')

  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.xlabel('Time (Days)', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  plt.title('Min and Max Br', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    h1=plt.plot(time_list[run]*tfac, np.abs(brmin_list[run]),color='red',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor='red')
  
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
  fig.suptitle(title, fontsize=20)
  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    plt.plot(time_list[run]*tfac, ax_dipole_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

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
  fig.suptitle(title, fontsize=20)
  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    plt.plot(time_list[run]*tfac, eq_dipole_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

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
    fig.suptitle(title, fontsize=20)

    for run in range(len(time_list)):
      LWu = LW-run*LWM
      plt.plot(time_list[run]*tfac, (1e5)*valerr_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
          markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

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
