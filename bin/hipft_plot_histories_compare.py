#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime, timezone
from sunpy.coordinates.sun import carrington_rotation_time, carrington_rotation_number
import os
import itertools

# Version 1.6.1

def argParsing():
  parser = argparse.ArgumentParser(description='HipFt History Plots.')

  parser.add_argument('-runtag',
    help='Add a run tag to the output plot file names.',
    dest='runtag',
    default='',
    required=False)

  parser.add_argument('-samples',
    help='Number of points to plot, this helps with larger files (default:all)',
    default=-1)

  parser.add_argument('-histfiles',
    help='A comma separated list of history files',
    type=str,
    required=False,
    default=' ')

  parser.add_argument('-labels',
    help='A comma separated list of labels to display in the plots for each run (default:"r1","r2",...)',
    type=str,
    default=' ')

  parser.add_argument('-val',
    help='Plot validation run errors.',
    dest='valrun',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-dpi',
    help='DPI for plots.',
    dest='dpi',
    default=120,
    type=int,
    required=False)

  parser.add_argument('-utstart',
    help='Start Date in UT: YYYY-MM-DDTHH:MM:SS',
    dest='utstart',
    required=False)
    
  parser.add_argument('-utstartxtick',
    help='UT date of first xtick: YYYY-MM-DDTHH:MM:SS',
    dest='utstartxtick',
    required=False)
    
  parser.add_argument('-tfac',
    help='Time factor to get time into units of hours (Default is 1.0).',
    dest='tfac',
    default='1.0',
    required=False)

  parser.add_argument('-xunits',
    help='Units of the x-axis (date, seconds, minutes, hours, days, weeks, cr, or years).',
    dest='xunits',
    default='hours',
    required=False)

  parser.add_argument('-xcrpos',
    help='Carrington rotation marker location (start, center, end).',
    dest='xcrpos',
    default='start',
    required=False)

  parser.add_argument('-xformat',
    help='Format for the date option where // is treated as a newline (example: %%H:%%M:%%S//%%Y/%%m/%%d).',
    dest='xformat',
    default='%H:%M:%S//%Y/%m/%d',
    required=False)

  parser.add_argument('-xcadence',
    help='Cadence of xc_units for xaxis tick marks (default is matplotlib auto-ticks cadence).',
    dest='xcadence',
    required=False)

  parser.add_argument('-xcadence_units',
    help='Units for x-axis tick cadence (default is xunits) (seconds, minutes, hours, days, weeks, cr, months, or years).',
    dest='xc_units',
    required=False)

  parser.add_argument('-xslantdeg',
    help='Flag to slant x-axis labels by a desired degree.',
    dest='slant',
    required=False)

  parser.add_argument('-xha',
    help='Horizontal alignment of slanted x-axis labels (left, center, or right) (default:center).',
    dest='ha',
    default='center',
    required=False)

  parser.add_argument('-xrmode', 
    action='store_true',
    help='Set text rotation mode to "anchor".',
    dest='rmode',
    required=False)

  parser.add_argument('-xma',
    help='Multi line alignment of slanted x-axis labels (left, center, or right) (default:center).',
    dest='ma',
    default='center',
    required=False)

  parser.add_argument('-fsize',
    help='Font size',
    dest='fsize',
    default=25,
    required=False)

  parser.add_argument('-fsize_xticklabels',
    help='Font size for xtick labels',
    dest='xlabelfsize',
    default=25,
    required=False)

  parser.add_argument('-xlabel',
    help='Label for x axis',
    dest='xlabel',
    required=False)

  parser.add_argument('-lgwidth',
    help='Legend box width (default:1)',
    dest='lgwidth',
    default=1,
    type=float,
    required=False)

  parser.add_argument('-lgyoffset',
    help='Legend box height offset (default:-0.1)',
    dest='lgyoffset',
    default=-0.1,
    type=float,
    required=False)

  parser.add_argument('-lgfsize',
    help='Legend font size',
    dest='lgfsize',
    default=25,
    type=float,
    required=False)

  parser.add_argument('-lw',
    help='Line width',
    dest='lw',
    type=float,
    required=False)

  parser.add_argument('-ms',
    help='Marker size',
    dest='ms',
    type=float,
    default=6.0,
    required=False)

  return parser.parse_args()


def stats(data):

  data_stats = np.zeros(4)

  data_stats[0] = np.min(data)
  data_stats[1] = np.max(data)
  data_stats[2] = np.mean(data)
  data_stats[3] = np.std(data)

  return data_stats


def run(args):  

  flux_fac=1e-21

  arg_dict = vars(args)
  hist_list = arg_dict['histfiles'].split(',')
  label_list = arg_dict['labels'].split(',')

  if arg_dict['histfiles'] == ' ':
    hist_list=[]
    wDir=os.getcwd()
    for file in os.listdir(wDir):
      if "hipft_history_sol_r" in file and file.endswith(".out"):
        hist_list.append(wDir+'/'+file)
    hist_list = sorted(hist_list)

  # Build the default list of labels based on the number runs the user entered
  # if no label list is entered
  if arg_dict['labels'] == ' ':
      label_list = []
      def_label = 'r'
      for i,dire in enumerate(hist_list):
          label_list.append(def_label+str(i+1))

  LABEL_LEN = len(label_list)
  # Validate the list arguments:
  if not len(hist_list) == LABEL_LEN:
      print('ERROR: Number of runs, dirs, and labels must match. Use -h for more information.')
      sys.exit()

  if LABEL_LEN > 256:
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

  ###### PLOTTING ######

  width = 0.3
  fsize = args.fsize
  lgfsize = args.lgfsize
  MS = args.ms
  LW = 1.0
  LWM = LW/LABEL_LEN
  fc = 'w'
  tc = 'k'

  ###### TIME ######

  if '/' in args.tfac:
    num, denom = args.tfac.split('/')
    tfac = float(num) / float(denom)
  else:  
    tfac = float(args.tfac)

  ######################

  for i,dire in enumerate(hist_list):
    h_file_name = dire
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


  cmap = plt.get_cmap('turbo',LABEL_LEN)
  COLORS = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
  MARKERS = ['o','v','^','<','>','8','s','p','*','h','H','D','d','P','X']
  MARKERS=MARKERS*int(np.ceil(LABEL_LEN/len(MARKERS)))

#
# ****** Total flux imbalance.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    if (args.lw):
      LWu=args.lw
    plt.plot(time_list[run]*tfac, flux_imb_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.title('Relative Flux Imbalance', {'fontsize': fsize, 'color': tc})
  
  init_locs = plt.xticks()[0]
  locs, labels, utstartSecs = get_xticks(args,xmn,xmx,init_locs)
  xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)

  plt.ylabel('%', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.grid(zorder=0)

  renderer = fig.canvas.get_renderer()
  for ncol in range(1,LABEL_LEN+1):
    lg = fig.legend(label_list,loc='outside lower center', ncol=ncol,fontsize=lgfsize)
    fig.canvas.draw()
    lgbbox = lg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    lg.remove()
    if lgbbox.width > args.lgwidth:
      ncol -= 1   
      break

  lg = fig.legend(label_list,loc='outside lower center', bbox_to_anchor=(0.5,args.lgyoffset), ncol=ncol,fontsize=lgfsize)       
  fig.tight_layout()
  fig.savefig('history_'+args.runtag+'_flux_imb_pm.png', bbox_inches='tight', dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Total (+) and (-) flux.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  firstRun=True
  for run in range(len(time_list)):
    LWu = LW-run*LWM
    if (args.lw):
      LWu=args.lw
    if firstRun:
      h=plt.plot(time_list[run]*tfac, flux_fac*np.abs(fluxm_list[run]),color='Blue',linewidth=LWu,marker=MARKERS[run],
          markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Blue')
      h2=plt.plot(time_list[run]*tfac, flux_fac*fluxp_list[run],color='Red',linewidth=LWu,marker=MARKERS[run],
          markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Red',label='_nolegend_')
      firstRun=False
    else:
      plt.plot(time_list[run]*tfac, flux_fac*np.abs(fluxm_list[run]),color='Blue',linewidth=LWu,marker=MARKERS[run],
          markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Blue')
      plt.plot(time_list[run]*tfac, flux_fac*fluxp_list[run],color='Red',linewidth=LWu,marker=MARKERS[run],
          markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Red',label='_nolegend_')

  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.title('Total Positive and Negative Flux', {'fontsize': fsize, 'color': tc})
  xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.grid(zorder=0)
  
  legend1 = plt.legend([h[0],h2[0]],["|Flux (-)|","Flux (+)"],loc='upper right',fontsize=fsize)
  renderer = fig.canvas.get_renderer()

  for ncol in range(1,LABEL_LEN+1):
    lg = fig.legend(label_list,loc='outside lower center', ncol=ncol,fontsize=lgfsize)
    fig.canvas.draw()
    lgbbox = lg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    lg.remove()
    if lgbbox.width > args.lgwidth:
      ncol -= 1      
      break

  lg = fig.legend(label_list,loc='outside lower center', bbox_to_anchor=(0.5,args.lgyoffset),ncol=ncol,fontsize=lgfsize)
  fig.tight_layout()
  ax.add_artist(legend1)  
  
  fig.savefig('history_'+args.runtag+'_flux_pm.png', bbox_inches="tight", dpi=args.dpi, \
              facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Total unsigned flux.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    if (args.lw):
      LWu=args.lw
    plt.plot(time_list[run]*tfac, flux_fac*flux_tot_un_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.title('Total Unsigned Flux', {'fontsize': fsize, 'color': tc})
  xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.grid(zorder=0)

  renderer = fig.canvas.get_renderer()
  for ncol in range(1,LABEL_LEN+1):
    lg = fig.legend(label_list,loc='outside lower center', ncol=ncol,fontsize=lgfsize)
    fig.canvas.draw()
    lgbbox = lg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    lg.remove()
    if lgbbox.width > args.lgwidth:
      ncol -= 1        
      break

  lg = fig.legend(label_list,loc='outside lower center', bbox_to_anchor=(0.5,args.lgyoffset),ncol=ncol,fontsize=lgfsize)
  fig.tight_layout()
  fig.savefig('history_'+args.runtag+'_flux_total_unsigned.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Total signed flux.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    if (args.lw):
      LWu=args.lw
    plt.plot(time_list[run]*tfac, flux_fac*flux_tot_s_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)
  plt.xlim(xmin=xmn,xmax=xmx)
  plt.title('Total Signed Flux', {'fontsize': fsize, 'color': tc})
  xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.grid(zorder=0)

  renderer = fig.canvas.get_renderer()
  for ncol in range(1,LABEL_LEN+1):
    lg = fig.legend(label_list,loc='outside lower center', ncol=ncol,fontsize=lgfsize)
    fig.canvas.draw()
    lgbbox = lg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    lg.remove()
    if lgbbox.width > args.lgwidth:
      ncol -= 1       
      break

  lg = fig.legend(label_list,loc='outside lower center', bbox_to_anchor=(0.5,args.lgyoffset),ncol=ncol,fontsize=lgfsize)
  fig.tight_layout()
  fig.savefig('history_'+args.runtag+'_flux_total_signed.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')  
#
# ****** Polar +/- fluxes.
#
  ymax = 1.1*np.amax([np.amax(np.abs(flux_fac*np.array(fluxp_pn_list))),np.amax(np.abs(flux_fac*np.array(fluxm_pn_list))),\
                      np.amax(np.abs(flux_fac*np.array(fluxp_ps_list))),np.amax(np.abs(flux_fac*np.array(fluxm_ps_list)))])
  ymin = -ymax 
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)

  firstRun=True
  for run in range(len(time_list)):
    LWu = LW-run*LWM
    if (args.lw):
      LWu=args.lw
    if firstRun:
      h=plt.plot(time_list[run]*tfac, flux_fac*fluxp_pn_list[run],color='Red',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Red')
      h2=plt.plot(time_list[run]*tfac, flux_fac*fluxm_pn_list[run],color='Blue',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Blue',label='_nolegend_')
      h3=plt.plot(time_list[run]*tfac, flux_fac*fluxp_ps_list[run],color='firebrick',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='firebrick',label='_nolegend_')
      h4=plt.plot(time_list[run]*tfac, flux_fac*fluxm_ps_list[run],color='navy',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.3,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='navy',label='_nolegend_')
      firstRun=False
    else:
      plt.plot(time_list[run]*tfac, flux_fac*fluxp_pn_list[run],color='Red',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Red')
      plt.plot(time_list[run]*tfac, flux_fac*fluxm_pn_list[run],color='Blue',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Blue',label='_nolegend_')
      plt.plot(time_list[run]*tfac, flux_fac*fluxp_ps_list[run],color='firebrick',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='firebrick',label='_nolegend_')
      plt.plot(time_list[run]*tfac, flux_fac*fluxm_ps_list[run],color='navy',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='navy',label='_nolegend_')

  plt.xlim(xmin=xmn,xmax=xmx)
  plt.ylim(ymin=ymin,ymax=ymax)
  xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
  plt.title('Polar Flux (within 30 degrees of poles)', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.grid(zorder=0)

  legend1 = plt.legend([h[0],h2[0],h3[0],h4[0]],["N (+)","N (-)","S (+)","S (-)"],loc='upper right',fontsize=fsize)
  renderer = fig.canvas.get_renderer()

  for ncol in range(1,LABEL_LEN+1):
    lg = fig.legend(label_list,loc='outside lower center', ncol=ncol,fontsize=lgfsize)
    fig.canvas.draw()
    lgbbox = lg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    lg.remove()
    if lgbbox.width > args.lgwidth:
      ncol -= 1        
      break

  lg = fig.legend(label_list,loc='outside lower center', bbox_to_anchor=(0.5,args.lgyoffset),ncol=ncol,fontsize=lgfsize)
  ax.add_artist(legend1)  

  fig.tight_layout()
  fig.savefig('history_'+args.runtag+'_flux_poles_30.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Polar average field strengths.
#
  ymax = 1.1*np.amax([np.amax(np.abs(np.array(pole_n_avg_field_list))),np.amax(np.abs(np.array(pole_s_avg_field_list)))])
  ymin = -ymax 
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)

  firstRun=True
  for run in range(len(time_list)):
    LWu = LW-run*LWM
    if (args.lw):
      LWu=args.lw
    if firstRun:
      h=plt.plot(time_list[run]*tfac, pole_n_avg_field_list[run],color='Black',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Black')
      h2=plt.plot(time_list[run]*tfac, pole_s_avg_field_list[run],color='Blue',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Blue',label='_nolegend_')
      firstRun=False
    else:
      plt.plot(time_list[run]*tfac, pole_n_avg_field_list[run],color='Black',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Black')
      plt.plot(time_list[run]*tfac, pole_s_avg_field_list[run],color='Blue',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='Blue',label='_nolegend_')

  plt.xlim(xmin=xmn,xmax=xmx)
  plt.ylim(ymin=ymin,ymax=ymax)
  xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  plt.title('Polar Average Field (within 30 degrees of poles)', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.grid(zorder=0)
  
  legend1 = plt.legend([h[0],h2[0]],["North","South"],loc='upper right',fontsize=fsize)
  renderer = fig.canvas.get_renderer()
  for ncol in range(1,LABEL_LEN+1):
    lg = fig.legend(label_list,loc='outside lower center', ncol=ncol,fontsize=lgfsize)
    fig.canvas.draw()
    lgbbox = lg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    lg.remove()
    if lgbbox.width > args.lgwidth:
      ncol -= 1
      break

  lg = fig.legend(label_list,loc='outside lower center', bbox_to_anchor=(0.5,args.lgyoffset), ncol=ncol,fontsize=lgfsize)  
  ax.add_artist(legend1)  

  fig.tight_layout()
  fig.savefig('history_'+args.runtag+'_field_poles_30.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Min and max field.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  firstRun=True
  for run in range(len(time_list)):
    LWu = LW-run*LWM
    if (args.lw):
      LWu=args.lw
    if firstRun:
      h=plt.plot(time_list[run]*tfac, brmax_list[run],color='blue',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='blue')
      firstRun=False
    else:
      plt.plot(time_list[run]*tfac, brmax_list[run],color='blue',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=COLORS[run],markeredgecolor='blue')

  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)
  plt.xlim(xmin=xmn,xmax=xmx)
  xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  plt.title('Min and Max Br', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)

  firstRun=True
  for run in range(len(time_list)):
    LWu = LW-run*LWM
    if (args.lw):
      LWu=args.lw
    if firstRun:
      h1=plt.plot(time_list[run]*tfac, np.abs(brmin_list[run]),color='red',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full', markerfacecolor=COLORS[run],markeredgecolor='red',label='_nolegend_')
      firstRun=False
    else:
      plt.plot(time_list[run]*tfac, np.abs(brmin_list[run]),color='red',linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full', markerfacecolor=COLORS[run],markeredgecolor='red',label='_nolegend_')
  
  ax.grid(zorder=0)
  
  legend1 = plt.legend([h[0],h1[0]],["max(Br)","|min(Br)|"],loc='upper right',fontsize=fsize)
  renderer = fig.canvas.get_renderer()

  for ncol in range(1,LABEL_LEN+1):
    lg = fig.legend(label_list,loc='outside lower center', ncol=ncol,fontsize=lgfsize)
    fig.canvas.draw()
    lgbbox = lg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    lg.remove()
    if lgbbox.width > args.lgwidth:
      ncol -= 1       
      break

  lg = fig.legend(label_list,loc='outside lower center', bbox_to_anchor=(0.5,args.lgyoffset), ncol=ncol,fontsize=lgfsize) 
  ax.add_artist(legend1)  

  fig.tight_layout()
  fig.savefig('history_'+args.runtag+'_br.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Axial Dipole strength
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()
  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    if (args.lw):
      LWu=args.lw
    plt.plot(time_list[run]*tfac, ax_dipole_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

  plt.xlim(xmin=xmn,xmax=xmx)
  xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)
  plt.title('Axial Dipole Strength', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.grid(zorder=0)

  renderer = fig.canvas.get_renderer()
  for ncol in range(1,LABEL_LEN+1):
    lg = fig.legend(label_list,loc='outside lower center', ncol=ncol,fontsize=lgfsize)
    fig.canvas.draw()
    lgbbox = lg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    lg.remove()
    if lgbbox.width > args.lgwidth:
      ncol -= 1      
      break

  lg = fig.legend(label_list,loc='outside lower center', bbox_to_anchor=(0.5,args.lgyoffset), ncol=ncol,fontsize=lgfsize) 
  fig.tight_layout()
  fig.savefig('history_'+args.runtag+'_dipole_axial.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Equatorial Dipole strength
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc, frameon=True)
  ax = plt.gca()
  xmn = np.min(time_list[0]*tfac)
  xmx = np.max(time_list[0]*tfac)

  for run in range(len(time_list)):
    LWu = LW-run*LWM
    if (args.lw):
      LWu=args.lw
    plt.plot(time_list[run]*tfac, eq_dipole_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
        markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

  plt.xlim(xmin=xmn,xmax=xmx)
  xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)
  plt.title('Equatorial Dipole Strength', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  ax.tick_params(axis='y',labelsize=fsize)
  ax.grid(zorder=0)

  renderer = fig.canvas.get_renderer()
  for ncol in range(1,LABEL_LEN+1):
    lg = fig.legend(label_list,loc='outside lower center', ncol=ncol,fontsize=lgfsize)
    fig.canvas.draw()
    lgbbox = lg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    lg.remove()
    if lgbbox.width > args.lgwidth:
      ncol -= 1         
      break

  lg = fig.legend(label_list,loc='outside lower center', bbox_to_anchor=(0.5,args.lgyoffset), ncol=ncol,fontsize=lgfsize) 
  fig.tight_layout()
  fig.savefig('history_'+args.runtag+'_dipole_eq.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
  
#
# ****** Validation
#
  if (args.valrun):
    fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
    ax = plt.gca()

    for run in range(len(time_list)):
      LWu = LW-run*LWM
      if (args.lw):
        LWu=args.lw
      plt.plot(time_list[run]*tfac, (1e5)*valerr_list[run],color=COLORS[run],linewidth=LWu,marker=MARKERS[run],
          markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=COLORS[run])

    plt.title('Validation Error', {'fontsize': fsize, 'color': tc})
    xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)
    plt.ylabel('(CV)RMSD ($10^{-5}$)', {'fontsize': fsize, 'color': tc})
    ax.tick_params(axis='y',labelsize=fsize)
    ax.grid(zorder=0)

    renderer = fig.canvas._get_renderer()
    for ncol in range(1,LABEL_LEN+1):
      lg = fig.legend(label_list,loc='outside lower center', ncol=ncol,fontsize=lgfsize)
      fig.canvas.draw()
      lgbbox = lg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
      lg.remove()
      if lgbbox.width > args.lgwidth:
        ncol -= 1     
        break

    lg = fig.legend(label_list,loc='outside lower center', bbox_to_anchor=(0.5,args.lgyoffset), ncol=ncol,fontsize=lgfsize) 
    fig.tight_layout()  
    fig.savefig('history_'+args.runtag+'_val.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  

#
# ****** x-axis label and ticks
#  


def get_xticks(args,xmn,xmx,init_locs):
  seconds = 1
  minutes = 60
  hours = 3600
  days = 86400
  weeks = 604800
  months = -100
  cr = 2356586
  years = 31556952
  default = hours
  
  xcUnitsSec = default
  if (args.xc_units):
    if args.xc_units == 'years':
      xcUnitsSec = years
    elif args.xc_units == 'cr':
      xcUnitsSec = cr
    elif args.xc_units == 'months':
      if (args.xunits != 'date'):
        raise Exception("Can only use 'months' with xunits of date")
    elif args.xc_units == 'weeks':
      xcUnitsSec = weeks
    elif args.xc_units == 'days':
      xcUnitsSec = days
    elif args.xc_units == 'hours':
      xcUnitsSec = hours
    elif args.xc_units == 'minutes':
      xcUnitsSec = minutes
    elif args.xc_units == 'seconds':
      xcUnitsSec = seconds
  else:
    if args.xunits == 'years':
      xcUnitsSec = years
    elif args.xunits == 'cr':
      xcUnitsSec = cr
    elif args.xunits == 'weeks':
      xcUnitsSec = weeks
    elif args.xunits == 'days':
      xcUnitsSec = days
    elif args.xunits == 'hours':
      xcUnitsSec = hours
    elif args.xunits == 'minutes':
      xcUnitsSec = minutes
    elif args.xunits == 'seconds':
      xcUnitsSec = seconds
      
  if (args.utstart):
    sDate = datetime.strptime(args.utstart, '%Y-%m-%dT%H:%M:%S')
    utstartSecs = int(sDate.replace(tzinfo=timezone.utc).timestamp())
  else:
    utstartSecs = 0
  initLocs_uttime = init_locs*hours+utstartSecs
  xmn_uttime = xmn*hours+utstartSecs
  xmx_uttime = xmx*hours+utstartSecs
  if args.xunits == "date":
    if (args.utstart):
      locs, labels = date_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime)
    else:
      raise Exception("Did not specify a utstart")
  elif args.xunits == "seconds":
    locs, labels = since_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime,seconds)
  elif args.xunits == "minutes":
    locs, labels = since_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime,minutes)
  elif args.xunits == "hours":
    locs, labels = since_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime,hours)
  elif args.xunits == "days":
    locs, labels = since_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime,days)
  elif args.xunits == "weeks":
    locs, labels = since_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime,weeks)
  elif args.xunits == "cr":
    if (args.utstart): 
      locs, labels = cr_xticks(args,xcUnitsSec,xmn_uttime,xmx_uttime)
    else:
      locs, labels = since_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime,cr)
  elif args.xunits == "years":
    locs, labels = since_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime,years)
  locs = (np.array(locs)-utstartSecs)/3600 
  labels = [label.replace('//','\n') for label in labels]
  return locs, labels, utstartSecs


def date_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime):
  locs = []
  labels = []
  if (args.xc_units == 'months'):
    xformat = '%m/%Y/%d/%H/%M/%S'
    if (args.xcadence):
      cadence = int(args.xcadence)
    else:
      tempArray = np.array(initLocs_uttime)
      cadence = int(np.average(np.diff(tempArray))/2678400)
      if cadence == 0:
        cadence = 1 
    currDate = datetime.utcfromtimestamp(xmn_uttime).strftime(xformat)
    locs.append(xmn_uttime)
    currDate=currDate.split('/')
    currDate[0]=str(int(currDate[0])+cadence)
    if int(currDate[0]) > 12:
      tempM = int(currDate[0])
      currDate[0]=str(tempM%12)
      currDate[1]=str(int(currDate[1])+int(tempM/12))
    currDate='/'.join(currDate)
    currDateSeconds = datetime.strptime(currDate,xformat).replace(tzinfo=timezone.utc).timestamp()
    while currDateSeconds <= xmx_uttime:
      locs.append(currDateSeconds)
      currDate=currDate.split('/')
      currDate[0]=str(int(currDate[0])+cadence)
      if int(currDate[0]) > 12:
        tempM = int(currDate[0])
        currDate[0]=str(tempM%12)
        currDate[1]=str(int(currDate[1])+int(tempM/12))
      currDate='/'.join(currDate)
      currDateSeconds = datetime.strptime(currDate,xformat).replace(tzinfo=timezone.utc).timestamp()
  elif (args.xunits == 'years' or args.xc_units == 'years'):
    xformat = '%Y/%m/%d/%H/%M/%S'
    if (args.xcadence):
      cadence = int(args.xcadence)
    else:
      tempArray = np.array(initLocs_uttime)
      cadence = int(np.average(np.diff(tempArray))/31556952)
      if cadence == 0:
        cadence = 1 
    currDate = datetime.utcfromtimestamp(xmn_uttime).strftime(xformat)
    locs.append(xmn_uttime)
    currDate=currDate.split('/')
    currDate[0]=str(int(currDate[0])+cadence)
    currDate='/'.join(currDate)
    currDateSeconds = datetime.strptime(currDate,xformat).replace(tzinfo=timezone.utc).timestamp()
    while currDateSeconds <= xmx_uttime:
      locs.append(currDateSeconds)
      currDate=currDate.split('/')
      currDate[0]=str(int(currDate[0])+cadence)
      currDate='/'.join(currDate)
      currDateSeconds = datetime.strptime(currDate,xformat).replace(tzinfo=timezone.utc).timestamp()
  else:
    if (args.xcadence):
      cadence = xcUnitsSec*int(args.xcadence)
    else:
      tempArray = np.array(initLocs_uttime)
      cadence = np.average(np.diff(tempArray)) 
    if (args.utstartxtick):
      currDate = datetime.strptime(args.utstartxtick,'%Y-%m-%dT%H:%M:%S').replace(tzinfo=timezone.utc).timestamp()
    else:
      currDate = xmn_uttime
    locs.append(currDate)
    skip = int(cadence)
    currDate = currDate + skip
    while currDate <= xmx_uttime:
      locs.append(currDate)
      currDate = currDate + skip
  for loc in locs:
    labels.append(datetime.utcfromtimestamp(loc).strftime(args.xformat))
  return locs, labels


def since_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime,secTimeUnit):
  locs = []
  labels = []
  if (args.xcadence):
    cadence = xcUnitsSec*int(args.xcadence)
  else:
    tempArray = np.array(initLocs_uttime)
    cadence = np.average(np.diff(tempArray)) 
  currDate = xmn_uttime
  locs.append(currDate)
  skip = int(cadence)
  currDate = currDate + skip
  while currDate <= xmx_uttime:
    locs.append(currDate)
    currDate = currDate + skip
  labelOffSet = int(locs[0]/secTimeUnit)
  for loc in locs:
      labels.append(str(int(loc/secTimeUnit)-labelOffSet))
  return locs, labels


def cr_xticks(args,xcUnitsSec,xmn_uttime,xmx_uttime):
  locs = []
  labels = []
  xformat = '%Y-%m-%dT%H:%M:%S'

  if (args.utstartxtick):
    currDateSec = datetime.strptime(args.utstartxtick,'%Y-%m-%dT%H:%M:%S').replace(tzinfo=timezone.utc).timestamp()
  else:
    currDateSec = xmn_uttime

  if (args.xcadence):
    if (args.xc_units):
      if (args.xc_units == 'cr'):
        endOffset = 0 
        if (args.xcrpos == 'end'):
          endOffset = 1
        cr_num = int(carrington_rotation_number(np.datetime64(int(currDateSec),'s')))
        if (args.xcrpos == 'center'):
          cr_num = cr_num-0.5
          currDate = int(carrington_rotation_time(cr_num).unix)
          if currDate < xmn_uttime:
            cr_num = cr_num+1
            currDate = int(carrington_rotation_time(cr_num).unix)
        else:
          currDate = int(carrington_rotation_time(cr_num).unix)
        if currDate<xmn_uttime:
          cr_num = cr_num+1
          currDate =  int(carrington_rotation_time(cr_num).unix)
        locs.append(currDate)
        labels.append(str(int(cr_num)-endOffset))
        cr_num = cr_num+int(args.xcadence)
        currDate = int(carrington_rotation_time(cr_num).unix)
        while currDate <= xmx_uttime:
          locs.append(currDate)
          labels.append(str(int(cr_num)-endOffset))
          cr_num = cr_num+int(args.xcadence)
          currDate = int(carrington_rotation_time(cr_num).unix)
      else:
        endOffset = 0 
        if (args.xcrpos == 'end'):
          endOffset = 1
        cr_num = int(carrington_rotation_number(np.datetime64(int(currDateSec),'s')))
        if (args.xcrpos == 'center'):
          cr_num = cr_num-0.5
          currDate = int(carrington_rotation_time(cr_num).unix)
          if currDate < xmn_uttime:
            cr_num = cr_num+1
            currDate = int(carrington_rotation_time(cr_num).unix)
        else:
          currDate = int(carrington_rotation_time(cr_num).unix)
        if currDate<xmn_uttime:
          cr_num = cr_num+1
          currDate =  int(carrington_rotation_time(cr_num).unix)
        locs.append(currDate)
        labels.append(str(int(cr_num)-endOffset))
        cr_num = float(carrington_rotation_number(np.datetime64(int(currDate+xcUnitsSec*int(args.xcadence)),'s')))
        currDate = int(carrington_rotation_time(cr_num).unix)
        while currDate <= xmx_uttime:
          locs.append(currDate)
          labels.append(str(int(cr_num)-endOffset))
          cr_num = float(carrington_rotation_number(np.datetime64(int(currDate+xcUnitsSec*int(args.xcadence)),'s')))
          currDate = int(carrington_rotation_time(cr_num).unix)
    else:
      endOffset = 0 
      if (args.xcrpos == 'end'):
        endOffset = 1
      cr_num = int(carrington_rotation_number(np.datetime64(int(currDateSec),'s')))
      if (args.xcrpos == 'center'):
        cr_num = cr_num-0.5
        currDate = int(carrington_rotation_time(cr_num).unix)
        if currDate < xmn_uttime:
          cr_num = cr_num+1
          currDate = int(carrington_rotation_time(cr_num).unix)
      else:
        currDate = int(carrington_rotation_time(cr_num).unix)
      if currDate<xmn_uttime:
        cr_num = cr_num+1
        currDate =  int(carrington_rotation_time(cr_num).unix)
      locs.append(currDate)
      labels.append(str(int(cr_num)-endOffset))
      cr_num = cr_num+int(args.xcadence)
      currDate = int(carrington_rotation_time(cr_num).unix)
      while currDate <= xmx_uttime:
        locs.append(currDate)
        labels.append(str(int(cr_num)-endOffset))
        cr_num = cr_num+int(args.xcadence)
        currDate = int(carrington_rotation_time(cr_num).unix)
  else:
    endOffset = 0 
    if (args.xcrpos == 'end'):
      endOffset = 1
    cr_num = int(carrington_rotation_number(np.datetime64(int(currDateSec),'s')))
    if (args.xcrpos == 'center'):
      cr_num = cr_num-0.5
      currDate = int(carrington_rotation_time(cr_num).unix)
      if currDate < xmn_uttime:
        cr_num = cr_num+1
        currDate =  int(carrington_rotation_time(cr_num).unix)
    else:
      currDate = int(carrington_rotation_time(cr_num).unix)
    if currDate<xmn_uttime:
      cr_num = cr_num+1
      currDate =  int(carrington_rotation_time(cr_num).unix)
    locs.append(currDate)
    labels.append(str(int(cr_num)-endOffset))
    cr_num = cr_num+1
    currDate = int(carrington_rotation_time(cr_num).unix)
    while currDate <= xmx_uttime:
      locs.append(currDate)
      labels.append(str(int(cr_num)-endOffset))
      cr_num = cr_num+1
      currDate = int(carrington_rotation_time(cr_num).unix)
  return locs, labels


def xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs):
  if (args.rmode):
    rmode = 'anchor'
  else:
    rmode = 'default'
  if args.xunits == "date":
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
  elif args.xunits == "seconds":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    elif (args.utstart): 
      plt.xlabel('Seconds since '+ datetime.utcfromtimestamp(utstartSecs).strftime('UT-%Y-%m-%dT%H:%M:%S'), {'fontsize': args.fsize, 'color': tc})
    else:
      plt.xlabel('Seconds', {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)   
  elif args.xunits == "minutes":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    elif (args.utstart): 
      plt.xlabel('Minutes since '+ datetime.utcfromtimestamp(utstartSecs).strftime('UT-%Y-%m-%dT%H:%M:%S'), {'fontsize': args.fsize, 'color': tc})
    else:
      plt.xlabel('Minutes', {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)  
  elif args.xunits == "hours":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    elif (args.utstart): 
      plt.xlabel('Hours since '+ datetime.utcfromtimestamp(utstartSecs).strftime('UT-%Y-%m-%dT%H:%M:%S'), {'fontsize': args.fsize, 'color': tc})
    else:
      plt.xlabel('Hours', {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)  
  elif args.xunits == "days":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    elif (args.utstart): 
      plt.xlabel('Days since '+ datetime.utcfromtimestamp(utstartSecs).strftime('UT-%Y-%m-%dT%H:%M:%S'), {'fontsize': args.fsize, 'color': tc})
    else:
      plt.xlabel('Days', {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)  
  elif args.xunits == "weeks":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    elif (args.utstart): 
      plt.xlabel('Weeks since '+ datetime.utcfromtimestamp(utstartSecs).strftime('UT-%Y-%m-%dT%H:%M:%S'), {'fontsize': args.fsize, 'color': tc})
    else:
      plt.xlabel('Weeks', {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)  
  elif args.xunits == "cr":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    elif (args.utstart): 
      if (args.xcrpos == 'start'):
        plt.xlabel('Carrington Rotation', {'fontsize': args.fsize, 'color': tc})
      else: 
        plt.xlabel('Carrington Rotation ('+args.xcrpos+')', {'fontsize': args.fsize, 'color': tc})
    else:
      if (args.xcrpos == 'start'):
        plt.xlabel('Carrington Rotation', {'fontsize': args.fsize, 'color': tc})
      else: 
        plt.xlabel('Carrington Rotation ('+args.xcrpos+')', {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize) 
  elif args.xunits == "years":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    elif (args.utstart): 
      plt.xlabel('Years since '+ datetime.utcfromtimestamp(utstartSecs).strftime('UT-%Y-%m-%dT%H:%M:%S'), {'fontsize': args.fsize, 'color': tc})
    else:
      plt.xlabel('Years', {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)  
  else:
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    else:
      plt.xlabel('Hours', {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)

def rearr(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()

