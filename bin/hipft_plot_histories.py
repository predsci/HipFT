#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.time import Time
from sunpy.coordinates.sun import carrington_rotation_time, carrington_rotation_number
import os
import itertools

# Version 1.8.3

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

  parser.add_argument('-rlist',
    help='A comma separated list of history file numbers to include',
    type=str,
    required=False,
    default=' ')

  parser.add_argument('-rexclude',
    help='A comma separated list of history file numbers to ignore (Overrides rlist)',
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
    help='Start Date of full data assimilation csv file (not start date of run) in UT: YYYY-MM-DDTHH:MM:SS',
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
    default=-0.025,
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
    default=1.0,
    required=False)

  parser.add_argument('-flw',
    help='Line width of fill area (default is -lw)',
    dest='flw',
    type=float,
    default=0.00001,
    required=False)

  parser.add_argument('-ms',
    help='Marker size',
    dest='ms',
    type=float,
    default=6.0,
    required=False)

  parser.add_argument('-lms',
    help='Lend marker size (default is -ms)',
    dest='lms',
    type=float,
    default=10.0,
    required=False)

  parser.add_argument('-summary',
    help='Summary mode that plots a mean line and upper-lower bound area.',
    dest='summary',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-tai',
    help='Use tai which ignores leap seconds.',
    dest='tai',
    action='store_true',
    default=False,
    required=False)

  return parser.parse_args()


def run(args):  

  flux_fac=1e-21

  arg_dict = vars(args)
  temphist_list = arg_dict['histfiles'].split(',')
  rexclude_list = arg_dict['rexclude'].split(',')
  rlist_list = arg_dict['rlist'].split(',')
  label_list = arg_dict['labels'].split(',')

  if rlist_list[0] == ' ':
    rlist_list='all'

  rList = []

  if arg_dict['histfiles'] == ' ':
    temphist_list=[]
    wDir=os.getcwd()
    for file in os.listdir(wDir):
      if "hipft_history_sol_r" in file and file.endswith(".out"):
        temphist_list.append(wDir+'/'+file)
    temphist_list = sorted(temphist_list)


  hist_list=[]

  for file in temphist_list:
    if arg_dict['histfiles'] == ' ':
      r=int((file.split('/')[-1]).replace('hipft_history_sol_r','').replace('.out',''))
      if str(r) in rexclude_list:
        continue
      elif str(r) in rlist_list or 'all' in rlist_list:
        hist_list.append(file)
        rList.append(r)
    else:
      hist_list.append(file)
      rList.append(1)

  NOTindividual=True
  if len(hist_list) == 1:
    NOTindividual=False 

  # Build the default list of labels based on the number runs the user entered
  # if no label list is entered
  if arg_dict['labels'] == ' ':
      label_list = []
      def_label = 'r'
      for i in rList:
        label_list.append(def_label+str(i))

  LABEL_LEN = len(label_list)
  # Validate the list arguments:
  if not len(hist_list) == LABEL_LEN:
      print('ERROR: Number of runs, dirs, and labels must match. Use -h for more information.')
      quit()

  if LABEL_LEN > 256:
      print('ERROR: Can only compare a maximum of 256 runs.')
      quit()

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
  LMS = args.lms

  LW = args.lw
  FLW = args.flw

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

    number_of_data_points = len(hist_sol_full)
    
    #thin out data in hist_sol
    if samples > 1:
      indices = np.linspace(0, number_of_data_points-1, samples, endpoint=True, dtype=int)
      hist_sol = hist_sol_full.iloc[indices]
    else:
      hist_sol = hist_sol_full

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
# ****** Create needed parameters and lists.
#  
  time_tfac = np.array(time_list,dtype=object)*tfac
  xmn = np.amin([np.amin(arr) for arr in time_tfac])
  xmx = np.amax([np.amax(arr) for arr in time_tfac])
  flux_tot_un_FF = np.array(flux_tot_un_list,dtype=object)*flux_fac
  fluxm_FF = np.abs(np.array(fluxm_list,dtype=object))*flux_fac
  fluxp_FF = np.array(fluxp_list,dtype=object)*flux_fac

  flux_tot_s_FF = np.array(flux_tot_s_list,dtype=object)*flux_fac
  fluxp_pn_FF = np.array(fluxp_pn_list,dtype=object)*flux_fac
  fluxm_pn_FF = np.array(fluxm_pn_list,dtype=object)*flux_fac
  fluxp_ps_FF = np.array(fluxp_ps_list,dtype=object)*flux_fac
  fluxm_ps_FF = np.array(fluxm_ps_list,dtype=object)*flux_fac

  valerr=np.array(valerr_list,dtype=object)*(1e5)

#
# ****** Total flux imbalance.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  if args.summary:
    summaryMode(flux_imb_list,time_tfac[0],LW,FLW,fsize,plt)
  else:
    normalMode(plt,ax,fig,args,flux_imb_list,time_tfac,COLORS,MARKERS,LW,MS,LMS,NOTindividual,LABEL_LEN,label_list,lgfsize,'k-')

  ymax = np.amax([np.amax(np.abs(arr)) for arr in flux_imb_list])
  ymin = -ymax 

  plt.title('Relative Flux Imbalance', {'fontsize': fsize, 'color': tc})
  plt.ylabel('%', {'fontsize': fsize, 'color': tc})
  init_locs = plt.xticks()[0]
  locs, labels, utstartSecs = get_xticks(args,xmn,xmx,init_locs)
  
  makeAxes(args,locs,labels,tc,ax,fig,plt,utstartSecs,xmn,xmx,ymin,ymax,fsize)
  fig.savefig('history_'+args.runtag+'_flux_imb_pm.png', bbox_inches='tight', dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Total (+) and (-) flux.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  if args.summary:
    legend1 = summaryMode2(fluxm_FF,fluxp_FF,time_tfac[0],LW,FLW,fsize,plt,"Blue","Red","|Flux (-)|","Flux (+)")
  else:
    legend1 = normalMode2(plt,ax,fig,args,fluxm_FF,fluxp_FF,time_tfac,COLORS,MARKERS,LW,MS,LMS,fsize,\
      NOTindividual,LABEL_LEN,label_list,lgfsize,'Blue','Red',["|Flux (-)|","Flux (+)"])

  ymin=0.0#np.amin([np.amin(fluxm_FF),np.amin(fluxp_FF)])
  ymax = max(np.amax(np.abs(arr)) for arr in fluxm_FF + fluxp_FF)
  plt.title('Total Positive and Negative Flux', {'fontsize': fsize, 'color': tc})
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})  
  ax.add_artist(legend1)  

  makeAxes(args,locs,labels,tc,ax,fig,plt,utstartSecs,xmn,xmx,ymin,ymax,fsize)
  plt.ylim(ymin=0.0)
  
  fig.savefig('history_'+args.runtag+'_flux_pm.png', bbox_inches="tight", dpi=args.dpi, \
              facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Total unsigned flux.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  if args.summary:
    summaryMode(flux_tot_un_FF,time_tfac[0],LW,FLW,fsize,plt)
  else:
    normalMode(plt,ax,fig,args,flux_tot_un_FF,time_tfac,COLORS,MARKERS,LW,MS,LMS,NOTindividual,LABEL_LEN,label_list,lgfsize,'k-')
  
  plt.title('Total Unsigned Flux', {'fontsize': fsize, 'color': tc})
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})

  ymin=0.0 #np.amin(flux_tot_un_FF)
  ymax = np.amax([np.amax(arr) for arr in flux_tot_un_FF])
  makeAxes(args,locs,labels,tc,ax,fig,plt,utstartSecs,xmn,xmx,ymin,ymax,fsize)
  plt.ylim(ymin=0.0)
  fig.savefig('history_'+args.runtag+'_flux_total_unsigned.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Total signed flux.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  if args.summary:
    summaryMode(flux_tot_s_FF,time_tfac[0],LW,FLW,fsize,plt)
  else:
    normalMode(plt,ax,fig,args,flux_tot_s_FF,time_tfac,COLORS,MARKERS,LW,MS,LMS,NOTindividual,LABEL_LEN,label_list,lgfsize,'b-')

  plt.title('Total Signed Flux', {'fontsize': fsize, 'color': tc})
  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})

  #ymin=np.amin(flux_tot_s_FF)
  #ymax=np.amax(flux_tot_s_FF)
  ymax = np.amax([np.amax(np.abs(arr)) for arr in flux_tot_s_FF])
  ymin = -ymax 
  
  makeAxes(args,locs,labels,tc,ax,fig,plt,utstartSecs,xmn,xmx,ymin,ymax,fsize)
  fig.savefig('history_'+args.runtag+'_flux_total_signed.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')  
#
# ****** Polar +/- fluxes.
#
  ymax = max(np.amax(np.abs(arr)) for arr in fluxp_pn_FF + fluxm_pn_FF + fluxp_ps_FF + fluxm_ps_FF)
  ymin = -ymax 
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  if args.summary:
    legend1 = summaryMode4(fluxp_pn_FF,fluxm_pn_FF,fluxp_ps_FF,fluxm_ps_FF,time_tfac[0],\
      LW,FLW,fsize,plt,"Red","Blue","firebrick","navy","N (+)","N (-)","S (+)","S (-)")
  else:
    legend1 = normalMode4(plt,ax,fig,args,fluxp_pn_FF,fluxm_pn_FF,fluxp_ps_FF,fluxm_ps_FF,time_tfac,COLORS,MARKERS,LW,MS,LMS,\
        fsize,NOTindividual,LABEL_LEN,label_list,lgfsize,"Red","Blue","firebrick","navy",["N (+)","N (-)","S (+)","S (-)"])

  plt.ylabel('$10^{21}$ Mx', {'fontsize': fsize, 'color': tc})
  plt.title('Polar Flux (within 30 degrees of poles)', {'fontsize': fsize, 'color': tc})

  ax.add_artist(legend1)  
  makeAxes(args,locs,labels,tc,ax,fig,plt,utstartSecs,xmn,xmx,ymin,ymax,fsize)
  fig.savefig('history_'+args.runtag+'_flux_poles_30.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Polar average field strengths.
#
  ymax = max(np.amax(np.abs(arr)) for arr in np.array(pole_n_avg_field_list,dtype=object) + np.array(pole_s_avg_field_list,dtype=object))
  ymin = -ymax 
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  if args.summary:
    legend1 = summaryMode2(pole_n_avg_field_list,pole_s_avg_field_list,time_tfac[0],LW,FLW,fsize,plt,"Black","Blue","North","South")
  else:
    legend1 = normalMode2(plt,ax,fig,args,pole_n_avg_field_list,pole_s_avg_field_list,time_tfac,COLORS,MARKERS,LW,MS,LMS,fsize,\
      NOTindividual,LABEL_LEN,label_list,lgfsize,'Black','Blue',["North","South"])

  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  plt.title('Polar Average Field (within 30 degrees of poles)', {'fontsize': fsize, 'color': tc})
  
  ax.add_artist(legend1)  
  makeAxes(args,locs,labels,tc,ax,fig,plt,utstartSecs,xmn,xmx,ymin,ymax,fsize)
  fig.savefig('history_'+args.runtag+'_field_poles_30.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Min and max field.
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  if args.summary:
    legend1 = summaryMode2(brmax_list,np.abs(np.array(brmin_list,dtype=object)),time_tfac[0],LW,FLW,fsize,plt,"blue","red","max(Br)","|min(Br)|")
  else:
    legend1 = normalMode2(plt,ax,fig,args,brmax_list,np.abs(np.array(brmin_list,dtype=object)),time_tfac,COLORS,MARKERS,LW,MS,LMS,fsize,\
      NOTindividual,LABEL_LEN,label_list,lgfsize,'blue','red',["max(Br)","|min(Br)|"])
  
  ymax = max(np.amax(np.abs(arr)) for arr in brmax_list + brmin_list)
  ymin = 0.0 
  
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})
  plt.title('Min and Max Br', {'fontsize': fsize, 'color': tc})

  ax.add_artist(legend1)  
  makeAxes(args,locs,labels,tc,ax,fig,plt,utstartSecs,xmn,xmx,ymin,ymax,fsize)

  if not args.summary:
      plt.ylim(ymin=0.0)
      
  fig.savefig('history_'+args.runtag+'_br.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Axial Dipole strength
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
  ax = plt.gca()

  if args.summary:
    summaryMode(ax_dipole_list,time_tfac[0],LW,FLW,fsize,plt)
  else:
    normalMode(plt,ax,fig,args,ax_dipole_list,time_tfac,COLORS,MARKERS,LW,MS,LMS,NOTindividual,LABEL_LEN,label_list,lgfsize,'k-')

  plt.title('Axial Dipole Strength', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})

  #ymin=np.amin(ax_dipole_list)
  #ymax=np.amax(ax_dipole_list)
  
  ymax = np.amax([np.amax(np.abs(arr)) for arr in ax_dipole_list])
  ymin = -ymax 
  
  
  makeAxes(args,locs,labels,tc,ax,fig,plt,utstartSecs,xmn,xmx,ymin,ymax,fsize)
  fig.savefig('history_'+args.runtag+'_dipole_axial.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
#
# ****** Equatorial Dipole strength
#
  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc, frameon=True)
  ax = plt.gca()

  if args.summary:
    summaryMode(eq_dipole_list,time_tfac[0],LW,FLW,fsize,plt)
  else:
    normalMode(plt,ax,fig,args,eq_dipole_list,time_tfac,COLORS,MARKERS,LW,MS,LMS,NOTindividual,LABEL_LEN,label_list,lgfsize,'k-')
 
  ymin = np.amin([np.amin(arr) for arr in eq_dipole_list])
  ymax = np.amax([np.amax(arr) for arr in eq_dipole_list])
  plt.title('Equatorial Dipole Strength', {'fontsize': fsize, 'color': tc})
  plt.ylabel('Gauss', {'fontsize': fsize, 'color': tc})

  makeAxes(args,locs,labels,tc,ax,fig,plt,utstartSecs,xmn,xmx,ymin,ymax,fsize)
  fig.savefig('history_'+args.runtag+'_dipole_eq.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  plt.close('all')
  
#
# ****** Validation
#
  if (args.valrun):
    fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=True)
    ax = plt.gca()

    if args.summary:
      summaryMode(valerr,time_tfac[0],LW,FLW,fsize,plt)
    else:
      normalMode(plt,ax,fig,args,valerr,time_tfac,COLORS,MARKERS,LW,MS,LMS,NOTindividual,LABEL_LEN,label_list,lgfsize,'k-')

    plt.title('Validation Error', {'fontsize': fsize, 'color': tc})
    xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)
    plt.ylabel('(CV)RMSD ($10^{-5}$)', {'fontsize': fsize, 'color': tc})
    ax.tick_params(axis='y',labelsize=fsize)
    ax.grid(zorder=0)

    fig.tight_layout()  
    fig.savefig('history_'+args.runtag+'_val.png', bbox_inches="tight", dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)
  
#
# ****** Helper Functions
#


def makeAxes(args,locs,labels,tc,ax,fig,plt,utstartSecs,xmn,xmx,ymin,ymax,fsize):
  plt.xlim(xmin=xmn,xmax=xmx)
  deltay = np.abs(ymax-ymin)
  ybuf = 0.175*deltay
  plt.ylim(ymin=ymin-ybuf,ymax=ymax+ybuf)
  xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)
  ax.tick_params(axis='y',labelsize=fsize)
  ax.grid(zorder=0)
  fig.tight_layout()


def summaryMode(llist,xlist,LW,FLW,fsize,plt):
  summaryHelper1(llist,xlist,LW,FLW,plt,"k",'blue',"$\mu$","$\sigma$","[Min,Max]")
  plt.legend(loc='upper left',fontsize=fsize, ncol=3)


def summaryMode2(llist1,llist2,xlist,LW,FLW,fsize,plt,CM1,CM2,LT1,LT2):
  summaryHelper1(llist1,xlist,LW,FLW,plt,CM1,CM1,LT1,'_nolegend_','_nolegend_')
  summaryHelper1(llist2,xlist,LW,FLW,plt,CM2,CM2,LT2,'_nolegend_','_nolegend_')
  legend1 = plt.legend(loc='upper right',fontsize=fsize, ncol=2)
  addLegendGeneric(LW,FLW,fsize,plt,"k","gray","$\mu$","$\sigma$","[Min,Max]")
  return legend1


def summaryMode4(llist1,llist2,llist3,llist4,xlist,LW,FLW,fsize,plt,CM1,CM2,CM3,CM4,LT1,LT2,LT3,LT4):
  summaryHelper1(llist1,xlist,LW,FLW,plt,CM1,CM1,LT1,'_nolegend_','_nolegend_')
  summaryHelper1(llist2,xlist,LW,FLW,plt,CM2,CM2,LT2,'_nolegend_','_nolegend_')
  summaryHelper1(llist3,xlist,LW,FLW,plt,CM3,CM3,LT3,'_nolegend_','_nolegend_')
  summaryHelper1(llist4,xlist,LW,FLW,plt,CM4,CM4,LT4,'_nolegend_','_nolegend_')
  legend1 = plt.legend(loc='upper right',fontsize=fsize, ncol=4)
  addLegendGeneric(LW,FLW,fsize,plt,"k","gray","$\mu$","$\sigma$","[Min,Max]")
  return legend1


def summaryHelper1(llist,xlist,LW,FLW,plt,CM,CF,LT1,LT2,LT3):
  upperBound=np.max(llist,axis=0)
  lowerBound=np.min(llist,axis=0)
  middleVal=np.mean(llist,axis=0)
  stdVal=np.std(llist,axis=0)
  stdVal0=middleVal - stdVal
  stdVal1=middleVal + stdVal
  plt.plot(xlist, middleVal,color=CM,linewidth=LW,label=LT1)
  plt.fill_between(xlist,stdVal0,stdVal1,linewidth=FLW,alpha=0.5,color=CF,label=LT2)
  plt.fill_between(xlist,lowerBound,upperBound,linewidth=FLW,alpha=0.25,color=CF,label=LT3)


def addLegendGeneric(LW,FLW,fsize,plt,CM,CF,LT1,LT2,LT3):
  h1=plt.plot(np.NaN, np.NaN,color=CM,linewidth=LW,label='_nolegend_')
  h2=plt.fill_between([np.NaN],[0],[0],linewidth=FLW,alpha=0.5,color=CF,label='_nolegend_')
  h3=plt.fill_between([np.NaN],[0],[0],linewidth=FLW,alpha=0.25,color=CF,label='_nolegend_')
  plt.legend([h1[0],h2,h3],[LT1,LT2,LT3],loc='lower left',fontsize=fsize, ncol=3)


def normalMode(plt,ax,fig,args,llist,xlist,COLORS,MARKERS,LW,MS,LMS,NOTindividual,LABEL_LEN,label_list,lgfsize,CL):
  if NOTindividual:
    for run in range(len(xlist)):
      plothelper1(xlist[run],llist[run],COLORS[run],MARKERS[run],LW,MS)
    addLegend(args,LABEL_LEN,LMS,label_list,lgfsize,ax,fig)
  else:
    plt.plot(xlist[0],llist[0],CL,linewidth=LW,markersize=MS)


def normalMode2(plt,ax,fig,args,llist1,llist2,xlist,COLORS,MARKERS,LW,MS,LMS,fsize,NOTindividual,\
                LABEL_LEN,label_list,lgfsize,CL1,CL2,LegT):
  if NOTindividual:
    firstRun=True
    for run in range(len(xlist)):
      if firstRun:
        h1=plothelper2(xlist[run],llist1[run],CL1,COLORS[run],MARKERS[run],LW,MS)
        h2=plothelper4(xlist[run],llist2[run],CL2,COLORS[run],MARKERS[run],LW,MS)
        firstRun=False
      else:
        plothelper2(xlist[run],llist1[run],CL1,COLORS[run],MARKERS[run],LW,MS)
        plothelper3(xlist[run],llist2[run],CL2,COLORS[run],MARKERS[run],LW,MS)
    legend1 = plt.legend([h1[0],h2[0]],LegT,loc='upper right',fontsize=fsize, ncol=2)
    addLegend(args,LABEL_LEN,LMS,label_list,lgfsize,ax,fig)
    return legend1
  else:
    h1 = plt.plot(xlist[0],llist1[0],color=CL1,linewidth=LW,markersize=MS)
    h2 = plt.plot(xlist[0],np.abs(llist2[0]),color=CL2,linewidth=LW,markersize=MS)
    return plt.legend([h1[0],h2[0]],LegT,loc='upper right',fontsize=fsize, ncol=2)

def normalMode4(plt,ax,fig,args,llist1,llist2,llist3,llist4,xlist,COLORS,MARKERS,LW,MS,LMS,fsize,NOTindividual,\
                LABEL_LEN,label_list,lgfsize,CL1,CL2,CL3,CL4,LegT):
  if NOTindividual:
    firstRun=True
    for run in range(len(xlist)):
      if firstRun:
        h1=plothelper2(xlist[run],llist1[run],CL1,COLORS[run],MARKERS[run],LW,MS)
        h2=plothelper4(xlist[run],llist2[run],CL2,COLORS[run],MARKERS[run],LW,MS)
        h3=plothelper4(xlist[run],llist3[run],CL3,COLORS[run],MARKERS[run],LW,MS)
        h4=plothelper4(xlist[run],llist4[run],CL4,COLORS[run],MARKERS[run],LW,MS)
        firstRun=False
      else:
        plothelper2(xlist[run],llist1[run],CL1,COLORS[run],MARKERS[run],LW,MS)
        plothelper3(xlist[run],llist2[run],CL2,COLORS[run],MARKERS[run],LW,MS)
        plothelper3(xlist[run],llist3[run],CL3,COLORS[run],MARKERS[run],LW,MS)
        plothelper3(xlist[run],llist4[run],CL4,COLORS[run],MARKERS[run],LW,MS)
    legend1 = plt.legend([h1[0],h2[0],h3[0],h4[0]],LegT,loc='upper right',fontsize=fsize, ncol=4)
    addLegend(args,LABEL_LEN,LMS,label_list,lgfsize,ax,fig)
    return legend1
  else:
    h1 = plt.plot(xlist[0],llist1[0],color=CL1,linewidth=LW,markersize=MS)
    h2 = plt.plot(xlist[0],llist2[0],color=CL2,linewidth=LW,markersize=MS) 
    h3 = plt.plot(xlist[0],llist3[0],color=CL3,linewidth=LW,markersize=MS)
    h4 = plt.plot(xlist[0],llist4[0],color=CL4,linewidth=LW,markersize=MS)
    return plt.legend([h1[0],h2[0],h3[0],h4[0]],LegT,loc='upper right',fontsize=fsize, ncol=4)


def plothelper1(x,l,c,m,LW,MS):
  plt.plot(x,l,color=c,linewidth=LW,marker=m,markersize=MS,markeredgewidth=0.0,fillstyle='full',markeredgecolor=c)

def plothelper2(x,l,c1,c2,m,LW,MS):
  return plt.plot(x,l,color=c1,linewidth=LW,marker=m,markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=c2,markeredgecolor=c1)

def plothelper3(x,l,c1,c2,m,LW,MS):
  plt.plot(x,l,color=c1,linewidth=LW,marker=m,markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=c2,markeredgecolor=c1,label='_nolegend_')

def plothelper4(x,l,c1,c2,m,LW,MS):
  return plt.plot(x,l,color=c1,linewidth=LW,marker=m,markersize=MS,markeredgewidth=0.0,fillstyle='full',markerfacecolor=c2,markeredgecolor=c1,label='_nolegend_')


def addLegend(args,LABEL_LEN,LMS,label_list,lgfsize,ax,fig):
  renderer = fig.canvas.get_renderer()
  for ncol in range(1,LABEL_LEN+1):
    lg = fig.legend(label_list,loc='outside lower center', ncol=ncol,fontsize=lgfsize)
    fig.canvas.draw()
    lgbbox = lg.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    lg.remove()
    if lgbbox.width > args.lgwidth:
      ncol -= 1
      break
  lg = fig.legend(label_list,ncol=ncol,fontsize=lgfsize) 
  lgh = args.lgyoffset - lg.get_window_extent().height/(7*args.dpi)
  lg.remove()
  lg = fig.legend(label_list,loc='outside lower center', bbox_to_anchor=(0.5,lgh), ncol=ncol,fontsize=lgfsize) 
  
  for handle in lg.legend_handles:
    handle.set_markersize(LMS)
    handle.set_linestyle("")

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
    utstartSecs = getSec(args.tai,args.utstart,'%Y-%m-%dT%H:%M:%S')
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
    t = getTimeObj(args.tai,xmn_uttime)
    currDate = t.strftime(xformat)
    locs.append(xmn_uttime)
    currDate=currDate.split('/')
    currDate[0]=str(int(currDate[0])+cadence)
    if int(currDate[0]) > 12:
      tempM = int(currDate[0])
      currDate[0]=str(tempM%12)
      currDate[1]=str(int(currDate[1])+int(tempM/12))
    currDate='/'.join(currDate)
    currDateSeconds = getSec(args.tai,currDate,xformat)
    while currDateSeconds <= xmx_uttime:
      locs.append(currDateSeconds)
      currDate=currDate.split('/')
      currDate[0]=str(int(currDate[0])+cadence)
      if int(currDate[0]) > 12:
        tempM = int(currDate[0])
        currDate[0]=str(tempM%12)
        currDate[1]=str(int(currDate[1])+int(tempM/12))
      currDate='/'.join(currDate)
      currDateSeconds = getSec(args.tai,currDate,xformat)
  elif (args.xunits == 'years' or args.xc_units == 'years'):
    xformat = '%Y/%m/%d/%H/%M/%S'
    if (args.xcadence):
      cadence = int(args.xcadence)
    else:
      tempArray = np.array(initLocs_uttime)
      cadence = int(np.average(np.diff(tempArray))/31556952)
      if cadence == 0:
        cadence = 1 
    t = getTimeObj(args.tai,xmn_uttime)
    currDate = t.strftime(xformat)
    locs.append(xmn_uttime)
    currDate=currDate.split('/')
    currDate[0]=str(int(currDate[0])+cadence)
    currDate='/'.join(currDate)
    currDateSeconds = getSec(args.tai,currDate,xformat)
    while currDateSeconds <= xmx_uttime:
      locs.append(currDateSeconds)
      currDate=currDate.split('/')
      currDate[0]=str(int(currDate[0])+cadence)
      currDate='/'.join(currDate)
      currDateSeconds = getSec(args.tai,currDate,xformat)
  else:
    if (args.xcadence):
      cadence = xcUnitsSec*int(args.xcadence)
    else:
      tempArray = np.array(initLocs_uttime)
      cadence = np.average(np.diff(tempArray)) 
    if (args.utstartxtick):
      currDate = getSec(args.tai,args.utstartxtick,'%Y-%m-%dT%H:%M:%S')
    else:
      currDate = xmn_uttime
    locs.append(currDate)
    skip = int(cadence)
    currDate = currDate + skip
    while currDate <= xmx_uttime:
      locs.append(currDate)
      currDate = currDate + skip
  if args.tai:
    for loc in locs:
      labels.append(Time(loc, format='unix_tai').strftime(args.xformat))
  else:
    for loc in locs:
      labels.append(Time(loc, format='unix').strftime(args.xformat))
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
    currDateSeconds = getSec(args.tai,args.utstartxtick,xformat)
  else:
    currDateSec = xmn_uttime

  if (args.xcadence):
    if (args.xc_units):
      if (args.xc_units == 'cr'):
        endOffset = 0 
        if (args.xcrpos == 'end'):
          endOffset = 1
        if args.tai:
          cr_num = int(carrington_rotation_number(Time(currDateSec, format='unix_tai')))
        else:
          cr_num = int(carrington_rotation_number(Time(currDateSec, format='unix')))
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
        if args.tai:
          cr_num = int(carrington_rotation_number(Time(currDateSec, format='unix_tai')))
        else:
          cr_num = int(carrington_rotation_number(Time(currDateSec, format='unix')))
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
        if args.tai:
          cr_num = float(carrington_rotation_number(Time(currDate+xcUnitsSec*int(args.xcadence), format='unix_tai')))
          currDate = int(carrington_rotation_time(cr_num).unix_tai)
        else:
          cr_num = float(carrington_rotation_number(Time(currDate+xcUnitsSec*int(args.xcadence), format='unix')))
          currDate = int(carrington_rotation_time(cr_num).unix)
        while currDate <= xmx_uttime:
          locs.append(currDate)
          labels.append(str(int(cr_num)-endOffset))
          if args.tai:
            cr_num = float(carrington_rotation_number(Time(currDate+xcUnitsSec*int(args.xcadence), format='unix_tai')))
            currDate = int(carrington_rotation_time(cr_num).unix_tai)
          else:
            cr_num = float(carrington_rotation_number(Time(currDate+xcUnitsSec*int(args.xcadence), format='unix')))
            currDate = int(carrington_rotation_time(cr_num).unix)
    else:
      endOffset = 0 
      if (args.xcrpos == 'end'):
        endOffset = 1
      if args.tai:
        cr_num = int(carrington_rotation_number(Time(currDateSec, format='unix_tai')))
      else:
        cr_num = int(carrington_rotation_number(Time(currDateSec, format='unix')))
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
    if args.tai:
      cr_num = int(carrington_rotation_number(Time(currDateSec, format='unix_tai')))
    else:
      cr_num = int(carrington_rotation_number(Time(currDateSec, format='unix')))
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
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Seconds since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
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
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Minutes since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
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
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Hours since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
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
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Days since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
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
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Weeks since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
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
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Years since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
    else:
      plt.xlabel('Years', {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)  
  else:
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    else:
      plt.xlabel('Hours', {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)


def gettLabel(atai,x):
  if atai:
    t = Time(x,format='unix_tai')
  else:
    t = Time(x,format='unix')
  return t.strftime('UT-%Y-%m-%dT%H:%M:%S')


def getTimeObj(atai,x):
    if atai:
      t = Time(x,format='unix_tai')
    else:  
      t = Time(x,format='unix')
    return t


def getSec(atai,x,xformat):
  if atai:
    t = Time.strptime(x,xformat,scale='tai')
    currDateSeconds = t.unix_tai
  else:
    t = Time.strptime(x,xformat,scale='utc')
    currDateSeconds = t.unix
  return currDateSeconds


def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()
