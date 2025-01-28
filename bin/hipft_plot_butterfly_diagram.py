#!/usr/bin/env python3
import argparse
import signal
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from sunpy.coordinates.sun import carrington_rotation_time, carrington_rotation_number
import os
import psihdf as ps
import psipals
import psimath
import multiprocessing as mp

# Version 1.3.1

def signal_handler(signal, frame):
  print('You pressed Ctrl+C! Stopping!')
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)


def argParsing():
  parser = argparse.ArgumentParser(description='plot_butterfly_diagram: This tool saves a plot of a 2D h5 butterfly diagram.')

  parser.add_argument(
    'iFile',
    help='Name of 2D hdf/h5 file')

  parser.add_argument('-o',
    help='Name of output png',
    dest='oFile',
    required=False)

  parser.add_argument('-cmin',
    help='Colormap Minimum',
    dest='cmin',
    type=float,
    required=False)

  parser.add_argument('-cmax',
    help='Colormap Maximum',
    dest='cmax',
    type=float,
    required=False)

  parser.add_argument('-csym',
    help='Make colormap center at 0 if there is negative data.',
    dest='csym',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-unit_label',
    help='Label for colorbar (units)',
    dest='unit_label',
    default='Gauss',
    required=False)

  parser.add_argument('-unit_fac',
    help='Unit factor to multiply data',
    dest='unit_fac',
    default='1.0',
    type=float,
    required=False)

  parser.add_argument('-title',
    help='Title to show on plot',
    dest='title',
    default=' ',
    required=False)

  parser.add_argument('-k',
    help='Plot with black background',
    dest='k',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-nogrid',
    help='Do not plot grid lines',
    dest='nogrid',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-finegrid',
    help='Make the grid more fine.',
    dest='finegrid',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-noax',
    help='Turn off axis and display image flush with border.',
    dest='noax',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-noscales',
    help='Turn off axis scales.',
    dest='noscales',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-nocb',
    help='Turn off colorbar.',
    dest='nocb',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-cmap',
    help='String of colormap name.  Can be any PSI colormap or matplotlib colormap. Run plot_psi_colormaps to view all PSI-specific maps.',
    dest='cmap',
    default='psi_blue_red',
    required=False)

  parser.add_argument('-smooth',
    help='Use smooth shading for plot.',
    dest='smooth',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-gridlines',
    help='Plot gridlines.',
    dest='gridlines',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-uniform',
    help='Ignore scales and plot uniform',
    dest='uscales',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-xmin',
    help='Set xmin.',
    dest='xmin',
    type=float,
    required=False)

  parser.add_argument('-xmax',
    help='Set xmax.',
    dest='xmax',
    type=float,
    required=False)

  parser.add_argument('-ymin',
    help='Set ymin.',
    dest='ymin',
    type=float,
    required=False)

  parser.add_argument('-ymax',
    help='Set ymax.',
    dest='ymax',
    type=float,
    required=False)

  parser.add_argument('-dpi',
    help='DPI for the resulting image.',
    dest='dpi',
    type=int,
    default=200,
    required=False)

  parser.add_argument('-x_cbar',
    help='Place the colorbar on the x-axis',
    dest='x_cbar',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-bgtrans',
    help='Make background outside of plot transparent.',
    dest='bgtrans',
    action='store_false',
    default=True,
    required=False)

  parser.add_argument('-v',
    help='Print details as script runs.',
    dest='verbose',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-c',
    help='Cadence of how many pixels to skip.',
    dest='cadence',
    default=0,
    required=False)

  parser.add_argument('-ignore_data_uttime',
    action='store_true',
    help='Ignore the UT time data',
    dest='ignore_data_uttime',
    required=False)
    
  parser.add_argument('-utstartxtick',
    help='UT date of first xtick: YYYY-MM-DDTHH:MM:SS',
    dest='utstartxtick',
    required=False)
    
  parser.add_argument('-xunits',
    help='Units of the x-axis (date, seconds, minutes, hours, days, weeks, cr, or years).',
    dest='xunits',
    required=False)

  parser.add_argument('-xcrpos',
    help='Carrington rotation marker location (start, center, end).',
    dest='xcrpos',
    default='start',
    required=False)

  parser.add_argument('-xformat',
    help="Format for the date option where // is treated as a newline (example: %%H:%%M:%%S//%%Y/%%m/%%d).",
    dest='xformat',
    default="%H:%M:%S//%Y/%m/%d",
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
    default=18,
    required=False)

  parser.add_argument('-fsize_xticklabels',
    help='Font size for xtick labels',
    dest='xlabelfsize',
    default=18,
    required=False)

  parser.add_argument('-3d',
    action='store_true',
    help='Flag to denote 3d h5 file.',
    dest='dim3',
    required=False)

  parser.add_argument('-slice',
    help='Slice to use for averaging if 3d.',
    dest='slice', 
    type=int,
    default=1,
    required=False)

  parser.add_argument('-slices',
    help='Slices to use for averaging if 3d in format: 1 2 3',
    dest='slices', 
    type=int, 
    nargs='+',
    required=False)

  parser.add_argument('-sall',
    action='store_true',
    help='Do all slices',
    dest='sall', 
    required=False)

  parser.add_argument('-xlabel',
    help='Label for x axis',
    dest='xlabel',
    required=False)

  parser.add_argument('-tai',
    help='Use tai which ignores leap seconds.',
    dest='tai',
    action='store_true',
    default=False,
    required=False)

  parser.add_argument('-np',
    help='Number of threads to use for movies.',
    dest='np',
    type=int,
    required=False)
  
  return parser.parse_args()


def run(args):

  if not args.np:
    args.np = int(os.getenv('OMP_NUM_THREADS', 1))

  #Check that file exists.
  if not os.path.exists(args.iFile):
    print('ERROR!  Butterfly h5 file not found:  '+args.iFile)
    exit(1)

  # Load colormaps:
  psipals.load()

  if (args.oFile is None):
    if (str(args.iFile).endswith('h5')):
      n = 3
    else:
      n = 4
    oFile = args.iFile[0:len(args.iFile) - n] + ".png"
  else:
    oFile = args.oFile
    if '.png' not in oFile:
      oFile = oFile + '.png'

  if bool(args.dim3):
    if (args.sall):
      xvec, yvec, zvec, data_in = ps.rdhdf_3d(args.iFile)
      with mp.Pool(processes=args.np) as pool:
            pool.starmap(process_file, [(args, islice, oFile, xvec, yvec, data_in, zvec[islice]) for islice in range(len(zvec))])
    elif (args.slices):
      xvec, yvec, zvec, data_in = ps.rdhdf_3d(args.iFile)
      with mp.Pool(processes=args.np) as pool:
            pool.starmap(process_file, [(args, islice-1, oFile, xvec, yvec, data_in, zvec[islice-1]) for islice in args.slices])
    else:
      xvec, yvec, zvec, data_in = ps.rdhdf_3d(args.iFile)
      data = np.squeeze(data_in[int(args.slice)-1,:,:])
      plot(args, xvec, yvec, data, oFile)
      matplotlib.pyplot.close()
  else:
    # Read the data:
    xvec, yvec, data = ps.rdhdf_2d(args.iFile)
    plot(args, xvec, yvec, data, oFile)
    matplotlib.pyplot.close()


def process_file(args, islice, oFile, xvec, yvec, data_in, zvec_i):
  oFileNew = oFile.replace('.png','_r'+str(int(zvec_i)).zfill(6)+'.png')
  data = np.squeeze(data_in[islice,:,:])
  plot(args, xvec, yvec, data, oFileNew)


def plot(args, xvec, yvec, data, oFile):

  #thin out data
  skip = int(args.cadence)+1
  xvec = xvec[::skip]
  data = data[:,::skip]

  if (len(xvec) == 0):
    xvec = np.array(range(0, len(data[0, :])))
  if (len(yvec) == 0):
    yvec = np.array(range(0, len(data[:, 0])))

  yvec = np.array([psimath.Math.t_rad_lat_deg(y_pt) for y_pt in yvec])

  # Set up data units:
  data = np.multiply(data, float(args.unit_fac))
  cbstr = args.unit_label

  # Set up uniform grid scales if requested:
  if (args.uscales):
    xvec_plot = np.array(range(0, len(xvec)))/3600
    yvec_plot = np.array(range(0, len(yvec)))
  else:
    xvec_plot = xvec/3600
    yvec_plot = yvec

  # Modify scales to be on a half-mesh of size N+1
  # for proper plotting with pcolormesh.

  # First get original limits (for plotting axis):
  xmin = np.min(xvec_plot)
  xmax = np.max(xvec_plot)
  ymin = np.min(yvec_plot)
  ymax = np.max(yvec_plot)

  if args.xmin is not None:
    xmin = float(args.xmin)
  if args.xmax is not None:
    xmax = float(args.xmax)
  if args.ymin is not None:
    ymin = float(args.ymin)
  if args.ymax is not None:
    ymax = float(args.ymax)

  if (not args.smooth):
    xvec_plot2 = np.zeros(len(xvec_plot) + 1)
    xvec_plot2[1:-1] = (xvec_plot[0:-1] + xvec_plot[1:]) / 2.0
    xvec_plot2[0] = xvec_plot2[1] - (xvec_plot[1] - xvec_plot[0])
    xvec_plot2[-1] = xvec_plot2[-2] + (xvec_plot[-1] - xvec_plot[-2])

    yvec_plot2 = np.zeros(len(yvec_plot) + 1)
    yvec_plot2[1:-1] = (yvec_plot[0:-1] + yvec_plot[1:]) / 2.0
    yvec_plot2[0] = yvec_plot2[1] - (yvec_plot[1] - yvec_plot[0])
    yvec_plot2[-1] = yvec_plot2[-2] + (yvec_plot[-1] - yvec_plot[-2])

    xvec_plot = xvec_plot2
    yvec_plot = yvec_plot2

  # Set up colormap scales:
  cmin = np.min(data)
  cmax = np.max(data)
  if args.cmin is not None:
    cmin = float(args.cmin)
  if args.cmax is not None:
    cmax = float(args.cmax)
  if (args.csym):
    cmin = min(cmin, -np.abs(cmax))
    cmax = max(cmax, np.abs(cmin))

  fsize = args.fsize

  if (args.k):
    fc = 'k'
    tc = 'w'
  else:
    fc = 'w'
    tc = 'k'

  if (args.gridlines):
    ecol = 'k'
  else:
    ecol = "None"

  fig = plt.figure(num=None, figsize=(14, 7), dpi=args.dpi, facecolor=fc,frameon=args.bgtrans)
  ax = plt.gca()

  try:
    if (args.smooth):
      plot_h = plt.pcolormesh(xvec_plot, yvec_plot, data, edgecolor=ecol, linewidth=0.01, shading='gouraud')
    else:
      plot_h = plt.pcolormesh(xvec_plot, yvec_plot, data, edgecolor=ecol, linewidth=0.01)

  except:
    print(xvec_plot)
    print(yvec_plot)
    print(data.shape)
    raise

  plt.set_cmap(args.cmap)
  plt.clim([cmin, cmax])

  xmn = np.min(xvec)
  xmx = np.max(xvec)


  init_locs = plt.xticks()[0]
  locs, labels, utstartSecs = get_xticks(args,xmn,xmx,init_locs)
  xaxis_TicksLabel(args,locs,labels,tc,ax,utstartSecs)

  plt.xlim(xmin=xmin, xmax=xmax)
  plt.ylim(ymin=ymin, ymax=ymax)

  xdtmi = np.abs(xmax-xmin)/15.0

  if args.ymin is None:
    plt.ylim(ymin=-90)
  if args.ymax is None:
    plt.ylim(ymax=90)
  if args.ymin is None and args.ymax is None:
    plt.yticks((-90,-45, 0, 45, 90),
      ('-90', '-45', '0', '45', '90')
    )
    ydtmi = 5
  else:
    locs, _ = plt.yticks()
    labels = np.round(np.rad2deg(locs))
    labels = [int(x) for x in labels]
    locs = np.deg2rad(labels)
    labels = [str(x) for x in labels]
    plt.yticks(locs,
      labels
    )
    ydtmi = np.abs(ymax-ymin)/15.0

  plt.ylabel(r'Latitude ($^\circ$)', {'fontsize': fsize, 'color': tc})

  if (args.finegrid):
      minor_xticks = np.arange(xmin, xmax, xdtmi)
      minor_yticks = np.arange(ymin, ymax, ydtmi)
      ax.set_xticks(minor_xticks, minor=True)
      ax.set_yticks(minor_yticks, minor=True)
      ax.grid(which='minor', alpha=0.5)


  ax_divider = make_axes_locatable(ax)
  if (args.x_cbar):
    cb_ax = ax_divider.append_axes("bottom", size="5%", pad=.8)
    cb = plt.colorbar(
      plot_h, cax=cb_ax, orientation="horizontal"
    )
    cb_ax.set_xlabel(cbstr, fontsize=fsize, color=tc)
  else:
    cb_ax = ax_divider.append_axes("right", size="2.5%", pad=.2)
    cb = plt.colorbar(
      plot_h, cax=cb_ax, orientation="vertical"
    )
    cb_ax.set_ylabel(cbstr, fontsize=fsize, color=tc)

  cb.outline.set_edgecolor(tc)
  cb_ax.tick_params(axis='both', color=tc, labelsize=fsize, length=6)

  cb.ax.tick_params(axis='y',colors=tc)

  if (args.uscales):
    plt.xlabel('Grid Index', {'fontsize': fsize, 'color': tc})
    plt.ylabel('Grid Index', {'fontsize': fsize, 'color': tc})

  ax.tick_params(axis='both', colors=tc, labelsize=fsize)

  ax.set_facecolor(fc)
  ax.spines['bottom'].set_color(tc)
  ax.spines['top'].set_color(tc)
  ax.spines['right'].set_color(tc)
  ax.spines['left'].set_color(tc)

  if (args.title is not None):
    ax.set_title(args.title, color=tc, size=fsize)

  if (not args.nogrid):
    ax.grid()

  if (args.noscales):
    plot_h.axes.get_xaxis().set_visible(False)
    plot_h.axes.get_yaxis().set_visible(False)

  if (args.noax):
    ax.set_frame_on(False)
    cb.remove()
    plot_h.axes.get_xaxis().set_visible(False)
    plot_h.axes.get_yaxis().set_visible(False)
    plt.axis('off')

  if (args.nocb):
    cb.remove()

  if args.verbose:
    print('{}'.format(oFile))
  fig.savefig(oFile, bbox_inches="tight", pad_inches=0, dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)

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

  if not args.xunits:
    total_time = xmx - xmn
    if total_time < 2*minutes:
      args.xunits = 'seconds'
    elif total_time < 2*hours:
      args.xunits = 'minutes'
    elif total_time < 7*days:
      args.xunits = 'hours'
    elif total_time < years:
      args.xunits = 'days'
    elif total_time < 2*years:
      args.xunits = 'cr'
      if not args.xcadence:
        args.xcadence = 2
    elif total_time < 2*years:
        args.xunits = 'days'
    else:
      args.xunits = 'years'
  
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
      
  if (args.ignore_data_uttime):
    utstartSecs = 0
  else:
    utstartSecs = xmn
  initLocs_uttime = init_locs*3600
  xmn_uttime = xmn
  xmx_uttime = xmx
  if args.xunits == "date":
    if (args.ignore_data_uttime):
      raise Exception("Did not specify a utstart")
    else:
      locs, labels = date_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime)
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
    if (args.ignore_data_uttime): 
      locs, labels = since_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime,cr)
    else:
      locs, labels = cr_xticks(args,xcUnitsSec,xmn_uttime,xmx_uttime)
  elif args.xunits == "years":
    locs, labels = since_xticks(args,xcUnitsSec,initLocs_uttime,xmn_uttime,xmx_uttime,years)
  locs = np.array(locs)/3600
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
    else: 
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Seconds since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)   
  elif args.xunits == "minutes":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    else: 
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Minutes since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)  
  elif args.xunits == "hours":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    else: 
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Hours since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)  
  elif args.xunits == "days":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    else: 
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Days since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)  
  elif args.xunits == "weeks":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    else: 
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Weeks since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
    ax.tick_params(axis='x',labelsize=args.xlabelfsize)  
  elif args.xunits == "cr":
    if (args.slant):
      plt.xticks(locs,labels,rotation=int(args.slant), ha=args.ha, rotation_mode=rmode, ma=args.ma)
    else:
      plt.xticks(locs,labels, ha=args.ha, ma=args.ma) 
    if (args.xlabel):
      plt.xlabel(args.xlabel, {'fontsize': args.fsize, 'color': tc})
    elif not (args.ignore_data_uttime): 
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
    else: 
      tlabel = gettLabel(args.tai,utstartSecs)
      plt.xlabel('Years since '+ tlabel, {'fontsize': args.fsize, 'color': tc})
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
