#!/usr/bin/env python3
import argparse
import signal
import sys
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import psihdf as ps
import psipals
import psimath

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
            default=' ',
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

  parser.add_argument('-xunits',
            help='Units of the x-axis.',
            dest='xunits',
            default='years',
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

  parser.add_argument('-xc',
            help='Cadence of xunits to skip (Default ever year).',
            dest='xcadence',
            default=1,
            required=False)

  return parser.parse_args()

def run(args):
  # Load colormaps:
  psipals.load()

  # Read the data:
  xvec, yvec, data = ps.rdhdf_2d(args.iFile)

  #thin out data
  skip = int(args.cadence)+1
  xvec = xvec[::skip]
  data = data[:,::skip]

  if (len(xvec) == 0):
    xvec = np.array(range(0, len(data[0, :])))
  if (len(yvec) == 0):
    yvec = np.array(range(0, len(data[:, 0])))

  # Set up data units:
  data = np.multiply(data, float(args.unit_fac))
  cbstr = args.unit_label

  # Set up uniform grid scales if requested:
  if (args.uscales):
    xvec_plot = np.array(range(0, len(xvec)))
    yvec_plot = np.array(range(0, len(yvec)))
  else:
    xvec_plot = xvec
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
  else:
    ymin = 0
  if args.ymax is not None:
    ymax = float(args.ymax)
  else:
    ymax = np.pi

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

  fsize = 18

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

  plt.xlim(xmin=xmin, xmax=xmax)
  plt.ylim(ymin=ymin, ymax=ymax)


  if args.xunits == "years":
    startTimeYear = datetime.fromtimestamp(xvec_plot[0]).strftime('%Y')
    if xvec_plot[0] > datetime.strptime(startTimeYear,'%Y').timestamp():
      startTimeYear=int(startTimeYear)+1
    else:
      startTimeYear=int(startTimeYear)
    endTimeYear = datetime.fromtimestamp(xvec_plot[-1]).strftime('%Y')
    if xvec_plot[-1] < datetime.strptime(endTimeYear,'%Y').timestamp():
      endTimeYear=int(endTimeYear)-1
    else:
      endTimeYear=int(endTimeYear)

    locsYear = list(range(startTimeYear, endTimeYear+1))
    locs = [datetime.strptime(str(jtime),'%Y').timestamp() for jtime in locsYear]

    labels = [str(itime) for itime in locsYear]

    xskip = int(args.xcadence)
    locs = locs[::xskip]
    labels = labels[::xskip]

    plt.xticks(locs,labels)
    plt.xlabel(r'Years', {'fontsize': fsize, 'color': tc})

  xdtmi = np.abs(xmax-xmin)/15.0

  if args.ymin is None:
    plt.ylim(ymin=0)
  if args.ymax is None:
    plt.ylim(ymax=np.pi)
  if args.ymin is None and args.ymax is None:
    plt.yticks((0, np.pi*0.25, np.pi*0.5, np.pi*0.75, np.pi),
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

  if (args.oFile is None):
    if (str(args.iFile).endswith('h5')):
      n = 3
    else:
      n = 4
    oFile = args.iFile[0:len(args.iFile) - n] + ".png"
  else:
    oFile = args.oFile

  if args.verbose:
    print('{}'.format(oFile))
  fig.savefig(oFile, bbox_inches="tight", pad_inches=0, dpi=args.dpi, facecolor=fig.get_facecolor(), edgecolor=None)


def main():
  ## Get input agruments:
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()

