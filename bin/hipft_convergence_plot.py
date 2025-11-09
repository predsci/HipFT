#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')

def argParsing():
  parser = argparse.ArgumentParser(description='HipFt Convergence Plotting.')

  parser.add_argument('-cfile',
    help='Convergence file',
    dest='cfile',
    type=str,
    required=True)

  parser.add_argument('-ofile',
    help='Output plot file name',
    dest='ofile',
    type=str,
    required=False)

  parser.add_argument('-errtype',
    help='Error to use for plotting: 1: RMSD, 2: MAXABS, 3: CVRMSD, 4: MAPE, 5: MAXAPE, 6: HHABS (Default: 6)',
    dest='errtype',
    type=int,
    required=False,
    default=6) 

  parser.add_argument('-lgl',
    help='A comma separated list of legend labels',
    dest='lgl',
    type=str,
    required=False) 

  parser.add_argument('-dpi',
    help='DPI for plots (default 120).',
    dest='dpi',
    default=120,
    type=int,
    required=False)

  parser.add_argument('-figsize',
    help='Comma separated figure size (Default: 96,48)',
    dest='fgsize',
    type=str,
    default='96,48',
    required=False)

  parser.add_argument('-fsize',
    help='Font size',
    dest='fsize',
    default=50,
    required=False)

  parser.add_argument('-xlabel',
    help='Label for x axis',
    dest='xlabel',
    required=False)

  parser.add_argument('-ylabel',
    help='Label for y axis',
    dest='ylabel',
    required=False)

  parser.add_argument('-ymax',
    help='Y axis max',
    dest='ymax',
    type=float,
    required=False)

  parser.add_argument('-ms',
    help='Marker size',
    dest='ms',
    type=float,
    default=4000,
    required=False)

  parser.add_argument('-lms',
    help='Legend marker size',
    dest='lms',
    type=float,
    default=3000,
    required=False)

  parser.add_argument('-markers',
    help='Comma separated marker list',
    dest='markers',
    type=str,
    required=False)

  parser.add_argument('-colors',
    help='Comma separated color list',
    dest='colors',
    type=str,
    required=False)

  parser.add_argument('-orders',
    help='Comma separated list of order lines to plot (can use [] to skip some)',
    dest='orders',
    type=str,
    required=False)

  parser.add_argument('-order_tethers',
    help='Comma separated list of indices to tether order lines to',
    dest='order_tethers',
    type=str,
    required=False)

  parser.add_argument('-order_labels',
    help='Comma separated list of legend labels for orders',
    dest='order_labels',
    type=str,
    required=False)
    
  parser.add_argument('-lw',
    help='Line width',
    dest='lw',
    type=float,
    default=10,
    required=False)    

  return parser.parse_args()


def run(args):
  arg_dict = vars(args)
  cfile_list = arg_dict['cfile'].split(',')
  fgsize = arg_dict['fgsize'].split(',')

  LABEL_LEN = len(cfile_list)

  create_lgl=False
  if args.lgl:
    lgl_list=arg_dict['lgl'].split(',')
    if len(lgl_list) != LABEL_LEN:
      print('Not enough legend labels provided')
      return
    elif len(lgl_list) != LABEL_LEN:
      print('Too many legend labels provided')
      return    
  else:
    create_lgl=True

  if args.markers:
    MARKERS=arg_dict['markers'].split(',')
    if len(MARKERS) != LABEL_LEN:
      print('Not enough markers provided')
      return
  else:
    MARKERS = ['o','v','^','<','>','8','s','p','*','h','H','D','d','P','X']

  if args.colors:
    COLORS=arg_dict['colors'].split(',')
    if len(COLORS) != LABEL_LEN:
      print('Not enough colors provided')
      return
  else:
    cmap = plt.get_cmap('rainbow',LABEL_LEN)
    COLORS = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]  

  if args.orders:
    ORDERS=arg_dict['orders'].split(',')

  if args.order_tethers:
    ORDER_TETHERS=arg_dict['order_tethers'].split(',')
  elif(args.orders):
    ORDER_TETHERS = [-1 for i in range(len(ORDERS))]

  if args.order_labels:
    ORDER_LABELS=arg_dict['order_labels'].split(',')
  elif(args.orders):
    ORDER_LABELS = ['$O(h^'+str(ORDERS[i])+')$' for i in range(len(ORDERS))]

  lw=args.lw

  etype = args.errtype+1

  fsize=args.fsize
  ms=args.ms
  fig, ax = plt.subplots(num=None, figsize=(int(fgsize[0]),int(fgsize[1])), dpi=args.dpi, facecolor='w')

  ox0=[]
  ox1=[]
  oy0=[]    
  oy1=[] 

  for (idx,cfile) in enumerate(cfile_list):
    nxl=[]
    err=[]   
    testtype=None
    with open(cfile, 'r') as f:
      testtype=f.readline().replace('\n','')
      if create_lgl:
        if 'fnm1' in testtype:
          leglabel='FE+Upwind'
        elif 'fnm2' in testtype:
          leglabel='RK3TVD+Upwind'
        elif 'fnm3' in testtype:
          leglabel='RK3TVD+WENO3'
        elif 'fnm4' in testtype:
          leglabel='SSPRK(4,3)+WENO3'
        else:
          leglabel=''
        if 'dnm1' in testtype:
          leglabel+=' dnm1'
        elif 'dnm2' in testtype:
          leglabel+=' dnm2'
        elif 'dnm3' in testtype:
          leglabel+=' dnm3'
        if 'strang' in testtype:
          leglabel+=' strang'
      else:
        leglabel=lgl_list[idx]

      if 'AdvThG' in testtype:
        xindx = 1
        dsize = np.pi
      else:
        xindx = 0
        dsize = 2*np.pi

      f.readline()
      f.readline()
      f.readline()
      for line in f:
        if '])' in line:
          break
        line = line.replace('[','').replace('], \\','')
        data = line.split(',')
        nxl.append(int(data[xindx]))
        err.append(float(data[etype]))

    xvec = dsize/np.array(nxl)
    yvec = np.array(err)

    plt.scatter(xvec,yvec,s=ms,c=COLORS[idx],edgecolors=COLORS[idx],zorder=3,marker=MARKERS[idx],label=leglabel)
    
    if args.orders and idx <= len(ORDERS)-1:
      xvec = dsize/np.array(nxl)
      yvec = np.array(err)    
      xtether = xvec[int(ORDER_TETHERS[idx])]
      ytether = yvec[int(ORDER_TETHERS[idx])]
      xmin = xvec[0]
      xmax = xvec[-1]
      ox0.append(float(xmin))
      oy0.append(float(10.0**(np.log10(ytether)-float(ORDERS[idx])*(np.log10(xtether)-np.log10(xmin)))))
      ox1.append(float(xmax))
      oy1.append(float(10.0**(np.log10(ytether)+float(ORDERS[idx])*(np.log10(xmax)-np.log10(xtether)))))


  if args.orders:
    for idx in range(len(ORDERS)):
      leglabel = '$O(h^'+str(ORDERS[idx])+')$'
      ax.plot([ox0[idx],ox1[idx]],[oy0[idx],oy1[idx]],'--',c=COLORS[idx],linewidth=lw,label=leglabel)


  if args.xlabel:
    xlabel=args.xlabel
  else:
    if 'AdvPhG' in testtype:
      xlabel=r'$\Delta \phi$'
    elif 'AdvPhSc' in testtype:
      xlabel=r'$\Delta \phi$'
    elif 'AdvThG' in testtype:
      xlabel=r'$\Delta \theta$'
    elif 'AdvPhThG' in testtype:
      xlabel=r'$\Delta \theta,\Delta \phi$'
    elif 'DiffSc' in testtype:
      xlabel=r'$\Delta \theta,\Delta \phi$'
    elif 'DiffAdvPhSc' in testtype:
      xlabel=r'$\Delta \theta,\Delta \phi$'

  if args.ylabel:
    ylabel=args.ylabel
  else:
    if etype == 2:
      ylabel='RMSD'
    elif etype == 3:
      ylabel='MAXABS'
    elif etype == 4:
      ylabel='(CV)RMSD'
    elif etype == 5:
      ylabel='MAPE'
    elif etype == 6:
      ylabel='MAXAPE'
    elif etype == 7:
      ylabel='HH$_{||}$'


  ax.set_ylabel(ylabel,fontsize=fsize)
  ax.set_xlabel(xlabel,fontsize=fsize)
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.set_xscale("log")
  ax.set_yscale("log")
  if (args.ymax):
      ax.set_ylim(ymax=args.ymax)
  ax.set_aspect('equal')
  ax.grid(zorder=0)

  xfirst=int(1/ ((dsize/nxl[0])/np.pi) )
  np_list=[r'$\frac{\pi}{'+str(xfirst * 2 ** (n - 1))+'}$' for n in range(1, len(nxl) + 1)]

  plt.xticks(xvec,np_list)
  plt.minorticks_off()

  lg = plt.legend(fontsize=fsize,loc='lower right')

  for handle in lg.legend_handles:
    handle._sizes= [args.lms]

  fig.tight_layout()

  ofile = args.ofile if args.ofile else testtype+'.png'
  fig.savefig(ofile,bbox_inches='tight',dpi=args.dpi)

  ofile = args.ofile if args.ofile else testtype+'.pdf'
  fig.savefig(ofile,bbox_inches='tight',dpi=args.dpi)

  plt.close('all')


def main():
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()


