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
    help='Error to use for plotting: 1: RMSD, 2: MAXABS, 3: CVRMSD, 4: MAPE, 5: MAXAPE, 6: HHABS (Default: 3)',
    dest='errtype',
    type=int,
    required=False,
    default=3) 

  parser.add_argument('-lgl',
    help='A comma separated list of legend labels',
    dest='lgl',
    type=str,
    required=False) 

  parser.add_argument('-dpi',
    help='DPI for plots.',
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

  parser.add_argument('-ms',
    help='Marker size',
    dest='ms',
    type=float,
    default=2000,
    required=False)

  parser.add_argument('-lms',
    help='Legend marker size',
    dest='lms',
    type=float,
    default=2000,
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
    if len(lgl_list) != LABEL_LEN:
      print('Not enough markers provided')
      return
  else:
    MARKERS = ['o','v','^','<','>','8','s','p','*','h','H','D','d','P','X']

  if args.colors:
    COLORS=arg_dict['colors'].split(',')
    if len(lgl_list) != LABEL_LEN:
      print('Not enough colors provided')
      return
  else:
    cmap = plt.get_cmap('rainbow',LABEL_LEN)
    COLORS = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]  

  etype = args.errtype+1

  fsize=args.fsize
  ms=args.ms
  fig, ax = plt.subplots(num=None, figsize=(int(fgsize[0]),int(fgsize[1])), dpi=args.dpi, facecolor='w')

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

  if args.xlabel:
    xlabel=args.xlabel
  else:
    if 'AdvPhG' in testtype:
      xlabel='$\Delta \phi$'
    elif 'AdvPhSc' in testtype:
      xlabel='$\Delta \phi$'
    elif 'AdvThG' in testtype:
      xlabel='$\Delta \\theta$'
    elif 'AdvPhThG' in testtype:
      xlabel='$\Delta \\theta,\Delta \phi$'
    elif 'DiffSc' in testtype:
      xlabel='$\Delta \\theta,\Delta \phi$'
    elif 'DiffAdvPhSc' in testtype:
      xlabel='$\Delta \\theta,\Delta \phi$'

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
      ylabel='HHABS'


  ax.set_ylabel(ylabel,fontsize=fsize)
  ax.set_xlabel(xlabel,fontsize=fsize)
  ax.tick_params(axis='y',labelsize=fsize)
  ax.tick_params(axis='x',labelsize=fsize)
  ax.set_xscale("log")
  ax.set_yscale("log")
  ax.set_aspect('equal')
  ax.grid(zorder=0)

  xfirst=int(1/ ((dsize/nxl[0])/np.pi) )
  np_list=['$\\frac{\pi}{'+str(xfirst * 2 ** (n - 1))+'}$' for n in range(1, len(nxl) + 1)]

  plt.xticks(xvec,np_list)
  plt.minorticks_off()

  lg = plt.legend(fontsize=fsize,loc='best')

  for handle in lg.legend_handles:
    handle._sizes= [args.lms]

  ofile = args.ofile if args.ofile else testtype+'.png'

  fig.tight_layout()                               
  fig.savefig(ofile,bbox_inches='tight',dpi=args.dpi)
  plt.close('all')


def main():
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()


