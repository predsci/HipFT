#!/usr/bin/env python3

import glob
import re
import os
import argparse
import signal

def handle_int(signum, frame):
    print('You pressed Ctrl+C! Stopping!')
    signal.signal(signum, signal.SIG_DFL)
    os.kill(os.getpid(), signum)
    sys.exit(1)
    
signal.signal(signal.SIGINT, handle_int)

def argParsing():
  parser = argparse.ArgumentParser(description='hipft_post_process_run:  This tool makes all plots for a run.')

  parser.add_argument('-rundir',
    help='Path to the directory where hipft was run.  Default is current directory.',
    dest='rundir',
    required=False)

  parser.add_argument('-outpath',
    help='Path to folder where output data is. Default: args.rundir/output_maps',
    dest='outpath',
    required=False)

  parser.add_argument('-basefn',
    help='Base file name for output maps.',
    dest='basefn',
    required=False)

  parser.add_argument('-utstart',
    help='Start Date in UT: YYYY-MM-DDTHH:MM:SS',
    dest='utstart',
    required=False)

  parser.add_argument('-maplist',
    help='Full path to the map file list.  Default: args.rundir/hipft_output_map_list.out',
    dest='maplist',
    required=False)

  parser.add_argument('-t0',
    help='Sequence start index.',
    dest='t0',
    required=False)

  parser.add_argument('-tf',
    help='Sequence stop index. If not specified, all maps in mapfile are used.',
    dest='tf',
    required=False)

  parser.add_argument('-tfac',
    help='Time factor to get time into units of hours (HipFT time is in hours, so default tfac=1).',
    dest='tfac',
    required=False)
 
  parser.add_argument('-xunits',
    help='Units of the x-axis (date, seconds, minutes, hours, days, weeks, cr, or years).',
    dest='xunits',
    default='hours',
    required=False)

  parser.add_argument('-nomovie',
    help='Do not make mov movie.',
    dest='nomovie',
    default=False,
    action='store_true',
    required=False)

  parser.add_argument('-hipftbin',
    help='HipFT bin location',
    dest='hipft',
    type=str,
    required=False)

  parser.add_argument('-samples',
    help='Number of points to plot, this helps with larger files (default:300)',
    default=300)

  parser.add_argument('-cores',
    help='Number of cores to use for movies.',
    dest='cores',
    type=int,
    required=False)

  parser.add_argument('-butterfly',
    help='Create butterfly plots only',
    dest='butterfly',
    default=False,
    action='store_true',
    required=False)

  parser.add_argument('-movies',
    help='Create movies only',
    dest='movies',
    default=False,
    action='store_true',
    required=False)

  parser.add_argument('-history',
    help='Create history plots only, Options{ 1: All, 2: individual, 3: together, 4: summary}',
    dest='history',
    type=int,
    required=False)

  parser.add_argument('-force',
    help='Force remaking the plots',
    dest='force',
    default=False,
    action='store_true',
    required=False)

  return parser.parse_args()


def run(args):

  if not args.cores:
    args.cores = int(os.getenv('OMP_NUM_THREADS', 1))

  run_all = not (args.butterfly or args.movies or args.history)

  args.hipftbin = args.hipftbin if args.hipft else os.popen('which hipft').read().replace('\n','').replace('bin/hipft','bin/')
  args.rundir = args.rundir or os.getcwd()
  cur_dir = os.getcwd()

  ## ~~~~~~ Get number of realizations
  pattern = os.path.join(args.rundir, 'hipft_history_sol*.out')
  file_list = glob.glob(pattern)
  check_exist = [os.path.join(cur_dir,os.path.basename(file)) for file in file_list]
  is3d = len(file_list) > 1

  ## ~~~~~~ Plot individual histories 
  if run_all or args.history == 1 or args.history == 2:
    if args.force:
      hist_ind = file_list
    else:
      hist_ind = []
      for check, file in zip(check_exist, file_list):
        pattern = check.replace('hipft_history_sol', 'history').replace('.out', '*.png')
        found_list = glob.glob(pattern)
        if len(found_list) != 9:
          hist_ind.append(file)

    if hist_ind:
        plot_individual_history(args, hist_ind)
    else:
      print("All individual history plots exist")

  ## ~~~~~~ Plot all histories together
  if (run_all or args.history == 1 or args.history == 3) and is3d:
    hist_ind_flag = args.force or len(glob.glob(os.path.join(cur_dir, 'history_together*.png'))) != 9

    if hist_ind_flag:
      plot_together_history(args, file_list)
    else:
      print("Together history plots exist")

  ## ~~~~~~ Plot histories summary
  if (run_all or args.history == 1 or args.history == 4) and is3d:
    hist_ind_flag = args.force or len(glob.glob(os.path.join(cur_dir, 'history_summary*.png'))) != 9
    
    if hist_ind_flag:
      plot_summary_history(args, file_list)
    else:
      print("Summary history plots exist")

  ## ~~~~~~ Make butterfly h5
  if run_all or args.butterfly:
    if args.force or not os.path.exists(os.path.join(cur_dir, 'butterfly.h5')):
      make_butterfly(args, is3d)
    else:
      print("butterfly h5 exists")

  ## ~~~~~~ Plot butterfly h5
    butterfly_plots = []
    if not args.force:
      if is3d:
        for file in check_exist:
          butterfly_file = file.replace('hipft_history_sol', 'butterfly').replace('.out', '.png')
          if not os.path.exists(butterfly_file):
            match = re.search(r'r(\d+)', butterfly_file)
            r = int(match.group(1)) if match else -1
            butterfly_plots.append(str(r-1))
      else:
        if not os.path.exists(os.path.join(cur_dir, 'butterfly.png')):
            butterfly_plots.append('-1')


    if butterfly_plots or args.force:
      plot_butterfly(args, is3d, None if len(butterfly_plots) == len(file_list) else butterfly_plots)
    else:
      print("All butterfly plots exist")

  ## ~~~~~~ Make movies
  if run_all or args.movies:
    movie_ind = []

    if not args.outpath:
      args.outpath = os.path.join(args.rundir, 'output_maps')

    if not args.force:
      if is3d:
        for file in check_exist:
          pattern = file.replace('hipft_history_sol', 'hipft_movie').replace('.out', '*.mov')
          found_list = glob.glob(pattern)
          pattern = os.path.join(args.outpath, 'hipft_brmap_idx*.h5')
          num_frames = glob.glob(pattern)
          match = re.search(r'r(\d+)', file)
          r = match.group() if match else ""
          pattern = os.path.join(args.outpath, 'plots', f'hipft_brmap_idx*{r}.png')
          num_plots = glob.glob(pattern)

          if not found_list or len(num_plots) != len(num_frames):
            r = int(match.group(1)) if match else -1
            movie_ind.append(str(r-1))
      else:
        if not os.path.exists(os.path.join(cur_dir, 'hipft_movie.mov')):
          butterfly_plots.append('1')

    if movie_ind or args.force:
      make_movies(args, None if len(movie_ind) == len(file_list) else movie_ind)
    else:
      print("All movies exist")


def plot_individual_history(args, hist_ind):
  print("Plotting histories individual...")

  bc_plotHistories = os.path.join(args.hipftbin, 'hipft_plot_histories.py')
  bc_plotHistories += f' -samples {args.samples}'

  if args.tfac:
    bc_plotHistories += f' -tfac {args.tfac}'
  if args.utstart:
    bc_plotHistories += f' -utstart {args.utstart}'

  for file in hist_ind:
    match = re.search(r'r(\d+)', file)
    realization = match.group() if match else 'individual'
    os.system(f'{bc_plotHistories} -histfiles {file} -runtag {realization}')


def plot_together_history(args, file_list):
  print("Plotting histories together...")

  bc_plotHistories = os.path.join(args.hipftbin, 'hipft_plot_histories.py')
  bc_plotHistories += f' -samples {args.samples}'

  if args.tfac:
    bc_plotHistories += f' -tfac {args.tfac}'
  if args.utstart:
    bc_plotHistories += f' -utstart {args.utstart}'

  histfiles = ','.join(file_list)
  os.system(f'{bc_plotHistories} -histfiles {histfiles} -runtag together')  


def plot_summary_history(args, file_list):
  print("Plotting histories summary...")

  bc_plotHistories = os.path.join(args.hipftbin, 'hipft_plot_histories.py')
  bc_plotHistories += f' -samples {args.samples}'

  if args.tfac:
    bc_plotHistories += f' -tfac {args.tfac}'
  if args.utstart:
    bc_plotHistories += f' -utstart {args.utstart}'

  histfiles = ','.join(file_list)
  os.system(f'{bc_plotHistories} -summary -histfiles {histfiles} -runtag summary')


def make_butterfly(args, is3d):
  print("Make butterfly h5...")

  bc_makeButterfly = os.path.join(args.hipftbin, 'hipft_make_butterfly_diagram.py')
  if args.rundir:
    bc_makeButterfly += f' -rundir {args.rundir}'
  if args.outpath:
    bc_makeButterfly += f' -outpath {args.outpath}'
  if args.maplist:
    bc_makeButterfly += f' -maplist {args.maplist}'
  if args.basefn:
    bc_makeButterfly += f' -basefn {args.basefn}'
  if args.utstart:
    bc_makeButterfly += f' -utstart {args.utstart}'
  if args.t0:
    bc_makeButterfly += f' -t0 {args.t0}'
  if args.tf:
    bc_makeButterfly += f' -tf {args.tf}'
  if args.tfac:
    bc_makeButterfly += f' -tfac {args.tfac}'
  if is3d:
    bc_makeButterfly += ' -3d -sall'
  
  os.system(bc_makeButterfly)


def plot_butterfly(args, is3d, butterfly_plots):
  print("Plot butterfly h5...")

  bc_plotButterfly = os.path.join(args.hipftbin, 'hipft_plot_butterfly_diagram.py')
  bc_plotButterfly += ' butterfly.h5'

  if args.xunits:
    bc_plotButterfly += f' -xunits {args.xunits}'
  if (args.cores):
    bc_plotButterfly += f' -cores {args.cores}'
  if is3d:
    bc_plotButterfly += ' -3d'
    if butterfly_plots:
      bc_plotButterfly += f' -slices {" ".join(butterfly_plots)}'
    else:
      bc_plotButterfly += ' -sall'

  os.system(bc_plotButterfly)


def make_movies(args, movie_ind):
  print("Making movies...")

  bc_makeMovie = os.path.join(args.hipftbin, 'hipft_make_plots_and_movies.py')
  bc_makeMovie += f' -dir {args.outpath}'

  if (args.nomovie):
    bc_makeMovie += ' -nomovie'
  if (args.cores):
    bc_makeMovie += f' -cores {args.cores}'

  if movie_ind:
    for r in movie_ind: 
      os.system(f'{bc_makeMovie} -slice {r}')
  else:
    os.system(bc_makeMovie)


def main():
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()
