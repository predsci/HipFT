#!/usr/bin/env python3

import glob
import re
import os
import argparse
import signal
import numpy as np
import pandas as pd
import shutil

# Version 1.3.0

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

# This may not be needed....
  parser.add_argument('-maplist',
    help='Full path to the map file list.  Default: args.rundir/hipft_output_map_list.out',
    dest='maplist',
    required=False)

# These will be changed to UT start and stop later...
  parser.add_argument('-t0',
    help='Sequence start index.',
    dest='t0',
    required=False)
  parser.add_argument('-tf',
    help='Sequence stop index. If not specified, all maps in mapfile are used.',
    dest='tf',
    required=False)

  parser.add_argument('-xunits',
    help='Units of the x-axis (date, seconds, minutes, hours, days, weeks, cr, or years).',
    dest='xunits',
    required=False)

  parser.add_argument('-plot_axis_options',
    help='String of plot options for the axises.',
    dest='plot_axis_options',
    required=False)

  parser.add_argument('-map_plot_options',
    help='String of plot options for the history plots.',
    dest='map_plot_options',
    required=False)

  parser.add_argument('-movie_plot_options',
    help='String of plot options for the images and movies plotted.',
    dest='movie_plot_options',
    required=False)
  
  parser.add_argument('-butterfly_plot_options',
    help='String of plot options for the butterfly plots.',
    dest='butterfly_plot_options',
    default="-cmin -10 -cmax 10",
    required=False)

  parser.add_argument('-hipft_home',
    help='HipFT bin location',
    dest='hipft_home',
    type=str,
    required=False)

  parser.add_argument('-history_plot_samples',
    help='Number of points to plot, this helps with larger files (default: 300) (Can set to "all" to plot all points)',
    dest='history_plot_samples',
    default=300,
    required=False)

  parser.add_argument('-history_plot_samples_markers',
    help='Number of marker points to plot (default: value of -history_plot_samples) (Can set to "all" to plot all points)',
    dest='history_plot_samples_markers',
    required=False)
  
  parser.add_argument('-serial',
    help='Set the number of threads to use to 1 instead of using the environment variable OMP_NUM_THREADS.',
    dest='serial',
    default=False,
    action='store_true',
    required=False)

  parser.add_argument('-butterfly',
    help='Create butterfly plots',
    dest='butterfly',
    default=False,
    action='store_true',
    required=False)

  parser.add_argument('-movies',
    help='Create movies',
    dest='movies',
    default=False,
    action='store_true',
    required=False)

  parser.add_argument('-history',
    nargs='?',
    const=1,
    help='Create history plots, Options{ 1: Plot all histories (2,3,4), \
                                         2: Plot individual realizations, \
                                         3: Plot all realizations together, \
                                         4: Plot realization summary (mean, stddev, etc)}, \
                                         If not provided, defaults to 1.',
    dest='history',
    type=int,
    required=False)

  parser.add_argument('-hist_dir',
    help='Location of history files (default is the run folder)',
    dest='hist_dir',
    type=str,
    required=False)

  parser.add_argument('-overwrite',
    help='Overwrite selected processing',
    dest='overwrite',
    default=False,
    action='store_true',
    required=False)

  parser.add_argument('-output_dir',
    help='Output directory of script (Default: current_folder/post_processing)',
    dest='output_dir',
    type=str,
    required=False)
  
  parser.add_argument('-tai_off',
    help='Turn off using tai which ignores leap seconds',
    dest='tai_off',
    default=False,
    action='store_true',
    required=False)

  return parser.parse_args()


def run(args):

  print('')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')
  print('                         OFT Post Processing                          ')
  print('')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

  args.num_threads = 1 if args.serial else int(os.getenv('OMP_NUM_THREADS', os.cpu_count() or 1))

  args.history_plot_samples_markers = args.history_plot_samples if not args.history_plot_samples_markers else args.history_plot_samples_markers

  run_all = not (args.butterfly or args.movies or args.history)

  args.hipft_home = args.hipft_home if args.hipft_home else os.path.dirname(os.path.normpath(os.popen('which hipft').read().strip()))

  print(f'HipFT bin found at :')
  print(f'\t {args.hipft_home}')
  print(f'Number of threads running on : {args.num_threads}')

  args.rundir = args.rundir or os.getcwd()
  if not args.output_dir:
    args.output_dir = os.path.join(os.getcwd(), 'post_processing')
  os.makedirs(args.output_dir, exist_ok=True)
  print(f'==> Created post processing output folder at :')
  print(f'\t {args.output_dir}')

  if run_all or args.history:
    history_folder = os.path.join(args.output_dir, "histories")
    os.makedirs(history_folder, exist_ok=True)
    print(f'==> Created map histories output folder at :')
    print(f'\t {history_folder}')

  ## ~~~~~~ Get number of realizations
  if args.hist_dir is not None:
    hdir = args.hist_dir
  else:
    hdir = args.rundir
  pattern = os.path.join(hdir, 'hipft_history_sol*.out')
  file_list = sorted(glob.glob(pattern))
  check_exist = [os.path.basename(file) for file in file_list]
  is3d = len(file_list) > 1

  ## ~~~~~~ Get xunits to use if not set
  if not args.xunits and (run_all or args.history or args.butterfly):
    hist_sol = pd.read_table(file_list[0],header=0,sep='\\s+')
    time_list = np.array(hist_sol['TIME'])
    total_time = (time_list[-1] - time_list[0])*3600
    minutes = 60
    hours = 3600
    days = 86400
    years = 31556952
    if total_time < 2*minutes:
      args.xunits = 'seconds'
    elif total_time < 2*hours:
      args.xunits = 'minutes'
    elif total_time < 7*days:
      args.xunits = 'hours'
    elif total_time < years:
      args.xunits = 'days'
    elif args.utstart and total_time < 2*years:
      args.xunits = 'cr'
      if not args.plot_axis_options or 'xcadence' not in args.plot_axis_options:
        args.plot_axis_options = " -xcadence 2"
    elif total_time < 2*years:
        args.xunits = 'days'
    else:
      args.xunits = 'years'

## ~~~~~~ Make realization parameters table
  if is3d:
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('==> Starting to make realization parameters table')
    if args.overwrite:
        print('==> Overwriting any existing realization parameters table!')

    table_flag = args.overwrite or not os.path.exists(os.path.join(args.output_dir, "hipft_realization_parameters_table.pdf"))
    if table_flag:
        make_realization_table(args)
    else:
        print("==> Realization parameters table already exists")

  ## ~~~~~~ Plot individual histories 
  if run_all or args.history == 1 or args.history == 2:
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('==> Starting to make individual history plots')
    
    if args.overwrite:
      print('==> Overwriting any existing individual history plots!')
      hist_ind = file_list
    else:
      hist_ind = []
      histories_directory = os.path.join(args.output_dir, "histories")
      for check, file in zip(check_exist, file_list):
        match = re.search(r'r(\d+)', file)
        r = match.group() if match else "r000001"
        check_file = check.replace('hipft_history_sol', 'history').replace('.out', '*.png')
        pattern = os.path.join(histories_directory, r, check_file)
        found_list = glob.glob(pattern)
        if len(found_list) != 9:
          hist_ind.append(file)
        else:
          print(f'==>        Found existing individual plots for {r}')

    if hist_ind:
        plot_individual_history(args, hist_ind)
    else:
      print("==> All individual history plots already exist")

  ## ~~~~~~ Plot all histories together
  if (run_all or args.history == 1 or args.history == 3) and is3d:
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('==> Starting to make plotted together history plots')

    if args.overwrite:
      print('==> Overwriting any existing plotted together history plots!')

    hist_ind_flag = args.overwrite or len(glob.glob(os.path.join(os.path.join(args.output_dir, "histories"), 'history_together*.png'))) != 9

    if hist_ind_flag:
      plot_together_history(args, file_list)
    else:
      print("==> All plotted together history plots already exist")

  ## ~~~~~~ Plot histories summary
  if (run_all or args.history == 1 or args.history == 4) and is3d:
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('==> Starting to make plotted summary history plots')

    if args.overwrite:
      print('==> Overwriting any existing summary history plots!')
    
    hist_ind_flag = args.overwrite or len(glob.glob(os.path.join(os.path.join(args.output_dir, "histories"), 'history_summary*.png'))) != 9
    
    if hist_ind_flag:
      plot_summary_history(args, file_list)
    else:
      print("==> All summary history plots already exist")

  ## ~~~~~~ Make butterfly h5
  if run_all or args.butterfly:
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('==> Starting to make butterfly plots')

    butterfly_folder = os.path.join(args.output_dir, "butterfly")
    os.makedirs(butterfly_folder, exist_ok=True)
    print(f'==> Created butterfly output folder at :')
    print(f'\t {butterfly_folder}')

    if args.overwrite:
      print('==> Overwriting any existing butterfly h5 file!')  

    if args.overwrite or not os.path.exists(os.path.join(os.path.join(args.output_dir, "butterfly"), 'butterfly.h5')):
      make_butterfly(args, is3d)
    else:
      print("==> All butterfly h5 files already exists")

  ## ~~~~~~ Plot butterfly h5
    butterfly_plots = []
    if not args.overwrite:
      if is3d:
        for file in check_exist:
          check_file = file.replace('hipft_history_sol', 'butterfly').replace('.out', '.png')
          butterfly_file = os.path.join(args.output_dir, "butterfly", check_file)
          match = re.search(r'r(\d+)', butterfly_file)
          if not os.path.exists(butterfly_file):
            r = int(match.group(1)) if match else -1
            butterfly_plots.append(str(r))
          else:
            print(f'==>        Found existing butterfly plots for {match.group()}')
      else:
        if not os.path.exists(os.path.join(args.output_dir, "butterfly", 'butterfly.png')):
            butterfly_plots.append('-1')
        else:
            print(f'==>         Found existing butterfly plots for r000001')

    if args.overwrite:
      print('==> Overwriting any existing butterfly plots!')  

    if butterfly_plots or args.overwrite:
      plot_butterfly(args, is3d, None if len(butterfly_plots) == len(file_list) else butterfly_plots)
    else:
      print("==> All butterfly plots already exist")

  ## ~~~~~~ Make movies


  if run_all or args.movies:
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('==> Starting to make movies and images')

    if args.overwrite:
      print('==> Overwriting any existing movies and images!')  

    movies_folder = os.path.join(args.output_dir, "map_plotting")
    os.makedirs(movies_folder, exist_ok=True)
    print(f'==> Created map plotting output folder at :')
    print(f'\t {movies_folder}')

  if run_all or args.movies:
    movie_ind = []

    if not args.outpath:
      args.outpath = os.path.join(args.rundir, 'output_maps')

    if not args.overwrite:
      if is3d:
        movies_directory = os.path.join(args.output_dir, "map_plotting")
        for file in check_exist:
          match = re.search(r'r(\d+)', file)
          r = match.group() if match else "r000001"
          check_file = file.replace('hipft_history_sol', 'hipft_movie').replace('.out', '*.mov')
          pattern = os.path.join(movies_directory, 'movies', check_file)
          found_list = glob.glob(pattern)
          pattern = os.path.join(args.outpath, 'hipft_brmap_idx*.h5')
          num_frames = len(glob.glob(pattern))
          pattern = os.path.join(movies_directory, 'images', r, f'hipft_brmap_idx*{r}.png')
          num_plots = glob.glob(pattern)
          if not found_list or len(num_plots) != num_frames:
            r = int(match.group(1)) if match else -1
            movie_ind.append(str(r))
          else:
            print(f'==>        Found existing movies and images for {r}')
      else:
        if not os.path.exists(os.path.join(args.output_dir, "map_plotting",  'movies', 'hipft_movie.mov')):
          movie_ind.append('individual')
        else:
            print(f'==>        Found existing movies and images for r000001')
      isSubset = len(movie_ind) != len(file_list)
    else:
      if is3d:
        movies_directory = os.path.join(args.output_dir, "map_plotting")
        for file in check_exist:
          match = re.search(r'r(\d+)', file)
          r = int(match.group(1)) if match else -1
          movie_ind.append(str(r))
      else:
          movie_ind.append('individual')
      isSubset = False

    if movie_ind or args.overwrite:
      make_movies(args, movie_ind, isSubset)
    else:
      print("==> All movies and images already exist")


def plot_individual_history(args, hist_ind):
  print("==> Plotting all individual history plots...")

  bc_plotHistories = os.path.join(args.hipft_home, 'hipft_plot_histories.py')
  bc_plotHistories += f' -tai' if not args.tai_off else ''
  bc_plotHistories += f" -samples {args.history_plot_samples}" if args.history_plot_samples != 'all' else ''
  bc_plotHistories += f" -samples_markers {args.history_plot_samples_markers}" if args.history_plot_samples_markers != 'all' else ''
  bc_plotHistories += f" -utstart {args.utstart}" if args.utstart else ''
  bc_plotHistories += f" {args.map_plot_options}" if args.map_plot_options else ''
  bc_plotHistories += f" -xunits {args.xunits}" if args.xunits else ''
  bc_plotHistories += f" {args.plot_axis_options}" if args.plot_axis_options else ''

  histories_directory = os.path.join(args.output_dir, "histories")
  for file in hist_ind:
    match = re.search(r'r(\d+)', file)
    realization = match.group() if match else 'individual'
    os.chdir(histories_directory)
    os.system(f'{bc_plotHistories} -histfiles {file} -runtag {realization}')
    os.chdir(args.output_dir)

  for file in hist_ind:
    match = re.search(r'r(\d+)', file)
    r = match.group() if match else "r000001"
    histories_sub_directory = os.path.join(histories_directory, r)
    os.makedirs(histories_sub_directory, exist_ok=True)
    os.system(f'mv {os.path.join(histories_directory,f"*{r}*.png")} {histories_sub_directory}')


def plot_together_history(args, file_list):
  print("==> Plotting all together history plots...")

  bc_plotHistories = os.path.join(args.hipft_home, 'hipft_plot_histories.py')
  bc_plotHistories += f' -tai' if not args.tai_off else ''
  bc_plotHistories += f" -samples {args.history_plot_samples}" if args.history_plot_samples != 'all' else ''
  bc_plotHistories += f" -samples_markers {args.history_plot_samples_markers}" if args.history_plot_samples_markers != 'all' else ''
  bc_plotHistories += f" -utstart {args.utstart}" if args.utstart else ''
  bc_plotHistories += f" {args.map_plot_options}" if args.map_plot_options else ''
  bc_plotHistories += f" -xunits {args.xunits}" if args.xunits else ''
  bc_plotHistories += f" {args.plot_axis_options}" if args.plot_axis_options else ''

  histfiles = ','.join(file_list)
  os.chdir(os.path.join(args.output_dir, "histories"))
  os.system(f'{bc_plotHistories} -histfiles {histfiles} -runtag together')
  os.chdir(args.output_dir)


def make_realization_table(args):
  print("==> Making realization parameters table...")
  bc_plotTable = os.path.join(args.hipft_home, 'hipft_make_realization_table.py')
  rpfile = os.path.join(args.rundir, 'hipft_realization_parameters.out')
  os.chdir(args.output_dir)
  os.system(f'{bc_plotTable} -rpfile {rpfile}')


def plot_summary_history(args, file_list):
  print("==> Plotting all summary history plots...")

  bc_plotHistories = os.path.join(args.hipft_home, 'hipft_plot_histories.py')
  bc_plotHistories += f' -tai' if not args.tai_off else ''
  bc_plotHistories += f" -samples {args.history_plot_samples}" if args.history_plot_samples != 'all' else ''
  bc_plotHistories += f" -samples_markers {args.history_plot_samples_markers}" if args.history_plot_samples_markers != 'all' else ''
  bc_plotHistories += f" -utstart {args.utstart}" if args.utstart else ''
  bc_plotHistories += f" {args.map_plot_options}" if args.map_plot_options else ''
  bc_plotHistories += f" -xunits {args.xunits}" if args.xunits else ''
  bc_plotHistories += f" {args.plot_axis_options}" if args.plot_axis_options else ''

  histfiles = ','.join(file_list)
  os.chdir(os.path.join(args.output_dir, "histories"))
  os.system(f'{bc_plotHistories} -summary -histfiles {histfiles} -runtag summary')
  os.chdir(args.output_dir)


def make_butterfly(args, is3d):
  print("==> Creating butterfly h5 file...")

  bc_makeButterfly = os.path.join(args.hipft_home, 'hipft_make_butterfly_diagram.py')
  bc_makeButterfly += f' -tai' if not args.tai_off else ''
  bc_makeButterfly += f" -rundir {args.rundir}" if args.rundir else ''
  bc_makeButterfly += f" -outpath {args.outpath}" if args.outpath else ''
  bc_makeButterfly += f" -maplist {args.maplist}" if args.maplist else ''
  bc_makeButterfly += f" -basefn {args.basefn}" if args.basefn else ''
  bc_makeButterfly += f" -utstart {args.utstart}" if args.utstart else ''
  bc_makeButterfly += f" -t0 {args.t0}" if args.t0 else ''
  bc_makeButterfly += f" -tf {args.tf}" if args.tf else ''
  bc_makeButterfly += f" -3d -sall" if is3d else ''

  mldfile = None

  if args.utstart and not args.maplist:
    make_mldfile = os.path.join(args.hipft_home, 'hipft_add_dates_to_map_output_list.py')
    make_mldfile += f' {args.utstart}'
    make_mldfile += f' -maplist {os.path.join(args.rundir, "hipft_output_map_list.out")}'
    make_mldfile += f' -o {os.path.join(args.output_dir, "hipft_output_map_list")}'
    os.system(make_mldfile)
    map_list_out = 'hipft_output_map_list_utc.out' if args.tai_off else 'hipft_output_map_list_tai.out'
    mldfile = os.path.join(args.output_dir, map_list_out)
    args.maplist = mldfile
  elif args.maplist:
    mldfile = args.maplist

  if mldfile:
    bc_makeButterfly += f' -maplist {mldfile}'

  os.chdir(os.path.join(args.output_dir, "butterfly"))
  os.system(bc_makeButterfly)
  os.chdir(args.output_dir)


def plot_butterfly(args, is3d, butterfly_plots):
  print("==> Plotting butterfly plots...")

  bc_plotButterfly = os.path.join(args.hipft_home, 'hipft_plot_butterfly_diagram.py')
  bc_plotButterfly += f' -tai' if not args.tai_off else ''
  bc_plotButterfly += ' butterfly.h5'

  bc_plotButterfly += f" -ignore_data_uttime" if not args.utstart else ''
  bc_plotButterfly += f" -xunits {args.xunits}" if args.xunits else ''
  bc_plotButterfly += f' -np {args.num_threads}'
  if is3d:
    bc_plotButterfly += ' -3d'
    bc_plotButterfly += f' -slices {" ".join(butterfly_plots)}' if butterfly_plots else ' -sall'
  bc_plotButterfly += f' {args.butterfly_plot_options}'
  bc_plotButterfly += f" {args.plot_axis_options}" if args.plot_axis_options else ''

  os.chdir(os.path.join(args.output_dir, "butterfly"))
  os.system(bc_plotButterfly)
  os.chdir(args.output_dir)


def make_movies(args, movie_ind, isSubset):
  print("==> Making movies and images...")

  mldfile = None

  if args.utstart and not args.maplist:
    make_mldfile = os.path.join(args.hipft_home, 'hipft_add_dates_to_map_output_list.py')
    make_mldfile += f' {args.utstart}'
    make_mldfile += f' -maplist {os.path.join(args.rundir, "hipft_output_map_list.out")}'
    make_mldfile += f' -o {os.path.join(args.output_dir, "hipft_output_map_list")}'
    os.system(make_mldfile)
    map_list_out = 'hipft_output_map_list_utc.out' if args.tai_off else 'hipft_output_map_list_tai.out'
    mldfile = os.path.join(args.output_dir, map_list_out)
  elif args.maplist:
    mldfile = args.maplist

  bc_makeMovie = os.path.join(args.hipft_home, 'hipft_make_plots_and_movies.py')
  bc_makeMovie += f' -dir {args.outpath}'
  bc_makeMovie += f' -np {args.num_threads}'
  if mldfile:
    bc_makeMovie += f' -mldfile {mldfile}'
  bc_makeMovie += f' {args.movie_plot_options}' if args.movie_plot_options else ''

  movies_directory = os.path.join(args.output_dir, "map_plotting", "movies")
  os.makedirs(movies_directory, exist_ok=True)
  os.chdir(movies_directory)
  if isSubset:
    print(f'==>        Making missing movies and images')
    os.system(f'{bc_makeMovie} -slices {" ".join(movie_ind)}')
  else:
    print(f'==>        Making all movies and images')
    os.system(bc_makeMovie)

  images_directory = os.path.join(args.output_dir, "map_plotting", "images")
  os.makedirs(images_directory, exist_ok=True)

  if movie_ind[0] == 'individual':
    images_sub_directory = os.path.join(images_directory, 'r000001')
    os.makedirs(images_sub_directory, exist_ok=True)
    os.system(f'mv {os.path.join(args.outpath, "plots", f"*.png")} {images_sub_directory}')
  else:
    for idx in movie_ind:
      r = f'r{str(int(idx)).zfill(6)}'
      images_sub_directory = os.path.join(images_directory, r)
      os.makedirs(images_sub_directory, exist_ok=True)
      os.system(f'mv {os.path.join(args.outpath, "plots", f"*{r}*.png")} {images_sub_directory}')
  shutil.rmtree(os.path.join(args.outpath, "plots"))

  os.chdir(args.output_dir)


def main():
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()
  
  
  
# NOTES
# More and pretty verbosity on steps it is doing.
# OFT Post Processing
# ==> Loading.....
# ==> Plotting butterfly diagram using <N> processes...
# ==> Done!  Results are in [DIRECTORY]

# Improvements TODO
#   encode movies to take up less space??
#   script option to keep images or not
