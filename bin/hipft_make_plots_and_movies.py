#!/usr/bin/env python3
import argparse
import os
import h5py
import re
import signal
import sys
from pathlib import Path
import concurrent.futures
import glob

# Version 2.4.3

def handle_int(signum, frame):
    print('You pressed Ctrl+C! Stopping!')
    signal.signal(signum, signal.SIG_DFL)
    os.kill(os.getpid(), signum)
    sys.exit(1)
    
signal.signal(signal.SIGINT, handle_int)

def argParsing():
  parser = argparse.ArgumentParser(description='HipFT Movie Maker.')

  parser.add_argument('-dir',
    help='Directory of output files (full path!)',
    dest='datadir',
    type=str,
    required=True)

  parser.add_argument('-dpi',
    help='DPI for plots.',
    dest='dpi',
    default=128,
    type=int,
    required=False)

  parser.add_argument('-odir',
    help='Output directory of movie (Default: current directory).',
    dest='odir',
    type=str,
    required=False)

  parser.add_argument('-o',
    help='Base filename for output mov movies.',
    dest='outfile',
    default='hipft_movie',
    type=str,
    required=False)

  parser.add_argument('-slice',
    help='Slice index (if not specified all slices wil be plotted).',
    dest='s',
    type=int,
    required=False)
  
  parser.add_argument('-slices',
    help='Slice indices (if not specified all slices wil be plotted)',
    dest='slices', 
    type=int, 
    nargs='+',
    required=False)

  parser.add_argument('-cmin',
    help='Colormap Minimum',
    dest='cmin',
    default='-25',
    type=float,
    required=False)

  parser.add_argument('-cmax',
    help='Colormap Maximum',
    dest='cmax',
    default='25',
    type=float,
    required=False)

  parser.add_argument('-label',
    help='Colorbar Label',
    dest='label',
    default='Gauss',
    type=str,
    required=False)

  parser.add_argument('-mldfile',
    default='hipft_output_map_list_tai.out',
    required=False,
    help='Name of the HipFT map text file with dates.')

  parser.add_argument('-nomovie',
    help='Do not make mov movie.',
    dest='nomovie',
    default=False,
    action='store_true',
    required=False)

  parser.add_argument('-np',
    help='Number of threads to use (Default 1 or the environment variable OMP_NUM_THREADS if set).',
    dest='np',
    type=int,
    required=False)

  return parser.parse_args()


def run(args):

  if not args.np:
    args.np = int(os.getenv('OMP_NUM_THREADS', 1))

  if args.odir:
    odir=str(Path(args.odir).resolve())
  else:
    odir=os.getcwd()

  args.datadir = str(Path(args.datadir).resolve())

  title_str=[]

  date_fmt = ''

  h5_files_tmp = [file for file in os.listdir(args.datadir) if file.endswith('.h5')]

  sorted_listdir = sorted(h5_files_tmp)

  if len(sorted_listdir) < 1:
    print('No .h5 files found.')
    exit()

  if os.path.exists(args.mldfile):
    with open(args.mldfile, "r") as ftmp:
      next(ftmp)
      for line in ftmp:
        title_str.append(str(line.split()[4]))
    with open(args.mldfile, "r") as ftmp:
        first_line = ftmp.readline().strip('\n').split()
    date_fmt=first_line[4]+': '
  elif os.path.exists('hipft_output_map_list.out'):
    with open('hipft_output_map_list.out', "r") as ftmp:
      next(ftmp)
      for line in ftmp:
        title_str.append(str(float(line.split()[1]))+' hours')
  else:
    idx=0
    for filetmp in sorted_listdir:
      idx+=1
      title_str.append("File number: {:06d}".format(idx))

  if not os.path.exists(args.datadir+'/plots'):
    os.makedirs(args.datadir+'/plots')

  os.chdir(args.datadir+'/plots')

  def process_file(filetmp, idx):
    dim3 = None
    file=args.datadir+'/'+filetmp
    if filetmp.endswith('.h5'):
      re_idx = re.search("idx[0-9]+",filetmp)
      if not re_idx:
        print("Expecting the file "+ filetmp  +" to contain an idx[0-9]+ number.")
        return dim3
      idxx = int(re.search("[0-9]+",re.search("idx[0-9]+",filetmp).group()).group())
      TITLE = date_fmt+title_str[idx-1]
      with h5py.File(file,'r') as f1:
        twoD = True
        if 'dim3' in list(f1.keys()):
          dim3 = f1['dim3'].shape[0]
          twoD = False
      if twoD:
        fileOut = os.path.basename(file).replace('.h5','.png')
        plot2d(TITLE,args,file,fileOut)
      elif (args.s):
        if (args.s > dim3):
          print("Slice requested is outside range defaulting to last slice : "+ str(dim3))
          extractANDplot(TITLE,args,file, dim3, idx)
        else:
          extractANDplot(TITLE,args,file, args.s, idx)
      elif (args.slices):
        for slice in args.slices:
          if (slice > dim3):
            print(f"Slice {slice} requested is outside range defaulting to last slice : {dim3}")
            extractANDplot(TITLE, args, file, dim3, idx)
        else:
          extractANDplot(TITLE, args, file, slice, idx)
      else:
        for i in range(1,dim3+1):
          extractANDplot(TITLE, args, file, i, idx)
      return dim3
    return dim3

  with concurrent.futures.ThreadPoolExecutor(max_workers=args.np) as executor:
    futures = [executor.submit(process_file, filetmp, idx + 1) for idx, filetmp in enumerate(sorted_listdir)]

  twoD = True
  for future in concurrent.futures.as_completed(futures):
    dim3 = future.result()
    if dim3 is not None:
        twoD = False 

  for file in glob.glob("tmp_file*.h5"):
    os.remove(file)

  if not os.path.exists(args.datadir+'/plots/tmp'):
    os.makedirs(args.datadir+'/plots/tmp')

  os.chdir(args.datadir+'/plots/tmp')

  if not args.nomovie:
    if twoD is None:
      print('No files found to make movie from.')
      exit()
    elif twoD:
      for filetmp in sorted(os.listdir(args.datadir+'/plots')):
        file=args.datadir+'/plots/'+filetmp
        if filetmp.endswith('.png'):
          linkFile(file,filetmp)
      ffmpefMovie(args,odir,args.outfile)
    elif (args.s):
      if (args.s > dim3):
        makeMovie(args,odir,dim3)
      else:
        makeMovie(args,odir,args.s)
    elif (args.slices):
      for slice in args.slices:
        if (slice > dim3):
          makeMovie(args,odir,dim3)
        else:
          makeMovie(args,odir,slice)
    else:
      for i in range(1,dim3+1):
        makeMovie(args,odir,i)

  for filetmp in sorted(os.listdir(args.datadir+'/plots/tmp')):
    os.remove(filetmp)

  os.chdir(args.datadir+'/plots')
  os.rmdir(args.datadir+'/plots/tmp')
  for file in glob.glob("tmp_file*.h5"):
    os.remove(file)


def makeMovie(args,odir,r):
  fslice="%06d" % r
  for filetmp in sorted(os.listdir(args.datadir+'/plots')):
    file=args.datadir+'/plots/'+filetmp
    if filetmp.endswith('_r'+ fslice +'.png'):
      linkFile(file,filetmp)
  name=args.outfile+'_r'+fslice
  ffmpefMovie(args,odir,name)


def linkFile(file,filetmp):
  idxx=int(re.search("[0-9]+",re.search("idx[0-9]+",filetmp).group()).group())
  os.system('ln -sf '+file+' movie'+ str(idxx) +'.png')


def ffmpefMovie(args,odir,name):
  if os.path.dirname(os.path.normpath(os.popen('which ffmpeg').read().strip())) == '':
    print("Warning - ffmpeg not detected, not making movie file.")
  else:
    os.system('ffmpeg -hide_banner -loglevel error -y -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "'+name+'.mov"')
    os.system('cp '+args.datadir+'/plots/tmp/'+name+'.mov '+odir+'/'+name+'.mov')


def extractANDplot(TITLE, args, file, r, idx):
  tmp_file = f'tmp_file{idx}.h5'
  os.system(f'hipft_extract_realization.py {file} -r {r} -o {tmp_file}')
  fileOut = os.path.basename(file).replace('.h5', '') + f'_r{r:06d}.png'
  plot2d(TITLE, args, tmp_file, fileOut)


def plot2d(TITLE,args,file,fileOut):
  os.system('plot2d -title "'+ TITLE +'" -cmin '+ str(args.cmin) +' -cmax '+ str(args.cmax) +' -dpi ' + str(args.dpi) \
    +' -tp -ll -finegrid -unit_label '+ args.label + " "+ file +' -o "'+ str(fileOut) +'"')


def main():
  ## Get input agruments:
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()
