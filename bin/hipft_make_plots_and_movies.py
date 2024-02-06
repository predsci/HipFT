#!/usr/bin/env python3
import argparse
import os
import h5py

def argParsing():
  parser = argparse.ArgumentParser(description='HipFt Movie Maker.')

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

  parser.add_argument('-mlut',
    default='hipft_output_map_list_ut.out',
    required=False,
    help='Name of the HipFT map text file with ut dates.')

  return parser.parse_args()


def run(args):
  if args.odir:
    odir=args.odir
  else:
    odir=os.getcwd()

  ut_dates=[]

  if os.path.exists(args.mlut):
    with open(args.mlut, "r") as ftmp:
      next(ftmp)
      for line in ftmp:
        ut_dates.append(line.split()[4])
  elif os.path.exists('hipft_output_map_list.out'):
    with open('hipft_output_map_list.out', "r") as ftmp:
      next(ftmp)
      for line in ftmp:
        ut_dates.append(str(float(line.split()[1]))+' hours')

  if not os.path.exists(args.datadir+'/plots'):
    os.makedirs(args.datadir+'/plots')

  os.chdir(args.datadir+'/plots')

  idx=0
  for filetmp in os.listdir(args.datadir):
    file=args.datadir+'/'+filetmp
    if filetmp.endswith('.h5'):
      idxx=int(filetmp[-9:-3])

      TITLE=ut_dates[idxx-1]
      idx+=1
      with h5py.File(file,'r') as f1:
        twoD=True
        if 'dim3' in list(f1.keys()):
          dim3 = f1['dim3'].shape[0]
          twoD=False

      if twoD:
        os.system('plot2d -title "'+ TITLE +'" -cmin '+ str(args.cmin) +' -cmax '+ str(args.cmax) +' -dpi ' + str(args.dpi) \
          +' -tp -ll -finegrid -unit_label '+ args.label + " "+ file +' -o "'+ os.path.basename(file).replace('.h5','.png'+'"'))
      elif (args.s == "all"):  
        for i in range(1,dim3):
          os.system('hipft_extract_realization.py '+ file +' -r '+ i +' -o tmp_file.h5')
          fslice="%06d" % (i)
          os.system('plot2d -title "'+ TITLE +'" -cmin '+ str(args.cmin) +' -cmax '+ str(args.cmax) +' -dpi ' + str(args.dpi) \
            +' -tp -ll -finegrid -unit_label '+ args.label +' tmp_file.h5 -o "'+ os.path.basename(file).replace('.h5','') +'_r' + str(fslice) +'.png"')
      elif (args.s):
        if (args.s > dim3):
          print("Slice requested is outside range defaulting to last slice : "+ str(dim3))
          os.system('hipft_extract_realization.py '+ file +' -r '+ str(dim3) +' -o tmp_file.h5')
          fslice="%06d" % (dim3)
          os.system('plot2d -title "'+ TITLE +'" -cmin '+ str(args.cmin) +' -cmax '+ str(args.cmax) +' -dpi ' + str(args.dpi) \
            +' -tp -ll -finegrid -unit_label '+ args.label +' tmp_file.h5 -o "'+ os.path.basename(file).replace('.h5','') +'_r' + str(fslice) +'.png"')
        else:
          os.system('hipft_extract_realization.py '+ file +' -r '+ str(args.s) +' -o tmp_file.h5')
          fslice="%06d" % (args.s)
          os.system('plot2d -title "'+ TITLE +'" -cmin '+ str(args.cmin) +' -cmax '+ str(args.cmax) +' -dpi ' + str(args.dpi) \
            +' -tp -ll -finegrid -unit_label '+ args.label +' tmp_file.h5 -o "'+ os.path.basename(file).replace('.h5','') +'_r' + str(fslice) +'.png"')
      else:
        os.system('hipft_extract_realization.py '+ file +' -r 1 -o tmp_file.h5')
        os.system('plot2d -title "'+ TITLE +'" -cmin '+ str(args.cmin) +' -cmax '+ str(args.cmax) +' -dpi ' + str(args.dpi) \
            +' -tp -ll -finegrid -unit_label '+ args.label +' tmp_file.h5 -o "'+ os.path.basename(file).replace('.h5','') +'_r000001.png"')
  if os.path.exists("tmp_file.h5"):
    os.remove("tmp_file.h5")

  if not os.path.exists(args.datadir+'/plots/tmp'):
    os.makedirs(args.datadir+'/plots/tmp')

  os.chdir(args.datadir+'/plots/tmp')


  if twoD:
    for filetmp in os.listdir(args.datadir+'/plots'):
      file=args.datadir+'/plots/'+filetmp
      if filetmp.endswith('.png'):
        idxx=int(filetmp[-10:-4])
        os.system('ln -s '+file+' movie'+ str(idxx) +'.png')
    os.system('ffmpeg -y -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "movie.mov"')
    os.system('cp '+args.datadir+'/plots/tmp/movie.mov '+odir+'/'+args.outfile+'.mov')
  elif (args.s == "all"):
    for j in range(1,dim3):
      fslice="%06d" % (i)
      for filetmp in os.listdir(args.datadir+'/plots'):
        file=args.datadir+'/plots/'+filetmp
        if filetmp.endswith('_r'+ fslice +'.png'):
          idxx=int(filetmp[-10:-4])
          os.system('ln -s '+file+' movie'+ str(idxx) +'.png')
      os.system('ffmpeg -y -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "movie.mov"')
      os.system('cp '+args.datadir+'/plots/tmp/movie.mov '+odir+'/'+args.outfile+'_r'+fslice+'.mov')
      if os.path.exists("movie.mov"):
        os.remove('movie.mov')
  elif (args.s):
    if (args.s > dim3):
      fslice="%06d" % (dim3)
      for filetmp in os.listdir(args.datadir+'/plots'):
        file=args.datadir+'/plots/'+filetmp
        if filetmp.endswith('_r'+ fslice +'.png'):
          idxx=int(filetmp[-10:-4])
          os.system('ln -s '+file+' movie'+ str(idxx) +'.png')
      os.system('ffmpeg -y -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "movie.mov"')
      os.system('cp '+args.datadir+'/plots/tmp/movie.mov '+odir+'/'+args.outfile+'_r'+fslice+'.mov')
    else:
      fslice="%06d" % (args.s)
      for filetmp in os.listdir(args.datadir+'/plots'):
        file=args.datadir+'/plots/'+filetmp
        if filetmp.endswith('_r'+ fslice +'.png'):
          idxx=int(filetmp[-10:-4])
          os.system('ln -s '+file+' movie'+ str(idxx) +'.png')
      os.system('ffmpeg -y -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "movie.mov"')
      os.system('cp '+args.datadir+'/plots/tmp/movie.mov '+odir+'/'+args.outfile+'_r'+fslice+'.mov')
  else:
    for filetmp in os.listdir(args.datadir+'/plots'):
      file=args.datadir+'/plots/'+filetmp
      if filetmp.endswith('.png') and 'idx' in file:
#      if filetmp.endswith('.png') and file.contains('_r'):
        idxx=int(filetmp[-10:-4])
        os.system('ln -s '+file+' movie'+ str(idxx) +'.png')
    os.system('ffmpeg -y -framerate 15 -i "movie%d.png" -pix_fmt yuv420p -c:a copy -crf 20 -r 15 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -codec:v libx264 "movie.mov"')
    os.system('cp '+args.datadir+'/plots/tmp/movie.mov '+odir+'/'+args.outfile+'_r000001.mov')

  for filetmp in os.listdir(args.datadir+'/plots/tmp'):
    os.remove(filetmp)
  os.chdir(args.datadir+'/plots')
  os.rmdir(args.datadir+'/plots/tmp')



def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()
