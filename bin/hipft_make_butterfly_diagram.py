#!/usr/bin/env python3
import argparse
import psihdf as ps
import numpy as np
from astropy.time import Time

# Version 1.1.0

def argParsing():
  parser = argparse.ArgumentParser(description='hipft_make_butterfly_diagram:  This tool makes the data for a butterfly diagram by averaging each hipft map along longitude, and then taking a running average of the results over the files.  The result is a single 2D file where each column is the running average of logitudinal averages.  By specifying a UT start time, the x-axis scales is set correctly for easy plotting.')

  parser.add_argument('-rundir',
    help='Path to the directory where hipft was run.  Default is current directory.',
    dest='rundir',
    required=False)

  parser.add_argument('-outpath',
    help='Path to folder where output data is. Default: args.rundir/output_maps',
    dest='outpath',
    required=False)
#   Could grep for output_map_directory in hipft.in, if not there,
#   use default.

  parser.add_argument('-basefn',
    help='Base file name for output maps.',
    dest='bfile',
    default='hipft_brmap',
    required=False)
#   Could grep for output_map_root_filename in hipft.in, and use default otherwise.

  parser.add_argument('-utstart',
    help='Start Date in UT: YYYY-MM-DDTHH:MM:SS',
    dest='utstart',
    default='0000-00-00T00:00:00',
    required=False)
#   Could extract from data assimilation csv file combined with start index...  too much?    

  parser.add_argument('-maplist',
    help='Full path to the map file list.  Default: args.rundir/hipft_output_map_list.out',
    dest='mapfile',
    required=False)

  parser.add_argument('-t0',
    help='Sequence start index.',
    dest='t0',
    default=1,
    required=False)

  parser.add_argument('-tf',
    help='Sequence stop index. If not specified, all maps in mapfile are used.',
    dest='tf',
    required=False)

  parser.add_argument('-o',
    help='Name of output h5 file.',
    dest='oFile',
    default='butterfly.h5',
    required=False)

  parser.add_argument('-tfac',
    help='Time factor to get time into units of hours (HipFT time is in hours, so default tfac=1).',
    dest='tfac',
    default='1',
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

  parser.add_argument('-al',
    help='Time length to average over in hours.',
    dest='al',
    default=0,
    required=False)

  parser.add_argument('-c',
    help='Cadence of how many files to skip.',
    dest='cadence',
    default=0,
    required=False)

  parser.add_argument('-tai',
    help='Use tai which ignores leap seconds.',
    dest='tai',
    action='store_true',
    default=False,
    required=False)

  return parser.parse_args()

def run(args):

  if (args.rundir is None):
    args.rundir = os.getcwd()
  if (args.outpath is None):
    args.outpath = args.rundir+'/output_maps'
  if (args.mapfile is None):
    args.mapfile = args.rundir+'/hipft_output_map_list.out'
  
  #Check that mapfile exists.
  if not os.path.exists(args.mapfile):
    print('ERROR!  Map list file not found:  '+args.mapfile)
    exit(1)
  
  ###### TIME ######

  if '/' in args.tfac:
    num, denom = args.tfac.split('/')
    tfac = float(num) / float(denom)
  else:  
    tfac = float(args.tfac)

  ######################

  i = 0
  time = []
  count = 0
  with open(args.mapfile, 'r') as infile:
    FirstLine=True
    for line in infile:
        if FirstLine:
            FirstLine=False
        else:
            time.append(int(float(line.split()[1]))*3600*tfac)
            count = count + 1
  time = np.array(time)
  skip = int(args.cadence)+1
  time = time[::skip]

  if (args.tf is None):
    args.tf = count

  iRange = list(range(int(args.t0),int(args.tf)+1))
  
  #thin out range
  iRange = iRange[::skip]
  
  liRange=len(iRange)

  firstPass = True

  for idx in iRange:
    filename=args.outpath+"/"+args.bfile+"_idx"+"{:06d}".format(idx)+".h5"
    # Read the data:
    if bool(args.dim3):
      if (args.sall):
        xvec, yvec, zvec, data = ps.rdhdf_3d(filename)
        dim3 = len(zvec)
        if firstPass:
          new_data=np.zeros((dim3,len(yvec), len(range(int(args.t0),int(args.tf)+1))))
          firstPass = False
        ave_data = np.average(data, axis=2)
        new_data[:,:,i] = ave_data
        i+=1
      elif (args.slices):
        xvec, yvec, zvec, data_in = ps.rdhdf_3d(filename)
        dim3 = len(args.slices)
        if firstPass:
          new_data=np.zeros((dim3,len(yvec), len(range(int(args.t0),int(args.tf)+1))))
          firstPass = False
        j = 0
        for islice in args.slices:
          data = np.squeeze(data_in[islice-1,:,:])
          ave_data = np.average(data, axis=1)
          new_data[j,:,i] = ave_data
          j+=1
        i+=1
      else:
        xvec, yvec, zvec, data_in = ps.rdhdf_3d(filename)
        if firstPass:
          new_data=np.zeros((len(yvec), len(range(int(args.t0),int(args.tf)+1))))
          firstPass = False
        data = np.squeeze(data_in[int(args.slice)-1,:,:])
        ave_data = np.average(data, axis=1)
        new_data[:,i] = ave_data
        i+=1
    else:
      xvec, yvec, data = ps.rdhdf_2d(filename)
      if firstPass:
        new_data=np.zeros((len(yvec), len(range(int(args.t0),int(args.tf)+1))))
        firstPass = False
      ave_data = np.average(data, axis=1)
      new_data[:,i] = ave_data
      i+=1

  if args.utstart == "0000-00-00T00:00:00":
    sDateSec = 0
  else:
    sDateSec = getSec(args.tai,args.utstart,'%Y-%m-%dT%H:%M:%S')
  uttime = time+sDateSec

  if (args.sall):
    writeOut3d(args,liRange,time,zvec,new_data,uttime,yvec,dim3)
  elif (args.slices):
    writeOut3d(args,liRange,time,args.slices,new_data,uttime,yvec,dim3)
  else:
    writeOut2d(args,liRange,time,new_data,uttime,yvec)

def getSec(atai,x,xformat):
  if atai:
    t = Time.strptime(x,xformat,scale='tai')
    currDateSeconds = t.unix_tai
  else:
    t = Time.strptime(x,xformat,scale='utc')
    currDateSeconds = t.unix
  return currDateSeconds

def writeOut3d(args,liRange,time,zvec,newData,uttime,yvec,dim3):
  aveOut = np.zeros((dim3,len(yvec), liRange))
  aveIds = []
  for i in range(len(time)):
    for j in range (i, len(time)):
        if time[j]-time[i] <= (int(args.al)*3600):
            if j not in aveIds:
                aveIds.append(j)
        else:
            break

    while len(aveIds) > 0:
        if abs(time[aveIds[0]]-time[i]) > (int(args.al)*3600):
            del aveIds[0]
        else:
            break

    aveOut[:,:,i] = np.average(np.take(newData, aveIds, axis=2), axis=2)

  ps.wrhdf_3d(args.oFile,uttime,yvec,zvec,aveOut)


def writeOut2d(args,liRange,time,newData,uttime,yvec):
  aveOut = np.zeros((len(yvec), liRange))
  aveIds = []
  for i in range(len(time)):
    for j in range (i, len(time)):
        if time[j]-time[i] <= (int(args.al)*3600):
            if j not in aveIds:
                aveIds.append(j)
        else:
            break

    while len(aveIds) > 0:
        if abs(time[aveIds[0]]-time[i]) > (int(args.al)*3600):
            del aveIds[0]
        else:
            break

    aveOut[:,i] = np.average(np.take(newData, aveIds, axis=1), axis=1)

  ps.wrhdf_2d(args.oFile,uttime,yvec,aveOut)


def main():
  ## Get input agruments:
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()
