#!/usr/bin/env python3
import argparse
import psihdf as ps
import numpy as np
from datetime import datetime, timezone

def argParsing():
  parser = argparse.ArgumentParser(description='hipft_make_butterfly_diagram:  This tool makes the data for a butterfly diagram by averaging each hipft map along longitude, and then taking a running average of the results over the files.  The result is a single 2D file where each column is the running average of logitudinal averages.  By specifying a UT start time, the x-axis scales is set correctly for easy plotting.')

  parser.add_argument('-folder',
    help='Name of folder.',
    dest='folder',
    required=True)

  parser.add_argument('-bfile',
    help='Base file name.',
    dest='bfile',
    required=True)

  parser.add_argument('-utstart',
    help='Start Date in UT: YYYYMMDDTHH:MM:SS',
    dest='utstart',
    default='00000000T00:00:00',
    required=False)

  parser.add_argument('-mapfile',
    help='File path to the hipft_output_map_list file.',
    dest='mapfile',
    required=False)

  parser.add_argument('-t0',
    help='Sequence start index.',
    dest='t0',
    required=True)

  parser.add_argument('-tf',
    help='Sequence stop index.',
    dest='tf',
    required=True)

  parser.add_argument('-o',
    help='Name of output h5 file.',
    dest='oFile',
    required=True)

  parser.add_argument('-tfac',
    help='Time factor to get time into units of hours (Default is days so tfac=1/24.0).',
    dest='tfac',
    default='1/24',
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

  return parser.parse_args()

def run(args):

  ###### TIME ######

  if '/' in args.tfac:
    num, denom = args.tfac.split('/')
    tfac = float(num) / float(denom)
  else:  
    tfac = float(args.tfac)

  ######################

  i = 0

  iRange = list(range(int(args.t0),int(args.tf)+1))
  #thin out range
  skip = int(args.cadence)+1
  iRange = iRange[::skip]

  firstPass = True

  for idx in iRange:
    filename=args.folder+"/"+args.bfile+"_idx"+"{:06d}".format(idx)+".h5"
    # Read the data:
    if bool(args.dim3):
      if (args.sall):
        xvec, yvec, zvec, data = ps.rdhdf_3d(filename)
        dim3 = len(zvec)
        if firstPass:
          new_data=np.zeros((dim3,512, len(range(int(args.t0),int(args.tf)+1))))
          firstPass = False
        ave_data = np.average(data, axis=2)
        new_data[:,:,i] = ave_data
        i+=1
      elif (args.slices):
        dim3 = len(args.slices)
        if firstPass:
          new_data=np.zeros((dim3,512, len(range(int(args.t0),int(args.tf)+1))))
          firstPass = False
        j = 0
        xvec, yvec, zvec, data_in = ps.rdhdf_3d(filename)
        for islice in args.slices:
          data = np.squeeze(data_in[islice-1,:,:])
          ave_data = np.average(data, axis=1)
          new_data[j,:,i] = ave_data
          j+=1
        i+=1
      else:
        if firstPass:
          new_data=np.zeros((512, len(range(int(args.t0),int(args.tf)+1))))
          firstPass = False
        xvec, yvec, zvec, data_in = ps.rdhdf_3d(filename)
        data = np.squeeze(data_in[int(args.slice)-1,:,:])
        ave_data = np.average(data, axis=1)
        new_data[:,i] = ave_data
        i+=1
    else:
      if firstPass:
        new_data=np.zeros((512, len(range(int(args.t0),int(args.tf)+1))))
        firstPass = False
      xvec, yvec, data = ps.rdhdf_2d(filename)
      ave_data = np.average(data, axis=1)
      new_data[:,i] = ave_data
      i+=1

  if args.mapfile is None:
      mapfile = args.folder+"/hipft_output_map_list.out"
  else:
      mapfile=args.mapfile

  time = []
  with open(mapfile, 'r') as infile:
    FirstLine=True
    for line in infile:
        if FirstLine:
            FirstLine=False
        else:
            time.append(int(float(line.split()[1]))*3600*tfac)
  time = np.array(time)

  time = time[::skip]

  liRange=len(iRange)

  if args.utstart == "00000000T00:00:00":
    sDateSec = 0
  else:
    sDate = datetime.strptime(args.utstart, '%Y%m%dT%H:%M:%S')
    sDateSec = int(sDate.replace(tzinfo=timezone.utc).timestamp())
  uttime = time+sDateSec

  if (args.sall):
    writeOut3d(args,liRange,time,zvec,new_data,uttime,yvec,dim3)
  elif (args.slices):
    writeOut3d(args,liRange,time,args.slices,new_data,uttime,yvec,dim3)
  else:
    writeOut2d(args,liRange,time,new_data,uttime,yvec)


def writeOut3d(args,liRange,time,zvec,newData,uttime,yvec,dim3):
  aveOut = np.zeros((dim3,512, liRange))
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
  aveOut = np.zeros((512, liRange))
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
