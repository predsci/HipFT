#!/usr/bin/env python3
import argparse
import psihdf as ps
import numpy as np
from datetime import datetime

def argParsing():
  parser = argparse.ArgumentParser(description='make_radial_average_over_time:  This tool converts each 2D h5 file into a 1D avergaed latitude array and creates a new h5 file that is a time average in the latitude direction.')

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

  parser.add_argument('-3d',
                      help='Flag to denote 3d h5 file.',
                      dest='dim3',
                      default=False,
                      required=False)

  parser.add_argument('-slice',
                      help='Slice to use for averaging if 3d.',
                      dest='slice',
                      default=1,
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
  new_xvec=[]
  new_data=np.zeros((512, len(range(int(args.t0),int(args.tf)+1))))
  i = 0

  iRange = list(range(int(args.t0),int(args.tf)+1))
  #thin out range
  skip = int(args.cadence)+1
  iRange = iRange[::skip]

  for idx in iRange:
    filename=args.folder+"/"+args.bfile+"_idx"+"{:06d}".format(idx)+".h5"
    # Read the data:
    if bool(args.dim3):
        xvec, yvec, zvec, data_in = ps.rdhdf_3d(filename)
        data = np.squeeze(data_in[int(args.slice)-1,:,:])
    else:
        xvec, yvec, data = ps.rdhdf_2d(filename)
    ave_data = np.average(data, axis=1)
    new_xvec.append(idx)
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
            time.append(int(float(line.split()[1])*3600))
  time = np.array(time)

  time = time[::skip]

  aveOut = np.zeros((512, len(iRange)))
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

    aveOut[:,i] = np.average(np.take(new_data, aveIds, axis=1), axis=1)
  if args.utstart == "00000000T00:00:00":
    sDateSec = 0
  else:
    sDate = datetime.strptime(args.utstart, '%Y%m%dT%H:%M:%S')
    sDateSec = int(sdate.timestamp())

  uttime = time+sDateSec

  oFile = args.oFile+"_al"+str(args.al)+".h5"

  ps.wrhdf_2d(oFile,uttime,yvec,aveOut)


def main():
  ## Get input agruments:
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()
