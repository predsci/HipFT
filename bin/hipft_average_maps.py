#!/usr/bin/env python3
import argparse
import psihdf as ps
import numpy as np
from astropy.time import Time
import os

def argParsing():
  parser = argparse.ArgumentParser(description='hipft_average_maps:  This tool makes and average map based on a sequence of hipft maps.')

  parser.add_argument('-rundir',
    help='Path to the directory where hipft was run.  Default is current directory.',
    dest='rundir',
    required=True)

  parser.add_argument('-i0',
    help='Sequence start index.',
    dest='i0',
    required=True)

  parser.add_argument('-i1',
    help='Sequence stop index.',
    dest='i1',
    required=True)

  parser.add_argument('-o',
    help='Output map name.',
    dest='outfile',
    default='hipft_brmap_average.h5',
    required=False)


  return parser.parse_args()

def run(args):

  #Check that rundir exists.
  if not os.path.exists(args.rundir):
    print('ERROR!  Rundir not found:  '+args.rundir)
    exit(1)

  index_range = list(range(int(args.i0),int(args.i1)+1))
  num_files = len(index_range)

  firstPass = True

  for idx in index_range:

    filename = args.rundir+"/hipft_brmap_idx"+"{:06d}".format(idx)+".h5"

    if firstPass:
      xvec, yvec, data = ps.rdhdf_2d(filename)
      avg_data = data.copy()
      firstPass = False
    else:
      data = []
      _, _, data = ps.rdhdf_2d(filename)
      
      avg_data[:,:] = avg_data[:,:] + data[:,:]

    
  avg_data[:,:] = avg_data[:,:] / np.float64(num_files)

  ps.wrhdf_2d(args.outfile, xvec, yvec, avg_data)


def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()
