#!/usr/bin/env python3
import argparse
import signal
import sys
from datetime import datetime
import numpy as np


def signal_handler(signal, frame):
  print('You pressed Ctrl+C! Stopping!')
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)


def argParsing():
  parser = argparse.ArgumentParser(description='plot_2d:  This tool saves a png plot of a 2D hdf/h5 file.')

  parser.add_argument('-maplist',
                      dest='maplist',
                      help='Name of map text file',
                      required=True)

  parser.add_argument('-utstart',
                      help='Start Date in UT: YYYYMMDDTHH:MM:SS',
                      dest='utstart',
                      default='00000000T00:00:00',
                      required=True)


  return parser.parse_args()

def run(args):

  mapfile=args.maplist

  time = []
  with open(mapfile, 'r') as infile:
    FirstLine=True
    for line in infile:
        if FirstLine:
            FirstLine=False
        else:
            time.append(int(float(line.split()[1])*3600))

  time = np.array(time)

  if args.utstart == "00000000T00:00:00":
    sDateSec = 0
  else:
    sDate = datetime.strptime(args.utstart, '%Y%m%dT%H:%M:%S')
    sDateSec = int(sDate.timestamp())

  uttime = time + sDateSec

  # Create map file with ut time column appended to end:

  outfile = open('map_list_with_ut.txt','w')
  outfile = open('map_list_with_ut.txt','a')
  i = 0

  with open(mapfile, 'r') as infile:
    FirstLine=True
    for line in infile:
        if FirstLine:
            outfile.write(line[:-1] + ' UT(Sec) UT(Str)\n')
            FirstLine=False
        else:
            outfile.write(line[:-1] + ' '+str(uttime[i])+ ' ' + datetime.strftime(datetime.fromtimestamp(uttime[i]),'%Y%m%dT%H:%M:%S')+'\n')
            i = i+1


def main():
  ## Get input agruments:
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()

