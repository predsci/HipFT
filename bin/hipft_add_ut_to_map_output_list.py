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

# Version 2.0

def argParsing():
  parser = argparse.ArgumentParser(description='hipft_add_ut_to_map_output_list:  This tool converts a hipft map output text file to one with UT dates based on a chosen start UT time (YYYY-MM-DDTHH:MM:SS).')

  parser.add_argument('utstart',
                      help='Start Date in UT: YYYY-MM-DDTHH:MM:SS',
                      default='0000-00-00T00:00:00')

  parser.add_argument('-maplist',
                      default='hipft_output_map_list.out',
                      required=False,
                      help='Name of the HipFT output map list text file')

  parser.add_argument('-o',
                      dest='outfile',
                      help='Choose name of output map text file with ut dates.',
                      default='hipft_output_map_list_ut.out',
                      required=False)

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

  if args.utstart == "0000-00-00T00:00:00":
    sDateSec = 0
  else:
    sDate = datetime.strptime(args.utstart, '%Y-%m-%dT%H:%M:%S')
    sDateSec = int(sDate.timestamp())

  uttime = time + sDateSec

  # Create map file with ut time column appended to end:
  outfile = open(args.outfile,'w')
  outfile = open(args.outfile,'a')
  i = 0

  with open(mapfile, 'r') as infile:
    FirstLine=True
    for line in infile:
        if FirstLine:
            outfile.write(line[:-1] + ' UT(Sec) UT(Str)\n')
            FirstLine=False
        else:
            outfile.write(line[:-1] + ' '+str(uttime[i])+ ' ' + datetime.strftime(datetime.fromtimestamp(uttime[i]),'%Y-%m-%dT%H:%M:%S')+'\n')
            i = i+1


def main():
  ## Get input agruments:
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()

