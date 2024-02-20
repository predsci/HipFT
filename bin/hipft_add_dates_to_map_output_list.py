#!/usr/bin/env python3
import argparse
import signal
import sys
#from datetime import datetime
import numpy as np
from astropy.time import Time
from astropy.time import TimeDelta

def signal_handler(signal, frame):
  print('You pressed Ctrl+C! Stopping!')
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

# Version 2.0

def argParsing():
  parser = argparse.ArgumentParser(description='hipft_add_dates_to_map_output_list:  This tool converts a hipft map output text file to one with dates based on a chosen start datetime (YYYY-MM-DDTHH:MM:SS).')

  parser.add_argument('utstart',
                      help='Start Date of first HipFT dump (regardless of hipft code time)in UTC: YYYY-MM-DDTHH:MM:SS')

  parser.add_argument('-maplist',
                      default='hipft_output_map_list.out',
                      required=False,
                      help='Name of the HipFT output map list text file')

  parser.add_argument('-outfmt',
                      required=False,
                      default='tai',
                      help='Output format for dates.  Default is TAI.  Options are "tai" or "utc".')

  parser.add_argument('-o',
                      dest='outfile',
                      help='Choose name of output map text file.',
                      default='hipft_output_map_list_tai.out',
                      required=False)

  return parser.parse_args()

def run(args):

  mapfile=args.maplist

  times = []
  with open(mapfile, 'r') as infile:
    FirstLine=True
    for line in infile:
        if FirstLine:
            FirstLine=False
        else:
            times.append(float(line.split()[1]))

  times = np.asarray(times)

  #Reset to start at 0.
  tstart = times[0]
  times = times - tstart

  #Assume start date is in UTC and store it.
  start_date_utc = Time(args.utstart, format='isot', scale='utc', out_subfmt='date_hms')
  start_date_tai = start_date_utc.tai

  if args.outfmt == "utc":
      start_date = start_date_utc
  else:
      start_date = start_date_tai

  out_times = []
  for time in times:
     out_times.append(start_date + TimeDelta(time*3600.0, format='sec') )

  # Create new map file with time columns added:
  outfile = open(args.outfile,'w')

  i = 0
  with open(mapfile, 'r') as infile:
    FirstLine=True
    for line in infile:
        if FirstLine:
            if args.outfmt == 'tai':
                outfile.write(line[:-1] + ' TAI(sec) TAI(str)\n')
            else:
                outfile.write(line[:-1] + ' UTC(sec) UTC(str)\n')
            FirstLine=False
        else:
            outfile.write(line[:-1] + ' '+\
                          "{:.3f}".format(24.0*3600.0*out_times[i].to_value('jd'))+' ' +\
                          str(out_times[i])+'\n')
            i = i + 1


def main():
  ## Get input agruments:
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()

