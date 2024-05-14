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

# Version 2.0.1

def argParsing():
  parser = argparse.ArgumentParser(description='hipft_add_dates_to_map_output_list:  This tool converts a hipft map output text file to one with dates based on a chosen start datetime (YYYY-MM-DDTHH:MM:SS).  It outputs both a file using UTC and one using TAI.')

  parser.add_argument('utstart',
                      help='Start date of first HipFT dump (regardless of HipFT''s start_time) in UTC: YYYY-MM-DDTHH:MM:SS')

  parser.add_argument('-maplist',
                      default='hipft_output_map_list.out',
                      required=False,
                      help='Name of the HipFT output map list text file')

  parser.add_argument('-o',
                      dest='outfile',
                      help='Choose base name of output map text file.',
                      default='hipft_output_map_list',
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

  out_times_utc = []
  out_times_tai = []
  for time in times:
     out_times_utc.append(start_date_utc + TimeDelta(time*3600.0, format='sec') )
     out_times_tai.append(start_date_tai + TimeDelta(time*3600.0, format='sec') )

  # Create new map file with time columns added:
  outfile_utc = open(str(args.outfile)+'_utc.out','w')
  outfile_tai = open(str(args.outfile)+'_tai.out','w')

  i = 0
  with open(mapfile, 'r') as infile:
    FirstLine=True
    for line in infile:
        if FirstLine:
            outfile_tai.write(line[:-1] + ' TAI(sec) TAI\n')
            outfile_utc.write(line[:-1] + ' UTC(sec) UTC\n')
            FirstLine=False
        else:
            outfile_utc.write(line[:-1] + ' '+\
                          "{:.3f}".format(24.0*3600.0*out_times_utc[i].to_value('jd'))+' ' +\
                          str(out_times_utc[i])+'\n')
            outfile_tai.write(line[:-1] + ' '+\
                          "{:.3f}".format(24.0*3600.0*out_times_tai[i].to_value('jd'))+' ' +\
                          str(out_times_tai[i])+'\n')
            i = i + 1

  outfile_utc.close()
  outfile_tai.close()

def main():
  ## Get input agruments:
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()

