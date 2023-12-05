#!/usr/bin/env python3
import argparse
import numpy as np
import sys

# This script needs error checking!!

def argParsing():
  parser = argparse.ArgumentParser(description='HipFt Compare Run Diagnostics.')

  parser.add_argument('run1dir',
    help='Directory of first run.',
    type=str)

  parser.add_argument('run2dir',
    help='Directory of 2nd run.',
    type=str)

  parser.add_argument('-p',
    help='Set precision (in terms of decimal places) to check.',
    type=int,
    dest='p',
    required=False,
    default=6)

  parser.add_argument('-v',
    help='Verbose mode',
    dest='v',
    required=False,
    default=False,
    action='store_true')    

  return parser.parse_args()


def compare_files(file1, file2, precision, verbose):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        last_line1 = list(f1)[-1]
        last_line2 = list(f2)[-1]
    with open(file1, 'r') as f1:
        first_line1 = list(f1)[0]

# Loop over quantities we want to check:
    PASS=True
    for i in (2,3,4,5,6,15):

      num1=last_line1.split()[i]
      num1p=num1[0:precision+3]+num1[-5:]

      num2=last_line2.split()[i]
      num2p=num2[0:precision+3]+num1[-5:]

      if (num1p != num2p):
        if (float(num1) > 1e-64):
          print('FAIL in '+first_line1.split()[i]+':')
          print('        '+num1)
          print('        '+num2)
          PASS=False
        else:
          print('PASSING A FAIL in tiny number in '+first_line1.split()[i]+':')
          print('        '+num1)
          print('        '+num2)

      if (verbose):
        print(first_line1.split()[i]+':')
        print('        '+num1)
        print('        '+num2)

    return PASS


def run(args):

# Get number of realizations from input file (assume hipft.in)
  nr = 1
  with open(args.run1dir+'/hipft.in','r') as file:
    for line in file:
      words = line.split()
      if 'n_realizations' in words:
        nr = int(words[2])
        break

  STATUS = 0

# Loop over all realizations.
  for i in range(nr):
  
    if (compare_files(args.run1dir+'/hipft_history_sol_r00000'+str(i+1)+'.out', 
                      args.run2dir+'/hipft_history_sol_r00000'+str(i+1)+'.out', args.p, args.v)):
      STATUS=0
    else:
      STATUS=1

# Return STATUS to the calling shell and print it.

  print(STATUS)
  sys.exit(STATUS)

def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()
