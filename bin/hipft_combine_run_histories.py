#!/usr/bin/env python3

import os
import shutil
import argparse

def argParsing():
  parser = argparse.ArgumentParser(description='HipFt History Plots.')

  parser.add_argument('-odir',
    help='Output directory (Default: current directory).',
    dest='odir',
    type=str,
    required=False)


  parser.add_argument('-rundirs',
    help='A comma separated list of run directories',
    type=str,
    required=False,
    default=' ')

  return parser.parse_args()


def run(args):
  arg_dict = vars(args)
  rundir_list = arg_dict['rundirs'].split(',')


  if not (args.odir):
    dirname = os.getcwd()
  else:
    dirname = args.odir
  print(dirname)

  f=open("hipft_combined_reference.txt", "w")
  f.write("New file  ,  Original file \n")

  r=1
  for dire in rundir_list:
    for file in sorted(os.listdir(dire)):
      if "hipft_history_sol_r" in file and file.endswith(".out"):
        realization="%06d" % r
        ifile=dire+'/'+file
        nfile=dirname+'/hipft_history_sol_r'+realization+'.out'
        print('Copied '+ifile)
        shutil.copy2(ifile,nfile)
        f.write('hipft_history_sol_r'+realization+'.out'+'  ,  '+ifile+'\n')
        r=r+1

  f.close()

def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()
