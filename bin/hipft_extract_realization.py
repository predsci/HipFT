#!/usr/bin/env python3
import psihdf as ps
import numpy as np
import argparse
from pathlib import Path

def argParsing():

    parser = argparse.ArgumentParser(description='Extract realization slices from an output 3D hipft map file and write it out in 2D files.  If no slice indecies are specified, it extracts all of them.')

    parser.add_argument(
            'hipft_3d_file',
            help='Name of 3D hipft output map file.')

    parser.add_argument('-r',
            help='Realization slices to extract (1-indexed)',
            dest='rlist',
            type=str,
            required=False)

    parser.add_argument('-o',
           help='Output file name (when using 1 realization)',
           dest='ofile',
           type=str,
           required=False)
            
    args = parser.parse_args()

    return parser.parse_args()

def main():

    ## Get input agruments:
    args = argParsing()

    pvec, tvec, rvec, data = ps.rdhdf_3d(args.hipft_3d_file)

    NR = rvec.size
    NT = tvec.size
    NP = pvec.size

    arg_dict = vars(args)

    if (args.rlist is not None):
        rlist = np.array(arg_dict['rlist'].split(','))
        rlist = rlist.astype(float)
    else:
        rlist = np.array(rvec)

    for i in rlist:
        j = int(i)
        idx = format(j, '#06d')
        if (args.ofile is not None):
            fname=args.ofile
        else:
            fname=Path(args.hipft_3d_file).stem
            fname=fname+'_r'+idx+'.h5'
        ps.wrhdf_2d(fname,pvec,tvec,data[j-1,:,:])

    

if __name__ == '__main__':
    main()


