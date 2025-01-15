#!/usr/bin/env python3
import psihdf as ps
import re
import numpy as np
import argparse
from pathlib import Path

# Version 1.0.1

def argParsing():

    parser = argparse.ArgumentParser(description='Pack realization slices from 2D files into a 3D hipft map file.')

    parser.add_argument(
            'hipft_2d_file_list',
            help='A comma separated list of 2d hipft files')

    parser.add_argument('-o',
           help='Output file name',
           dest='ofile',
           type=str,
           required=False)

    return parser.parse_args()

def main():

    ## Get input agruments:
    args = argParsing()

    arg_dict = vars(args)
    hipft_2d_file_list = arg_dict['hipft_2d_file_list'].split(',')

    data_arr = []
    rvec=[]
    pvec_first, tvec_first = None, None

    for file2d in hipft_2d_file_list:
        match = re.search(r'r\d{6}', file2d)
        print(match)
        if match:
            idx = int(match.group(0).replace('r',''))
        else:
            raise ValueError(f"Index not found in filename: {file2d}")


        pvec, tvec, data = ps.rdhdf_2d(file2d)

        if pvec_first is None and tvec_first is None:
            pvec_first = pvec.copy()
            tvec_first = tvec.copy()
        else:
            if not np.array_equal(tvec, tvec_first):
                print(f"This has a different tvec than the first file: {file2d}")
            if not np.array_equal(pvec, pvec_first):
                print(f"This has a different pvec than the first file: {file2d}")

        data_arr.append(data[:, :])
        rvec.append(idx)

    if (args.ofile is not None):
        fname=args.ofile
    else:
        fname = Path(file2d).stem
        fname = re.sub(r'_r\d{6}', '', fname)+'.h5'
    print(rvec)
    ps.wrhdf_3d(fname, pvec, tvec, rvec, data_arr)


if __name__ == '__main__':
    main()
