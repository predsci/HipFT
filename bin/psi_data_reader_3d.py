#!/usr/bin/env python3
import psihdf as ps
import numpy as np
import argparse

def argParsing():

    parser = argparse.ArgumentParser(description='Read PSI 3D hdf data.')

    parser.add_argument("psi_3d_hdf_file_name",
        help='Name of 3D PSI HDF or H5 file (e.g. rho001.hdf or rho001.h5).')

    args = parser.parse_args()

    return parser.parse_args()

def main():

    ## Get input agruments:
    args = argParsing()

    rvec, tvec, pvec, data = ps.rdhdf_3d(args.psi_3d_hdf_file_name)

    NR = rvec.size
    NT = tvec.size
    NP = pvec.size

    print('Opened file:   '+args.psi_3d_hdf_file_name)
    print(' ')
    print('Shape of data:'+str(data.shape))
    print(' ')
    print('Nr:         '+str(NR))
    print('Ntheta:     '+str(NT))
    print('Nphi:       '+str(NP))
    print(' ')
    print('min(r): '+str(np.min(rvec)))
    print('max(r): '+str(np.max(rvec)))
    print('min(theta): '+str(np.min(tvec)))
    print('max(theta): '+str(np.max(tvec)))
    print('min(phi):   '+str(np.min(pvec)))
    print('max(phi):   '+str(np.max(pvec)))
    print(' ')
    print('min(data):  '+str(np.min(data)))
    print('max(data):  '+str(np.max(data)))
    print('mean(data): '+str(np.mean(data)))
    print(' ')
    print('Example data point at data[5,4,3]:')
    print('{}\t{}\t{}\t{}'.format('r[3]', 'theta[4]', 'phi[5]', 'value'))
    print('{}\t{}\t{}\t{}'.format(rvec[3], tvec[4], pvec[5], data[5,4,3]))


if __name__ == '__main__':
    main()


