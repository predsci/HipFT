#!/usr/bin/env python3
import psihdf as ps
import numpy as np
import argparse

def argParsing():

    parser = argparse.ArgumentParser(description='Read PSI 2D hdf data.')

    parser.add_argument("psi_2d_hdf_file_name",
        help='Name of 2D PSI HDF file (e.g. br_photo.h5).')

    args = parser.parse_args()

    return parser.parse_args()

def main():

    ## Get input agruments:
    args = argParsing()

    xvec, yvec, data = ps.rdhdf_2d(args.psi_2d_hdf_file_name)

    # If the data is in tp format, transpose to pt:
    # This makes plotting better/easier.
    if (np.max(yvec) > 3.5):
      tmpt = xvec
      tmpp = yvec
      yvec = tmpt
      xvec = tmpp
      data = np.transpose(data)

    # If the data is in sinlat, convert scale to theta:
    # No interpolation!  Rather, data is associated with nonuniform theta coordinates)
    if (np.min(yvec) < -0.5):
      yvec = np.arccos(yvec)

    tvec = yvec
    pvec = xvec

    tmin = np.min(tvec)
    tmax = np.max(tvec)
    pmin = np.min(pvec)
    pmax = np.max(pvec)

    NT = tvec.size
    NP = pvec.size

    print('Opened file:   '+args.psi_2d_hdf_file_name)
    print(' ')
    print('Shape of data:'+str(data.shape))
    print(' ')
    print('Ntheta:     '+str(NT))
    print('Nphi:       '+str(NP))
    print(' ')
    print('min(theta): '+str(np.min(tvec)))
    print('max(theta): '+str(np.max(tvec)))
    print('min(phi):   '+str(np.min(pvec)))
    print('max(phi):   '+str(np.max(pvec)))
    print(' ')
    print('min(data):  '+str(np.min(data)))
    print('max(data):  '+str(np.max(data)))
    print('mean(data): '+str(np.mean(data)))
    print('Example data point at data[5,3]:')
    print('{}\t{}\t{}'.format('theta[5]', 'phi[3]', 'value'))
    print('{}\t{}\t{}'.format(tvec[5], pvec[3], data[5,3]))


if __name__ == '__main__':
    main()


