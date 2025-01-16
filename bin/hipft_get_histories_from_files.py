#!/usr/bin/env python3
import numpy as np
import pandas as pd
import argparse
import h5py
import glob
from datetime import datetime, timezone

# Version 1.1.0

def argParsing():

    parser = argparse.ArgumentParser(description='HIPFT: Generate history files from sequence of hipft-formatted output maps.')

    parser.add_argument('-folder',
        help='Name of folder.',
        dest='folder',
        required=True)

    parser.add_argument('-bfile',
        help='Base file name.',
        dest='bfile',
        default='hipft_brmap_idx',
        required=False)

    parser.add_argument('-t0',
        help='Sequence start index.',
        dest='t0',
        default=1,
        required=False)

    parser.add_argument('-tf',
        help='Sequence stop index.',
        dest='tf',
        required=True)

    parser.add_argument('-cadence',
        help='Cadence for time in hours (default: 1).',
        dest='cadence',
        default=1.0,
        required=False)

    parser.add_argument('-o',
        help='Output file name (default: hipft_history_sol.out).',
        dest='oname',
        type=str,
        default='hipft_history_sol.out',
        required=False)

    
    return parser.parse_args()
    

# hipft_analysis_step
#
# Compute analysis of HIPFT run.
# Takes in a folder name, base file name, and sequence start and stop indexes.
# Produces the HIPFT history (analysis) file from the maps.


def run(args):

    time=0.0
    cadence = float(args.cadence)

    #Calculation constants
    io_hist_sol_filename = args.oname
    pole_flux_lat_limit = 30.  # Make this an input parameter!
    d2r = 0.017453292519943295
    pi = 3.1415926535897932
    pi_i = 0.3183098861837907
    rsun_cm2 = 4.84416e21

    #Header
    header_txt = '      STEP                   TIME                 BR_MIN                 BR_MAX         \
    BR_ABS_MIN          FLUX_POSITIVE          FLUX_NEGATIVE    NPOLE_FLUX_POSITIVE    NPOLE_FLUX_NEGATIVE\
    SPOLE_FLUX_POSITIVE    SPOLE_FLUX_NEGATIVE             NPOLE_AREA             SPOLE_AREA          \
    EQ_DIPOLE              AX_DIPOLE  VALIDATION_ERR_CVRMSD \n'
    f_out = open(io_hist_sol_filename, "w")
    f_out.write(header_txt)

    #Count for writing to history file
    first_file = True

    #Loop through files in folder from start to stop time.
    for idx in range(int(args.t0),int(args.tf)+1):
        filename=args.folder+"/"+args.bfile+"{:06d}".format(idx)+".h5"
    
        #Initialize variables
        h_minbr = 3.40282347e38
        h_maxbr = 1.17549435e-38
        h_minabsbr = 3.40282347e38
        h_valerr = 0.
        sumfval2 = 0.

        file = h5py.File(filename, "r")
        f = file['Data'][()]
        if first_file:
            ptemp = file['dim1'][()]
            npm = np.size(ptemp)+1
            t = file['dim2'][()]
            ntm = np.size(t)
            p = np.zeros(npm)
            p[0:-1]=ptemp
            p[npm-1]=p[-2]+(p[-2]-p[-3])

            dth = np.zeros(ntm+1)
            th = np.zeros(ntm+1)
            for i in range(1,ntm):
                th[i] = 0.5*(t[i] + t[i-1])
                dth[i] = t[i]-t[i-1]
            th[0] = th[1] - dth[1]
            th[ntm] = th[ntm-1] + dth[ntm-1]
            dth[0] = dth[1]
            dth[ntm] = dth[ntm-1]

            dt = np.zeros(ntm)
            for i in range(0,ntm):
                dt[i] = th[i+1]-th[i]

            ph = np.zeros(npm)
            dph = np.zeros(npm)
            for i in range(1,npm):
                ph[i] = 0.5*(p[i] + p[i-1])
                dph[i] = p[i]-p[i-1]
            ph[0] = ph[npm-2] - 2*pi
            dph[0] = dph[npm-2]

            dp = np.zeros(npm)
            for i in range(0,npm-1):
                dp[i] = ph[i+1] - ph[i]
            dp[npm-1] = dp[1]
            
            first_file = False
            print('First step done - grid calculated.')
        else:
            print('Time '+str(time)+' IDX: '+str(idx)+' completed.')
            time = time + cadence

        ntime = int(filename[-9:-3])-1

        #Get max and min metrics
        for j in range(0,npm-2):
            for i in range(0,ntm):
                h_minbr = min(f[i,j],h_minbr)
                h_maxbr = max(f[i,j],h_maxbr)
                h_minabsbr = min(abs(f[i,j]),h_minabsbr)


        #Get integrated metrics
        h_fluxp = 0.
        h_fluxm = 0.
        h_fluxp_pn = 0.
        h_fluxm_pn = 0.
        h_fluxp_ps = 0.
        h_fluxm_ps = 0.
        h_area_pn = 0.
        h_area_ps = 0.
        eqd1 = 0.
        eqd2 = 0.
        h_eq_dipole = 0.
        h_ax_dipole = 0.

        for i in range(0,npm-1):
            for j in range(0,ntm):
                if j==0:
                    tav=0.5*(t[0]+t[1])
                    sn_t = np.sin(tav)
                    cs_t = np.cos(tav)
                    d_t = dth[1]
                    da_t=0.25*sn_t*d_t
                elif j==ntm-1:
                    tav=0.5*(t[ntm-1]+t[ntm-2])
                    sn_t = np.sin(tav)
                    cs_t = np.cos(tav)
                    d_t = dth[ntm-1]
                    da_t=0.25*sn_t*d_t    
                else:
                    sn_t = np.sin(t[j])
                    cs_t = np.cos(t[j])
                    d_t = dt[j]
                    da_t = sn_t*d_t    

                if i==0:
                    da_p=0.5*dph[0]
                elif i==npm-2:
                    da_p=0.5*dph[npm-2]
                else:
                    da_p=dp[i]  

                cs_p = np.cos(p[i])
                sn_p = np.sin(p[i])

                fv = f[j,i]*da_t*da_p

                #Fluxes
                if fv>0:
                    h_fluxp = h_fluxp + fv
                else:
                    h_fluxm = h_fluxm + fv

                #Polar fluxes and areas
                if t[j]<pole_flux_lat_limit*d2r:
                    if fv>0:
                        h_fluxp_pn = h_fluxp_pn + fv
                    else:
                        h_fluxm_pn = h_fluxm_pn + fv
                    h_area_pn = h_area_pn + da_t*da_p

                if t[j]>pi-pole_flux_lat_limit*d2r:
                    if fv>0:
                        h_fluxp_ps = h_fluxp_ps + fv
                    else:
                        h_fluxm_ps = h_fluxm_ps + fv
                    h_area_ps = h_area_ps + da_t*da_p

                #Dipoles
                eqd1 = eqd1 + f[j,i]*sn_t*cs_p*da_t*da_p
                eqd2 = eqd2 + f[j,i]*sn_t*sn_p*da_t*da_p

                h_ax_dipole = h_ax_dipole + f[j,i]*cs_t*da_t*da_p

        #Set equatorial dipole strength
        eqd1 = 0.75*pi_i*eqd1
        eqd2 = 0.75*pi_i*eqd2
        h_eq_dipole = np.sqrt(eqd1**2 + eqd2**2)

        #Set axial dipole strength
        h_ax_dipole = 0.75*pi_i*h_ax_dipole
        
        #Set fluxes to be in units of Mx
        h_fluxp    = rsun_cm2*h_fluxp
        h_fluxm    = rsun_cm2*h_fluxm
        h_fluxp_pn = rsun_cm2*h_fluxp_pn
        h_fluxm_pn = rsun_cm2*h_fluxm_pn
        h_fluxp_ps = rsun_cm2*h_fluxp_ps
        h_fluxm_ps = rsun_cm2*h_fluxm_ps

        #Set polar areas to be in units of cm
        h_area_pn = rsun_cm2*h_area_pn
        h_area_ps = rsun_cm2*h_area_ps

        #Write out histories
        hist_sol='         %.0f %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e \n' %(\
                ntime, time, h_minbr, h_maxbr, h_minabsbr, h_fluxp, h_fluxm, h_fluxp_pn, \
                h_fluxm_pn, h_fluxp_ps, h_fluxm_ps, h_area_pn, h_area_ps, h_eq_dipole,h_ax_dipole, h_valerr)
                
        f_out.write(hist_sol)
        
    f_out.close()
    


def main():
    ## Get input agruments:
    args = argParsing()
    run(args)

if __name__ == '__main__':
    main()
            


