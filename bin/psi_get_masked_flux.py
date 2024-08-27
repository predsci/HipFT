#!/usr/bin/env python3
import numpy as np
import sys
import signal
import ntpath
import argparse
from scipy import ndimage
from scipy.integrate import simps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
#
import psipals
import psihdf as ps

# GET_MASKED_FLUX: Version 1.1.0
#
# Get flux results of mask overlaid on BR map.
# Mask can be floating point or binary.
# If using a mask file, both the mask map and BR map must be same resolution with periodic points included.
# Signed Flux in each detected mask region is computed, and totals.
#
# 1.0.3: Changed per-CH area calculation to include floating point mask.
# 1.0.4: Fixed issue in detection of staggered phi.  It was detecting non-staggered phi as staggered!
# 1.1.0: Added mask box mode, quiet mode now prints all quantites.

def signal_handler(signal, frame):
        print('You pressed Ctrl+C! Stopping!')
        sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)

def argParsing():

    parser = argparse.ArgumentParser(description='get_masked_flux')

    parser.add_argument('br_file',help='Name of 2D Br hdf file in phi-theta coordinates.')

    parser.add_argument('-br_unit_fac',
                        help='Unit factor for Br such that br_unit_fac*Br is in Gauss. For MAS code data, you can use the special keyword "mas".  If Br is already in Gauss, you can use the keyword "gauss" (the default).',
                        dest='br_unit_fac',
                        default='gauss',
                        required=False)

    parser.add_argument('-mask_file',
                        help='Name of 2D mask hdf file in phi-theta coordinates.  Mask can be 0,1 or -1,0,1 and is modded. Default is full map.',
                        dest='mask_file',
                        required=False)

    parser.add_argument('-mask_box',
                        help='Use a box for the mask.  Set this to p1,p2,t1,t2 representing the corners of the phi-theta box in DEGREES.',
                        dest='mask_box',
                        type=str,
                        required=False)

    parser.add_argument('-r',
                        help='Assumed radii of maps in Rs units (default 1).',
                        dest='r0fac',
                        default=1.0,
                        required=False)

    parser.add_argument('-q',
                        help='Output values in one line.',
                        dest='quiet',
                        action='store_true',
                        default=False,
                        required=False)

    parser.add_argument('-plot',
                        help='Output plots. Names based on br and mask file names.',
                        dest='plot',
                        action='store_true',
                        default=False,
                        required=False)

    parser.add_argument('-k',
                        help='Plot with black background',
                        dest='k',
                        action='store_true',
                        default=False,
                        required=False)

    parser.add_argument('-sinlat',
                        help='Plot in Sin-latitute',
                        dest='sinlat',
                        action='store_true',
                        default=False,
                        required=False)

    parser.add_argument('-crange',
                        help='Colormap limits of magnetic field plots (-crange,crange)',
                        dest='crange',
                        default=15,
                        required=False)

    parser.add_argument('-max_num_ch_plot',
                        help='Maximum number of mask regions to plot flux percentage',
                        dest='max_num_ch_plot',
                        default=5,
                        required=False)

    parser.add_argument('-v',
                        help='Verbose output listing flux of each CH',
                        dest='verb',
                        action='store_true',
                        default=False,
                        required=False)


    return parser.parse_args()

r_s=np.double(6.957e10)
r1au=np.double(1.496e13)
gauss2nt=np.double(1e5)
flux_fac=np.double(1e-21)

## Get input agruments:
args = argParsing()

#Check unit factor:
if (args.br_unit_fac=='mas'):
    br_unit_fac=np.double(2.20689140)
elif (args.br_unit_fac=='Gauss'):
    br_unit_fac=np.double(1.0)
elif (args.br_unit_fac=='gauss'):
    br_unit_fac=np.double(1.0)
else:
    br_unit_fac=np.double(args.br_unit_fac)

max_num_holes_to_plot = int(args.max_num_ch_plot)

#Get radii:
r0 = r_s*np.double(args.r0fac)

#Read the Br data and apply unit fac:
pvec,tvec,br_data = ps.rdhdf_2d(args.br_file)
br_data = np.multiply(br_data,br_unit_fac)

if (args.mask_file is not None):
   #Read the mask data:
   pvec2,tvec2,mask_data = ps.rdhdf_2d(args.mask_file)
elif (args.mask_box is not None):
   pvec2=pvec
   tvec2=tvec
   T,P = np.meshgrid(tvec2,pvec2,indexing='ij')
   mask_data=np.zeros(br_data.shape)
   arg_dict = vars(args)
   box_coords = np.array(arg_dict['mask_box'].split(','))
   p1 = (np.pi/180.0)*np.float64(box_coords[0])
   p2 = (np.pi/180.0)*np.float64(box_coords[1])
   t1 = (np.pi/180.0)*np.float64(box_coords[2])
   t2 = (np.pi/180.0)*np.float64(box_coords[3])
   boxarea = ((T>t1) & (T<t2) & (P>p1) & (P<p2))
   mask_data[boxarea] = 1.0
   args.mask_file = 'box.h5'
else:
   pvec2=pvec
   tvec2=tvec 
   mask_data=np.ones(br_data.shape)

#If the data is in pt format, transpose to tp:
if (np.max(tvec)>3.5):
    tmpt=tvec
    tmpp=pvec
    tvec=tmpp
    pvec=tmpt
    br_data=np.transpose(br_data)

#If the data is in pt format, transpose to tp:
if (np.max(tvec2)>3.5):
    tmpt=tvec2
    tmpp=pvec2
    tvec2=tmpp
    pvec2=tmpt
    mask_data=np.transpose(mask_data)

eps=1e-6

if (not np.allclose(tvec,tvec2,rtol=eps) or not np.allclose(pvec,pvec2,rtol=eps)):
   print ("Error! Br map and Mask do not have the same scales!")
   quit()

#Check for sinlat and convert to theta:
if (np.min(tvec)<-0.5):
    tvec=np.acos(tvec)

#Make sure |mask| <=1
#and set boundaries for color plots:
if (np.max(np.abs(mask_data))>1+eps):
   print ("Error!  Mask has values |x|>1!")
   quit()

cmax_mask=1
if (np.min(mask_data)<-eps):
     cmin_mask = -1
     mask_data_pm = mask_data
else:
     cmin_mask = 0

mask_data = np.abs(mask_data)

#Check for staggering past domain bounderies.
#If staggered in theta, replace end points with interpolation and reset theta points.
#If staggered in phi, remove one duplicate point to get proper integrals.

if (tvec[0]<0 and tvec[1]>0 and tvec[-2]<np.pi and tvec[-1]>np.pi):
    print("Theta is staggered!")
    # Find interpolated values at 0 and pi using: y = y0 + (x-x0)/(x1-x0) * (y1-y0)
    br_data[ 0,:]  = br_data[ 0,:] + ((  0.0-tvec[ 0])/(tvec[ 1]-tvec[ 0]))*(br_data[ 1,:]-br_data[ 0,:])
    br_data[-1,:]  = br_data[-1,:] + ((np.pi-tvec[-1])/(tvec[-2]-tvec[-1]))*(br_data[-2,:]-br_data[-1,:])
    mask_data[ 0,:]  = mask_data[ 0,:] + ((  0.0-tvec[ 0])/(tvec[ 1]-tvec[ 0]))*(mask_data[ 1,:]-mask_data[ 0,:])
    mask_data[-1,:]  = mask_data[-1,:] + ((np.pi-tvec[-1])/(tvec[-2]-tvec[-1]))*(mask_data[-2,:]-mask_data[-1,:])
    #Reset scale values to 0 and pi:
    tvec[0] = 0.0
    tvec[-1] = np.pi
elif (tvec[0]>eps and tvec[-1]<np.pi-eps):
    # Inner half-mesh, need to create polar condition.
    # Need to ADD theta boundary points.
    # Weighted average Sum of br. set value to boundary.
    print ("Error!  Inner theta half-mesh not implemented yet!")
    quit()    
elif (not (abs(tvec[0])<=eps and abs(tvec[-1]-np.pi)<=eps)):
    print ("Error!  Full theta domain not covered!")
    quit()

if ((pvec[-1]-pvec[0])>2.0*np.pi+eps):
    print("Phi is staggered!")
    br_data = np.delete(br_data,[0],axis=1)
    mask_data = np.delete(mask_data,[0],axis=1)
    pvec = np.delete(pvec,[0])
elif (not (abs(pvec[-1]-pvec[0]-2.0*np.pi)<=eps)):
    print ("Error!  Full phi domain not covered!")
    quit()

#Create coordinates:
T,P = np.meshgrid(tvec,pvec,indexing='ij')
r2SinT = (r0**2)*np.sin(T)

#Get total area of Br to test integration:
int_tmp = r2SinT
int_tmp = simps(int_tmp,x=tvec,axis=0)
total_area = simps(int_tmp,x=pvec)

#Get total positive flux of Br:
int_tmp = (br_data>0)*br_data*r2SinT
int_tmp = simps(int_tmp,x=tvec,axis=0)
total_positive_flux_br = simps(int_tmp,x=pvec)

#Get total negative flux of Br:
int_tmp = (br_data<0)*br_data*r2SinT
int_tmp = simps(int_tmp,x=tvec,axis=0)
total_negative_flux_br = simps(int_tmp,x=pvec)

#Get total unsigned flux of Br:
total_unsigned_flux_br = np.abs(total_positive_flux_br)+np.abs(total_negative_flux_br)

#Get total signed flux of Br:
total_signed_flux_br = total_positive_flux_br+total_negative_flux_br

#Get total area of masked region:
int_tmp = mask_data*r2SinT
int_tmp = simps(int_tmp,x=tvec,axis=0)
total_masked_area = simps(int_tmp,x=pvec)

#Get total positive flux of masked Br:
int_tmp = (mask_data*br_data>0)*mask_data*br_data*r2SinT
int_tmp = simps(int_tmp,x=tvec,axis=0)
total_positive_flux = simps(int_tmp,x=pvec)

#Get total negative flux of masked Br:
int_tmp = (mask_data*br_data<0)*mask_data*br_data*r2SinT
int_tmp = simps(int_tmp,x=tvec,axis=0)
total_negative_flux = simps(int_tmp,x=pvec)

#Get total unsigned flux of Br:
total_unsigned_flux = np.abs(total_positive_flux)+np.abs(total_negative_flux)

#Get total signed flux of masked Br:
total_signed_flux = total_positive_flux+total_negative_flux

#Get individual mask areas:
CHS, num_CHS  = ndimage.label(mask_data)

if (num_CHS>0):

    signed_flux_per_hole   = np.zeros(num_CHS)
    area_per_hole          = np.zeros(num_CHS)
    unsigned_flux_per_hole = np.zeros(num_CHS)

#i=0 is the background.  Multiply by mask_data in case one is using a floating-point mask.
    for i in range(num_CHS):
        CH_IDX=(CHS==i+1)

        int_tmp = CH_IDX*mask_data*br_data*r2SinT
        int_tmp = simps(int_tmp,x=tvec,axis=0)
        signed_flux_per_hole[i] = simps(int_tmp,x=pvec)

        int_tmp = CH_IDX*mask_data*np.abs(br_data)*r2SinT
        int_tmp = simps(int_tmp,x=tvec,axis=0)
        unsigned_flux_per_hole[i] = simps(int_tmp,x=pvec)

        int_tmp = CH_IDX*mask_data*r2SinT
        int_tmp = simps(int_tmp,x=tvec,axis=0)
        area_per_hole[i] = simps(int_tmp,x=pvec)

    total_open_flux_allholes     = np.sum(np.abs(signed_flux_per_hole))
    total_signed_flux_allholes   = np.sum(signed_flux_per_hole)
    total_unsigned_flux_allholes = np.sum(unsigned_flux_per_hole)
    total_area_allholes          = np.sum(area_per_hole)
    Br1au = gauss2nt*total_open_flux_allholes/(4.0*np.pi*r1au**2)

    if (args.verb):
        print(" ")
        print("#### FLUX PER HOLE DATA SORTED BY ABS(SIGNED_FLUX) ####")
        print("CH#  AREA(1e21 cm^2)  SIGNED_FLUX (1e21 Mx)  UNSIGNED_FLUX (1e21 MX)")

        sf_idx_sort=np.argsort(np.abs(signed_flux_per_hole))
        sf_idx_sort=np.flip(sf_idx_sort,axis=0)

        
        for i in range(num_CHS):
            j=sf_idx_sort[i]
            print("%3d  %8.4f  %8.4f  %8.4f" % (i+1,(flux_fac*area_per_hole[j]),(flux_fac*signed_flux_per_hole[j]),(flux_fac*unsigned_flux_per_hole[j])))
        print(" ")
        print("#################### SUMMARY DATA ##########################")


    if (args.plot):
        relative_flux = np.zeros(num_CHS)
        relative_area = np.zeros(num_CHS)
        results_vec   = np.zeros([num_CHS,6])

        num_holes_to_plot=np.min([max_num_holes_to_plot,num_CHS])

        for chm_idx in range(num_CHS):
            relative_flux[chm_idx] = np.abs(signed_flux_per_hole[chm_idx])/total_open_flux_allholes
            relative_area[chm_idx] = np.abs(area_per_hole[chm_idx])/total_area_allholes
            results_vec[chm_idx,:] = np.array([chm_idx+1, area_per_hole[chm_idx], signed_flux_per_hole[chm_idx], unsigned_flux_per_hole[chm_idx], relative_area[chm_idx], relative_flux[chm_idx]])

        results_tmp = results_vec
        results_tmp[:,2] = np.abs(results_tmp[:,2])
        results_sorted_vec = results_tmp[results_tmp[:,2].argsort(),]
        results_sorted_vec = np.flipud(results_sorted_vec)
        results_sorted_vec[:,2] = results_vec[np.asarray(results_sorted_vec[:,0]-1,dtype=int),2]

        CH_sorted = 0*CHS
        for ii in range(num_CHS):
            CH_sorted[CHS==results_sorted_vec[ii,0]]=ii+1
        CH_sorted[CH_sorted>num_holes_to_plot]=num_holes_to_plot

        fluxtemp=0
        #Get percentage of flux from num_holes_to_plot holes:
        for ii in range(num_holes_to_plot-1):
           fluxtemp = fluxtemp + np.abs(results_sorted_vec[ii,2])
        open_flux_for_colored_holes = fluxtemp

        fluxtemp=0
        #Get percentage of flux from num_holes_to_plot holes:
        for ii in range(num_holes_to_plot-1,num_CHS):
            fluxtemp = fluxtemp + np.abs(results_sorted_vec[ii,2])
        open_flux_for_rest_of_holes = fluxtemp
        remaining_per = open_flux_for_rest_of_holes/total_open_flux_allholes
else:
    total_open_flux_allholes     = -9999
    total_signed_flux_allholes   = -9999
    total_unsigned_flux_allholes = -9999
    total_area_allholes          = -9999
    Br1au = gauss2nt*total_unsigned_flux_br/(4.0*np.pi*r1au**2)

if (not args.quiet):
    print("")
    print("Integration test:")
    print("  Total surface area (true):       %12.8f (1e21) cm^2" % (flux_fac*4.0*np.pi*(r0**2)))
    print("  Total surface area (int):        %12.8f (1e21) cm^2" % (flux_fac*total_area))
    print("")
    print("Full map results:")
    print("  Total unsigned flux of BR:       %12.8f (1e21) Mx" % (flux_fac*total_unsigned_flux_br))
    print("  Total positive flux of BR:       %12.8f (1e21) Mx" % (flux_fac*total_positive_flux_br))
    print("  Total negative flux of BR:       %12.8f (1e21) Mx" % (flux_fac*total_negative_flux_br))
    print("  Total   signed flux of BR:       %12.8f (1e21) Mx" % (flux_fac*total_signed_flux_br))
    favg=0.5*total_unsigned_flux_br
    print("  Fractional flux imbalance of BR:           %12.8e" % (total_signed_flux_br/favg))
    print(" ")
    if (args.mask_file is not None or args.mask_box is not None):
        print("Mask results:")
        print("  Masked surface area:             %8.4f (1e21) cm^2" % (flux_fac*total_masked_area))
        print("  Total unsigned flux of MASK*BR:  %8.4f (1e21) Mx" % (flux_fac*total_unsigned_flux))
        print("  Total positive flux of MASK*BR:  %8.4f (1e21) Mx" % (flux_fac*total_positive_flux))
        print("  Total negative flux of MASK*BR:  %8.4f (1e21) Mx" % (flux_fac*total_negative_flux))
        print("  Total   signed flux of MASK*BR:  %8.4f (1e21) Mx" % (flux_fac*total_signed_flux))
        favg=0.5*total_unsigned_flux
        print("  Fractional flux imbalance of MASK*BR:     %12.8e" % (total_signed_flux/favg))
        print(" ")
        if (num_CHS>1):
            print("Number of individual masked regions (MR): %8i" % (num_CHS))
            print(" ")
            print("Total  surface area of all MRs:  %8.4f (1e21) cm^2" % (flux_fac*total_area_allholes))
            print("Total   signed flux of all MRs:  %8.4f (1e21) Mx" % (flux_fac*total_signed_flux_allholes))
            print("Total unsigned flux of all MRs:  %8.4f (1e21) Mx" % (flux_fac*total_unsigned_flux_allholes))
            print("Total     open flux of all MRs:  %8.4f (1e21) Mx" % (flux_fac*total_open_flux_allholes))
            print(" ")
    print("Extrapolated IMF at 1AU (assuming mask is open field map):  %8.4f nT" % (Br1au))
else:
    print("%16.8f,%16.8f,%16.8f,%16.8f,%16.8f,%16.8f,%16.8f,%16.8f,%16.8f,%16.8f,%16.8f,%16.8f,%16.8f,%16.8f" % \
    (np.double(args.r0fac),\
    (flux_fac*total_area),(flux_fac*total_unsigned_flux_br),(flux_fac*total_positive_flux_br),(flux_fac*total_negative_flux_br),\
    (flux_fac*total_signed_flux_br),(total_signed_flux_br/(0.5*total_unsigned_flux_br)),\
    (flux_fac*total_masked_area),(flux_fac*total_unsigned_flux),(flux_fac*total_positive_flux),(flux_fac*total_negative_flux),\
    (flux_fac*total_signed_flux),(total_signed_flux/(0.5*total_unsigned_flux)),Br1au))

#print(" ")
#code_fac=flux_fac*br_unit_fac*r0**2
#print("Flux unit conversion factor:           "+str(code_fac))
#print("Total pre-unit positive flux of BR:   %20.16f" % (flux_fac*total_positive_flux_br/code_fac))
#print("Total pre-unit negative flux of BR:   %20.16f" % (flux_fac*total_negative_flux_br/code_fac))

if (args.plot):
    if (not args.quiet):
        print(" ")
        print("Plotting...")

# Make base filenames
    if (str(args.br_file).endswith('h5')):
       br_ofile = args.br_file[0:len(args.br_file)-3]
    else:
       br_ofile = args.br_file[0:len(args.br_file)-4]

    if (args.mask_file is not None):
        if (str(args.mask_file).endswith('h5')):
           mask_ofile = args.mask_file[0:len(args.mask_file)-3]
        else:
           mask_ofile = args.mask_file[0:len(args.mask_file)-4]

    br_ofile = ntpath.basename(br_ofile)
    if (args.mask_file is not None): mask_ofile = ntpath.basename(mask_ofile)
    if (args.mask_file is not None): ofile_base = mask_ofile+'_OVER_'+br_ofile

    cmin=-float(args.crange)
    cmax=float(args.crange)

    fsize=18

    if (args.k):
        fc='k'
        tc='w'
    else:
        fc='w'
        tc='k'

# Load color maps:
    psipals.load()

    if (args.sinlat):
        tvec=np.cos(tvec)

# Set up scales for plots:
    xvec_plot2=np.zeros(len(pvec)+1)
    xvec_plot2[1:-1] = (pvec[0:-1] + pvec[1:])/2.0
    xvec_plot2[0]    = pvec[1]  - (pvec[1]  - pvec[0])
    xvec_plot2[-1]   = pvec[-2] + (pvec[-1] - pvec[-2])

    yvec_plot2=np.zeros(len(tvec)+1)
    yvec_plot2[1:-1] = (tvec[0:-1] + tvec[1:])/2.0
    yvec_plot2[0]    = yvec_plot2[1]  - (tvec[1]  - tvec[0])
    yvec_plot2[-1]   = yvec_plot2[-2] + (tvec[-1] - tvec[-2])

    pvec_plot = xvec_plot2
    tvec_plot = yvec_plot2

#  Plot Br
    fig,ax = plt.subplots(num=None, figsize=[14,7], dpi=200, facecolor=fc, frameon=False)

    plot_mag = plt.pcolormesh(pvec_plot, tvec_plot, br_data, zorder=1)
    plt.set_cmap('psi_blue_red')
    plt.clim([cmin,cmax])
    plt.xlim(xmin=0,xmax=2*np.pi)
    plt.xlabel('Longitude',{'fontsize': fsize, 'color': tc})
    plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi),(0,90,180,270,360), color=tc, size=fsize)
    if(args.sinlat):
        plt.ylim(ymin=-1,ymax=1)
        plt.ylabel('Sin-Latitude',{'fontsize': fsize, 'color': tc})
        plt.yticks((-1,-1.0/2,0,1.0/2,1),('$-1$','$-1/2$','$0$','$1/2$','$1$'))
        cb=plt.colorbar(plot_mag,fraction=0.015, pad=0.02, aspect=20)
    else:
        plt.ylim(ymin=0,ymax=np.pi)
        plt.ylabel('Colatitude',{'fontsize': fsize, 'color': tc})
        plt.yticks((0,np.pi/4,np.pi/2,3*np.pi/4,np.pi),('0','45','90','135','180'))
        plt.gca().invert_yaxis()
        cb=plt.colorbar(plot_mag,fraction=0.024, pad=0.02, aspect=20)

    plt.tick_params(axis='both', color=tc, labelsize=fsize, size=fsize)
    ax.tick_params(axis='x', colors=tc, which='both')
    ax.tick_params(axis='y', colors=tc, which='both')

    ax.set_aspect('equal')
    ax.set_facecolor(fc)
    ax.spines['bottom'].set_color(tc)
    ax.spines['top'].set_color(tc)
    ax.spines['right'].set_color(tc)
    ax.spines['left'].set_color(tc)
    ax.grid(zorder=2)

    cbstr='Gauss'

    cb.outline.set_edgecolor(tc)
    cb.set_label(cbstr,fontsize=fsize,color=tc)
    cb.ax.yaxis.set_tick_params(labelsize=fsize,color=tc)
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=tc)
 
    if args.mask_box is not None:
        points = [(p1,t1),(p1,t2),(p2,t2),(p2,t1)]
        line = plt.Polygon(points, closed=True, fill=None, edgecolor='g', linestyle='--',linewidth=2)
        ax.add_patch(line)

    fig.savefig("BR_"+br_ofile+".png",bbox_inches='tight',dpi=200,facecolor=fig.get_facecolor(), edgecolor='none')

    plt.close('all')

    if (args.mask_file is None):
       if (not args.quiet):
           print("Done!")
       quit()

#  Plot Mask
    fig,ax = plt.subplots(num=None, figsize=[14,7], dpi=200, facecolor=fc)
    if (cmin_mask==-1):
        plot_mask = plt.pcolormesh(pvec_plot, tvec_plot, mask_data_pm, zorder=1)
        plt.set_cmap('psi_blue_red')
        cbticks=[-1, 0, 1]
        cbtick_labels=['CH (-)', 'No CH', 'CH (+)']
    else:
        plot_mask = plt.pcolormesh(pvec_plot, tvec_plot, mask_data, zorder=1)
        plt.set_cmap('psi_inverted_grayscale')
        cbticks=[0, 1]
        cbtick_labels=['No Mask','Mask']
    plt.clim([cmin_mask,cmax_mask])
    ax.spines['bottom'].set_color(tc)
    ax.spines['top'].set_color(tc)
    ax.spines['right'].set_color(tc)
    ax.spines['left'].set_color(tc)
    plt.tick_params(axis='both', color=tc, labelsize=fsize, size=fsize)
    ax.tick_params(axis='x', colors=tc, which='both')
    ax.tick_params(axis='y', colors=tc, which='both')
    ax.set_aspect('equal')
    ax.set_facecolor(fc)
    plt.xlim(xmin=0,xmax=2*np.pi)
    plt.xlabel('Longitude',{'fontsize': fsize, 'color': tc})
    plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi),(0,90,180,270,360), color=tc, size=fsize)
    if(args.sinlat):
        plt.ylim(ymin=-1,ymax=1)
        plt.ylabel('Sin-Latitude',{'fontsize': fsize, 'color': tc})
        plt.yticks((-1,-1.0/2,0,1.0/2,1),('$-1$','$-1/2$','$0$','$1/2$','$1$'))
        cb=plt.colorbar(fraction=0.015, pad=0.02, aspect=20, ticks=cbticks)
    else:
        plt.ylabel('Colatitude',{'fontsize': fsize, 'color': tc})
        plt.ylim(ymin=0,ymax=np.pi)
        plt.yticks((0,np.pi/4,np.pi/2,3*np.pi/4,np.pi),('0','45','90','135','180'))
        cb=plt.colorbar(fraction=0.024, pad=0.02, aspect=20, ticks=cbticks)
        plt.gca().invert_yaxis()
    ax.grid(zorder=2)
    cb.outline.set_edgecolor(tc)
    cb.ax.set_yticklabels(cbtick_labels)
    cb.ax.yaxis.set_tick_params(labelsize=fsize,color=tc)
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=tc)
    plt.contour(pvec, tvec, np.asarray(mask_data>0.1), levels=[0.9], colors='black',linewidths=[1],zorder=3)
    fig.savefig("MASK_"+mask_ofile+".png",bbox_inches='tight',dpi=200,facecolor=fig.get_facecolor(), edgecolor='none')

    plt.close('all')

#  Plot Br+Mask

    fig,ax = plt.subplots(num=None, figsize=[14,7], dpi=200, facecolor=fc)
    plot_maskedbr = plt.pcolormesh(pvec_plot, tvec_plot, mask_data*br_data, zorder=1)
    plt.set_cmap('psi_blue_red')
    plt.clim([cmin,cmax])
    plt.xlim(xmin=0,xmax=2*np.pi)
    plt.xlabel('Longitude',{'fontsize': fsize, 'color': tc})
    plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi),(0,90,180,270,360), color=tc, size=fsize)
    if(args.sinlat):
        plt.ylim(ymin=-1,ymax=1)
        plt.ylabel('Sin-Latitude',{'fontsize': fsize, 'color': tc})
        plt.yticks((-1,-1.0/2,0,1.0/2,1),('$-1$','$-1/2$','$0$','$1/2$','$1$'))
        cb=plt.colorbar(plot_maskedbr,fraction=0.015, pad=0.02, aspect=20)
    else:
        plt.ylim(ymin=0,ymax=np.pi)
        plt.ylabel('Colatitude',{'fontsize': fsize, 'color': tc})
        plt.yticks((0,np.pi/4,np.pi/2,3*np.pi/4,np.pi),('0','45','90','135','180'))
        plt.gca().invert_yaxis()
        cb=plt.colorbar(plot_maskedbr,fraction=0.024, pad=0.02, aspect=20)

    plt.tick_params(axis='both', color=tc, labelsize=fsize, size=fsize)
    ax.tick_params(axis='x', colors=tc, which='both')
    ax.tick_params(axis='y', colors=tc, which='both')

    ax.set_aspect('equal')
    ax.set_facecolor(fc)
    ax.spines['bottom'].set_color(tc)
    ax.spines['top'].set_color(tc)
    ax.spines['right'].set_color(tc)
    ax.spines['left'].set_color(tc)
    ax.grid(zorder=2)

    cbstr='Gauss'

    cb.outline.set_edgecolor(tc)
    cb.set_label(cbstr,fontsize=fsize,color=tc)
    cb.ax.yaxis.set_tick_params(labelsize=fsize,color=tc)
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=tc)
    plt.contour(pvec, tvec, np.asarray(mask_data>0.1), levels=[0.9], colors='black',linewidths=[1],zorder=3)
    fig.savefig("MASKED_BR_"+mask_ofile+"_"+br_ofile+".png",bbox_inches='tight',dpi=200,facecolor=fig.get_facecolor(), edgecolor='none')

    plt.close('all')

#  Plot Flux per hole

    colormapfull=[(1,1,1)]
    tmp=1
    for ii in range(50):
        colormaptmp=[(0,tmp,tmp), (tmp,0,tmp), (0,tmp,0), (tmp,tmp,0), (0,0,0), (0,0,tmp), (tmp,0,0)]
        colormapfull=colormapfull+colormaptmp
        tmp=tmp*0.75
    colormap_used=colormapfull[0:num_holes_to_plot+1]

    cmap_name = 'flux'
    cmap_flux = ListedColormap(colormap_used, name=cmap_name, N=None)
    matplotlib.colormaps.register(name=cmap_name, cmap=cmap_flux)
#    cm.register_cmap(name=cmap_name, cmap=cmap_flux)

    cbticks=np.asarray(range(num_holes_to_plot+1))+0.5
    cbtick_labels = ["%5.1f" % x for x in 100*results_sorted_vec[range(0,num_holes_to_plot-1),5]]
    cbtick_labels = [item+'%' for item in cbtick_labels]
    cbtick_labels=["     "]+cbtick_labels
    strtmp="%5.1f" % (100*remaining_per)
    cbtick_labels=cbtick_labels+[strtmp+"%"]

    fig,ax = plt.subplots(num=None, figsize=[14,7], dpi=200, facecolor=fc)
    plot_mag = plt.pcolormesh(pvec_plot, tvec_plot, CH_sorted, zorder=1, cmap=cmap_name)
    cbstr='Percentage of total open flux'
    plt.clim([0,num_holes_to_plot+1])

    plt.xlim(xmin=0,xmax=2*np.pi)
    plt.xlabel('Longitude',{'fontsize': fsize, 'color': tc})
    plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi),(0,90,180,270,360), color=tc, size=fsize)
    if(args.sinlat):
        plt.ylim(ymin=-1,ymax=1)
        plt.ylabel('Sin-Latitude',{'fontsize': fsize, 'color': tc})
        plt.yticks((-1,-1.0/2,0,1.0/2,1),('$-1$','$-1/2$','$0$','$1/2$','$1$'))
        cb=plt.colorbar(plot_mag,fraction=0.015, pad=0.02, aspect=20, ticks=cbticks)
    else:
        plt.ylim(ymin=0,ymax=np.pi)
        plt.ylabel('Colatitude',{'fontsize': fsize, 'color': tc})
        plt.yticks((0,np.pi/4,np.pi/2,3*np.pi/4,np.pi),('0','45','90','135','180'))
        plt.gca().invert_yaxis()
        cb=plt.colorbar(plot_mag,fraction=0.024, pad=0.02, aspect=20, ticks=cbticks)

    plt.tick_params(axis='both', color=tc, labelsize=fsize, size=fsize)
    ax.tick_params(axis='x', colors=tc, which='both')
    ax.tick_params(axis='y', colors=tc, which='both')

    ax.set_aspect('equal')
    ax.set_facecolor(fc)
    ax.spines['bottom'].set_color(tc)
    ax.spines['top'].set_color(tc)
    ax.spines['right'].set_color(tc)
    ax.spines['left'].set_color(tc)
    ax.grid(zorder=2)

  #  test = ma.masked_array(CH_sorted, mask=mask_data)

    cb.outline.set_edgecolor(tc)
    cb.set_label(cbstr,fontsize=fsize,color=tc)
    cb.ax.yaxis.set_tick_params(labelsize=fsize,color=tc)
    cb.ax.set_yticklabels(cbtick_labels)
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=tc)
    plt.contour(pvec, tvec, np.asarray(mask_data>0.1), levels=[0.9], colors='black',linewidths=[1],zorder=3)

    fig.savefig("FLUX_PERCENT_"+mask_ofile+"_"+br_ofile+".png", bbox_inches='tight',dpi=200,facecolor=fig.get_facecolor(), edgecolor='none')

#    for i,j in zip(plot_mag.get_facecolors(),mask_data.flatten()):
#        i[3] = j # Set the alpha value of the RGBA tuple using maskdata
#    fig.savefig("FLUX_PERCENT_"+mask_ofile+"_"+br_ofile+"_2.png", bbox_inches='tight',dpi=200,facecolor=fig.get_facecolor(), edgecolor='none')

    if (not args.quiet):
        print("Done!")

