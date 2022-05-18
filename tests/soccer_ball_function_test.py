#!/usr/bin/env python
import numpy as np
from scipy.special import sph_harm
#
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm, colors
#from matplotlib.colors import LightSource

import psihdf as ps
import psipals

psipals.load()
cmap2use = cm.get_cmap('psi_blue_red', 1024)

# Create 1D grid arrays for theta and phi:
n_theta_vec = np.array([16,32,64,128,256,512,1024,2048,4096,8192])
n_phi_vec = 2.0*n_theta_vec

for i in range(len(n_theta_vec)):

    n_theta = int(n_theta_vec[i])
    n_phi = int(n_phi_vec[i])

    theta_1d = np.linspace(0, np.pi, n_theta, dtype=np.float64)
    phi_1d = np.linspace(0, 2.0*np.pi, n_phi, dtype=np.float64)

    # Create 3D grid arrays from the 1D grid arrays:
    theta_3d, phi_3d = np.meshgrid(theta_1d, phi_1d)

    # Create the 3D Cartesian grid (unit sphere, so r=1)
    x_3d = np.sin(theta_3d) * np.cos(phi_3d)
    y_3d = np.sin(theta_3d) * np.sin(phi_3d)
    z_3d = np.cos(theta_3d)

    # Calculate the spherical harmonic soccer ball solution:

    f1 = sph_harm(0, 6, phi_3d, theta_3d).real
    f2 = np.sqrt(2)*sph_harm(5, 6, phi_3d, theta_3d).real
    f = f1 + np.sqrt(14.0/11.0)*f2

    #ps.wrhdf_2d('y_0_6_scipy_pt.h5',phi_1d,theta_1d,f1.transpose())
    #ps.wrhdf_2d('y_5_6_scipy_pt.h5',phi_1d,theta_1d,f2.transpose())


    #fmax, fmin = f.max(), f.min()
    #f = (f - fmin)/(fmax - fmin)
    #f = -1000.0 + 2000.0*f

    f = f*1000.0

    fmax, fmin = f.max(), f.min()

    #absfmnmx=np.max([np.abs(fmax),np.abs(fmin)])

    ps.wrhdf_2d('soccer_ball_scipy_pt_t'+str(n_theta)+'_p'+str(n_phi)+'.h5',phi_1d,theta_1d,f.transpose())

    ### PLOTTING ###
    continue

    fig = plt.figure(figsize=(6, 6), dpi=120, facecolor='k',frameon=False)
    fig.add_subplot(projection='3d')
    ax = plt.gca()
    ax.set_facecolor('k')

    # Set the lightsource location:
    #lm = LightSource(azdeg=0., altdeg=0., hsv_min_val=0, hsv_max_val=1, hsv_min_sat=1, hsv_max_sat=0)

    # Plot the surface data:
    norm = matplotlib.colors.Normalize(vmin=f.min(), vmax=f.max())
    plot_handle = ax.plot_surface(x_3d, y_3d, z_3d, rstride=1, cstride=1, \
                              facecolors = cmap2use(norm(f)), \
                              shade=False, linewidth=0, vmin=0, vmax=1, antialiased=False)
    sz=0.63
    ax.axes.set_xlim3d(left=-sz, right=sz)
    ax.axes.set_ylim3d(bottom=-sz, top=sz) 
    ax.axes.set_zlim3d(bottom=-sz, top=sz) 
                              
    # Set the view:
    ax.view_init(azim=0.,elev=0.)

    # Turn off the axis planes:
    ax.set_axis_off()

    # Set limits of plot and aspect ratio:
    limits = np.array([getattr(ax, f'get_{axis}lim')() for axis in 'xyz'])
    ax.set_box_aspect(np.ptp(limits, axis = 1))

    # Output plot:
    fig.savefig('soccer_ball_function_t'+str(n_theta)+'_p'+str(n_phi)+'_ball.png', bbox_inches="tight", pad_inches=0, dpi=120, facecolor=fig.get_facecolor(), edgecolor=None)


