# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile

# Import modules to interact with OS
import os
import time

# List galaxies
galaxies = ['M33', 'MW', 'M31']

# Run with all three galaxies centered on the Milky Way
for i in range(802):
    # Define the time step
    t = 14.2857*i

    # Start plot
    fig = plt.figure(figsize=(13,10))
    ax = plt.subplot(111)

    # Store the COM of the Milky Way to center on
    COMD = CenterOfMass("MW_{:0003d}.txt".format(i),2)
    COMP = COMD.COM_P(0.1, 4.0)

    # Run throught each galaxy and add to plot
    for gal in galaxies:

        # Get name of snapshot for the galaxy
        name = "{:s}_{:0003d}".format(gal,i)
        print(name)
        
        COMD = CenterOfMass("{}.txt".format(name),2)

        # Adjust to the Milky Way Center of Frame
        xD = COMD.x - float(COMP[0]/u.kpc)
        yD = COMD.y - float(COMP[1]/u.kpc)
        zD = COMD.z - float(COMP[2]/u.kpc)
    
        # Add galaxy to plot with normalized logmarithmic color scheme
        # between 1 and 1000 particles per bin
        #plt.hist2d(xD, yD, bins=250, norm=LogNorm(vmin=1e0,vmax=1e3), cmap='magma')
        # Create plot of norma
        plt.plot(xD, yD, '.', markersize=1, alpha=0.075)
        plt.title('t={:0.003f} Myr'.format(t), fontsize=22)
        
        plt.title('t={:0.003f} Myr'.format(t), fontsize=22)
        #plt.colorbar()

        # Add axis labels
        plt.xlabel('x (kpc)', fontsize=22)
        plt.ylabel('y (kpc)', fontsize=22)

        #set axis limits
        plt.ylim(-100,100)
        plt.xlim(-100,100)

        #adjust tick label font size
        label_size = 14
        matplotlib.rcParams['xtick.labelsize'] = label_size 
        matplotlib.rcParams['ytick.labelsize'] = label_size

    # Save to a file
    ax.set_rasterized(True)
    plt.savefig(("all_{:0003d}.png".format(i)))
    plt.close()

# Convert the pictures to a gif
print('converting...')
os.system('convert -loop 0 -delay 10 all_*.png all.gif')

# Create individual gifs for each galaxy
for gal in galaxies:
    for i in range( 802):
        # Get name of snapshot for the galaxy
        name = '{:s}_{:0003d}'.format(gal,i)
        print(name)

        # Find time step values
        t = i*14.2857

        # Create plot
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Find COM position
        COMD = CenterOfMass("{}.txt".format(name),2)
        COMP = COMD.COM_P(0.1, 4.0)

        # Recenter to COM position
        xD = COMD.x - float(COMP[0]/u.kpc)
        yD = COMD.y - float(COMP[1]/u.kpc)
        zD = COMD.z - float(COMP[2]/u.kpc)
    
        # Create plot of norma
        plt.plot(xD, yD, 'k.', markersize=1, alpha=0.075)
        plt.title('t={:0.003f} Myr'.format(t), fontsize=22)
                
        # Add axis labels
        plt.xlabel('x (kpc)', fontsize=14)
        plt.ylabel('y (kpc)', fontsize=14)

        #set axis limits
        plt.ylim(-100,100)
        plt.xlim(-100,100)

        #adjust tick label font size
        label_size = 12
        matplotlib.rcParams['xtick.labelsize'] = label_size 
        matplotlib.rcParams['ytick.labelsize'] = label_size

        # Save to a file
        ax.set_rasterized(True)
        plt.savefig(("{:s}.png".format(name)))
        plt.close()

    print('converting...')
    os.system('convert -loop 0 -delay 10 {}_*.png {}.gif'.format(gal, gal))
    print('quick break')
    time.sleep(60)
