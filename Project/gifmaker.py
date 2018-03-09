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

import os
import time

galaxies = ['M33', 'M31', 'MW']

'''
for i in range(802):
    t = 14.2857*i
    
    fig = plt.figure(figsize=(13,10))
    ax = plt.subplot(111)

    COMD = CenterOfMass("MW_{:0003d}.txt".format(i),2)
    COMP = COMD.COM_P(0.1, 4.0)
    
    for gal in galaxies:
        name = "{:s}_{:0003d}".format(gal,i)
        print(name)
        
        COMD = CenterOfMass("{}.txt".format(name),2)
    
        xD = COMD.x - float(COMP[0]/u.kpc)
        yD = COMD.y - float(COMP[1]/u.kpc)
        zD = COMD.z - float(COMP[2]/u.kpc)
    
        #### PLACE FUNCTION HERE ####
        plt.hist2d(xD, yD, bins=250, norm=LogNorm(vmin=0,vmax=3), cmap='magma')
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

print('converting...')
os.system('convert -loop 0 -delay 10 all_*.png all.gif')
'''

for gal in galaxies:
    for i in range(802):
        name = '{:s}_{:0003d}'.format(gal,i)
        print(name)

        t = i*14.2857
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        COMD = CenterOfMass("{}.txt".format(name),2)
        COMP = COMD.COM_P(0.1, 4.0)
        
        xD = COMD.x - float(COMP[0]/u.kpc)
        yD = COMD.y - float(COMP[1]/u.kpc)
        zD = COMD.z - float(COMP[2]/u.kpc)
    
        #### PLACE FUNCTION HERE ####
        plt.hist2d(xD, yD, bins=25, norm=LogNorm(vmin=1e0,vmax=1e3.5), cmap='magma')
        plt.title('t={:0.003f} Myr'.format(t), fontsize=22)
        plt.colorbar()
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
        plt.savefig(("{:s}.png".format(name)))
        plt.close()

    print('converting...')
    os.system('convert -loop 0 -delay 10 {}_*.png {}.gif'.format(gal, gal))
    print('quick break')
    time.sleep(60)
