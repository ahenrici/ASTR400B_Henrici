#!/usr/bin/env python3

'''
File name:  GalaxyMass.py
Author: Andrew Henrici
Date Creaated:  1/28/2018
Date last modified: 1/30/2018
Description:  Calculate the total mass for each component
of the galaxies in the Local Group
'''

from ReadFile import read
from astropy import units as u
from numpy import around, where, array
from astropy.table import Table, Column



# Function to find the total mass of a given component
def ComponentMass(FileName, PType):
    # Read the data from the file
    time, N_particles, data = read(FileName)

    # Selects data for particles of chosen type
    index = where(around(data['type'], 3) == PType)
    Particle = data[index]

    # Finds only the mass information
    # converts to proper units and sums
    mass = Particle['m']*10**10*u.Msun
    Mass = around(sum(mass), 3)
    
    return Mass

# Dictionary for the particle types
types = {'Dark Matter': 1.000, 'Disk': 2.000, 'Bulge': 3.000}

# List of the galaxies that will be looked at
Galaxies = ['MW', 'M31', 'M33']

# Set up list to collect all data for output
Table_arr = []

# Cycles through all of the galaxies in the list
for Gal in Galaxies:
    
    # Set up list to collect all data per galaxy
    Gal_Array = []

    # Set up variables to 
    Total = 0
    Baryon_tot = 0
    
    # Cycles through all of the particle types
    # and finds the totla mass
    for t in types:
        
        # Set up file name
        fnum = '000'
        fname = '{}_{}'.format(Gal, fnum)

        # Find the mass for selected component
        Mass = ComponentMass('{}.txt'.format(fname), types[t])

        # Collect the component information for later output
        Gal_Array.append([Gal, t, Mass])
        
        # Find the total mass for the galaxy
        Total += Mass
        

        # Find the total baryonic mass og the galaxy
        if t != 'Dark Matter':
            Baryon_tot += Mass

    # Find the baryonic ratio of the galaxy
    fbar = Baryon_tot/Total

    # Add galaxy information to output list
    print(Gal_Array, Total, fbar)


