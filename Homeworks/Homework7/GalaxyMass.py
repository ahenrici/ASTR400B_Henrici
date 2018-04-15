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

