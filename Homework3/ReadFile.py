#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
File name:  ReadFile.py
Author: Andrew Henrici
Date Creaated:  1/16/2018
Date last modified: 4/25/2018
Description:  Create function to read in data files.
'''

# import modules and packages to run the script
import numpy as np
from astropy import units as u

# function to read file and grab information
def read(filename):
    file = open(filename, 'r')
    
    # Read in the first line and grab value for the time
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*10.0*u.Myr

    # Read in the second line and grab value for the number of particles
    line2 = file.readline()
    label, value = line2.split()
    N_particles = float(value)*10*u.Msun
    file.close()

    # Read in the rest of the file and stores data for each particle
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)

    return time, N_particles, data
