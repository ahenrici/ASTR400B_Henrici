import numpy as np
from astropy import units as u


def Read(filename):
    # Read in the first line and grab value for the time
    file = open(filename, 'r')
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
