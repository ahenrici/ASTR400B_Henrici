#!/usr/bin/env python3

'''
File name:  VLowRes.py
Author: Andrew Henrici
Date Created:  2/17/2018
Date last modified: 2/6/2018
Description:  Calculate the center of mass
position and velocity vectors for chosen galaxies
'''
import numpy as np
import astropy.units as u
from ReadFile import read
from CenterOfMass import CenterOfMass
import matplotlib.pyplot as plt


def OrbitCOM(galaxy, start, end, n):
    fileout = "Orbit_{}.txt".format(galaxy)

    Orbit = np.zeros((int(end/n + 1),7))

    # Define the values for COM
    delta = 0.1
    VolDec = 2

    for i in np.arange(start, end+n, n):
        # add a string of the filenumber to the value “000”
        ilbl = '000' + str(i)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        
        # create the filename
        filename = "{:s}_{:s}.txt".format(galaxy, ilbl)


        # Create a Center of Mass object for the Galaxy
        COM = CenterOfMass(filename, ptype = 2)
        Gal_COM = COM.COM_P(delta, VolDec)
        Gal_COMV = COM.COM_V(Gal_COM)
        
        # Get the time, position, and velocities of the particles
        time  = float(COM.time/u.Myr)/1000.0
        X = float(Gal_COM[0]/u.kpc)
        Y = float(Gal_COM[1]/u.kpc)
        Z = float(Gal_COM[2]/u.kpc)
        Vx = float(Gal_COMV[0]/(u.km/u.s))
        Vy = float(Gal_COMV[1]/(u.km/u.s))
        Vz = float(Gal_COMV[2]/(u.km/u.s))

        # Store the values in the Orbit array
        Orbit[int(i/n),] = np.array([time, X, Y, Z, Vx, Vy, Vz])

        print("{} {}".format(galaxy, i))

    np.savetxt(fileout, Orbit, header='t x y z vx vy vz', comments='# ',\
               fmt=['%.2f ', '%.2f ','%.2f ','%.2f ','%.2f ','%.2f ','%.2f '])
            
    return None


for Galaxy in ['M33', 'M31', 'MW']:
    OrbitCOM(Galaxy, 0, 800, 5)

M33 = np.genfromtxt('Orbit_M33.txt', dtype = float, names=True)
M31 = np.genfromtxt('Orbit_M31.txt', dtype = float, names=True)
MW = np.genfromtxt('Orbit_MW.txt', dtype = float, names=True)


fig = plt.figure(0)
plt.suptitle("Center of Mass Position")

plt.subplot(313)
plt.xlabel('Time')
plt.ylabel('X position')
plt.plot(MW['t'], MW['x'], label='MW', color='red')
plt.plot(M33['t'], M33['x'], label='M33', color='green')
plt.plot(M31['t'], M31['x'], label='M31', color='blue')

plt.subplot(312)
plt.xlabel('Time')
plt.ylabel('Y position')
plt.plot(MW['t'], MW['y'], label='MW', color='red')
plt.plot(M33['t'], M33['y'], label='M33', color='green')
plt.plot(M31['t'], M31['y'], label='M31', color='blue')

plt.subplot(311)
plt.xlabel('Time')
plt.ylabel('X position')
plt.plot(MW['t'], MW['z'], label='MW', color='red')
plt.plot(M33['t'], M33['z'], label='M33', color='green')
plt.plot(M31['t'], M31['z'], label='M31', color='blue')
plt.legend()

plt.show()
