#!/usr/bin/env python3

'''
File name:  Homework7.py
Author: Andrew Henrici
Date Created:  3/22/2018
Date last modified: 3/28/2018
Description:  Integrate the orbit of M33 around M31
'''

import numpy as np
import astropy.units as u
import  astropy.constants as const
from ReadFile import read
from CenterOfMass import CenterOfMass
import matplotlib.pyplot as plt
from GalaxyMass import ComponentMass


class M33AnalyticOrbit:

    def __init__(self, filename):
        # Set the name of the file for the output
        self.filename = filename

        # Define the Gravitational constant
        self.G = 4.498768e-6*u.kpc**3/u.Msun/u.Gyr**2

        # Create a Center of Mass object for M31
        M31_COM = CenterOfMass("M31_000.txt", ptype=2)
        M31_COMP = M31_COM.COM_P(0.1, 2)
        M31_COMV = M31_COM.COM_V(M31_COMP)

        # Create a Center of Mass object for M33
        M33_COM = CenterOfMass("M33_000.txt", ptype=2)
        M33_COMP = M33_COM.COM_P(0.1, 2)
        M33_COMV = M33_COM.COM_V(M33_COMP)

        # Recenter to COM position of M33 centered on M31
        self.x = M31_COMP[0] - M33_COMP[0]
        self.y = M31_COMP[1] - M33_COMP[1]
        self.z = M31_COMP[2] - M33_COMP[2]

        # Recenter to COM velocity with respect to M31
        self.vx = M31_COMV[0] - M33_COMV[0]
        self.vy = M31_COMV[1] - M33_COMV[1]
        self.vz = M31_COMV[2] - M33_COMV[2]

        # Set the scale lengths of the disk, bulge, and halo of M31
        self.rd = 5.0*u.kpc
        self.rbulge = 1.0*u.kpc
        self.rhalo = 63.0*u.kpc

        # Set Masses for the disk, bulge, and halo
        self.Md = ComponentMass('M31_000.txt', 2)
        self.Mbulge = ComponentMass('M31_000.txt', 3)
        self.Mhalo = ComponentMass('M31_000.txt', 1)

        return None


    # Calculate acceleration for halo and bulge particles
    # in a given direction using the Hernquist profile
    def HernquistAccel(self, M, ra, x, y, z, direction):
        # Takes in the compnenet mass, scale length of the component,
        # position of M33, and chosen direction and returns
        # the acceleration of that direction

        r = np.sqrt(x**2 + y**2 + z**2)

        # Use the Hernquist potential acceleration
        a = -self.G*M/(r*(ra + r)**2)

        # Calculate accelaration of particles in chosen direction
        if direction == 0:
            a_comp = a*x
        elif direction == 1:
            a_comp = a*y
        elif direction == 2:
            a_comp = a*z

        return a_comp


    # Calculate the acceleration from the disk particles
    # in a given direction using Miaymoto-Nagai profile
    def MiyamotoNagaiAccel(self, M, rd, x, y, z, direction):
    	# Takes in the mass of the disk, the scale length, the 
    	# position coordinates, and chosen direction and returns
    	# the acceleration in that direction

        r = np.sqrt(x**2 + y**2 + z**2)
        zd = rd/5.0
        B = rd + np.sqrt(z**2 + zd**2)

        # Use the Miyamoto Nagai potential acceleration
        a = -self.G*M/(r**2 + B**2)**1.5

        # Calculate accelaration of particles in chosen direction
        if direction == 0:
        	a_comp = a*x
        elif direction == 1:
            a_comp = a*y
        elif direction == 2:
            a_comp = a*B*z/np.sqrt(z**2 + zd**2)

        return a_comp


    # Sums the acceleration terms for each particle
    def M31Accel(self, x, y, z, direction):
    	# Take in the position coordinates and find the 
    	# acceleration for the given direction

    	# Calculate the acceleration for each particle type
        Disk_accel = self.MiyamotoNagaiAccel(self.Md, self.rd, x, y, z, direction)
        Bulge_accel = self.HernquistAccel(self.Mbulge, self.rbulge, x, y, z, direction)
        Halo_accel = self.HernquistAccel(self.Mhalo, self.rhalo, x, y, z, direction)

        # Add together all of the accelerations to one term
        Total_accel = (Disk_accel + Bulge_accel + Halo_accel).to(u.km/u.s**2)

        return Total_accel


    # Use Leap Frog method to calculate delta pos and delta vel
    def LeapFrog(self, dt, x, y, z, vx, vy, vz):
    	# Take in the step size, the position coordinates, and
    	# velocity components and return the updated position
    	# and velocities

        # Do a half step guess of the new positions
        x_half = x + (vx*dt/2.0).to(u.kpc)
        y_half = y + (vy*dt/2.0).to(u.kpc)
        z_half = z + (vz*dt/2.0).to(u.kpc)

        # Do a guess of the acceleration from the half step
        ax_half = self.M31Accel(x_half, y_half, z_half, 0)
        ay_half = self.M31Accel(x_half, y_half, z_half, 1)
        az_half = self.M31Accel(x_half, y_half, z_half, 2)

        # Update the velocity
        vx_step = vx + (ax_half*dt).to(u.km/u.s)
        vy_step = vy + (ay_half*dt).to(u.km/u.s)
        vz_step = vz + (az_half*dt).to(u.km/u.s)

        # Update the Position
        x += ((vx + vx_step)*dt/2.0).to(u.kpc)
        y += ((vy + vy_step)*dt/2.0).to(u.kpc)
        z += ((vz + vz_step)*dt/2.0).to(u.kpc)

        # Store the position and velocities for the return statement
        pos = [x, y, z]
        vel = [vx_step, vy_step, vz_step]

        return pos, vel


    # Use LeapFrog function to integrate over a time period between
    # ti(start) and tf(end) with steps dt
    def OrbitIntegrator(self, ti, tf, dt):
    	# Takes in the start, finish, and step times and saves the 
    	# output to a .txt file

    	# Set up the position values so they can be updated
        x = self.x
        y = self.y
        z = self.z

        # Set up the velocity values so they can be updated
        vx = self.vx
        vy = self.vy
        vz = self.vz

        # Initialize array to store the time, postion, and velocity values
        storage = np.zeros((int((tf/dt)+2), 7))
        storage[0,:] = np.array([ti/u.Gyr, x/u.kpc, y/u.kpc, z/u.kpc, 
        	vx*u.s/u.km, vy*u.s/u.km, vz*u.s/u.km])

        t = ti  # initialize the time 
        c = 1   # initialize the counter

        # Loop through the time period integrate the orbit of M33 around M31
        while t <= tf:
        	print(t)

        	# Update the position and velocities
        	pos, vel = self.LeapFrog(dt, x, y, z, vx, vy, vz)
        	
        	# Reset the variables to the updated values
        	x, y, z = pos[0], pos[1], pos[2]
        	vx, vy, vz = vel[0], vel[1], vel[2]

            # Update the time 
        	t = np.around(t + dt, 3)
        	
        	# Store the updated values for future
        	storage[c,] = np.array([t/u.Gyr, x/u.kpc, y/u.kpc, z/u.kpc, 
        	                vx*u.s/u.km, vy*u.s/u.km, vz*u.s/u.km])
        	
        	# Update the counter
        	c += 1
 
        # Save the output of the orbit to
        np.savetxt(self.filename, storage, header='t x y z vx vy vz', comments='# ',\
			fmt=['%.2f ', '%.2f ','%.2f ','%.2f ','%.2f ','%.2f ','%.2f '])

        return None


# Run the code for now to 10 billion years into the future with
# a time step of 100 million years
Orbit = M33AnalyticOrbit("M31M33_Orbit.txt")
Orbit.OrbitIntegrator(0.0*u.Gyr, 10.0*u.Gyr, 0.1*u.Gyr)

# Read in the results of the program
OrbitInt = np.genfromtxt("M31M33_Orbit.txt", dtype = float, names=True)
M33 = np.genfromtxt('Orbit_M33.txt', dtype = float, names=True)
M31 = np.genfromtxt('Orbit_M31.txt', dtype = float, names=True)

# Compute the distance between the galaxies from this homework
r = np.sqrt(OrbitInt['x']**2 + OrbitInt['y']**2 + OrbitInt['z']**2)

# Compute the distances between the galaxies from the previous homework
R = np.sqrt((M31['x']-M33['x'])**2 + (M31['y']-M33['y'])**2+(M31['z']-M33['z'])**2)

# Plot the two different methods on the same graph
plt.plot(M31['t'], R, label='HW6')
plt.plot(OrbitInt['t'], r, label='HW7')

# Set up plot
plt.title('Distance of M33 from M31')
plt.xlabel('Time [Gyr]')
plt.ylabel('Distance [kpc]')
plt.legend()

# Show plot
plt.show()