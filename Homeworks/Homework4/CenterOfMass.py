#!/usr/bin/env python3

'''
File name:  CenterOfMass.py
Author: Andrew Henrici
Date Creaated:  2/6/2018
Date last modified: 2/6/2018
Description:  Calculate the center of mass
position and velocity vectors for chosen galaxies
'''
import numpy as np
import astropy.units as u
from ReadFile import read


class CenterOfMass:

    def __init__(self, filename, ptype):
        # Read in the file and particle type
        self.time, self.total, self.data = read(filename)

        # Create an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)

        # Store the mass of only the particles of the given type
        self.m = self.data['m'][self.index]*1e10*u.Msun

        # Store the position values of only the particles of the given type
        self.X = self.data['x'][self.index]*u.kpc
        self.Y = self.data['y'][self.index]*u.kpc
        self.Z = self.data['z'][self.index]*u.kpc

        # Store the velocity values of only the particles of the given type
        self.Vx = self.data['vx'][self.index]*u.km/u.s
        self.Vy = self.data['vy'][self.index]*u.km/u.s
        self.Vz = self.data['vz'][self.index]*u.km/u.s
        
        return None

    def COMdefine(self, M, x, y, z):
        # Finds the center of mass for the x,y,z coordinates with the particle
        # Input the list of the particles values of mass and 3D vector
        COM = [sum(M*x)/sum(M),
               sum(M*y)/sum(M),
               sum(M*z)/sum(M)]
        return COM

    def COM_P(self, delta):
        # Find the magnitude of the position vector for the center of mass
        # that has converges below given tolerance values (delta)

        # Declare initial variables 
        xi = self.X
        yi = self.Y
        zi = self.Z
        mi = self.m

        # Find the initial center of mass for allparticles
        COM = self.COMdefine(mi, xi, yi, zi)

        # Calculate the magnitude of the initialguess Center of Mass position
        RCOM = np.sqrt(COM[0]**2 + COM[1]**2 + COM[2]**2)

        # Change the coordinates to be around the new center of mass
        X_new = xi - COM[0]
        Y_new = yi - COM[1]
        Z_new = zi - COM[2]
        M_new = mi
        
        # New magnitude of the paricle positions
        RNEW = np.sqrt(X_new**2 + Y_new**2 + Z_new**2)

        # Find half of the distance of the farthest particle
        RMAX = max(RNEW)/2

        # Declare initial difference for loop
        diff = 1e5*u.kpc

        # While loop to hone in on the center of mass
        while diff > delta*u.kpc:
            # Find the list of particles in the new limit
            new_list_index = np.where(RNEW < RMAX)

            # Limit the particles to new list
            M_new = self.m[new_list_index]
            X_new = self.X[new_list_index]
            Y_new = self.Y[new_list_index]
            Z_new = self.Z[new_list_index]

            # Find the new center of mass
            COM2 = self.COMdefine(M_new, X_new, Y_new, Z_new)

            # Find the new magnitude of the center of mass
            RCOM2 = np.sqrt(COM2[0]**2 + COM2[1]**2 + COM2[2]**2)
            
            # Shift the position to new center of mass
            X_new = self.X - COM2[0]
            Y_new = self.Y - COM2[1]
            Z_new = self.Z - COM2[2]

            # New magnitude of the particle positions
            RNEW = np.sqrt(X_new**2 + Y_new**2 + Z_new**2)

            # Find half of the distance of the farthest particle
            RMAX = RMAX/2
            
            # Calculate the difference between the new and old center of masses
            diff = RCOM - RCOM2
            
            # Store COM and RCOM 
            RCOM = RCOM2
            COM = COM2
            
        return COM


    def COM_V(self, COM):
        # Finds the COM velocity vector after taking in the COM position coordinates

        # Shift the particles to the COM pos
        X_new = self.X - COM[0]
        Y_new = self.Y - COM[1]
        Z_new = self.Z - COM[2]

        # Find the magnitue of the distance of the particles from the COM
        RNEW = np.sqrt(X_new**2 + Y_new**2 + Z_new**2)

        # Narrow particles to within 15 kpc from COM
        new_list_index = np.where(RNEW < 15*u.kpc)
        
        # Select the velocity vectors and masses of particles within 15 kpc
        # from the COM
        Vxi = self.Vx[new_list_index]
        Vyi = self.Vy[new_list_index]
        Vzi = self.Vz[new_list_index]
        mi = self.m[new_list_index]
        
        # Find the initial center of mass for allparticles
        COMV = self.COMdefine(mi, Vxi, Vyi, Vzi)

        return COMV

# Output:
# Finds the position and velocity vectors for input galaxies


# Create a Center of mass object for the MW
MWCOM = CenterOfMass("MW_000.txt", 2)

# Calculate COM position for MW data
MW_COM = MWCOM.COM_P(1e-2)
print('\nThe COM position of the Milky Way is:\n', MW_COM)

# Calculate COM velolcity for the MW data
MW_COMV = MWCOM.COM_V(MW_COM)
print('The COM velocity of the Milky Way is:\n', MW_COMV)


# Create a Center of mass object for the M31
M31COM = CenterOfMass("M31_000.txt", 2)

# Calculate COM position for M31 data
M31_COM = M31COM.COM_P(1e-2)
print('\nThe COM position of Andromeda is:\n', M31_COM)

# Calculate COM velocity for M31 data
M31_COMV = M31COM.COM_V(M31_COM)
print('The COM velocity of Andromeda is:\n', M31_COMV)


# Create a Center of mass object for the M33
M33COM = CenterOfMass("M33_000.txt", 2)

# Calculate COM position for M33 data
M33_COM = M33COM.COM_P(1e-2)
print('\nThe COM position of M33 is:\n', M33_COM)

# Calculate COM velocity for M33 data
M33_COMV = M33COM.COM_V(M33_COM)
print('The COM velocity of M33 is:\n', M33_COMV)

# Find the distance between the MW and M31
R = np.sqrt((MW_COM[0] - M31_COM[0])**2\
            + (MW_COM[1] - M31_COM[1])**2\
            + (MW_COM[2] - M31_COM[2])**2)
print('\nThe distance between MW and M31 is:  ', R)

# Find the magnitude of the velocity between MW and M31
V = np.sqrt((MW_COMV[0] - M31_COMV[0])**2\
            + (MW_COMV[1] - M31_COMV[1])**2\
            + (MW_COMV[2] - M31_COMV[2])**2)
print('The speed between MW and M31 is:  ', V)

# Find the distance between the M31 and M33
R2 = np.sqrt((M31_COM[0] - M33_COM[0])**2\
            + (M31_COM[1] - M33_COM[1])**2\
            + (M31_COM[2] - M33_COM[2])**2)
print('\nThe distance between M31 and M33 is:  ', R2)

# Find the magnitude of the velocity between M31 and M33
V2 = np.sqrt((M31_COMV[0] - M33_COMV[0])**2\
            + (M31_COMV[1] - M33_COMV[1])**2\
            + (M31_COMV[2] - M33_COMV[2])**2)
print('The speed between M31 and M33 is:  ', V2)
