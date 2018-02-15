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

