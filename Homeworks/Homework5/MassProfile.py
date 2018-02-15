#!/usr/bin/env python3

'''
File name:  CenterOfMass.py
Author: Andrew Henrici
Date Creaated:  2/6/2018
Date last modified: 2/6/2018
Description:  Calculate the center of mass
position and velocity vectors for chosen galaxies
'''
from CenterOfMass import CenterOfMass
from ReadFile import read
import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt


class MassProfile:
    
    def __init__(self, galaxy, snap):
        # Set the galaxy name as a global variable
        self.gname = galaxy

        # Define Gravitational constant in our units
        self.G = c.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        # add a string of the filenumber to the value “000”
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]

        # create the filename as a global variable
        self.filename = "{:s}_{:s}.txt".format(self.gname, ilbl)
        
        # Read in the file and particle type
        self.time, self.total, self.data = read(self.filename)

        # Store the mass of only the particles of the given type
        self.m = self.data['m']

        # Store the position values of only the particles of the given type
        self.X = self.data['x']*u.kpc
        self.Y = self.data['y']*u.kpc
        self.Z = self.data['z']*u.kpc

        return None


    def MassEnclosed(self, ptype, radii):
        # Finds the mass enclosed in each radius from the list of radii,
        # for a given particle type

        # Find the center of mass for given particle type
        COM = CenterOfMass(self.filename, ptype)
        Gal_COM = COM.COM_P(1)

        # Create an array to store indexes of particles of desired Ptype
        index = np.where(self.data['type'] == ptype)

        # Define the new coordinates and masses
        m_new = self.m[index]
        x_new = self.X[index] - Gal_COM[0]
        y_new = self.Y[index] - Gal_COM[1]
        z_new = self.Z[index] - Gal_COM[2]

        # The magnitude of the distance of the particle from COM
        R = np.sqrt(x_new**2 + y_new**2 + z_new**2)
        
        # Initialize the array that will contain the masses
        mass = np.array([])

        # Loop throught list of radii and find mass inside that radius
        for r in radii:
            # Select all particles within radius r
            list_index = np.where(R < r*u.kpc)

            # Find the total mass within the radius
            m = sum(m_new[list_index])*1e10*u.Msun

            # Add the total mass to list
            mass = np.append(mass, m)

        return mass


    def MassEnclosedTotal(self, radii):
        # Find the mass within the list of radii for each type of particle

        # Initialize array to contain masses
        Mass = np.array([])

        # Variable to sum the masses for each particle
        mass = 0
        
        # Loop through each particle type and add the masses to list
        for p in range(1,4):
            if not (bool(self.gname == 'M33') and bool(p == 3)):    # Ignore the bulge in M33
                mass += self.MassEnclosed(p, radii)
                print(p)
        Mass = np.append(Mass, mass)
        
    
        return Mass


    def HernquistMass(self, radius, a, Mhalo):
        # Calculate the theoretical Hernquist Mass Profile of a halo
        # for given radius, scale factor, and halo mass.
        M = Mhalo*radius**2/(a + radius)**2

        return M
    

    def CircularVelocity(self, ptype, radii):
        # Find the circular velocity at different radii from the mass
        # enclosed of a certain particle type assuming spherical symmetry
        
        # Collect the mass enclosed at different radii
        M = np.array(self.MassEnclosed(ptype, radii))
        
        
        # The circular velocity for mass enclosed
        vel = np.sqrt(self.G*M*u.Msun/(radii*u.kpc)).to(u.km/u.s)
        
        return vel


    def CircularVelocityTotal(self, radii):
        # Find the circular velocity at different radii for all
        # particle types
       
        # Get the total mass profile
        M = self.MassEnclosedTotal(radii)

        # Find the rotation curve for the total mass profile
        Vel = np.sqrt(self.G*M*u.Msun/(radii*u.kpc)).to(u.km/u.s)
       
        return Vel


    def HernquistVCirc(self, radius, a, Mhalo):
        #Find the Hernquist rotation curve for given radius 
        
        # Find the circular velocities 
        M = self.HernquistMass(radius, a, Mhalo)

        # The circular velocity for mass enclosed
        V = np.sqrt(self.G*M*u.Msun/(radius*u.kpc))
        
        return V


# Define array of radii that will be observed
radii = np.linspace(0.1, 30, 600)

# Define dictionary of particle types that will be used
types = {'Halo': 1, 'Disk': 2, 'Bulge': 3}

# Define dectionary of galaxies that will be examined
# and the Hernquist profile scale factors for each galaxy
Galaxies = {'MW': 63, 'M31': 63, 'M33': 25}

# Define dictionary for the total halo masses of each galaxy
Mhalos = {'MW': 1.975, 'M31': 1.921, 'M33': 0.187}

# Set counter for plots
count = 0

print('Mass Profiles:\n')
# Loop throught each galaxy in the list and plot the mass profiles
for Galaxy in ['M33', 'M31', 'MW']:
    # Define the galaxy being used
    Profile = MassProfile(Galaxy, 0)
    # Set up the figure for the galaxy
    plt.figure(count)
    plt.title("{} Mass Profile".format(Galaxy))

    # Loop throught the particle types
    for p in ['Halo', 'Disk', 'Bulge']:
        if not (bool(Galaxy == 'M33') and bool(p == 'Bulge')):
            print(Galaxy, p)
            # Get the mass profile 
            m = Profile.MassEnclosed(types[p], radii)
            
            # Add the mass profile to the plot
            plt.semilogy(radii, m,  label = p)

    # Plot the total mass profile
    print('{} Total mass'.format(Galaxy))
    plt.semilogy(radii, Profile.MassEnclosedTotal(radii), color='k', label = 'Total')

    # Plot the Hernquist profile
    Mhalo = Mhalos[Galaxy]*1e12*u.Msun
    HM_Profile = Profile.HernquistMass(radii, Galaxies[Galaxy], Mhalo)
    plt.semilogy(radii, HM_Profile, color='darkblue', linestyle='--',\
                 label='Hernquist: a={}'.format(Galaxies[Galaxy]))

    # Add legend to bottom right of graph
    plt.legend(loc = 4)

    # Add labels to each axis
    plt.xlabel(r'Radius [$kpc$]')
    plt.ylabel(r'Mass [M$\odot$]')
    
    count += 1
   

print('\n\nRotation Curves:\n')
# Loop throught each galaxy in the list and find the rotation curves
for Galaxy in ['M33', 'M31', 'MW']:
    # Define the galaxy being used
    Profile = MassProfile(Galaxy, 0)

    # Set up the figure for the galaxy
    plt.figure(count)
    plt.title("{} Rotation Curve".format(Galaxy))

    # Loop throught the particle types
    for p in ['Disk', 'Bulge','Halo']:
        if not (bool(Galaxy == 'M33') and bool(p == 'Bulge')):
            print(Galaxy, p)
            # Get the rotation curve for the particle
            v = Profile.CircularVelocity(types[p], radii)

            # Add the rotation curve to the plot
            plt.semilogy(radii, v,  label = p)
 

    # Plot the total mass rotation curve
    print('{} Total Velocity curve'.format(Galaxy))
    Tot_Vel = Profile.CircularVelocityTotal(radii)
    plt.semilogy(radii, Tot_Vel, color='k', label='Total')

    # Calculate the Hernquist rotation curve
    Mhalo = Mhalos[Galaxy]*1e12*u.Msun
    H_Vel = Profile.HernquistVCirc(radii, Galaxies[Galaxy], Mhalo)

    # Plot the Hernquist rotation curve
    plt.semilogy(radii, H_Vel, color='darkblue', linestyle='--',\
                 label='Hernquist: a={}'.format(Galaxies[Galaxy]))

    # Add legend to bottom right of plot
    plt.legend(loc = 4)

    # Add labels to each axis
    plt.xlabel(r'Radius [$kpc$]')
    plt.ylabel(r'Velocity [$m/s$]')
    
    count += 1


plt.show()

