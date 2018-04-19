'''
Author:  Andrew Henrici
Date Made:
'''
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
from RotateFrame import RotateFrame
# Grab useful functions
from numpy import pi, sqrt


class binmaker:
    
    def __init__(self, filename):
        # Store the COM of the Milky Way to center on
        COMD = CenterOfMass(filename,2)
        COMP = COMD.COM_P(0.1, 4.0)

        # Recenter to COM position
        self.xD = COMD.x - float(COMP[0]/u.kpc)
        self.yD = COMD.y - float(COMP[1]/u.kpc)
        self.zD = COMD.z - float(COMP[2]/u.kpc)
        
    # Transform coordinates from cartesian to spherical
    def cart2sph(self, x, y, z):
        hxy = np.hypot(x, y)
        r = np.hypot(hxy, z)
        phi = np.arctan2(z, hxy)
        theta = np.arctan2(y, x)
        return r, theta, phi

       
    def cart2cyl(self, x, y, z):
        rho = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y,x)
        z = z
        return rho, theta, z

    def rmax(self):
        return np.amax(self.rS)
        
    def Radial_Bins(self, N):
        # Number of bins
        N_radial = N

        # Convert COM to spherical
        self.rS, self.thetaS, self.phiS = self.cart2cyl(self.xD, self.yD, self.zD)
        
        r_max = self.rmax()
        
        # Create list of radii that will divide the bins
        radii = np.linspace(0, r_max, N_radial+1)
        radial_bins = []

        
        # Sort the particles into each bin
        for r in range(N_radial):
            radial_bins.append(np.where((self.rS <= radii[r+1]) & (self.rS > radii[r])))

        return radial_bins, r_max


    def Angular_Bins(self, N):
        # Number of bins + 1
        N_angular = N + 1
        # Create list of angles that will divide the bins
        angles = np.linspace(-pi, pi, N_angular)
        angular_bins = []

        # Convert COM to spherical
        self.rS, self.thetaS, self.phiS = self.cart2cyl(self.xD, self.yD, self.zD)
        
        # Sort the particles into each bin
        for a in range(N_angular-1):
            angular_bins.append(np.where((self.thetaS <= angles[a+1]) & (self.thetaS > angles[a])))

     