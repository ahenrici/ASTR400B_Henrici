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
from astropy.coordinates import cartesian_to_spherical, spherical_to_cartesian

def cart2sph(x, y, z):
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    phi = np.arctan2(z, hxy)
    theta = np.arctan2(y, x)
    return r, theta, phi


def sph2cart(r, theta, phi):
    rcos_theta = r * np.cos(phi)
    x = rcos_theta * np.cos(theta)
    y = rcos_theta * np.sin(theta)
    z = r * np.sin(phi)
    return x, y, z

def cart2cyl(x, y, z):
    rho = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)
    z = z
    return rho, theta, z

def cyl2cart(rho, theta, z):
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = z
    return x, y, z


# Store the COM of the Milky Way to center on
COMD = CenterOfMass("MW_000.txt",2)
COMP = COMD.COM_P(0.1, 4.0)

# Recenter to COM position
xD = COMD.x - float(COMP[0]/u.kpc)
yD = COMD.y - float(COMP[1]/u.kpc)
zD = COMD.z - float(COMP[2]/u.kpc)

# Convert COM to spherical
#rD, thetaD, phiD = cart2sph(xD, yD, zD)
rD, thetaD, zD = cart2cyl(xD, yD, zD)

# Number of bins + 1
N_angular = 5
# Create list of angles that will divide the bins
angles = np.linspace(-pi, pi, N_angular)
angular_bins = []

# Sort the particles into each bin
for a in range(N_angular-1):
    angular_bins.append(np.where((thetaD <= angles[a+1]) & (thetaD > angles[a])))

# Number of bins
N_radial = 10
# Create list of radii that will divide the bins
radii = np.linspace(0, 40, N_radial+1)
radial_bins = []

# Sort the particles into each bin
for r in range(N_radial):
    radial_bins.append(np.where((rD <= radii[r+1]) & (rD > radii[r])))


fig = plt.figure()

for nr in range(N_radial):

    rbin=radial_bins[nr]
    R = rD[rbin]
    Theta = thetaD[rbin]
    Z = zD[rbin]

    for na in range(N_angular-1):
        abin = np.where((Theta <= angles[na+1]) & (Theta > angles[na]))

        r = R[abin]
        theta = Theta[abin]
        z = Z[abin]

        # Plot 
        x, y, z = cyl2cart(r, theta, z)
        plt.plot(x,y,'.', markersize=1)

plt.show()

