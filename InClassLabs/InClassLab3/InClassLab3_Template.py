#In Class Lab 3
# Template For Lab
# Feb 6 2018

# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib


# LOAD DATA
#****************
 
# Comes from   http://stellar.dartmouth.edu/models/isolf_new.html

# files have been modified from download.  ( M/Mo --> M;   Log L/Lo --> L)
# removed #'s from all lines except column heading

# NOTE SETTINGS USED:  Y = 0.245 default   [Fe/H] = -2.0  alpha/Fe = -0.2
# These could all be changed and it would generate a different isochrone


# Filename for data with Isochrone fit for 1 Gyr
filename1="Isochrone1.txt"


# READ IN DATA
# "dtype=None" means line is split using white spaces
# "skip_header=8"  skipping the first 8 lines 
# the flag "names=True" creates arrays to store the date
#       with the column headers given in line 8 

# Read in data for 1 Gyr
data1 = np.genfromtxt(filename1,dtype=None,names=True,skip_header=8)


# Plot Isochrones 
# For Carina
################

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# PLOT ISOCHRONES

### EDIT THIS ! 
plt.plot(data1['B']-data1['R'], data1['R'], color='blue', linewidth=5, label='1 Gyr')


# Add axis labels
plt.xlabel('B-R', fontsize=22)
plt.ylabel(r'R', fontsize=22)

#set axis limits
plt.xlim(-0.5,2)
plt.ylim(5,-2.5)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.figtext(0.6, 0.15, 'Isochrone Carina', fontsize=22)


# Save to a file
ax.set_rasterized(True)
plt.savefig('IsochroneLabCarina.eps', rasterized=True, dpi=350)





