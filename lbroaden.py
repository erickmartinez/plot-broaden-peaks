# -*- coding: utf-8 -*-
"""
Created on Mon May  1 09:27:35 2017

@author: erickrmartinez@gmail.com
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
#from scipy import constants
import re
import os

# Enter the file with the peak data points
inputFileName = r'proof-tddft.dat'

FWHM = 2.5
DX = 0.01
MIN_FACTOR = 2.0


# Get the file name without extention
baseFileName = os.path.splitext(inputFileName)[0]
# Get the current path
cwd = os.path.dirname(os.path.abspath(__file__))
# Get the full path to the data file:
path2DataFile = os.path.join(cwd,inputFileName)
# Create a CSV file
outputFileName = baseFileName + '.csv'
# Get the full path to the output file
pathToOutput = os.path.join(cwd,outputFileName)
# Get the full path for the output file
pathToFigure = os.path.join(cwd,baseFileName)
# Define a data type which will be used to import the peaks
lbType = np.dtype([('energy','d'),('lambda','d'),('peak_y','d')])
# Define a data type to export the spectrum
spType = np.dtype([('x','d'),('y','d')])

outputHeader = "Wavelength (nm), Intensity (a.u.)"

peak_wl = []
peak_y = []


# Construct a regular expression to match the data in the file
# Regular expression detail:
# Group 1 matches the description of the state
reg_group = [None]*7
reg_group[0] = '(^\s*Excited State\s*\d*:\s*[a-zA-Z-]*\s*)'
# Group 2 matches the energy of the state
reg_group[1] = '(\d*\.\d*)'
# Group 3 matches the units of the energy column
reg_group[2] = '(\s*eV\s*)'
# Group 4 matches the wavelegth
reg_group[3] = '(\d*\.\d*)'
# Group 5 matches the units of the peak_wl and 'f='
reg_group[4] = '(\s*nm\s*f\=)'
# Group 6 matches the peak_y of the peak
reg_group[5] = '(\d*\.\d*)'
# Group 7 matches any character at after this
reg_group[6] = '(.*$)'
# Join the groups in a single pattern
pattern = ''.join(reg_group)
# Read the file
with open(path2DataFile, 'r') as infile:
    for line in infile:
        m = re.match(pattern,line)
        if m:
            peak_wl.append(float(m.group(4)))
            peak_y.append(float(m.group(6)))
#print("File closed? %s " % infile.closed)

peak_wl = np.array(peak_wl)
peak_y = np.array(peak_y)

LMIN = np.amin(peak_wl) - MIN_FACTOR*FWHM
LMAX = np.amax(peak_wl) + MIN_FACTOR*FWHM
steps = int((LMAX - LMIN)/DX)
gamma = 0.5 * FWHM
# Area normalization factor (can be set to 1 to normalize to height)
factor = 1.0 / (gamma * np.pi)
peaks = len(peak_wl)
print("Peaks in file: %d" % peaks)
print("dx : %.5f" % DX)
print("LMIN : %.3f" % LMIN)
print("LMAX : %.3f" % LMAX)
print("FWHM : %.3f" % FWHM)
print("gamma : %.3f" % gamma)
    
spectrum = np.array([],dtype=spType)

# Iterate over all the steps to get the y at each x
for i in range(steps):
    x = i*DX + LMIN
    y = 0.0
    # iterate over all the peaks and sum the contribution from each one
    for j in range(peaks):
        squarex = x - peak_wl[j]
        squarex = squarex / gamma
        squarex = squarex*squarex
        sumi = 1 + squarex
        y += peak_y[j] / sumi
    # Add the (x,y) data point to the spectrum
    row = np.array((x,y),dtype=spType)
    spectrum = np.append(spectrum,row)

# Save the broaden spectrum
np.savetxt(pathToOutput,
            spectrum,
            delimiter=',',
            fmt=('%.3f',
                    '%.4e'),
            header=outputHeader,
            comments='',)

# Plot the results
golden_ratio = (np.sqrt(5)+1)/2
asect_ratio = golden_ratio**2
# Plot style parameters
plotStyle = {'font.size': 14,
            'legend.fontsize': 10,
            'mathtext.fontset': 'custom',
            'mathtext.rm': 'Arial',
            'mathtext.it': 'Arial:italic',
            'mathtext.bf': 'Arial:bold',
            'xtick.major.size' : 6,
            'xtick.major.width' : 1.5,
            'ytick.major.size' : 6,
            'ytick.major.width' : 1.5,
            'xtick.minor.size' : 3,
            'xtick.minor.width' : 1.0,
            'ytick.minor.size' : 3,
            'ytick.minor.width' : 1.0,
            'figure.dpi': 100}
mpl.rc('font',family='Arial')
mpl.rc('axes', linewidth=2.5)
mpl.rcParams.update(plotStyle)
# Width of the spectrum
linew = 2.0


plt.close('all')
fig = plt.figure()
fig.set_size_inches(8.1,4.5,forward=True)
ax1 = plt.subplot2grid((1,1), (0,0))
ax1.autoscale_view(True,True,True)

# Plot the broaden spectrum
broaden_plot = ax1.plot(spectrum['x'],spectrum['y'],label="Broaden",color='C0',lw=linew)
# Plot the peaks
peaks_plot = ax1.bar(peak_wl,peak_y,label="Peaks",color='C1',lw=linew)

xlabel = "Wavelength (nm)"
ylabel = "Intensity (a.u.)"
ax1.set_xlabel(xlabel,fontweight='bold',fontsize=18)
ax1.set_ylabel(ylabel,fontweight='bold',fontsize=18)

# Remove ticks from the y axis
ax1.get_yaxis().set_ticks([])
# Eyecandy for the x axis
ax1.xaxis.set_major_locator(mticker.MaxNLocator(8,prune='lower'))
ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))

#leg = ax1.legend(loc='upper right',prop={'family':'Arial','size':14},frameon=False)

# Change the aspect ratio of the plot
x1,x2 = ax1.axes.get_xlim()
y1,y2 = ax1.axes.get_ylim()
ratio = (x2 - x1) / (y2 - y1)
ax1.set_aspect(ratio/asect_ratio)
            
plt.tight_layout()
plt.show()

# Save the plot in three different formats
fig.savefig(pathToFigure + '.png' , dpi=300)
fig.savefig(pathToFigure + '.eps', format='eps', dpi=600)
fig.savefig(pathToFigure + '.pdf', format='pdf', dpi=600)