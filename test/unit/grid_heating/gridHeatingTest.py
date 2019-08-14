#!usr/bin/env python3

'''  This script tests the heating performance of VPIC.  It assumes you are
using one of two decks, gridHeatingTest1D.cxx or gridHeatingTest2D.cxx, and
relies on some very specifically formated output therein.

VPIC does not have the best heating rates because of low particle shapes, and
this script checks if it matches what is expected, not if it is close enough to
physical reality.

Written by Scott V. Luedtke, XCP-6, August 1, 2019'''

import sys
import numpy as np
import scipy.constants as const
import os

base_dir = ""
if len(sys.argv) > 1:
    base_dir = sys.argv[1]

try:
    # TODO: presumably it only makes sense to use this if it has value 1-3?
    # We should do some checking on it
    dim = int(sys.argv[2])
except:
    dim = 2

# Read in simulation specific parameters from a file
params_path = os.path.join(base_dir, "params_" + str(dim) + "d.txt")
print("Param path is " + params_path)
params = open(params_path, 'r')

params.readline()
dt = float(params.readline().split()[0])
lamb = float(params.readline().split()[0])
ne_over_nc = float(params.readline().split()[0])
nx = int(params.readline().split()[0])
ny = int(params.readline().split()[0])
nz = int(params.readline().split()[0])
box_size_x = float(params.readline().split()[0])
box_size_y = float(params.readline().split()[0])
box_size_z = float(params.readline().split()[0])
field_interval = int(params.readline().split()[0])
tracer_interval = int(params.readline().split()[0])
numtrace = int(params.readline().split()[0])
numvars = int(params.readline().split()[0])
numspecs = int(params.readline().split()[0])
nstep_total = int(params.readline().split()[0])
nprocs = int(params.readline().split()[0])
maxEelec = float(params.readline().split()[0])*const.c**2*const.m_e/const.value('electron volt')*1e-6
maxEI2 = float(params.readline().split()[0])*const.c**2*const.m_p*12./const.value('electron volt')*1e-6
numbins = int(params.readline().split()[0])
spec_interval = int(params.readline().split()[0])

# Get density normalization
omega = 2.*np.pi*const.c/lamb
ncr = const.epsilon_0 * const.m_e * omega**2 / const.e**2
delta = lamb/(np.sqrt(ne_over_nc)*(2.*np.pi))
num_physical_per_macro = delta**3 * ne_over_nc*ncr

# Convert to SI, and don't blame me if some of these are wrong
wpe = 2.*np.pi*const.c*np.sqrt(ne_over_nc)/lamb # Plasma frequency in Hz
timeToSI = 1./wpe
lengthToSI = lamb/(2*np.pi*np.sqrt(ne_over_nc))
efieldToSI = const.m_e*np.sqrt(ne_over_nc)*2.*np.pi*const.c*const.c/(lamb*const.e)
energyToSI = efieldToSI**2 * lengthToSI**3 * const.epsilon_0

filename = os.path.join(base_dir,'rundata_' + str(dim) + 'd/energies')

print("Reading " + filename)

d = np.loadtxt(filename, comments='%')
x = d[:,0]*dt*timeToSI
ex = d[:,1]*energyToSI
ey = d[:,2]*energyToSI
ez = d[:,3]*energyToSI
bx = d[:,4]*energyToSI
by = d[:,5]*energyToSI
bz = d[:,6]*energyToSI
I2 = d[:,7]*energyToSI
elec = d[:,8]*energyToSI

fields = ex + ey + ez + bx + by + bz

# Set this to 1 if you want visual output
draw_plot = 0
if draw_plot:
    import matplotlib.pyplot as plt
    plt.plot(x,fields/elec[0], 'r-', label='Fields')
    plt.plot(x,elec/elec[0], 'b-', label='Electrons')
    plt.plot(x,I2/elec[0], 'g-', label='Protons')
    #plt.plot(x,elec+I2, 'k-', label='Everything')
    plt.legend(loc='upper left')
    #title = "$\Delta x = 0.5c/\omega_p$ 64ppc"
    title = "Grid Heating"
    ylabel = "Normalized Energy"
    xlabel = "Time (s)"
    plt.title(title, fontsize=25)
    plt.ylabel(ylabel, fontsize=20)
    plt.xlabel(xlabel, fontsize=20)
    #plt.yscale('log', nonposy='clip')
    plt.show()

volume = box_size_x*box_size_y*box_size_z*1e-6**3 # SI
start = 400
stdpass = 5

coefs = np.polyfit(x[start:],elec[start:], 1)
#print("Electron coefs are ", coefs)
elecrate = coefs[0]/volume
if dim==2:
    avg = 1.34586133288e+28
    std = 1.819809542e+26
elif dim==1:
    avg = 4.06512504887e+27
    std = 6.42447884639e+25
else:
    print("Warning, no comparison data found; using fake data")
    avg = 1.
    std = 1.
err = elecrate-avg
errstd = abs(err/std)
print("The electron heating rate is", elecrate, "Joules per second per cubic meter.")
print("Five runs during the creation of this test had an average of", avg, "and standard deviation of", std, "J/(s m^3)")
print("This run is off by", err, "which is", errstd, "standard deviations.")
print("This test will pass if within", stdpass, "standard deviations.\n")
if errstd < stdpass:
    elecpass = True
    print("Electron heating rate test PASS.\n")
else:
    elecpass = False
    print("Electron heating rate test FAIL.\n")
    sys.exit(1);


print("The ions in this test do not get into a linear heating regime, and vary much more in their heating rate than the electrons.  This test will always pass.")
coefs = np.polyfit(x[start:],I2[start:], 1)
protrate = coefs[0]/volume
if dim == 2:
    pavg = 2.88450472127e+25
    pstd = 1.70730720784e+24
elif dim==1:
    pavg = 6.22504867595e+25
    pstd = 4.21707647083e+24
else:
    print("Warning, no comparison data found; using fake data")
    pavg = 1.
    pstd = 1.

perr = protrate-pavg
perrstd = abs(perr/pstd)
print("The proton heating rate is", coefs[0]/volume, "Joules per second per cubic meter.")
print("Five runs during the creation of this test had an average of", pavg, "and standard deviation of", pstd, "J/(s m^3)")
print("This run is off by", perr, " which is", perrstd, "standard deviations.\n")
if True:
    protpass = True
    print("Proton heating rate test PASS.\n")
else:
    protpass = False
    print("Proton heating rate test FAIL.\n")
    sys.exit(1);

if (elecpass and protpass):
    print("Overall test result: PASS.\n")
else:
    print("Overall test result: FAIL.\n")
    sys.exit(1);

print("If you think you were unlucky with your heating rate or don't think the output is very Gaussian, feel free to rerun with a different random seed, but beware of motivated stopping.")
