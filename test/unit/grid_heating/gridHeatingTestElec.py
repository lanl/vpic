#!usr/bin/env python3

# This script tests the heating performance of VPIC.  It assumes you are
# using gridHeatingTestElec.cxx, and relies on some very specifically formated
# output therein.

# VPIC does not have the best heating rates because of low particle shapes, and
# this script checks if it matches what is expected, not if it is close enough to
# physical reality.

# Written by Scott V. Luedtke, XCP-6, August 14, 2019'''

import sys
import numpy as np
import scipy.constants as const
import os

base_dir = ""
if len(sys.argv) > 1:
    base_dir = sys.argv[1]

# Read in simulation specific parameters from a file
params_path = os.path.join(base_dir, 'heating_params.txt')
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

# Convert to SI, and don't blame me if some of these are wrong
wpe = 2. * np.pi * const.c * np.sqrt(ne_over_nc) / lamb  # Plasma frequency in Hz
timeToSI = 1. / wpe
lengthToSI = lamb / (2 * np.pi * np.sqrt(ne_over_nc))
efieldToSI = const.m_e * np.sqrt(ne_over_nc) * 2. * np.pi * const.c * const.c / (lamb * const.e)
energyToSI = efieldToSI ** 2 * lengthToSI ** 3 * const.epsilon_0

filename = os.path.join(base_dir, 'heating_rundata', 'energies')

d = np.loadtxt(filename, comments='%')
x = d[:, 0] * dt * timeToSI
# ex = d[:,1]*energyToSI
# ey = d[:,2]*energyToSI
# ez = d[:,3]*energyToSI
# bx = d[:,4]*energyToSI
# by = d[:,5]*energyToSI
# bz = d[:,6]*energyToSI
elec = d[:, 7] * energyToSI

# fields = ex + ey + ez + bx + by + bz

# To plot, uncomment this block
# plt.plot(x,fields/elec[0], 'r-', label='Fields')
# plt.plot(x,elec/elec[0], 'b-', label='Electrons')
# #plt.plot(x,elec+I2, 'k-', label='Everything')
# plt.legend(loc='upper left')
# #title = "$\Delta x = 0.5c/\omega_p$ 64ppc"
# title = "Grid Heating"
# ylabel = "Normalized Energy"
# xlabel = "Time (s)"
# plt.title(title, fontsize=25)
# plt.ylabel(ylabel, fontsize=20)
# plt.xlabel(xlabel, fontsize=20)
# plt.yscale('log', nonposy='clip')
# plt.show()

volume = box_size_x * box_size_y * box_size_z * 1e-6 ** 3  # SI
start = 100
stdpass = 5

coefs = np.polyfit(x[start:], elec[start:], 1)
# print("Electron coefs are ", coefs)
elecrate = coefs[0] / volume

# Computed from 5 runs on my desktop with different random seeds
avg = 1.24560647508e+28
std = 6.98124068576e+26

err = elecrate - avg
errstd = abs(err / std)
print("The electron heating rate is", elecrate, "Joules per second per cubic meter.")
print("Five runs during the creation of this test had an average of", avg, "and standard deviation of",
      std, "J/(s m^3)")
print("This run is off by", err, "which is", errstd, "standard deviations.")
print("This test will pass if within", stdpass, "standard deviations.\n")
if errstd < stdpass:
    elecpass = True
    print("Electron heating rate test PASS.\n")
else:
    elecpass = False
    print("Electron heating rate test FAIL.\n")

print("If you think you were unlucky with your heating rate or don't think the output is very Gaussian, "
      "feel free to rerun with a different random seed, but beware of motivated stopping.")
