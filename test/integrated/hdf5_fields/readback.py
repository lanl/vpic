#!/usr/bin/env python
import h5py
import numpy as np
import os.path
import sys
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
if len(sys.argv) != 2:
    sys.stderr.write("Usage: "+str(sys.argv[0])+" rundir\n")
    sys.exit(1)

rundir = sys.argv[1]
print(rundir)
filename = rundir+"/field_hdf5/T.2/fields_2.h5"

if not os.path.isfile(filename):
    print("FAIL: "+filename+" is missing")
    sys.exit(1)

infile = h5py.File(filename, 'r')
datagroup = infile["Timestep_2"]


cbx = np.array(datagroup["cbx"])
if np.all(cbx == 17) == False:
    print('cbx does not contain all 17.')
    sys.exit(1)

cby = np.array(datagroup["cby"])
if np.all(cby == 18) == False:
    print('cby does not contain all 17.')
    sys.exit(1)

cbz = np.array(datagroup["cbz"])
if np.all(cbz == 19) == False:
    print('cbz does not contain all 17.')
    sys.exit(1)

ex = np.array(datagroup["ex"])
if np.all(ex == 5) == False:
    print('ex does not contain all 17.')
    sys.exit(1)

ey = np.array(datagroup["ey"])
if np.all(ey == 6) == False:
    print('ey does not contain all 17.')
    sys.exit(1)

ez = np.array(datagroup["ez"])
if np.all(ez == 7) == False:
    print('ez does not contain all 17.')
    sys.exit(1)

other_fields_name = ["cmat", "div_b_err", "div_e_err", "ematx", "ematy", "ematz",
                     "fmatx", "fmaty", "fmatz", "jfx",  "jfy", "jfz", "nmat", "rhob", "rhof", "tcax", "tcay", "tcaz"]

for field_name in other_fields_name:
    field_data = np.array(datagroup[field_name])
    if np.all(field_data == 0) == False:
        print(field_name, 'does not contain all 0.')
        sys.exit(1)

infile.close()
