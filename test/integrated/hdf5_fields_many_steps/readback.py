#!/usr/bin/env python
import h5py
import numpy as np
import os.path
import sys

if len(sys.argv) != 2:
    sys.stderr.write("Usage: "+str(sys.argv[0])+" rundir\n")
    sys.exit(1)

rundir = sys.argv[1]
print(rundir)

for t in range(1, 5):
    file_name = rundir+"/field_hdf5/T."+str(t)+"/fields_"+str(t)+".h5"
    group_name = "Timestep_"+str(t)
    # print(file_name)
    # print(group_name)
    ##file_name = rundir+"/field_hdf5/T.2/fields_2.h5"

    if not os.path.isfile(file_name):
        print("FAIL: "+file_name+" is missing")
        sys.exit(1)

    infile = h5py.File(file_name, 'r')
    datagroup = infile[group_name]

    cbx = np.array(datagroup["cbx"])
    if np.array_equal(cbx, [3, 4, 3, 4]):
        print('cbx does not contain [3, 4, 3, 4].')
        sys.exit(1)

    cby = np.array(datagroup["cby"])
    if np.array_equal(cby, [3, 4, 3, 4]):
        print('cby does not contain all 17.')
        sys.exit(1)

    cbz = np.array(datagroup["cbz"])
    if np.array_equal(cbz, [3, 4, 3, 4]):
        print('cbz does not contain all 17.')
        sys.exit(1)

    ex = np.array(datagroup["ex"])
    if np.array_equal(ex, [3, 4, 3, 4]):
        print('ex does not contain all 17.')
        sys.exit(1)

    ey = np.array(datagroup["ey"])
    if np.array_equal(ey, [3, 4, 3, 4]):
        print('ey does not contain all 17.')
        sys.exit(1)

    ez = np.array(datagroup["ez"])
    if np.array_equal(ez, [3, 4, 3, 4]):
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
