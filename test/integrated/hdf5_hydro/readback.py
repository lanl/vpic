#!/usr/bin/env python
import h5py
import numpy as np
import os.path
import sys

if len(sys.argv) != 2:
    sys.stderr.write("Usage: "+str(sys.argv[0])+" rundir\n")
    sys.exit(1)


# print(dd)
# print(dd.shape)
# exit(0)
rundir = sys.argv[1]
# print(rundir)


step_names = ["0",  "1"]

for step_name in step_names:
    filename = rundir+"/hydro_hdf5/T."+step_name+"/hydro_electron_" + step_name+".h5"

    if not os.path.isfile(filename):
        print("FAIL: "+filename+" is missing")
        sys.exit(1)

    infile = h5py.File(filename, 'r')
    datagroup = infile["Timestep_"+step_name+""]

    hydro_names = ["jx",  "jy", "jz", "rho"]

    for hydro_name in hydro_names:
        bin_filename = "step"+step_name+"_rank0_"+hydro_name+".bin"
        hydro_data_bi = np.fromfile(bin_filename, dtype=np.float32)
        # print(hydro_data_bi.shape)
        hydro_data_h5 = np.array(datagroup[hydro_name]).flatten()
        #print("hydro_data_bi=", hydro_data_bi[1:10])
        #print("hydro_data_h5=", hydro_data_h5[1:10])
        # print(hydro_data_h5.shape)
        if np.allclose(hydro_data_bi, hydro_data_h5, atol=0.0000000001) == False:
            print(hydro_name, 'in ', 'hydro_data_h5', " does not contain same value as ", bin_filename)
            print_max_element = 0
            print("     Binary Output", ",  ", "HDF5 Output",  ",  ", "Difference")
            for i in range(0, hydro_data_h5.shape[0]):
                if hydro_data_bi[i] - hydro_data_h5[i] > 0.0000000001:
                    print(i, hydro_data_bi[i], ",  ", hydro_data_h5[i], ",  ", hydro_data_bi[i] - hydro_data_h5[i])
                    print_max_element = print_max_element + 1
                if print_max_element == 10:
                    break
            sys.exit(1)

    infile.close()
