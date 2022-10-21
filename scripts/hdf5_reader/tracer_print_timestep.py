#!/usr/bin/env python
import copy
import h5py
import numpy as np
import os
import sys
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
if len(sys.argv) < 3:
	sys.stderr.write("Usage: "+str(sys.argv[0])+" timestep tracers.h5p [other.h5p...]\n")
	sys.exit(1)

timestep=int(sys.argv[1])

#iterate through input files
for infilename in sys.argv[2:]:
	infile=h5py.File(infilename, 'r')

	# iterate through timestep groups
	for groupname in infile:
		if not groupname.startswith("Step#"):
			continue
		group = infile[groupname]

		# iterate through species
		for speciesname in group:
			if speciesname.startswith("grid_meta_data"):
				continue
			subgroup=group[speciesname]

			# dataset names
			dsetnames = []
			for datasetname in subgroup:
				if datasetname != "timestep":
					dsetnames.append(datasetname)
			print("#timestep "+" ".join(map(str, dsetnames)))

			ts=np.array(subgroup["timestep"])

			# iterate over particles
			for n in range(len(ts)):
				data = []
				if ts[n] == timestep:
					data.append(timestep)
					for dset in dsetnames:
						data.append(subgroup[dset][n])
					print(" ".join(map(str, data)))

	infile.close()
