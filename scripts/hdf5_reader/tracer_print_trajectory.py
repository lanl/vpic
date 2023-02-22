#!/usr/bin/env python
import copy
import h5py
import numpy as np
import os
import sys
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
if len(sys.argv) < 3:
	sys.stderr.write("Usage: "+str(sys.argv[0])+" id tracers.h5p [other.h5p...]\n")
	sys.exit(1)

id=int(sys.argv[1])

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
				if datasetname != "q":
					dsetnames.append(datasetname)
			print("#id "+" ".join(map(str, dsetnames)))

			ids=np.array(subgroup["q"])

			# iterate over particles
			for n in range(len(ids)):
				data = []
				if ids[n] == id:
					data.append(id)
					for dset in dsetnames:
						data.append(subgroup[dset][n])
					print(" ".join(map(str, data)))

	infile.close()
