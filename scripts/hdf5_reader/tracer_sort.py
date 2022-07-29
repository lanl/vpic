#!/usr/bin/env python
import copy
import h5py
import numpy as np
import os
import sys
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
if len(sys.argv) != 2:
	sys.stderr.write("Usage: "+str(sys.argv[0])+" tracers.h5p\n")
	sys.exit(1)

infilename=sys.argv[1]
dirname=os.path.dirname(infilename)
outfilename=dirname+"/tracers_sorted.h5p"

infile=h5py.File(infilename, 'r')
outfile=h5py.File(outfilename, 'w')


for groupname in infile:
	if not groupname.startswith("Step#"):
		continue
	group = infile[groupname]
	outgroup = outfile.create_group(groupname)

	# load meta data
	meta=group["grid_meta_data"]
	nx,ny,nz=np.array(meta["cells"])
	x0,y0,z0=np.array(meta["offset"])
	dx,dy,dz=np.array(meta["resolution"])
	print("cells = "+str(nx)+" "+str(ny)+" "+str(nz))
	print("offset = "+str(x0)+" "+str(y0)+" "+str(z0))
	print("resolution = "+str(dx)+" "+str(dy)+" "+str(dz))

	# copy meta data to output file
	metaoutgroup=outgroup.create_group("grid_meta_data")
	for datasetname in meta:
		data = meta[datasetname]
		metaoutgroup.create_dataset(datasetname, data=data)

	# iterate through data
	for speciesname in group:
		if speciesname.startswith("grid_meta_data"):
			continue
		subgroup=group[speciesname]
		suboutgroup=outgroup.create_group(speciesname)
		ts=np.array(subgroup["timestep"])
		cond = ts!=-1
		id=np.array(subgroup["q"])[cond]
		shuffle=np.argsort(id, kind="mergesort") #stable sort algorithm is important to make sure timesteps stay in order

		# write sorted version of the data to output file
		for datasetname in subgroup:
			print(datasetname)
			data=subgroup[datasetname][cond]
			shuffleddata=np.empty_like(data)
			for i in range(len(shuffle)):
				shuffleddata[i] = data[shuffle[i]]
			suboutgroup.create_dataset(datasetname, data=shuffleddata)

		# add computed x/y/z position
		ii = np.array(subgroup["i"])[cond]

		if(np.any(ii<0) or np.any(ii>=nx*ny*nz)):
			sys.stderr.write("Invalid voxel index, not computing x/y/z\n")
			continue

		dX = np.array(subgroup["dX"])[cond]
		dY = np.array(subgroup["dY"])[cond]
		dZ = np.array(subgroup["dZ"])[cond]

		x = np.empty_like(dX)
		y = np.empty_like(dY)
		z = np.empty_like(dZ)

		for i in range(len(ii)):
			voxel = ii[shuffle[i]]
			xypart = voxel % (nx*ny)
			iz = (voxel - xypart) // (nx*ny)
			ix = xypart % nx
			iy = (xypart - ix) // nx

			x[i] = x0 + ix*dx + 0.5*(1.+dX[shuffle[i]])*dx
			y[i] = y0 + iy*dy + 0.5*(1.+dY[shuffle[i]])*dy
			z[i] = z0 + iz*dz + 0.5*(1.+dZ[shuffle[i]])*dz

		suboutgroup.create_dataset("x", data=x)
		suboutgroup.create_dataset("y", data=y)
		suboutgroup.create_dataset("z", data=z)

infile.close()
outfile.close()
