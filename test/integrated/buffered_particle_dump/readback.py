#!/usr/bin/env python
import h5py
import numpy as np
import os.path
import sys

if len(sys.argv) != 2:
    sys.stderr.write("Usage: "+str(sys.argv[0])+" rundir\n")
    sys.exit(1)

rundir=sys.argv[1]
filename=rundir+"/tracers/T.2/tracers.h5p"

if not os.path.isfile(filename):
    print("FAIL: "+filename+" is missing")
    sys.exit(1)

infile=h5py.File(filename, 'r')

timestepgroup=infile["Step#2"]
metagroup=timestepgroup["grid_meta_data"]

# cells counts should match the global 2x2x1 domain size
cells=np.array(metagroup["cells"])
if not np.array_equal(cells, np.array([2,2,1])):
    print("FAIL: wrong global cells counts")
    sys.exit(1)

# both ranks should have two particles in each of the timesteps
np_local=np.array(metagroup["np_local_electron_tracer"])
if not np.array_equal(np_local, np.array([[2,2],[2,2],[2,2]])):
    print("FAIL: some rank didn't have the expected two particles")
    sys.exit(1)

# the global domain starts at 0,0,0
offset=np.array(metagroup["offset"])
if not np.array_equal(offset, np.array([0,0,0])):
    print("FAIL: domain doesn't start at coordinate origin")
    sys.exit(1)

# resolution should be 1/8 d_e in all directions
resolution=np.array(metagroup["resolution"])
if not np.array_equal(resolution, np.array([0.125,0.125,0.125])):
    print("FAIL: resolution unexpectedly deviates from 0.125 d_e")
    sys.exit(1)

datagroup=timestepgroup["electron_tracer"]
Ux=np.array(datagroup["Ux"])
Uy=np.array(datagroup["Uy"])
Uz=np.array(datagroup["Uz"])
dX=np.array(datagroup["dX"])
dY=np.array(datagroup["dY"])
dZ=np.array(datagroup["dZ"])
i=np.array(datagroup["i"])
w=np.array(datagroup["w"])
timestep=np.array(datagroup["timestep"])

if np.any(w!=0.):
    print("FAIL: tracer particle with non-zero weight")
    sys.exit(1)

for ts in range(3):
    cond = timestep==ts
    if np.sum(cond) != 4:
        print("FAIL: "+str(numpy.sum(cond))+" particles in step "+str(ts)+" instead of 4")
        sys.exit(1)

    if np.any (Ux[cond]!=0.5):
        printf("FAIL: Ux changed")
        sys.exit(1)

    if np.any (Uy[cond]!=0.0):
        printf("FAIL: Uy changed")
        sys.exit(1)

    if np.any (Uz[cond]!=0.0):
        printf("FAIL: Uz changed")
        sys.exit(1)

#timestep 0
# particle 1
cond_t   = timestep==0
cond_i   = i==0
cond_dx  = dX==0.
cond_dyl = 0.33<dY
cond_dyu =      dY<0.34
cond_dzl = 0.33<dZ
cond_dzu =      dZ<0.34
cond_p1 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),cond_dx), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p1) != 1:
    printf("FAIL: can't find the particle in cell 0, offset 0.,0.33,0.33 in timestep 0")
    sys.exit(1)

#particle 2
cond_t   = timestep==0
cond_i   = i==0
cond_dx  = dX==0.
cond_dyl = 0.66<dY
cond_dyu =      dY<0.67
cond_dzl = 0.33<dZ
cond_dzu =      dZ<0.34
cond_p2 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),cond_dx), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p2) != 1:
    printf("FAIL: can't find the particle in cell 0, offset 0.,0.66,0.33 in timestep 0")
    sys.exit(1)

#particle 3
cond_t   = timestep==0
cond_i   = i==1
cond_dx  = dX==0.
cond_dyl = 0.33<dY
cond_dyu =      dY<0.34
cond_dzl = 0.66<dZ
cond_dzu =      dZ<0.67
cond_p3 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),cond_dx), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p3) != 1:
    printf("FAIL: can't find the particle in cell 1, offset 0.,0.33,0.66 in timestep 0")
    sys.exit(1)

#particle 3
cond_t   = timestep==0
cond_i   = i==1
cond_dx  = dX==0.
cond_dyl = 0.66<dY
cond_dyu =      dY<0.67
cond_dzl = 0.66<dZ
cond_dzu =      dZ<0.67
cond_p4 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),cond_dx), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p4) != 1:
    printf("FAIL: can't find the particle in cell 1, offset 0.,0.66,0.66 in timestep 0")
    sys.exit(1)

#timestep 1
# particle 1
cond_t   = timestep==1
cond_i   = i==0
cond_dxl = 0.49<dX
cond_dxu =      dX<0.50
cond_dyl = 0.33<dY
cond_dyu =      dY<0.34
cond_dzl = 0.33<dZ
cond_dzu =      dZ<0.34
cond_p1 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),np.logical_and(cond_dxl,cond_dxu)), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p1) != 1:
    print("FAIL: can't find the particle in cell 0, offset 0.49,0.33,0.33 in timestep 1")
    sys.exit(1)

# particle 2
cond_t   = timestep==1
cond_i   = i==0
cond_dxl = 0.49<dX
cond_dxu =      dX<0.50
cond_dyl = 0.66<dY
cond_dyu =      dY<0.67
cond_dzl = 0.33<dZ
cond_dzu =      dZ<0.34
cond_p2 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),np.logical_and(cond_dxl,cond_dxu)), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p2) != 1:
    print("FAIL: can't find the particle in cell 0, offset 0.49,0.66,0.33 in timestep 1")
    sys.exit(1)

# particle 3
cond_t   = timestep==1
cond_i   = i==1
cond_dxl = 0.49<dX
cond_dxu =      dX<0.50
cond_dyl = 0.33<dY
cond_dyu =      dY<0.34
cond_dzl = 0.66<dZ
cond_dzu =      dZ<0.67
cond_p3 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),np.logical_and(cond_dxl,cond_dxu)), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p3) != 1:
    print("FAIL: can't find the particle in cell 1, offset 0.49,0.33,0.66 in timestep 1")
    sys.exit(1)

# particle 4
cond_t   = timestep==1
cond_i   = i==1
cond_dxl = 0.49<dX
cond_dxu =      dX<0.50
cond_dyl = 0.66<dY
cond_dyu =      dY<0.67
cond_dzl = 0.66<dZ
cond_dzu =      dZ<0.67
cond_p4 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),np.logical_and(cond_dxl,cond_dxu)), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p4) != 1:
    print("FAIL: can't find the particle in cell 1, offset 0.49,0.66,0.66 in timestep 1")
    sys.exit(1)



#timestep 2
# particle 1
cond_t   = timestep==2
cond_i   = i==0
cond_dxl = 0.98<dX
cond_dxu =      dX<0.99
cond_dyl = 0.33<dY
cond_dyu =      dY<0.34
cond_dzl = 0.33<dZ
cond_dzu =      dZ<0.34
cond_p1 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),np.logical_and(cond_dxl,cond_dxu)), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p1) != 1:
    print("FAIL: can't find the particle in cell 0, offset 0.98,0.33,0.33 in timestep 2")
    sys.exit(1)

# particle 2
cond_t   = timestep==2
cond_i   = i==0
cond_dxl = 0.98<dX
cond_dxu =      dX<0.99
cond_dyl = 0.66<dY
cond_dyu =      dY<0.67
cond_dzl = 0.33<dZ
cond_dzu =      dZ<0.34
cond_p2 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),np.logical_and(cond_dxl,cond_dxu)), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p2) != 1:
    print("FAIL: can't find the particle in cell 0, offset 0.98,0.66,0.33 in timestep 2")
    sys.exit(1)

# particle 3
cond_t   = timestep==2
cond_i   = i==1
cond_dxl = 0.98<dX
cond_dxu =      dX<0.99
cond_dyl = 0.33<dY
cond_dyu =      dY<0.34
cond_dzl = 0.66<dZ
cond_dzu =      dZ<0.67
cond_p3 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),np.logical_and(cond_dxl,cond_dxu)), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p3) != 1:
    print("FAIL: can't find the particle in cell 1, offset 0.98,0.33,0.66 in timestep 2")
    sys.exit(1)

# particle 4
cond_t   = timestep==2
cond_i   = i==1
cond_dxl = 0.98<dX
cond_dxu =      dX<0.99
cond_dyl = 0.66<dY
cond_dyu =      dY<0.67
cond_dzl = 0.66<dZ
cond_dzu =      dZ<0.67
cond_p4 = np.logical_and(np.logical_and(np.logical_and(cond_t,cond_i),np.logical_and(cond_dxl,cond_dxu)), np.logical_and(np.logical_and(cond_dyl,cond_dyu),np.logical_and(cond_dzl,cond_dzu)))
if np.sum(cond_p4) != 1:
    print("FAIL: can't find the particle in cell 1, offset 0.98,0.66,0.66 in timestep 2")
    sys.exit(1)

infile.close()
