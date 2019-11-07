#!/usr/bin/env python
import struct
import sys
import numpy
from vpic_binary_reader import *

#### Sample Usage
# Read fields
# python read_file.py ./vpic/build/field/T.2/field.2.0
# Read hydro
# python read_file.py ./vpic/build/hydro/T.0/dihydro.0.0
# Read particles
# python read_file.py ./vpic/build/eparticle.168.0

if len(sys.argv) < 2:
    sys.stderr.write("Usage: "+sys.argv[0]+" filename [filename...]\n")
    sys.exit(1)

filename = sys.argv[1]

f = open(filename, "rb")
header = read_header(f)
print "HEADER: "
for k,v in sorted(header.items()):
    print str(k)+": "+str(v)
print ""

array_header = read_array_header(f, header)
print "ARRAY HEADER: "
for k,v in sorted(array_header.items()):
    print str(k)+": "+str(v)
print ""

cur_byte = f.tell()
print str(header["end_byte"] - cur_byte)+" bytes data"

if header["dump_type"] == 1: # fields
    print("Num Fields %s" % array_header["no_fields"])

    #WARNING: Users may have to edit this line by hand if they want to unpack different fields
    ex,ey,ez,bx,by,bz = read_data(f, array_header["dim"], header["floatfmt"], array_header["no_fields"])

    cur_byte = f.tell()
    print str(header["end_byte"] - cur_byte)+" bytes unread"

    # print ex
    minimum = numpy.nanmin(ex)
    maximum = numpy.nanmax(ex)
    print "Ex: "+str(minimum)+"..."+str(maximum)

    minimum = numpy.nanmin(ey)
    maximum = numpy.nanmax(ey)
    print "Ey: "+str(minimum)+"..."+str(maximum)

    minimum = numpy.nanmin(ez)
    maximum = numpy.nanmax(ez)
    print "Ez: "+str(minimum)+"..."+str(maximum)

    minimum = numpy.nanmin(bx)
    maximum = numpy.nanmax(bx)
    print "Bx: "+str(minimum)+"..."+str(maximum)

    minimum = numpy.nanmin(by)
    maximum = numpy.nanmax(by)
    print "By: "+str(minimum)+"..."+str(maximum)

    minimum = numpy.nanmin(bz)
    maximum = numpy.nanmax(bz)
    print "Bz: "+str(minimum)+"..."+str(maximum)

elif header["dump_type"] == 2: # hydro
    jx,jy,jz,rho,txx,tyy,tzz,txy,tyz,tzy = read_data(f, array_header["dim"], header["floatfmt"], array_header["no_fields"])
    cur_byte = f.tell()
    print str(header["end_byte"] - cur_byte)+" bytes unread"

    # print jx
    minimum = numpy.nanmin(jx)
    maximum = numpy.nanmax(jx)
    print "jx: "+str(minimum)+"..."+str(maximum)

elif header["dump_type"] == 3: #particles
    sys.stderr.write("Dealing with a particle file\n")
    x,y,z, ux,uy,uz, w, voxel, ix,iy,iz, frac_x,frac_y,frac_z = get_next_particle(f, header)
    print "cell: "+str(ix)+" "+str(iy)+" "+str(iz)
    print "frac: "+str(frac_x)+" "+str(frac_y)+" "+str(frac_z)
    print "pos:  "+str(x)+" "+str(y)+" "+str(z)
    print "u:    "+str(ux)+" "+str(uy)+" "+str(uz)
    print "w:    "+str(w)
    print "voxel:"+str(voxel)

    # exhaust the file
    #for i in range(numpy.prod(array_header["dim"])-1):
    #    x,y,z, ux,uy,uz, w, voxel, ix,iy,iz, frac_x,frac_y,frac_z = get_next_particle(f, header)

f.close()
print ""

collection_header = match_files(sys.argv[1:])
print "COLLECTION HEADER: "
for k,v in sorted(collection_header.items()):
    print str(k)+": "+str(v)
print ""
