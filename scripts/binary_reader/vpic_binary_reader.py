import struct
import sys
import numpy
from collections import Counter

def read_header(inputfile):
    f = inputfile
    header = {}

    f.seek(0,2) # go to end of file
    end_byte = f.tell()
    f.seek(0,0) # go to beginning of file
    start_byte = f.tell()

    header["end_byte"] = end_byte
    header["file_bytes"] = end_byte - start_byte

    char_bits = struct.unpack('=1b', f.read(1))[0]
    if char_bits != 8:
        sys.stderr.write("Header has "+str(char_bits)+" per byte, not eight. That is cool but scary\n")
        sys.exit(1)
    else:
        header["CHAR_BITS"] = char_bits

    shortint_len = struct.unpack('=1b', f.read(1))[0]
    if shortint_len != 2:
        sys.stderr.write("sizeof(short int)="+str(shortint_len)+", not two. That is possible but scary\n")
        sys.exit(1)
    else:
        header["sizeof(short int)"] = shortint_len

    int_len = struct.unpack('=1b', f.read(1))[0]
    if int_len != 4:
        sys.stderr.write("sizeof(int)="+str(int_len)+", not four. That is possible but scary\n")
        sys.exit(1)
    else:
        header["sizeof(int)"] = int_len

    float_len = struct.unpack('=1b', f.read(1))[0]
    if float_len == 4:
        header["sizeof(float)"] = float_len
        floatfmt = "f"
        header["floatfmt"] = floatfmt
    elif float_len == 8:
        sys.stderr.write("WARNING: this file was produced by the code running in double precission\n")
        header["sizeof(float)"] = float_len
        floatfmt = "d"
        header["floatfmt"] = floatfmt
    else:
        sys.stderr.write("sizeof(float)="+str(float_len)+", not four. That is possible but scary\n")
        sys.exit(1)

    double_len = struct.unpack('=1b', f.read(1))[0]
    if double_len != 8:
        sys.stderr.write("sizeof(double)="+str(double_len)+", not eight. That is possible but scary\n")
        sys.exit(1)
    else:
        header["sizeof(double)"] = double_len

    magic1 = struct.unpack('=1H', f.read(shortint_len))[0] # technically H is unsigned short int, but that is the decoding that int( ,16) uses
    if magic1 != int("0xcafe", 16):
        sys.stderr.write("magic1="+str(magic1)+", not 0xcafe. Something is wrong\n")
        sys.exit(1)
    else:
        header["magic1"] = magic1

    magic2 = struct.unpack('=1I', f.read(int_len))[0]
    if magic2 != int("0xdeadbeef", 16):
        sys.stderr.write("magic2="+str(magic2)+", not 0xdeadbeef. Something is wrong\n")
        sys.exit(1)
    else:
        header["magic2"] = magic2

    float_const = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    if float_const != 1.:
        sys.stderr.write("float_const="+str(float_const)+", not 1.0. Something is wrong\n")
        sys.exit(1)
    else:
        header["float_const"] = float_const

    double_const = struct.unpack('=1d', f.read(double_len))[0]
    if double_const != 1.:
        sys.stderr.write("double_const="+str(double_const)+", not 1.0. Something is wrong\n")
        sys.exit(1)
    else:
        header["double_const"] = double_const

    version = struct.unpack('=1i', f.read(int_len))[0]
    if version != 0:
        sys.stderr.write("version="+str(version)+", not 0. Something is wrong\n")
        sys.exit(1)
    else:
        header["version"] = version

    dump_type = struct.unpack('=1i', f.read(int_len))[0]
    if dump_type not in [0, 1, 2, 3, 4, 5]:
        sys.stderr.write("dump_type="+str(dump_type)+", not a known type. Something is wrong\n")
        sys.exit(1)
    else:
        header["dump_type"] = dump_type

    header["step()"] = struct.unpack('=1i', f.read(int_len))[0]
    header["nxout"] = struct.unpack('=1i', f.read(int_len))[0]
    header["nyout"] = struct.unpack('=1i', f.read(int_len))[0]
    header["nzout"] = struct.unpack('=1i', f.read(int_len))[0]
    header["grid->dt"] = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    header["dxout"] = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    header["dyout"] = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    header["dzout"] = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    header["grid->x0"] = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    header["grid->y0"] = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    header["grid->z0"] = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    header["grid->cvac"] = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    header["grid->eps0"] = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    header["damp"] = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    header["rank()"] = struct.unpack('=1i', f.read(int_len))[0]
    header["nproc()"] = struct.unpack('=1i', f.read(int_len))[0]
    header["sp_id"] = struct.unpack('=1i', f.read(int_len))[0]
    header["q_m"] = struct.unpack('=1'+floatfmt, f.read(float_len))[0]

    header["offset_x"] = int(round(header["grid->x0"] / header["dxout"]))
    header["offset_y"] = int(round(header["grid->y0"] / header["dyout"]))
    header["offset_z"] = int(round(header["grid->z0"] / header["dzout"]))

    return header

def read_array_header(inputfile, header):
    f = inputfile
    array_header = {}

    array_header["sizeof(p[0])"] = struct.unpack('=1i', f.read(header["sizeof(int)"]))[0]
    ndim = struct.unpack('=1i', f.read(header["sizeof(int)"]))[0]
    array_header["ndim"] = ndim
    dim = struct.unpack('='+str(ndim)+'i', f.read(ndim*int(header["sizeof(int)"])))[0:ndim]
    array_header["dim"] = dim
    elements = numpy.prod(dim)

    cur_byte = f.tell()
    unread_byte = header["end_byte"] - cur_byte
    no_fields = unread_byte / (elements * header["sizeof(float)"])

    array_header["no_fields"] = no_fields

    return array_header

def match_files(filelist):
    min_x = sys.maxint
    max_x =-sys.maxint-1
    min_y = sys.maxint
    max_y =-sys.maxint-1
    min_z = sys.maxint
    max_z =-sys.maxint-1
    no_fields = []
    no_files = 0

    for filename in filelist:
        f = open(filename, "rb")
        header = read_header(f)
        array_header = read_array_header(f, header)
        no_fields.append(array_header["no_fields"])

        min_x = min(min_x, header["offset_x"])
        min_y = min(min_y, header["offset_y"])
        min_z = min(min_z, header["offset_z"])

        max_x = max(max_x, header["offset_x"]+header["nxout"])
        max_y = max(max_y, header["offset_y"]+header["nyout"])
        max_z = max(max_z, header["offset_z"]+header["nzout"])

        f.close()
        no_files += 1

    size_x = max_x - min_x
    size_y = max_y - min_y
    size_z = max_z - min_z

    header = {}
    header["min_x"] = min_x
    header["min_y"] = min_y
    header["min_z"] = min_z
    header["max_x"] = max_x
    header["max_y"] = max_y
    header["max_z"] = max_z
    header["size_x"] = size_x
    header["size_y"] = size_y
    header["size_z"] = size_z

    if no_files == 0:
        sys.stderr.write("No files of dumptype 1 or 2 found\n")
        sys.exit(1)

    if no_fields.count(min(no_fields)) != len(no_fields):
        cnt = Counter(no_fields)
        sys.stderr.write("Not all files seemed to contain the same number of fields\n")
        sys.stderr.write("The most common numbers were "+str(cnt)+"\n")
        sys.exit(1)

    header["no_fields"] = min(no_fields)
    header["no_files"] = no_files

    return header

def read_data(inputfile, dim, floatfmt, no_fields):
    f = inputfile
    volume = no_fields*dim[0]*dim[1]*dim[2]

    #data1d = struct.unpack('='+str(volume)+'f', f.read(volume*4))
    data1d = numpy.fromfile(f,dtype=numpy.single,count=volume)

    #dim4 = [dim[0], dim[1], dim[2], no_fields]
    dim4 = [dim[2], dim[1], dim[0], no_fields]
    data4d = numpy.reshape(data1d, dim4, order='C')

    ret_list = []
    for i in range(no_fields):
        ret_list.append(data4d[:,:,:, i])

    return ret_list

def read_field(inputfile, dim, floatfmt, field_i, no_fields):
    f = inputfile
    volume = no_fields*dim[0]*dim[1]*dim[2]

    if field_i<0 or field_i >= no_fields:
        sys.stderr.write("Invalid field #"+str(field_i)+" out of "+str(no_fields)+" fields available\n")
        sys.exit(1)

    data1d = struct.unpack('='+str(volume)+'f', f.read(volume*4))

    #dim4 = [dim[0], dim[1], dim[2], no_fields]
    dim4 = [dim[2], dim[1], dim[0], no_fields]
    data4d = numpy.reshape(data1d, dim4, order='C')

    return data4d[:,:,:, field_i]

def get_next_particle(inputfile, header):
    f = inputfile
    int_len   = header["sizeof(int)"]
    float_len = header["sizeof(float)"]
    floatfmt  = header["floatfmt"]

    frac_x = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    frac_y = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    frac_z = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    voxel  = struct.unpack('=1I', f.read(int_len))[0]
    ux     = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    uy     = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    uz     = struct.unpack('=1'+floatfmt, f.read(float_len))[0]
    w      = struct.unpack('=1'+floatfmt, f.read(float_len))[0]

    # from vpic/src/grid/ops.c
    sx = 1
    sy = (header["nxout"]+2)*sx
    sz = (header["nyout"]+2)*sy
    nv = (header["nzout"]+2)*sz

    if voxel < 0 or voxel >= nv:
        sys.stderr.write("Invalid voxel index "+str(voxel)+"\n")
        sys.stderr.write("sx = "+str(sx)+", sy = "+str(sy)+", sz = "+str(sz)+"\n")
        sys.stderr.write("nv = "+str(nv)+"\n")

    xypart = voxel % sz
    iz = (voxel - xypart) // sz
    ix = xypart % sy
    iy = (xypart - ix) // sy

    x = ix + 0.5*(1.+frac_x)
    y = iy + 0.5*(1.+frac_y)
    z = iz + 0.5*(1.+frac_z)

    return x,y,z, ux,uy,uz, w, voxel, ix,iy,iz, frac_x,frac_y,frac_z

def read_global_field(filelist, i_field):
    collection_header = match_files(filelist)

    global_field = numpy.full( (collection_header["size_x"],collection_header["size_y"],collection_header["size_z"]), numpy.nan)
    covered = numpy.full( (collection_header["size_x"],collection_header["size_y"],collection_header["size_z"]), False, dtype="bool")

    for filename in filelist:
        f = open(filename, "rb")
        header = read_header(f)
        array_header = read_array_header(f, header)
        local_field = read_field(f, array_header["dim"], header["floatfmt"], i_field, collection_header["no_fields"])

        # the data in the file has shape
        dim = array_header["dim"]
        # and we want to start reading at the second element in each direction (removing 1 cell of ghost zone)
        file_start = [1, 1, 1]
        # and stop 1 element short of the upper end
        file_stop  = [array_header["dim"][0]-1, array_header["dim"][1]-1, array_header["dim"][2]-1]

        # in memory we want to start writing at the correct offset
        mem_start = [header["offset_x"]-collection_header["min_x"], header["offset_y"]-collection_header["min_y"], header["offset_z"]-collection_header["min_z"]]
        # and stop at the correct end
        mem_stop  = [mem_start[0]+header["nxout"], mem_start[1]+header["nyout"], mem_start[2]+header["nzout"]]

        global_field[mem_start[0]:mem_stop[0], mem_start[1]:mem_stop[1], mem_start[2]:mem_stop[2]] = local_field[file_start[0]:file_stop[0], file_start[1]:file_stop[1], file_start[2]:file_stop[2] ]
        covered[mem_start[0]:mem_stop[0], mem_start[1]:mem_stop[1], mem_start[2]:mem_stop[2]] = True

        f.close()

    if not numpy.all(covered):
        sys.stderr.write("Not all parts of the global data have been recovered\n")
        sys.exit(1)
    elif numpy.any(numpy.isnan(global_field)):
        sys.stderr.write("field data contains NaN\n")
        sys.exit(1)

    return global_field
