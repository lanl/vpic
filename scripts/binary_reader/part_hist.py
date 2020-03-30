#!/usr/bin/env python

#### Sample Usage
# python part_hist.py --ybins 100 --xvar voxel --plot $HOME/vpic/build/eparticle.168.0

import argparse
import glob
import numpy
import os.path
import struct
import sys
import sys
import textwrap
from collections import Counter
from vpic_binary_reader import *

parser = argparse.ArgumentParser(description="Make 1d or 2d histogramms of particle quantites", formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("filenames", nargs='+', metavar='filename')
parser.add_argument("--xmin", type=float, help="lower limit of the plot on the x axis")
parser.add_argument("--xmax", type=float, help="upper limit of the plot on the x axis")
parser.add_argument("--xvar", choices=["x","y","z", "ux","uy","uz","u", "gamma","gamma-1", "mu", "w", "voxel", "ix","iy","iz", "frac_x","frac_y","frac_z"], help=textwrap.dedent('''partice quantity to plot on the x axis
 x y z:                position of the particle
 ux uy uz:             components of u = \\beta * \\gamma
 u:                    magnitude of \\vec{u}
 gamma:                Lorentz factor \\gamma = sqrt(1. + u^2)
 gamma-1:              excess Lorentz factor \\gamma - 1
 mu:                   cosine of the pitch angle
 voxel:                cell index in VPIC
 ix iy iz:             integer cell coordinate
 frac_x frac_y frac_z: fractional position in cell [-1:1]'''))
parser.add_argument("--xbins", type=int, help="number of bins on the x axis. if in doubt use 100")
parser.add_argument("--ymin", type=float, help="lower limit of the plot on the y axis")
parser.add_argument("--ymax", type=float, help="upper limit of the plot on the y axis")
parser.add_argument("--yvar", help="particle quantity to plot on the y axis", choices=["x","y","z", "ux","uy","uz","u", "gamma", "gamma-1", "mu", "w", "voxel", "ix","iy","iz", "frac_x","frac_y","frac_z"])
parser.add_argument("--ybins", type=int, help="number of bins on the y axis. if in doubt use 100")
parser.add_argument("--species", action='append', type=int, help="accepted species IDs")
parser.add_argument("--pos-xmin", type=float, help="only accept particles right of POS_XMIN")
parser.add_argument("--pos-xmax", type=float, help="only accept particles left of POS_XMAX")
parser.add_argument("--pos-ymin", type=float, help="only accept particles right of POS_YMIN")
parser.add_argument("--pos-ymax", type=float, help="only accept particles left of POS_YMAX")
parser.add_argument("--pos-zmin", type=float, help="only accept particles right of POS_ZMIN")
parser.add_argument("--pos-zmax", type=float, help="only accept particles left of POS_ZMAX")
parser.add_argument("--first", type=int, help="only use the FIRST number of particles and skip the rest")
parser.add_argument("--random", type=int, help="read als particle but use only RANDOM of them for the hisotgram")
parser.add_argument("--out", metavar='filename', type=str, help="write histogram to the supplied filename")
parser.add_argument("--plot", action='store_true', help="plot the histogram")
parser.add_argument("--fields", metavar="DIR", help="load fields from this directory")
args = parser.parse_args()

if args.xvar is None:
    sys.stderr.write("You need to decide on a x variable\n")
    sys.exit(1)

if (args.ymin or args.ymax or args.ybins) and args.yvar is None:
    sys.stderr.write("You need to decide on a y variable\n")
    sys.exit(1)

if args.out is None and args.plot is False:
    sys.stderr.write("If you don't want to store nor plot the histogram, I'd rather be lazy\n")
    sys.stderr.write("Specify --out filename to save the histogram to a text file, or\n")
    sys.stderr.write("specify --plot to view the histogram interactivly\n")
    sys.exit(1)

collection_header = match_files(args.filenames)

if args.fields is not None:
    if args.xvar == "mu" or args.yvar == "mu":
        fieldfiles = glob.glob(args.fields+"/fields.*")
        bx = read_global_field(fieldfiles, 3)
        by = read_global_field(fieldfiles, 4)
        bz = read_global_field(fieldfiles, 5)
    else:
        sys.stderr.write("No need to load field files for these quantities\n")
        args.fields = None

X = []
Y = []
no_files = 0
no_part = 0
w_part = 0.
species_seen = []

for filename in args.filenames:
    f = open(filename, "rb")
    header = read_header(f)

    if header["dump_type"] != 3: #particles
        sys.stderr.write("Skipping "+filename+" (not a particle file)\n")
        continue

    if args.species is None:
        species_seen.append(header["sp_id"])
    else:
        if not (header["sp_id"] in args.species):
            sys.stderr.write("Skipping "+filename+" (wrong species)\n")
            continue

    array_header = read_array_header(f, header)

    no_files += 1

    for i in range(numpy.prod(array_header["dim"])):
        x,y,z, ux,uy,uz, w, voxel, ix,iy,iz, frac_x,frac_y,frac_z = get_next_particle(f, header)
        x += header["grid->x0"]
        y += header["grid->y0"]
        z += header["grid->z0"]
        ix += header["grid->x0"]
        iy += header["grid->y0"]
        iz += header["grid->z0"]
        u = numpy.sqrt(ux**2 + uy**2 + uz**2)
        gamma = numpy.sqrt(1. + u**2)
        if args.fields is not None:
            bxlocal = bx[ix-1,iy-1,iz-1]
            bylocal = by[ix-1,iy-1,iz-1]
            bzlocal = bz[ix-1,iy-1,iz-1]
            b = numpy.sqrt(bxlocal**2 + bylocal**2 + bzlocal**2)
            mu = (ux*bxlocal + uy*bylocal + uz*bzlocal)/(u * b)

        if args.pos_xmin is not None and x < args.pos_xmin:
            continue
        if args.pos_xmax is not None and x > args.pos_xmax:
            continue
        if args.pos_ymin is not None and y < args.pos_ymin:
            continue
        if args.pos_ymax is not None and y > args.pos_ymax:
            continue
        if args.pos_zmin is not None and z < args.pos_zmin:
            continue
        if args.pos_zmax is not None and z > args.pos_zmax:
            continue

        if args.xvar == "x":
            xval = float(x)
        elif args.xvar == "y":
            xval = float(y)
        elif args.xvar == "z":
            xval = float(z)
        elif args.xvar == "ux":
            xval = float(ux)
        elif args.xvar == "uy":
            xval = float(uy)
        elif args.xvar == "uz":
            xval = float(uz)
        elif args.xvar == "u":
            xval = float(u)
        elif args.xvar == "gamma":
            xval = float(gamma)
        elif args.xvar == "gamma-1":
            xval = float(gamma) - 1.
        elif args.xvar == "mu":
            xval = mu
        elif args.xvar == "w":
            xval = float(w)
        elif args.xvar == "voxel":
            xval = float(voxel)
        elif args.xvar == "ix":
            xval = float(ix)
        elif args.xvar == "iy":
            xval = float(iy)
        elif args.xvar == "iz":
            xval = float(iz)
        elif args.xvar == "frac_x":
            xval = float(frac_x)
        elif args.xvar == "frac_y":
            xval = float(frac_y)
        elif args.xvar == "frac_z":
            xval = float(frac_z)
        else:
            sys.stderr.write("Unknown args.xvar "+str(args.xvar)+"\n")
            sys.exit(1)
        if (args.xmin is not None and xval < args.xmin) or (args.xmax is not None and xval > args.xmax):
            continue

        if args.yvar is not None:
            if args.yvar == "x":
                yval = float(x)
            elif args.yvar == "y":
                yval = float(y)
            elif args.yvar == "z":
                yval = float(z)
            elif args.yvar == "ux":
                yval = float(ux)
            elif args.yvar == "uy":
                yval = float(uy)
            elif args.yvar == "uz":
                yval = float(uz)
            elif args.yvar == "u":
                yval = float(u)
            elif args.yvar == "gamma":
                yval = float(gamma)
            elif args.yvar == "gamma-1":
                yval = float(gamma) - 1.
            elif args.yvar == "mu":
                yval = mu
            elif args.yvar == "w":
                yval = float(w)
            elif args.yvar == "voxel":
                yval = float(voxel)
            elif args.yvar == "ix":
                yval = float(ix)
            elif args.yvar == "iy":
                yval = float(iy)
            elif args.yvar == "iz":
                yval = float(iz)
            elif args.yvar == "frac_x":
                yval = float(frac_x)
            elif args.yvar == "frac_y":
                yval = float(frac_y)
            elif args.yvar == "frac_z":
                yval = float(frac_z)
            else:
                sys.stderr.write("Unknown args.yvar "+str(args.yvar)+"\n")
                sys.exit(1)
            if (args.ymin is not None and yval < args.ymin) or (args.ymax is not None and yval > args.ymax):
                continue

        no_part += 1
        w_part += w

        if args.random is None or no_part <= args.random:
            X.append(xval)
            if args.yvar is not None:
                Y.append(yval)
        else:
            i = numpy.random.randint(0,args.random)
            X[i] = xval
            if args.yvar is not None:
                Y[i] = yval

        if args.first is not None and no_part >= args.first:
            break

    f.close()

if no_files == 0:
    sys.stderr.write("No files loaded. Abort.\n")
    sys.exit(1)

sys.stderr.write(str(no_part)+" particles loaded, eqivalent to "+str(w_part)+" physical particles\n")
sys.stderr.write(str(len(X))+" particles in histogram\n")

if args.xmin is None:
    args.xmin = min(X)
if args.xmax is None:
    args.xmax = max(X)
if args.yvar is not None:
    if args.ymin is None:
        args.ymin = min(Y)
    if args.ymax is None:
        args.ymax = max(Y)

if args.yvar is None:
    if args.xbins is None:
        args.xbins = 'auto'
    h, xedges = numpy.histogram(X, bins=args.xbins, range=(args.xmin,args.xmax))
else:
    if args.xbins is None:
        args.xbins = numpy.sqrt(len(X))
    if args.ybins is None:
        args.ybins = numpy.sqrt(len(X))

    h, xedges, yedges = numpy.histogram2d(X, Y, bins=[args.xbins, args.ybins], range=[[args.xmin,args.xmax],[args.ymin,args.ymax]])

# how to reproduce the plot
cmdstr = sys.argv[0]+" "

if args.out is not None:
    cmdstr += "--out "+str(args.out)+" "

if args.plot:
    cmdstr += "--plot "

if args.xvar is not None:
    cmdstr += "--xvar "+str(args.xvar)+" "
    cmdstr += "--xmin "+str(args.xmin)+" "
    cmdstr += "--xmax "+str(args.xmax)+" "
    cmdstr += "--xbins "+str(len(xedges)-1)+" "

if args.yvar is not None:
    cmdstr += "--yvar "+str(args.yvar)+" "
    cmdstr += "--ymin "+str(args.ymin)+" "
    cmdstr += "--ymax "+str(args.ymax)+" "
    cmdstr += "--ybins "+str(len(yedges)-1)+" "

if args.pos_xmin is not None:
    cmdstr += "--pos-xmin +"+str(args.pos_xmin)+" "
if args.pos_xmax is not None:
    cmdstr += "--pos-xmax +"+str(args.pos_xmax)+" "
if args.pos_ymin is not None:
    cmdstr += "--pos-ymin +"+str(args.pos_ymin)+" "
if args.pos_ymax is not None:
    cmdstr += "--pos-ymax +"+str(args.pos_ymax)+" "
if args.pos_zmin is not None:
    cmdstr += "--pos-zmin +"+str(args.pos_zmin)+" "
if args.pos_zmax is not None:
    cmdstr += "--pos-zmax +"+str(args.pos_zmax)+" "

if args.first is not None:
    cmdstr += "--first "+str(args.first)+" "

if args.species is not None:
    for s in args.species:
        cmdstr += "--species "+str(s)+" "
else:
    cnt = Counter(species_seen)
    for s in cnt.keys():
        cmdstr += "--species "+str(s)+" "

if args.fields is not None:
    cmdstr += "--fields "+str(args.fields)+" "

cmdstr += " ".join(args.filenames)

# save to file
if args.out is not None:
    with open(args.out, 'w') as outfile:
        if args.yvar is None:
            for i in range(len(xedges)-1):
                outfile.write(str(xedges[i])+" "+str(xedges[i+1])+" "+str(h[i])+"\n")
        else:
            for i in range(len(xedges)-1):
                for j in range(len(yedges)-1):
                    outfile.write(str(xedges[i])+" "+str(xedges[i+1])+" "+str(yedges[j])+" "+str(yedges[j+1])+" "+str(h[i][j])+"\n")
        outfile.write("# "+cmdstr+"\n")

# interactive plot
if args.plot:
    print cmdstr
    import matplotlib.pyplot as plt
    if args.yvar is None:
        left,right = xedges[:-1],xedges[1:]
        X = numpy.array([left,right]).T.flatten()
        Y = numpy.array([h,h]).T.flatten()
        plt.plot(X,Y)
        plt.xlabel(args.xvar)
        plt.ylabel("f("+str(args.xvar)+") (a.u.)")
        plt.show()
    else:
        Xgrid, Ygrid = numpy.meshgrid(xedges, yedges)
        plt.pcolormesh(Xgrid, Ygrid, h.T)
        plt.xlabel(args.xvar)
        plt.ylabel(args.yvar)
        plt.title("f("+str(args.xvar)+","+str(args.yvar)+") (a.u.)")
        plt.xlim(xedges[0],xedges[-1])
        plt.ylim(yedges[0],yedges[-1])
        plt.show()

