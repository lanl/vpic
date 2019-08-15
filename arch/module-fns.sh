#!/usr/bin/env bash
#
# Copyright (c) 2019 Carnegie Mellon University,
# Copyright (c) 2019 Triad National Security, LLC, as operator of
#     Los Alamos National Laboratory.
#
# All rights reserved.
#
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file. See the AUTHORS file for names of contributors.
#

#
# module-fns.sh  functions to interact with the module(1) system
# 18-Jul-2019  chuck@ece.cmu.edu
#

#
# input:
#  - we assume the set of loaded modules is in ${LOADEDMODULES}
# output variables:
#  - VCOM: currently loaded compiler (GNU, INT, CCE, PGI)
#  - VMPI: currently loaded MPI (CMPI, IMPI, OMPI)
#  - VCPU: currently selected CPU (broadwell, haswell, knl)
#  - mod_X: if 'X' module is loaded, set to X's version (or 1 if versionless)
#

if [ x${LOADEDMODULES} = x ]; then
    echo "ERROR: No modules loaded.  A compiler, MPI, and cmake are required."
    exit 1
fi
if [ "`printenv | fgrep CRAY_PRGENV`" ]; then
    craype=1
else
    craype=0
fi

#
# probe environment for list of loaded modules and their versions
#
lmods=`echo ${LOADEDMODULES} | sed -e 's/:/ /g'`
lmodsidx=0
for lm in ${lmods}
do
    name=`echo $lm | sed -e 's@/.*@@'`
    vers=`echo $lm | sed -e 's@.*/@@'`
    if [ x${name} = x${vers} ]; then
        vers=1      # fake version number for modules without one
    fi
    # convert - and . to underline so we can use it as a shell variable
    vname=`echo ${name} | sed -e 's/-/_/g' -e 's/\./_/g'`
    eval mod_${vname}=${vers}
    eval modpos_${vname}=${lmodsidx}
    lmodsidx=`expr ${lmodsidx} + 1`
done

#
# look for obvious problems...
#
if [ x${mod_cmake} = x ]; then
    echo "ERROR: A cmake module must be loaded to compile VPIC."
    echo "Try 'module load cmake'"
    exit 1
fi
if [ x${mod_craype_hugepages2M} != x ]; then
    echo "ERROR: craype-hugepages2M causes VPIC memory usage issues"
    echo "Try 'module unload craype-hugepages2M'"
    exit 1
fi

#
# determine which compiler is loaded and set VCOM
#
if [ x${mod_gcc} != x ]; then
    VCOM="GNU"
    VCOMVER=${mod_gcc}
    VCOMPOS=${modpos_gcc}
elif [ x${mod_intel} != x ]; then
    VCOM="INT"
    VCOMVER=${mod_intel}
    VCOMPOS=${modpos_intel}
elif [ x${mod_cce} != x ]; then
    VCOM="CCE"
    VCOMVER=${mod_cce}
    VCOMPOS=${modpos_cce}
elif [ x${mod_pgi} != x ]; then
    VCOM="PGI"
    VCOMVER=${mod_pgi}
    VCOMPOS=${modpos_pgi}
else
    echo "ERROR: no compiler module loaded."
    if [ ${craype} = 1 ]; then
        echo "Try 'module load PrgEnv-gnu' (or PrgEnv-intel, PrgEnv-cray, etc.)"
    else
        echo "Try 'module load gcc' (or intel, pgi, etc.)"
    fi
    exit 1
fi

#
# determine which MPI is loaded
#
if [ x${mod_cray_mpich} != x ]; then
    VMPI="CMPI"
    VMPIVER=${mod_cray_mpich}
elif [ x${mod_intel_mpi} != x ]; then
    VMPI="IMPI"
    VMPIVER=${mod_intel_mpi}
elif [ x${mod_openmpi} != x ]; then
    VMPI="OMPI"
    VMPIVER=${mod_openmpi}
    if [ x$craype = x0 -a ${modpos_openmpi} -lt ${VCOMPOS} ]; then
        #
        # XXX: there is one complete compiled version of openmpi for
        # every version of every compiled installed on the system.
        # when you "module load openmpi" it looks to see what compiler
        # has been loaded and uses that to select which version of
        # openmpi to add to your path.   if you load openmpi before
        # loading the compiler, the module system cannot determine
        # which compiled version of openmpi to put into your path.
        # the openmpi module should raise an error in this case, but
        # instead it picks up the non-module gcc from /usr/bin and
        # tries to use that, and that doesn't work as intented...
        # so we catch it here and die.
        #
        echo "ERROR: openmpi module added before a compiler module."
        echo "Try 'module purge' and load a compiler first, then openmpi."
        exit 1
    fi
else
    echo "ERROR: no MPI module loaded."
    echo "Try 'module load openmpi' (or intel-mpi, cray-mpich, etc.)"
    exit 1
fi

#
# heuristic to choose a CPU target
#
if [ x${mod_craype_haswell} != x ]; then
    VCPU="haswell"
elif [ x${mod_craype_mic_knl} != x ]; then
    VCPU="knl"
else
    VCPU="broadwell"      # XXX best guess
fi


echo "Compiler loaded: ${VCOM} version ${VCOMVER}"
echo "MPI loaded: ${VMPI} version ${VMPIVER}"
echo "CPU target: ${VCPU}"
