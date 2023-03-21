#!/bin/bash

#SBATCH -N 1
#SBATCH -t 4:00:00 
#SBATCH --output=vpic-build-log.txt

src_dir=/users/matsekh/VPIC/vpic-master
builds=/lustre/scratch5/.mdt0/matsekh/VPIC/test-builds
mkdir -p $builds
openmp=$builds/openmp
pthreads=$builds/pthreads

### Build weak and strong scaling lpi_2d_F6_test input decks ###
function lpi_deck(){
    $1/bin/vpic $2/lpi-deck/lpi_2d_F6_test        # strong scaling deck
    $1/bin/vpic $2/lpi-deck/lpi_2d_F6_test-nx12   # weak scaling decks 
    $1/bin/vpic $2/lpi-deck/lpi_2d_F6_test-nx24
    $1/bin/vpic $2/lpi-deck/lpi_2d_F6_test-nx48
    $1/bin/vpic $2/lpi-deck/lpi_2d_F6_test-nx96
    $1/bin/vpic $2/lpi-deck/lpi_2d_F6_test-nx288
    $1/bin/vpic $2/lpi-deck/lpi_2d_F6_test-nx864
}

### Configure VPIC cmake ###
function vpic_cmake(){
cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_INTEGRATED_TESTS=OFF \
  -DUSE_V4_AVX=ON \
  -DUSE_V4_AVX2=ON \
  -DUSE_V8_AVX=ON \
  -DUSE_V8_AVX2=ON \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_C_FLAGS="-O3 -dynamic -craype-verbose" \
  -DCMAKE_CXX_FLAGS="-O3 -dynamic -craype-verbose" \
  -DCMAKE_EXE_LINKER_FLAGS="-dynamic" \
  -DUSE_PTHREADS=$2 \
  -DUSE_OPENMP=$3 \
  -DUSE_LEGACY_SORT=$4 \
  $1
}

# PrgEnv-cray, OpenMP, TSORT
build_dir=$openmp/cray/tsort 
mkdir -p $build_dir 
cd $build_dir 
cmake_pthreads="OFF"
cmake_openmp="ON" 
cmake_sort="TSORT" 
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort 
make
lpi_deck $build_dir $src_dir

# PrgEnv-cray, OpenMP, LSORT 
build_dir=$openmp/cray/lsort
mkdir -p $build_dir
cd $build_dir
cmake_pthreads="OFF"
cmake_openmp="ON"
cmake_sort="LSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-gnu, OpenMP, TSORT 
build_dir=$openmp/gnu/tsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-gnu
cmake_pthreads="OFF"
cmake_openmp="ON"
cmake_sort="TSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-gnu, OpenMP, LSORT 
build_dir=$openmp/gnu/lsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-gnu
cmake_pthreads="OFF"
cmake_openmp="ON"
cmake_sort="LSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-aocc, OpenMP, TSORT 
build_dir=$openmp/aocc/tsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-aocc
cmake_pthreads="OFF"
cmake_openmp="ON"
cmake_sort="TSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-aocc, OpenMP, LSORT 
build_dir=$openmp/aocc/lsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-aocc
cmake_pthreads="OFF"
cmake_openmp="ON"
cmake_sort="LSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-intel, OpenMP, TSORT 
build_dir=$openmp/intel/tsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-intel
cmake_pthreads="OFF"
cmake_openmp="ON"
cmake_sort="TSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-intel, OpenMP, LSORT
build_dir=$openmp/intel/lsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-intel
cmake_pthreads="OFF"
cmake_openmp="ON"
cmake_sort="LSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-cray, PTHREADS, TSORT
build_dir=$pthreads/cray/tsort 
mkdir -p $build_dir
cd $build_dir
cmake_pthreads="ON"
cmake_openmp="OFF"
cmake_sort="TSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-cray, PTHREADS, LSORT
build_dir=$pthreads/cray/lsort
mkdir -p $build_dir
cd $build_dir
cmake_pthreads="ON"
cmake_openmp="OFF"
cmake_sort="LSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-gnu, PTHREADS, TSORT
build_dir=$pthreads/gnu/tsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-gnu
cmake_pthreads="ON"
cmake_openmp="OFF"
cmake_sort="TSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-gnu, PTHREADS, LSORT
build_dir=$pthreads/gnu/lsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-gnu
cmake_pthreads="ON"
cmake_openmp="OFF"
cmake_sort="LSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-aocc, PTHREADS, TSORT 
build_dir=$pthreads/aocc/tsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-aocc
cmake_pthreads="ON"
cmake_openmp="OFF"
cmake_sort="TSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-aocc, PTHREADS, LSORT
build_dir=$pthreads/aocc/lsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-aocc
cmake_pthreads="ON"
cmake_openmp="OFF"
cmake_sort="LSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-intel, PTHREADS, TSORT
build_dir=$pthreads/intel/tsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-intel
cmake_pthreads="ON"
cmake_openmp="OFF"
cmake_sort="TSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

# PrgEnv-intel, PTHREADS, LSORT
build_dir=$pthreads/intel/lsort
mkdir -p $build_dir
cd $build_dir
module swap PrgEnv-cray PrgEnv-intel
cmake_pthreads="ON"
cmake_openmp="OFF"
cmake_sort="LSORT"
vpic_cmake $src_dir $cmake_pthreads $cmake_openmp $cmake_sort
make
lpi_deck $build_dir $src_dir

