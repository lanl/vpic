#!/bin/bash

master=$HOME/VPIC/vpic-master
scripts=$master/NewTestScripts
scratch=/lustre/scratch5/.mdt0/matsekh/VPIC
builds=$scratch/test-builds
results=$scratch/NewTestResults

function vpic_benchmark(){

    prg_env=$1 
    dir_path=$2

    # set up directory tree 
    test_dir=$results/$dir_path
    mkdir -p $test_dir
    cd $test_dir
    mkdir -p strong
    mkdir -p thread
    mkdir -p weak
    build_dir=$builds/$dir_path

    # set up strong scaling test
    cd $test_dir/strong
    deck=$build_dir/lpi_2d_F6_test.Linux
    # launch strong scaling test
    sbatch -N 1 -t 10:00:00 $scripts/vpic_strong_scaling.sh $prg_env $deck

    # set up thread scaling test (const cpu)
    cd $test_dir/thread
    deck=$build_dir/lpi_2d_F6_test.Linux
    # launch thread scaling test (const cpu) 
    sbatch -N 1 -t 10:00:00 $scripts/vpic_thread_scaling.sh $prg_env $deck

    # set up weak scaling test
    cd $test_dir/weak
    deck1=$build_dir/lpi_2d_F6_test-nx864.Linux
    deck2=$build_dir/lpi_2d_F6_test-nx288.Linux
    deck3=$build_dir/lpi_2d_F6_test-nx96.Linux
    deck4=$build_dir/lpi_2d_F6_test-nx48.Linux
    deck5=$build_dir/lpi_2d_F6_test-nx24.Linux
    deck6=$build_dir/lpi_2d_F6_test-nx12.Linux
    # launch weak scaling test
    sbatch -N 1 -t 10:00:00 $scripts/vpic_weak_scaling.sh $prg_env $deck1 $deck2 $deck3 $deck4 $deck5 $deck6
}

 prog_env="cray"
 bench_path=openmp/cray/tsort
 vpic_benchmark $prog_env $bench_path  
 
 prog_env="cray"
 bench_path=openmp/cray/lsort
 vpic_benchmark $prog_env $bench_path  
 
 prog_env="aocc"
 bench_path=openmp/aocc/tsort
 vpic_benchmark $prog_env $bench_path  
 
 prog_env="aocc"
 bench_path=openmp/aocc/lsort
 vpic_benchmark $prog_env $bench_path  

 prog_env="gnu"
 bench_path=openmp/gnu/tsort
 vpic_benchmark $prog_env $bench_path  

 prog_env="gnu"
 bench_path=openmp/gnu/lsort
 vpic_benchmark $prog_env $bench_path  

 prog_env="intel"
 bench_path=openmp/intel/tsort
 vpic_benchmark $prog_env $bench_path  
 
 prog_env="intel"
 bench_path=openmp/intel/lsort
 vpic_benchmark $prog_env $bench_path  
 
 prog_env="cray"
 bench_path=pthreads/cray/tsort
 vpic_benchmark $prog_env $bench_path  
 
 prog_env="cray"
 bench_path=pthreads/cray/lsort
 vpic_benchmark $prog_env $bench_path  
 
 prog_env="aocc"
 bench_path=pthreads/aocc/tsort
 vpic_benchmark $prog_env $bench_path  
 
 prog_env="aocc"
 bench_path=pthreads/aocc/lsort
 vpic_benchmark $prog_env $bench_path  

 prog_env="gnu"
 bench_path=pthreads/gnu/tsort
 vpic_benchmark $prog_env $bench_path  

 prog_env="gnu"
 bench_path=pthreads/gnu/lsort
 vpic_benchmark $prog_env $bench_path  

 prog_env="intel"
 bench_path=pthreads/intel/tsort
 vpic_benchmark $prog_env $bench_path  
 
 prog_env="intel"
 bench_path=pthreads/intel/lsort
 vpic_benchmark $prog_env $bench_path  
