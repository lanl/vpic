#!/bin/bash 

prog_env=$1       # programming environment

nn=1              # number of nodes
threads=8         # number of cpu threads
p_scale=2         # scaling constant
log=./log-weak    # log file

let np=32\*$nn     # number of processes / ranks

if [ $prog_env = "aocc" ] 
then    
    module swap PrgEnv-cray PrgEnv-aocc
elif [ $prog_env = "gnu" ] 
then
    module swap PrgEnv-cray PrgEnv-gnu
elif [ $prog_env = "intel" ]
then
    module swap PrgEnv-cray PrgEnv-intel
fi

module list
cc --version
CC --version 
#lscpu
#env 

function weak_run(){

    srun -n $1 -c $2 --cpu-bind=cores $3 --tpp $2 >> $log
}

code=$2             # executable
let np=$np/$p_scale
weak_run $np $threads $code

code=$3             # executable
let np=$np/$p_scale
weak_run $np $threads $code

code=$4             # executable 
let np=$np/$p_scale
weak_run $np $threads $code

code=$5             # executable
let np=$np/$p_scale
weak_run $np $threads $code

code=$6             # executable
let np=$np/$p_scale
weak_run $np $threads $code

code=$7             # executable
let np=$np/$p_scale
weak_run $np $threads $code

