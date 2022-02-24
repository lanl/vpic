#!/bin/bash 

prog_env=$1       # programming environment
code=$2           # executable

nn=1              # number of nodes
threads=8         # number of cpu threads
p_scale=2         # scaling constant
log=./log-strong  # log file

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

let np=32\*$nn       # number of processes / ranks
while [ $np -ge 1 ]
do   
    srun -n $np -c $threads --cpu-bind=cores $code --tpp $threads >> $log
    
    let np=$np/$p_scale
done
