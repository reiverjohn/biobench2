#!/bin/sh -login
#PBS -q main
#PBS -j oe
#PBS -m abe
#PBS -M johnj@msu.edu
#PBS -l nodes=1:ppn=1,walltime=04:00:00,mem=8gb

### Set environment and modules here  MODIFY AS NEEDED!!!
module swap intelcc gnu
module load gcc
module unload mkl mvapich 

cd /mnt/home/john1545/BioBench2
perl ./runBench.pl
