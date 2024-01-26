#!/bin/bash
#PBS -N P
#PBS -q luna
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -l walltime=10:59:59
#PBS -l place=group=infiniband
##PBS -o stdout
##PBS -j oe
##PBS -m ae

module add python-3.6.2-gcc

export OMP_NUM_THREADS=$TORQUE_RESC_TOTAL_PROCS

 cd $PBS_O_WORKDIR
echo ""
echo "OMP_NUM_THREADS:"
echo $OMP_NUM_THREADS 
echo ""

echo "Beginning the test for PPSTM_test.py script"
python PPSTM_simpleUpdated2.py > out.txt

echo "Everything seems to be done, so... bye!"

