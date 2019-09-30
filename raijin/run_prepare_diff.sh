#!/bin/bash -l
#PBS -N "lssa_2"
#PBS -l walltime=1:00:00
#PBS -l ncpus=16
#PBS -l mem=10GB
#PBS -l jobfs=0MB
#PBS -l wd
#PBS -l software=CHIPS_lssa
#PBS -q normal
#PBS -P ru04
module load gcc/6.2.0 openmpi/3.0.0
source /home/563/nb9897/MWA/chips_2019/scripts/raijin_env_variables.sh
echo $PBS_JOBID
export OMP_NUM_THREADS=16
printenv
cd $CODEDIR
./prepare_diff _FHD_19_minusone 236 0 'xx' _FHD_19_minusone_redshift7 1 -c 80000.00000 -p 2.000 -n 168515000.00000
./prepare_diff _FHD_1_plusone 384 0 'yy' _FHD_1_plusone 1 -c 80000.00000 -p 2.000 -n 167075000.00000

./prepare_diff _RTS_${num}_minusone 236 0 'yy' _RTS_${num}_minusone_redshift7 1 -c 80000.00000 -p 8.000 -n 168475000.00000
./prepare_diff _RTS_${num}_minustwo 384 0 'xx' _RTS_${num}_minustwo 1 -c 80000.00000 -p 8.000 -n 167035000.00000
