#!/bin/bash -l
#PBS -N "lssa_xx"
#PBS -l walltime=7:00:00
#PBS -l ncpus=16
#PBS -l mem=20GB
#PBS -l jobfs=0MB
#PBS -l wd
#PBS -l software=CHIPS_lssa
#PBS -q express
#PBS -P ru04
module load gcc/6.2.0 openmpi/3.0.0
source /home/563/nb9897/MWA/chips_2019/scripts/raijin_env_variables_RTS.sh
echo $PBS_JOBID
export OMP_NUM_THREADS=16
printenv
cd $CODEDIR


./lssa_fg_thermal _RTS_all_678_redshift7 236 80 'xx' 300. _RTS_all_678_redshift7 1 1 0 -c 80000.00000 -p 8.000
./lssa_fg_thermal _RTS_all_plusone 384 80 'yy' 300. _RTS_all_plusone 0 1 0 -c 80000.00000 -p 8.000


#extension, Nchan, nbins, pol, maxu, date, bias_mode, band, bandwidth
