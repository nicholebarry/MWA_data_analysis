#!/bin/bash -l
#PBS -r y
#PBS -N "combine_xx"
#PBS -l walltime=03:00:00
#PBS -l ncpus=16
#PBS -l mem=20GB
#PBS -l wd
#PBS -l software=CHIPS_lssa
#PBS -P ru04
#PBS -q express

module load gcc/6.2.0 openmpi/3.0.0
source /home/563/nb9897/MWA/chips_2019/scripts/raijin_env_variables_RTS.sh
echo $PBS_JOBID
export OMP_NUM_THREADS=16
printenv
cd /home/563/nb9897/MWA/chips_2019/bin/

./combine_data '/home/563/nb9897/MWA/chips_2019/scripts/obs_list/combine_lists/RTS/all_plusone_yy.txt' 384 yy._RTS_all_plusone 0
./combine_data '/home/563/nb9897/MWA/chips_2019/scripts/obs_list/combine_lists/RTS_redshift7/678_yy.txt' 236 yy._RTS_all_678_redshift7 0

#usage: ./combine_data "input text file of extensions to be added" Nchan output_extension 0=add/1=subtract
