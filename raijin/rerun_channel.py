#!/usr/bin/env python
import os
import time
import glob
from subprocess import check_output, call


if __name__ == '__main__':
	#check outputs for errored logs
        coarse_chans = []
	OUTPUTDIR = '/g/data1b/ru04/NB_CHIPS_OUT/'
	coarse_channel=['12']
	output_tag='FHD_1_minusone'
        freq_key = open(OUTPUTDIR + '/freq_key_' + coarse_channel +'.txt')
        for freq_coarse in freq_key:
        	os.remove(OUTPUTDIR + '/bv_freq' + freq_coarse + '_*._' + output_tag + '.dat')
                os.remove(OUTPUTDIR + '/bvdiff_freq' + freq_coarse + '_*._' + output_tag + '.dat')
                os.remove(OUTPUTDIR + '/noise_freq' + freq_coarse + '_*._' + output_tag + '.dat')
                os.remove(OUTPUTDIR + '/noisediff_freq' + freq_coarse + '_*._' + output_tag + '.dat')
                os.remove(OUTPUTDIR + '/weights_freq' + freq_coarse + '_*._' + output_tag + '.dat')

	obs_list = open('/home/563/nb9897/MWA/chips_2019/scripts/obs_list/1_minusone.txt','r').read().split('\n')
        obs_list = array([obs for obs in obs_list if obs != '' and obs != ' ' ])

	for obs in obs_list:
		stripped_grid_job='/short/ru04/nb9897/logs/1_minusone/run_grid_'+obs+'_'+output_tag
		grid_job_message = check_output('qsub %s_%s.sh' %(stripped_grid_job,coarse_channel),shell=True)
		grid_job_ID = grid_job_message.split()[-1].strip('\n')
		grid_jobs_ID_grep.append(' -e ')
		grid_jobs_ID_grep.append(grid_job_ID)
		grid_jobs_ID_grep_join = ''.join(grid_jobs_ID_grep))
		wait_job = check_output('qstat | grep -e %s | wc -l' %(grid_jobs_ID_grep_join),shell=True)
                while int(wait_job) > 1:
                	time.sleep(30)
                	wait_job = check_output('qstat | grep -e %s | wc -l' %(grid_jobs_ID_grep_join),shell=True)

