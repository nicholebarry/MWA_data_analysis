import string

fname='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_presidelobe/grid_out/firstpass.o70326.'
ending = (range(95))[1:]
fname_array=[fname+str(end) for end in ending]

for fname_array_element in fname_array:

	with open(fname_array_element) as f:
    		lines = f.readlines()

		for line in lines:
			#data0 = line.find('OBSID')
			#if data0 != -1:
			#	print line 
			data1 = line.find('from low signal to noise after')
			data2 = line.find('sources detected with maximum signal-to-noise of')
			data3 = line.find('Full pipeline time (minutes):')

			if data1 != -1:
				data_start = line.find('noise after')+12
				data_end = line.find('seconds')
				decon_time = line[data_start:data_end]
				data_start = line.find('with')+4
				data_end = line.find('sources')
				decon_components = line[data_start:data_end]
				
			if data2 != -1:
				data_end = line.find('sources')
				sources = line[0:data_end]
				data_start = line.find('of')+2
				snr = line[data_start:-1]

			if data3 != -1:
				data_start = line.find(':')+1
				pipe_time = line[data_start:-1]

		try:
  			pipe_time
		except NameError:
  			pipe_time = 'Undefined'

		print decon_time + ',' + decon_components + ',' + sources + ',' + snr + ',' + pipe_time

