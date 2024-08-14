pro extra_rfi_checker

  dir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_gd_woden_redo/'
  cal_dir = 'calibration/'
  obs_dir = 'metadata/'
  vis_dir = 'vis_data/'
  pol_name=['XX','YY']

  ; Parameters for the observation id files
  dir_obsids = '/home/nbarry/MWA/pipeline_scripts/bash_scripts/ozstar/obs_list/2014/'
  day = ['10','11','13','14','15','16']
  n_days = N_elements(day)
  pointingnames=['_minustwo_ssins','_minusone_ssins','_zenith_ssins','_plusone_ssins','_plustwo_ssins']
  n_pointings = N_elements(pointingnames)
  obs_ptr=PTRARR(n_days,n_pointings,/allocate_heap)

  max_n_obs=0
  ; Get observation ids and put them in a pointing pointer
  FOR poi_i=0,n_pointings-1 DO BEGIN
    max_n=0
    FOR day_i=0,n_days-1 DO BEGIN
      readcol, dir_obsids + day[day_i] + pointingnames[poi_i]+'.txt', obs_temp, format='A', /silent
      *obs_ptr[day_i,poi_i]=obs_temp
      max_n_obs = max([N_elements(obs_temp),max_n_obs])
    ENDFOR
  ENDFOR

  for day_i=0, n_days-1 do begin
    for poi_i=0, n_pointings-1 do begin

      obs_ids = *obs_ptr[day_i,poi_i]
      n_obs = N_elements(obs_ids)

      for obs_i=0,n_obs-1 do begin

        if ~file_test(dir + obs_dir + obs_ids[obs_i] + '_params.sav') then continue
        if file_test('/fred/oz048/MWA/CODE/FHD/fhd_nb_data_gd_woden_redo_redo/vis_data/' + obs_ids[obs_i] + '_flags.sav') then continue

        restore, dir + vis_dir + obs_ids[obs_i] + '_flags.sav'
        restore, dir + obs_dir + obs_ids[obs_i] + '_params.sav'

        for pol_i=0, 1 do begin

          restore, dir + vis_dir + obs_ids[obs_i] + '_vis_' + pol_name[pol_i] + '.sav'

          vis_weights_use=0>*vis_weights[pol_i]<1
          data = abs(*vis_ptr * vis_weights_use)

          inds_nan = where(data EQ 0,n_count)
          data[inds_nan] = !Values.F_NAN

          data_zeromean = data - rebin(reform(mean(data,dimension=1,/nan),1,obs.nbaselines*obs.n_time),obs.n_freq,obs.nbaselines*obs.n_time)
          standard_dev =  stddev(data_zeromean,dimension=2,/nan)
          standard_dev = rebin(reform(standard_dev,obs.n_freq,1),obs.n_freq,obs.nbaselines*obs.n_time)
          standard_dev[inds_nan]=0
          data_zeromean[inds_nan]=0

          inds_out = where(data_zeromean GT standard_dev*6,n_count)
          if n_count gt 0 then begin
            row = inds_out / obs.n_freq
            uniq_row = row[uniq(row,sort(row))]
            (*vis_weights[pol_i])[*,uniq_row] = -1
          endif

        endfor

        save, vis_weights, filename='/fred/oz048/MWA/CODE/FHD/fhd_nb_data_gd_woden_redo_redo/vis_data/' + obs_ids[obs_i] + '_flags.sav'
      endfor
    endfor
  endfor

end