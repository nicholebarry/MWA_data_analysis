pro flag_stats

  dir_file = '/fred/oz048/MWA/data/2014/van_vleck_corrected/coarse_corr_no_ao/obs_lists_extra/'
  file_name = ['26']

calc_stats=0
  if keyword_set(calc_stats) then begin

  for file_i=0,N_elements(file_name)-1 do begin

    obsid = dir_file+file_name[file_i]+'.txt'

    if file_test(obsid) then begin
      readcol, obsid, obsids, FORMAT='(A)'
      obsid = obsids
    endif

    n_obs = N_elements(obsid)
    percent_unflagged = FLTARR(n_obs)

    for obs_i=0, n_obs-1 do begin
      if file_test('/fred/oz048/MWA/data/2014/van_vleck_corrected/coarse_corr_no_ao/'+obsid[obs_i]+'.uvfits') then begin
        uvfits_read,hdr,params,layout,vis_arr,vis_weights,file_path_vis='/fred/oz048/MWA/data/2014/van_vleck_corrected/coarse_corr_no_ao/'+obsid[obs_i]+'.uvfits',silent=1
        percent_unflagged[obs_i] = double(total(*vis_weights[0] > 0)) / double(N_elements(*vis_weights[0]))
        print, obsid[obs_i], percent_unflagged[obs_i], format='(A,F10)'
        undefine_fhd, hdr, params, layout, vis_arr, vis_weights
      endif
    endfor
    save, percent_unflagged, obsid, filename='/fred/oz048/MWA/data/2014/van_vleck_corrected/coarse_corr_no_ao/flag_stats_'+file_name[file_i]+'.sav'

  endfor
stop
  endif

  days_to_convert = ['25','26']
  pointings_to_convert = ['_minustwo','_minusone','_zenith','_plusone','_plustwo']
  n_days = N_elements(days_to_convert)
  n_poi = N_elements(pointings_to_convert)

  for day_i=0, n_days-1 do begin

    restore,'/fred/oz048/MWA/data/2014/van_vleck_corrected/coarse_corr_no_ao/flag_stats_'+days_to_convert[day_i]+'.sav'
    flagged_inds = where(percent_unflagged LT .6)
    unflagged_inds = where(percent_unflagged GT .6)
    obs = obsid[unflagged_inds]

    readcol,dir_file + days_to_convert[day_i]+'.txt',obs_from_txt,format='(A)'
    match,obs,obs_from_txt,suba,subb
    print, N_elements(subb)/float(N_elements(obs_from_txt))
    openw,21,dir_file + days_to_convert[day_i]+'_ssins.txt'
    printf,21, transpose(obs_from_txt[subb]),format='(A)'
    close, 21

    for poi_i=0, n_poi-1 do begin
      readcol,dir_file + days_to_convert[day_i]+pointings_to_convert[poi_i]+'.txt',obs_from_txt,format='(A)'
      match,obs,obs_from_txt,suba,subb
      print, N_elements(subb)/float(N_elements(obs_from_txt))
      openw,21,dir_file + days_to_convert[day_i]+pointings_to_convert[poi_i]+'_ssins.txt'
      printf,21, transpose(obs_from_txt[subb]),format='(A)'
      close, 21
    endfor

  endfor

end
