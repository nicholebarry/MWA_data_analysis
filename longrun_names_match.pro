pro longrun_names_match,obs_names=obs_names

  ;Long run directory
  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  
  ;Use this as a proxy to get the obsids that were run
  ;vis_files=findfile(dir+'vis_data/106*_vis_XX.sav')
  ;obsids=file_basename(vis_files,'_vis_XX.sav')
  
  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/beardsley_thesis_list.txt'
  readcol, filename, obsids, format='A', /silent
  
  ;obs_names=STRARR(2441,3)
  obs_names=STRARR(1029,3)
  obs_names[*,0]=obsids
  
  day_name_array=['Aug23','Aug27','Sep02','Sep04','Sep06','Sep09','Sep11','Sep13','Sep17','Sep19','Sep30','Oct02','Oct04','Oct08','Oct10','Oct15','Oct23','Oct25','Oct28','Oct29','Oct31','Nov17','Nov18','Nov29','Nov30']
  pointing_name_array=['minusfour','minusthree','minustwo','minusone','zenith','plusone','plustwo','plusthree','plusfour']
  
  for day_i=0, N_elements(day_name_array)-1 do begin
    for poi_i=0, N_elements(pointing_name_array)-1 do begin
    
      filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/uncut/'+day_name_array[day_i]+'_'+pointing_name_array[poi_i]+'.txt'
      obsids='empty'
      readcol, filename, obsids, format='A', /silent
      ;If day_name_array[day_i] EQ 'Nov18' then stop
      if obsids[0] EQ 'empty' then continue
      
      for obs_i=0, N_elements(obsids)-1 do begin
      
        obs_index = WHERE(STRMATCH(obs_names[*,0],obsids[obs_i]) EQ 1, n_count)
        ;If day_name_array[day_i] EQ 'Nov18' then stop
        If n_count GT 0 then begin
          obs_names[obs_index,1] = day_name_array[day_i]
          obs_names[obs_index,2] = pointing_name_array[poi_i]
          ;If pointing_name_array[poi_i] EQ 'minusfour' OR pointing_name_array[poi_i] EQ 'minusthree' then begin
          ;  minus43 = [minus43,obs_index]
        endif
        
      endfor
      
    endfor
  endfor
  
  return
end