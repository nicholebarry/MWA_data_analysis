pro pfb_edge_compare

  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/beardsley_thesis_list.txt'
  readcol, filename, obsids, format='A', /silent
  
  longrun_names_match, obs_names=obs_names
  day_name_array=['Aug23','Aug27','Sep02','Sep04','Sep06','Sep09','Sep11','Sep13','Sep17','Sep19','Sep30','Oct02','Oct04','Oct08','Oct10','Oct15','Oct23','Oct25','Oct28','Oct29','Oct31','Nov17','Nov18','Nov29','Nov30']
  pointing_name_array=['minustwo','minusone','zenith','plusone','plustwo']
  pfb_set = FLTARR(N_elements(day_name_array),2,336,6)
  
  amp_degree=2
  
  obs = getvar_savefile('/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/metadata/1061316296_obs.sav','obs')
  freq_use=where((*obs.baseline_info).freq_use,nf_use)
  anti_freq_use=where((*obs.baseline_info).freq_use EQ 0,nf_use)
  tile_names=ULONG((*obs.baseline_info).tile_names)
  
  ;Using preexisting file to extract information about which tiles have which cable length
  mode_filepath=filepath('mwa_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  ;Taking tile information and cross-matching it with the nonflagged tiles array, resulting in nonflagged tile arrays
  ;grouped by cable length
  cable_length_ref=cable_len[Uniq(cable_len,Sort(cable_len))]
  n_cable=N_Elements(cable_length_ref)
  tile_use_arr=Ptrarr(n_cable)
  tile_arr=FINDGEN(128)
  FOR cable_i=0,n_cable-1 DO tile_use_arr[cable_i]=Ptr_new(where(cable_len EQ cable_length_ref[cable_i]))
  tile_cable_use = [(*tile_use_arr[0])[0],(*tile_use_arr[1])[0],(*tile_use_arr[2])[0],(*tile_use_arr[3])[0],(*tile_use_arr[4])[0],(*tile_use_arr[5])[0]]
  
  subband_low =  where((findgen(384) mod 16) EQ 0) + 1
  subband_high =  where((findgen(385) mod 16) EQ 0) -2
  
  ;****************************Read in and parse temperature data
  
  ;Find the smaller temperature data file with just the right obs to get temperature data.
  ;To get this file, run this code with make_obs_to_query set and then run query.sh
  filename='/nfs/eor-00/h1/nbarry/longrun_poly_bftemp/obs_queried_v2.txt'
  
  ;Read out temperature data in the form of the string due to the funky format
  textfast,data_array,/read,file_path=filename,string=1
  
  temperature_array=FLTARR(N_elements(day_name_array),128) ;Want day x tile
  
  tile_temp=STRARR(8)
  
  for day_i=0,N_elements(day_name_array)-1 do begin
    day_inds = where((obs_names[*,1] EQ day_name_array[day_i]) AND (obs_names[*,2] EQ 'zenith'),n_count)
    if n_count EQ 0 then continue
    temp_index=-1
    for inday_i=0, n_count-1 do begin
      day_index=where(strmatch(data_array,'*'+obs_names[day_inds[inday_i],0]+'*') EQ 1,n_count2)
      if n_count2 GT 0 then temp_index = day_index
    endfor
    if temp_index[0] LT 0 then continue
    split_temp_array=STRARR(N_elements(temp_index),3)
    for i=0,N_elements(temp_index)-1 do begin
      split_temp_array[i,*]=strsplit(data_array[temp_index[i]], '|',/EXTRACT)  ;Split line in data file into receiver, obsid full, and temps per tile
      split_temp_array[i,1]=(strsplit(split_temp_array[i,1], '|',/EXTRACT))[0]  ;Remove extra ms from end of obsid
      split_temp_array[i,0]=STRING(ULONG(split_temp_array[i,0])*10)  ;Add zero to end of receiver for getting tile names easier
      split_temp_array[i,2]=strmid(split_temp_array[i,2],1) ;Remove '{' from beginning of temp data
      str_length=strlen(split_temp_array[i,2]) ;Find length of string for next command
      split_temp_array[i,2]=strmid(split_temp_array[i,2],0,str_length-1) ;Remove '}' from end of temp data
      tile_temp[*]=strsplit(split_temp_array[i,2], ',',/EXTRACT) ;Split up temp data into 8 discrete temperatures per receiver
      
      
      
      for tile_i=0,7 do begin
        tile_index=where(tile_names EQ (ULONG(split_temp_array[i,0])+tile_i+1))
        tile_temp_float=Convert_To_Type(tile_temp[tile_i],4)
        temperature_array[day_i,tile_index]=tile_temp_float ;Add per tile data into array as a float
      endfor
    endfor
    
  endfor
  
  ;****************************End of read in and parse temperature data
  stop
  
  dig_jump_edges = FLTARR(2,N_elements(day_name_array),n_cable)
  edges = FLTARR(2,N_elements(day_name_array),24,n_cable)
  
  for day_i=0, N_elements(day_name_array)-1 do begin
    print, day_i
    obsid_inds = where((obs_names[*,1] EQ day_name_array[day_i]) AND (obs_names[*,2] EQ 'zenith'),n_count)
    if n_count GT 0 then obsids = reform(obs_names[obsid_inds,0]) else continue
    
    ;filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/Oct10_zenith.txt'
    ;readcol, filename, obsids, format='A', /silent
    n_obs = N_elements(obsids)
    
    
    all_subbands = FLTARR(2,n_obs,24,14,6) ;pol,obs,subband num, num of freq channels in subband, cable
    
    gain = (FLTARR(2,n_obs,384,128))
    ave_gain = (FLTARR(2,n_obs,384,n_cable))
    for pol_i=0,0 do begin
      for obs_i=0, n_obs-1 do begin
        cal = getvar_savefile('/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/calibration/'+strtrim(obsids[obs_i],2)+'_cal.sav','cal')
        (*cal.gain[pol_i]) = (*cal.gain[pol_i]) + (*cal.gain_residual[pol_i])
        
        ;        for tile_j=0, N_elements(*tile_use_arr[1])-1 do begin
        ;          tile_i = (*tile_use_arr[1])[tile_j]
        ;          if (cal.mode_params[pol_i,tile_i] NE !NULL) then begin
        ;            mode_i = (*cal.mode_params[pol_i,tile_i])[0]
        ;            amp_use = (*cal.mode_params[pol_i,tile_i])[1]
        ;            phase_use = (*cal.mode_params[pol_i,tile_i])[2]
        ;            gain_mode_fit=amp_use*exp(-Complex(0,1)*2.*!Pi*(mode_i*findgen(384)/384)+Complex(0,1)*phase_use)
        ;            (*cal.gain[pol_i])[*,tile_i] -= gain_mode_fit
        ;          endif
        ;        endfor
        
        gain_fit = FLTARR(384,128)
      endfor
    endfor
  endfor
  
end