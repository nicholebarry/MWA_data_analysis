pro cable_pointing_cal_maker

  ;**************setup
  ;set to the default if day is not set
  If ~keyword_set(day) then day='beardsley_thesis_list_'
  
  ;parsednames=[day+'minusthree',day+'minustwo',day+'minusone',day+'zenith',day+'plusone',day+'plustwo',day+'plusthree']
  parsednames=[day+'minustwo',day+'minusone',day+'zenith',day+'plusone',day+'plustwo']
  
  ;obs ptr is where obs of a pointing are stored. It can be made bigger than needed to account for the pointings of the longrun
  obs_ptr=PTRARR(10,/allocate_heap)
  
  ;Different pointing ranges given if it was the longrun or not
  ;pointing_num=[-3,-2,-1,0,1,2,3]
  pointing_num=[-2,-1,0,1,2]
  
  ;For each pointing, get the obs from the correct text file and read it to an obs_ptr array. Fill variable with 'empty' string first to create a tag if it is unfilled by
  ;the text file. Save the number of obs in that pointing.
  FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
  
  
    filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/' + parsednames[j] + '.txt'
    obs_temp='empty'
    readcol, filename, obs_temp, format='A', /silent
    If j NE 0 then parsednumbers=[parsednumbers,N_elements(obs_temp)] else parsednumbers=N_elements(obs_temp)
    *obs_ptr[j]=obs_temp
    
  ENDFOR
  ;*************end of setup
  
  obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/'+(*obs_ptr[0])[0]+'_obs.sav','obs')
  cal = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/calibration/'+(*obs_ptr[0])[0]+'_cal.sav','cal')
  mode_filepath=filepath(obs.instrument+'_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  
  ;Taking tile information and cross-matching it with the nonflagged tiles array, resulting in nonflagged tile arrays
  ;grouped by cable length
  cable_length_ref=cable_len[Uniq(cable_len,Sort(cable_len))]
  n_cable=N_Elements(cable_length_ref)
  
  n_pol=cal.n_pol
  n_freq=cal.n_freq
  n_tile=cal.n_tile
  freq_arr=cal.freq
  freq_use = where((*obs.baseline_info).freq_use)
  nf_use=N_Elements(freq_use)
  
  FOR j=0, (size(pointing_num))[1]-1 DO BEGIN
    obsid=(*obs_ptr[j])
    tile_use = PTRARR(N_elements(obsid),/allocate)
    raw_gain = PTRARR(2,N_elements(obsid),/allocate)
    FOR obs_i=0, N_elements(obsid)-1 do begin
      cal = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/calibration/'+obsid[obs_i]+'_cal.sav','cal')
      obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/'+obsid[obs_i]+'_obs.sav','obs')
      
      *tile_use[obs_i]=where((*obs.baseline_info).tile_use)
      FOR pol_i=0, 1 do begin
        *raw_gain[pol_i,obs_i] = *cal.gain_residual[pol_i]+*cal.gain[pol_i]
      endfor
    ENDFOR
    
    bandpass_arr=Fltarr(N_elements(obsid),(n_pol)*n_cable+1,n_freq)
    
    FOR obs_i=0, N_elements(obsid)-1 do begin
      tile_use_arr=Ptrarr(n_cable)
      FOR cable_i=0,n_cable-1 DO tile_use_arr[cable_i]=Ptr_new(where(*tile_use[obs_i] AND cable_len EQ cable_length_ref[cable_i]))
      
      
      bandpass_arr[obs_i,0,*]=freq_arr
      bandpass_col_count=0
      
      ;Main gain calculation loop
      FOR cable_i=0,n_cable-1 DO BEGIN
        tile_use_cable=*tile_use_arr[cable_i]
        
        nt_use_cable=N_Elements(tile_use_cable)
        
        FOR pol_i=0,n_pol-1 DO BEGIN
          gain=*raw_gain[pol_i, obs_i] ;n_freq x n_tile element complex array
          
          ;gain2 is a temporary variable used in place of the gain array for an added layer of safety
          IF cable_i EQ 0 AND pol_i EQ 0 THEN gain2=Complexarr(n_pol,(size(gain))[1],(size(gain))[2])
          
          ;Only use gains from unflagged tiles and frequencies, and calculate the amplitude and phase
          gain_use=extract_subarray(gain,freq_use,tile_use_cable)
          amp=Abs(gain_use)
          phase=Atan(gain_use,/phase)
          
          ;amp2 is a temporary variable used in place of the amp array for an added layer of safety
          amp2=fltarr(nf_use,nt_use_cable)
          
          ;This is the normalization loop for each tile. If the mean of gain amplitudes over all frequencies is nonzero, then divide
          ;the gain amplitudes by that number, otherwise make the gain amplitudes zero.
          FOR tile_i=0,nt_use_cable-1 DO BEGIN
            resistant_mean,amp[*,tile_i],2,res_mean
            IF res_mean NE 0 THEN amp2[*,tile_i]=amp[*,tile_i]/res_mean ELSE amp2[*,tile_i]=0.
          ENDFOR
          
          ;This finds the normalized gain amplitude mean per frequency over all tiles, which is the final bandpass per cable group.
          
          bandpass_single=Fltarr(nf_use)
          FOR f_i=0L,nf_use-1 DO BEGIN
            resistant_mean,amp2[f_i,*],2,res_mean
            bandpass_single[f_i]=res_mean
          ENDFOR
          
          ;Want iterative to start at 1 (to not overwrite freq) and store final bandpass per cable group.
          bandpass_col_count += 1
          bandpass_arr[obs_i,bandpass_col_count,freq_use]=bandpass_single
          
        endfor
      endfor
    endfor
    
    bandpass_arr2=Fltarr((n_pol)*n_cable+1,n_freq)
    bandpass_arr2[0,*]=freq_arr
    bandpass_col_count=0
    
    for col_i=1, (n_pol)*n_cable do begin
      bandpass_single=Fltarr(nf_use)
      FOR f_i=0L,nf_use-1 DO BEGIN
        resistant_mean,bandpass_arr[*,col_i,freq_use[f_i]],2,res_mean
        bandpass_single[f_i]=res_mean
      ENDFOR
      
      ;n_freq x 13 array. columns are frequency, 90m xx, 90m yy, 150m xx, 150m yy, 230m xx, 230m yy, 320m xx, 320m yy, 400m xx, 400m yy, 524m xx, 524m yy
      bandpass_col_count += 1
      bandpass_arr2[bandpass_col_count,freq_use]=bandpass_single
      
      textfast, bandpass_arr2,/write,append=0,file_path='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/FHD/instrument_config/'+strtrim(pointing_num[j],2)+'_bandpass_2013longrun_Jan2017.txt'
      
    endfor
    
  ENDFOR
  
end
