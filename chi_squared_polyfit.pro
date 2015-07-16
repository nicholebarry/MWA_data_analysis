pro chi_squared_polyfit, data_filename, param_num,split_jump=split_jump,nopointing=nopointing, unreduced=unreduced
  ;Routine for calculating the reduced chi-squared statistic of polyfits and modefits to the cross gains for each tile.  Outputs are printed min,
  ;max, and mean reduced chi-squared statistic over all tiles, with a sav file including reduced chi-squared statistics for each tile per polarization
  ;and per pointing in the directory specified by data_filename.
  ;data_filename is the location of the polyfits and modefits. param_num is the number of degrees of freedom of the fit per tile. split_jump is
  ;an indication that the fit as a split in the files to accomadate the digital gain jump. mode_calc tells the program to pick out the tiles with
  ;150m cables to calculate the chi-squared statistic. nopointing tells the program to modify the degrees of freedom to accomated runs where the
  ;modefit was not calculated by pointing.


  ;**************setup, stolen from mkbp_pointings
  ;set to the default if day is not set
  If ~keyword_set(day) then day='Aug23'
  
  ;I accidently named the longrun txt files full of correct obs slightly differently
  If keyword_set(longrun) then parsednames=[day+'_minusfour',day+'_minusthree',day+'_minustwo',day+'_minusone',day+'_zenith',day+'_plusone',day+'_plustwo',day+'_plusthree',day+'_plusfour'] else $
    parsednames=[day+'minustwo',day+'minusone',day+'zenith',day+'plusone',day+'plustwo',day+'plusthree']
    
  ;obs ptr is where obs of a pointing are stored. It can be made bigger than needed to account for the pointings of the longrun
  obs_ptr=PTRARR(10,/allocate)
  
  ;Different pointing ranges given if it was the longrun or not
  If keyword_set(longrun) then pointing_num=[-4,-3,-2,-1,0,1,2,3,4] else pointing_num=[-2,-1,0,1,2,3]
  
  ;For each pointing, get the obs from the correct text file and read it to an obs_ptr array. Fill variable with 'empty' string first to create a tag if it is unfilled by
  ;the text file. Save the number of obs in that pointing.
  FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
  
    If keyword_set(longrun) then filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/' + parsednames[j] + '.txt' else $
      filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/' + parsednames[j] + '.txt'
    obs_temp='empty'
    readcol, filename, obs_temp, format='A', /silent
    If j NE 0 then parsednumbers=[parsednumbers,N_elements(obs_temp)] else parsednumbers=N_elements(obs_temp)
    
    If j NE 0 then all_obs=[all_obs,obs_temp] else all_obs=obs_temp
    
    *obs_ptr[j]=obs_temp
    
  ENDFOR
  ;*************end of setup
  
  ;getting cable lengths
  mode_filepath=filepath('mwa_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  
  ;Defining arrays for later
  tile_use_ptr=PTRARR(94,/allocate)
  freq_use=PTRARR(94,/allocate)
  residuals=FLTARR(94,2,384,128)

  
  ;***************Residual calculation over all obs loop
  for obs_i=0, N_elements(all_obs)-1 do begin
  
    restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_pointing_May2015/metadata/'+strtrim(string(all_obs[obs_i]),2)+'_obs.sav'
    *tile_use_ptr[obs_i]=(*obs.baseline_info).tile_use
    
    tile_index=FINDGEN(128)
    freq_use=where((*obs.baseline_info).freq_use)
    raw_temp=FLTARR(2,384,128)
    fit_temp=FLTARR(2,384,128)
    
    
    ;saved_run location
    ;print, 'Saved run activated!'
    
    IF ulong(all_obs[obs_i]) LE 1061313496 THEN pointing='-2'
    IF (ulong(all_obs[obs_i]) LE 1061315320) AND (ulong(all_obs[obs_i]) GE 1061313616) THEN pointing='-1'
    IF (ulong(all_obs[obs_i]) LE 1061317152) AND (ulong(all_obs[obs_i]) GE 1061315448) THEN pointing='0'
    IF (ulong(all_obs[obs_i]) LE 1061318984) AND (ulong(all_obs[obs_i]) GE 1061317272) THEN pointing='1'
    IF (ulong(all_obs[obs_i]) LE 1061320816) AND (ulong(all_obs[obs_i]) GE 1061319104) THEN pointing='2'
    IF (ulong(all_obs[obs_i]) LE 1061323008) AND (ulong(all_obs[obs_i]) GE 1061320936) THEN pointing='3'
    
    filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_pointing_gainrescor_plots/fhd_nb_reg_Aug23_nopoly_fromnocable/'+pointing+'_bandpass.txt
    readcol, filename, freq_arr_input, cable90xx, cable90yy, cable150xx, cable150yy, cable230xx, cable230yy, cable320xx, cable320yy, cable400xx, cable400yy, cable524xx, cable524yy, /silent
    
    
    restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_pointing_May2015/calibration/'+strtrim(string(all_obs[obs_i]),2)+'_cal.sav'
    for pol_i=0,1 do begin
      *cal.gain[pol_i]=*cal.gain[pol_i]+*cal.gain_residual[pol_i]
      
      for tile_i=0,127 do begin
      
        if (cable_len[tile_i] EQ 90) AND (pol_i EQ 0) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable90xx
        if (cable_len[tile_i] EQ 90) AND (pol_i EQ 1) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable90yy
        if (cable_len[tile_i] EQ 150) AND (pol_i EQ 0) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable150xx
        if (cable_len[tile_i] EQ 150) AND (pol_i EQ 1) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable150yy
        if (cable_len[tile_i] EQ 230) AND (pol_i EQ 0) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable230xx
        if (cable_len[tile_i] EQ 230) AND (pol_i EQ 1) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable230yy
        if (cable_len[tile_i] EQ 320) AND (pol_i EQ 0) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable320xx
        if (cable_len[tile_i] EQ 320) AND (pol_i EQ 1) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable320yy
        if (cable_len[tile_i] EQ 400) AND (pol_i EQ 0) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable400xx
        if (cable_len[tile_i] EQ 400) AND (pol_i EQ 1) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable400yy
        if (cable_len[tile_i] EQ 524) AND (pol_i EQ 0) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable524xx
        if (cable_len[tile_i] EQ 524) AND (pol_i EQ 1) then (*cal.gain[pol_i])[*,tile_i]=(*cal.gain[pol_i])[*,tile_i]/cable524yy
        
      endfor
      
      raw_temp[pol_i,*,*]=abs(*cal.gain[pol_i])
    endfor
    
    if keyword_set(split_jump) then begin
      restore, data_filename+'/'+strtrim(string(all_obs[obs_i]),2)+'pre_cal.sav'
      ;for pol_i=0,1 do fit_temp[pol_i,0:255,*]=abs((*cal.gain[pol_i])[0:255,*])
      for pol_i=0,1 do fit_temp[pol_i,0:255,*]=abs((*cal_pre.gain[pol_i])[0:255,*])
      restore, data_filename+'/'+strtrim(string(all_obs[obs_i]),2)+'post_cal.sav'
      ;for pol_i=0,1 do fit_temp[pol_i,256:383,*]=abs((*cal.gain[pol_i])[256:383,*])
      for pol_i=0,1 do fit_temp[pol_i,256:383,*]=abs((*cal_post.gain[pol_i])[256:383,*])
    endif else begin
      restore, data_filename+'/'+strtrim(string(all_obs[obs_i]),2)+'_cal.sav'
      for pol_i=0,1 do fit_temp[pol_i,*,*]=abs(*cal.gain[pol_i])
    endelse
    
    for pol_i=0,1 do residuals[obs_i,pol_i,*,*]=raw_temp[pol_i,*,*]-fit_temp[pol_i,*,*]
    
  endfor
  ;*********************End of residual calculation over all obs loop
  
  ;*********************Calculate reduced chi-squared block
  

  tile_use=INTARR(94)
  tile_use[*]=1
  for obs_j=0, 93 do tile_use= *tile_use_ptr[obs_j]*tile_use

  tile_use = where(tile_use)

  
  reduced_chi=FLTARR(6,2,128)
  for poi_i=0,5 do begin
  
    ;Choose the right set of obsids for the pointing.
    if poi_i EQ 0 then poi_range=[0,parsednumbers[0]-1]
    if poi_i EQ 1 then poi_range=[parsednumbers[0],total(parsednumbers[0:1])-1]
    if poi_i EQ 2 then poi_range=[total(parsednumbers[0:1]),total(parsednumbers[0:2])-1]
    if poi_i EQ 3 then poi_range=[total(parsednumbers[0:2]),total(parsednumbers[0:3])-1]
    if poi_i EQ 4 then poi_range=[total(parsednumbers[0:3]),total(parsednumbers[0:4])-1]
    if poi_i EQ 5 then poi_range=[total(parsednumbers[0:4]),total(parsednumbers[0:5])-1]
    
    
    for pol_i=0,1 do begin
    
      for tile_i=0, N_elements(tile_use)-1 do begin
      
        ;Initialize a new temp frequency variance array for each tile and pol
        freq_variance=FLTARR(N_elements(freq_use))
        
        for freq_i=0, N_elements(freq_use)-1 do begin
        
          ;Find the mean per frequency channel of the calculated residual
          freq_mean=mean(residuals[poi_range[0]:poi_range[1],pol_i,freq_use[freq_i],tile_use[tile_i]])
          ;Then find the variance (Sum of the squared difference between the residual and mean, divided by the number of observations which went into the calculation of the fit).
          freq_variance[freq_i]=Total((residuals[poi_range[0]:poi_range[1],pol_i,freq_use[freq_i],tile_use[tile_i]]-freq_mean)^2.)/(poi_range[1]-poi_range[0])
          
        endfor
        ;Find the unreduced chi squared for by tile, pol, and pointing
        indiv_chi=Total(residuals[poi_range[0]:poi_range[1],pol_i,freq_use,tile_use[tile_i]]^2./freq_variance)
        
        ;Find the reduced chi squared by tile, pol, and pointing
        If ~keyword_set(unreduced) then begin
          If keyword_set(nopointing) then reduced_chi[poi_i,pol_i,tile_i]=indiv_chi/((N_elements(freq_use)*(poi_range[1]-poi_range[0]))-param_num*(poi_range[1]-poi_range[0])) else $
            reduced_chi[poi_i,pol_i,tile_i]=indiv_chi/((N_elements(freq_use)*(poi_range[1]-poi_range[0]))-param_num)
        endif else begin
          reduced_chi[poi_i,pol_i,tile_i]=indiv_chi
        endelse
        
        if cable_len[tile_use[tile_i]] EQ 90 then begin
          if keyword_set(tile90) then tile90=[tile90,tile_use[tile_i]] else tile90=tile_use[tile_i]
        endif
        if cable_len[tile_use[tile_i]] EQ 150 then begin
          if keyword_set(tile150) then tile150=[tile150,tile_use[tile_i]] else tile150=tile_use[tile_i]
        endif
        if cable_len[tile_use[tile_i]] EQ 230 then begin
          if keyword_set(tile230) then tile230=[tile230,tile_use[tile_i]] else tile230=tile_use[tile_i]
        endif
        if cable_len[tile_use[tile_i]] EQ 320 then begin
          if keyword_set(tile320) then tile320=[tile320,tile_use[tile_i]] else tile320=tile_use[tile_i]
        endif
        if cable_len[tile_use[tile_i]] EQ 400 then begin
          if keyword_set(tile400) then tile400=[tile400,tile_use[tile_i]] else tile400=tile_use[tile_i]
        endif
        if cable_len[tile_use[tile_i]] EQ 524 then begin
          if keyword_set(tile524) then tile524=[tile524,tile_use[tile_i]] else tile524=tile_use[tile_i]
        endif
        
      endfor
    endfor
  endfor
  ;********************End of calculated reduced chi squared block
  ;if keyword_set(unreduced) then name='unreduced_chi' else name='reduced_chi'
  ;If keyword_set(mode_calc) then save, reduced_chi, FILENAME=data_filename+'/'+name+'_150mode.sav' else save, reduced_chi, FILENAME=data_filename+'/'+name+'.sav'
  
  print, '90m'
  print, minmax(reduced_chi[tile90])
  print, mean(reduced_chi[tile90],/NAN)
  
  print, '150m'
  print, minmax(reduced_chi[tile150])
  print, mean(reduced_chi[tile150],/NAN)
  
  print, '230m'
  print, minmax(reduced_chi[tile230])
  print, mean(reduced_chi[tile230],/NAN)
  
  print, '320m'
  print, minmax(reduced_chi[tile320])
  print, mean(reduced_chi[tile320],/NAN)
  
  print, '400m'
  print, minmax(reduced_chi[tile400])
  print, mean(reduced_chi[tile400],/NAN)
  
  print, '524m'
  print, minmax(reduced_chi[tile524])
  print, mean(reduced_chi[tile524],/NAN)
  
  
end