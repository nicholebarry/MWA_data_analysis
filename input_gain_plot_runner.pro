pro input_gain_plot_runner, day
  ;day options are Aug27, Aug23, both, or both_ratio
  ;Aug23: input cross-calc gains, unscaled auto gains
  ;Aug27: input cross-calc gains, unscaled auto gains
  ;both: input cross-calc gains Aug23, unscaled auto gains Aug23, unscaled auto gains Aug27
  ;both_ratio: ratio of input cross-calc gains Aug23 and unscaled auto gains Aug23, input cross-calc gains Aug27 and unscaled auto gains Aug27

  IF day EQ 'both_ratio' then both_ratio=1
  IF day EQ 'both_ratio' then day='both'
  
  n_pol=2
  
  ;Either find the obsids of Aug23 or Aug27, or find the obsids for both
  If day NE 'both' then begin
    If day EQ 'Aug27' then filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/obs_id_27.txt'
    If day EQ 'Aug23' then filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/obs_id_23.txt'
    readcol, filename, obsid, format='A'
  endif else begin
    filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/obs_id_27.txt'
    readcol, filename, obsidAug27, format='A'
    filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/obs_id_23.txt'
    readcol, filename, obsidAug23, format='A'
    ;Set obsid to get the number of elements right
    obsid=obsidAug23
  endelse
  obs_amount=N_elements(obsid)
  
  ;Define temp gain holders
  gain_x_temp=FLTARR(384,obs_amount,128)
  gain_y_temp=FLTARR(384,obs_amount,128)
  gain_autox_temp=FLTARR(384,obs_amount,128)
  gain_autoy_temp=FLTARR(384,obs_amount,128)
  If day EQ 'both' then gain_autoxAug27_temp=FLTARR(384,obs_amount,128) & If day EQ 'both' then gain_autoyAug27_temp=FLTARR(384,obs_amount,128)
  If keyword_set(both_ratio) then gain_x_temp_ratio=FLTARR(384,obs_amount,128) & If keyword_set(both_ratio) then gain_y_temp_ratio=FLTARR(384,obs_amount,128)
  
  ;restore_path_Aug27='/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/'
  restore_path_Aug27='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_missingobs_Aug27_May2015/calibration/'
  restore_path_Aug23='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_std_May2015/calibration/'
  restore_path_Aug27_auto='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_autogainsonly_Aug27_May2015/calibration/'
  restore_path_Aug23_auto='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_autogainsonly_May2015/calibration/'
  If day EQ 'Aug27' then restore_path=[restore_path_Aug27,restore_path_Aug27_auto]
  If day EQ 'Aug23' then restore_path=[restore_path_Aug23,restore_path_Aug23_auto]
  
  ;At the moment, I think I only want to compare Xgains, AutoAug23, AutoAug27
  If day EQ 'both' then restore_path=[restore_path_Aug23,restore_path_Aug23_auto,restore_path_Aug27_auto]
  If keyword_set(both_ratio) then restore_path=[restore_path_Aug23,restore_path_Aug23_auto,restore_path_Aug27_auto,restore_path_Aug27] ;order specific for below
  
  ;auto_params_Aug23_full=FLTARR(obs_amount,2,2,128) ;obs_amount x pol x linfit x tile
  ;auto_params_Aug27_full=FLTARR(obs_amount,2,2,128)
  n_freq=384
  for obs_i=0, obs_amount-1 do begin
  
    ;Begin cross restore and store
    restore, restore_path[0] + obsid[obs_i] + '_cal.sav'
    for pol_i=0,1 do *cal.gain[pol_i]=*cal.gain[pol_i]+*cal.gain_residual[pol_i]
    
    gain_x_temp[*,obs_i,*]=(*cal.gain[0])[*,*]
    gain_y_temp[*,obs_i,*]=(*cal.gain[1])[*,*]
    
    print, obs_i
    
    ;Begin auto restore and store
    restore, restore_path[1] + obsid[obs_i] + '_cal.sav'
    
    ;gain_autox_temp[*,obs_i,*]=(*cal.gain[0])[*,*]
    ;gain_autoy_temp[*,obs_i,*]=(*cal.gain[1])[*,*]
    for freq_i=0,n_freq-1 do gain_autox_temp[freq_i,obs_i,*]=(abs((*cal.gain[0])[freq_i,*])-(*cal.auto_params[0])[0,*])/(*cal.auto_params[0])[1,*]
    for freq_i=0,n_freq-1 do gain_autoy_temp[freq_i,obs_i,*]=(abs((*cal.gain[1])[freq_i,*])-(*cal.auto_params[1])[0,*])/(*cal.auto_params[1])[1,*]
    ;auto_params_Aug23_full[obs_i,0,*,*]=(*cal.auto_params[0])[*,*]
    ;auto_params_Aug23_full[obs_i,1,*,*]=(*cal.auto_params[1])[*,*]
    ;for freq_i=0,n_freq-1 do gain_autox_temp[freq_i,obs_i,*]=((*cal.gain[0])[freq_i,*])/(*cal.auto_params[0])[1,*]
    ;for freq_i=0,n_freq-1 do gain_autoy_temp[freq_i,obs_i,*]=((*cal.gain[1])[freq_i,*])/(*cal.auto_params[1])[1,*]
    
    If day EQ 'both' then begin
      obsid=obsidAug27
      restore, restore_path[2] + obsid[obs_i] + '_cal.sav'
      for freq_i=0,n_freq-1 do gain_autoxAug27_temp[freq_i,obs_i,*]=(abs((*cal.gain[0])[freq_i,*])-(*cal.auto_params[0])[0,*])/(*cal.auto_params[0])[1,*]
      for freq_i=0,n_freq-1 do gain_autoyAug27_temp[freq_i,obs_i,*]=(abs((*cal.gain[1])[freq_i,*])-(*cal.auto_params[1])[0,*])/(*cal.auto_params[1])[1,*]
      ;auto_params_Aug27_full[obs_i,0,*,*]=(*cal.auto_params[0])[*,*]
      ;auto_params_Aug27_full[obs_i,1,*,*]=(*cal.auto_params[1])[*,*]
      ;for freq_i=0,n_freq-1 do gain_autoxAug27_temp[freq_i,obs_i,*]=((*cal.gain[0])[freq_i,*])/(*cal.auto_params[0])[1,*]
      ;for freq_i=0,n_freq-1 do gain_autoyAug27_temp[freq_i,obs_i,*]=((*cal.gain[1])[freq_i,*])/(*cal.auto_params[1])[1,*]
      obsid=obsidAug23
    endif
    
    If keyword_Set(both_ratio) then begin
      obsid=obsidAug27
      restore, restore_path[3] + obsid[obs_i] + '_cal.sav'
      for pol_i=0,1 do *cal.gain[pol_i]=*cal.gain[pol_i]+*cal.gain_residual[pol_i]
      gain_x_temp_ratio[*,obs_i,*]=(*cal.gain[0])[*,*]
      gain_y_temp_ratio[*,obs_i,*]=(*cal.gain[1])[*,*]
      obsid=obsidAug23
    endif
    
  endfor ;end cal restore loop
  
  ;auto_params_Aug23=mean(auto_params_Aug23_full,dimension=1)
  ;auto_params_Aug27=mean(auto_params_Aug27_full,dimension=1)
  
  ;At the moment, not taking into account flagging (just for ease)
  obs_use_binary_full=INTARR(obs_amount,128)
  obs_use_binary_full[*,*]=1
  
  ;for obs_i=0, obs_amount-1 do begin
  ;
  ;   If keyword_set(Aug27) then restore_path='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_missingobs_Aug27_May2015/metadata/' + obsid[obs_i] + '_obs.sav' $
  ;   else restore_path='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_std_May2015/metadata/' + obsid[obs_i] + '_obs.sav'
  ;
  ;   restore, restore_path
  ;   obs_use_binary_full[obs_i,*]=(*obs.baseline_info).tile_use[*]
  ;
  ; endfor
  
  auto_input=PTRARR(2,/allocate)
  If day EQ 'both' then autoAug27=PTRARR(2,/allocate)
  If keyword_set(both_ratio) then gainAug27_forratio=PTRARR(2,/allocate)
  
  forbidden_tiles=[76,103]
  If day EQ 'Aug27' then forbidden_tiles=[forbidden_tiles,80,17]
  If day EQ 'both' then forbidden_tiles=[forbidden_tiles,80,17]
  
  ;Set tile beginning and end to get a range of plots
  for tile_i=0,127 do begin
    tile_index=tile_i
    tile_name=strtrim(string(cal.tile_names[tile_index]),2)
    
    IF where(tile_i EQ forbidden_tiles) EQ -1 then begin
      ;IF tile_i NE 17 AND tile_i NE 76 AND tile_i NE 80 AND tile_i NE 103 then begin
    
      *cal.gain[0]=gain_x_temp[*,*,tile_index]
      *cal.gain[1]=gain_y_temp[*,*,tile_index]
      *auto_input[0]=gain_autox_temp[*,*,tile_index]
      *auto_input[1]=gain_autoy_temp[*,*,tile_index]
      If day EQ 'both' then begin
        *autoAug27[0]=gain_autoxAug27_temp[*,*,tile_index]
        *autoAug27[1]=gain_autoyAug27_temp[*,*,tile_index]
      endif
      If keyword_set(both_ratio) then begin
        *gainAug27_forratio[0]=gain_x_temp_ratio[*,*,tile_index]
        *gainAug27_forratio[1]=gain_y_temp_ratio[*,*,tile_index]
      endif
      obs_use_binary=obs_use_binary_full[*,tile_index]
      print, 'On tile name '+ tile_name+' index '+strtrim(string(tile_i),2)
      bp=input_gain_plotter(cal,obs,output_file_path='/nfs/eor-00/h1/nbarry/Aug23_Aug27_pertile_ratio/',$
        tile_index=tile_index,tile_name=tile_name,obs_use_binary=obs_use_binary,auto_input=auto_input,$
        autoAug27=autoAug27,gainAug27_forratio=gainAug27_forratio)
    ;auto_params_Aug23=auto_params_Aug23[*,*,tile_i], auto_params_Aug27=auto_params_Aug27[*,*,tile_i])
        
    endif
  endfor
  
  
END