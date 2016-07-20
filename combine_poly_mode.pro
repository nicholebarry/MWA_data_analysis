pro combine_poly_mode, remove_poly=remove_poly, seperate_phase=seperate_phase

  day='Aug23'
  n_pol=2
  n_freq=384
  parsednames=[day+'minustwo',day+'minusone',day+'zenith',day+'plusone',day+'plustwo',day+'plusthree']
  pointing_num=[-2,-1,0,1,2,3]
  
  FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
  
  
    filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/' + parsednames[j] + '.txt'
    obs_temp='empty'
    readcol, filename, obs_temp, format='A', /silent
    If j NE 0 then parsednumbers=[parsednumbers,N_elements(obs_temp)] else parsednumbers=N_elements(obs_temp)
    if j EQ 0 then obs_name_arr=obs_temp else obs_name_arr=[obs_name_arr,obs_temp]
    
  ENDFOR
  
  
  For obs_i=0, N_elements(obs_name_arr)-1 do begin
  
    restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/'+obs_name_arr[obs_i]+'_obs.sav'
    ;restore, '/nfs/eor-00/h1/nbarry/Aug23_cross_polyfit_quad_data/'+obs_name_arr[obs_i]+'pre_cal.sav'
    ;restore, '/nfs/eor-00/h1/nbarry/Aug23_twopolyquad_polyonly_constrainedphase/'+obs_name_arr[obs_i]+'_cal.sav'
    restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/calibration/'+obs_name_arr[obs_i]+'_cal.sav'
    ;restore, '/nfs/eor-00/h1/nbarry/Aug23_onequad_polyscaled_90150/cross_polyfit_data/'+obs_name_arr[obs_i]+'post_cal.sav'
    ;restore,'/nfs/eor-00/h1/nbarry/Aug23_polyratio_quadsplit/'+obs_name_arr[obs_i]+'_cal.sav'
    
    nt_use=N_Elements(where((*obs.baseline_info).tile_use EQ 1))
    gain_arr=complex(DBLARR(n_pol,n_freq,128))
    gain_arr[*,*,*]=1.
    
    tile_use=where((*obs.baseline_info).tile_use,nt_use_temp)
    
    cross_poly=complex(FLTARR(384,128,2))
    cross_poly[*,*,*]=1.
    
    ;cross_poly[256:383,tile_use,0]=((*cal_post.gain[0])[256:383,tile_use])
    ;cross_poly[0:255,tile_use,0]=((*cal_pre.gain[0])[0:255,tile_use])
    
    ;cross_poly[256:383,tile_use,1]=((*cal_post.gain[1])[256:383,tile_use])
    ;cross_poly[0:255,tile_use,1]=((*cal_pre.gain[1])[0:255,tile_use])
    
    cross_poly[*,tile_use,0]=((*cal.gain[0])[*,tile_use])
    cross_poly[*,tile_use,1]=((*cal.gain[1])[*,tile_use])
    
    if keyword_set(remove_poly) then begin
    
      gain_arr_fit=complex(DBLARR(n_pol,n_freq,nt_use))
      FOR pol_i=0, 1 DO BEGIN
        FOR tile_i=0,nt_use-1 DO BEGIN
        
          mode_params=DBLARR(3)
          If (cal.mode_params[pol_i,tile_use[tile_i]]) NE !NULL THEN mode_params=*(cal.mode_params[pol_i,tile_use[tile_i]])
          
          gain_arr_fit[pol_i,*,tile_i]=mode_params[1]*Exp(-Complex(0,1)*2.*!Pi*mode_params[0]*findgen(n_freq)/n_freq+Complex(0,1)*(mode_params[2]))
          
        ENDFOR
        
        cross_poly[*,tile_use,pol_i]=((*cal.gain[pol_i])[*,tile_use]-gain_arr_fit[pol_i,*,*])
      ENDFOR
      
    endif
    
    If keyword_set(seperate_phase) then begin
      mode_amp = getvar_savefile('/nfs/eor-00/h1/nbarry/Aug23_twopolyquad_x2fix_updatecompare_amponly/'+obs_name_arr[obs_i]+'_cal.sav','cal')
      phases = getvar_savefile('/nfs/eor-00/h1/nbarry/Aug23_twopolyquad_x2fix_updatecompare_phaseonly/'+obs_name_arr[obs_i]+'_cal.sav','cal')
      for i =0,127 do if mode_amp.mode_params[0,i] NE !NULL then mode_amp.mode_params[0,i]=PTR_NEW([(*mode_amp.mode_params[0,i])[0],(*mode_amp.mode_params[0,i])[1],(*phases.mode_params[0,i])[2]])
      for i =0,127 do if mode_amp.mode_params[1,i] NE !NULL then mode_amp.mode_params[1,i]=PTR_NEW([(*mode_amp.mode_params[1,i])[0],(*mode_amp.mode_params[1,i])[1],(*phases.mode_params[1,i])[2]])
      cal=mode_amp
    endif else restore,'/nfs/eor-00/h1/nbarry/Aug23_twopolyquad_autocheck/'+obs_name_arr[obs_i]+'_cal.sav'
    
    
    ;for ruby's plotting program
    ;GET_LUN, lun
    ;OPENU, lun, '~/mode_params_scaled_autos_even_more_flagging.txt'
    ;printf, lun,'tile ','pol ','mode ','amplitude ','phase '
    
    
    FOR pol_i=0, n_pol-1 DO BEGIN
      FOR tile_i=0,nt_use-1 DO BEGIN
      
        mode_params=DBLARR(3)
        
        If (cal.mode_params[pol_i,tile_use[tile_i]]) NE !NULL THEN mode_params=*(cal.mode_params[pol_i,tile_use[tile_i]])
        
        gain_arr[pol_i,*,tile_use[tile_i]]=cross_poly[*,tile_use[tile_i],pol_i]+mode_params[1]*Exp(-Complex(0,1)*2.*!Pi*mode_params[0]*findgen(n_freq)/n_freq+Complex(0,1)*(mode_params[2]))
        
      ;for ruby's plotting program
      ;mode_params_good_units=mode_params
      ;mode_params_good_units[0]=mode_params_good_units[0]/(3.072*10^7.)
      ;printf,lun,tile_use[tile_i],pol_i,mode_params_good_units
        
      ENDFOR
      (*cal.gain[pol_i])[*,*]=gain_arr[pol_i,*,*]
    ENDFOR
    
    ;free_lun,lun
    
    save, cal, Filename='/nfs/eor-00/h1/nbarry/Aug23_twopolyquad_autocheck/stdpoly/'+obs_name_arr[obs_i]+'_cal.sav'
  endfor
end