pro mkbp_pointings_autoinput, dir_name,day=day, advanced_plotting=advanced_plotting, longrun=longrun, split_jump=split_jump, save_autos=save_autos,one_mode=one_mode
  ;set longrun to use obs used in the longrun

  ;List of dirnames I've been using
  ;dir_name='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_pointing_gainrescor_plots/fhd_nb_longrun_predigjump/fhd_nb_longrun_'+day    -- longrun dirname
  ;dir_name='/nfs/eor-00/h1/nbarry/Aug23_autos_modesparams_plusmodepointing_frombp_unsplit'   -- Aug23 autos that have polyfit/modefit by pointing, restored from recent run, no digital split
  ;dir_name='/nfs/eor-00/h1/nbarry/Aug23_pointing_nodigjump_v2_plusmodepointing_frombp'    --- Aug23 crosses that have polyfit/modefit by pointing, restored from recent run, digital split
  ;dir_name='/nfs/eor-00/h1/nbarry/Aug23_pointing_plusmodepointing_frombp'    --- Aug23 cross that have polyfit/modefit by pointing, resotred from recent run, no digital split

  ;**************setup
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
    *obs_ptr[j]=obs_temp
    
  ENDFOR
  ;*************end of setup
  
  
  ;*************Analysis by pointing
  FOR j=0, (size(pointing_num))[1]-1 DO BEGIN
    obsid=(*obs_ptr[j])[0]
    
    ;If there are obs in the current pointing, then do the restore. 'empty' was a filler tag placed above
    IF obsid NE 'empty' THEN BEGIN
    
      ;****Restore loop, per obs in a pointing
      For i=0, parsednumbers[j]-1 DO BEGIN
        obsid=(*obs_ptr[j])[i]
        obsid_count=parsednumbers[j]
        
        ;I think this makes sure the obs actually exists, as there might be obs in the text file that passed tests but didn't get run? I think this is superfluous.
        IF obsid NE 'empty' THEN BEGIN
        
          ;Restore obs save file, either from a regular run or a longrun (different obs potentially)
          If keyword_set(longrun) then restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/' + obsid + '_obs.sav' else $
            ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/metadata/' + obsid + '_obs.sav'
            restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_pointing_May2015/metadata/' + obsid + '_obs.sav'
            
          ;Setup the obs structure array and fill it on successive loops
          If (i eq 0) Then obs_array = replicate(obs,parsednumbers[j]) ELSE obs_array[i]=obs
          
          ;Restore params save file, either from regular run or a longrun
          If keyword_set(longrun) then restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/' + obsid + '_params.sav' else $
            ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/metadata/' + obsid + '_params.sav'
            restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_pointing_May2015/metadata/' + obsid + '_params.sav'
            
          ;Setup the params structure array and fill it on successive loops
          If (i eq 0) Then params_array = replicate(params,parsednumbers[j]) ELSE params_array[i]=params
          
          ;Restore cal structure from longrun or regular run
          If keyword_set(longrun) then restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/' + obsid + '_cal.sav' else $
            ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/calibration/' + obsid + '_cal.sav'
            restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_pointing_May2015/calibration/' + obsid + '_cal.sav'
            
          ;If the autos keyword is set, restore the cal structure from a regular run (Aug23 obsids). Only an option for nonlongruns
          ;If keyword_set(save_autos) then restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_autogainsonly_May2015/calibration/' + obsid + '_cal.sav'
            
          ;Make an array of structures.  First make a new cal structure and fill it with the value restored.  There are various points in the init code that set the below values to
          ;zero, so they have to be refilled individually.
          cal2=fhd_struct_init_cal(obs, params, n_pol=cal.n_pol,n_freq=cal.n_freq,n_tile=cal.n_tile,freq=cal.freq,gain=cal.gain,gain_residual=cal.gain_residual, $
            n_vis_cal=cal.n_vis_cal, amp_params=cal.amp_params, phase_params=cal.phase_params,mode_params=cal.mode_params)
          cal2.mode_params = cal.mode_params
          cal2.amp_params=cal.amp_params
          cal2.phase_params=cal.phase_params
          cal2.gain_residual=cal.gain_residual
          
          ;Make and array of structures from the new cal structure made just previously
          If (i eq 0) Then cal_array = replicate(cal2,parsednumbers[j]) ELSE cal_array[i]=cal2
          
          If (i eq 0) Then n_pol = intarr(parsednumbers[j])
          n_pol[i]=cal_array[i].n_pol
          for pol_i=0,n_pol[i]-1 do *cal_array[i].gain[pol_i]=*cal_array[i].gain[pol_i]+*cal_array[i].gain_residual[pol_i]
          
          
          ;If the autos keyword is set, restore the cal structure from a regular run (Aug23 obsids). Only an option for nonlongruns
          restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_autogainsonly_May2015/calibration/' + obsid + '_cal.sav'
          
          cal2=fhd_struct_init_cal(obs, params, n_pol=cal.n_pol,n_freq=cal.n_freq,n_tile=cal.n_tile,freq=cal.freq,gain=cal.gain,gain_residual=cal.gain_residual, $
            n_vis_cal=cal.n_vis_cal, amp_params=cal.amp_params, phase_params=cal.phase_params,mode_params=cal.mode_params)
          cal2.mode_params = cal.mode_params
          cal2.amp_params=cal.amp_params
          cal2.phase_params=cal.phase_params
          cal2.gain_residual=cal.gain_residual
          cal2.auto_params=cal.auto_params
          
          ;Make and array of structures from the new cal structure made just previously
          If (i eq 0) Then begin
            cal_array_autos = cal2
            cal_array_scaled=cal2
            cal_array_unscaled=cal2
            cal_array_unscaled.gain=POINTER_COPY(cal2.gain)
            cal_array_scaled.gain=POINTER_COPY(cal2.gain)
            cal_array_unscaled.mode_params=POINTER_COPY(cal2.mode_params)
            cal_array_scaled.mode_params=POINTER_COPY(cal2.mode_params)
            
          endif ELSE begin
            cal_array_autos=[cal_array_autos,cal2]
            cal_array_scaled=[cal_array_scaled,cal2]
            cal_array_unscaled=[cal_array_unscaled,cal2]
            
            cal_array_unscaled[i].gain=POINTER_COPY(cal2.gain)
            cal_array_scaled[i].gain=POINTER_COPY(cal2.gain)
            cal_array_unscaled[i].mode_params=POINTER_COPY(cal2.mode_params)
            cal_array_scaled[i].mode_params=POINTER_COPY(cal2.mode_params)
          endelse
          
          
          ;Note that the abs is there. shouldn't effect bp, but poly for sure since autos have no phase info.
          ;unscaled_temp_xx=complex(FLTARR(384,128))
          ;unscaled_temp_yy=complex(FLTARR(384,128))
          for freq_i=0,cal_array_autos[i].n_freq-1 do begin
            (*cal_array_unscaled[i].gain[0])[freq_i,*]=(abs((*cal_array_autos[i].gain[0])[freq_i,*])-(*cal_array_autos[i].auto_params[0])[0,*])/(*cal_array_autos[i].auto_params[0])[1,*]
            ; unscaled_temp_xx[freq_i,*]=(abs((*cal_array_autos[i].gain[0])[freq_i,*])-(*cal_array_autos[i].auto_params[0])[0,*])/(*cal_array_autos[i].auto_params[0])[1,*]
            (*cal_array_unscaled[i].gain[1])[freq_i,*]=(abs((*cal_array_autos[i].gain[1])[freq_i,*])-(*cal_array_autos[i].auto_params[1])[0,*])/(*cal_array_autos[i].auto_params[1])[1,*]
          ;unscaled_temp_yy[freq_i,*]=(abs((*cal_array_autos[i].gain[1])[freq_i,*])-(*cal_array_autos[i].auto_params[1])[0,*])/(*cal_array_autos[i].auto_params[1])[1,*]
          endfor
          
          cal_bandpass=vis_cal_bandpass(cal_array[i],obs_array[i],cal_remainder=cal_remainder,saved_run_bp=1,cable_bandpass_fit=1)
          for pol_i=0,n_pol[i]-1 do *cal_array[i].gain[pol_i]=*cal_remainder.gain[pol_i]
          
          cal_bandpass=vis_cal_bandpass(cal_array_scaled[i],obs_array[i],cal_remainder=cal_remainder,saved_run_bp=1,cable_bandpass_fit=1)
          for pol_i=0,n_pol[i]-1 do *cal_array_scaled[i].gain[pol_i]=*cal_remainder.gain[pol_i]
          
          cal_bandpass=vis_cal_bandpass(cal_array_unscaled[i],obs_array[i],cal_remainder=cal_remainder,saved_run_bp=1,cable_bandpass_fit=1)
          for pol_i=0,n_pol[i]-1 do *cal_array_unscaled[i].gain[pol_i]=*cal_remainder.gain[pol_i]

        ENDIF
        
      ENDFOR
      ;****End of restore loop
      
      file_mkdir, dir_name, /NOEXPAND_PATH
      ;bp=vis_cal_bandpass_with_cable_pointingv2(cal_array,obs_array,parsednumbers[j],pointing_num[j],cal_remainder=cal_remainder,file_path_fhd=dir_name+'/',advanced_plotting=advance_plotting)
      
      ;bp=vis_cal_polyfit_pointing(cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=2,phase_degree=1,$
      ;  file_path=dir_name+'/')
      
      If keyword_set(split_jump) then begin
        ;Only a bandpass within the specified frequencies are calculated, so a only pre dig jump looks truncated
      
      
        ;*************************Split pointing calc block
        ;Fit parameters are only calculated in the specified frequencies, but the final fit is over all frequencies.
        ;Pre dig jump calcs should be fine, even though the product looks like it's done over all freq
        ;Just fitting linear trend, bc don't want to overfit and it will be interesting to see if degree=2 was just trying to capture dig jump
      
        cal_array_cross=cal_array
        
        for iter_i=0,2 do begin
        
          if iter_i EQ 1 then cal_array=cal_array_unscaled
          if iter_i EQ 2 then cal_array=cal_array_scaled
          
          original_freq_use=(*obs_array[0].baseline_info).freq_use
          ;flag frequencies after the digital gain jump
          for obs_i=0, parsednumbers[j]-1 do (*obs_array[obs_i].baseline_info).freq_use[256:383]=0
          ;Calc the pre dig jump poly+mode fit over the pointing to try to avoid sidelobe effects
          cal_polyfit_pre=vis_cal_polyfit_pointing(cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=1,phase_degree=1,$
            file_path=dir_name+'/',cal_cable_reflection_mode_fit=1,cal_cable_reflection_fit=cal_cable_reflection_fit,cal_mode_fit=cal_mode_fit)
          ;Restore the original freq_use to erase frequency flagging from before
          for obs_i=0, parsednumbers[j]-1 do (*obs_array[obs_i].baseline_info).freq_use=original_freq_use
          ;flag the frequencies before the digital gain jump
          for obs_i=0, parsednumbers[j]-1 do (*obs_array[obs_i].baseline_info).freq_use[0:255]=0
          ;Calc the pre dig jump poly+mode fit over the pointing to try to avoid sidelobe effects
          cal_polyfit_post=vis_cal_polyfit_pointing(cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=1,phase_degree=1,$
            file_path=dir_name+'/',cal_cable_reflection_mode_fit=1,cal_cable_reflection_fit=cal_cable_reflection_fit,cal_mode_fit=cal_mode_fit)
          ;Restore the original freq_use to erase frequency flaggin
          for obs_i=0, parsednumbers[j]-1 do (*obs_array[obs_i].baseline_info).freq_use=original_freq_use
          
          
          if iter_i EQ 0 then begin
            cal_polyfit_post_cross=cal_polyfit_post
            cal_polyfit_pre_cross=cal_polyfit_pre
          endif
          if iter_i EQ 1 then begin
            cal_polyfit_post_auto_unscaled=cal_polyfit_post
            cal_polyfit_pre_auto_unscaled=cal_polyfit_pre
          endif
          if iter_i EQ 2 then begin
            cal_polyfit_post_auto_scaled=cal_polyfit_post
            cal_polyfit_pre_auto_scaled=cal_polyfit_pre
          endif
          
        endfor
        
        cal_array=cal_array_cross
        cal_polyfit_post=cal_polyfit_post_cross
        cal_polyfit_pre=cal_polyfit_pre_cross
        
        
        poly_input_gains=complex(FLTARR(384,128,2,parsednumbers[j]))
        poly_input_gains[*,*,*]=1.
        
        for obs_i=0, parsednumbers[j]-1 do begin
        
          tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
          
          poly_input_gains[256:383,tile_use,0,obs_i]=((*cal_polyfit_post[obs_i].gain[0])[256:383,tile_use])
          poly_input_gains[0:255,tile_use,0,obs_i]=((*cal_polyfit_pre[obs_i].gain[0])[0:255,tile_use])
          
          poly_input_gains[256:383,tile_use,1,obs_i]=((*cal_polyfit_post[obs_i].gain[1])[256:383,tile_use])
          poly_input_gains[0:255,tile_use,1,obs_i]=((*cal_polyfit_pre[obs_i].gain[1])[0:255,tile_use])
          
        endfor
        
        poly_input_gains_unscaled=complex(FLTARR(384,128,2,parsednumbers[j]))
        poly_input_gains_unscaled[*,*,*]=1.
        for obs_i=0, parsednumbers[j]-1 do begin
        
          tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
          
          poly_input_gains_unscaled[256:383,tile_use,0,obs_i]=((*cal_polyfit_post_auto_unscaled[obs_i].gain[0])[256:383,tile_use])
          poly_input_gains_unscaled[0:255,tile_use,0,obs_i]=((*cal_polyfit_pre_auto_unscaled[obs_i].gain[0])[0:255,tile_use])
          
          poly_input_gains_unscaled[256:383,tile_use,1,obs_i]=((*cal_polyfit_post_auto_unscaled[obs_i].gain[1])[256:383,tile_use])
          poly_input_gains_unscaled[0:255,tile_use,1,obs_i]=((*cal_polyfit_pre_auto_unscaled[obs_i].gain[1])[0:255,tile_use])
          
        endfor
        poly_input_gains_scaled=complex(FLTARR(384,128,2,parsednumbers[j]))
        poly_input_gains_scaled[*,*,*]=1.
        for obs_i=0, parsednumbers[j]-1 do begin
        
          tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
          
          poly_input_gains_scaled[256:383,tile_use,0,obs_i]=((*cal_polyfit_post_auto_scaled[obs_i].gain[0])[256:383,tile_use])
          poly_input_gains_scaled[0:255,tile_use,0,obs_i]=((*cal_polyfit_pre_auto_scaled[obs_i].gain[0])[0:255,tile_use])
          
          poly_input_gains_scaled[256:383,tile_use,1,obs_i]=((*cal_polyfit_post_auto_scaled[obs_i].gain[1])[256:383,tile_use])
          poly_input_gains_scaled[0:255,tile_use,1,obs_i]=((*cal_polyfit_pre_auto_scaled[obs_i].gain[1])[0:255,tile_use])
          
        endfor
        
        mode_input=PTRARR(obsid_count,2,128)
        
        cal_polyfit_total=vis_cal_modefit_only_pointing(poly_input_gains_unscaled,cal_array_unscaled,obs_array,parsednumbers[j],pointing_num[j],degree=1,phase_degree=1,$
          file_path=dir_name+'/',cal_cable_reflection_mode_fit=1)
        ;         for pol_i, n_pol-1 do begin
        ;          for tile_i,
        for obs_j=0,obsid_count-1 do mode_input[obs_j,*,*]=((cal_polyfit_total[obs_j].mode_params))
        ;     endfor

        cal_polyfit_total=vis_cal_modefit_only_pointing(poly_input_gains_scaled,cal_array_scaled,obs_array,parsednumbers[j],pointing_num[j],degree=1,phase_degree=1,$
          file_path=dir_name+'/',cal_cable_reflection_mode_fit=1, mode_input=mode_input)

        for pol_i=0,1 do begin
          for obs_i=0,obsid_count-1 do begin
            gain_arr=*cal_array[obs_i].gain[pol_i]
            gain_arr_fit=poly_input_gains[*,*,pol_i,obs_i]
            
            IF N_Elements(obs_array[obs_i]) GT 0 THEN tile_use_temp=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
            
            FOR ti=0L,nt_use_temp-1 DO BEGIN
              tile_i=tile_use_temp[ti]
              
              IF mode_input[obs_i,pol_i,tile_i] NE !NULL then mode_use=(*mode_input[obs_i,pol_i,tile_i])[0] else mode_use=0
              If mode_use EQ 0 then continue
              
              amp_use=(*cal_polyfit_total[obs_i].mode_params[pol_i,tile_i])[1]
              phase_use=(*cal_polyfit_total[obs_i].mode_params[pol_i,tile_i])[2]
              
              gain_mode_fit=amp_use*exp(-Complex(0,1)*2.*!Pi*(mode_use*findgen(384)/384.)+Complex(0,1)*phase_use)
              gain_arr_fit[*,tile_i]+=gain_mode_fit
              cal_array_cross[obs_i].mode_params[pol_i,tile_i]=Ptr_new([mode_use,amp_use,phase_use])
              
            endfor
            *cal_array_cross[obs_i].gain[pol_i]=gain_arr_fit
            
          endfor
        endfor
        ;*************************Split pointing calc block
        
        
        ;for saving the cal structures individually for restoring, especially with adams fourth line plots
        
        For i=0, parsednumbers[j]-1 DO BEGIN
          obsid=(*obs_ptr[j])[i]
          cal=cal_array_cross[i]
          save, cal, FILENAME=dir_name+'/'+obsid+'_cal.sav
        endfor
        
      endif else begin
      
      
        ;*************************Unsplit pointing calc block
        cal_polyfit=vis_cal_polyfit_pointing(cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=2,phase_degree=1,$
          file_path=dir_name+'/',cal_cable_reflection_mode_fit=1)
          
        For i=0, parsednumbers[j]-1 DO BEGIN
          obsid=(*obs_ptr[j])[i]
          cal=cal_polyfit[i]
          save, cal, FILENAME=dir_name+'/'+obsid+'_cal.sav
        endfor
      ;*************************end of unsplit pointing calc block
        
      endelse
      
      
    ENDIF
  ENDFOR
;*************end of analysis by pointing
  
END