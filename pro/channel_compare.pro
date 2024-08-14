pro channel_compare, dir_name,day=day, advanced_plotting=advanced_plotting, longrun=longrun, split_jump=split_jump, $
    save_autos=save_autos,one_mode=one_mode,phase_transfer=phase_transfer,compare=compare,digjump_calc=digjump_calc
  ;set longrun to use obs used in the longrun
    
  ;List of dirnames I've been using
  ;dir_name='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_pointing_gainrescor_plots/fhd_nb_longrun_predigjump/fhd_nb_longrun_'+day    -- longrun dirname
  ;dir_name='/nfs/eor-00/h1/nbarry/Aug23_autos_modesparams_plusmodepointing_frombp_unsplit'   -- Aug23 autos that have polyfit/modefit by pointing, restored from recent run, no digital split
  ;dir_name='/nfs/eor-00/h1/nbarry/Aug23_pointing_nodigjump_v2_plusmodepointing_frombp'    --- Aug23 crosses that have polyfit/modefit by pointing, restored from recent run, digital split
  ;dir_name='/nfs/eor-00/h1/nbarry/Aug23_pointing_plusmodepointing_frombp'    --- Aug23 cross that have polyfit/modefit by pointing, resotred from recent run, no digital split
    
  ;**************setup
  edges_array=[1,14,17,30,33,46,49,62,65,78,81,94,97,110,113,126,129,142,145,158,161,174,177,190,193,206,209,222,225,238,241,254];,257,286,289,318,321,350,353,382]
  ;center_array=[15,15,47,47,79,79,111,111,143,143,175,175,207,207,239,239];,271,271,303,303,335,335,367,367]
  center_array=[7,7,23,23,39,39,55,55,71,71,87,87,103,103,119,119,135,135,151,151,167,167,183,183,199,199,215,215,231,231,247,247]
  
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
  
  ;***
  ;cgplot, [1,2,3],[1,2,3], psym=2,symsize=0.2, yrange=[4,7],xrange=[5,10],/NODATA
  ;cgplot, [1,2,3],[1,2,3], psym=2,symsize=0.2, yrange=[0,2],xrange=[0,2],/NODATA
  ;***
  
  
  ;*************Analysis by pointing
  FOR j=0, (size(pointing_num))[1]-1 DO BEGIN
    obsid=(*obs_ptr[j])[0]
    
    ;If there are obs in the current pointing, then do the restore. 'empty' was a filler tag placed above
    IF obsid NE 'empty' THEN BEGIN
    
      for ii=0,1 do begin
        if ii EQ 0 then save_autos=1
        if ii EQ 1 then undefine, save_autos
        ;****Restore loop, per obs in a pointing
        For i=0, parsednumbers[j]-1 DO BEGIN
          obsid=(*obs_ptr[j])[i]
          
          ;I think this makes sure the obs actually exists, as there might be obs in the text file that passed tests but didn't get run? I think this is superfluous.
          IF obsid NE 'empty' THEN BEGIN
          
            ;Restore obs save file, either from a regular run or a longrun (different obs potentially)
            If keyword_set(longrun) then restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/' + obsid + '_obs.sav' else $
              ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/metadata/' + obsid + '_obs.sav'
              restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/' + obsid + '_obs.sav'
              
            ;Setup the obs structure array and fill it on successive loops
            If (i eq 0) Then obs_array = replicate(obs,parsednumbers[j]) ELSE obs_array[i]=obs
            
            ;Restore params save file, either from regular run or a longrun
            If keyword_set(longrun) then restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/' + obsid + '_params.sav' else $
              ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/metadata/' + obsid + '_params.sav'
              restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/' + obsid + '_params.sav'
              
            ;Setup the params structure array and fill it on successive loops
            If (i eq 0) Then params_array = replicate(params,parsednumbers[j]) ELSE params_array[i]=params
            
            ;Restore cal structure from longrun or regular run
            If keyword_set(longrun) then restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/' + obsid + '_cal.sav' else $
              ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/calibration/' + obsid + '_cal.sav'
              restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/calibration/' + obsid + '_cal.sav'
              
            if keyword_set(phase_transfer) then begin
              phase_transfer=FLTARR(2,384,128)
              phase_transfer[0,*,*]=atan(*cal.gain[0],/phase)
              phase_transfer[1,*,*]=atan(*cal.gain[1],/phase)
            endif
            
            ;If the autos keyword is set, restore the cal structure from a regular run (Aug23 obsids). Only an option for nonlongruns
            If keyword_set(save_autos) then restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_autogainsonly_May2015/calibration/' + obsid + '_cal.sav'
            
            ;Make an array of structures.  First make a new cal structure and fill it with the value restored.  There are various points in the init code that set the below values to
            ;zero, so they have to be refilled individually.
            cal2=fhd_struct_init_cal(obs, params, n_pol=cal.n_pol,n_freq=cal.n_freq,n_tile=cal.n_tile,freq=cal.freq,gain=cal.gain,gain_residual=cal.gain_residual, $
              n_vis_cal=cal.n_vis_cal, amp_params=cal.amp_params, phase_params=cal.phase_params,mode_params=cal.mode_params)
            cal2.mode_params = cal.mode_params
            cal2.amp_params=cal.amp_params
            cal2.phase_params=cal.phase_params
            cal2.gain_residual=cal.gain_residual
            if keyword_set(save_autos) then cal2.auto_params=cal.auto_params
            
            ;Make and array of structures from the new cal structure made just previously
            If (i eq 0) Then cal_array = replicate(cal2,parsednumbers[j]) ELSE cal_array[i]=cal2
            
            If (i eq 0) Then n_pol = intarr(parsednumbers[j])
            n_pol[i]=cal_array[i].n_pol
            If ~keyword_set(save_autos) then begin
              for pol_i=0,n_pol[i]-1 do *cal_array[i].gain[pol_i]=*cal_array[i].gain[pol_i]+*cal_array[i].gain_residual[pol_i]
            endif else begin
              ;Note that the abs is there. shouldn't effect bp, but poly for sure since autos have no phase info.
              for freq_i=0,cal_array[i].n_freq-1 do (*cal_array[i].gain[0])[freq_i,*]=(abs((*cal_array[i].gain[0])[freq_i,*])-(*cal_array[i].auto_params[0])[0,*])/(*cal_array[i].auto_params[0])[1,*]
              for freq_i=0,cal_array[i].n_freq-1 do (*cal_array[i].gain[1])[freq_i,*]=(abs((*cal_array[i].gain[1])[freq_i,*])-(*cal_array[i].auto_params[1])[0,*])/(*cal_array[i].auto_params[1])[1,*]
              if keyword_set(phase_transfer) then begin
                *cal_array[i].gain[0]=*cal_array[i].gain[0]*exp(Complex(0,1)*phase_transfer[0,*,*])
                *cal_array[i].gain[1]=*cal_array[i].gain[1]*exp(Complex(0,1)*phase_transfer[1,*,*])
              endif
            endelse
            cal_array_orig=cal_array
            
            if i EQ 0 AND ii EQ 0 then cal_array_autos=cal_array
            if i EQ 0 AND ii EQ 1 then cal_array_cross=cal_array
            
            if ii EQ 0 then cal_array_autos[i].gain=Pointer_copy(cal_array[i].gain)
            if ii EQ 1 then cal_array_cross[i].gain=Pointer_copy(cal_array[i].gain)
          ;cal_bandpass=vis_cal_bandpass(cal_array[i],obs_array[i],cal_remainder=cal_remainder,saved_run_bp=1,cable_bandpass_fit=1)
          ;for pol_i=0,n_pol[i]-1 do *cal_array[i].gain[pol_i]=*cal_remainder.gain[pol_i]
            
            
          ENDIF
        ENDFOR
      ;****End of restore loop
        
        
      endfor
      
      mode_filepath=filepath(obs_array[0].instrument+'_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
      textfast,data_array,/read,file_path=mode_filepath,first_line=1
      cable_len=Reform(data_array[2,*])
      tile_use_cable=PTRARR(6,/allocate)
      *tile_use_cable[0]=where((*obs.baseline_info).tile_use AND cable_len EQ 90)
      *tile_use_cable[1]=where((*obs.baseline_info).tile_use AND cable_len EQ 150)
      *tile_use_cable[2]=where((*obs.baseline_info).tile_use AND cable_len EQ 230)
      *tile_use_cable[3]=where((*obs.baseline_info).tile_use AND cable_len EQ 320)
      *tile_use_cable[4]=where((*obs.baseline_info).tile_use AND cable_len EQ 400)
      *tile_use_cable[5]=where((*obs.baseline_info).tile_use AND cable_len EQ 524)
      
      
      freq_use=where((*obs.baseline_info).freq_use)
      ratio_array=FLTARR(2,parsednumbers[j],N_elements(edges_array),128)
      
      
      for cable_i=0, 5 do begin
        if cable_i EQ 0 then cablecolor='black'
        if cable_i EQ 1 then cablecolor='blue'
        if cable_i EQ 2 then cablecolor='green'
        if cable_i EQ 3 then cablecolor='red'
        if cable_i EQ 4 then cablecolor='purple'
        if cable_i EQ 5 then cablecolor='orange'
        
        For i=0, parsednumbers[j]-1 DO BEGIN
        ;For i=0, 0 DO BEGIN
          ratio_array[0,i,*,*]=abs((*cal_array_autos[i].gain[0])[edges_array,*])/abs((*cal_array_cross[i].gain[0])[edges_array,*])
          ratio_array[1,i,*,*]=abs((*cal_array_autos[i].gain[1])[edges_array,*])/abs((*cal_array_cross[i].gain[1])[edges_array,*])
          if i EQ 0 AND cable_i EQ 0 then cgplot, abs((*cal_array_cross[0].gain[0])[center_array,(*tile_use_cable[0])[0]]),ratio_array[0,0,*,(*tile_use_cable[0])[0]],psym=2,symsize=0.2, color=cablecolor, $
            xtitle='amp of center of channel, crosses', ytitle='amp of edge of channel, crosses',/NODATA
          ;if i EQ 0 AND cable_i EQ 0 then $
          ;  cgplot, abs((*cal_array_autos[0].gain[0])[center_array,(*tile_use_cable[0])[0]]),$
          ;    abs((*cal_array_autos[0].gain[0])[edges_array,(*tile_use_cable[0])[0]])-abs((*cal_array_autos[0].gain[0])[center_array,(*tile_use_cable[0])[0]]),$
          ;    psym=2,symsize=0.2, color=cablecolor, /NODATA

          for tile_i=0,N_elements(*tile_use_cable[cable_i])-1 do begin
          ;for tile_i=0,0 do begin
            cgoplot, abs((*cal_array_autos[i].gain[0])[center_array,(*tile_use_cable[cable_i])[tile_i]]),ratio_array[0,i,*,(*tile_use_cable[cable_i])[tile_i]],psym=2,symsize=0.2, color=cablecolor
          ;  cgoplot, abs((*cal_array_autos[i].gain[0])[center_array,(*tile_use_cable[cable_i])[tile_i]]),$
          ;    abs((*cal_array_autos[i].gain[0])[edges_array,(*tile_use_cable[cable_i])[tile_i]])-abs((*cal_array_autos[i].gain[0])[center_array,(*tile_use_cable[cable_i])[tile_i]]),$
          ;    psym=2,symsize=0.2, color=cablecolor

          endfor ; end tile for

        endfor ;end obs for
        stop
      endfor ;end cable for
      
      stop
    endif
  endfor
end