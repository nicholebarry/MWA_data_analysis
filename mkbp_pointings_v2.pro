pro mkbp_pointings_v2, dir_name,day=day, advanced_plotting=advanced_plotting, longrun=longrun, split_jump=split_jump, $
    save_autos=save_autos,one_mode=one_mode,phase_transfer=phase_transfer,compare=compare,digjump_calc=digjump_calc, $
    split_compare=split_compare, mode_input=mode_input
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
  ;If keyword_set(longrun) then parsednames=[day+'_minusfour',day+'_minusthree',day+'_minustwo',day+'_minusone',day+'_zenith',day+'_plusone',day+'_plustwo',day+'_plusthree',day+'_plusfour'] else $
    parsednames=[day+'minustwo',day+'minusone',day+'zenith',day+'plusone',day+'plustwo',day+'plusthree']
  ;  parsednames=[day+'_minustwo',day+'_minusone',day+'_zenith',day+'_plusone',day+'_plustwo',day+'_plusthree']
    
  ;obs ptr is where obs of a pointing are stored. It can be made bigger than needed to account for the pointings of the longrun
  obs_ptr=PTRARR(10,/allocate)
  
  ;Different pointing ranges given if it was the longrun or not
  ;If keyword_set(longrun) then pointing_num=[-4,-3,-2,-1,0,1,2,3,4] else 
  pointing_num=[-2,-1,0,1,2,3]
  
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
            ;Crosses are saved such that solutions must be added to residuals to get back the original unfit gain
            for pol_i=0,n_pol[i]-1 do *cal_array[i].gain[pol_i]=*cal_array[i].gain[pol_i]+*cal_array[i].gain_residual[pol_i]
          endif else begin
          
          ;This is a reconstruction of the auto original gain, since window power scaled autos have had problems in the past
          ;Note that the abs is there. shouldn't effect bp, but poly for sure since autos have no phase info.
            for freq_i=0,cal_array[i].n_freq-1 do (*cal_array[i].gain[0])[freq_i,*]=(abs((*cal_array[i].gain[0])[freq_i,*])-(*cal_array[i].auto_params[0])[0,*])/(*cal_array[i].auto_params[0])[1,*]
            for freq_i=0,cal_array[i].n_freq-1 do (*cal_array[i].gain[1])[freq_i,*]=(abs((*cal_array[i].gain[1])[freq_i,*])-(*cal_array[i].auto_params[1])[0,*])/(*cal_array[i].auto_params[1])[1,*]
          
          ;Adding the cross phases to the amplitude of the auto gains to give a complex gain
            if keyword_set(phase_transfer) then begin
              *cal_array[i].gain[0]=*cal_array[i].gain[0]*exp(Complex(0,1)*phase_transfer[0,*,*])
              *cal_array[i].gain[1]=*cal_array[i].gain[1]*exp(Complex(0,1)*phase_transfer[1,*,*])
            endif
          
          endelse
          cal_array_orig=cal_array
          
          ;Do a simple bandpass (saved run std by pointing) to prep the unfit gains to go through various polyfits/modefits
          cal_bandpass=vis_cal_bandpass(cal_array[i],obs_array[i],cal_remainder=cal_remainder,saved_run_bp=1,cable_bandpass_fit=1)
          for pol_i=0,n_pol[i]-1 do *cal_array[i].gain[pol_i]=*cal_remainder.gain[pol_i]
          
        ENDIF
      ENDFOR
      ;*************End of restore loop
      
      
      
      
      file_mkdir, dir_name, /NOEXPAND_PATH
      
      if keyword_set(one_mode) then pointingmode=one_mode
      
      ;*******************Split the analysis about the digital gain jump
      If keyword_set(split_jump) then begin
      
        ;*********Read in and setup of a polyfit that will be used in a scaling ratio later
        if keyword_set(compare) then begin
        
          compare_boo=file_test(compare,/Directory)
          If compare_boo EQ 0 then begin
            print, "Compare keyword must be a directory to the data to compare with"
            stop
          endif
          
          cross_poly=complex(FLTARR(384,128,2,parsednumbers[j]))
          cross_poly[*,*,*]=1.
          
          for obs_i=0, parsednumbers[j]-1 do begin
          
            ;CAN CHANGE!!!!!!
            ;restore,'/nfs/eor-00/h1/nbarry/Aug23_onequad_polyscaled_90150/cross_polyfit_data/'+obs_array[obs_i].obsname+'pre_cal.sav'
            ;restore,'/nfs/eor-00/h1/nbarry/Aug23_onequad_polyscaled_90150/cross_polyfit_data/'+obs_array[obs_i].obsname+'post_cal.sav'
          
            if keyword_set(split_compare) then begin
              restore,compare+obs_array[obs_i].obsname+'pre_cal.sav'
              restore,compare+obs_array[obs_i].obsname+'post_cal.sav'
              
              tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
              
              cross_poly[256:383,tile_use,0,obs_i]=((*cal_post.gain[0])[256:383,tile_use])
              cross_poly[0:255,tile_use,0,obs_i]=((*cal_pre.gain[0])[0:255,tile_use])
              cross_poly[256:383,tile_use,1,obs_i]=((*cal_post.gain[1])[256:383,tile_use])
              cross_poly[0:255,tile_use,1,obs_i]=((*cal_pre.gain[1])[0:255,tile_use])
            endif else begin
              restore,compare+obs_array[obs_i].obsname+'_cal.sav'
              
              tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
              
              cross_poly[*,tile_use,0,obs_i]=((*cal.gain[0])[*,tile_use])
              cross_poly[*,tile_use,1,obs_i]=((*cal.gain[1])[*,tile_use])
              
            endelse
            
          endfor
          
        endif
        ;*********End of read in and setup of a polyfit that will be used in a scaling ratio later
        
        
        
        ;*************Split pointing calc block
        ;Fit parameters are only calculated in the specified frequencies, but the final fit is over all frequencies.
        ;Pre dig jump calcs should be fine, even though the product looks like it's done over all freq
        ;Just fitting linear trend, bc don't want to overfit and it will be interesting to see if degree=2 was just trying to capture dig jump
        cal_mode_fit=1
        cal_cable_reflection_fit=150
        degree=2
        undefine, cal_mode_fit, cal_cable_reflection_fit
        
        ;undefine, cal_mode_fit, cal_cable_reflection_fit
        
        original_freq_use=(*obs_array[0].baseline_info).freq_use
        ;flag frequencies after the digital gain jump
        for obs_i=0, parsednumbers[j]-1 do (*obs_array[obs_i].baseline_info).freq_use[256:383]=0
        ;Calc the pre dig jump poly+mode fit over the pointing to try to avoid sidelobe effects
        cal_polyfit_pre=vis_cal_polyfit_pointing(cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=degree,phase_degree=1,$
          file_path=dir_name+'/',cal_cable_reflection_mode_fit=1,cal_cable_reflection_fit=cal_cable_reflection_fit,cal_mode_fit=cal_mode_fit,pointingmode=pointingmode)
        ;Restore the original freq_use to erase frequency flagging from before
        for obs_i=0, parsednumbers[j]-1 do (*obs_array[obs_i].baseline_info).freq_use=original_freq_use
        ;flag the frequencies before the digital gain jump
        for obs_i=0, parsednumbers[j]-1 do (*obs_array[obs_i].baseline_info).freq_use[0:255]=0
        ;Calc the pre dig jump poly+mode fit over the pointing to try to avoid sidelobe effects
        cal_polyfit_post=vis_cal_polyfit_pointing(cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=degree,phase_degree=1,$
          file_path=dir_name+'/',cal_cable_reflection_mode_fit=1,cal_cable_reflection_fit=cal_cable_reflection_fit,cal_mode_fit=cal_mode_fit,pointingmode=pointingmode)
        ;Restore the original freq_use to erase frequency flaggin
        for obs_i=0, parsednumbers[j]-1 do (*obs_array[obs_i].baseline_info).freq_use=original_freq_use
        
        if ~keyword_set(pointingmode) then begin
          poly_input_gains=complex(FLTARR(384,128,2,parsednumbers[j]))
          poly_input_gains[*,*,*]=1.
          
          for obs_i=0, parsednumbers[j]-1 do begin
          
            tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
            
            poly_input_gains[256:383,tile_use,0,obs_i]=((*cal_polyfit_post[obs_i].gain[0])[256:383,tile_use])
            poly_input_gains[0:255,tile_use,0,obs_i]=((*cal_polyfit_pre[obs_i].gain[0])[0:255,tile_use])
            
            poly_input_gains[256:383,tile_use,1,obs_i]=((*cal_polyfit_post[obs_i].gain[1])[256:383,tile_use])
            poly_input_gains[0:255,tile_use,1,obs_i]=((*cal_polyfit_pre[obs_i].gain[1])[0:255,tile_use])
            
          endfor
          
        ;cal_polyfit_total=vis_cal_modefit_only_pointing(poly_input_gains,cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=1,phase_degree=1,$
        ;  file_path=dir_name+'/',cal_cable_reflection_mode_fit=1,cal_cable_reflection_fit=150)
          
        endif
        
        freq_use=where((*obs_array[0].baseline_info).freq_use)
        freq_dig_index=where(freq_use GT 256)
        freq_dig_num=N_elements(freq_use[freq_dig_index])
        
        
        if keyword_set(digjump_calc) then begin
        
          digjump=(FLTARR(parsednumbers[j],2,128))
          cal_nodig=cal_array
          for obs_i=0,parsednumbers[j]-1 do begin
            cal_nodig[obs_i].gain=Pointer_copy(cal_array[obs_i].gain)
            for pol_i=0,1 do begin
              ;The correct difference is between the two nonflagged channels surrounding the digital gain jump (255,256)
              digjump[obs_i,pol_i,*]=reform(abs((*cal_polyfit_pre[obs_i].gain[pol_i])[254,*])-abs((*cal_polyfit_post[obs_i].gain[pol_i])[257,*]))
              rebinned_digjump=FLTARR(freq_dig_num,128) ;freq by tiles
              for freqbin_i=0,freq_dig_num-1 do rebinned_digjump[freqbin_i,*]=digjump[obs_i,pol_i,*]
              (*cal_nodig[obs_i].gain[pol_i])[freq_use[freq_dig_index],*]=abs((*cal_array[obs_i].gain[pol_i])[freq_use[freq_dig_index],*])+rebinned_digjump
              
            endfor
          endfor
          
          ;This can only capture the amplitude fit since the digital gain jump is an amplitude correction
          cal_polyfit_nodigjump=vis_cal_polyfit_pointing(cal_nodig,obs_array,parsednumbers[j],pointing_num[j],degree=2,phase_degree=0,$
            file_path=dir_name+'/',cal_cable_reflection_mode_fit=1,cal_cable_reflection_fit=cal_cable_reflection_fit,cal_mode_fit=cal_mode_fit,pointingmode=pointingmode)
            
          for obs_i=0,parsednumbers[j]-1 do begin
            for pol_i=0,1 do begin
              rebinned_digjump=FLTARR(127,128)
              for freqbin_i=0,126 do rebinned_digjump[freqbin_i,*]=digjump[obs_i,pol_i,*]
              (*cal_polyfit_nodigjump[obs_i].gain[pol_i])[257:383,*]=(*cal_polyfit_nodigjump[obs_i].gain[pol_i])[257:383,*]-rebinned_digjump
            endfor
          endfor
        endif
        
        ;if keyword_set(one_mode) then begin
        poly_input_gains=complex(FLTARR(384,128,2,parsednumbers[j]))
        poly_input_gains[*,*,*]=1.
        
        for obs_i=0, parsednumbers[j]-1 do begin
        
          tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
          
          if ~keyword_set(digjump_calc) then begin
            poly_input_gains[256:383,tile_use,0,obs_i]=((*cal_polyfit_post[obs_i].gain[0])[256:383,tile_use])
            poly_input_gains[0:255,tile_use,0,obs_i]=((*cal_polyfit_pre[obs_i].gain[0])[0:255,tile_use])
            
            poly_input_gains[256:383,tile_use,1,obs_i]=((*cal_polyfit_post[obs_i].gain[1])[256:383,tile_use])
            poly_input_gains[0:255,tile_use,1,obs_i]=((*cal_polyfit_pre[obs_i].gain[1])[0:255,tile_use])
          endif else begin
          
            ;In order to capture the phase polyfit for the digital gain jump split option
            cal_polyfit_phase=vis_cal_polyfit_pointing(cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=2,phase_degree=1,$
              file_path=dir_name+'/',cal_cable_reflection_mode_fit=1,cal_cable_reflection_fit=cal_cable_reflection_fit,cal_mode_fit=cal_mode_fit,pointingmode=pointingmode)
              
            phase_fit=fltarr(384,2,128)
            for pol_i=0,1 do begin
              for tile_i=0,127 do begin
                phase_params=(*cal_polyfit_phase[0].phase_params[pol_i,tile_i])
                FOR di=0L,1 DO phase_fit[*,pol_i,tile_i]+=phase_params[di]*findgen(384)^di
              endfor
            endfor
            
            ;gain_arr[*,tile_i]=gain_fit*Exp(i_comp*phase_fit)
            
            poly_input_gains[*,tile_use,0,obs_i]=((*cal_polyfit_nodigjump[obs_i].gain[0])[*,tile_use])*Exp(Complex(0,1)*phase_fit[*,0,tile_use])
            poly_input_gains[*,tile_use,1,obs_i]=((*cal_polyfit_nodigjump[obs_i].gain[1])[*,tile_use])*Exp(Complex(0,1)*phase_fit[*,1,tile_use])
          endelse
          
          
        endfor
        
        if keyword_set(compare) then begin
          compare_ratio_array=abs(poly_input_gains)/abs(cross_poly)
          compare_ratio=FLTARR(128,2)
          for pol_j=0,1 do begin
            for tile_i=0,127 do begin
              resistant_mean,compare_ratio_array[*,tile_i,pol_j,*],2,res_mean
              compare_ratio[tile_i,pol_j]=res_mean
            endfor
          endfor
        endif
        
        
        if keyword_set(mode_input) then begin
        
          mode_input_arr=PTRARR(parsednumbers[j],2,128,/allocate)
          
          for obs_i=0, parsednumbers[j]-1 do begin
          
            restore,mode_input+obs_array[obs_i].obsname+'_cal.sav'
            
            tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
            
            for tile_i=0, N_elements(tile_use)-1 do begin
            
              mode_input_arr[obs_i,0,tile_use[tile_i]]=(cal.mode_params[0,tile_use[tile_i]])
              mode_input_arr[obs_i,1,tile_use[tile_i]]=(cal.mode_params[1,tile_use[tile_i]])
              
            endfor
            
          endfor
          
        endif
        
        cal_polyfit_total=cal_array
        for obs_j=0,parsednumbers[j]-1 do begin
          for pol_j=0,1 do begin
            (cal_polyfit_total[obs_j].gain[pol_j])=Pointer_copy(cal_array[obs_j].gain[pol_j])
            (*cal_polyfit_total[obs_j].gain[pol_j])[*,*]=poly_input_gains[*,*,pol_j,obs_j]
          endfor
        endfor
        cal_polyfit_total=vis_cal_modefit_only_pointing(poly_input_gains,cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=1,phase_degree=1,$
          file_path=dir_name+'/',cal_cable_reflection_mode_fit=1,compare_ratio=compare_ratio,cal_cable_reflection_both=1,pointingmode=pointingmode,mode_input=mode_input_arr)
        ;cal_polyfit_total=vis_cal_modefit_only_pointing(poly_input_gains,cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=1,phase_degree=1,$
        ;  file_path=dir_name+'/',cal_cable_reflection_mode_fit=1,cal_cable_reflection_fit=150, compare_ratio=compare_ratio, pointingmode=pointingmode)
        ;endif
        ;*************plit pointing calc block
        
        
        ;for saving the cal structures individually for restoring, especially with adams fourth line plots
        ;if ~keyword_set(one_mode) then begin
        ;  For i=0, parsednumbers[j]-1 DO BEGIN
        ;    obsid=(*obs_ptr[j])[i]
        ;    cal_pre=cal_polyfit_pre[i]
        ;    save, cal_pre, FILENAME=dir_name+'/'+obsid+'pre_cal.sav
        ;    cal_post=cal_polyfit_post[i]
        ;    save, cal_post, FILENAME=dir_name+'/'+obsid+'post_cal.sav
        ;  endfor
        ;endif else begin
        For i=0, parsednumbers[j]-1 DO BEGIN
          obsid=(*obs_ptr[j])[i]
          cal=cal_polyfit_total[i]
          save, cal, FILENAME=dir_name+'/'+obsid+'_cal.sav
        endfor
      ;endelse
        
      endif else begin
        ;*******************End of split analysis across digital gain jump
      
      
      
        if keyword_set(compare) then begin
        
          cross_poly=complex(FLTARR(384,128,2,parsednumbers[j]))
          cross_poly[*,*,*,*]=1.
          
          for obs_i=0, parsednumbers[j]-1 do begin
          
            ;CAN CHANGE!!!!!!
            restore,'/nfs/eor-00/h1/nbarry/Aug23_onequad_polyscaled_90150/cross_polyfit_data/'+obs_array[obs_i].obsname+'pre_cal.sav'
            restore,'/nfs/eor-00/h1/nbarry/Aug23_onequad_polyscaled_90150/cross_polyfit_data/'+obs_array[obs_i].obsname+'post_cal.sav'
            
            tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
            
            cross_poly[256:383,tile_use,0,obs_i]=((*cal_post.gain[0])[256:383,tile_use])
            cross_poly[0:255,tile_use,0,obs_i]=((*cal_pre.gain[0])[0:255,tile_use])
            
            cross_poly[256:383,tile_use,1,obs_i]=((*cal_post.gain[1])[256:383,tile_use])
            cross_poly[0:255,tile_use,1,obs_i]=((*cal_pre.gain[1])[0:255,tile_use])
            
          endfor
          
        endif
        
        ;*************************Unsplit pointing calc block
        cal_mode_fit=1
        cal_cable_reflection_fit=150
        if keyword_set(one_mode) then undefine, cal_mode_fit, cal_cable_reflection_fit
        
        cal_polyfit=vis_cal_polyfit_pointing(cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=2,phase_degree=1,$
          file_path=dir_name+'/',cal_cable_reflection_fit=cal_cable_reflection_fit,cal_mode_fit=cal_mode_fit,pointingmode=pointingmode)
          
        if keyword_set(one_mode) then begin
          poly_input_gains=complex(FLTARR(384,128,2,parsednumbers[j]))
          poly_input_gains[*,*,*,*]=1.
          
          for obs_i=0, parsednumbers[j]-1 do begin
          
            tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
            
            poly_input_gains[*,tile_use,0,obs_i]=((*cal_polyfit[obs_i].gain[0])[*,tile_use])
            poly_input_gains[*,tile_use,1,obs_i]=((*cal_polyfit[obs_i].gain[1])[*,tile_use])
            
          endfor
          
          if keyword_set(compare) then begin
            compare_ratio_array=abs(poly_input_gains)/abs(cross_poly)
            compare_ratio=FLTARR(128,2)
            for pol_j=0,1 do begin
              for tile_i=0,127 do begin
                resistant_mean,compare_ratio_array[*,tile_i,pol_j,*],2,res_mean
                compare_ratio[tile_i,pol_j]=res_mean
              endfor
            endfor
          endif
          
          cal_polyfit=vis_cal_modefit_only_pointing(poly_input_gains,cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=2,phase_degree=1,$
            file_path=dir_name+'/',cal_cable_reflection_mode_fit=1,compare_ratio=compare_ratio,cal_cable_reflection_fit=150,cal_cable_reflection_both=1,pointingmode=pointingmode)
            
        endif
        
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