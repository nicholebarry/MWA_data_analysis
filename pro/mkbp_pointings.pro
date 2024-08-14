pro mkbp_pointings, dir_name,day=day, advanced_plotting=advanced_plotting, longrun=longrun, split_jump=split_jump, save_autos=save_autos,one_mode=one_mode,phase_transfer=phase_transfer
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
        
        ;I think this makes sure the obs actually exists, as there might be obs in the text file that passed tests but didn't get run? I think this is superfluous.
        IF obsid NE 'empty' THEN BEGIN
        
          ;Restore obs save file, either from a regular run or a longrun (different obs potentially)
          If keyword_set(longrun) then restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/' + obsid + '_obs.sav' else $
            ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/metadata/' + obsid + '_obs.sav'
            restore, '/nfs/eor-03/r1/EoR2013/fhd_nb_catalog_July2015/metadata/' + obsid + '_obs.sav'
            
          ;Setup the obs structure array and fill it on successive loops
          If (i eq 0) Then obs_array = replicate(obs,parsednumbers[j]) ELSE obs_array[i]=obs
          
          ;Restore params save file, either from regular run or a longrun
          If keyword_set(longrun) then restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/' + obsid + '_params.sav' else $
            ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/metadata/' + obsid + '_params.sav'
            restore, '/nfs/eor-03/r1/EoR2013/fhd_nb_catalog_July2015/metadata/1061311664_params.sav'
            
          ;Setup the params structure array and fill it on successive loops
          ;If (i eq 0) Then params_array = replicate(params,parsednumbers[j]) ELSE params_array[i]=params
          
          ;Restore cal structure from longrun or regular run
          If keyword_set(longrun) then restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/' + obsid + '_cal.sav' else $
            ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/calibration/' + obsid + '_cal.sav'
            restore, '/nfs/eor-03/r1/EoR2013/fhd_nb_catalog_July2015/calibration/' + obsid + '_cal.sav'
            
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
          
          if keyword_set(split_jump) then begin
            cal_bandpass=vis_cal_bandpass(cal_array[i],obs_array[i],cal_remainder=cal_remainder,saved_run_bp=1,cable_bandpass_fit=1)
            for pol_i=0,n_pol[i]-1 do *cal_array[i].gain[pol_i]=*cal_remainder.gain[pol_i]

          endif
        ENDIF
      ENDFOR
      ;****End of restore loop
      
      file_mkdir, dir_name, /NOEXPAND_PATH
      bp=vis_cal_bandpass_with_cable_pointingv2(cal_array,obs_array,parsednumbers[j],pointing_num[j],cal_remainder=cal_remainder,file_path_fhd=dir_name+'/',advanced_plotting=advance_plotting)
      
      ;bp=vis_cal_polyfit_pointing(cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=2,phase_degree=1,$
      ;  file_path=dir_name+'/')
      
      If keyword_set(split_jump) then begin
        ;Only a bandpass within the specified frequencies are calculated, so a only pre dig jump looks truncated
      
      
        ;*************************Split pointing calc block
        ;Fit parameters are only calculated in the specified frequencies, but the final fit is over all frequencies.
        ;Pre dig jump calcs should be fine, even though the product looks like it's done over all freq
        ;Just fitting linear trend, bc don't want to overfit and it will be interesting to see if degree=2 was just trying to capture dig jump
        cal_mode_fit=1
        cal_cable_reflection_fit=150
        if keyword_set(one_mode) then undefine, cal_mode_fit, cal_cable_reflection_fit
        
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
  
        if keyword_set(one_mode) then begin
          poly_input_gains=complex(FLTARR(384,128,2,parsednumbers[j]))
          poly_input_gains[*,*,*]=1.
          
          for obs_i=0, parsednumbers[j]-1 do begin
          
            tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
            
            poly_input_gains[256:383,tile_use,0,obs_i]=((*cal_polyfit_post[obs_i].gain[0])[256:383,tile_use])
            poly_input_gains[0:255,tile_use,0,obs_i]=((*cal_polyfit_pre[obs_i].gain[0])[0:255,tile_use])
            
            poly_input_gains[256:383,tile_use,1,obs_i]=((*cal_polyfit_post[obs_i].gain[1])[256:383,tile_use])
            poly_input_gains[0:255,tile_use,1,obs_i]=((*cal_polyfit_pre[obs_i].gain[1])[0:255,tile_use])
            
          endfor
          
          cal_polyfit_total=vis_cal_modefit_only_pointing(poly_input_gains,cal_array,obs_array,parsednumbers[j],pointing_num[j],degree=1,phase_degree=1,$
            file_path=dir_name+'/',cal_cable_reflection_mode_fit=1)
        endif
        ;*************************Split pointing calc block
        
        
        ;for saving the cal structures individually for restoring, especially with adams fourth line plots
        if ~keyword_set(one_mode) then begin
          For i=0, parsednumbers[j]-1 DO BEGIN
            obsid=(*obs_ptr[j])[i]
            cal=cal_polyfit_pre[i]
            save, cal, FILENAME=dir_name+'/'+obsid+'pre_cal.sav
            cal=cal_polyfit_post[i]
            save, cal, FILENAME=dir_name+'/'+obsid+'post_cal.sav
          endfor
        endif else begin
          For i=0, parsednumbers[j]-1 DO BEGIN
            obsid=(*obs_ptr[j])[i]
            cal=cal_polyfit_total[i]
            save, cal, FILENAME=dir_name+'/'+obsid+'_cal.sav
          endfor
        endelse
        
        
        
        mode_params_pre=cal_polyfit_pre.mode_params
        mode_params_post=cal_polyfit_post.mode_params
        
        ;Restoring this obs structure just to get the cable lengths of the tiles, and storing the indices of the 150m cables
        restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/metadata/1061316296_obs.sav'
        mode_filepath=filepath('mwa_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
        textfast,data_array,/read,file_path=mode_filepath,first_line=1
        cable_len=Reform(data_array[2,*])
        tile_use_cable=where((*obs.baseline_info).tile_use AND cable_len EQ 150)
        
        ;Initialize complex arrays
        post_modefit=complex(FLTARR(384,N_elements(tile_use_cable)))
        pre_modefit=complex(FLTARR(384,N_elements(tile_use_cable)))
        amp_ratio_data=FLTARR(N_elements(tile_use_cable))
        n_freq=384.
        
        ;If a histogram plot is to be made of the ratio of amplitudes, this loop will give you that data. Serves no other purpose
        for tile_i=0,N_elements(tile_use_cable)-1 do begin
          mode_params_tile_post=*mode_params_post[0,tile_use_cable[tile_i],0]
          mode_params_tile_pre=*mode_params_pre[0,tile_use_cable[tile_i],0]
          post_modefit[*,tile_i]=mode_params_tile_post[1]*Exp(-Complex(0,1)*2.*!Pi*mode_params_tile_post[0]*findgen(384)/n_freq+Complex(0,1)*(mode_params_tile_post[2]))
          pre_modefit[*,tile_i]=mode_params_tile_pre[1]*Exp(-Complex(0,1)*2.*!Pi*mode_params_tile_pre[0]*findgen(384)/n_freq+Complex(0,1)*(mode_params_tile_pre[2]))
          amp_ratio_data[tile_i]=(minmax(abs(post_modefit[*,tile_i])))[0]/(minmax(abs(pre_modefit[*,tile_i])))[0]
        endfor
        ;cgHistoplot, amp_ratio_data, BINSIZE=.1, /FILL
        
        
        
        freq_arr=cal_array[0].freq
        
        ;Initialize arrays for final cal solutions calculated.  Fill with 1 as filler, due to flagged tiles needing a placeholder
        poly_input_gains=complex(FLTARR(384,128,2))
        poly_input_gains[*,*,*]=1.
        
        ;freq x tiles x obs x pol
        for obs_i=0, parsednumbers[j]-1 do begin
        
          tile_use=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp)
          
          poly_input_gains[256:383,tile_use,0]=((*cal_polyfit_post[obs_i].gain[0])[256:383,tile_use])
          poly_input_gains[0:255,tile_use,0]=((*cal_polyfit_pre[obs_i].gain[0])[0:255,tile_use])
          poly_input_gains[256:383,tile_use,1]=((*cal_polyfit_post[obs_i].gain[1])[256:383,tile_use])
          poly_input_gains[0:255,tile_use,1]=((*cal_polyfit_pre[obs_i].gain[1])[0:255,tile_use])
          
        endfor
        
        ;plot_split_jump=1
        if keyword_set(plot_split_jump) then begin
          tile_names=(*obs_array[0].baseline_info).tile_names
          obsname=strtrim(string(pointing_num[j]),2)
          filename=dir_name+'/'+obsname+'_poly.png'
          save, poly_input_gains, filename=dir_name+'/'+obsname+'_poly_input_gains.sav'
          ;tile_use=where((*obs_array[obs_i].baseline_info).tile_use)
          tile_use=(*obs_array[0].baseline_info).tile_use ; just first for now
          tile_exist=(*obs_array[0].baseline_info).tile_use
          
          plot_pos=calculate_plot_positions(.8, nplots=128, /no_colorbar, ncol=17, nrow=9, plot_margin = [.05, .05, .02, .25]); use 1 extra row/column to leave room for title and axes
          plot_pos=reform(plot_pos.plot_pos,17,9,4)
          ; Now remove unwanted row/column
          plot_pos = plot_pos[1:*,*,*]
          plot_pos = plot_pos[*,1:*,*]
          plot_pos = reform(plot_pos,128,4)
          ; Shift a half a width over and up
          width = plot_pos[1,0]-plot_pos[0,0]
          height = abs(plot_pos[16,1]-plot_pos[0,1])
          plot_pos[*,0] -= width/2
          plot_pos[*,1] += height/2
          plot_pos[*,2] -= width/2
          plot_pos[*,3] += height/2
          cal_plot_charsize=.5
          
          plot_cals_sub_quickanddirty_ploy,freq_arr,poly_input_gains[*,*,0],poly_input_gains[*,*,0],filename=filename,phase=phase,real_vs_imaginary=real_vs_imaginary,$
            tile_A=tile_A,tile_B=tile_B,tile_use=tile_use,tile_exist=tile_exist,tile_names=tile_names,yrange=yrange,$
            obsname=obsname,plot_pos=plot_pos,cal_plot_charsize=cal_plot_charsize,cal_plot_symsize=cal_plot_symsize,cal_plot_resize=cal_plot_resize
        endif
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