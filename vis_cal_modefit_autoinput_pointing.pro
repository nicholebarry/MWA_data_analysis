FUNCTION vis_cal_modefit_autoinput_pointing,polyfit_gains,cal_array,obs_array,obsid_count,pointing_num,degree=degree,phase_degree=phase_degree,$
    file_path=file_path, cal_cable_reflection_mode_fit=cal_cable_reflection_mode_fit
  ;My program to get a polyfit over a pointing for a tile and not by obs id.  This is to avoid sidelobe effects we may be seeing.
  ;polyfit_gains is freq x tile x pol x obs
    
  cal_cable_reflection_fit=150
  
  cal_mode_fit=1
  tile_use=PTRARR(obsid_count,/allocate)
  nt_use=INTARR(obsid_count)
  
  ;These should be const across all obs.
  n_pol=cal_array[0].n_pol
  n_freq=cal_array[0].n_freq
  n_tile=cal_array[0].n_tile
  freq_arr=cal_array[0].freq
  freq_use=where((*obs_array[0].baseline_info).freq_use,nf_use)
  gain_arr_ptr=PTRARR(obsid_count,n_pol,/allocate)
  cal_return=cal_array
  
  For obs_i=0, obsid_count-1 do begin
  
    IF N_Elements(obs_array[obs_i]) GT 0 THEN *tile_use[obs_i]=where((*obs_array[obs_i].baseline_info).tile_use,nt_use_temp) ELSE tile_use=lindgen(n_tile)
    nt_use[obs_i]=nt_use_temp
    
    FOR pol_i=0,n_pol-1 DO BEGIN
      cal_return[obs_i].gain[pol_i]=Ptr_new(*cal_array[obs_i].gain[pol_i]) ;essential, to keep original cal gains from being overwritten!
      *gain_arr_ptr[obs_i,pol_i]=((*cal_array[obs_i].gain[pol_i])[*,*])    ;freq x tile
    endfor
    
  endfor
  
  c_light=299792458.
  i_comp=Complex(0,1)
  
  gain_residual=ptrarr(n_pol,n_tile,obsid_count)
  
  undefine, gain_arr
  
  IF Keyword_Set(cal_mode_fit) THEN BEGIN
    CASE 1 OF
      Keyword_Set(cal_cable_reflection_fit): BEGIN
      
        cable_filepath=filepath(obs_array[0].instrument+'_cable_length.txt',root=rootdir('FHD'),subdir='instrument_config')
        textfast,data_array,/read,file_path=cable_filepath,first_line=1
        tile_i_file=Reform(data_array[0,*])
        tile_name_file=Reform(data_array[1,*])
        cable_len=Reform(data_array[2,*])
        cable_vf=Reform(data_array[3,*])
        tile_ref_flag=0>Reform(data_array[4,*])<1
        
        IF cal_cable_reflection_fit GT 1 THEN BEGIN
          cable_cut_i=where(cable_len NE cal_cable_reflection_fit,n_cable_cut)
          IF n_cable_cut GT 0 THEN tile_ref_flag[cable_cut_i]=0
        ENDIF ELSE IF cal_cable_reflection_fit LT -1 THEN BEGIN
          cable_cut_i=where(cable_len EQ Abs(cal_cable_reflection_fit),n_cable_cut)
          IF n_cable_cut GT 0 THEN tile_ref_flag[cable_cut_i]=0
        ENDIF
        
        reflect_time=2.*cable_len/(c_light*cable_vf)
        bandwidth=(Max(freq_arr)-Min(freq_arr))*n_freq/(n_freq-1)
        mode_i_arr=Fltarr(n_pol,n_tile)
        FOR pol_i=0,n_pol-1 DO mode_i_arr[pol_i,*]=bandwidth*reflect_time*tile_ref_flag
      END
      (cal_mode_fit EQ -1): BEGIN
        spec_mask=fltarr(n_freq)
        spec_mask[freq_use]=1
        freq_cut=where(spec_mask EQ 0,n_mask)
        spec_psf=(Abs(FFT(spec_mask)))
        spec_inds=lindgen(n_freq/2)
        spec_psf=spec_psf[spec_inds]
        spectrum=fltarr(n_freq/2)
        FOR pol_i=0,n_pol-1 DO BEGIN
          FOR ti=0L,nt_use-1 DO BEGIN
            tile_i=tile_use[ti]
            spec0=Abs(FFT(*gain_residual[pol_i,tile_i]))
            spectrum+=spec0[spec_inds]
          ENDFOR
        ENDFOR
        mode_test=spectrum
        psf_mask=fltarr(n_freq/2)
        IF n_mask GT 0 THEN BEGIN
          psf_mask[where(spec_psf GT Max(spec_psf)/1E3)]=1
          psf_mask=smooth(psf_mask,5,/edge_truncate)
          mask_i=where(psf_mask,n_mask2)
          IF n_mask2 GT 0 THEN mode_test[mask_i]=0
        ENDIF
        mode_max=Max(mode_test,mode_i)
        mode_i_arr=Fltarr(n_pol,n_tile)+mode_i
      END
      ELSE: mode_i_arr=Fltarr(n_pol,n_tile)+cal_mode_fit
    ENDCASE
    
    
    mode_params=PTRARR(128,obsid_count)
    FOR pol_i=0,n_pol-1 DO BEGIN
    
      for obs_i=0,obsid_count-1 do begin
        restore,'/nfs/eor-00/h1/nbarry/Aug23_autos_onemode/'+obs_array[obs_i].obsname+'_cal.sav'
        for tile_i=0,127 do begin
          if (cal_array[obs_i].mode_params[pol_i,tile_i]) NE !NULL then mode_params[tile_i,obs_i]=PTR_NEW((*cal_array[obs_i].mode_params[pol_i,tile_i]))
        endfor
        
        
        gain_arr=*cal_array[obs_i].gain[pol_i]
        ;        gain_amp=Abs(gain_arr)
        ;        gain_phase=Atan(gain_arr,/phase)
        gain_arr_fit=polyfit_gains[*,*,pol_i,obs_i]
        gain_arr-=gain_arr_fit ; Subtract the polyfit outright so they don't talk to one another
        FOR ti=0L,nt_use[obs_i]-1 DO BEGIN
          tile_i=(*tile_use[obs_i])[ti]
          If mode_params[tile_i,obs_i] NE !NULL then mode_i=(*mode_params[tile_i,obs_i])[0] else mode_i=0
          IF mode_i EQ 0 THEN CONTINUE
          IF Keyword_Set(cal_cable_reflection_mode_fit) THEN BEGIN
            ; We are going to fit the actual mode to subtract.
            mode_use=(*mode_params[tile_i,obs_i])[0]
            phase_use=(*mode_params[tile_i,obs_i])[2]
            
            period=n_freq/mode_use
            minmax_ind=INTARR(floor(n_freq/period),2)
            minmax_ind_use=INTARR(floor(n_freq/period),2)
            
            for period_i=0,floor(n_freq/period)-1 do begin
              minmax_values=minmax(real_part(exp(-2.*!Pi*i_comp*mode_use*(findgen(floor(period))+floor(period_i*period))/n_freq)*exp(+i_comp*phase_use)),minmax_ind_temp)
              minmax_ind[period_i,*]=minmax_ind_temp+floor(period_i*period)

              if where(minmax_ind[period_i,0] EQ freq_use) NE -1 then minmax_ind_use[period_i,0]=1 else minmax_ind_use[period_i,0]=0
              if where(minmax_ind[period_i,1] EQ freq_use) NE -1 then minmax_ind_use[period_i,1]=1 else minmax_ind_use[period_i,1]=0
              
            endfor
            
            unflagged_ind=where(minmax_ind_use EQ 1, n_unflagged_ind)
            
            resistant_mean, abs(gain_arr[freq_use,tile_i]), 2,gain_mean
            resistant_mean,abs(real_part((gain_arr[minmax_ind[unflagged_ind],tile_i]-gain_mean)*exp(2.*!Pi*i_comp*mode_i*minmax_ind[unflagged_ind]/n_freq)*exp(-i_comp*phase_use))),2,amp_use
            stop
          ENDIF ELSE IF Keyword_Set(amp_arr) OR Keyword_Set(phase_arr) THEN BEGIN
            amp_use=amp_arr[pol_i,tile_i]
            phase_use=phase_arr[pol_i,tile_i]
          ENDIF ELSE BEGIN
            mode_fit=Total(exp(i_comp*2.*!Pi/n_freq*(mode_i)*freq_use)*Reform(gain_arr[freq_use,tile_i]))
            amp_use=abs(mode_fit)/nf_use ;why factor of 2? Check FFT normalization
            phase_use=atan(mode_fit,/phase)
          ENDELSE
          gain_mode_fit=amp_use*exp(-i_comp*2.*!Pi*(mode_i*findgen(n_freq)/n_freq)+i_comp*phase_use)
          gain_arr_fit[*,tile_i]+=gain_mode_fit
          cal_return.mode_params[pol_i,tile_i]=Ptr_new([mode_i,amp_use,phase_use])
          debug=1
        ENDFOR
        *cal_return.gain[pol_i]=gain_arr_fit
      endfor ;end obs for
      
      
    endfor ;end pol for
    
  ENDIF
  
  undefine_fhd,gain_residual
  RETURN,cal_return
END
