FUNCTION vis_cal_polyfit_pointing,cal_array,obs_array,obsid_count,pointing_num,degree=degree,phase_degree=phase_degree,$
    file_path=file_path, cal_cable_reflection_mode_fit=cal_cable_reflection_mode_fit,cal_cable_reflection_fit=cal_cable_reflection_fit, $
    cal_mode_fit=cal_mode_fit,pointingmode=pointingmode
  ;My program to get a polyfit over a pointing for a tile and not by obs id.  This is to avoid sidelobe effects we may be seeing.
    
  IF N_Elements(degree) EQ 0 THEN degree=2 ELSE BEGIN
    IF degree LE 0 THEN BEGIN
      amp_free_fit=1
      degree=Round(abs(degree))>1
    ENDIF ELSE BEGIN
      amp_free_fit=0
      degree=Round(degree)>1
    ENDELSE
  ENDELSE
  IF N_Elements(phase_degree) EQ 0 THEN phase_degree=1.
  
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
  FOR pol_i=0,n_pol-1 DO BEGIN
  
  
  
  
    gain_arr_ptr_per_tile=PTRARR(n_tile,/allocate)
    phase_arr_ptr_per_tile=PTRARR(n_tile,/allocate)
    freq_use_ptr=PTRARR(n_tile,/allocate)
    undefine, gain_arr
    
    For tile_i=0, n_tile-1 do begin
      freq_use_count=0
      
      FOR obs_i=0, obsid_count-1 do begin
      
        IF where(tile_i EQ *tile_use[obs_i]) NE -1 then begin
        
          freq_use_count=freq_use_count+1
          
          If ~keyword_set(gain_arr) then begin
            gain_arr=reform((*gain_arr_ptr[obs_i,pol_i])[freq_use[0],tile_i])
            phase_arr=PhUnwrap(Atan(reform((*gain_arr_ptr[obs_i,pol_i])[freq_use,tile_i]),/phase))
            for freq_i=1, nf_use-1 do gain_arr=[gain_arr,reform((*gain_arr_ptr[obs_i,pol_i])[freq_use[freq_i],tile_i])]
          endif else begin
            for freq_i=0, nf_use-1 do gain_arr=[gain_arr,reform((*gain_arr_ptr[obs_i,pol_i])[freq_use[freq_i],tile_i])]
            phase_arr=[phase_arr,PhUnwrap(Atan(reform((*gain_arr_ptr[obs_i,pol_i])[freq_use,tile_i]),/phase))]
          endelse
          
        endif
        
      endfor ;end obs loop
      
      If keyword_set(gain_arr) then *gain_arr_ptr_per_tile[tile_i]=gain_arr
      If keyword_set(phase_arr) then *phase_arr_ptr_per_tile[tile_i]=phase_arr
      If freq_use_count GT 0 then begin
        *freq_use_ptr[tile_i]=freq_use
        for iter_i=1,freq_use_count-1 do *freq_use_ptr[tile_i]=[*freq_use_ptr[tile_i],freq_use]
      endif
      undefine,gain_arr,phase_arr
      
    endfor ;end tile loop
    
    
    fake_gain=complex(FLTARR(nf_use))
    fake_gain[*]=1.
    gain_arr=complex(FLTARR(n_freq,n_tile))
    FOR tile_i=0L,n_tile-1 DO BEGIN
      IF *gain_arr_ptr_per_tile[tile_i] EQ !Null then *gain_arr_ptr_per_tile[tile_i]=fake_gain
      IF *phase_arr_ptr_per_tile[tile_i] EQ !Null then *phase_arr_ptr_per_tile[tile_i]=Atan(fake_gain,/phase)
      IF *freq_use_ptr[tile_i] EQ !Null then *freq_use_ptr[tile_i]=freq_use
      gain_amp=Abs(*gain_arr_ptr_per_tile[tile_i])
      ;gain_phase=Atan(*gain_arr_ptr_per_tile[tile_i],/phase)
      fit_params=poly_fit(*freq_use_ptr[tile_i],gain_amp,degree) ;gain needs to be a 1D array of a tile's gain for the obs of the pointing for all freq. freq_use needs to be a repeated 1D array to match gain
      cal_return.amp_params[pol_i,tile_i]=Ptr_new(fit_params)
      
      gain_fit=fltarr(n_freq)
      FOR di=0L,degree DO gain_fit+=fit_params[di]*findgen(n_freq)^di
      
      ;gain_residual[pol_i,tile_i]=Ptr_new(Reform(gain_amp[*,tile_i])-gain_fit)
      
      IF phase_degree GT 0 THEN BEGIN
        phase_use=*phase_arr_ptr_per_tile[tile_i]
        ;phase_use=PhUnwrap(reform(gain_phase[freq_use,tile_i]))
        phase_params=poly_fit(*freq_use_ptr[tile_i],phase_use,phase_degree,yfit=phase_fit)
        cal_return.phase_params[pol_i,tile_i]=Ptr_new(phase_params)
        phase_fit=fltarr(n_freq)
        FOR di=0L,phase_degree DO phase_fit+=phase_params[di]*findgen(n_freq)^di
        
        
        gain_arr[*,tile_i]=gain_fit*Exp(i_comp*phase_fit)

      ;stop ;no failsafes for skipped tiles yet
      ENDIF ELSE gain_arr[*,tile_i]=gain_fit;*weight_invert(gain_amp[*,tile_i]) ;this preserves the original phase
    ENDFOR ;end tile loop
    For obs_i=0, obsid_count-1 do *cal_return[obs_i].gain[pol_i]=gain_arr
  ENDFOR
  
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
    
    
    
    ;******************Mode fit by obs calc block
    mode_params=PTRARR(128,obsid_count)
    FOR pol_i=0,n_pol-1 DO BEGIN
    
      For obs_i=0, obsid_count-1 do begin
      
        gain_arr=*cal_array[obs_i].gain[pol_i]
        gain_arr_fit=*cal_return[obs_i].gain[pol_i]
        gain_arr-=gain_arr_fit ; Subtract the polyfit outright so they don't talk to one another
        FOR ti=0L,nt_use[obs_i]-1 DO BEGIN
          tile_i=(*tile_use[obs_i])[ti]
          mode_i=mode_i_arr[pol_i,tile_i]
          IF mode_i EQ 0 THEN CONTINUE
          IF Keyword_Set(cal_cable_reflection_mode_fit) THEN BEGIN
            ; We are going to fit the actual mode to subtract.
            mode0=mode_i ; start with nominal cable length
            dmode=0.05 ; pretty fine
            nmodes=101 ; range around the central mode to test
            modes=(dindgen(nmodes)-nmodes/2)*dmode+mode0 ; array of modes to try
            modes=rebin(modes,nmodes,nf_use) ; hopefully this is right...
            gainr=rebin(transpose(reform(real_part(gain_arr[freq_use,tile_i]))),nmodes,nf_use)
            gaini=rebin(transpose(reform(imaginary(gain_arr[freq_use,tile_i]))),nmodes,nf_use) ; and this...
            gain_temp=gainr+i_comp*gaini ; for some reason I cant rebin complex numbers
            freq_mat=rebin(transpose(freq_use),nmodes,nf_use) ; this too...
            test_fits=Total(exp(i_comp*2.*!Pi/n_freq*modes*freq_mat)*gain_temp,2)
            amp_use=max(abs(test_fits),mode_ind)/nf_use
            phase_use=atan(test_fits[mode_ind],/phase)
            mode_i=modes[mode_ind,0]
            
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
          mode_params[tile_i,obs_i]=PTR_NEW([mode_i,amp_use,phase_use])
        ENDFOR
        if ~keyword_set(pointingmode) then *cal_return[obs_i].gain[pol_i]=gain_arr_fit
      ENDFOR
      ;******************end of modefit by obs
      
      if keyword_set(pointingmode) then begin
        ;******************Calc median mode block
        mode_median=FLTARR(128)
        mode_temp=FLTARR(obsid_count)
        
        FOR ti=0L,127 DO BEGIN
          ;fill a temporary array with all the modes of a single tile over the pointings obsids
          for obs_i=0,obsid_count-1 do begin
            if mode_params[ti,obs_i] NE !Null then mode_temp[obs_i]=(*mode_params[ti,obs_i])[0] else mode_temp[obs_i]=0
          endfor
          
          ;Begin a quartile analysis of the mode array over obsids. This is the same as a boxplot analysis, and assumes no prior distribution to the data
          ;Start with the median value over obsid
          medianVal=Median(mode_temp)
          
          ;If the median value is 0, then it is a tile without a mode, and should fill the median array with a 0
          IF medianVal NE 0 then begin
          
            ;Find the indices that correspond to the 25% and 75% quartile, and count how many are in each
            qtr_25th_index = Where(mode_temp LT medianVal, countlowerhalf)
            qtr_75th_index = Where(mode_temp GT medianVal, countupperhalf)
            
            ;Unless the amount in the quartile is 0, find the median within that quartile, which corresponds to the 25% and 75% value.
            If countlowerhalf EQ 0 then qtr_25th=medianVal else qtr_25th=Median(mode_temp[qtr_25th_index])
            If countupperhalf EQ 0 then qtr_75th=medianVal else qtr_75th=Median(mode_temp[qtr_75th_index])
            
            ;Find the bottom and top outlier cutoff value given the inner quartile range.
            bottomoutlier=qtr_25th-1.*(qtr_75th-qtr_25th)
            topoutlier=qtr_75th+1.*(qtr_75th-qtr_25th)
            
            ;Find the indicies that are above the bottom outlier cutoff and below the top outlier cutoff. Calculate the mean from that subset
            topindex=where(mode_temp LE topoutlier,meancount1)
            If meancount1 GT 0 then indexformean=where(mode_temp[topindex] GE bottomoutlier, meancount2) else meancount2=0
            If (meancount1 GT 0) AND (meancount2 GT 0) then box_mean=mean(mode_temp[indexformean]) else box_mean=0
            
          endif else box_mean=medianVal
          mode_median[ti]=box_mean
        ;if (ti EQ 14) and (pol_i EQ 0) then stop
        endfor
        ;****************end of calc median mode block
        
        
        ;****************run through with a fixed mode
        For obs_i=0, obsid_count-1 do begin
        
          ;Fill the arrays with the original gain (bp subtracted) and the polyfit
          gain_arr=*cal_array[obs_i].gain[pol_i]
          gain_arr_fit=*cal_return[obs_i].gain[pol_i]
          
          ; Subtract the polyfit outright so the modefit and polyfit don't talk to one another
          gain_arr-=gain_arr_fit
          
          FOR ti=0L,nt_use[obs_i]-1 DO BEGIN
          
            tile_i=(*tile_use[obs_i])[ti]
            ;Fix the mode to the median value calculated above
            mode_i=mode_median[tile_i]
            
            ;If the mode for that tile is 0, restart the tile loop and advance a tile index. This avoids the calculation of a modefit
            ;if the tile doesn't have the specified reflection. This also leaves a !Null pointer in the mode_params pointer array,
            ;which will be used extensively later to pinpoint tiles with reflection fits.
            IF mode_i EQ 0 THEN CONTINUE
            mode_fit=Total(exp(i_comp*2.*!Pi/n_freq*(mode_i)*freq_use)*Reform(gain_arr[freq_use,tile_i]))
            
            ;Find the real and imaginary part of the modefit (instead of the amp and phase) to avoid phase wrapping issues
            ;in taking the mean over the pointing later. Also divide by the number of frequencies used to normalize
            ;the IDL convention for the fourier transform.
            real_use=real_part(mode_fit)/nf_use
            imag_use=imaginary(mode_fit)/nf_use
            
            *mode_params[tile_i,obs_i]=[mode_i,real_use,imag_use]
          ENDFOR
          
        ENDFOR
        ;***************end of fixed mode
        
        ;******************Calc mean amp and phase block
        real_mean=FLTARR(128)
        imag_mean=FLTARR(128)
        amp_mean=FLTARR(128)
        phase_mean=FLTARR(128)
        real_temp=FLTARR(obsid_count)
        imag_temp=FLTARR(obsid_count)
        
        FOR ti=0L,127 DO BEGIN
          ;fill a temporary array with all the amp and phase of a single tile over the pointings obsids
          for obs_i=0,obsid_count-1 do begin
            if mode_params[ti,obs_i] NE !Null then real_temp[obs_i]=(*mode_params[ti,obs_i])[1] else real_temp[obs_i]=0
            if mode_params[ti,obs_i] NE !Null then imag_temp[obs_i]=(*mode_params[ti,obs_i])[2] else imag_temp[obs_i]=0
          endfor
          
          ;Begin a quartile analysis of the real array over obsids. This is the same as a boxplot analysis, and assumes no prior distribution to the data
          ;Start with the median value over obsid
          medianVal=Median(real_temp)
          
          ;If the median value is 0, then it is a tile without a reflection fit, and should fill the median array with a 0
          IF medianVal NE 0 then begin
            ;Find the indices that correspond to the 25% and 75% quartile, and count how many are in each
            qtr_25th_index = Where(real_temp LT medianVal, countlowerhalf)
            qtr_75th_index = Where(real_temp GT medianVal, countupperhalf)
            
            ;Unless the amount in the quartile is 0, find the median within that quartile, which corresponds to the 25% and 75% value.
            If countlowerhalf EQ 0 then qtr_25th=medianVal else qtr_25th=Median(real_temp[qtr_25th_index])
            If countupperhalf EQ 0 then qtr_75th=medianVal else qtr_75th=Median(real_temp[qtr_75th_index])
            
            ;Find the bottom and top outlier cutoff value given the inner quartile range.
            bottomoutlier=qtr_25th-1.*(qtr_75th-qtr_25th)
            topoutlier=qtr_75th+1.*(qtr_75th-qtr_25th)
            
            ;Find the indicies that are above the bottom outlier cutoff and below the top outlier cutoff. Calculate the mean from that subset
            topindex=where(real_temp LE topoutlier,meancount1)
            If meancount1 GT 0 then indexformean=where(real_temp[topindex] GE bottomoutlier, meancount2) else meancount2=0
            If (meancount1 GT 0) AND (meancount2 GT 0) then box_mean=mean(real_temp[indexformean]) else box_mean=0
            
          endif else box_mean=medianVal
          real_mean[ti]=box_mean
          
          ;Begin a quartile analysis of the imaginary array over obsids. This is the same as a boxplot analysis, and assumes no prior distribution to the data
          ;Start with the median value over obsid
          medianVal=Median(imag_temp)
          
          ;If the median value is 0, then it is a tile without a reflection fit, and should fill the median array with a 0
          IF medianVal NE 0 then begin
            ;Find the indices that correspond to the 25% and 75% quartile, and count how many are in each
            qtr_25th_index = Where(imag_temp LT medianVal, countlowerhalf)
            qtr_75th_index = Where(imag_temp GT medianVal, countupperhalf)
            
            ;Unless the amount in the quartile is 0, find the median within that quartile, which corresponds to the 25% and 75% value.
            If countlowerhalf EQ 0 then qtr_25th=medianVal else qtr_25th=Median(imag_temp[qtr_25th_index])
            If countupperhalf EQ 0 then qtr_75th=medianVal else qtr_75th=Median(imag_temp[qtr_75th_index])
            
            ;Find the bottom and top outlier cutoff value given the inner quartile range.
            bottomoutlier=qtr_25th-1.*(qtr_75th-qtr_25th)
            topoutlier=qtr_75th+1.*(qtr_75th-qtr_25th)
            
            ;Find the indicies that are above the bottom outlier cutoff and below the top outlier cutoff. Calculate the mean from that subset
            topindex=where(imag_temp LE topoutlier,meancount1)
            If meancount1 GT 0 then indexformean=where(imag_temp[topindex] GE bottomoutlier, meancount2) else meancount2=0
            If (meancount1 GT 0) AND (meancount2 GT 0) then box_mean=mean(imag_temp[indexformean]) else box_mean=0
            
          endif else box_mean=medianVal
          imag_mean[ti]=box_mean
          
          ;Calculate the amp and phase from the mean of the real and imaginary parts of the fit
          amp_mean[ti]=sqrt(real_mean[ti]^2.+imag_mean[ti]^2.)
          phase_mean[ti]=atan(real_mean[ti]+i_comp*imag_mean[ti],/phase)
        endfor
        ;****************end of calc mean amp and phase block
        
        for obs_i=0, obsid_count-1 do begin
          gain_arr=*cal_array[obs_i].gain[pol_i]
          gain_arr_fit=*cal_return[obs_i].gain[pol_i]
          gain_arr-=gain_arr_fit ; Subtract the polyfit outright so they don't talk to one another
          FOR ti=0L,nt_use[obs_i]-1 DO BEGIN
            tile_i=(*tile_use[obs_i])[ti]
            mode_i=mode_median[tile_i]
            IF mode_i EQ 0 THEN CONTINUE
            
            gain_mode_fit=amp_mean[tile_i]*exp(-i_comp*2.*!Pi*(mode_median[tile_i]*findgen(n_freq)/n_freq)+i_comp*phase_mean[tile_i])
            gain_arr_fit[*,tile_i]+=gain_mode_fit
            cal_return[obs_i].mode_params[pol_i,tile_i]=Ptr_new([mode_median[tile_i],amp_mean[tile_i],phase_mean[tile_i]])
            
          ENDFOR
          *cal_return[obs_i].gain[pol_i]=gain_arr_fit
        endfor
      endif
    endfor ;end pol for
    
  ENDIF
  
  undefine_fhd,gain_residual
  RETURN,cal_return
END
