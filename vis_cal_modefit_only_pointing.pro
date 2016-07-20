FUNCTION vis_cal_modefit_only_pointing,polyfit_gains,cal_array,obs_array,obsid_count,pointing_num,degree=degree,phase_degree=phase_degree,$
    file_path=file_path, cal_cable_reflection_mode_fit=cal_cable_reflection_mode_fit,mode_input=mode_input,compare_ratio=compare_ratio,$
    cal_cable_reflection_fit=cal_cable_reflection_fit,cal_cable_reflection_both=cal_cable_reflection_both,pointingmode=pointingmode,longrun_mode=longrun_mode, $
    orig_mode_flagging=orig_mode_flagging
  ;My program to get a polyfit over a pointing for a tile and not by obs id.  This is to avoid sidelobe effects we may be seeing.
  ;polyfit_gains is freq x tile x pol x obs
    
  ;************Setup
  ;
  ;Fit at least the 150m
  If ~keyword_set(cal_cable_reflection_fit) then cal_cable_reflection_fit=150
  
  ;Force fitting of the mode
  cal_mode_fit=1
  tile_use=PTRARR(obsid_count,/allocate)
  nt_use=INTARR(obsid_count)
  
  ;These should be const across all obs.
  n_pol=cal_array[0].n_pol
  n_freq=cal_array[0].n_freq
  n_tile=cal_array[0].n_tile
  freq_arr=cal_array[0].freq
  flagged_index=where((*obs_array[0].baseline_info).freq_use EQ 0)
  
  ;Optionally flag one extra channel on either side of course channels in order to remove contamination from van Vleck corrections
  If ~keyword_set(orig_mode_flagging) then begin
    for fl_i=0, (size(flagged_index))[1]-1 do begin
      If fl_i GT 0 then begin
        (*obs_array[0].baseline_info).freq_use[flagged_index[fl_i]-1]=0
      endif
      If fl_i LT (size(flagged_index))[1]-1 then begin
        (*obs_array[0].baseline_info).freq_use[flagged_index[fl_i]+1]=0
      endif
    endfor
  endif
  
  freq_use=where((*obs_array[0].baseline_info).freq_use,nf_use)
  gain_arr_ptr=PTRARR(obsid_count,n_pol,/allocate)
  cal_return=cal_array
  
  ;Create new pointers for the return structure, and find tiles to use for each observation
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
  
  ;************End of Setup
  
  
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
          
          If keyword_set(cal_cable_reflection_both) then begin
            cable_90=where(cable_len EQ 90,n_cable_cut)
            IF n_cable_cut GT 0 THEN tile_ref_flag[cable_90]=1
            
          endif
          
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
      
        If ~keyword_set(longrun_mode) then begin
          gain_arr=*cal_array[obs_i].gain[pol_i]
          gain_arr_fit=polyfit_gains[*,*,pol_i,obs_i]
          
          ; Subtract the polyfit outright so the modefit and polyfit don't talk to one another
          gain_arr-=gain_arr_fit        
        endif else begin
          
          gain_arr=polyfit_gains[*,*,pol_i,obs_i]
          gain_arr_fit=complex(FLTARR(384,128))
          gain_arr_fit[*,*]=1.

        endelse

        ;Extra removal of leftover gain TEST
        ;For obs_j=0, obsid_count-1 do begin
        ;  for tile_j=0,127 do begin
        ;    extra_removal=mean(abs(gain_arr_res[freq_use,tile_j,obs_j])
        ;    gain_arr_res[freq_use,tile_j,obs_j]-=extra_removal
        ;  endfor
        ;endfor

        for tile_j=0,127 do begin
          extra_removal=mean(real_part(gain_arr[freq_use,tile_j]))
          gain_arr[freq_use,tile_j]-=extra_removal
        endfor
        
        FOR ti=0L,nt_use[obs_i]-1 DO BEGIN
          tile_i=(*tile_use[obs_i])[ti]
          mode_i=mode_i_arr[pol_i,tile_i]
          
          
          IF mode_i EQ 0 THEN CONTINUE
          IF Keyword_Set(cal_cable_reflection_mode_fit) THEN BEGIN
            ; We are going to fit the actual mode to subtract.
          
            freq_not_use=where((*obs_array[0].baseline_info).freq_use EQ 0)
            gain_arr[freq_not_use,tile_i]=0
            tile_name=strtrim((*obs_array[0].baseline_info).tile_names[tile_i],2)
            
            If ~keyword_set(orig_mode_flagging) then begin
              for fl_i=0, (size(freq_not_use))[1]-1 do begin
                If (fl_i GT 0) AND (fl_i LT (size(freq_not_use))[1]-3) then begin
                  ;If (fl_i GT 0) AND (fl_i LT (size(freq_not_use))[1]-1) then begin
                
                  ;IF freq_not_use[fl_i+1] EQ (freq_not_use[fl_i]+1) then begin
                  IF freq_not_use[fl_i+3] EQ (freq_not_use[fl_i]+3) then begin
                  
                    slope_real=(real_part(gain_arr[freq_not_use[fl_i]+4,tile_i]) - real_part(gain_arr[freq_not_use[fl_i]-1,tile_i]))/5
                    intercept_real=real_part(gain_arr[freq_not_use[fl_i]-1,tile_i])
                    
                    slope_imag=(imaginary(gain_arr[freq_not_use[fl_i]+4,tile_i]) - imaginary(gain_arr[freq_not_use[fl_i]-1,tile_i]))/5
                    intercept_imag=imaginary(gain_arr[freq_not_use[fl_i]-1,tile_i])
                    
                    ;TEST
                    ;slope_amp=(abs(gain_arr[freq_not_use[fl_i]+4,tile_i]) - abs(gain_arr[freq_not_use[fl_i]-1,tile_i]))/5
                    ;intercept_amp=abs(gain_arr[freq_not_use[fl_i]-1,tile_i])
                    
                    ;slope_phase=(atan(gain_arr[freq_not_use[fl_i]+4,tile_i],/phase) - atan(gain_arr[freq_not_use[fl_i]-1,tile_i],/phase))/5
                    ;intercept_phase=atan(gain_arr[freq_not_use[fl_i]-1,tile_i],/phase)
                    ;END OF TEST
                    
                    ;slope_real=(real_part(gain_arr[freq_not_use[fl_i]+2,tile_i]) - real_part(gain_arr[freq_not_use[fl_i]-1,tile_i]))/3
                    ;intercept_real=real_part(gain_arr[freq_not_use[fl_i]-1,tile_i])
                    
                    ;slope_imag=(imaginary(gain_arr[freq_not_use[fl_i]+2,tile_i]) - imaginary(gain_arr[freq_not_use[fl_i]-1,tile_i]))/3
                    ;intercept_imag=imaginary(gain_arr[freq_not_use[fl_i]-1,tile_i])
                    
                    
                    gain_arr[freq_not_use[fl_i],tile_i]=(slope_real*1.+intercept_real)+i_comp*(slope_imag*1.+intercept_imag)
                    gain_arr[freq_not_use[fl_i+1],tile_i]=(slope_real*2.+intercept_real)+i_comp*(slope_imag*2.+intercept_imag)
                    gain_arr[freq_not_use[fl_i+2],tile_i]=(slope_real*3.+intercept_real)+i_comp*(slope_imag*3.+intercept_imag)
                    gain_arr[freq_not_use[fl_i+3],tile_i]=(slope_real*4.+intercept_real)+i_comp*(slope_imag*4.+intercept_imag)
                    
                  ;TEST
                  ;gain_arr[freq_not_use[fl_i],tile_i]=abs((slope_amp)*1.+(intercept_amp))*exp(i_comp*(slope_phase*1.+intercept_phase))
                  ;gain_arr[freq_not_use[fl_i+1],tile_i]=abs((slope_amp)*2.+(intercept_amp))*exp(i_comp*(slope_phase*2.+intercept_phase))
                  ;gain_arr[freq_not_use[fl_i+2],tile_i]=abs((slope_amp)*3.+(intercept_amp))*exp(i_comp*(slope_phase*3.+intercept_phase))
                  ;gain_arr[freq_not_use[fl_i+3],tile_i]=abs((slope_amp)*4.+(intercept_amp))*exp(i_comp*(slope_phase*4.+intercept_phase))
                  ;                  gain_arr[freq_not_use[fl_i],tile_i]=abs((slope_amp)*1.+(intercept_amp))
                  ;gain_arr[freq_not_use[fl_i+1],tile_i]=abs((slope_amp)*2.+(intercept_amp))
                  ;gain_arr[freq_not_use[fl_i+2],tile_i]=abs((slope_amp)*3.+(intercept_amp))
                  ;gain_arr[freq_not_use[fl_i+3],tile_i]=abs((slope_amp)*4.+(intercept_amp))
                  ;END OF TEST
                    
                  endif
                endif
              endfor
            endif
            ;**************
            
            gain_arr[0:1,tile_i]=0
            gain_arr[382:383,tile_i]=0
            
            wh_freq_nonzero=where(gain_arr[*,tile_i] NE 0, nonzero_freq_count)
            
            mode0=mode_i ;start with nominal cable length
            dmode=0.025 ;mode delta amount
            nmodes=71 ;range around the central mode to test
            modes=(dindgen(nmodes)-nmodes/2)*dmode+mode0 ; array of modes to try

            modes=rebin(modes,nmodes,nonzero_freq_count)

            gainr=rebin(transpose(reform(real_part(gain_arr[wh_freq_nonzero,tile_i]))),nmodes,nonzero_freq_count)
            gaini=rebin(transpose(reform(imaginary(gain_arr[wh_freq_nonzero,tile_i]))),nmodes,nonzero_freq_count)

            gain_temp=gainr+i_comp*gaini 

            freq_mat=rebin(transpose(Findgen(nonzero_freq_count)+2),nmodes,nonzero_freq_count)

            test_fits=Total(exp(i_comp*2.*!Pi/n_freq*modes*freq_mat)*gain_temp,2)
           
            ;test_fits2=Total(exp(i_comp*2.*!Pi*(modes/bandwidth)*(freq_mat/n_freq*bandwidth+min(freq_arr)))*gain_temp,2)
            
            
            tile_name=strtrim((*obs_array[0].baseline_info).tile_names[tile_i],2)
            
            ;If (tile_name EQ '33') AND (pol_i EQ 1) then tile33=abs(test_fits)
            
            ;If (tile_name EQ '24') and pol_i EQ 0 then stop
            ;    plot_title=strtrim(string(obs_array[obs_i].obsname),2)+'_'+tile_name+'_2interp_only'
            ;cgPS_Open,'/nfs/eor-00/h1/nbarry/Aug23_90m_trys/fakechannel_compare/'+plot_title+'.png',/quiet,/nomatch
            ; cgplot, modes[*,0],abs(test_fits),yrange=[0,20], title=plot_title
            ; cgoplot, [mode_i,mode_i], [0,10], linestyle=2
            
            ;  temp=max(abs(test_fits),mode_ind)
            ;  cgoplot, [modes[mode_ind,0],modes[mode_ind,0]], [0,10]
            
            ; cgLegend, Title=['abs(DFT)','Theory mode'],$
            ;    Color=['black','black'],Location=[0.65,0.75], charsize=1,VSpace=1.1,Length=0.05, $
            ;    linestyle=[0,2]
            ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
            ;endif
            
            ;temp=max(abs(test_fits),mode_ind)
            ;If (tile_name EQ '12') AND (strtrim(string(obs_array[obs_i].obsname),2) EQ '1061314832') then stop
            
            ;If (tile_name EQ '21') OR (tile_name EQ '25') then stop
            
            
            ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
            
            If ~keyword_set(mode_input) then begin
              ;amp_use=2*max(abs(test_fits),mode_ind)/380
              amp_use=max(abs(test_fits),mode_ind)/nonzero_freq_count
              
              If mode_ind EQ 0 then begin
                amp_use=max(abs(test_fits[round(nmodes/2):nmodes-1]),mode_ind)/nonzero_freq_count
                mode_ind=mode_ind+round(nmodes/2)
                IF mode_ind EQ nmodes-1 then amp_use=max(abs(test_fits),mode_ind)/nonzero_freq_count
              endif

              if keyword_set(compare_ratio) then amp_use=amp_use/compare_ratio[tile_i,pol_i]
              phase_use=atan(test_fits[mode_ind],/phase)
              
            endif else begin
              mode_ind_temp=min(abs((*mode_input[obs_i,pol_i,tile_i])[0] - modes[*,0]),mode_ind)
              amp_use=abs(test_fits[mode_ind])/nonzero_freq_count
              phase_use=atan(test_fits[mode_ind],/phase)
            endelse
            
            mode_i=modes[mode_ind,0]
            
          ENDIF ELSE IF Keyword_Set(amp_arr) OR Keyword_Set(phase_arr) THEN BEGIN
            amp_use=amp_arr[pol_i,tile_i]
            phase_use=phase_arr[pol_i,tile_i]
          ENDIF ELSE BEGIN
            mode_fit=Total(exp(i_comp*2.*!Pi/n_freq*(mode_i)*freq_use)*Reform(gain_arr[freq_use,tile_i]))
            amp_use=abs(mode_fit)/nf_use ;why factor of 2? Check FFT normalization
            phase_use=atan(mode_fit,/phase)
          ENDELSE
          
          if keyword_set(pointingmode) then begin
            real_use=amp_use*cos(phase_use)
            imag_use=amp_use*sin(phase_use)
          endif
          
          gain_mode_fit=amp_use*exp(-i_comp*2.*!Pi*(mode_i*findgen(n_freq)/n_freq)+i_comp*phase_use)

          gain_arr_fit[*,tile_i]+=gain_mode_fit
          
          if ~keyword_set(pointingmode) then cal_return[obs_i].mode_params[pol_i,tile_i]=PTR_NEW([mode_i,amp_use,phase_use]) else $
            cal_return[obs_i].mode_params[pol_i,tile_i]=PTR_NEW([mode_i,real_use,imag_use])
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
            if cal_return[obs_i].mode_params[pol_i,ti] NE !Null then mode_temp[obs_i]=(*cal_return[obs_i].mode_params[pol_i,ti])[0] else mode_temp[obs_i]=0
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
            If (meancount1 GT 0) AND (meancount2 GT 0) then box_mean=median(mode_temp[indexformean]) else box_mean=0
            
          endif else box_mean=medianVal
          mode_median[ti]=box_mean
          
        endfor
        ;****************end of calc median mode block
        
        
        ;****************run through with a fixed mode
        For obs_i=0, obsid_count-1 do begin
        
          ;Fill the arrays with the original gain (bp subtracted) and the polyfit
          gain_arr=*cal_array[obs_i].gain[pol_i]
          gain_arr_fit=polyfit_gains[*,*,pol_i,obs_i]
          
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
            
            if keyword_set(compare_ratio) then begin
              amp_use=abs(mode_fit)/nf_use
              amp_use=amp_use/compare_ratio[tile_i,pol_i]
              
              phase_use=atan(mode_fit,/phase)
              
              real_use=amp_use*cos(phase_use)
              imag_use=amp_use*sin(phase_use)
            endif
            
          ;(*cal_return[obs_i].mode_params[pol_i,tile_i])=[mode_i,real_use,imag_use]
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
            if cal_return[obs_i].mode_params[pol_i,ti] NE !Null then real_temp[obs_i]=(*cal_return[obs_i].mode_params[pol_i,ti])[1] else real_temp[obs_i]=0
            if cal_return[obs_i].mode_params[pol_i,ti] NE !Null then imag_temp[obs_i]=(*cal_return[obs_i].mode_params[pol_i,ti])[2] else imag_temp[obs_i]=0
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
            If (meancount1 GT 0) AND (meancount2 GT 0) then box_mean=median(real_temp[indexformean]) else box_mean=0
            
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
            If (meancount1 GT 0) AND (meancount2 GT 0) then box_mean=median(imag_temp[indexformean]) else box_mean=0
            
          endif else box_mean=medianVal
          imag_mean[ti]=box_mean
          
          ;Calculate the amp and phase from the mean of the real and imaginary parts of the fit
          amp_mean[ti]=sqrt(real_mean[ti]^2.+imag_mean[ti]^2.)
          phase_mean[ti]=atan(real_mean[ti]+i_comp*imag_mean[ti],/phase)
        endfor
        ;****************end of calc mean amp and phase block
        
        for obs_i=0, obsid_count-1 do begin
          gain_arr=*cal_array[obs_i].gain[pol_i]
          gain_arr_fit=polyfit_gains[*,*,pol_i,obs_i]
          gain_arr-=gain_arr_fit
          
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
    
  ENDIF else begin
    for obs_i=0, obsid_count-1 do begin
      for pol_i=0,1 do begin
        *cal_return[obs_i].gain[pol_i]=polyfit_gains[*,*,pol_i,obs_i]
      endfor
    endfor
  endelse
  
  undefine_fhd,gain_residual
  RETURN,cal_return
END
