pro poly_comparisons_bftemp_polybylongrun, longrun=longrun, lst_compare=lst_compare, raw=raw, mean_tile=mean_tile, quick_check=quick_check, line_plots=line_plots, spread=spread, $
    make_obs_to_query=make_obs_to_query, dir_name=dir_name,restore_for_mode=restore_for_mode
  ;Program for plotting polyfit by pointings solutions over the longrun to check for stability. Longrun must be set at this time.
  ;lst_compare is not implemented at this time. Raw underplots the raw (minus bandpass by cable by pointing) for each obsid of the
  ;pointing or day and pointing depending on the next keyword.  by_pointing plots each all the pointings over all the days
  ;for a specific tile.  If that is not set, then plots are made for each pointing for each day.
  ;quick_check reads in only the -2 pointing to make plots quicker
  ;line_plots makes the plots described above
  ;spread is a new option to plot spreads of the zeroth order amp params by obs , overplotted with by pointing
  bypointingdiff=1
  cal_savefile = '/nfs/eor-00/h1/nbarry/Aug23_longrunpoly_onephase/'
  ;days of the long run to use
  ;Oct 15 bad
  ;Oct 04 phase is bad for tile 27
  day=['Aug23','Aug27','Sep02','Sep04','Sep06','Sep09','Sep11','Sep13','Sep17','Sep19','Sep30','Oct02','Oct04','Oct08','Oct10','Oct23'];,'Oct25','Oct29']
  
  day_num=N_elements(day)
  rgbcolors=[[69,190,207],[144,23,3],[240,87,249],[45,165,30],[42,25,77],[1,53,8],[167,138,28],[251,193,236],[52,86,193],[246,229,194],[148,33,145],[253,88,89],[239,233,110],$
    [14,128,63],[96,50,14],[176,247,108],[74,126,157],[47,14,36],[212,120,250],[202,29,168]]
    
  ;Directory to print plots -- CAN CHANGE
  outdir='/nfs/eor-00/h1/nbarry/longrun_poly_bftemp_polybylongrun/'
  
  pointing_num=[-2,-1,0,1,2,3]
  ;pointing_num=[3]
  ;pointing_num=[-5,-4,-3,-2,-1,0,1,2,3,4]
  pointing_name=['-2','-1','0','1','2','3']
  ;pointing_name=['3']
  ;pointing_name=['-5','4','-3','-2','-1','0','1','2','3','4']
  If keyword_set(quick_check) then pointing_num=[-2]
  
  parsednames=STRARR(day_num,N_elements(pointing_num))
  parsednames_Aug23=STRARR(day_num,N_elements(pointing_num))
  for day_i=0,day_num-1 do begin
    If keyword_set(longrun) AND ~keyword_set(quick_check) then parsednames[day_i,*]=[day[day_i]+'_minustwo',day[day_i]+'_minusone',day[day_i]+'_zenith',day[day_i]+'_plusone',day[day_i]+'_plustwo',day[day_i]+'_plusthree']
    ;If keyword_set(longrun) AND ~keyword_set(quick_check) then parsednames[day_i,*]=[day[day_i]+'_plusthree']
    If keyword_set(quick_check) then parsednames[day_i,*]=[day[day_i]+'_minustwo']
    parsednames_Aug23[day_i,*]=[day[day_i]+'minustwo',day[day_i]+'minusone',day[day_i]+'zenith',day[day_i]+'plusone',day[day_i]+'plustwo',day[day_i]+'plusthree']
  ;parsednames[day_i,*]=[day[day_i]+'_minusfive',day[day_i]+'_minusfour',day[day_i]+'_minusthree',day[day_i]+'_minustwo',$
  ;  day[day_i]+'_minusone',day[day_i]+'_zenith',day[day_i]+'_plusone',day[day_i]+'_plustwo',day[day_i]+'_plusthree',day[day_i]+'_plusfour']
  endfor
  
  If ~keyword_set(restore_for_mode) then begin
  
    restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/1061316296_obs.sav'
    tile_names=ULONG((*obs.baseline_info).tile_names)
    
    ;***************Loop to get the obsids of each day/pointing of the longrun
    obsid_day_pointing=PTRARR((size(pointing_num))[1],day_num,/allocate)
    obsid_count=INTARR((size(pointing_num))[1])
    beginning_obsid_count=INTARR((size(pointing_num))[1])
    bftemp_pointing=FLTARR((size(pointing_num))[1],day_num)
    
    if keyword_set(make_obs_to_query) then begin
      GET_LUN,lun
      OPENW,lun,outdir+'obs_to_query_v2.txt'
    endif
    
    FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
      undefine, obsid
      undefine, beginning_obsid
      
      FOR day_i=0,day_num-1 DO BEGIN
      
        If keyword_set(longrun) then filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/' + parsednames[day_i,j] + '.txt'
        ;If parsednames[day_i,j] EQ 'Aug23_plusthree' then filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23plusthree.txt' ;TEMP to get Aug23 +3
        If day_i EQ 0 then filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/' + parsednames_Aug23[day_i,j] + '.txt' ;TEMP to get Aug23 +3
        
        obs_array='empty'
        readcol, filename, obs_array, format='A', /silent
        
        If obs_array[0] NE 'empty' Then begin
          textfast,obs_array,/read,file_path=filename,string=1
          
          If keyword_set(obsid) then obsid=[obsid,ULONG(obs_array[*])] else obsid=ULONG(obs_array[*])
          If keyword_set(beginning_obsid) then beginning_obsid=[beginning_obsid,ULONG(obs_array[0])] else beginning_obsid=ULONG(obs_array[0])
          *obsid_day_pointing[j, day_i]=(obs_array[*])
          if keyword_set(make_obs_to_query) then PRINTF, lun,(*obsid_day_pointing[j,day_i])[0]
          
        endif
        
      ENDFOR ; end day for
      
      obsid_count[j]=N_elements(obsid)
      beginning_obsid_count[j]=N_elements(beginning_obsid)
      
    endfor ;end pointing for
    ;***************End of loop to get the obsids of each day/pointing of the longrun
    
    if keyword_set(make_obs_to_query) then begin
      FREE_LUN,lun
      stop
    endif
    
    ;****************************Read in the polyfit calibration solutions, and make the raw gains minus bp sol if specified
    gain=complex(DBLARR(2,day_num,384,128,(size(pointing_num))[1],20))
    
    ;scalefactor=(DBLARR(2,day_num,128,(size(pointing_num))[1],20))
    scalefactorpre=(DBLARR(2,day_num,128,(size(pointing_num))[1],20))
    scalefactorpost=(DBLARR(2,day_num,128,(size(pointing_num))[1],20))
    
    If keyword_set(raw) then gain_raw=complex(DBLARR(2,day_num,384,128,(size(pointing_num))[1],20))
    
    for j=0,(size(pointing_num))[1]-1 do begin
      print, 'Reading in pointing ' + pointing_name[j]
      for day_i=0,day_num-1 do begin
      
        If *obsid_day_pointing[j, day_i] NE !NULL then begin
        
          for obs_i=0, N_elements(*obsid_day_pointing[j, day_i])-1 do begin
          
          
            if day[day_i] NE 'Aug23' then begin
              filename='/nfs/eor-00/h1/nbarry/longrun_std_test_twopolyquad_nomode_onephase/'+day[day_i]+'/'+(*obsid_day_pointing[j, day_i])[obs_i]+'_cal.sav'
              restore, filename
            endif else begin
              ;TEMP to get Aug23 +3
              ;if day[day_i] EQ 'Aug23' then begin
              filename='/nfs/eor-00/h1/nbarry/Aug23_std_test_twopolyquad_polyonly/'+(*obsid_day_pointing[j, day_i])[obs_i]+'_cal.sav'
              restore, filename
            endelse
            
            gain[0,day_i,*,*,j,obs_i]=*cal.gain[0]
            gain[1,day_i,*,*,j,obs_i]=*cal.gain[1]
            
            for pol_j=0,1 do begin
              for tile_i=0,127 do begin
                gain_temp=reform(gain[pol_j,day_i,0:255,tile_i,j,obs_i])
                avestep='prepostsep'
                
                ;Potentially doing nothing, but check for flags nevertheless
                ;ind_full=where(abs(gain_temp) NE 1)
                
                ;scale the amplitude
                ;scalefactor[pol_j,day_i,tile_i,j,obs_i]= mean(abs(gain_temp[ind_full]))
                
                ;scale the amplitude, pre
                scalefactorpre[pol_j,day_i,tile_i,j,obs_i]= mean(abs(reform(gain[pol_j,day_i,0:255,tile_i,j,obs_i])))
                
                ;scale the amplitude, post
                scalefactorpost[pol_j,day_i,tile_i,j,obs_i]= mean(abs(reform(gain[pol_j,day_i,256:383,tile_i,j,obs_i])))
              endfor
            endfor
          endfor
          
          If keyword_set(raw) then begin
            for obs_i=0, N_elements(*obsid_day_pointing[j, day_i])-1 do begin
            
              if day[day_i] NE 'Aug23' then begin
                filename= '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/'+strtrim((*obsid_day_pointing[j, day_i])[obs_i],2)+'_cal.sav'
                restore,filename
                
              endif else begin
                filename='/nfs/eor-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/calibration/'+(*obsid_day_pointing[j, day_i])[obs_i]+'_cal.sav'
                restore, filename
              endelse
              
              If ~keyword_set(spread) AND keyword_set(line_plots) then begin
              
                if day[day_i] NE 'Aug23' then begin
                  ;Reinstate the raw gains from the longrun
                  filename= '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/'+strtrim((*obsid_day_pointing[j, day_i])[obs_i],2)+'_obs.sav'
                  restore,filename
                  
                endif else begin
                  filename='/nfs/eor-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/'+(*obsid_day_pointing[j, day_i])[obs_i]+'_obs.sav'
                  restore, filename
                endelse
                
                for pol_i=0,1 do *cal.gain[pol_i]=*cal.gain[pol_i]+*cal.gain_residual[pol_i]
                
                ;Take the bandpass out of the raw gains. Just added the longrun bp to see effects (smoother freq dependence is what I predict)
                cal_bandpass=vis_cal_bandpass(cal,obs,cal_remainder=cal_remainder,cable_bandpass_fit=1, saved_run_bp=1)
                for pol_i=0,1 do *cal.gain[pol_i]=*cal_remainder.gain[pol_i]
                
                gain_raw[0,day_i,*,*,j,obs_i]=*cal.gain[0]  ;pol x day x freq x tile x pointing x obs
                gain_raw[1,day_i,*,*,j,obs_i]=*cal.gain[1]
                
              endif
              
            endfor
          endif
          
        endif else begin
          gain[0,day_i,*,*,j,*]=1
          gain[1,day_i,*,*,j,*]=1
        endelse
        
      endfor
    endfor ;end pointing for
    ;save, gain_obs_amp, filename=outdir+'gain_obs_amp.sav'
    ;save, gain_poi_amp, filename=outdir+'gain_poi_amp.sav'
    
    flagged_ind=where(abs(gain_raw) EQ 0,flagged_count)
    
    for freq_i=0,383 do begin
    
      ;Scale the amplitudes by the calculated scale factor
      ;gain_raw[*,*,freq_i,*,*,*]=abs(gain_raw[*,*,freq_i,*,*,*])/(scalefactor)*exp(Complex(0,1)*atan(gain_raw[*,*,freq_i,*,*,*],/phase))
    
      ;Scale the amplitudes by the calculated scale factor
      If freq_i LT 256 then gain_raw[*,*,freq_i,*,*,*]=abs(gain_raw[*,*,freq_i,*,*,*])/(scalefactorpre)*exp(Complex(0,1)*atan(gain_raw[*,*,freq_i,*,*,*],/phase))
      
      ;Scale the amplitudes by the calculated scale factor
      If freq_i GT 255 then gain_raw[*,*,freq_i,*,*,*]=abs(gain_raw[*,*,freq_i,*,*,*])/(scalefactorpost)*exp(Complex(0,1)*atan(gain_raw[*,*,freq_i,*,*,*],/phase))
      
    ;scale the poly amplitudes by the calculated scale factor just for fun
    ;gain[*,*,freq_i,*]=abs(gain[*,*,freq_i,*])/(scalefactor)*exp(Complex(0,1)*atan(gain[*,*,freq_i,*],/phase))
      
    endfor
    
    
    
    ;**********Section for determining two quadratic polyfits for firstpass cal
    
    freq_use=where((*obs.baseline_info).freq_use,nf_use)
    
    ;obs_elements_temp=where(abs(gain_raw[0,0,0,0,0,*]) NE 0,obs_elements)
    ;summed_elements=FLOAT(obs_elements*(size(gain_raw))[2])
    ;gain_raw_sub=gain_raw[*,*,freq_use,*,*,*]
    gain_raw_ave=complex(FLTARR(2,384,128,(size(pointing_num))[1]))
    summed_elements=FLTARR(2,384,128,(size(pointing_num))[1])
    gain_fit=FLTARR(2,day_num,384,128,(size(pointing_num))[1])
    gain_fit_Aug23=FLTARR(384,128,2,94);polyfit_gains is freq x tile x pol x obs
    phase_fit_Aug23=FLTARR(384,128,2,94);polyfit_gains is freq x tile x pol x obs
    unscaled_gain_fit=FLTARR(2,day_num,384,128,(size(pointing_num))[1])
    modeadd_gain_fit=complex(FLTARR(2,day_num,384,128,(size(pointing_num))[1]))
    phase_fit=FLTARR(2,day_num,384,128,(size(pointing_num))[1])
    
    
    FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
      filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/' + parsednames[0,j] + '.txt'
      obs_temp='empty'
      readcol, filename, obs_temp, format='A', /silent
      If j NE 0 then parsednumbers=[parsednumbers,N_elements(obs_temp)] else parsednumbers=N_elements(obs_temp)
      If j NE 0 then parsednumbers_tot=[parsednumbers_tot,N_elements(obs_temp)+parsednumbers_tot[j-1]] else parsednumbers_tot=N_elements(obs_temp)
      *obs_ptr[j]=obs_temp
    ENDFOR
    
    ;for pol_i=0, 1 do begin
    ;for tile_i=0,127 do begin
    for j=0,(size(pointing_num))[1]-1 do begin
      for tile_i=0,127 do begin
        for pol_i=0, 1 do begin
        
          gain_raw_temp_pre=reform(abs(gain_raw[pol_i,0,freq_use[0:223],tile_i,j,0]))
          x_arr_pre=freq_use[0:223]
          gain_raw_temp_post=reform(abs(gain_raw[pol_i,0,freq_use[224:335],tile_i,j,0]))
          x_arr_post=freq_use[224:335]
          fit_params_phase=FLTARR(2,day_num)
          
          ;To get the modefit through the ave for just Aug23 for testing          ;gain_raw: pol x day x freq x tile x pointing x obs
          If j EQ 0 then gain_raw_ave[*,tile_i,pol_i,0:parsednumbers_tot[j]-1]=mean(reform(abs(gain_raw[pol_i,0,*,tile_i,j,*])),dimension=2)*mean(atan((gain_raw[pol_i,0,freq_use[0:223],tile_i,j,*]),/phase),dimension=2) else $
            gain_raw_ave[*,tile_i,pol_i,parsednumbers_tot[j-1]-1:parsednumbers_tot[j]-1]=mean(reform(abs(gain_raw[pol_i,0,*,tile_i,j,*])),dimension=2)*mean(atan((gain_raw[pol_i,0,freq_use[0:223],tile_i,j,*]),/phase),dimension=2)

          ;End of modefit through the ave
            
            
          for day_i=0,day_num-1 do begin
          
            x_arr_phase=freq_use
            gain_raw_temp_phase=reform(atan(gain_raw[pol_i,day_i,freq_use,tile_i,j,0],/phase))
            
            for obs_i=0,19 do begin
            
              IF (obs_i NE 0) AND (day_i NE 0) then begin
                gain_raw_temp_pre=[gain_raw_temp_pre,reform(abs(gain_raw[pol_i,day_i,freq_use[0:223],tile_i,j,obs_i]))]
                x_arr_pre=[x_arr_pre,freq_use[0:223]]
                
                gain_raw_temp_post=[gain_raw_temp_post,reform(abs(gain_raw[pol_i,day_i,freq_use[224:335],tile_i,j,obs_i]))]
                x_arr_post=[x_arr_post,freq_use[224:335]]
              endif
              IF (obs_i NE 0) then begin
                gain_raw_temp_phase=[gain_raw_temp_phase,reform(atan(gain_raw[pol_i,day_i,freq_use,tile_i,j,obs_i],/phase))]
                x_arr_phase=[x_arr_phase,freq_use]
              endif
              
            endfor
            
            ;phase fitting
            fit_ind=where(FINITE(gain_raw_temp_phase),count)
            
            If count NE 0 then begin
              fit_params_phase[*,day_i]=poly_fit(x_arr_phase[fit_ind],gain_raw_temp_phase[fit_ind],1) ;gain needs to be a 1D array of a tile's gain for the obs of the pointing for all freq. freq_use needs to be a repeated 1D array to match gain
            endif else begin
              fit_params_phase[*,day_i]=[0,0]
            endelse
            
            
          endfor
          
          fit_ind=where(FINITE(gain_raw_temp_pre),count)
          
          If count NE 0 then begin
            fit_params_pre=poly_fit(x_arr_pre[fit_ind],gain_raw_temp_pre[fit_ind],2) ;gain needs to be a 1D array of a tile's gain for the obs of the pointing for all freq. freq_use needs to be a repeated 1D array to match gain
            
            fit_ind=where(FINITE(gain_raw_temp_post),count)
            fit_params_post=poly_fit(x_arr_post[fit_ind],gain_raw_temp_post[fit_ind],2) ;gain needs to be a 1D array of a tile's gain for the obs of the pointing for all freq. freq_use needs to be a repeated 1D array to match gain
          endif else begin
            fit_params_post=[0,0,0]
            fit_params_pre=[0,0,0]
          endelse
          
          
          ;cal_return.amp_params[pol_i,tile_i]=Ptr_new(fit_params)
          n_freq=384
          
          gain_fit_post=fltarr(n_freq)
          FOR di=0L,2 DO gain_fit_post+=fit_params_post[di]*findgen(n_freq)^di
          gain_fit_pre=fltarr(n_freq)
          FOR di=0L,2 DO gain_fit_pre+=fit_params_pre[di]*findgen(n_freq)^di
          
          gain_fit_phase=fltarr(n_freq,day_num)
          for day_j=0, day_num-1 do begin
            FOR di=0L,1 DO gain_fit_phase[*,day_j]+=fit_params_phase[di,day_j]*findgen(n_freq)^di
          endfor
          
          for day_i=0,day_num-1 do begin
            gain_fit[pol_i,day_i,*,tile_i,j]=[gain_fit_pre[0:255]*scalefactorpre[pol_i,day_i,tile_i,j,0],gain_fit_post[256:383]*scalefactorpost[pol_i,day_i,tile_i,j,0]]
            If day_i EQ 0 then begin
              If j EQ 0 then begin
                gain_fit_Aug23[*,tile_i,pol_i,0:parsednumbers_tot[j]-1]=[gain_fit_pre[0:255],gain_fit_post[256:383]]
                phase_fit_Aug23[*,tile_i,pol_i,0:parsednumbers_tot[j]-1]=gain_fit_phase[*,day_i]
              endif else begin
                gain_fit_Aug23[*,tile_i,pol_i,parsednumbers_tot[j-1]-1:parsednumbers_tot[j]-1]=[gain_fit_pre[0:255],gain_fit_post[256:383]]
                phase_fit_Aug23[*,tile_i,pol_i,parsednumbers_tot[j-1]-1:parsednumbers_tot[j]-1]=gain_fit_phase[*,day_i]
              endelse
            endif
            unscaled_gain_fit[pol_i,day_i,*,tile_i,j]=[gain_fit_pre[0:255],gain_fit_post[256:383]]
            phase_fit[pol_i,day_i,*,tile_i,j]=gain_fit_phase[*,day_i]
          endfor
          
          
        endfor ;end pol for
        
        
      endfor
    endfor
    
    
    ;polyfit_gains is freq x tile x pol x obs
    
    
    cal_polyfit_total=vis_cal_modefit_only_pointing(gain_raw_ave-gain_fit_Aug23*exp(Complex(0,1)*phase_fit_Aug23),cal_array_input,obs_array,parsednumbers[j],pointing_num[j],degree=1,phase_degree=1,$
      file_path=dir_name+'/',cal_cable_reflection_mode_fit=1,cal_cable_reflection_both=1)
      
      
    ;Recent addition to save the resulting gain fits for reading in (after combining a mode) with mkbp_pointings_v2
    for obs_i=0, N_elements(*obsid_day_pointing[j, 0])-1 do begin
      ;filename= '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/'+strtrim((*obsid_day_pointing[j, 0])[obs_i],2)+'_cal.sav'
      filename= '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/calibration/'+strtrim((*obsid_day_pointing[j, 0])[obs_i],2)+'_cal.sav'
      restore,filename
      (*cal.gain[0])[*,*]=gain_fit[0,0,*,*,j]*exp(Complex(0,1)*phase_fit[0,0,*,*,j])
      (*cal.gain[1])[*,*]=gain_fit[1,0,*,*,j]*exp(Complex(0,1)*phase_fit[1,0,*,*,j])
      save, cal, filename=cal_savefile+strtrim((*obsid_day_pointing[j, 0])[obs_i],2)+'_cal.sav'
    endfor
    
    
    save,gain_fit, filename='/nfs/eor-00/h1/nbarry/longrun_polyfit_v2.sav' ;pol x day x freq x tile x pointing
    
    
    ;**********End section for determining two quadratic polyfits for firstpass cal
    stop
    
    
    
    
    
    
    
    
    
    freq_use=where((*obs.baseline_info).freq_use,nf_use)
    
    ;obs_elements_temp=where(abs(gain_raw[0,0,0,0,0,*]) NE 0,obs_elements)
    ;summed_elements=FLOAT(obs_elements*(size(gain_raw))[2])
    ;gain_raw_sub=gain_raw[*,*,freq_use,*,*,*]
    gain_raw_ave=complex(FLTARR(2,384,128,(size(pointing_num))[1]))
    summed_elements=FLTARR(2,384,128,(size(pointing_num))[1])
    
    for pol_i=0, 1 do begin
      for freq_i=0, 383 do begin
        for tile_i=0,127 do begin
          for j=0,(size(pointing_num))[1]-1 do begin
          
            gain_raw_temp=gain_raw[pol_i,0,freq_i,tile_i,j,*]
            for day_i=1,day_num-1 do begin
              gain_raw_temp=[gain_raw_temp,gain_raw[pol_i,day_i,freq_i,tile_i,j,*]]
            endfor
            
            resistant_mean, abs(gain_raw_temp),2,res_mean_real
            ;resistant_mean, imaginary(gain_raw_temp),2,res_mean_imag
            
            gain_raw_ave[pol_i,freq_i,tile_i,j]=res_mean_real;*scalefactorpre[pol_i,day_i,tile_i,j,0];+Complex(0,1)*res_mean_imag
            
          ;gain_raw_ave[pol_i,freq_i,tile_i,j]=total(real_part(gain_raw[pol_i,*,freq_i,tile_i,j,*]),/NAN)+Complex(0,1)*total(imaginary(gain_raw[pol_i,*,freq_i,tile_i,j,*]),/NAN);*(FLOAT((size(gain_raw))[2])/summed_elements)
            
          ;obs_elements_temp=where(FINITE(abs(gain_raw[pol_i,*,freq_i,tile_i,j,*])),obs_elements)
          ;summed_elements[pol_i,freq_i,tile_i,j]=FLOAT(obs_elements)
            
          ;+ mean(real_part(gain_raw),dimension=6,/NAN)*(FLOAT(obs_elements)/summed_elements) + $
          ;  Complex(0,1)*(mean(imaginary(gain_raw),dimension=2,/NAN)*(FLOAT((size(gain_raw))[2])/summed_elements) + mean(imaginary(gain_raw),dimension=6,/NAN)*(FLOAT(obs_elements)/summed_elements))
          endfor
        endfor
      endfor
    endfor
    
    save,gain_raw_ave, filename='/nfs/eor-00/h1/nbarry/longrun_polyfit_average.sav' ;pol x freq x tile x pointing
    
    ;polyfit_input=complex(FLTARR(384,128,2,1))
    
    for pol_i=0, 1 do begin
      for tile_i=0,127 do begin
        for j=0,(size(pointing_num))[1]-1 do begin
        
          for day_i=0,day_num-1 do begin
            ;gain_fit[pol_i,day_i,*,tile_i,j]=[gain_fit_pre[0:255]*scalefactorpre[pol_i,day_i,tile_i,j,0],gain_fit_post[256:383]*scalefactorpost[pol_i,day_i,tile_i,j,0]]
            modeadd_gain_fit[pol_i,day_i,*,tile_i,j]=[reform((unscaled_gain_fit[pol_i,day_i,0:255,tile_i,j]*gain_raw_ave[pol_i,0:255,tile_i,j])*scalefactorpre[pol_i,day_i,tile_i,j,0]),$
              reform((unscaled_gain_fit[pol_i,day_i,256:383,tile_i,j]*gain_raw_ave[pol_i,256:383,tile_i,j])*scalefactorpost[pol_i,day_i,tile_i,j,0])]
            modeadd_gain_fit[pol_i,day_i,*,tile_i,j]=modeadd_gain_fit[pol_i,day_i,*,tile_i,j]*exp(Complex(0,1)*phase_fit[pol_i,day_i,*,tile_i,j])
            
          endfor
        ;polyfit_input[*,tile_i,pol_i,0:94,j]=gain_fit[pol_i,0,*,tile_i,j]
          
        endfor
      endfor
    endfor
    
    parsednames=['Aug23minustwo','Aug23minusone','Aug23zenith','Aug23plusone','Aug23plustwo','Aug23plusthree']
    
    ;obs ptr is where obs of a pointing are stored. It can be made bigger than needed to account for the pointings of the longrun
    obs_ptr=PTRARR(10,/allocate_heap)
    
    ;Different pointing ranges given if it was the longrun or not
    pointing_num=[-2,-1,0,1,2,3]
    
    ;For each pointing, get the obs from the correct text file and read it to an obs_ptr array. Fill variable with 'empty' string first to create a tag if it is unfilled by
    ;the text file. Save the number of obs in that pointing.
    FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
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
      
      polyfit_input=complex(FLTARR(384,128,2,parsednumbers[j]))
      
      ;****Restore loop, per obs in a pointing
      For i=0, parsednumbers[j]-1 DO BEGIN
        obsid=(*obs_ptr[j])[i]
        
        ;Restore cal structure from longrun or regular run
        restore, '/nfs/eor-00/h1/nbarry/Aug23_longrunpoly_afterpolfix/polyonly/' + obsid + '_cal.sav'
        restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/' + obsid + '_params.sav'
        restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/' + obsid + '_obs.sav'
        
        ;Make an array of structures.  First make a new cal structure and fill it with the value restored.  There are various points in the init code that set the below values to
        ;zero, so they have to be refilled individually.
        cal2=fhd_struct_init_cal(obs, params, n_pol=cal.n_pol,n_freq=cal.n_freq,n_tile=cal.n_tile,freq=cal.freq,gain=cal.gain,gain_residual=cal.gain_residual, $
          n_vis_cal=cal.n_vis_cal, amp_params=cal.amp_params, phase_params=cal.phase_params,mode_params=cal.mode_params)
        ; cal2=fhd_struct_init_cal_noparams_July2015(obs, n_pol=cal.n_pol,n_freq=cal.n_freq,n_tile=cal.n_tile,freq=cal.freq,gain=cal.gain,gain_residual=cal.gain_residual, $
        ;   n_vis_cal=cal.n_vis_cal, amp_params=cal.amp_params, phase_params=cal.phase_params,mode_params=cal.mode_params)
        cal2.mode_params = cal.mode_params
        cal2.amp_params=cal.amp_params
        cal2.phase_params=cal.phase_params
        cal2.gain_residual=cal.gain_residual
        
        ;Make and array of structures from the new cal structure made just previously
        If (i eq 0) Then cal_array = replicate(cal2,parsednumbers[j]) ELSE cal_array[i]=cal2
        
        If (i eq 0) Then n_pol = intarr(parsednumbers[j])
        n_pol[i]=cal_array[i].n_pol
        
        polyfit_input[*,*,0,i]=*cal_array[i].gain[0]
        polyfit_input[*,*,1,i]=*cal_array[i].gain[1]
        
        ;*cal_array[i].gain[0]=reform(modeadd_gain_fit[0,0,*,*,j])
        ;*cal_array[i].gain[1]=reform(modeadd_gain_fit[1,0,*,*,j])
        
        ;Setup the obs structure array and fill it on successive loops
        If (i eq 0) Then obs_array = replicate(obs,parsednumbers[j]) ELSE obs_array[i]=obs
        
      ENDfor
      
      If ~keyword_set(quick_check) then begin
        For i=0, parsednumbers[j]-1 DO BEGIN
          *cal_array[i].gain[0]=reform(modeadd_gain_fit[0,0,*,*,j])
          *cal_array[i].gain[1]=reform(modeadd_gain_fit[1,0,*,*,j])
        endfor
      endif else begin
      
        *cal_array[i].gain[0]=reform(modeadd_gain_fit[0,0,*,*])
        *cal_array[i].gain[1]=reform(modeadd_gain_fit[1,0,*,*])
        
      endelse
      
      save,polyfit_input, FILENAME=dir_name+'/' + 'polyfit_input_'+pointing_name[j]+'.sav'
      save,cal_array, FILENAME=dir_name+'/' + 'cal_array_'+pointing_name[j]+'.sav'
      save,obs_array, FILENAME=dir_name+'/' + 'obs_array_'+pointing_name[j]+'.sav'
      
      cal_polyfit=vis_cal_modefit_only_pointing(polyfit_input,cal_array,obs_array,parsednumbers[j],file_path=dir_name+'/',cal_cable_reflection_both=1,cal_cable_reflection_mode_fit=1,cal_cable_reflection_fit=150)
      
      For i=0, parsednumbers[j]-1 DO BEGIN
        obsid=(*obs_ptr[j])[i]
        cal=cal_polyfit[i]
        save, cal, FILENAME=dir_name+'/'+obsid+'_cal.sav
      endfor
      
    ENDFOR
  ;*************End of restore loop
    
  endif else begin
  
    parsednames=['Aug23minustwo','Aug23minusone','Aug23zenith','Aug23plusone','Aug23plustwo','Aug23plusthree']
    
    ;obs ptr is where obs of a pointing are stored. It can be made bigger than needed to account for the pointings of the longrun
    obs_ptr=PTRARR(10,/allocate_heap)
    
    ;Different pointing ranges given if it was the longrun or not
    pointing_num=[-2,-1,0,1,2,3]
    
    ;For each pointing, get the obs from the correct text file and read it to an obs_ptr array. Fill variable with 'empty' string first to create a tag if it is unfilled by
    ;the text file. Save the number of obs in that pointing.
    FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
      filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/' + parsednames[j] + '.txt'
      obs_temp='empty'
      readcol, filename, obs_temp, format='A', /silent
      If j NE 0 then parsednumbers=[parsednumbers,N_elements(obs_temp)] else parsednumbers=N_elements(obs_temp)
      *obs_ptr[j]=obs_temp
    ENDFOR
    ;*************end of setup
    
    FOR j=0, (size(pointing_num))[1]-1 DO BEGIN
    
      restore,dir_name+'/' + 'polyfit_input_'+pointing_name[j]+'.sav'
      restore,dir_name+'/' + 'cal_array_'+pointing_name[j]+'.sav'
      restore,dir_name+'/' + 'obs_array_'+pointing_name[j]+'.sav'
      
      
      cal_polyfit=vis_cal_modefit_only_pointing(polyfit_input,cal_array,obs_array,parsednumbers[j],file_path=dir_name+'/',cal_cable_reflection_both=1,cal_cable_reflection_mode_fit=1,cal_cable_reflection_fit=150)
      
      
      For i=0, parsednumbers[j]-1 DO BEGIN
        obsid=(*obs_ptr[j])[i]
        cal=cal_polyfit[i]
        save, cal, FILENAME=dir_name+'/'+obsid+'_cal.sav
      endfor
      
    endfor
    
  endelse
  
  
  stop
  
  ;gain_raw_ave=gain_raw_ave/summed_elements
  
  ;stop
  
  ;linfit(atan(gain_raw[0,0,freq_use,tile_i,0,0],/phase)
  
  color_array=['black','blue','green','purple','red','brown']
  
  for tile_i=0,127 do begin
  
    cgPS_Open,'/nfs/eor-00/h1/nbarry/scaled_poly_ave/scaled_poly_resave/gain_mean_plusallave_'+avestep+'_tile'+strtrim(tile_names[tile_i],2)+'.png',/quiet,/nomatch
    ;cgplot, abs(gain_raw[0,0,freq_use,tile_i,0,0]), yrange=[.8,1.2],xrange=[0,336], title='Scaled longrun raw gains by '+avestep+ ', tile '+strtrim(tile_names[tile_i],2)+' xx -2pointing',$
    ;  xtitle='frequency index', ytitle='scaled gain', charsize=1,/NODATA
    
    cgplot, abs(gain_raw[0,0,freq_use,tile_i,0,0]), yrange=[.8,1.2],xrange=[0,336], title='Scaled longrun mean gain amplitudes (pre/post sep), tile '+strtrim(tile_names[tile_i],2)+' xx',$
      xtitle='frequency index', ytitle='scaled mean gain amplitude', charsize=1,/NODATA
    ;cgplot, atan(gain_raw[0,0,freq_use,tile_i,0,0],/phase), yrange=[-3.14,3.14],xrange=[0,336], title='Longrun raw phase, tile '+strtrim(tile_names[tile_i],2)+' xx -2pointing', xtitle='frequency index', ytitle='phase', charsize=1
    ;Device, Decomposed=0
    ;for day_i=0,(size(gain_raw))[2]-1 do begin
    ;  for obs_i=0,(size(gain_raw))[6]-1 do begin
    ;cgoplot,abs(gain_raw[0,day_i,freq_use,tile_i,0,obs_i])
    ;TVLCT, rgbcolors[0,day_i], rgbcolors[1,day_i], rgbcolors[2,day_i],100
    ;    cgoplot,abs(gain_raw[0,day_i,freq_use,tile_i,0,obs_i]);,/phase),color=100B
      
      
    ;  endfor
    ;endfor
    ;Device, Decomposed=1
    for j=0,(size(pointing_num))[1]-1 do begin
      cgoplot, abs(gain_raw_ave[0,freq_use,tile_i,j]),color=color_array[j]
    endfor
    
    gain_raw_ave_p=mean(abs(gain_raw_ave[0,freq_use,tile_i,*]),dimension=4)
    cgoplot, abs(gain_raw_ave_p), color='brown',thick=4
    
    
    cgLegend, Title=[pointing_name,'Mean'], color=color_array, $
      Location=[0.7,0.8],charsize=1,VSpace=.9,thick=3
      
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    
  endfor
  
  ;cgPS_Open,'/nfs/eor-00/h1/nbarry/scaled_gain_raw_ave_'+avestep+'_tile13.png',/quiet,/nomatch
  ;cgplot, abs(gain_raw[0,0,freq_use,2,0,0]), yrange=[1,2],xrange=[0,336], title='Scaled longrun raw gains by '+avestep+ ', tile 13 xx -2pointing', xtitle='frequency index', ytitle='scaled gain', charsize=1
  ;for day_i=0,(size(gain_raw))[2]-1 do begin
  ;  for obs_i=0,(size(gain_raw))[6]-1 do begin
  ;    cgoplot,abs(gain_raw[0,day_i,freq_use,2,0,obs_i])
  ;  endfor
  ;endfor
  ;cgoplot, abs(gain_raw_ave[0,freq_use,2,0]),color='green'
  ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  
  stop
  ;****************************End of read in the polyfit calibration solutions, and make the raw gains minus bp sol if specified
  
  
  
  ;****************************Read in and parse temperature data
  
  ;Find the smaller temperature data file with just the right obs to get temperature data.
  ;To get this file, run this code with make_obs_to_query set and then run query.sh
  filename='/nfs/eor-00/h1/nbarry/longrun_poly_bftemp/obs_queried_v2.txt'
  
  ;Read out temperature data in the form of the string due to the funky format
  textfast,data_array,/read,file_path=filename,string=1
  
  temperature_array=FLTARR(day_num,128,8) ;Want day x tile x pointing
  
  tile_temp=STRARR(8)
  
  for j=0,(size(pointing_num))[1]-1 do begin
    for day_i=0,day_num-1 do begin
    
      If *obsid_day_pointing[j, day_i] NE !NULL then begin
      
        temp_index=where(strmatch(data_array,'*'+(*obsid_day_pointing[j, day_i])[0]+'*') EQ 1)
        split_temp_array=STRARR(N_elements(temp_index),3)
        for i=0,N_elements(temp_index)-1 do begin
          split_temp_array[i,*]=strsplit(data_array[temp_index[i]], '|',/EXTRACT)  ;Split line in data file into receiver, obsid full, and temps per tile
          split_temp_array[i,1]=(strsplit(split_temp_array[i,1], '|',/EXTRACT))[0]  ;Remove extra ms from end of obsid
          split_temp_array[i,0]=STRING(ULONG(split_temp_array[i,0])*10)  ;Add zero to end of receiver for getting tile names easier
          split_temp_array[i,2]=strmid(split_temp_array[i,2],1) ;Remove '{' from beginning of temp data
          str_length=strlen(split_temp_array[i,2]) ;Find length of string for next command
          split_temp_array[i,2]=strmid(split_temp_array[i,2],0,str_length-1) ;Remove '}' from end of temp data
          tile_temp[*]=strsplit(split_temp_array[i,2], ',',/EXTRACT) ;Split up temp data into 8 discrete temperatures per receiver
          
          
          
          for tile_i=0,7 do begin
            tile_index=where(tile_names EQ (ULONG(split_temp_array[i,0])+tile_i+1))
            tile_temp_float=Convert_To_Type(tile_temp[tile_i],4)
            temperature_array[day_i,tile_index,j]=tile_temp_float ;Add per tile data into array as a float
          endfor
        endfor
        
      endif
      
    endfor
    
  endfor
  
  ;****************************End of read in and parse temperature data
  
  ;GET_LUN,lun
  ;OPENW,lun,outdir+'linear_fits_days_fixed.txt',width=600
  ;header=['pointing','tile_name','pol','Group_1_days','Group_1_y-inter','Group_1_slope','Group_1_chi^2','Group_2_days','Group_2_y-inter','Group_2_slope','Group_2_chi^2']
  
  ;Using preexisting file to extract information about which tiles have which cable length
  mode_filepath=filepath('mwa_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  
  ;Taking tile information and cross-matching it with the nonflagged tiles array, resulting in nonflagged tile arrays
  ;grouped by cable length
  cable_length_ref=cable_len[Uniq(cable_len,Sort(cable_len))]
  n_cable=N_Elements(cable_length_ref)
  tile_use_arr=Ptrarr(n_cable)
  FOR cable_i=0,n_cable-1 DO tile_use_arr[cable_i]=Ptr_new(where((cable_len EQ cable_length_ref[cable_i])))
  
  slopes_per_cable=PTRARR(n_cable,/allocate)
  slopes_per_cable1=PTRARR(n_cable,/allocate)
  slopes_per_cable2=PTRARR(n_cable,/allocate)
  slope_diff=PTRARR(4,/allocate)
  
  textfast, data_array,read=1,string=1,file_path=outdir+'linear_fits.txt',delimiter=' '
  
  
  for pol_i=0,1 do begin
    if pol_i EQ 0 then pol_name='xx'
    if pol_i EQ 1 then pol_name='yy'
    
    
    poi_name=['-2','-1','0','1','2','3']
    If keyword_set(line_plots) then begin
    
      for tile_i=0,127 do begin
        input_str=STRARR(11,((size(pointing_num))[1]))
        
        savelocation=outdir+strtrim(string(tile_names[tile_i]),2)+'_scatter_'+pol_name+'.png'
        ;cgPS_Open,savelocation,/quiet,/nomatch
        title='Tile ' + strtrim(string(tile_names[tile_i]),2)+', gain temperature slope scatter, '+pol_name
        
        ;cgplot, [1,2,3],[1,2,3],xtitle='Average Gain (after bandpass removal)',ytitle='Beamformer Temperature (C)', $
        ;  title=title,charsize=1, psym=2, symsize=0.5,yrange=[0,50],xrange=[1,2],/NODATA
        
        ;for j=0,((size(pointing_num))[1])-1 do begin
        for j=2,2 do begin
          ;Get inputs which are not flagged
          input=mean(abs(gain[pol_i,*,0:255,tile_i,j]),dimension=3)
          input_index=where(input NE 1)
          ;result = LINFIT(reform(input[input_index]),reform(temperature_array[input_index,tile_i,j]),chisqr=chisqr)
          
          ;Sort unflagged inputs in order of coldest to hotest in anticipation that they will be grouped on temperature,
          ;though it is no requirement
          temperature_temp=reform(temperature_array[input_index,tile_i,j])
          sorted_index=sort(temperature_temp)
          temperature_sorted=temperature_temp[sorted_index]
          input_sorted=reform(input[input_index[sorted_index]])
          
          day_temp=day[input_index]
          day_sorted=day_temp[sorted_index]
          
          if pol_i EQ 0 then pol_con=0
          if pol_i EQ 1 then pol_con=640
          current_file_index=pol_con+1+5*tile_i+j
          
          file_days1=data_array[3,current_file_index]
          file_days2=data_array[7,current_file_index]
          file_days1_split=strsplit(file_days1,',',/EXTRACT)
          file_days2_split=strsplit(file_days2,',',/EXTRACT)
          
          If keyword_set(bycable) OR keyword_set(bycablebypol) then begin
            For group_i=5,9,4 do begin
              For cable_i=0,n_cable-1 do begin
                temp=where(tile_i EQ *tile_use_arr[cable_i],in_cable_group)
                
                If in_cable_group GT 0 then begin
                
                  IF abs(FLOAT(data_array[group_i,current_file_index])) GT 0 then begin
                  
                    If slopes_per_cable[cable_i] NE !NULL then *slopes_per_cable[cable_i]=[*slopes_per_cable[cable_i],FLOAT(data_array[group_i,current_file_index])] $
                    else *slopes_per_cable[cable_i]=FLOAT(data_array[group_i,current_file_index])
                    
                  endif
                endif
              endfor
            endfor
          endif
          
          If keyword_set(bypointingdiff) then begin
          
            IF abs(FLOAT(data_array[4,current_file_index])) GT 0 then begin
            
              diff_minusone=(FLOAT(data_array[4,current_file_index]))-(FLOAT(data_array[4,current_file_index-1]))
              diff_minustwo=(FLOAT(data_array[4,current_file_index]))-(FLOAT(data_array[4,current_file_index-2]))
              diff_plusone=(FLOAT(data_array[4,current_file_index]))-(FLOAT(data_array[4,current_file_index+1]))
              diff_plustwo=(FLOAT(data_array[4,current_file_index]))-(FLOAT(data_array[4,current_file_index+2]))
              temp_slopes=[diff_minustwo,diff_minusone,diff_plusone,diff_plustwo]
              
              For poi_i=0,3 do begin
                If slope_diff[poi_i] NE !NULL then *slope_diff[poi_i]=[*slope_diff[poi_i],temp_slopes[poi_i]] else *slope_diff[poi_i]=temp_slopes[poi_i]
              endfor
              
            endif
          endif
          
          
          
          ;For group_i=5,9,4 do begin
          For cable_i=0,n_cable-1 do begin
            temp=where(tile_i EQ *tile_use_arr[cable_i],in_cable_group)
            
            If in_cable_group GT 0 then begin
            
              group_i=5
              IF abs(FLOAT(data_array[group_i,current_file_index])) GT 0 then begin
                If slopes_per_cable1[cable_i] NE !NULL then *slopes_per_cable1[cable_i]=[*slopes_per_cable1[cable_i],FLOAT(data_array[group_i,current_file_index])] $
                else *slopes_per_cable1[cable_i]=FLOAT(data_array[group_i,current_file_index])
              endif
              
              group_i=9
              IF abs(FLOAT(data_array[group_i,current_file_index])) GT 0 then begin
                If slopes_per_cable2[cable_i] NE !NULL then *slopes_per_cable2[cable_i]=[*slopes_per_cable2[cable_i],FLOAT(data_array[group_i,current_file_index])] $
                else *slopes_per_cable2[cable_i]=FLOAT(data_array[group_i,current_file_index])
              endif
              
            endif
          endfor
        ;endfor
          
          
          
        ;input_str[0:2,j]=data_array[0:2,current_file_index]
        ;input_str[4:6,j]=data_array[4:6,current_file_index]
        ;input_str[8:10,j]=data_array[8:10,current_file_index]
          
        ;undefine, indexmap1, indexmap2
        ;for day_j=0, N_elements(file_days1_split)-1 do begin
        ;  index_temp=where(strmatch(day,file_days1_split[day_j]) EQ 1, n_count)
        ;  if keyword_set(indexmap1) then indexmap1=[indexmap1,index_temp] else indexmap1=index_temp
        ;
        ;  If file_days2_split[0] NE '0' then begin
        ;    If day_j LT N_elements(file_days2_split) then begin
        ;    index_temp=where(strmatch(day,file_days2_split[day_j]) EQ 1, n_count)
        ;    if keyword_set(indexmap2) then indexmap2=[indexmap2,index_temp] else indexmap2=index_temp
          
        ;    endif
        ;  endif
          
        ;endfor
          
        ;group1days=day_sorted[indexmap1]
        ;group1days_in=strjoin(group1days,',')
        ;If file_days2_split[0] NE '0' then begin
        ;  group2days=day_sorted[indexmap2]
        ;  group2days_in=strjoin(group2days,',')
        ;endif else group2days_in=0
          
        ;input_str[3,j]=group1days_in
        ;input_str[7,j]=group2days_in
          
        endfor
        
      ;printf,lun, input_str
      endfor
      
      
    endif ;end of line plots
    
    
  endfor
  ;FREE_LUN,lun
  
  If keyword_set(bycable) then begin
  
    y_arr_fill_ptr=PTRARR(n_cable,/allocate)
    x_arr_fill_ptr=PTRARR(n_cable,/allocate)
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/slope_hist_bycable_filled.png',/quiet,/nomatch
    color_array=['black','red','blue','purple','brown','orange']
    names=['90m','150m','230m','320m','400m','524m']
    thick_array=[2,2,2,2,2,2]
    cgplot, [1,2],[1,2],xrange=[-150,-30], yrange=[0,80],/NODATA, title='Slopes of bf temps vs. average gains per day, grouped by cable type', $
      xtitle='Fit slope value', ytitle='Density', charsize=1
      
    ;temp=plot( [1,2],[1,2],xrange=[-150,-40], yrange=[0,100],/NODATA)
      
    For cable_ii=0,n_cable-1 do begin
      IF cable_ii EQ 0 then cable_i=5
      IF cable_ii EQ 1 then cable_i=0
      IF cable_ii EQ 2 then cable_i=1
      IF cable_ii EQ 3 then cable_i=2
      IF cable_ii EQ 4 then cable_i=4
      IF cable_ii EQ 5 then cable_i=3
      
      ;Histogram slopes by cable length
      ;
      ;The binsize is arbitrary, but it does give a pretty plot. Could be increased. Must be even.
      binsize=4.
      
      ;Histogram the uv baselines to the binsize provided above, and record that in the reverse indices for applying
      ;to the visibility difference later. Convert the uv from light speed time to wavelengths.
      result=histogram(*slopes_per_cable[cable_i],binsize=binsize,reverse_indices=ri, /NAN, locations=locations, omax=omax)
      
      ;Calculate the visibility sigma (normalized due to real/imag) for each bin using reverse indices
      ;slopes_per_bin=FLTARR(N_elements(result))
      ;for i=0, N_elements(result)-1 do if result[i] GT 0 then slopes_per_bin[i]=(*slopes_per_cable[cable_i])[ri[ri[i]:ri[i+1]-1]]
      
      ;Make the x input and y input pretty for the plot
      y_arr=[result[0], result, result[N_elements(result)-1]]
      x_arr=[locations[0],locations+binsize/2,omax]
      
      undefine, y_arr_fill, x_arr_fill
      for i=0,N_elements(x_arr)-1 do begin
        IF i EQ 0 then begin
          y_arr_fill=0
          x_arr_fill=x_arr[i]-binsize/2
        endif
        y_arr_fill=[y_arr_fill,y_arr[i],y_arr[i]]
        x_arr_fill=[x_arr_fill,x_arr[i]-binsize/2,x_arr[i]+binsize/2]
        IF i EQ N_elements(x_arr)-1 then begin
          y_arr_fill=[y_arr_fill,0]
          x_arr_fill=[x_arr_fill,x_arr[i]+binsize/2]
        endif
      endfor
      *y_arr_fill_ptr[cable_i]=y_arr_fill
      *x_arr_fill_ptr[cable_i]=x_arr_fill
      
      ;Make
      ;cgPS_Open,'/nfs/eor-00/h1/nbarry/freq'+strtrim(freq,2)+'vissigma_inter_time'+strtrim(time_i,2)+'.png',/quiet,/nomatch
      cgoplot, x_arr, y_arr, psym=10, color=color_array[cable_i], thick=thick_array[cable_i]
      
    endfor
    stop
    For cable_ii=0,n_cable-1 do begin
      IF cable_ii EQ 0 then cable_i=5
      IF cable_ii EQ 1 then cable_i=0
      IF cable_ii EQ 2 then cable_i=1
      IF cable_ii EQ 3 then cable_i=2
      IF cable_ii EQ 4 then cable_i=4
      IF cable_ii EQ 5 then cable_i=3
      cgcolorfill, *x_arr_fill_ptr[cable_i], *y_arr_fill_ptr[cable_i], color=color_array[cable_i]
    endfor
    
    al_legend, ['90m','150m','230m','320m','400m','524m'],linestyle=[0,0,0,0,0,0],color=color_array, charsize=1
    
    
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif
  
  
  If keyword_set(bypointingdiff) then begin
  
    y_arr_fill_ptr=PTRARR(5,/allocate)
    x_arr_fill_ptr=PTRARR(5,/allocate)
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/slope_histdiff_bypointing.png',/quiet,/nomatch
    color_array=['black','red','blue','purple'];,'brown','orange']
    names=['Zen-MinusTwo','Zen-MinusOne','Zen-PlusOne','Zen-PlusTwo']
    thick_array=[3,3,3,3]
    cgplot, [1,2],[1,2],xrange=[-50,50], yrange=[0,90],/NODATA, title='Slopes of bf temps vs. ave gains per day, differenced from zenith', $
      xtitle='Fit slope value difference from zenith', ytitle='Density', charsize=1
      
    ;temp=plot( [1,2],[1,2],xrange=[-150,-40], yrange=[0,100],/NODATA)
      
    For poi_i=0,3 do begin
      ;IF cable_ii EQ 0 then cable_i=5
      ;IF cable_ii EQ 1 then cable_i=0
      ;IF cable_ii EQ 2 then cable_i=1
      ;IF cable_ii EQ 3 then cable_i=2
      ;IF cable_ii EQ 4 then cable_i=4
      ;IF cable_ii EQ 5 then cable_i=3
    
      ;Histogram slopes by cable length
      ;
      ;The binsize is arbitrary, but it does give a pretty plot. Could be increased. Must be even.
      binsize=3.
      
      ;Histogram the uv baselines to the binsize provided above, and record that in the reverse indices for applying
      ;to the visibility difference later. Convert the uv from light speed time to wavelengths.
      result=histogram(*slope_diff[poi_i],binsize=binsize,reverse_indices=ri, /NAN, locations=locations, omax=omax)
      
      ;Calculate the visibility sigma (normalized due to real/imag) for each bin using reverse indices
      ;slopes_per_bin=FLTARR(N_elements(result))
      ;for i=0, N_elements(result)-1 do if result[i] GT 0 then slopes_per_bin[i]=(*slopes_per_cable[cable_i])[ri[ri[i]:ri[i+1]-1]]
      
      ;Make the x input and y input pretty for the plot
      y_arr=[0, result, 0]
      x_arr=[locations[0],locations+binsize/2,omax]
      
      undefine, y_arr_fill, x_arr_fill
      for i=0,N_elements(x_arr)-1 do begin
        IF i EQ 0 then begin
          y_arr_fill=0
          x_arr_fill=x_arr[i]-binsize/2
        endif
        y_arr_fill=[y_arr_fill,y_arr[i],y_arr[i]]
        x_arr_fill=[x_arr_fill,x_arr[i]-binsize/2,x_arr[i]+binsize/2]
        IF i EQ N_elements(x_arr)-1 then begin
          y_arr_fill=[y_arr_fill,0]
          x_arr_fill=[x_arr_fill,x_arr[i]+binsize/2]
        endif
      endfor
      *y_arr_fill_ptr[poi_i]=y_arr_fill
      *x_arr_fill_ptr[poi_i]=x_arr_fill
      
      ;Make
      ;cgPS_Open,'/nfs/eor-00/h1/nbarry/freq'+strtrim(freq,2)+'vissigma_inter_time'+strtrim(time_i,2)+'.png',/quiet,/nomatch
      cgoplot, x_arr, y_arr, psym=10, color=color_array[poi_i], thick=thick_array[poi_i]
      
    endfor
    
    For poi_i=0,3 do begin
    ;IF cable_ii EQ 0 then cable_i=5
    ;IF cable_ii EQ 1 then cable_i=0
    ;IF cable_ii EQ 2 then cable_i=1
    ;IF cable_ii EQ 3 then cable_i=2
    ;IF cable_ii EQ 4 then cable_i=4
    ;IF cable_ii EQ 5 then cable_i=3
    ;cgcolorfill, *x_arr_fill_ptr[poi_i], *y_arr_fill_ptr[poi_i], color=color_array[poi_i]
    endfor
    
    al_legend, names,linestyle=[0,0,0,0],color=color_array, charsize=1
    
    
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif
  
  
  if keyword_set(bypolbycable) then begin
  
    y_arr_fill_ptr=PTRARR(n_cable,/allocate)
    x_arr_fill_ptr=PTRARR(n_cable,/allocate)
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/slope_hist_bycable_'+pol_name+'.png',/quiet,/nomatch
    color_array=['black','red','blue','purple','brown','orange']
    names=['90m','150m','230m','320m','400m','524m']
    thick_array=[2,2,2,2,2,2]
    cgplot, [1,2],[1,2],xrange=[-150,-30], yrange=[0,80],/NODATA, title='Slopes of bf temps vs. gains per day, grouped by cable type, '+pol_name, $
      xtitle='Fit slope value', ytitle='Density', charsize=1
      
    For cable_ii=0,n_cable-1 do begin
      IF cable_ii EQ 0 then cable_i=5
      IF cable_ii EQ 1 then cable_i=0
      IF cable_ii EQ 2 then cable_i=1
      IF cable_ii EQ 3 then cable_i=2
      IF cable_ii EQ 4 then cable_i=4
      IF cable_ii EQ 5 then cable_i=3
      
      ;The binsize is arbitrary, but it does give a pretty plot. Could be increased. Must be even.
      binsize=4.
      
      ;Histogram the uv baselines to the binsize provided above, and record that in the reverse indices for applying
      ;to the visibility difference later. Convert the uv from light speed time to wavelengths.
      result=histogram(*slopes_per_cable[cable_i],binsize=binsize,reverse_indices=ri, /NAN, locations=locations, omax=omax)
      
      ;Make the x input and y input pretty for the plot
      y_arr=[0, result, 0]
      x_arr=[locations[0],locations+binsize/2,omax]
      
      undefine, y_arr_fill, x_arr_fill
      for i=0,N_elements(x_arr)-1 do begin
        IF i EQ 0 then begin
          y_arr_fill=0
          x_arr_fill=x_arr[i]-binsize/2
        endif
        y_arr_fill=[y_arr_fill,y_arr[i],y_arr[i]]
        x_arr_fill=[x_arr_fill,x_arr[i]-binsize/2,x_arr[i]+binsize/2]
        IF i EQ N_elements(x_arr)-1 then begin
          y_arr_fill=[y_arr_fill,0]
          x_arr_fill=[x_arr_fill,x_arr[i]+binsize/2]
        endif
      endfor
      *y_arr_fill_ptr[cable_i]=y_arr_fill
      *x_arr_fill_ptr[cable_i]=x_arr_fill
      
      cgoplot, x_arr, y_arr, psym=10, color=color_array[cable_i], thick=thick_array[cable_i]
      
    endfor
    
    For cable_ii=0,n_cable-1 do begin
      IF cable_ii EQ 0 then cable_i=5
      IF cable_ii EQ 1 then cable_i=0
      IF cable_ii EQ 2 then cable_i=1
      IF cable_ii EQ 3 then cable_i=2
      IF cable_ii EQ 4 then cable_i=4
      IF cable_ii EQ 5 then cable_i=3
    ;cgcolorfill, *x_arr_fill_ptr[cable_i], *y_arr_fill_ptr[cable_i], color=color_array[cable_i]
    endfor
    
    al_legend, ['90m','150m','230m','320m','400m','524m'],linestyle=[0,0,0,0,0,0],color=color_array, charsize=1
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif
  
  
  
  stop
  
  y_arr_fill_ptr=PTRARR(n_cable,/allocate)
  x_arr_fill_ptr=PTRARR(n_cable,/allocate)
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/slope_hist_bycable_bygroup.png',/quiet,/nomatch
  color_array=['black','red','blue','purple','brown','orange']
  names=['90m','150m','230m','320m','400m','524m']
  thick_array=[2,2,2,2,2,2]
  cgplot, [1,2],[1,2],xrange=[-150,-30], yrange=[0,80],/NODATA, title='Slopes of bf temps vs. gains per day, grouped by cable type', $
    xtitle='Fit slope value', ytitle='Density', charsize=1
    
  For group_i=5,9,4 do begin
    If group_i EQ 5 then slopes_per_cable=slopes_per_cable1
    If group_i EQ 9 then slopes_per_cable=slopes_per_cable2
    For cable_ii=0,n_cable-1 do begin
      IF cable_ii EQ 0 then cable_i=5
      IF cable_ii EQ 1 then cable_i=0
      IF cable_ii EQ 2 then cable_i=1
      IF cable_ii EQ 3 then cable_i=2
      IF cable_ii EQ 4 then cable_i=4
      IF cable_ii EQ 5 then cable_i=3
      
      ;The binsize is arbitrary, but it does give a pretty plot. Could be increased. Must be even.
      binsize=4.
      
      ;Histogram the uv baselines to the binsize provided above, and record that in the reverse indices for applying
      ;to the visibility difference later. Convert the uv from light speed time to wavelengths.
      result=histogram(*slopes_per_cable[cable_i],binsize=binsize,reverse_indices=ri, /NAN, locations=locations, omax=omax)
      
      ;Make the x input and y input pretty for the plot
      y_arr=[0, result, 0]
      x_arr=[locations[0],locations+binsize/2,omax]
      
      undefine, y_arr_fill, x_arr_fill
      for i=0,N_elements(x_arr)-1 do begin
        IF i EQ 0 then begin
          y_arr_fill=0
          x_arr_fill=x_arr[i]-binsize/2
        endif
        y_arr_fill=[y_arr_fill,y_arr[i],y_arr[i]]
        x_arr_fill=[x_arr_fill,x_arr[i]-binsize/2,x_arr[i]+binsize/2]
        IF i EQ N_elements(x_arr)-1 then begin
          y_arr_fill=[y_arr_fill,0]
          x_arr_fill=[x_arr_fill,x_arr[i]+binsize/2]
        endif
      endfor
      *y_arr_fill_ptr[cable_i]=y_arr_fill
      *x_arr_fill_ptr[cable_i]=x_arr_fill
      
      cgoplot, x_arr, y_arr, psym=10, color=color_array[cable_i], thick=thick_array[cable_i]
      If group_i EQ 9 then cgcolorfill, *x_arr_fill_ptr[cable_i], *y_arr_fill_ptr[cable_i], color=color_array[cable_i]
      
    endfor
  endfor
  
  For cable_ii=0,n_cable-1 do begin
    IF cable_ii EQ 0 then cable_i=5
    IF cable_ii EQ 1 then cable_i=0
    IF cable_ii EQ 2 then cable_i=1
    IF cable_ii EQ 3 then cable_i=2
    IF cable_ii EQ 4 then cable_i=4
    IF cable_ii EQ 5 then cable_i=3
  ;cgcolorfill, *x_arr_fill_ptr[cable_i], *y_arr_fill_ptr[cable_i], color=color_array[cable_i]
  endfor
  
  al_legend, ['90m','150m','230m','320m','400m','524m'],linestyle=[0,0,0,0,0,0],color=color_array, charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  stop
end
