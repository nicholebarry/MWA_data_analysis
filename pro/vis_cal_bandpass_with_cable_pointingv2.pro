FUNCTION vis_cal_bandpass_with_cable_pointingv2,cal_array,obs_array,obsid_count,pointing_num,cal_remainder=cal_remainder,file_path_fhd=file_path_fhd,advanced_plotting=advanced_plotting

  ;split_jump=1
  ;longrun=1
  if keyword_set(longrun) then begin
    for obs_i=0, obsid_count-1 do begin
    
      ;save original frequencies to use
      original_freq_use=(*obs_array[obs_i].baseline_info).freq_use
      
      If ~keyword_set(split_jump) then begin
        ;Only a bandpass within the specified frequencies are calculated, so a only pre dig jump looks truncated
        cal_bandpass=vis_cal_bandpass(cal_array[obs_i],obs_array[obs_i],cal_remainder=cal_remainder,cable_bandpass_fit=1,_Extra=extra)
        cal_polyfit=vis_cal_polyfit(cal_remainder,obs_array[obs_i],degree=2,_Extra=extra)
        cal_array[obs_i].amp_params=cal_polyfit.amp_params
        cal_array[obs_i].phase_params=cal_polyfit.phase_params
        cal_array[obs_i].mode_params=cal_polyfit.mode_params
      endif else begin
        ;Only a bandpass within the specified frequencies are calculated, so a only pre dig jump looks truncated
        cal_bandpass=vis_cal_bandpass(cal_array[obs_i],obs_array[obs_i],cal_remainder=cal_remainder,cable_bandpass_fit=1,_Extra=extra)
        
        ;Fit parameters are only calculated in the specified frequencies, but the final fit is over all frequencies.
        ;Pre dig jump calcs should be fine, even though the product looks like it's done over all freq
        ;Just fitting linear trend, bc don't want to overfit and it will be interesting to see if degree=2 was just trying to capture dig jump
        (*obs_array[obs_i].baseline_info).freq_use[256:383]=0
        cal_polyfit_pre=vis_cal_polyfit(cal_remainder,obs_array[obs_i],degree=1,cal_cable_reflection_fit=150,_Extra=extra)
        (*obs_array[obs_i].baseline_info).freq_use=original_freq_use
        (*obs_array[obs_i].baseline_info).freq_use[0:255]=0
        cal_polyfit_post=vis_cal_polyfit(cal_remainder,obs_array[obs_i],degree=1,cal_cable_reflection_fit=150,_Extra=extra)
        (*obs_array[obs_i].baseline_info).freq_use=original_freq_use
        
        ;Have two sets of params now to use
        amp_params_pre=cal_polyfit_pre.amp_params
        phase_params_pre=cal_polyfit_pre.phase_params
        mode_params_pre=cal_polyfit_pre.mode_params
        amp_params_post=cal_polyfit_post.amp_params
        phase_params_post=cal_polyfit_post.phase_params
        mode_params_post=cal_polyfit_post.mode_params
        
        freq_arr=cal_array[obs_i].freq
        
        If obs_i EQ 0 then poly_input_gains=complex(FLTARR(384,128,obsid_count,2))
        
        ;freq x tiles x obs x pol
        poly_input_gains[256:383,*,obs_i,0]=((*cal_polyfit_post.gain[0])[256:383,*])
        poly_input_gains[0:255,*,obs_i,0]=((*cal_polyfit_pre.gain[0])[0:255,*])
        poly_input_gains[256:383,*,obs_i,1]=((*cal_polyfit_post.gain[1])[256:383,*])
        poly_input_gains[0:255,*,obs_i,1]=((*cal_polyfit_pre.gain[1])[0:255,*])
        
        
        ;plot_split_jump=1
        if keyword_set(plot_split_jump) then begin
          tile_names=(*obs_array[0].baseline_info).tile_names
          obsname=obs_array[obs_i].obsname
          filename='/nfs/eor-00/h1/nbarry/Aug23_nodigjump/poly_plots/'+obsname+'_poly.png'
          tile_use=(*obs_array[obs_i].baseline_info).tile_use
          tile_exist=(*obs_array[obs_i].baseline_info).tile_use
          
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
          
          plot_cals_sub_quickanddirty_ploy,freq_arr,poly_input_gains_xx,poly_input_gains_yy,filename=filename,phase=phase,real_vs_imaginary=real_vs_imaginary,$
            tile_A=tile_A,tile_B=tile_B,tile_use=tile_use,tile_exist=tile_exist,tile_names=tile_names,yrange=yrange,$
            obsname=obsname,plot_pos=plot_pos,cal_plot_charsize=cal_plot_charsize,cal_plot_symsize=cal_plot_symsize,cal_plot_resize=cal_plot_resize
        endif
        
      endelse
    endfor
  endif

  n_cable=6
  
  ;Get the cable lengths by accessing the cable reflection txt file.  This could change, so be aware.
  mode_filepath=filepath(obs_array[0].instrument+'_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  
  ;Initialize arrays for excluding flagged tiles
  tile_use_cable=PTRARR(n_cable,obsid_count, /allocate)
  nt_use_cable=INTARR(n_cable,obsid_count)
  t_use=PTRARR(obsid_count,/allocate)
  nt_use=INTARR(obsid_count)
  
  ;Setup loop. Define arrays dependent on obs number of the pointing. Different obs ids have different unflagged tiles.
  FOR obs_i=0,obsid_count-1 DO BEGIN
  
    gain_arr_ptr=cal_array[obs_i].gain
    n_pol=cal_array[obs_i].n_pol
    n_freq=cal_array[obs_i].n_freq
    n_tile=cal_array[obs_i].n_tile
    IF N_Elements(obs_array[obs_i]) GT 0 THEN freq_use=where((*obs_array[obs_i].baseline_info).freq_use) ELSE freq_use=lindgen(n_freq)
    
    ;only pre dig jump
    ;freq_use=((*obs_array[0].baseline_info).freq_use)
    
    nf_use=N_Elements(freq_use)
    freq_arr=cal_array[obs_i].freq
    nt_use[obs_i]=N_Elements(where((*obs_array[obs_i].baseline_info).tile_use EQ 1))
    *t_use[obs_i]=(where((*obs_array[obs_i].baseline_info).tile_use EQ 1))
    
    n_pol=N_Elements(gain_arr_ptr)
    IF obs_i eq 0 THEN gain_arr_ptr_fitted=PTRARR(obsid_count,n_pol,/allocate)
    
    ;gain_arr_ptr2=Ptrarr(n_pol,/allocate)
    ;gain_arr_ptr3=Ptrarr(n_pol,/allocate)
    
    bandpass_arr=(dblarr((n_pol)*n_cable+1,n_freq))
    ;!!!!!!!!!!!!!!!!!!!Modification for splitting up the 524m cable
    ;cable5split=1
    If keyword_set(cable5split) then bandpass_arr=(dblarr((n_pol)*n_cable+3,n_freq))
    ;!!!!!!!!!!!!!!!!!!!Modification for splitting up the 524m cable
    
    ;Setting up the bandpass_arr with a column for frequency
    bandpass_arr[0,*]=freq_arr
    bandpass_col_count=0
    
    ;Removing tile 113 clear up cable2
    ;((*obs_array[obs_i].baseline_info).tile_use[113])=0
    ;Removing anything that looks cruddy in the cable residuals to see if it clears up more, esp cable5
    ;((*obs_array[obs_i].baseline_info).tile_use[124])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[123])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[122])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[112])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[104])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[98])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[97])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[96])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[95])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[88])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[87])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[79])=0
    ;((*obs_array[obs_i].baseline_info).tile_use[78])=0
    
    
    ;Group cable lengths in arrays
    *tile_use_cable[0,obs_i]=where((*obs_array[obs_i].baseline_info).tile_use AND cable_len EQ 90)
    *tile_use_cable[1,obs_i]=where((*obs_array[obs_i].baseline_info).tile_use AND cable_len EQ 150)
    *tile_use_cable[2,obs_i]=where((*obs_array[obs_i].baseline_info).tile_use AND cable_len EQ 230)
    *tile_use_cable[3,obs_i]=where((*obs_array[obs_i].baseline_info).tile_use AND cable_len EQ 320)
    *tile_use_cable[4,obs_i]=where((*obs_array[obs_i].baseline_info).tile_use AND cable_len EQ 400)
    *tile_use_cable[5,obs_i]=where((*obs_array[obs_i].baseline_info).tile_use AND cable_len EQ 524)
    
    ;***************************Polyfit removal block
    ;In order to see if there is a benifit to rm fit/mode fit first, it has been taken out
    degree=2
    phase_degree=1
    gain_arr=complex(DBLARR(n_pol,n_freq,nt_use[obs_i]))
    
    FOR pol_i=0, n_pol-1 DO BEGIN
      FOR tile_i=0,nt_use[obs_i]-1 DO BEGIN
        gain_fit=DBLARR(n_freq)
        phase_fit=DBLARR(n_freq)
        mode_params=DBLARR(3)
        amp_params=*(cal_array[obs_i].amp_params[pol_i,(*t_use[obs_i])[tile_i]])
        phase_params=*(cal_array[obs_i].phase_params[pol_i,(*t_use[obs_i])[tile_i]])
        FOR di=0L,degree DO gain_fit+=amp_params[di]*findgen(n_freq)^di
        FOR di=0L,phase_degree DO phase_fit+=phase_params[di]*findgen(n_freq)^di
        If (cal_array[obs_i].mode_params[pol_i,(*t_use[obs_i])[tile_i]]) NE !NULL THEN mode_params=*(cal_array[obs_i].mode_params[pol_i,(*t_use[obs_i])[tile_i]])
        
        gain_arr[pol_i,*,tile_i]=(gain_fit*Exp(Complex(0,1)*phase_fit)+mode_params[1]*Exp(-Complex(0,1)*2.*!Pi*mode_params[0]*findgen(n_freq)/n_freq+Complex(0,1)*(mode_params[2])))
      ;if mode_params[1] NE 0 then stop
        
      ENDFOR
      
      ;*gain_arr_ptr_fitted[obs_i,pol_i]=((*cal_array[obs_i].gain[pol_i])[*,(*t_use[obs_i])]/gain_arr[pol_i,*,*])
      ;Uncomment below to not remove polyfit before bp calculation
      *gain_arr_ptr_fitted[obs_i,pol_i]=((*cal_array[obs_i].gain[pol_i])[*,(*t_use[obs_i])])
      IF keyword_set(split_jump) then begin
        *gain_arr_ptr_fitted[obs_i,pol_i]=((*cal_array[obs_i].gain[pol_i])[*,*]/poly_input_gains[*,*,obs_i,pol_i])
      endif
      
    ENDFOR
    ;***************************End of polyfit removal block
    
    
    ;Continue to set up, now per cable length type
    FOR cable_i=0,n_cable-1 DO BEGIN
      ;number of tiles to cycle through
      nt_use_cable[cable_i,obs_i]=N_Elements(*tile_use_cable[cable_i,obs_i])
    ENDFOR
  ENDFOR ;endfor for obs
  
  ;Remove unneeded arrays and set up an array for comparisons later
  undefine, gain_arr, gain_fit, phase_fit, mode_params, amp_params, phase_params
  orig_amp_arr = complex(dblarr(obsid_count,N_elements(freq_use),Max(nt_use),n_pol))
  
  
  ;******************************Pointers block
  FOR obs_i=0,obsid_count-1 DO BEGIN
    FOR cable_i=0,n_cable-1 DO BEGIN
      ;Final setup dependence, which is pol per cable per obs.  Populate a gain array and point to it.
      FOR pol_i=0,n_pol-1 DO BEGIN
      
        index_size=size(*tile_use_cable[cable_i,obs_i])
        index=INTARR(index_size[1])
        for ind_i=0, index_size[1]-1 do index[ind_i]=where(*t_use[obs_i] EQ (*tile_use_cable[cable_i,obs_i])[ind_i])
        
        gain=*gain_arr_ptr_fitted[obs_i,pol_i] ;n_freq x n_tile element complex array
        
        ;Only use gains that are associated with the right cables and useable gains. Then get the amplitude and phase
        gain_use=extract_subarray(gain,freq_use,index)
        amp=(gain_use)
        
        ;Setup amp pointer only once and setup amp per freq+tile+pol for every obs and cable
        IF obs_i eq 0 and pol_i eq 0 and cable_i eq 0 THEN amp_array_ptr = PTRARR(n_cable, /allocate)
        IF pol_i eq 0 THEN amp_array = complex(dblARR(obsid_count,N_elements(freq_use),MAX(nt_use_cable[cable_i,*]),n_pol))
        
        ;Make a pointer that will point to an array which includes gain amplitude for all freq of all tiles per pol, which can change in size (hence the pointer)
        For tile_i=0, nt_use_cable[cable_i,obs_i]-1 DO amp_array[obs_i,*,tile_i,pol_i] = amp[*,tile_i]
        
        IF (obs_i eq 0) and (pol_i eq 0) THEN *amp_array_ptr[cable_i] = complex(dblarr(obsid_count,N_elements(freq_use),Max(nt_use_cable[cable_i,*]),n_pol))
        IF pol_i eq 1 THEN *amp_array_ptr[cable_i]+=amp_array
        
        ;Set up arrays of original methods for comparison later
        if cable_i eq 0 then gain_orig_use=extract_subarray(gain,freq_use,findgen(nt_use[obs_i]))
        if cable_i eq 0 then orig_amp_arr[obs_i,*,*t_use[obs_i],pol_i]+=gain_orig_use
        
      ENDFOR ;end per pol
      
    ENDFOR ;end per cable
    
  ENDFOR ;end per obs
  ;******************************End of pointers block
  
  
  undefine, gain_orig_use, amp_array, amp
  
  
  ;******************************Final bp solution calculation block
  ;Label the number of seperate bp solutions (5 for default cable cal)
  IF keyword_set(cable5split) then n_cable_normal=5 else n_cable_normal=6
  
  ;Now for the final loop with the median/mean calculation
  For cable_i=0,n_cable_normal-1 DO BEGIN
    For pol_i=0,n_pol-1 DO BEGIN
    
      amp2=(dblarr(obsid_count,nf_use,Max(nt_use_cable[cable_i,*])))
      
      ;Fill amp2 with the gains of each obs, each freq, per tile for the pol of the for loop divided by the mean of the gains over all obs and freq, but per tile and per pol.
      FOR tile_i=0,Max(nt_use_cable[cable_i,*])-1 DO BEGIN
        ;Filling temporary variable with gains on a per tile basis.  Find indices where gain is nonzero for mean calculation
        tile_amp_temp=(*amp_array_ptr[cable_i])[*,*,tile_i,pol_i]
        indices_tile_amp = where(tile_amp_temp NE 0)
        
        ;Take the resistant mean of the ABSOLUTE value of the gains and with an inputting sigma for which to include in the mean
        resistant_mean, abs(tile_amp_temp[indices_tile_amp]),2,tile_amp_temp_mean

        amp2[*,*,tile_i]=(tile_amp_temp_mean EQ 0) ? 0:(abs((*amp_array_ptr[cable_i])[*,*,tile_i,pol_i])/tile_amp_temp_mean)
      ENDFOR
      
      ;Reinitialize a temporary variable
      bandpass_single = (DBLARR(nf_use))
      
      ;Begin to take the mean over all tiles in a cable type and pointing for a specific frequency to get the solution for that frequency.
      For freq_i=0,nf_use-1 DO begin
        amp3=amp2[*,freq_i,*]
        
        ;Make sure no zeroed amplitudes and included in the resistant mean
        indices = (where(amp3 NE 0))
        amp4=amp3[indices]
        
        ;Take the resistant mean of the absolute value with just values between one sigma.
        resistant_mean, abs(amp4), 2, amp4_mean

        bandpass_single[freq_i]=amp4_mean
      endfor
      bandpass_col_count= bandpass_col_count+1
      bandpass_arr[bandpass_col_count,freq_use]=bandpass_single
    ;If bandpass_col_count EQ 12 THEN stop
      
      
    ENDFOR
  ENDFOR
  ;******************************End of final bp solution calculation block

  
  ;******************************Final bp solution calculation for 524m cable block
  ;Now for the final loop with the median/mean calculation
  If keyword_set(cable5split) then begin
    For cable_i=5,n_cable-1 DO BEGIN
      For pol_i=0,n_pol-1 DO BEGIN
      
        outside_list=[78,79,80,87,88,95,96,97,98,104,105,112,122,123,124]
        ;outside_list_cablearr_index=outside_list
        inside_list=[36,50,75,82,84,85,86,89,92,93,103,107,108,114,115]
        ;inside_list_cablearr_index=inside_list
        
        outside_list_ptr=PTRARR(obsid_count,/allocate)
        inside_list_ptr=PTRARR(obsid_count,/allocate)
        nt_use_outside=INTARR(obsid_count)
        nt_use_inside=INTARR(obsid_count)
        
        cable5tile_index=where(cable_len EQ 524)
        
        for obs_j=0,obsid_count-1 do begin
          for out_i=0, (size(outside_list))[1]-1 do begin
            temp_index=where(outside_list[out_i] EQ *tile_use_cable[cable_i,obs_j])
            if temp_index NE -1 and keyword_set(*outside_list_ptr[obs_j]) then (*outside_list_ptr[obs_j])=[(*outside_list_ptr[obs_j]),temp_index] else if temp_index NE -1 then (*outside_list_ptr[obs_j])=[temp_index]
          ;outside_list_cablearr_index[out_i]=where(outside_list[out_i] EQ FIX((*obs_array[obs_j].baseline_info).tile_names[where(cable_len EQ 524)]))
          ;outside_list_cablearr_index[out_i]=where(outside_list[out_i] EQ *tile_use_cable[5,obs_i])
          endfor
          for in_i=0, (size(inside_list))[1]-1 do begin
            temp_index=where(inside_list[in_i] EQ *tile_use_cable[cable_i,obs_j])
            if temp_index NE -1 and keyword_set(*inside_list_ptr[obs_j]) then (*inside_list_ptr[obs_j])=[(*inside_list_ptr[obs_j]),temp_index] else if temp_index NE -1 then (*inside_list_ptr[obs_j])=[temp_index]
          ;inside_list_cablearr_index[in_i]=where(inside_list[in_i] EQ FIX((*obs_array[obs_j].baseline_info).tile_names[where(cable_len EQ 524)]))
          ;inside_list_cablearr_index[in_i]=where(inside_list[in_i] EQ *tile_use_cable[5,obs_i])
          endfor
          nt_use_outside[obs_j]=N_Elements(*outside_list_ptr[obs_j])
          nt_use_inside[obs_j]=N_elements(*inside_list_ptr[obs_j])
        endfor
        
        outside_list_cablearr_index=*outside_list_ptr[(where(nt_use_outside EQ Max(nt_use_outside[*])))[0]]
        inside_list_cablearr_index=*inside_list_ptr[(where(nt_use_inside EQ Max(nt_use_inside[*])))[0]]
        
        amp2_outside=(fltarr(obsid_count,nf_use,Max(nt_use_outside[*])))
        amp2_inside=(fltarr(obsid_count,nf_use,Max(nt_use_inside[*])))
        
        ;Fill amp2 with the gains of each obs, each freq, per tile for the pol of the for loop divided by the median of the gains over all obs and freq, but per tile and per pol.
        FOR tile_i=0,Max(nt_use_outside[*])-1 DO BEGIN
          tile_amp_temp=(*amp_array_ptr[cable_i])[*,*,outside_list_cablearr_index[tile_i],pol_i]
          indices_tile_amp = where(tile_amp_temp NE 0)
          ;amp2[indices[*,*],indices[1,*],tile_i]=(Median(abs((*amp_array_ptr[cable_i])[indices[0,*],indices[1,*],tile_i,pol_i]) EQ 0) ? 0:(abs((*amp_array_ptr[cable_i])[indices[0,*],indices[1,*],tile_i,pol_i])/Median(abs((*amp_array_ptr[cable_i])[indices[0,*],indices[1,*],tile_i,pol_i]))))
          ;amp2[*,*,tile_i]=(Median(abs(tile_amp_temp[indices_tile_amp])) EQ 0) ? 0:(abs((*amp_array_ptr[cable_i])[*,*,tile_i,pol_i])/Median(abs(tile_amp_temp[indices_tile_amp])))
          resistant_mean, abs(tile_amp_temp[indices_tile_amp]),2,tile_amp_temp_mean
          amp2_outside[*,*,tile_i]=(tile_amp_temp_mean EQ 0) ? 0:(abs((*amp_array_ptr[cable_i])[*,*,outside_list_cablearr_index[tile_i],pol_i])/tile_amp_temp_mean)
        ENDFOR
        
        FOR tile_i=0,Max(nt_use_inside[*])-1 DO BEGIN
          tile_amp_temp=(*amp_array_ptr[cable_i])[*,*,inside_list_cablearr_index[tile_i],pol_i]
          indices_tile_amp = where(tile_amp_temp NE 0)
          ;amp2[indices[*,*],indices[1,*],tile_i]=(Median(abs((*amp_array_ptr[cable_i])[indices[0,*],indices[1,*],tile_i,pol_i]) EQ 0) ? 0:(abs((*amp_array_ptr[cable_i])[indices[0,*],indices[1,*],tile_i,pol_i])/Median(abs((*amp_array_ptr[cable_i])[indices[0,*],indices[1,*],tile_i,pol_i]))))
          ;amp2[*,*,tile_i]=(Median(abs(tile_amp_temp[indices_tile_amp])) EQ 0) ? 0:(abs((*amp_array_ptr[cable_i])[*,*,tile_i,pol_i])/Median(abs(tile_amp_temp[indices_tile_amp])))
          resistant_mean, abs(tile_amp_temp[indices_tile_amp]),2,tile_amp_temp_mean
          amp2_inside[*,*,tile_i]=(tile_amp_temp_mean EQ 0) ? 0:(abs((*amp_array_ptr[cable_i])[*,*,inside_list_cablearr_index[tile_i],pol_i])/tile_amp_temp_mean)
        ENDFOR
        
        bandpass_single_inside = (DBLARR(nf_use))
        bandpass_single_outside = (DBLARR(nf_use))
        
        ;!!Need to make sure potential zeros are not included in the next mean, which could involve them.
        
        For freq_i=0,nf_use-1 DO begin
          amp3=amp2_outside[*,freq_i,*]
          indices = (where(amp3 NE 0))
          amp4=amp3[indices]
          resistant_mean, abs(amp4), 2, amp4_mean
          bandpass_single_outside[freq_i]=amp4_mean
        endfor
        For freq_i=0,nf_use-1 DO begin
          amp3=amp2_inside[*,freq_i,*]
          indices = (where(amp3 NE 0))
          amp4=amp3[indices]
          resistant_mean, abs(amp4), 2, amp4_mean
          bandpass_single_inside[freq_i]=amp4_mean
        endfor
        bandpass_col_count= bandpass_col_count+1
        bandpass_arr[bandpass_col_count,freq_use]=bandpass_single_outside
        bandpass_arr[bandpass_col_count+2,freq_use]=bandpass_single_inside
      ENDFOR
    ENDFOR
  endif
  ;******************************End of final bp solution calculation block
  
  
  
  undefine, amp3
  
  
  ;******************************No cable cal pointing comparison block
  ;This block of code is for making a no cable cal by pointing comparison.  Only relevant when that is the case.
  If keyword_set(no_cable_cal_compare) then begin
    bandpass_arr_orig=(dblarr((n_pol)+1,n_freq))
    bandpass_arr_orig[0,*]=freq_arr
    bandpass_col_count=0
    
    For pol_i=0,n_pol-1 DO BEGIN
    
      amp2=dblarr(obsid_count,nf_use,Max(nt_use))
      
      ;Fill amp2 with the gains of each obs, each freq, per tile for the pol of the for loop divided by the median of the gains over all obs and freq, but per tile and per pol.
      FOR tile_i=0,Max(nt_use)-1 DO BEGIN
        amp2[*,*,tile_i]=(Median(abs(orig_amp_arr[*,*,tile_i,pol_i])) EQ 0) ? 0:(abs(orig_amp_arr[*,*,tile_i,pol_i])/Median(abs(orig_amp_arr[*,*,tile_i,pol_i])))
      ENDFOR
      
      bandpass_single = DBLARR(nf_use)
      
      For freq_i=0,nf_use-1 DO begin
        amp3=amp2[*,freq_i,*]
        indices = (where(amp3 NE 0))
        amp4=amp3[indices]
        bandpass_single[freq_i]=Median(abs(amp4))
      endfor
      
      bandpass_col_count= bandpass_col_count+1
      bandpass_arr_orig[bandpass_col_count,freq_use]=bandpass_single
      
    ENDFOR
  endif
  ;******************************End of no cable cal pointing comparison block
  
  
  ;******************************Taking FFTs of bandpass fit block -- outdated
  IF Keyword_set(make_fft) THEN BEGIN
    image_path='/nfs/eor-00/h1/nbarry/FFT_pointing_gains.png'
    cgPS_Open,image_path,/quiet,/nomatch
    rippleft_xx = ABS(FFT(bandpass_arr[3,*]))
    rippleft_yy = ABS(FFT(bandpass_arr[5,*]))
    cgplot,rippleft_xx,color='blue',yrange=[0,.01],xrange=[0,50],title='FFT of bandpass gains',xtitle='index',ytitle='amplitude'
    cgoplot,rippleft_yy,color='red'
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  ENDIF
  ;******************************End of FFT of bandpass fit block  -- outdated
  
  ;undefine, file_path_fhd
  ;******************************Basic plotting block
  IF Keyword_Set(file_path_fhd) THEN BEGIN
  
    If ~keyword_set(cable5split) then begin
      ;Finding the total number of tiles with a cable type per pointing
      total_use_cable = FLTARR(n_cable)
      FOR cable_i=0, n_cable-1 DO BEGIN
        FOR obs_i=0, obsid_count-1 DO total_use_cable[cable_i] += nt_use_cable[cable_i,obs_i]
      ENDFOR
      
      ;Setting up file paths for saving the image and txt file
      basename=file_basename(file_path_fhd)
      dirpath=file_dirname(file_path_fhd)
      ;image_path=filepath(basename,root=dirpath,sub='output_images')
      image_path=file_path_fhd
      IF file_test(file_dirname(image_path),/directory) EQ 0 THEN file_mkdir,file_dirname(image_path)
      ;export_path=filepath(basename,root=dirpath,sub='calibration')
      export_path=file_path_fhd
      
      ;Export the txt file of the bandpass fit
      IF file_test(file_dirname(export_path),/directory) EQ 0 THEN file_mkdir,file_dirname(export_path)
      Textfast,bandpass_arr,/write,file_path=export_path+strtrim(string(pointing_num),1)+'_bandpass'
      
      ;Setting up plotting variables (range, color, ticks, etc)
      freq=freq_arr/1E6
      xtickv=[ceil(min(freq)/10)*10,floor(max(freq)/10)*10]
      xtickname=strtrim(round(xtickv),2)
      xrange=[min(freq)-(max(freq)-min(freq))/8,max(freq)+(max(freq)-min(freq))/8]
      yrange=[min(bandpass_arr[1:(n_pol*6),freq_use]),max(bandpass_arr[1:(n_pol*6),freq_use])]
      ytickv=[yrange[0],mean(yrange),yrange[1]]
      axiscolor='black'
      
      ;Begin bandpass fit plot for xx
      cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_xx_cc.png',/quiet,/nomatch
      
      ;Plot all of the different cable types
      cgplot,freq[freq_use],bandpass_arr[1,freq_use],color='blue',title=string(pointing_num) + ' xx',xtitle='Frequency [MHz]',ytitle='Gain',$
        yrange=yrange,xrange=xrange,/noerase,axiscolor=axiscolor,psym=2,symsize=0.2
      cgoplot,freq[freq_use],bandpass_arr[3,freq_use],color='red',psym=2,symsize=0.2
      cgoplot,freq[freq_use],bandpass_arr[5,freq_use],color='green',psym=2,symsize=0.2
      cgoplot,freq[freq_use],bandpass_arr[7,freq_use],color='purple',psym=2,symsize=0.2
      cgoplot,freq[freq_use],bandpass_arr[9,freq_use],color='yellow',psym=2,symsize=0.2
      cgoplot,freq[freq_use],bandpass_arr[11,freq_use],color='black',psym=2,symsize=0.2
      cgoplot,freq[freq_use],bandpass_arr[11,freq_use],color='black',psym=2,symsize=0.2
      
      cgLegend, Title=['90m cables ('+ Strtrim(String(UINT(total_use_cable[0])),1) +')','150m cables ('+ Strtrim(String(UINT(total_use_cable[1])),1) +')', $
        '230m cables ('+ Strtrim(String(UINT(total_use_cable[2])),1) +')', '320m cables ('+ Strtrim(String(UINT(total_use_cable[3])),1) +')', $
        '400m cables ('+ Strtrim(String(UINT(total_use_cable[4])),1) +')','524m cables ('+ Strtrim(String(UINT(total_use_cable[5])),1) +')'], $
        Color=['blue','red','green','purple','yellow','black'],Psym=[2,2,2,2,2,2], $
        Length=0.0,Location=[0.55,0.85]
        
      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
      ;Begin bandpass fit plot for yy
      cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_yy_cc.png',/quiet,/nomatch
      
      ;Plot the different cable types
      cgplot,freq[freq_use],bandpass_arr[2,freq_use],color='blue',title=string(pointing_num) + ' yy',xtitle='Frequency [MHz]',ytitle='Gain',$
        yrange=yrange,xrange=xrange,/noerase,axiscolor=axiscolor,psym=2,symsize=0.2
      cgoplot,freq[freq_use],bandpass_arr[4,freq_use],color='red',psym=2,symsize=0.2
      cgoplot,freq[freq_use],bandpass_arr[6,freq_use],color='green',psym=2,symsize=0.2
      cgoplot,freq[freq_use],bandpass_arr[8,freq_use],color='purple',psym=2,symsize=0.2
      cgoplot,freq[freq_use],bandpass_arr[10,freq_use],color='yellow',psym=2,symsize=0.2
      cgoplot,freq[freq_use],bandpass_arr[12,freq_use],color='black',psym=2,symsize=0.2
      cgLegend, Title=['90m cables ('+ Strtrim(String(UINT(total_use_cable[0])),1) +')','150m cables ('+ Strtrim(String(UINT(total_use_cable[1])),1) +')', $
        '230m cables ('+ Strtrim(String(UINT(total_use_cable[2])),1) +')', '320m cables ('+ Strtrim(String(UINT(total_use_cable[3])),1) +')', $
        '400m cables ('+ Strtrim(String(UINT(total_use_cable[4])),1) +')','524m cables ('+ Strtrim(String(UINT(total_use_cable[5])),1) +')'], $
        Color=['blue','red','green','purple','yellow','black'],Psym=[2,2,2,2,2,2], $
        Length=0.0,Location=[0.55,0.85]
        
      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
    endif
    
  ENDIF
  ;***********************************End of basic plotting block
  
  
  ;***********************************Advanced plotting block -- includes comparisons between fits and inputs
  ;advanced_plotting=1
  If keyword_set(advanced_plotting) then begin
  
    If ~keyword_set(cable5split) then begin
      ;Set up ptrs to store data for plotting
      bandpass_arr_change_xx=ptrarr(n_cable, /allocate)
      bandpass_arr_change_yy=ptrarr(n_cable, /allocate)
      
      For cable_i=0, n_cable-1 do begin
      
        ;Create temporary arrays for the bandpass change. Hacky solution is used to account for the fact that some tiles are flagged in an observation,
        ;so the array sizes change. Random value of 12345678 is inputted as an initialization, which will be weeded out during actual plotting and calculations
        ;bandpass_arr_change_xx_temp=fltarr(obsid_count,nf_use,Max(nt_use_cable[cable_i,*]),/NoZero)
        ;bandpass_arr_change_yy_temp=fltarr(obsid_count,nf_use,Max(nt_use_cable[cable_i,*]),/NoZero)
        bandpass_arr_change_xx_temp=fltarr(obsid_count,Max(nt_use_cable[cable_i,*]),/NoZero)
        bandpass_arr_change_yy_temp=fltarr(obsid_count,Max(nt_use_cable[cable_i,*]),/NoZero)
        bandpass_arr_change_xx_temp[*,*,*]=12345678.
        bandpass_arr_change_yy_temp[*,*,*]=12345678.
        
        For obs_i = 0, obsid_count-1 do begin
        
          ;index and index_size account for the fact that the tile array changes per observation due to flagging.
          index_size=size(*tile_use_cable[cable_i,obs_i])
          index=INTARR(index_size[1])
          
          ;Find where the used tiles per observation array holds the current tile that is used for the cable. Save that location in an array for later referencing.
          for ind_i=0, index_size[1]-1 do index[ind_i]=where(*t_use[obs_i] EQ (*tile_use_cable[cable_i,obs_i])[ind_i])
          
          for tile_i=0,nt_use_cable[cable_i,obs_i]-1 do begin
          
            for freq_i=3, 3 do begin;nf_use-1 do begin
            
              ;Take the fitted bandpass per pointing per cable and find the fractional difference of the inputted (unfitted) gains for that cable type.
              ;bandpass_arr_change_xx_temp[obs_i,freq_i,tile_i]=(abs(bandpass_arr[(cable_i*2+1),freq_use[freq_i]])-$
              ;  abs((*gain_arr_ptr_fitted[obs_i,0])[freq_use[freq_i],index[tile_i]])) / abs((*gain_arr_ptr_fitted[obs_i,0])[freq_use[freq_i],index[tile_i]])
              bandpass_arr_change_xx_temp[obs_i,tile_i]=(abs(bandpass_arr[(cable_i*2+1),freq_use[freq_i]])-$
                abs((*gain_arr_ptr_fitted[obs_i,0])[freq_use[freq_i],index[tile_i]])) / abs((*gain_arr_ptr_fitted[obs_i,0])[freq_use[freq_i],index[tile_i]])
              ;bandpass_arr_change_xx_temp[obs_i,freq_i,tile_i]=abs((*gain_arr_ptr_fitted[obs_i,0])[freq_use[freq_i],index[tile_i]])
                
              ;bandpass_arr_change_xx_temp[obs_i,freq_i,tile_i]=(abs((*gain_arr_ptr_fitted[obs_i,0])[freq_use[freq_i],index[tile_i]]))/(abs(bandpass_arr[(cable_i*2+1),freq_use[freq_i]]))
                
              ;bandpass_arr_change_yy_temp[obs_i,freq_i,tile_i]=((abs(bandpass_arr[(cable_i*2+2),freq_use[freq_i]])-$
              ;  abs((*gain_arr_ptr_fitted[obs_i,1])[freq_use[freq_i],index[tile_i]])))/abs((*gain_arr_ptr_fitted[obs_i,1])[freq_use[freq_i],index[tile_i]])
              bandpass_arr_change_yy_temp[obs_i,tile_i]=((abs(bandpass_arr[(cable_i*2+2),freq_use[freq_i]])-$
                abs((*gain_arr_ptr_fitted[obs_i,1])[freq_use[freq_i],index[tile_i]])))/abs((*gain_arr_ptr_fitted[obs_i,1])[freq_use[freq_i],index[tile_i]])
            ;bandpass_arr_change_yy_temp[obs_i,freq_i,tile_i]=abs((*gain_arr_ptr_fitted[obs_i,1])[freq_use[freq_i],index[tile_i]])
                
            ;bandpass_arr_change_yy_temp[obs_i,freq_i,tile_i]=abs((*gain_arr_ptr_fitted[obs_i,1])[freq_use[freq_i],index[tile_i]])/(abs(bandpass_arr[(cable_i*2+2),freq_use[freq_i]]))
                
                
                
            endfor ;end freq for loop
            
          endfor ;end tile for loop
          
        endfor ;end obs for loop
        
        
        ;Fill the bandpass change ptrs for each cable
        *bandpass_arr_change_yy[cable_i]=bandpass_arr_change_yy_temp
        *bandpass_arr_change_xx[cable_i]=bandpass_arr_change_xx_temp
        
      endfor ;end cable for loop
      
      ;Set up ptrs for gaussian fit variables. This is not used at the moment, but could be used if all of the gaussians are plotted together after cable loop
      yfit_ptr=Ptrarr(n_cable,n_pol,/allocate)
      binCenters_ptr=Ptrarr(n_cable,n_pol,/allocate)
      maxfit_ptr=Ptrarr(n_cable,n_pol,/allocate)
      centerfit_ptr=Ptrarr(n_cable,n_pol,/allocate)
      fwhm_ptr=Ptrarr(n_cable,n_pol,/allocate)
      sigma_ptr=Ptrarr(n_cable,n_pol,/allocate)
      ;sum_ptr=Ptrarr(n_cable,n_pol,/allocate)
      ind_ptr=PTRARR(20,n_cable,/allocate)
      
      
      ;Histograms of fit bandpass and input gain fractional changes
      for cable_i=0, n_cable-1 do begin
        image_path=file_path_fhd
        
        ;Binsize which seems to work for the log axis histograms
        binsize=.005
        
        
        ;Log axis histograms for xx. Flagged tile spots are removed via the where statement and the initialized value (12345678).
        ;****xx histo plot
        cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_histogram_xx_cable'+strtrim(string(cable_i),1)+'.png',/quiet,/nomatch
        
        ;title=string(pointing_num)+ ': (B!Icp!N-g!Iunfit!N)/g!Iunfit!N for cable'+strtrim(string(cable_i),1)+', xx'
        cgHistoplot, (*bandpass_arr_change_xx[cable_i])[where(*bandpass_arr_change_xx[cable_i] NE 12345678)], binsize=binsize, $
          /FILL,title=string(pointing_num)+ ': g!Ibp,input!N for cable'+strtrim(string(cable_i),1)+', xx',$
          Charsize=1, histdata=h, locations=loc, reverse_indices=ri;, ytickformat='(F8.2)',/LOG
          
        ;Data needed for the gaussing fit of the histogram
        binCenters = loc + (binsize / 2.0)
        
        ;Sum may be used to normalize if needed
        ;sum=total((*bandpass_arr_change_yy[cable_i])[where(*bandpass_arr_change_yy[cable_i] NE 12345678)])
        
        ;Fit the histogram with a gaussian.  coeff is an array of three variables used to calculate the gaussian fit
        yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
        
        ;Overplot the gaussian on the histogram plot
        cgPlot, binCenters, yfit, COLOR='navy', THICK=2, /OVERPLOT
        
        ;Put the gaussian variables in a printable format for graph text
        maxfit = String(coeff[0], FORMAT='(I0)')
        centerfit = String(coeff[1], FORMAT='(F0.4)')
        fwhm = String(2 * SQRT(2 * ALOG(2)) * coeff[2], FORMAT='(F0.2)')
        sigma = String(coeff[2], FORMAT='(F0.4)')
        
        ;Plot the two sigma line (can be changed) in order to represent the two sigma resistant mean used to calculate the bandpass fit. Plot text of sigma variables if wanted
        cgPlot, [coeff[2]*1.+coeff[1],coeff[2]*1.+coeff[1]],[0.000001,coeff[0]], Linestyle=2, /OVERPLOT
        cgPlot, [-coeff[2]*1.+coeff[1],-coeff[2]*1.+coeff[1]],[0.000001,coeff[0]], Linestyle=2, /OVERPLOT
        cgText, 0.6, 0.80, /NORMAL, 'Maximum: ' + maxfit, COLOR='navy', charsize=1
        cgText, 0.6, 0.75, /NORMAL, 'Expected Value: ' + centerfit, COLOR='navy', charsize=1
        cgText, 0.6, 0.70, /NORMAL, '$\sigma$: ' + sigma, COLOR='navy', charsize=1
        
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        ;****end of xx histo plot
        
        ;Remove the fitted gaussian
        ;****xx histo gauss removed plot
        cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_histogram_gaussremove_xx_cable'+strtrim(string(cable_i),1)+'.png',/quiet,/nomatch
        ;title=string(pointing_num)+ ': Gauss removed (B!Icp!N-g!Iunfit!N)/g!Iunfit!N, cable'+strtrim(string(cable_i),1)+', xx'
        cgPlot,binCenters,h-yfit,title=string(pointing_num)+ ': Gauss removed g!Iunfit!N, cable'+strtrim(string(cable_i),1)+', xx',charsize=1, psym=10
        cgPlot, [-.2,.2],[0,0], linestyle=2, /overplot
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        ;****end of xx histo gauss removed plot
        
        ;n_bins=(size(h))[1]
        ;for bin_i=0,9 do begin
        ;  if ri[(n_bins/2+10)+(bin_i+1)] ne ri[(n_bins/2+10)+bin_i] then begin
        ;    indices=ri[ri[(n_bins/2+10)+bin_i]:ri[(n_bins/2+10)+bin_i+1]-1]
        ;    ind_arr=INTARR((size(indices))[1],3)
        ;    for ind_i=0, (size(indices))[1]-1 do begin
        ;      ind_obs= indices[ind_i] mod obsid_count
        ;      ind_freq= (indices[ind_i] / obsid_count) mod nf_use;max(nt_use_cable[cable_i,*])
        ;      ind_tile= indices[ind_i] / (obsid_count*nf_use);*max(nt_use_cable[cable_i,*]))
        ;      ind_arr[ind_i,0]=ind_obs
        ;      ind_arr[ind_i,1]=ind_freq
        ;      ind_arr[ind_i,2]=ind_tile
        ;    endfor
        ;    *ind_ptr[bin_i,cable_i]=ind_arr
        ;  endif
        ;endfor
        
        
        
        ;binsize_int=1.
        
        ;ind_arr_obs=0
        ;ind_arr_freq=0
        ;ind_arr_tile=0
        ;for bin_i=0,19 do begin
        ;  if keyword_set((*ind_ptr[bin_i,cable_i])) then begin
        ;    ind_arr_obs=[ind_arr_obs,reform((*ind_ptr[bin_i,cable_i])[*,0])]
        ;    ind_arr_freq=[ind_arr_freq,reform((*ind_ptr[bin_i,cable_i])[*,1])]
        ;    ind_arr_tile=[ind_arr_tile,reform((*ind_ptr[bin_i,cable_i])[*,2])]
        ;  endif
        ;endfor
        
        ;cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_change_histogram_xx_badgains_obs'+strtrim(string(cable_i),1)+'.png',/quiet,/nomatch
        ;cgHistoplot, ind_arr_obs, binsize=binsize_int, $
        ;  /FILL,title=string(pointing_num)+ ': potential bad gains, obs, '+strtrim(string(cable_i),1)+', xx',$
        ;  Charsize=1
        ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        ;cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_change_histogram_xx_badgains_freq'+strtrim(string(cable_i),1)+'.png',/quiet,/nomatch
        ;cgHistoplot, ind_arr_freq, binsize=binsize_int, $
        ;  /FILL,title=string(pointing_num)+ ': potential bad gains, freq, '+strtrim(string(cable_i),1)+', xx',$
        ;  Charsize=1
        ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        ;cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_change_histogram_xx_badgains_tile'+strtrim(string(cable_i),1)+'.png',/quiet,/nomatch
        ;cgHistoplot, ind_arr_tile, binsize=binsize_int, $
        ;  /FILL,title=string(pointing_num)+ ': potential bad gains, tile, '+strtrim(string(cable_i),1)+', xx',$
        ;  Charsize=1
        ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
        *yfit_ptr[cable_i,0]=yfit
        *binCenters_ptr[cable_i,0]=binCenters
        *maxfit_ptr[cable_i,0]=maxfit
        *centerfit_ptr[cable_i,0]=centerfit
        *fwhm_ptr[cable_i,0]=fwhm
        *sigma_ptr[cable_i,0]=sigma
        ;*sum_ptr[cable_i,0]=sum
        
        cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_histogram_yy_cable'+strtrim(string(cable_i),1)+'.png',/quiet,/nomatch
        ;title=string(pointing_num)+ ': (B!Icp!N-g!Iunfit!N)/g!Iunfit!N for cable'+strtrim(string(cable_i),1)+', yy'
        cgHistoplot, (*bandpass_arr_change_yy[cable_i])[where(*bandpass_arr_change_yy[cable_i] NE 12345678)], binsize=binsize, $
          /FILL,title=string(pointing_num)+ ': g!Ibp,input!N for cable'+strtrim(string(cable_i),1)+', yy',$
          Charsize=1, histdata=h, locations=loc;,ytickformat='(F8.2)',/LOG
        sum=total((*bandpass_arr_change_yy[cable_i])[where(*bandpass_arr_change_yy[cable_i] NE 12345678)])
        binCenters = loc + (binsize / 2.0)
        ;weights=1.0/h^2.
        ;a=[10000,1,0.03,0,0,0]
        ;yfit= CurveFit(binCenters, h, weights, a)
        
        yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
        cgPlot, binCenters, yfit, COLOR='navy', THICK=2, /OVERPLOT
        maxfit = String(coeff[0], FORMAT='(I0)')
        centerfit = String(coeff[1], FORMAT='(F0.4)')
        fwhm = String(2 * SQRT(2 * ALOG(2)) * coeff[2], FORMAT='(F0.2)')
        sigma = String(coeff[2], FORMAT='(F0.4)')
        cgPlot, [coeff[2]*1.+coeff[1],coeff[2]*1.+coeff[1]],[0.000001,coeff[0]], Linestyle=2, /OVERPLOT
        cgPlot, [-coeff[2]*1.+coeff[1],-coeff[2]*1.+coeff[1]],[0.000001,coeff[0]], Linestyle=2, /OVERPLOT
        cgText, 0.6, 0.80, /NORMAL, 'Maximum: ' + maxfit, COLOR='navy', charsize=1
        cgText, 0.6, 0.75, /NORMAL, 'Expected Value: ' + centerfit, COLOR='navy', charsize=1
        cgText, 0.6, 0.70, /NORMAL, '$\sigma$: ' + sigma, COLOR='navy', charsize=1
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
        cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_histogram_gaussremove_yy_cable'+strtrim(string(cable_i),1)+'.png',/quiet,/nomatch
        ;title=string(pointing_num)+ ': Gauss removed (B!Icp!N-g!Iunfit!N)/g!Iunfit!N, cable'+strtrim(string(cable_i),1)+', yy'
        cgPlot,binCenters,h-yfit,title=string(pointing_num)+ ': Gauss removed g!Iunfit!N, cable'+strtrim(string(cable_i),1)+', yy',charsize=1, psym=10
        cgPlot, [-.2,.2],[0,0], linestyle=2, /overplot
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
        *yfit_ptr[cable_i,1]=yfit
        *binCenters_ptr[cable_i,1]=binCenters
        *maxfit_ptr[cable_i,1]=maxfit
        *centerfit_ptr[cable_i,1]=centerfit
        *fwhm_ptr[cable_i,1]=fwhm
        *sigma_ptr[cable_i,1]=sigma
      ;*sum_ptr[cable_i,1]=sum
        
      endfor
      
      ;cgplot,freq[freq_use],bandpass_arr_change_xx[*,freq_use,*],color='blue',title=string(pointing_num)+ ': (cable-orig)/orig xx',xtitle='Frequency [MHz]',ytitle='Gain',$
      ;  yrange=yrange,xrange=xrange,/noerase,axiscolor=axiscolor,psym=2,symsize=0.2
      ;endif
      
      ;undefine, bp_change
      ;bp_change=1
      IF Keyword_Set(bp_change) THEN BEGIN
      
        ;Set up the bandpass change array, which is the difference between pointing calibration by cable fit and regular calibration (no cable)
        ;that was calculated on the fly above in the comparison block. There is the option to make it a fractional change
        bandpass_arr_change=bandpass_arr
        FOR i=1, 6 DO bandpass_arr_change[i*2-1,*]=(bandpass_arr[i*2-1,*]-bandpass_arr_orig[1,*]);/bandpass_arr_orig[1,*]
        FOR i=1, 6 DO bandpass_arr_change[i*2,*]=(bandpass_arr[i*2,*]-bandpass_arr_orig[2,*]);/bandpass_arr_orig[2,*]
        
        basename=file_basename(file_path_fhd)
        dirpath=file_dirname(file_path_fhd)
        ;image_path=filepath(basename,root=dirpath,sub='output_images')
        ;image_path='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_pointing_samecheck/output_images/'
        image_path=file_path_fhd
        IF file_test(file_dirname(image_path),/directory) EQ 0 THEN file_mkdir,file_dirname(image_path)
        ;export_path=filepath(basename,root=dirpath,sub='calibration')
        export_path=file_path_fhd
        IF file_test(file_dirname(export_path),/directory) EQ 0 THEN file_mkdir,file_dirname(export_path)
        Textfast,bandpass_arr_change,/write,file_path=export_path+strtrim(string(pointing_num),1)+'_bandpass_change'
        
        freq=freq_arr/1E6
        xtickv=[ceil(min(freq)/10)*10,floor(max(freq)/10)*10]
        xtickname=strtrim(round(xtickv),2)
        xrange=[min(freq)-(max(freq)-min(freq))/8,max(freq)+(max(freq)-min(freq))/8]
        yrange=[min(bandpass_arr_change[1:(n_pol*6),freq_use]),max(bandpass_arr_change[1:(n_pol*6),freq_use])]
        ;yrange=[-0.00001,0.00001]
        ytickv=[yrange[0],mean(yrange),yrange[1]]
        axiscolor='black'
        cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_change_xx_cc.png',/quiet,/nomatch
        cgplot,freq[freq_use],bandpass_arr_change[1,freq_use],color='blue',title=string(pointing_num)+ ': (cable-orig) xx',xtitle='Frequency [MHz]',ytitle='Gain',$
          yrange=yrange,xrange=xrange,/noerase,axiscolor=axiscolor,psym=2,symsize=0.2,CHARSIZE=1
        ;cgplot,freq[freq_use],bandpass_arr_change[1,freq_use],color='blue',title=string(pointing_num)+ ': (cable-orig) xx',xtitle='Frequency [MHz]',ytitle='Gain',$
        ;    yrange=yrange,xrange=xrange,/noerase,axiscolor=axiscolor,psym=2,symsize=0.2
          
        cgoplot,freq[freq_use],bandpass_arr_change[3,freq_use],color='red',psym=2,symsize=0.2
        cgoplot,freq[freq_use],bandpass_arr_change[5,freq_use],color='green',psym=2,symsize=0.2
        cgoplot,freq[freq_use],bandpass_arr_change[7,freq_use],color='purple',psym=2,symsize=0.2
        cgoplot,freq[freq_use],bandpass_arr_change[9,freq_use],color='yellow',psym=2,symsize=0.2
        cgoplot,freq[freq_use],bandpass_arr_change[11,freq_use],color='black',psym=2,symsize=0.2
        ;cgLegend, Title=['90m cables ('+ Strtrim(String(N_elements(tile_use_90)),1) +')','150m cables ('+ Strtrim(String(N_elements(tile_use_150)),1) +')', $
        ;    '230m cables ('+ Strtrim(String(N_elements(tile_use_230)),1) +')', '320m cables ('+ Strtrim(String(N_elements(tile_use_320)),1) +')', $
        ;    '400m cables ('+ Strtrim(String(N_elements(tile_use_400)),1) +')','524m cables ('+ Strtrim(String(N_elements(tile_use_524)),1) +')'], $
        ;    Color=['blue','red','green','purple','yellow','black'],Psym=[2,2,2,2,2,2], $
        ;    Length=0.0,Location=[0.55,0.85]
        
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
        cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_change_yy_cc.png',/quiet,/nomatch
        cgplot,freq[freq_use],bandpass_arr_change[2,freq_use],color='blue',title=string(pointing_num) + ': (cable-orig) yy',xtitle='Frequency [MHz]',ytitle='Gain',$
          yrange=yrange,xrange=xrange,/noerase,axiscolor=axiscolor,psym=2,symsize=0.2,CHARSIZE=1
        ;cgplot,freq[freq_use],bandpass_arr_change[2,freq_use],color='blue',title=string(pointing_num) + ': (cable-orig) yy',xtitle='Frequency [MHz]',ytitle='Gain',$
        ;    yrange=yrange,xrange=xrange,/noerase,axiscolor=axiscolor,psym=2,symsize=0.2
          
        cgoplot,freq[freq_use],bandpass_arr_change[4,freq_use],color='red',psym=2,symsize=0.2
        cgoplot,freq[freq_use],bandpass_arr_change[6,freq_use],color='green',psym=2,symsize=0.2
        cgoplot,freq[freq_use],bandpass_arr_change[8,freq_use],color='purple',psym=2,symsize=0.2
        cgoplot,freq[freq_use],bandpass_arr_change[10,freq_use],color='yellow',psym=2,symsize=0.2
        cgoplot,freq[freq_use],bandpass_arr_change[12,freq_use],color='black',psym=2,symsize=0.2
        ;cgLegend, Title=['90m cables ('+ Strtrim(String(N_elements(tile_use_90)),1) +')','150m cables ('+ Strtrim(String(N_elements(tile_use_150)),1) +')', $
        ;    '230m cables ('+ Strtrim(String(N_elements(tile_use_230)),1) +')', '320m cables ('+ Strtrim(String(N_elements(tile_use_320)),1) +')', $
        ;    '400m cables ('+ Strtrim(String(N_elements(tile_use_400)),1) +')','524m cables ('+ Strtrim(String(N_elements(tile_use_524)),1) +')'], $
        ;    Color=['blue','red','green','purple','yellow','black'],Psym=[2,2,2,2,2,2], $
        ;    Length=0.0,Location=[0.55,0.85]
        
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      endif ;end of bpchange if
      
    endif else begin;cablesplit5 if end
    
      ;Set up ptrs to store data for plotting  -- can be changed if want to plot them seperately
      bandpass_arr_change_xx=ptrarr(2, /allocate)
      bandpass_arr_change_yy=ptrarr(2, /allocate)
      
      For cable_i=5, n_cable-1 do begin
      
        ;Create temporary arrays for the bandpass change. Hacky solution is used to account for the fact that some tiles are flagged in an observation,
        ;so the array sizes change. Random value of 12345678 is inputted as an initialization, which will be weeded out during actual plotting and calculations
        bandpass_arr_change_xx_temp=fltarr(obsid_count,nf_use,Max(nt_use_cable[cable_i,*]),/NoZero)
        bandpass_arr_change_yy_temp=fltarr(obsid_count,nf_use,Max(nt_use_cable[cable_i,*]),/NoZero)
        bandpass_arr_change_xx_temp[*,*,*]=12345678.
        bandpass_arr_change_yy_temp[*,*,*]=12345678.
        
        inside_list_cablearr_index_histo_ptr=PTRARR(obsid_count,/allocate)
        outside_list_cablearr_index_histo_ptr=PTRARR(obsid_count,/allocate)
        
        For obs_i = 0, obsid_count-1 do begin
        
          ;index and index_size account for the fact that the tile array changes per observation due to flagging.
          index_size=size(*tile_use_cable[cable_i,obs_i])
          index=INTARR(index_size[1])
          
          ;Find where the used tiles per observation array holds the current tile that is used for the cable. Save that location in an array for later referencing.
          for ind_i=0, index_size[1]-1 do index[ind_i]=where(*t_use[obs_i] EQ (*tile_use_cable[cable_i,obs_i])[ind_i])
          
          for tile_i=0,nt_use_cable[cable_i,obs_i]-1 do begin
          
            for freq_i=0, nf_use-1 do begin
            
              if where(tile_i eq inside_list_cablearr_index) NE -1 then begin
              
                ;Take the fitted bandpass per pointing per cable and find the fractional difference of the inputted (unfitted) gains for that cable type.
                bandpass_arr_change_xx_temp[obs_i,freq_i,tile_i]=((abs(bandpass_arr[(cable_i*2+3),freq_use[freq_i]])-$
                  abs((*gain_arr_ptr_fitted[obs_i,0])[freq_use[freq_i],index[tile_i]])))/abs((*gain_arr_ptr_fitted[obs_i,0])[freq_use[freq_i],index[tile_i]])
                  
                bandpass_arr_change_yy_temp[obs_i,freq_i,tile_i]=((abs(bandpass_arr[(cable_i*2+4),freq_use[freq_i]])-$
                  abs((*gain_arr_ptr_fitted[obs_i,1])[freq_use[freq_i],index[tile_i]])))/abs((*gain_arr_ptr_fitted[obs_i,1])[freq_use[freq_i],index[tile_i]])
                  
              ;Assuming yy is good enough
              ;If (bandpass_arr_change_yy_temp[obs_i,freq_i,tile_i] NE 12345678. AND freq_i EQ 0) Then inside_list_cablearr_index_histo_temp=[inside_list_cablearr_index_histo_temp,tile_i]
                  
              endif
              
              if where(tile_i eq outside_list_cablearr_index) NE -1 then begin
              
                ;Take the fitted bandpass per pointing per cable and find the fractional difference of the inputted (unfitted) gains for that cable type.
                bandpass_arr_change_xx_temp[obs_i,freq_i,tile_i]=((abs(bandpass_arr[(cable_i*2+1),freq_use[freq_i]])-$
                  abs((*gain_arr_ptr_fitted[obs_i,0])[freq_use[freq_i],index[tile_i]])))/abs((*gain_arr_ptr_fitted[obs_i,0])[freq_use[freq_i],index[tile_i]])
                  
                bandpass_arr_change_yy_temp[obs_i,freq_i,tile_i]=((abs(bandpass_arr[(cable_i*2+2),freq_use[freq_i]])-$
                  abs((*gain_arr_ptr_fitted[obs_i,1])[freq_use[freq_i],index[tile_i]])))/abs((*gain_arr_ptr_fitted[obs_i,1])[freq_use[freq_i],index[tile_i]])
                  
              ;If (bandpass_arr_change_yy_temp[obs_i,freq_i,tile_i] NE 12345678. AND freq_i EQ 0) Then outside_list_cablearr_index_histo_temp=[outside_list_cablearr_index_histo_temp,tile_i]
                  
              endif
              
            endfor ;end freq for loop
            
          endfor ;end tile for loop
          
        ;*inside_list_cablearr_index_histo_ptr[obs_i]=inside_list_cablearr_index_histo_temp
        ;*outside_list_cablearr_index_histo_ptr[obs_i]=outside_list_cablearr_index_histo_temp
          
        endfor ;end obs for loop
        
        ;Fill the bandpass change ptrs for each cable  -- can be changed if you want to plot them seperately
        *bandpass_arr_change_yy[1]=bandpass_arr_change_yy_temp
        *bandpass_arr_change_xx[0]=bandpass_arr_change_xx_temp
        
      endfor ;end cable for loop
      
      
      ;Histograms of fit bandpass and input gain fractional changes
      for cable_i=5, n_cable-1 do begin
        image_path=file_path_fhd
        
        ;Binsize which seems to work for the log axis histograms
        binsize=.005
        
        cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_change_histogram_log_xx_cable'+strtrim(string(cable_i),1)+'_inside.png',/quiet,/nomatch
        
        bandpass_arr_change_xx_histo_temp=(*bandpass_arr_change_xx[0])[*,*,inside_list_cablearr_index]
        cgHistoplot, bandpass_arr_change_xx_histo_temp[where(bandpass_arr_change_xx_histo_temp NE 12345678.)], binsize=binsize, $
          /FILL,title=string(pointing_num)+ ': (B!Icp!N-g!Iunfit!N)/g!Iunfit!N for cable'+strtrim(string(cable_i),1)+', xx',$
          Charsize=1, xrange=[-.2,.2], histdata=h_inside, locations=loc_inside, reverse_indices=ri, ytickformat='(F8.2)',/LOG
        binCenters_inside = loc_inside + (binsize / 2.0)
        yfit_inside = GaussFit(binCenters_inside, h_inside, coeff_inside, NTERMS=3)
        cgPlot, binCenters_inside, yfit_inside, COLOR='navy', THICK=2, /OVERPLOT
        
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
        
        cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_change_histogram_log_xx_cable'+strtrim(string(cable_i),1)+'_outside.png',/quiet,/nomatch
        
        bandpass_arr_change_xx_histo_temp=(*bandpass_arr_change_xx[0])[*,*,outside_list_cablearr_index]
        ;Could have 12345678 errors
        cgHistoplot, bandpass_arr_change_xx_histo_temp[where(bandpass_arr_change_xx_histo_temp NE 12345678.)], binsize=binsize, $
          /FILL,title=string(pointing_num)+ ': (B!Icp!N-g!Iunfit!N)/g!Iunfit!N for cable'+strtrim(string(cable_i),1)+', xx',$
          Charsize=1, xrange=[-.2,.2], histdata=h_outside, locations=loc_outside, reverse_indices=ri, ytickformat='(F8.2)',/LOG
        binCenters_outside = loc_outside + (binsize / 2.0)
        yfit_outside = GaussFit(binCenters_outside, h_outside, coeff_outside, NTERMS=3)
        cgPlot, binCenters_outside, yfit_outside, COLOR='navy', THICK=2, /OVERPLOT
        
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
        ;Log axis histograms for xx. Flagged tile spots are removed via the where statement and the initialized value (12345678).
        ;****xx histo plot
        cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_change_histogram_log_xx_cable'+strtrim(string(cable_i),1)+'.png',/quiet,/nomatch
        
        
        cgHistoplot, (*bandpass_arr_change_xx[0])[where(*bandpass_arr_change_xx[0] NE 12345678)], binsize=binsize, $
          /FILL,title=string(pointing_num)+ ': (B!Icp!N-g!Iunfit!N)/g!Iunfit!N for cable'+strtrim(string(cable_i),1)+', xx',$
          Charsize=1, xrange=[-.2,.2], histdata=h, locations=loc, reverse_indices=ri, ytickformat='(F8.2)',/LOG
          
        ;Data needed for the gaussing fit of the histogram
        ;binCenters = loc + (binsize / 2.0)
          
        ;Sum may be used to normalize if needed
        ;sum=total((*bandpass_arr_change_yy[cable_i])[where(*bandpass_arr_change_yy[cable_i] NE 12345678)])
          
        ;Fit the histogram with a gaussian.  coeff is an array of three variables used to calculate the gaussian fit
        ;yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
          
        ;Overplot the gaussian on the histogram plot
        ;cgPlot, binCenters, yfit, COLOR='navy', THICK=2, /OVERPLOT
        cgPlot, binCenters_inside, yfit_inside, COLOR='navy', THICK=2, /OVERPLOT
        cgPlot, binCenters_outside, yfit_outside, COLOR='navy', THICK=2, /OVERPLOT
        
        ;Put the gaussian variables in a printable format for graph text
        ;maxfit = String(coeff[0], FORMAT='(I0)')
        ;centerfit = String(coeff[1], FORMAT='(F0.4)')
        ;fwhm = String(2 * SQRT(2 * ALOG(2)) * coeff[2], FORMAT='(F0.2)')
        ;sigma = String(coeff[2], FORMAT='(F0.4)')
        
        ;Plot the two sigma line (can be changed) in order to represent the two sigma resistant mean used to calculate the bandpass fit. Plot text of sigma variables if wanted
        ;cgPlot, [coeff[2]*1.+coeff[1],coeff[2]*1.+coeff[1]],[0.000001,coeff[0]], Linestyle=2, /OVERPLOT
        ;cgPlot, [-coeff[2]*1.+coeff[1],-coeff[2]*1.+coeff[1]],[0.000001,coeff[0]], Linestyle=2, /OVERPLOT
        ;cgText, 0.6, 0.80, /NORMAL, 'Maximum: ' + maxfit, COLOR='navy', charsize=1
        ;cgText, 0.6, 0.75, /NORMAL, 'Expected Value: ' + centerfit, COLOR='navy', charsize=1
        ;cgText, 0.6, 0.70, /NORMAL, '$\sigma$: ' + sigma, COLOR='navy', charsize=1
        
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      ;****end of xx histo plot
        
      ;Remove the fitted gaussian
      ;****xx histo gauss removed plot
      ;cgPS_Open,image_path+strtrim(string(pointing_num),1)+'_bandpass_change_histogram_gaussremove_xx_cable'+strtrim(string(cable_i),1)+'.png',/quiet,/nomatch
      ;cgPlot,binCenters,h-yfit,title=string(pointing_num)+ ': Gauss removed (B!Icp!N-g!Iunfit!N)/g!Iunfit!N, cable'+strtrim(string(cable_i),1)+', xx',charsize=1, xrange=[-.2,.2], psym=10
      ;cgPlot, [-.2,.2],[0,0], linestyle=2, /overplot
      ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      ;****end of xx histo gauss removed plot
      endfor
      
      
      
    endelse ;end of cablesplit5 else
  ENDIF
  
  
  
  ;RETURN,cal_bandpass
  
  RETURN,bandpass_arr
END
