pro vis_cal_temp

  obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included/metadata/1061316176_obs.sav', 'obs')
  params = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included/metadata/1061316176_params.sav', 'params')
  cal = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/calibration/1061316176_cal.sav', 'cal')
  
  vis_model_ptr=PTRARR(2,/allocate)
  vis_model_ptr[0] = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included/vis_data/1061316176_vis_model_XX.sav', 'vis_model_ptr')
  vis_model_ptr[1] = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included/vis_data/1061316176_vis_model_YY.sav', 'vis_model_ptr')
  
  
  vis_weight_ptr = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included/vis_data/1061316176_flags.sav','flag_arr')
  
  
  
  IF N_Elements(cal) EQ 0 THEN cal=fhd_struct_init_cal(obs,params,_Extra=extra)
  reference_tile=cal.ref_antenna
  min_baseline=obs.min_baseline
  max_baseline=obs.max_baseline
  dimension=obs.dimension
  elements=obs.elements
  double_precision=0
  IF Tag_Exist(obs, 'double_precision') THEN double_precision=obs.double_precision
  
  min_cal_baseline=cal.min_cal_baseline
  max_cal_baseline=cal.max_cal_baseline
  min_cal_solutions=cal.min_solns ;minimum number of calibration equations needed to solve for the gain of one baseline
  time_average=cal.time_avg
  max_cal_iter=cal.max_iter
  IF max_cal_iter LT 5 THEN print,'Warning! At least 5 calibration iterations recommended. Using '+Strn(Floor(max_cal_iter))
  conv_thresh=cal.conv_thresh
  
  n_pol=cal.n_pol
  n_freq=cal.n_freq
  n_tile=cal.n_tile
  n_time=cal.n_time
  
  vis_weight_ptr_use=vis_weight_ptr ;weights WILL be over-written! (Only for NAN gain solutions)
  tile_A_i=cal.tile_A-1
  tile_B_i=cal.tile_B-1
  freq_arr=cal.freq
  bin_offset=cal.bin_offset
  n_baselines=obs.nbaselines
  IF Tag_exist(cal,'phase_iter') THEN phase_fit_iter=cal.phase_iter ELSE phase_fit_iter=Floor(max_cal_iter/4.)<4
  
  kbinsize=obs.kpix
  
  FOR pol_i=0,n_pol-1 DO BEGIN
  
    ;inds_A=Rebin(Lindgen(n_freq),n_freq,n_baselines,/sample)+Rebin(transpose(tile_A_i)*n_freq,n_freq,n_baselines)
    ;inds_B=Rebin(Lindgen(n_freq),n_freq,n_baselines,/sample)+Rebin(transpose(tile_B_i)*n_freq,n_freq,n_baselines)
  
    vis_weight_use=0>*vis_weight_ptr_use[pol_i]<1
    vis_model=Temporary(*vis_model_ptr[pol_i])
    ;vis_model=Temporary(vis_model)*vis_weight_use
    weight=Temporary(vis_weight_use)
    
    kx_arr=cal.uu/kbinsize
    ky_arr=cal.vv/kbinsize
    kr_arr=Sqrt(kx_arr^2.+ky_arr^2.)
    dist_arr=(freq_arr#Temporary(kr_arr))*kbinsize
    xcen=freq_arr#Abs(Temporary(kx_arr))
    ycen=freq_arr#Abs(Temporary(ky_arr))
    flag_dist_cut=where((Temporary(dist_arr) LT min_cal_baseline) OR (Temporary(xcen) GT dimension/2.) OR (Temporary(ycen) GT elements/2.),n_dist_cut)
    IF n_dist_cut GT 0 THEN weight[flag_dist_cut]=0.
    kx_arr=(ky_arr=(dist_arr=0))
    
    
    tile_use_flag=(*obs.baseline_info).tile_use
    tile_use_flag[76]=0
    freq_use_flag=(*obs.baseline_info).freq_use
    
    
    freq_weight=Total(weight,2)
    baseline_weight=Total(weight,1)
    freq_use=where((freq_weight GT 0) AND (freq_use_flag GT 0),n_freq_use)
    baseline_use=where(baseline_weight,n_baseline_use)
    hist_tile_A=histogram(tile_A_i[baseline_use],min=0,/bin,max=n_tile-1,reverse_ind=riA)
    hist_tile_B=histogram(tile_B_i[baseline_use],min=0,/bin,max=n_tile-1,reverse_ind=riB)
    tile_use=where(((hist_tile_A+hist_tile_B) GT 0) AND (tile_use_flag GT 0),n_tile_use)
    
    tile_A_i_use=Lonarr(n_baseline_use)
    tile_B_i_use=Lonarr(n_baseline_use)
    FOR tile_i=0L,n_tile_use-1 DO BEGIN
      IF hist_tile_A[tile_use[tile_i]] GT 0 THEN tile_A_i_use[riA[riA[tile_use[tile_i]]:riA[tile_use[tile_i]+1]-1]]=tile_i
      IF hist_tile_B[tile_use[tile_i]] GT 0 THEN tile_B_i_use[riB[riB[tile_use[tile_i]]:riB[tile_use[tile_i]+1]-1]]=tile_i
    ENDFOR
    
    nan_i=where(Finite(vis_model,/nan),n_nan)
    IF n_nan GT 0 THEN vis_model[nan_i]=0
    
    gain_arr = complex(FLTARR(n_freq,n_tile_use))
    gain_arr[*,*]=1.
    
    
    max_cal_iter=100
    phase_fit_iter=10
    min_cal_solutions=10
    
    gain_arr_tau = (*cal.gain[0])
    tile_flags = where(tile_use_flag EQ 0, n_flag_tiles)
    if n_flag_tiles GT 0 then gain_arr_tau[*,tile_flags]=0
    FOR fii=0L,n_freq_use-1 DO BEGIN
      fi=freq_use[fii]
      gain_curr=reform(gain_arr[fi,tile_use])
      
      b_i_use=where(weight[fi,baseline_use] GT 0,n_baseline_weight)
      A_ind=tile_A_i_use & A_ind=A_ind[b_i_use]
      B_ind=tile_B_i_use & B_ind=B_ind[b_i_use]
      
      ;gain_arr_tau = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/gain_arr_tau_notimeavg.sav','gain_arr')
      
      ;gain_arr_old = *cal.gain[0] +*cal.gain_residual[0]
      vis_gain_old=reform(gain_arr_tau[fi,A_ind]*conj(gain_arr_tau[fi,B_ind]))
      ;vis_cal_old=fft(fft(vis_model,dimension=1)*(vis_gain_old),dimension=1,/inverse)
      vis_cal_old=reform(vis_model[fi,baseline_use[b_i_use]]);*(vis_gain_old))
      ;vis_tau_model = vis_model;*vis_gain_old
      
      
      weight2=reform(weight[fi,baseline_use[b_i_use]])
      vis_model2=reform(vis_model[fi,baseline_use[b_i_use]]);*weight_invert(weight2)
      ;vis_gains_freq=vis_gains_freq[b_i_use]
      
      
      
      A_ind_arr=Ptrarr(n_tile_use,/allocate)
      n_arr=Fltarr(n_tile_use)
      FOR tile_i=0L,n_tile_use-1 DO BEGIN
        ;should be set up so that using where is okay
        inds=where(A_ind EQ tile_i,n1)
        IF n1 GT 1 THEN *A_ind_arr[tile_i]=Reform(inds,1,n1) ELSE *A_ind_arr[tile_i]=-1
        n_arr[tile_i]=n1 ;NEED SOMETHING MORE IN CASE INDIVIDUAL TILES ARE FLAGGED FOR ONLY A FEW FREQUENCIES!!
      ENDFOR
      
      
      gain_A = complex(FLTARR(n_tile_use))
      
      
      FOR i=0L,(max_cal_iter-1)>1 DO BEGIN
        if (i EQ 20) OR (i EQ 40) then stop
        
        ;vis_gain_update = gain_curr[inds_A]*conj(gain_curr[inds_B])
        vis_gain_update = conj(gain_curr[B_ind])*vis_model2
        
        for tile_i=0, n_tile_use-1 do begin
          ;inds1 = where(cal.tile_a EQ tile_i,n_count1)
          ;inds2 = where(cal.tile_b[inds1] EQ tile_i,n_count2)
        
        
          ;if n_count2 GT 0 then gains_freq[*,tile_i] = sqrt(mean(gains_freq_full[*,inds1[inds2]], dimension=2))
        
          ;b_i_use=where(weight[fi,*A_ind_arr[tile_i]] GT 0,n_baseline_use2)
          ;b_i_use_data = where(abs(vis_gain_update[fi,(*A_ind_arr[tile_i])[b_i_use]]) GT 0, n_baseline_use3)
          ;b_i_use_model = where(abs(vis_cal_old[fi,(*A_ind_arr[tile_i])[b_i_use[b_i_use_data]]]) GT 0, n_baseline_use4)
          ;weight_inds = (*A_ind_arr[tile_i])[b_i_use[b_i_use_data[b_i_use_model]]]
          ;weight_inds=where(weight2[*A_ind_arr[tile_i]] GT 0,n_baseline_use2)
          ;if (n_baseline_use2 GT 10) AND (n_baseline_use3 GT 0) AND (n_baseline_use4 GT 0) AND (n_arr[tile_i] GE min_cal_solutions) then $
          ;if (n_baseline_use2 GT 10) AND (n_arr[tile_i] GE min_cal_solutions) then $
          ;  gain_A[tile_i]=LA_Least_Squares( vis_cal_old[*A_ind_arr[tile_i]], vis_gain_update[*A_ind_arr[tile_i]],  method=2, double=double_precision)
        ;if abs(gain_A[fi,tile_i]) LT .990 then stop
        ;if max(abs(vis_cal_old[fi,weight_inds] - vis_gain_update[fi,weight_inds])) GT 20 then begin
        ;  te = where((abs(vis_cal_old[fi,weight_inds] - vis_gain_update[fi,weight_inds])) GT 20,n_count)
        ;  print, minmax(B_ind[weight_inds[te]])
        ;  endif
        
          IF n_arr[tile_i] GE min_cal_solutions THEN BEGIN
            weight_inds=(where(weight2[*A_ind_arr[tile_i]] GT 0,n_baseline_use2))
            if n_baseline_use2 GT 0 then $
              gain_A[tile_i]=LA_Least_Squares( vis_cal_old[transpose((*A_ind_arr[tile_i])[weight_inds])], vis_gain_update[transpose((*A_ind_arr[tile_i])[weight_inds])],  method=2, double=double_precision)
          ENDIF

            
            
        ;print, minmax(abs(vis_cal_old[*,weight_inds] - vis_gain_update[*,weight_inds]))
            
        endfor
        
        
        
        gain_old=gain_curr
        IF phase_fit_iter-i GT 0 THEN gain_A*=Abs(gain_old)*weight_invert(Abs(gain_A))
        
        gain_curr=(gain_A+gain_old)/2.
        
      endfor
      stop
      gain_arr[fi,tile_use]=gain_curr
    endfor
    
  endfor
  
end