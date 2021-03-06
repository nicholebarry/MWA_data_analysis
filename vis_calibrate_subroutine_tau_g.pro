;function poly_cal_freq_mode, mode, val, mask
;
;  dims = size(mask, /dim)
;
;  case n_elements(dims) of
;    1: begin
;      n_freq = dims[0]
;      n_val = 1
;    end
;
;    2: begin
;      n_freq = dims[0]
;      n_val = dims[1]
;    end
;    else: stop
;  end
;
;  if n_elements(val) ne n_val then stop
;
;  val_arr = matrix_multiply(dblarr(n_freq)+1, val)
;  x_arr = (dindgen(n_freq)/(n_freq-1))*2-1
;  if n_val gt 1 then x_arr = rebin(x_arr, n_freq, n_val,/sample)
;
;  case mode of
;    0: freq_arr = val_arr
;    1: freq_arr = val_arr * x_arr
;    2: freq_arr = val_arr * (3*x_arr^2. - 1)/2.
;  endcase
;
;  ;freq_arr = sqrt((2*mode+1)/2.)*freq_arr
;
;  return, freq_arr*mask
;end
;
;function poly_cal_gain, gain_modes, mode_types, mode_num, mask, amp=amp, phase=phase
;
;  n_modes = n_elements(mode_types)
;  dims_gain = size(gain_modes, /dimension)
;  dims_mask = size(mask, /dimension)
;
;  if n_elements(dims_gain) gt 1 then if dims_gain[1] ne dims_mask[1] then message, 'second dimension of gain_modes and mask must match'
;
;  if dims_gain[0] ne n_modes then message, 'first dimension of gain_modes must match length of mode_types'
;
;
;  amp = dblarr(dims_mask)
;  phase = dblarr(dims_mask)
;  for fi=0L,n_modes-1 do begin
;    if mode_types[fi] eq 'amp' then amp += poly_cal_freq_mode(mode_num[fi], reform(gain_modes[fi,*]), mask) $
;    else phase += poly_cal_freq_mode(mode_num[fi], reform(gain_modes[fi,*]), mask)
;  endfor
;
;  if total(abs(amp)) eq 0 then amp +=1. ;; for phase-only calibration
;  if total(abs(phase)) gt 0 then gain = amp * exp(complex(0,1)*phase) else gain = amp
;
;  return, gain
;end
;
;function poly_cal_chi2, params
;
;  common poly_cal_data, model_vis, data_vis, poly_cal_gain_B, f_arr, mode_types, mode_num, mask
;
;  dims = size(mask, /dimension)
;  if n_elements(dims) gt 1 then n_vis = dims[1] else n_vis = 1
;
;  if n_elements(params) eq 1 then params = [params[0]]
;  if n_vis gt 1 then p = rebin(params, n_elements(params), n_vis) else p=params
;  gain_A = poly_cal_gain(p, mode_types, mode_num, mask, amp=amp)
;
;  data_use = data_vis*mask
;  model_use = model_vis*mask
;
;  diff = model_use*Conj(poly_cal_gain_B)*gain_A - data_use
;
;  ;; discourage negative amplitudes
;  if min(amp) lt 0 then return, total(abs(diff)^2.)*product(dims) else return, total(abs(diff)^2.)
;
;end
;
;function poly_cal_grad, params
;
;  common poly_cal_data, model_vis, data_vis, poly_cal_gain_B, f_arr, mode_types, mode_num, mask
;
;  dims = size(mask, /dimension)
;  if n_elements(dims) gt 1 then n_vis = dims[1] else n_vis = 1
;
;  if n_vis gt 1 then p = rebin(params, n_elements(params), n_vis) else p=params
;  gain_A = poly_cal_gain(p, mode_types, mode_num, mask, phase=phase_A)
;
;  data_use = data_vis*mask
;  model_use = model_vis*mask
;
;  amp_grad_factor = 2.*abs(model_use)^2.*abs(poly_cal_gain_B)^2.*abs(gain_A) - 2*Real_part(model_use*Conj(poly_cal_gain_B)*Conj(data_use)*exp(complex(0,1)*phase_A))
;  phase_grad_factor = 2*Imaginary(model_use*gain_A*Conj(poly_cal_gain_B)*Conj(data_use))
;
;  grad = params*0
;  for i=0, n_elements(mode_types)-1 do begin
;    if mode_types[i] eq 'amp' then grad[i] = total(amp_grad_factor*f_arr^(mode_num[i])) $
;    else grad[i] = total(phase_grad_factor*f_arr^(mode_num[i]))
;
;  endfor
;
;  return, grad
;
;end

FUNCTION vis_calibrate_subroutine_tau_G,vis_ptr,vis_model_ptr,vis_weight_ptr,obs,params,cal,preserve_visibilities=preserve_visibilities,$
    calib_freq_func=calib_freq_func,calibration_weights=calibration_weights,_Extra=extra
    
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
  
  cal_return=cal
  FOR pol_i=0,n_pol-1 DO cal_return.gain[pol_i]=Ptr_new(*cal.gain[pol_i])
  
  FOR pol_i=0,n_pol-1 DO BEGIN
    convergence=Fltarr(n_freq,n_tile)
    convergence_phase=Fltarr(n_freq,n_tile)
    convergence_plot=Fltarr(max_cal_iter)
    convergence_plot_phase=Fltarr(max_cal_iter)
    gain_arr=*cal.gain[pol_i]
    
    IF Keyword_Set(time_average) THEN BEGIN
      ;average over time
      ;the visibilities have dimension nfreq x (n_baselines x n_time),
      ; which can be reformed to nfreq x n_baselines x n_time
      tile_A_i=tile_A_i[0:n_baselines-1]
      tile_B_i=tile_B_i[0:n_baselines-1]
      vis_weight_use=0>Reform(*vis_weight_ptr_use[pol_i],n_freq,n_baselines,n_time)<1
      IF Keyword_Set(preserve_visibilities) THEN vis_model=Reform(*vis_model_ptr[pol_i],n_freq,n_baselines,n_time) $
      ELSE vis_model=Reform(Temporary(*vis_model_ptr[pol_i]),n_freq,n_baselines,n_time)
      vis_model=Total(Temporary(vis_model)*vis_weight_use,3)
      vis_measured=Reform(*vis_ptr[pol_i],n_freq,n_baselines,n_time)
      vis_avg=Total(Temporary(vis_measured)*vis_weight_use,3)
      weight=Total(Temporary(vis_weight_use),3)
      
      kx_arr=cal.uu[0:n_baselines-1]/kbinsize ;ignore slight variation with time
      ky_arr=cal.vv[0:n_baselines-1]/kbinsize
      kr_arr=Sqrt(kx_arr^2.+ky_arr^2.)
      dist_arr=(freq_arr#Temporary(kr_arr))*kbinsize
      xcen=freq_arr#Abs(kx_arr)
      ycen=freq_arr#Abs(ky_arr)
      
      IF Keyword_Set(calibration_weights) THEN BEGIN
        flag_dist_cut=where((dist_arr LT min_baseline) OR (Temporary(xcen) GT dimension/2.) OR (Temporary(ycen) GT elements/2.),n_dist_cut)
        IF min_cal_baseline GT min_baseline THEN taper_min=((Sqrt(2.)*min_cal_baseline-dist_arr)/min_cal_baseline)>0. ELSE taper_min=0.
        IF max_cal_baseline LT max_baseline THEN taper_max=((dist_arr-max_cal_baseline)/min_cal_baseline)>0. ELSE taper_max=0.
        baseline_weights=(1.-(taper_min+taper_max)^2.)>0.
      ENDIF ELSE flag_dist_cut=where((dist_arr LT min_cal_baseline) OR (Temporary(xcen) GT dimension/2.) OR (Temporary(ycen) GT elements/2.),n_dist_cut)
    ENDIF ELSE BEGIN
      vis_weight_use=0>*vis_weight_ptr_use[pol_i]<1
      IF Keyword_Set(preserve_visibilities) THEN vis_model=*vis_model_ptr[pol_i] $
      ELSE vis_model=Temporary(*vis_model_ptr[pol_i])
      vis_model=Temporary(vis_model)*vis_weight_use
      vis_avg=*vis_ptr[pol_i]*vis_weight_use
      weight=Temporary(vis_weight_use)
      
      kx_arr=cal.uu/kbinsize
      ky_arr=cal.vv/kbinsize
      kr_arr=Sqrt(kx_arr^2.+ky_arr^2.)
      dist_arr=(freq_arr#Temporary(kr_arr))*kbinsize
      xcen=freq_arr#Abs(Temporary(kx_arr))
      ycen=freq_arr#Abs(Temporary(ky_arr))
      IF Keyword_Set(calibration_weights) THEN BEGIN
        flag_dist_cut=where((Temporary(dist_arr) LT min_baseline) OR (Temporary(xcen) GT dimension/2.) OR (Temporary(ycen) GT elements/2.),n_dist_cut)
        IF min_cal_baseline GT min_baseline THEN taper_min=((Sqrt(2.)*min_cal_baseline-dist_arr)/min_cal_baseline)>0. ELSE taper_min=0.
        IF max_cal_baseline LT max_baseline THEN taper_max=((dist_arr-max_cal_baseline)/min_cal_baseline)>0. ELSE taper_max=0.
        baseline_weights=(1.-(taper_min+taper_max)^2.)>0.
      ENDIF ELSE flag_dist_cut=where((Temporary(dist_arr) LT min_cal_baseline) OR (Temporary(xcen) GT dimension/2.) OR (Temporary(ycen) GT elements/2.),n_dist_cut)
    ENDELSE
    kx_arr=(ky_arr=(dist_arr=0))
    
    IF n_dist_cut GT 0 THEN weight[flag_dist_cut]=0.
    vis_avg*=weight_invert(weight)
    vis_model*=weight_invert(weight)
    
    tile_use_flag=(*obs.baseline_info).tile_use
    freq_use_flag=(*obs.baseline_info).freq_use
    
    freq_weight=Total(weight,2)
    ;baseline_weight=Total(weight,1)
    ;TEMP;;;;;;;;;;;;;;;;;;;;;;;;;;
    baseline_weight=Total(*vis_ptr[pol_i],1)
    freq_use=where((freq_weight GT 0) AND (freq_use_flag GT 0),n_freq_use)
    ;baseline_use=where(baseline_weight,n_baseline_use)
    baseline_use=where(abs(baseline_weight) GT 0 ,n_baseline_use)
    
    
    hist_tile_A=histogram(tile_A_i[baseline_use],min=0,/bin,max=n_tile-1,reverse_ind=riA)
    hist_tile_B=histogram(tile_B_i[baseline_use],min=0,/bin,max=n_tile-1,reverse_ind=riB)
    tile_use=where(((hist_tile_A+hist_tile_B) GT 0) AND (tile_use_flag GT 0),n_tile_use)
    
    tile_A_i_use=Lonarr(n_baseline_use)
    tile_B_i_use=Lonarr(n_baseline_use)
    FOR tile_i=0L,n_tile_use-1 DO BEGIN
      IF hist_tile_A[tile_use[tile_i]] GT 0 THEN tile_A_i_use[riA[riA[tile_use[tile_i]]:riA[tile_use[tile_i]+1]-1]]=tile_i
      IF hist_tile_B[tile_use[tile_i]] GT 0 THEN tile_B_i_use[riB[riB[tile_use[tile_i]]:riB[tile_use[tile_i]+1]-1]]=tile_i
    ENDFOR
    
    ref_tile_use=Min(where(reference_tile EQ tile_use))
    IF ref_tile_use EQ -1 THEN BEGIN
      ref_tile_use=0L
      cal.ref_antenna=tile_use[ref_tile_use]
      cal.ref_antenna_name=(*obs.baseline_info).tile_names[cal.ref_antenna]
    ENDIF
    nan_i=where(Finite(vis_avg,/nan),n_nan)
    IF n_nan GT 0 THEN vis_model[nan_i]=(vis_avg[nan_i]=0)
    
    conv_test=fltarr(n_freq_use,max_cal_iter)
    ;    if not keyword_set(calib_freq_func) then begin
    ;vis_avg = sqrt(n_freq)*fft(vis_avg,dimension=1)
    ;vis_model = sqrt(n_freq)*fft(vis_model,dimension=1)
    ;gain_arr =sqrt(n_freq)* fft(gain_arr,dimension=1)
    
    ;vis_avg = fft(fft(vis_avg,dimension=1),dimension=1) ; double fourier transform
    ;vis_model = fft(vis_model,dimension=1)
    ;gain_arr =fft(gain_arr,dimension=1)
    
    
    gain_curr=Reform(gain_arr[*,tile_use])
    vis_data2=Reform(vis_avg[*,baseline_use]); & vis_data2=[[vis_data2],[Conj(vis_data2)]]
    vis_model2=Reform(vis_model[*,baseline_use]); & vis_model2=[[vis_model2],[Conj(vis_model2)]]
    weight2=Reform(weight[*,baseline_use]); & weight2=[[weight2],[weight2]]
    ;IF Keyword_Set(calibration_weights) THEN BEGIN baseline_wts2=Reform(baseline_weights[fi,baseline_use]) & baseline_wts2=[baseline_wts2,baseline_wts2] & ENDIF
    ;IF Keyword_Set(calibration_weights) THEN BEGIN baseline_wts2=Reform(baseline_weights[*,baseline_use]) & ENDIF
    
    ; Apply window function
    ;window = spectral_window(n_freq, type='Blackman-Harris',/periodic)
    ;norm_factor = sqrt(n_freq/total(window^2.))
    ;window = window * norm_factor
    ;window_expand = rebin(reform(window,n_freq,1), n_freq, 2*N_elements(baseline_use),  /sample)
    ;vis_data2 = vis_data2
    ;vis_model2 = vis_model2
    c = 299792458.
    
    uu = params.uu * c ; in meters
    vv = params.vv * c ; in meters
    
    z0_freq = 1420.40E6 ;; Hz
    freq = (*obs.baseline_info).freq
    redshifts = z0_freq/freq - 1 ;; frequencies will be identical if kx, ky, kz match
    ;mean_redshift = mean(redshifts)
    
    cosmology_measures, redshifts, wedge_factor = wedge_factor, comoving_dist_los=comov_dist_los, Ez=Ez
    
    mean_comoving_dist_los = mean(comov_dist_los)
    
    kx = 2. * !PI * uu / mean_comoving_dist_los
    ky = 2. * !PI * vv / mean_comoving_dist_los
    kperp = sqrt(abs(kx[baseline_use])^2. + abs(ky[baseline_use])^2.)
    
    ;*****calculate wedge
    ;; assume 20 degrees from pointing center to first null
    source_dist = 20d * !dpi / 180d
    fov_amp = wedge_factor * source_dist
    max_theta= 10d ;random guess!
    
    ;; calculate angular distance to horizon
    horizon_amp = wedge_factor * ((max_theta+90d) * !dpi / 180d)
    
    wedge_amp = [[fov_amp], [horizon_amp]] ; to be multiplied by k perp
    wedge_names = ['fov', 'horizon']
    ;*****
    
    comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
    comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
    
    z_mpc_delta = float(mean(comov_los_diff))
    z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
    kz_mpc_range =  (2.*!pi) / (z_mpc_delta)
    kz_mpc_delta = (2.*!pi) / z_mpc_length
    kz = findgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta - kz_mpc_range/2.
    folded_kz = abs(kz)
    
    ;Filtering
    kperp_cut = kperp#transpose(wedge_amp[*,0])
    kperp_mask = FLTARR(n_freq, N_elements(kperp))
    kperp_mask[*,*] = 1.
    
    cal_k_wedge=-1
    if keyword_set(cal_k_wedge) then begin
      if cal_k_wedge GT 0 then begin
        for freq_i=0, n_freq-1 do begin
          cut_freq_i = where(kperp_cut[*,freq_i] LT folded_kz[freq_i],n_count)
          if n_count GT 0 then kperp_mask[freq_i,cut_freq_i] = 0
        endfor
      endif else begin
        for freq_i=0, n_freq-1 do begin
          cut_freq_i = where(kperp_cut[*,freq_i] GT folded_kz[freq_i],n_count)
          if n_count GT 0 then kperp_mask[freq_i,cut_freq_i] = 0
        endfor
      endelse
    endif
    
    if keyword_set(cal_k_coverage) then begin
      coverage_cut = where(kperp GT .1,n_count)
      if n_count GT 0 then kperp_mask[*,coverage_cut] = 0
    endif
    
    ;weight2 = weight2 * [[kperp_mask],[kperp_mask]]
    kperp_mask = [[kperp_mask],[kperp_mask]]
    ;for freq_i=0, n_freq-1 do *kperp_mask_inds[freq_i] = where(kperp_mask GT 0,n_count)
    
    ;vis_data2 = fft(fft(vis_data2,dimension=1)*kperp_mask,dimension=1) ; double fourier transform
    ;vis_model2 = fft(vis_model2,dimension=1)*kperp_mask
    ;gain_curr =fft(gain_curr,dimension=1)
    ;     vis_data2 = fft(fft(vis_data2,dimension=1),dimension=1) ; double fourier transform
    ;vis_model2 = fft(vis_model2,dimension=1)
    ;gain_curr =fft(gain_curr,dimension=1)
    vis_data2 = fft(vis_data2,dimension=1) ; double fourier transform
    vis_model2 = fft(vis_model2,dimension=1)
    ;gain_curr =fft(gain_curr,dimension=1)
    
    ;b_i_use = PTRARR(n_freq,/allocate)
    ;for fi=0, n_freq-1 do $
    ;b_i_use_baseline=where(weight2[fi,*] GT 0,n_baseline_use2) & *b_i_use[fi]=b_i_use_baseline
    
    ;weight2=weight2[]
    ;vis_data2=vis_data2[b_i_use];*weight_invert(weight2)
    ;vis_model2=vis_model2[b_i_use];*weight_invert(weight2)
    
    A_ind=[tile_A_i_use];,tile_B_i_use]; & A_ind=A_ind[b_i_use]
    B_ind=[tile_B_i_use];,tile_A_i_use]; & B_ind=B_ind[b_i_use]
    ;window_expand2 = rebin(reform(window,n_freq,1), n_freq, N_elements(B_ind),  /sample)
    ;A_ind=tile_A_i_use; & A_ind=A_ind[b_i_use]
    ;B_ind=tile_B_i_use; & B_ind=B_ind[b_i_use]
    
    A_ind_arr=Ptrarr(n_tile_use,/allocate)
    n_arr=Fltarr(n_tile_use)
    FOR tile_i=0L,n_tile_use-1 DO BEGIN
      ;should be set up so that using where is okay
      inds=where(A_ind EQ tile_i,n1)
      print, n1
      IF n1 GT 1 THEN *A_ind_arr[tile_i]=Reform(inds,1,n1) ELSE *A_ind_arr[tile_i]=-1
      n_arr[tile_i]=n1 ;NEED SOMETHING MORE IN CASE INDIVIDUAL TILES ARE FLAGGED FOR ONLY A FEW FREQUENCIES!!
    ENDFOR
    
    
    vis_use=vis_data2
    
    FOR i=0L,(max_cal_iter-1)>1 DO BEGIN
      ;gain_curr_conj= fft(conj(fft(gain_curr,/inverse,dimension=1)),dimension=1)
      gain_curr_conj= conj(gain_curr)
      
      ;vis_model_matrix= 1/sqrt(n_freq)*fft( sqrt(n_freq)*fft(gain_curr[*,A_ind],dimension=1)*$
      ;  sqrt(n_freq)*fft(gain_curr_conj[*,B_ind],dimension=1)*sqrt(n_freq)*fft(vis_model2,dimension=1),/inverse, dimension=1)
      ;vis_model_matrix= 1/sqrt(n_freq)*fft( $
      ;  sqrt(n_freq)*fft(gain_curr_conj[*,B_ind],dimension=1)*sqrt(n_freq)*fft(vis_model2,dimension=1)$
      ;  ,/inverse, dimension=1)
      ;vis_model_matrix=fft(gain_curr_conj[*,B_ind],dimension=1)*fft(vis_model2,dimension=1)
      vis_model_matrix=gain_curr_conj[*,B_ind]*vis_model2
      gain_new=Complexarr(n_freq,n_tile_use)
      FOR fii=0L,n_freq_use-1 DO BEGIN
        fi=freq_use[fii]
        
        FOR tile_i=0L,n_tile_use-1 DO BEGIN
          IF n_arr[tile_i] GE min_cal_solutions THEN BEGIN
            b_i_use=where(weight2[fi,*A_ind_arr[tile_i]] GT 0,n_baseline_use2)
            if n_baseline_use2 GT 0 then $
              gain_new[fi,tile_i]=LA_Least_Squares( vis_model_matrix[fi,(*A_ind_arr[tile_i])[b_i_use]], vis_use[fi,(*A_ind_arr[tile_i])[b_i_use]], method=2, double=double_precision)
          ENDIF
        ENDFOR
        
      endfor
      
      ;gain_new = gain_new ;gain as a func of tau
      ;gain_new = fft(gain_new,dimension=1,/inverse)
      ;print, gain_new[0,0]
      gain_old=gain_curr
      
      ;IF phase_fit_iter-i GT 0 THEN gain_new*=Abs(gain_old)*weight_invert(Abs(gain_new)) ;fit only phase at first
      ;print, gain_curr[0,0]
      gain_curr=(gain_new+gain_old)/2.
      
      dgain=Abs(gain_curr)*weight_invert(Abs(gain_old))
      diverge_i=where(dgain LT Abs(gain_old)/2.,n_diverge)
      IF n_diverge GT 0 THEN gain_curr[diverge_i]=(gain_new[diverge_i]+gain_old[diverge_i]*2.)/3.
      IF nan_test(gain_curr) GT 0 THEN gain_curr[where(Finite(gain_curr,/nan))]=gain_old[where(Finite(gain_curr,/nan))]
      ;gain_curr*=Conj(gain_curr[ref_tile_use])/Abs(gain_curr[ref_tile_use])
      ;for tile_i=0, n_tile_use-1 do gain_curr[*,tile_i]*=Conj(gain_curr[*,ref_tile_use])*weight_invert(Abs(gain_curr[*,ref_tile_use]))
      conv_test[*,i]=Max(Abs(gain_curr-gain_old)*weight_invert(Abs(gain_old)))
      ;IF i GT phase_fit_iter THEN IF conv_test[fii,i] LE conv_thresh THEN BREAK
      
      ;convergence[freq_use,tile_use]=Abs(gain_curr-gain_old[freq_use,tile_use])*weight_invert(Abs(gain_old[freq_use,tile_use]))
      ;convergence_phase[freq_use,tile_use]=(phunwrap(atan(gain_curr,/phase))-phunwrap(atan(gain_old[freq_use,tile_use],/phase)));*weight_invert(phunwrap(Atan(gain_old,/phase)))
      convergence_plot[i]=mean(convergence)
      convergence_plot_phase[i]=mean(convergence_phase)
      ;if convergence_plot_phase[i] GT .4 then stop
      ;Ptr_free,A_ind_arr
      print, gain_curr[0,0]
      gain_arr[*,tile_use]=gain_curr
      
    ENDFOR
    stop
    
    Ptr_free,A_ind_arr
    
    ;gain_arr_freq = fft(gain_arr,dimension=1,/inverse)/384.
    
    
    n_freq=cal.n_freq
    tile_A_i=cal.tile_A-1 ;tile numbering starts at 1
    tile_B_i=cal.tile_B-1 ;tile numbering starts at 1
    n_baselines=Long(N_Elements(tile_A_i))
    
    inds_A=Rebin(Lindgen(n_freq),n_freq,n_baselines,/sample)+Rebin(transpose(tile_A_i)*n_freq,n_freq,n_baselines)
    inds_B=Rebin(Lindgen(n_freq),n_freq,n_baselines,/sample)+Rebin(transpose(tile_B_i)*n_freq,n_freq,n_baselines)
    vis_gain=gain_arr[inds_A]*conj(gain_arr[inds_B])
    vis_cal = fft(fft(*vis_ptr[0],dimension=1)*weight_invert(vis_gain),dimension=1,/inverse)
    
    
    ;vis_XX = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/vis_data/1061316296_vis_XX.sav', 'vis_ptr')
    vis_XX = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/vis_data/1061316176_vis_XX.sav', 'vis_ptr')
    gain_arr_old = *cal.gain[0]+*cal.gain_residual[0]
    vis_gain_old=gain_arr_old[inds_A]*conj(gain_arr_old[inds_B])
    vis_cal_old=*vis_ptr[0]*weight_invert(vis_gain_old)
    
    vis_gains_freq = vis_cal_old*weight_invert(vis_cal)*vis_gain_old
    gains_freq = (*cal.gain[0])
    
    
     A_ind_arr=Ptrarr(n_tile_use,/allocate)
    n_arr=Fltarr(n_tile_use)
    FOR tile_i=0L,n_tile_use-1 DO BEGIN
      ;should be set up so that using where is okay
      inds=where(tile_A_i_use EQ tile_i,n1)
      IF n1 GT 1 THEN *A_ind_arr[tile_i]=Reform(inds,1,n1) ELSE *A_ind_arr[tile_i]=-1
      n_arr[tile_i]=n1 ;NEED SOMETHING MORE IN CASE INDIVIDUAL TILES ARE FLAGGED FOR ONLY A FEW FREQUENCIES!!
    ENDFOR
    weight2=Reform(weight[*,baseline_use])
    
    
    gain_B = complex(FLTARR(n_freq,n_tile_use))
    gain_A = complex(FLTARR(n_freq,n_tile_use))
    gain_B[*,*]=1.
    gain_curr[*,*]=1.
    
    FOR i=0L,(max_cal_iter-1)>1 DO BEGIN
    
      for tile_i=0, n_tile_use-1 do begin
        
        
        ;inds1 = where(cal.tile_a EQ tile_i,n_count1)
        ;inds2 = where(cal.tile_b[inds1] EQ tile_i,n_count2)
        vis_gain_update = gain_curr[inds_A]*conj(gain_curr[inds_B])
      
        ;if n_count2 GT 0 then gains_freq[*,tile_i] = sqrt(mean(gains_freq_full[*,inds1[inds2]], dimension=2))
        FOR fii=0L,n_freq_use-1 DO BEGIN
          fi=freq_use[fii]
          
          b_i_use=where(weight2[fi,*A_ind_arr[tile_i]] GT 0,n_baseline_use2)
          if n_baseline_use2 GT 0 then $
            gain_A[fi,tile_i]=LA_Least_Squares( vis_gains_freq[fi,(*A_ind_arr[tile_i])[b_i_use]], vis_gain_update[fi,(*A_ind_arr[tile_i])[b_i_use]], method=2, double=double_precision)
        endfor
        
        gain_old=gain_curr
      gain_curr=(gain_A+gain_old)/2.
        
      endfor
    endfor
    stop
    
    ;    endif else begin
    ;      ;; calibrate using polynomials in frequency rather than fitting each frequency independently
    ;
    ;      common poly_cal_data, model_vis, data_vis, poly_cal_gain_B, f_arr, mode_type_use, mode_num_use, mask
    ;
    ;      n_phase_modes = 2
    ;      n_amp_modes = 3
    ;      n_mode = n_phase_modes + n_amp_modes
    ;      if n_phase_modes gt 0 then begin
    ;        mode_type = [strarr(n_phase_modes) + 'phase', strarr(n_amp_modes) + 'amp']
    ;        mode_num = [indgen(n_phase_modes), indgen(n_amp_modes)]
    ;      endif else begin
    ;        mode_type = strarr(n_amp_modes) + 'amp'
    ;        mode_num = indgen(n_amp_modes)
    ;      endelse
    ;
    ;      amp_modes = where(mode_type eq 'amp')
    ;      phase_modes = where(mode_type eq 'phase')
    ;
    ;      ;; keep frequency direction around
    ;      vis_data1=vis_avg[*,baseline_use]
    ;      vis_model1=vis_model[*,baseline_use]
    ;      weight1=weight[*,baseline_use]
    ;
    ;      hist_A = histogram(tile_A_i_use, min=0, max = n_tile-1, reverse_indices = ri_A)
    ;      hist_B = histogram(tile_B_i_use, min=0, max = n_tile-1, reverse_indices = ri_B)
    ;      for i=0, n_tile-1 do begin
    ;        if hist_A[i] ne 0 then begin
    ;          inds = ri_A[ri_A[i]:ri_A[i+1]-1]
    ;          if hist_B[i] ne 0 then inds = [inds, ri_B[ri_B[i]:ri_B[i+1]-1]]
    ;        endif else if hist_B[i] ne 0 then inds = ri_B[ri_B[i]:ri_B[i+1]-1] else continue
    ;
    ;        if i eq 0 then vis_tile_inds = create_struct('t'+number_formatter(i), inds) $
    ;        else vis_tile_inds = create_struct(vis_tile_inds, 't'+number_formatter(i), inds)
    ;      endfor
    ;
    ;      tile_freq_flag = bytarr(n_freq, n_tile_use)
    ;      n_arr=Fltarr(n_tile_use)
    ;      FOR tile_i=0L,n_tile_use-1 DO BEGIN
    ;        inds=vis_tile_inds.(tile_i)
    ;        IF n_elements(inds) GT 1 THEN begin
    ;          tile_freq_flag[*,tile_i] = total(weight1[*, inds],2) gt 0
    ;        endif
    ;        n_arr[tile_i]=n_elements(inds) ;NEED SOMETHING MORE IN CASE INDIVIDUAL TILES ARE FLAGGED FOR ONLY A FEW FREQUENCIES!!
    ;      ENDFOR
    ;      weight_mask = weight1 gt 0
    ;
    ;      ;; get starting point by fitting gain_arr
    ;      gain_arr_use = gain_arr[*,tile_use]
    ;      temp = dblarr(n_mode, n_tile_use)
    ;      f_arr = (dindgen(n_freq)/(n_freq-1))*2-1
    ;      for tile_i=0L, n_tile_use-1 do begin
    ;        wh_f_use = where(tile_freq_flag[*,tile_i] gt 0, count_f_use)
    ;        if count_f_use eq 0 then continue
    ;        temp[amp_modes,tile_i] = poly_fit(f_arr[wh_f_use], abs(gain_arr_use[[wh_f_use],tile_i]), n_amp_modes-1)
    ;        temp[phase_modes,tile_i] = poly_fit(f_arr[wh_f_use], phunwrap(atan(gain_arr_use[[wh_f_use],tile_i],/phase)), n_phase_modes-1)
    ;      end
    ;
    ;      gain_arr_mode = temp
    ;      ;; get legendre polynomial coefficients
    ;      if n_amp_modes gt 2 then begin
    ;        gain_arr_mode[amp_modes[2],*] = temp[amp_modes[2],*] * 2/3
    ;        gain_arr_mode[amp_modes[0],*] = temp[amp_modes[2],*] * 1/3 + temp[amp_modes[0],*]
    ;      endif
    ;
    ;      gain_fit_mode=gain_arr_mode
    ;      fit_phases = reform(gain_fit_mode[phase_modes[0],*])
    ;      if max(fit_phases) gt !pi then fit_phases[where(fit_phases gt !pi)] -= 2*!pi
    ;      if min(fit_phases) lt -!pi then fit_phases[where(fit_phases lt -!pi)] += 2*!pi
    ;      gain_fit_mode[phase_modes[0],*] = fit_phases
    ;      gain_fit = poly_cal_gain(gain_fit_mode, mode_type, mode_num, tile_freq_flag)
    ;
    ;      gain_curr_mode = gain_fit_mode
    ;
    ;      ;;start from gain=1.
    ;      ;gain_curr_mode = dblarr(n_mode, n_tile_use)
    ;      ;gain_curr_mode[amp_modes[0],*] = 1.
    ;
    ;      ;; double arrays to include all visibilities twice -- once as A,B & once as B,A
    ;      vis_data2=[[vis_data1],[Conj(vis_data1)]]
    ;      vis_model2=[[vis_model1],[Conj(vis_model1)]]
    ;      weight2=[[weight1],[weight1]]
    ;      A_ind=[tile_A_i_use,tile_B_i_use]
    ;      B_ind=[tile_B_i_use,tile_A_i_use]
    ;
    ;      gain_track = dblarr(n_mode, n_tile_use, max_cal_iter)
    ;      ncalls_track = lonarr(n_tile_use, max_cal_iter)
    ;
    ;      chi2_track = fltarr(n_tile_use, max_cal_iter)
    ;      conv_test=fltarr(n_mode, max_cal_iter)
    ;      conv_test2=fltarr(max_cal_iter)
    ;
    ;      loop_times = fltarr(max_cal_iter)
    ;      n_vis_use = lonarr(n_tile_use)
    ;
    ;      phase_fit_iter = 4
    ;      phase_0_iter = 2
    ;
    ;      time0 = systime(1)
    ;      FOR i=0L,(max_cal_iter-1)>1 DO BEGIN
    ;        loop_t0 = systime(1)
    ;
    ;
    ;        gain_new_mode=dblarr(n_mode, n_tile_use)
    ;        gain_track[*, *, i] = gain_curr_mode
    ;
    ;        gain_curr = poly_cal_gain(gain_curr_mode, mode_type, mode_num, tile_freq_flag)
    ;
    ;        if i lt phase_fit_iter then begin
    ;          if i lt phase_0_iter then begin
    ;            mode_type_use = mode_type[phase_modes[0], *]
    ;            mode_num_use = mode_num[phase_modes[0], *]
    ;          endif else begin
    ;            mode_type_use = mode_type[phase_modes, *]
    ;            mode_num_use = mode_num[phase_modes, *]
    ;          endelse
    ;        endif else begin
    ;          mode_type_use = mode_type
    ;          mode_num_use = mode_num
    ;        endelse
    ;
    ;        FOR tile_i=0L,n_tile_use-1 DO begin
    ;          wh_A = where(A_ind eq tile_i, n_vis_A)
    ;
    ;          if n_vis_A eq 0 then continue
    ;
    ;          if i eq 0 then n_vis_use[tile_i] = n_vis_A
    ;
    ;          f_arr = rebin((dindgen(n_freq)/(n_freq-1))*2-1, n_freq, n_vis_A)
    ;          model_vis = vis_model2[*,wh_A]
    ;          data_vis = vis_data2[*, wh_A]
    ;          poly_cal_gain_B = gain_curr[*,B_ind[wh_A]]
    ;          mask = weight_mask[*, wh_A]
    ;
    ;          if i lt phase_0_iter then p=gain_curr_mode[phase_modes[0], tile_i] else $
    ;            if i lt phase_fit_iter then p=gain_curr_mode[phase_modes, tile_i] else p=gain_curr_mode[*, tile_i]
    ;
    ;          ;p2=p
    ;          dfpmin, p, 1.0e-3, fval, 'poly_cal_chi2', 'poly_cal_grad', iter=ncalls, itmax=5000, stepmax = !pi/180.
    ;
    ;          ;dfpmin2, p2, 1.0e-3, fval2, 'poly_cal_chi2', 'poly_cal_grad', iter=ncalls2, itmax=50, stepmax = !pi/180.
    ;
    ;          ;p=p2
    ;          ;ncalls=ncalls2
    ;          ;stop
    ;          if i lt phase_fit_iter then begin
    ;            if i lt phase_0_iter then begin
    ;              gain_new_mode[*, tile_i]=gain_curr_mode[*,tile_i]
    ;              gain_new_mode[phase_modes[0],tile_i] = p[phase_modes[0]]
    ;              chi2_track[tile_i, i] = poly_cal_chi2(gain_new_mode[phase_modes[0], tile_i])
    ;            endif else begin
    ;              gain_new_mode[phase_modes, tile_i]=p
    ;              gain_new_mode[amp_modes,tile_i] = gain_curr_mode[amp_modes,tile_i]
    ;              chi2_track[tile_i, i] = poly_cal_chi2(gain_new_mode[phase_modes, tile_i])
    ;            endelse
    ;          endif else begin
    ;            gain_new_mode[*, tile_i]=p
    ;            chi2_track[tile_i, i] = poly_cal_chi2(gain_new_mode[*, tile_i])
    ;          endelse
    ;
    ;          ncalls_track[tile_i, i] = ncalls
    ;        endfor
    ;
    ;        IF Total(Abs(gain_new_mode[amp_modes])) EQ 0 THEN BEGIN
    ;          gain_curr_mode=gain_new_mode
    ;          BREAK
    ;        ENDIF
    ;
    ;        gain_old_mode=gain_curr_mode
    ;        gain_curr_mode=(gain_new_mode+gain_old_mode)/2.
    ;
    ;        IF i lt phase_fit_iter then gain_curr_mode[amp_modes,*] = gain_old_mode[amp_modes,*]
    ;
    ;        if min(gain_curr_mode[amp_modes[0],*]) lt 0 then $
    ;          gain_curr_mode[amp_modes[0],*] = abs(gain_curr_mode[amp_modes[0],*])
    ;
    ;        gain_curr_mode[phase_modes, *] -= rebin(gain_curr_mode[phase_modes, ref_tile_use], n_phase_modes, n_tile_use)
    ;
    ;        fit_phases = reform(gain_curr_mode[phase_modes[0],*])
    ;        if max(fit_phases) gt !pi then fit_phases[where(fit_phases gt !pi)] -= 2*!pi
    ;        if min(fit_phases) lt -!pi then fit_phases[where(fit_phases lt -!pi)] += 2*!pi
    ;
    ;        if max(fit_phases) gt !pi or min(fit_phases) lt -!pi then stop
    ;        gain_curr_mode[phase_modes[0],*] = fit_phases
    ;
    ;        gain_old = gain_curr
    ;        gain_curr = poly_cal_gain(gain_curr_mode, mode_type, mode_num, tile_freq_flag, phase=phase, amp=amp)
    ;
    ;        if min(amp) lt 0 then begin
    ;          neg_tile = where(amp lt 0) / n_freq
    ;          neg_tile = neg_tile[uniq(neg_tile)]
    ;
    ;          for tile_i=0, n_elements(neg_tile)-1 do gain_curr_mode[amp_modes, neg_tile[tile_i]] = gain_old_mode[amp_modes, neg_tile[tile_i]]
    ;
    ;          gain_curr = poly_cal_gain(gain_curr_mode, mode_type, mode_num, tile_freq_flag, phase=phase, amp=amp)
    ;
    ;        endif
    ;
    ;        small_amp = where(gain_curr_mode[amp_modes[0],*] lt 0.25, n_small)
    ;        ;small_amp = where(gain_curr_mode[amp_modes[0],*] lt 0.25 or min(abs(gain_curr)+abs(tile_freq_flag-1), dimension=1) eq 0, n_small)
    ;        if n_small gt 0 then begin
    ;          stop
    ;          gain_curr_mode[amp_modes[0],small_amp] = 1.
    ;          gain_curr_mode[phase_modes[0],small_amp] = 0.
    ;        endif
    ;
    ;        if i gt phase_fit_iter and n_small eq 0 then begin
    ;          large_chi2 = where(chi2_track[*,i] gt (mean(chi2_track[*,i]) + 2*stddev(chi2_track[*,i])), n_large_chi2)
    ;          if n_large_chi2 gt 0 then begin
    ;
    ;            if n_elements(kicked_tiles) gt 0 then begin
    ;              for large_chi2_i=0, n_large_chi2-1 do begin
    ;                wh_kick = where(kicked_tiles eq large_chi2[large_chi2_i], n_kick)
    ;                last_kick = kick_iter[wh_kick[n_kick-1]]
    ;                if n_kick eq 0 then begin
    ;                  gain_curr_mode[amp_modes[0],large_chi2[large_chi2_i]] = 1.
    ;                  gain_curr_mode[phase_modes[0],large_chi2[large_chi2_i]] = 0.
    ;
    ;                  kicked_tiles = [kicked_tiles, large_chi2[large_chi2_i]]
    ;                  kick_iter = [kick_iter, i]
    ;                endif else if n_kick lt 3 and (i-last_kick) gt 1 then begin
    ;                  ; regular kick didn't work, try a large phase slope
    ;                  gain_curr_mode[amp_modes[0],large_chi2[large_chi2_i]] = 1.
    ;                  gain_curr_mode[phase_modes[0],large_chi2[large_chi2_i]] = 0.
    ;                  gain_curr_mode[phase_modes[1],large_chi2[large_chi2_i]] = !pi * (-1)^(n_kick mod 2)
    ;
    ;                  kicked_tiles = [kicked_tiles, large_chi2[large_chi2_i]]
    ;                  kick_iter = [kick_iter, i]
    ;                endif
    ;              endfor
    ;            endif else begin
    ;              gain_curr_mode[amp_modes[0],large_chi2] = 1.
    ;              gain_curr_mode[phase_modes[0],large_chi2] = 0.
    ;              if n_elements(kicked_tiles) eq 0 then begin
    ;                kicked_tiles = large_chi2
    ;                kick_iter = intarr(n_large_chi2)+i
    ;              endif else begin
    ;                kicked_tiles = [kicked_tiles, large_chi2]
    ;                kick_iter = [kick_iter, intarr(n_large_chi2)+i]
    ;              endelse
    ;            endelse
    ;
    ;          endif
    ;        endif
    ;
    ;        conv_test[*, i]=Max(Abs(gain_curr_mode-gain_old_mode), dimension=2)
    ;        if i gt 0 then conv_test2[i]=Max((chi2_track[*,i-1]-chi2_track[*,i])/chi2_track[*,i-1])
    ;
    ;        loop_t1 = systime(1)
    ;        loop_times[i] = loop_t1-loop_t0
    ;        print, i, mean(ncalls_track[*, i]), loop_times[i]
    ;        n_iter = i+1
    ;
    ;        IF i GE phase_fit_iter THEN IF max(conv_test[*, i]) LE conv_thresh THEN BREAK
    ;
    ;      ENDFOR
    ;      time1 = systime(1)
    ;      print, 'polynomial cal loop time (m): ', (time1-time0)/60.
    ;
    ;      if n_iter lt max_cal_iter then begin
    ;        gain_track = gain_track[*,*,0:n_iter-1]
    ;        ncalls_track = ncalls_track[*,0:n_iter-1]
    ;
    ;        conv_test=conv_test[*,0:n_iter-1]
    ;        conv_test2=conv_test2[0:n_iter-1]
    ;        chi2_track=chi2_track[*,0:n_iter-1]
    ;
    ;        loop_times = loop_times[0:n_iter-1]
    ;      endif
    ;
    ;      undefine, kicked_tiles
    ;
    ;      gain_arr_mode[*,tile_use]=gain_curr_mode
    ;
    ;      gain_arr[*,tile_use] = gain_curr
    ;
    ;    endelse
    
    nan_i=where(Finite(gain_arr,/nan),n_nan)
    IF n_nan GT 0 THEN BEGIN
      ;any gains with NANs -> all tiles for that freq will have NANs
      freq_nan_i=nan_i mod n_freq
      freq_nan_i=freq_nan_i[Uniq(freq_nan_i,Sort(freq_nan_i))]
      (*vis_weight_ptr_use[pol_i])[freq_nan_i,*]=0
      gain_arr[nan_i]=0.
      
    ENDIF
    *cal_return.gain[pol_i]=gain_arr
    cal_return.convergence[pol_i]=Ptr_new(convergence)
  ENDFOR
  
  vis_count_i=where(*vis_weight_ptr_use[0],n_vis_cal)
  cal_return.n_vis_cal=n_vis_cal
  
  RETURN,cal_return
END
