pro pfb_edge_brightness_corr


  n_tile=128
  n_freq=384
  n_pol=2
  pol_name=['XX','YY']

  lambda = ['750','100','150','200','250']

  ;dir_path = ['/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_int/vis_data/',$
  ;dir_path = ['/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_int/vis_data/',$
  dir_path = ['/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_SSINS_base_flag_int/vis_data/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max100_int/vis_data/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max150_int/vis_data/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max200_int/vis_data/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max250_int/vis_data/']

  obs_path = ['/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_int/metadata/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max100_int/metadata/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max150_int/metadata/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max200_int/metadata/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max250_int/metadata/']

  n_dir = N_elements(dir_path)

  ;obs_ids = ['1061317400'];,'1061317520','1061317640','1061317888','1061318128','1061318248','1061318376','1061318496'];,'1061318736','1061318864']
  obs_ids = ['1061317400'];,'1061318128','1061318248'];,'1061318376','1061318496','1061318736','1061318864']
;obs_ids = ['1061312152','1061312272','1061312640'];,'1061312760','1061312880','1061313496','1061313616','1061313736',$
;'1061313856','1061313984','1061314104','1061314344','1061314472','1061314592','1061314712','1061315080',$
;'1061315320','1061315448','1061316544','1061316664','1061316784','1061316912','1061317032','1061317152']

  n_obs = N_Elements(obs_ids)

  rgbcolors=$
      [[228,134,161],$
;      [214,63,91],$
;      [164,79,86],$
;      [200,79,48],$
;      [217,144,103],$
      [219,133,46],$
;      [137,106,42],$
;      [198,167,62],$
;      [153,181,52],$
;      [159,178,106],$
      [82,115,55],$
;      [73,147,61],$
;      [95,196,87],$
;      [81,183,144],$
;      [70,174,206],$
      [100,138,212],$
;      [99,101,214],$
;      [101,87,160],$
;      [186,126,223],$
;      [194,150,217],$
      [164,79,201],$
;      [210,85,181],$
;      [156,79,136],$
      [209,67,131]]
  rgbcolors=reverse(rgbcolors)
  color_num = [10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56]

  for n_i=0, N_elements(rgbcolors[0,*])-1 do $
    TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]
  color_byte = [10B,12B,14B,16B,18B,20B,22B,24B,26B,28B,30B,32B,34B,36B,38B,40B,42B,44B,46B,48B,50B,52B,54B,56B]

  obs = getvar_savefile('/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_int/metadata/' + obs_ids[0] + '_obs.sav','obs')
  freq_use=(*obs.baseline_info).freq_use
  freq_use=where(freq_use,nf_use)

  subband_low =  where((findgen(384) mod 16) EQ 0) + 2 ;minus the edge
  subband_high =  where((findgen(385) mod 16) EQ 0) -3 ;minus the edge
  subband_high = subband_high[1:24]
  subband_low_edge =  where((findgen(384) mod 16) EQ 0) + 1
  subband_high_edge =  where((findgen(385) mod 16) EQ 0) - 2
  subband_high_edge = subband_high_edge[1:24]
  subband_dc = where((findgen(384) mod 16) EQ 0) + 8

  params = getvar_savefile('/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_int/metadata/' + obs_ids[0] + '_params.sav','params')


  n_vis = N_elements(params.uu)
  mean_per_band = FLTARR(n_dir, n_pol, n_obs, 24,n_vis)
  mean_per_edge = FLTARR(n_dir, n_pol, n_obs, 24,n_vis)

  for dir_i=0, 0 do begin

print, dir_path[dir_i] 

    for obs_i=0,n_obs-1 do begin

      restore, dir_path[dir_i] + obs_ids[obs_i] + '_flags.sav'

n_vis_current = N_elements((*vis_weights[0])[0,*])
if (n_vis_current NE n_vis) then begin
  obs_current = getvar_savefile(obs_path[dir_i] + obs_ids[obs_i] + '_obs.sav','obs')
  diff = (*obs_current.baseline_info).tile_b[0:127] - (*obs.baseline_info).tile_b[0:127]
  zeros = where(diff NE 0,n_count)
  defined_a = where((*obs.baseline_info).tile_b NE zeros[0]+1)
  defined_b = where((*obs.baseline_info).tile_a[defined_a] NE zeros[0]+1)
;for pol_i=0, obs.n_pol-1 do *vis_model_arr[pol_i] = (*vis_model_arr[pol_i])[*,defined_a[defined_b]]
endif

      for pol_i=0, n_pol-1 do begin

        inds = where((*vis_weights[pol_i]) LT 0,n_count)    
        if n_count GT 0 then (*vis_weights[pol_i])[inds] = 0

        ;48 flagged edge channels, 24 half flagged DC channels, .01 is a buffer
        inds_partial_flags = where ( mean(*vis_weights[pol_i],/nan,dimension=1) LT ((384 - 48 - 24*.5) / 384. * max(*vis_weights[pol_i],/nan,dimension=1) -.01), n_count)
        if n_count GT 0 then (*vis_weights[pol_i])[*,inds_partial_flags] = 0

        restore, dir_path[dir_i] + obs_ids[obs_i] + '_vis_'+pol_name[pol_i] + '.sav'        


        for subband_i=0,23 do begin

          freq_inds = INDGEN( subband_high[subband_i] - subband_low[subband_i] +1) + subband_low[subband_i] 
          freq_inds = [freq_inds[0:5],freq_inds[7:11]] ; minus the dc freq

          if (n_vis_current EQ n_vis) then begin
            ;vis_band = abs((*vis_ptr)[subband_low[subband_i]:subband_high[subband_i],*])*(*vis_weights[pol_i])[subband_low[subband_i]:subband_high[subband_i],*]
            vis_band = abs((*vis_ptr)[freq_inds,*])*(*vis_weights[pol_i])[freq_inds,*]
            inds = where(vis_band EQ 0,n_count)
            if n_count GT 0 then vis_band[inds] = !Values.F_NAN
            mean_per_band[dir_i,pol_i,obs_i,subband_i,*] = mean(vis_band,dimension=1,/nan)
            mean_per_edge[dir_i,pol_i,obs_i,subband_i,*] = ( abs((*vis_ptr)[subband_low_edge[subband_i],*])*(*vis_weights[pol_i])[subband_low_edge[subband_i],*] + abs((*vis_ptr)[subband_high_edge[subband_i],*])*(*vis_weights[pol_i])[subband_high_edge[subband_i],*] ) /2
          endif else begin
            if n_vis_current LT n_vis then begin
            vis_band = abs((*vis_ptr)[subband_low[subband_i]:subband_high[subband_i],*])*(*vis_weights[pol_i])[subband_low[subband_i]:subband_high[subband_i],*]
            inds = where(vis_band EQ 0,n_count)
            if n_count GT 0 then vis_band[inds] = !Values.F_NAN
            mean_per_band[dir_i,pol_i,obs_i,subband_i,defined_a[defined_b]] = mean(vis_band,dimension=1,/nan)
            mean_per_edge[dir_i,pol_i,obs_i,subband_i,defined_a[defined_b]] = ( abs((*vis_ptr)[subband_low_edge[subband_i],*])*(*vis_weights[pol_i])[subband_low_edge[subband_i],*] + abs((*vis_ptr)[subband_high_edge[subband_i],*])*(*vis_weights[pol_i])[subband_high_edge[subband_i],*] ) /2
            endif else begin
            ;an extra time step for no reason. cut off the end
            vis_band = abs((*vis_ptr)[subband_low[subband_i]:subband_high[subband_i],0:n_vis-1])*(*vis_weights[pol_i])[subband_low[subband_i]:subband_high[subband_i],0:n_vis-1]
mean_per_band[dir_i,pol_i,obs_i,subband_i,*] = mean(vis_band,dimension=1,/nan)
            mean_per_edge[dir_i,pol_i,obs_i,subband_i,*] = ( abs((*vis_ptr)[subband_low_edge[subband_i],0:n_vis-1])*(*vis_weights[pol_i])[subband_low_edge[subband_i],0:n_vis-1] + abs((*vis_ptr)[subband_high_edge[subband_i],0:n_vis-1])*(*vis_weights[pol_i])[subband_high_edge[subband_i],0:n_vis-1] ) /2
            endelse
          endelse



        endfor 

      endfor

    endfor

    inds = where(mean_per_band EQ 0,n_count)
    IF n_count GT 0 then mean_per_band[inds] = !Values.F_NAN
    inds = where(mean_per_edge EQ 0,n_count)
    IF n_count GT 0 then mean_per_edge[inds] = !Values.F_NAN
    mean_per_band_over_obs = mean(mean_per_band,dimension=3,/nan)
    mean_per_edge_over_obs = mean(mean_per_edge,dimension=3,/nan)

    inds = where(mean_per_band_over_obs EQ 0,n_count)
    IF n_count GT 0 then mean_per_band_over_obs[inds] = !Values.F_NAN
    inds = where(mean_per_edge_over_obs EQ 0,n_count)
    IF n_count GT 0 then mean_per_edge_over_obs[inds] = !Values.F_NAN

    mean_per_edge_over_obs_orig = mean_per_edge_over_obs
    mean_per_band_over_obs_orig = mean_per_band_over_obs

    mean_per_band_over_obs = reform(mean_per_band_over_obs[dir_i,*,0:15,*])
    mean_per_edge_over_obs = reform(mean_per_edge_over_obs[dir_i,*,0:15,*])
    ;mean_per_band_over_obs = mean(reform(mean_per_band_over_obs[dir_i,*,0:15,*]),dimension=2)
    ;mean_per_edge_over_obs = mean(reform(mean_per_edge_over_obs[dir_i,*,0:15,*]),dimension=2)

    ;Select baselines that go to the power spectrum
    lambda_cut = sqrt(params.uu^2 + params.vv^2)*182E6
    inds = where(lambda_cut GT 100, n_count)
    for pol_i=0,n_pol-1 do mean_per_band_over_obs[pol_i,*,inds] = !Values.F_NAN
    for pol_i=0,n_pol-1 do mean_per_edge_over_obs[pol_i,*,inds] = !Values.F_NAN

    n_bins = 50
    range = max(mean_per_band_over_obs,/nan)
    steps = (range) / n_bins

    inds_ptr=PTRARR(n_bins,/allocate)
    for step_i=0, n_bins-1 do *inds_ptr[step_i] = where(mean_per_band_over_obs LT steps*(step_i+1) AND mean_per_band_over_obs GT steps*(step_i))
    mean_per_band_amp = FLTARR(n_bins)
    for step_i=0, n_bins-1 do mean_per_band_amp[step_i] = mean(mean_per_band_over_obs[*inds_ptr[step_i]],/nan)
    mean_per_edge_amp = FLTARR(n_bins)
    for step_i=0, n_bins-1 do mean_per_edge_amp[step_i] = mean(mean_per_edge_over_obs[*inds_ptr[step_i]],/nan)
    

    result = histogram(mean_per_band_over_obs,nbins=n_bins,min=0,/nan,reverse_indices=ri)
    x_bins = FINDGEN(n_bins)*steps + steps
   ; result = [0,result]
  ;  x_bins = [0,x_bins]


    ;lambda = sqrt(params.uu^2 + params.vv^2)*182E6


;[leftAxisLoc, bottomAxisLoc, rightAxesLoc, topAxesLoc]
stop

max_x = 170
output_path = dir_path[dir_i] + 'vis_cal_edge_problem3_15_SSINS_'+lambda[dir_i]+'.png'
cgps_open,output_path,/quiet,/nomatch
print, output_path

small_inds = where(mean_per_band_amp[*] LT max_x)
cgplot, mean_per_band_amp[small_inds],mean_per_edge_amp[small_inds] / mean_per_band_amp[small_inds],psym=14,ERR_YHigh=1/sqrt(result[small_inds]),ERR_YLow=1/sqrt(result[small_inds]),Ystyle=8, $
  ytitle='Normalized calibrated visibility edge channel', charsize=1,Position=[0.10, 0.4, 0.9, 0.90],xrange=[0,max_x],yrange=[0,5],title='Pointing 1, cal <'+lambda[dir_i]+' $\lambda$'
temp = max(result,max_ind)
cgoplot, [x_bins[max_ind],x_bins[max_ind]],[0,20],color='blue',Position=[0.10, 0.4, 0.9, 0.90]
;cgplot, mean_per_edge_amp[0,*],mean_per_edge_amp[0,*] / mean_per_band_amp[0,*],psym=14,Ystyle=8
;result = [0,result]
; x_bins = [0,x_bins]
cgaxis, yaxis=1, yrange=[1,max(result) + 10000],/ylog, /save, title='# vis in bin', charsize=1, color='blue'
cgoplot, x_bins, result, psym=10,/ylog,color='blue'
small_inds = where(reform(mean_per_edge_amp / mean_per_band_amp) LT 1.6)

smaller_inds = where(mean_per_band_amp[small_inds] LT max_x)
cgplot, mean_per_band_amp[small_inds[smaller_inds]],mean_per_edge_amp[small_inds[smaller_inds]] / mean_per_band_amp[small_inds[smaller_inds]],$
  psym=14,ERR_YHigh=1/sqrt(result[small_inds[smaller_inds]]),ERR_YLow=1/sqrt(result[small_inds[smaller_inds]]),$
  charsize=1,xtitle='Avg calibrated visibility amplitude',/NoErase,Position=[0.10, 0.1, 0.9, 0.35],yrange=[.8,1.6],xrange=[0,max_x]
cgoplot, [x_bins[max_ind],x_bins[max_ind]],[0,20],color='blue',Position=[0.10, 0.1, 0.9, 0.35]

cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage




;    wavelength_dist =  mean((*obs.baseline_info).freq) * sqrt(params.uu^2. + params.vv^2.)

;    inds10 = where(wavelength_dist LT 10)
;    inds20 = where(wavelength_dist LT 20 AND wavelength_dist GT 10)
;    inds30 = where(wavelength_dist LT 30 AND wavelength_dist GT 20)
;    inds40 = where(wavelength_dist LT 40 AND wavelength_dist GT 30)
;    inds50 = where(wavelength_dist LT 50 AND wavelength_dist GT 40)
;    inds70 = where(wavelength_dist LT 70 AND wavelength_dist GT 50)
;    inds100 = where(wavelength_dist LT 100 AND wavelength_dist GT 70)
;    inds150 = where(wavelength_dist LT 150 AND wavelength_dist GT 100)
;    inds200 = where(wavelength_dist LT 200 AND wavelength_dist GT 150)
;    inds = where(wavelength_dist GT 200)

;    mean_per_ring = FLTARR(n_dir,n_pol,24,10)
;    mean_per_ring[*,*,*,0] = mean(mean_per_band_over_obs[*,*,*,inds10],/nan,dimension=4)
;    mean_per_ring[*,*,*,1] = mean(mean_per_band_over_obs[*,*,*,inds20],/nan,dimension=4)
;    mean_per_ring[*,*,*,2] = mean(mean_per_band_over_obs[*,*,*,inds30],/nan,dimension=4)
;    mean_per_ring[*,*,*,3] = mean(mean_per_band_over_obs[*,*,*,inds40],/nan,dimension=4)
;    mean_per_ring[*,*,*,4] = mean(mean_per_band_over_obs[*,*,*,inds50],/nan,dimension=4)
;    mean_per_ring[*,*,*,5] = mean(mean_per_band_over_obs[*,*,*,inds70],/nan,dimension=4)
;    mean_per_ring[*,*,*,6] = mean(mean_per_band_over_obs[*,*,*,inds100],/nan,dimension=4)
;    mean_per_ring[*,*,*,7] = mean(mean_per_band_over_obs[*,*,*,inds150],/nan,dimension=4)
;    mean_per_ring[*,*,*,8] = mean(mean_per_band_over_obs[*,*,*,inds200],/nan,dimension=4)
;    mean_per_ring[*,*,*,9] = mean(mean_per_band_over_obs[*,*,*,inds],/nan,dimension=4)

;    mean_per_edge_ring = FLTARR(n_dir,n_pol,24,10)
;    mean_per_edge_ring[*,*,*,0] = mean(mean_per_edge_over_obs[*,*,*,inds10],/nan,dimension=4)
;    mean_per_edge_ring[*,*,*,1] = mean(mean_per_edge_over_obs[*,*,*,inds20],/nan,dimension=4)
;    mean_per_edge_ring[*,*,*,2] = mean(mean_per_edge_over_obs[*,*,*,inds30],/nan,dimension=4)
;    mean_per_edge_ring[*,*,*,3] = mean(mean_per_edge_over_obs[*,*,*,inds40],/nan,dimension=4)
;    mean_per_edge_ring[*,*,*,4] = mean(mean_per_edge_over_obs[*,*,*,inds50],/nan,dimension=4)
;    mean_per_edge_ring[*,*,*,5] = mean(mean_per_edge_over_obs[*,*,*,inds70],/nan,dimension=4)
;    mean_per_edge_ring[*,*,*,6] = mean(mean_per_edge_over_obs[*,*,*,inds100],/nan,dimension=4)
;    mean_per_edge_ring[*,*,*,7] = mean(mean_per_edge_over_obs[*,*,*,inds150],/nan,dimension=4)
;    mean_per_edge_ring[*,*,*,8] = mean(mean_per_edge_over_obs[*,*,*,inds200],/nan,dimension=4)
;    mean_per_edge_ring[*,*,*,9] = mean(mean_per_edge_over_obs[*,*,*,inds],/nan,dimension=4)
;stop

  endfor

stop




end
