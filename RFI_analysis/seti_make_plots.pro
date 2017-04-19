pro seti_make_plots, plot_data_file=plot_data_file

  ;Restore data for plotting
  if not keyword_set(plot_data_file) then plot_data_file = '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/center_freq_pointing/individual/seti_binned_diff_thesis_evenodd2_0_minustwo.sav'
  binned_diff_old = getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/center_freq_pointing/individual/seti_binned_diff_thesis_evenodd2_7_minustwo.sav','binned_diff')
  
  restore, plot_data_file
  If keyword_set(binned_diff_total) then binned_diff = binned_diff_total
  ;restore,'/nfs/mwa-00/h1/nbarry/seti_all_col_40sig_2000.sav'
  ;restore,'/nfs/mwa-00/h1/nbarry/seti_all_row_40sig_2000.sav'
  
  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  obs = GETVAR_SAVEFILE(dir+'metadata/1061316296_obs.sav', 'obs')
  freq_arr=(*obs.baseline_info).freq/10.^6.
  
  area=1E11+1.179E9
  s=28
  ;y_arr_dist=area*x_arr_dist*exp(-x_arr_dist^2./(2.*s^2.))/s^2.
  x_arr_dist=FINDGEN(N_elements(binned_diff))+.5
  
  ;>>> popt, pcov = optimize.curve_fit(f, x[3:200], y[3:200]/100000000)
  ;>>> print popt
  ;[ 988.37832724   28.16802621]
  
  ;>>> popt, pcov = optimize.curve_fit(f, x[3:150], y[3:150]/1000000)
  ;>>> print popt
  ;[  1.02816564e+05   2.84462205e+01]
  
  ;evenodd2
  ;5229.98088072    38.89738168
  
  ;evenodd2 - minustwo
  ;[ 1217.50858531*10000000    39.5904041 ]
  ;evenodd2 - minusone
  ;[ 1262.57881436*10000000    38.35548915]
  ;evenodd2 - zenith
  ;[ 1114.74515471*10000000    37.48823737]
  ;evenodd2 - plusone
  ;[ 728.23738851*10000000   36.90903868]
  ;evenodd2 - plusone
  ;[ 688.71176183*10000000    37.37008034]
  
  ;total, not split, no even-odd
  ;[ 10177.61541761*10000000     28.04886802]
  
  ;area=10177.61541761*10000000
  ;s=28.04886802
  area=53015.49709256*100000.
  s=54.08900427
  y_arr_dist=area*x_arr_dist*exp(-x_arr_dist^2./(2.*s^2.))/s^2.
  
  
  binned_diff=Long64(binned_diff)
  neg_values=where(binned_diff LT 0, n_count)
  binned_diff[neg_values]=binned_diff[neg_values]+2147483648*2
  x_arr_full=[0,FINDGEN(N_elements(binned_diff))+.5,N_elements(binned_diff)-.5]
  y_arr_full=[binned_diff[0],binned_diff,binned_diff[N_elements(binned_diff)-1]]
  
  binned_diff_old=Long64(binned_diff_old)
  neg_values=where(binned_diff_old LT 0, n_count)
  binned_diff_old[neg_values]=binned_diff_old[neg_values]+2147483648*2
  x_arr_full_old=[0,FINDGEN(N_elements(binned_diff_old))+.5,N_elements(binned_diff_old)-.5]
  y_arr_full_old=[binned_diff_old[0],binned_diff_old,binned_diff_old[N_elements(binned_diff_old)-1]]
  
  ;cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/center_freq_pointing/plots/vis_evenodd_res_thesis_loglog_DC_plustwo.png',/quiet,/nomatch
  cgplot, x_arr_full,y_arr_full,  yrange=[1.,10^8.], xrange=[2,1000],xtitle='Residual Even-Odd Visibility Amplitude (Jy)', ytitle='Binned Result', title='-2 Fall 2013 Even-Odd Residual Visibilties, DC Channel 1', charsize=1.0,/ylog,/xlog
  ;, yrange=[1.,10^10.]
  ;cgoplot, x_arr_dist, y_arr_dist, color='blue'
  cgoplot, x_arr_full_old,y_arr_full_old, color='blue'
  ;cglegend, title=['Flagged','Unflagged'], color=['black','blue'],location=[.65,.80],charsize=1
  ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  stop
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/antenna_check/hot_baselines/plots/seti_res_vis_thesis_evenodd.png',/quiet,/nomatch
  cgplot, x_arr_full,y_arr_full, xrange=[1,200], xtitle='Residual Visibility Amplitude (Jy)', ytitle='Binned Result', title='-3 Fall 2013 Semester Residual Even-Odd Visibilties', charsize=1
  cgoplot, x_arr_dist, y_arr_dist, color='blue'
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  stop
  
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/seti_res_vis_log_40sig_2000.png',/quiet,/nomatch
  cgplot, x_arr_full,y_arr_full, xrange=[2,30000], xtitle='Residual Visibility Amplitude (Jy)', ytitle='Binned Result', title='Fall 2013 Semester (Obs 2000) Residual Visibilties', charsize=1, /ylog, /xlog
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  min_index = min(where(all_col NE 0))
  max_index = max(where(all_col NE 0))
  binsize=1.
  x_arr=[freq_arr[min_index],freq_arr[min_index],freq_arr[min_index:max_index]+.080*binsize/2,freq_arr[max_index]+.080*binsize,freq_arr[max_index]+.080*binsize]
  y_arr=[0,all_col[min_index], all_col[min_index:max_index], all_col[max_index],0]
  cgPS_Open,'/nfs/eor-00/h1/nbarry/seti_outlier_freq_40sig_2000.png',/quiet,/nomatch
  cgplot, x_arr, y_arr,  xtitle='Frequency (MHz)', ytitle='Binned Result', title='Fall 2013 Semester (Obs 2000) 40$\sigma$ Visibility Outliers',psym=10,charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  ;*******Movie style plots
  longrun_names_match
  restore, '~/vis_res/thesis/seti_all_3D_20_thesis_evenodd2_total.sav'
  ;quick_image, all_3D_20[*,*,70],title=obs_names[70,0] + ' ' +obs_names[70,1] + ' ' +obs_names[70,2], xtitle = 'frequency index', ytitle = 'time index' ;;;;;test
  for obs_i=0, 1028 do if max(all_3D_20_total[*,*,obs_i]) GT 0 then quick_image, all_3D_20_total[*,0:27,obs_i],title=obs_names[obs_i,0] + ' ' +obs_names[obs_i,1] + ' ' +obs_names[obs_i,2] $
    + ' even-odd thesis cut, start of deviation', xtitle = 'frequency index', ytitle = 'time index',charsize = 1.2, data_range=[0,10], $
    savefile = '/nfs/mwa-00/h1/nbarry/vis_res/thesis/vis_res_movie/evenodd2/'+string(obs_i,format='(I4.4)')
  ;*******
    
  minustwo_obs = where(obs_names[*,2] EQ 'minustwo')
  minusone_obs = where(obs_names[*,2] EQ 'minusone')
  zenith_obs = where(obs_names[*,2] EQ 'zenith')
  plusone_obs = where(obs_names[*,2] EQ 'plusone')
  plustwo_obs = where(obs_names[*,2] EQ 'plustwo')
  
  ;*******Timing across pointings
  longrun_names_match,obs_names=obs_names
  restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/seti_all_3D_20_thesis_evenodd2_total.sav'
  all_3D_20 = ULONG64(INTARR(384,56))
  for obs_i=0,1028 do if (where(obs_i EQ minustwo_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then all_3D_20 = all_3D_20 + all_3D_20_total[*,*,obs_i]
  ;quick_image, all_3D_20[*,*,70],title=obs_names[70,0] + ' ' +obs_names[70,1] + ' ' +obs_names[70,2], xtitle = 'frequency index', ytitle = 'time index' ;;;;;test
  quick_image, all_3D_20[*,0:27],freq_arr, INDGEN(28),title='-2 Fall 2013 Residual Even-Odd Visibilities', $
    xtitle = 'Frequency (MHz)', ytitle = 'Time index (Even-Odd 2 second intervals)',charsize = 1.5,$
    savefile = '/nfs/mwa-00/h1/nbarry/vis_res/thesis/vis_res_movie/evenodd2/'+string(obs_i,format='(I4.4)')
    
  ;*******Freq plots over obsids
  longrun_names_match,obs_names=obs_names, n_obs=n_obs_per_day
  n_obs_per_day = total(n_obs_per_day, 2) ;obs_i=20 for Oct31
  n_obs_inds = where(INTARR(N_elements(n_obs_per_day)) NE 20)
  n_obs_per_day_tot = ULONG64(INTARR(N_elements(n_obs_inds)))
  for obs_i=0, N_elements(n_obs_inds)-1 do n_obs_per_day_tot[obs_i] = total(n_obs_per_day[0:n_obs_inds[obs_i]])  
  n_obs = 1029
  all_3D_20 = ULONG64(INTARR(384,n_obs))
  for obs_i=0,n_obs-1 do if (obs_names[obs_i,1] NE 'Oct31') then all_3D_20[*,obs_i] = all_3D_20[*,obs_i] + total(all_3D_20_total_yy[*,*,obs_i],2)
  ;quick_image, all_3D_20[*,*,70],title=obs_names[70,0] + ' ' +obs_names[70,1] + ' ' +obs_names[70,2], xtitle = 'frequency index', ytitle = 'time index' ;;;;;test
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/newlongrun/freq_day_xx.png',/quiet,/nomatch
  quick_image, all_3D_20,freq_arr, INDGEN(n_obs),title='Fall 2013 Residual Even-Odd Visibilities', $
    xtitle = 'Frequency (MHz)', ytitle = 'Day index of longrun',charsize = 1.5, data_range=[0,25], cb_title = 'Counts'
  for obs_i=0,N_elements(n_obs_inds)-1 do $
    cgoplot, [0,383], [n_obs_per_day_tot[obs_i], n_obs_per_day_tot[obs_i]], Linestyle=2
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    
  ;*******Timing across pointings
  longrun_names_match,obs_names=obs_names
  restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/residual/seti_all_3D_20_thesis_residual_total.sav'
  all_3D_20 = ULONG64(INTARR(384,56))
  for obs_i=0,1028 do if (where(obs_i EQ minustwo_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then all_3D_20 = all_3D_20 + all_3D_20_total[*,*,obs_i]
  quick_image, all_3D_20[*,*],freq_arr, INDGEN(56),title='-2 Fall 2013 Residual Visibilities', $
    xtitle = 'Frequency (MHz)', ytitle = 'Time index (2 second intervals)',charsize = 1.5,$
    savefile = '/nfs/mwa-00/h1/nbarry/vis_res/thesis/residual/timing_pointing_-2.pdf'
    
  all_3D_20 = ULONG64(INTARR(384,56))
  for obs_i=0,1028 do if (where(obs_i EQ minusone_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then all_3D_20 = all_3D_20 + all_3D_20_total[*,*,obs_i]
  quick_image, all_3D_20[*,*],freq_arr, INDGEN(56),title='-1 Fall 2013 Residual Visibilities', $
    xtitle = 'Frequency (MHz)', ytitle = 'Time index (2 second intervals)',charsize = 1.5, $
    savefile = '/nfs/mwa-00/h1/nbarry/vis_res/thesis/residual/timing_pointing_-1.pdf'
    
  all_3D_20 = ULONG64(INTARR(384,56))
  for obs_i=0,1028 do if (where(obs_i EQ zenith_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then all_3D_20 = all_3D_20 + all_3D_20_total[*,*,obs_i]
  quick_image, all_3D_20[*,*],freq_arr, INDGEN(56),title='0 Fall 2013 Residual Visibilities', $
    xtitle = 'Frequency (MHz)', ytitle = 'Time index (2 second intervals)',charsize = 1.5, $
    savefile = '/nfs/mwa-00/h1/nbarry/vis_res/thesis/residual/timing_pointing_0.pdf'
    
  all_3D_20 = ULONG64(INTARR(384,56))
  for obs_i=0,1028 do if (where(obs_i EQ plusone_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then all_3D_20 = all_3D_20 + all_3D_20_total[*,*,obs_i]
  quick_image, all_3D_20[*,*],freq_arr, INDGEN(56),title='+1 Fall 2013 Residual Visibilities', $
    xtitle = 'Frequency (MHz)', ytitle = 'Time index (2 second intervals)',charsize = 1.5, $
    savefile = '/nfs/mwa-00/h1/nbarry/vis_res/thesis/residual/timing_pointing_+1.pdf'
    
  all_3D_20 = ULONG64(INTARR(384,56))
  for obs_i=0,1028 do if (where(obs_i EQ plustwo_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then all_3D_20 = all_3D_20 + all_3D_20_total[*,*,obs_i]
  quick_image, all_3D_20[*,*],freq_arr, INDGEN(56),title='+2 Fall 2013 Residual Visibilities', $
    xtitle = 'Frequency (MHz)', ytitle = 'Time index (2 second intervals)',charsize = 1.5, $
    savefile = '/nfs/mwa-00/h1/nbarry/vis_res/thesis/residual/timing_pointing_+2.pdf'
    
    
  ;*******Antenna plots
  longrun_names_match, obs_names=obs_names
  ;For missing data...(plot as white in quick_image)
  ;and to force little interpolation between white and real data pixels (IDL is weird)
  
  restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/seti_tiles20_thesis_evenodd2_total.sav'
  tiles20_total=shift(tiles20_total,1,0,0) ;OOPS, locations_tile_a is indexed from 0, tile_b is indexed from 0. Saved data in July 2016 is off by 1 bin
  n=129
  i = REBIN(LINDGEN(n), n, n)
  j = REBIN(TRANSPOSE(LINDGEN(n)), n, n)
  mask = rebin(((j LT i) * (-20000)),129,129,3000)
  tiles20_total += mask
  ;tiles20_total=congrid(tiles20_total,n*50,n*50) put into quick_image call to preserve obsid axis
  x_array=congrid(LINDGEN(n),n*50,n*50)
  y_array=x_array
  
  for obs_i=0, 1028 do if where(obs_i EQ minustwo_obs) NE -1 then quick_image, congrid(tiles20_total[*,*,obs_i],n*50,n*50),x_array,y_array,title=obs_names[obs_i,0] + ' ' +obs_names[obs_i,1] + ' ' +obs_names[obs_i,2] $
    + ' even-odd thesis cut', xtitle = 'Antenna b', ytitle = 'Antenna a',charsize = 1.5, data_range=[0,30], missing_value=-20000, $
    savefile = '/nfs/mwa-00/h1/nbarry/vis_res/thesis/antenna/vis_res_movie/'+string(obs_i,format='(I4.4)')
    
  for obs_i=0, 1028 do if where(obs_i EQ minustwo_obs) NE -1 then quick_image, tiles40_total[*,*,obs_i],title=obs_names[obs_i,0] + ' ' +obs_names[obs_i,1] + ' ' +obs_names[obs_i,2] $
    + ' even-odd thesis cut', xtitle = 'Antenna b', ytitle = 'Antenna a',charsize = 1.2, data_range=[0,30], $
    savefile = '/nfs/mwa-00/h1/nbarry/vis_res/thesis/antenna/vis_res_movie40/'+string(obs_i,format='(I4.4)')
    
  restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/minusthree/no_cal_hot_removed/seti_tiles20_thesis_evenodd_total.sav'
  ;tiles20_total=shift(tiles20_total,1,0,0)
  n=129
  i = REBIN(LINDGEN(n), n, n)
  j = REBIN(TRANSPOSE(LINDGEN(n)), n, n)
  mask = rebin(((j LT i) * (-20000)),129,129)
  tiles_20_totalobs = ULONG64(INTARR(129,129))
  for obs_i=0, 1028 do if (obs_names[obs_i,1] NE 'Oct31') then tiles_20_totalobs = tiles_20_totalobs +tiles20_total[*,*,obs_i]
  tiles_20_totalobs = tiles_20_totalobs + mask
  quick_image, congrid(tiles_20_totalobs[*,*],n*50,n*50),x_array,y_array,data_range=[0,3000], missing_value=-20000, xtitle = 'Antenna b', ytitle = 'Antenna a', title='-3 Flagged, Cal removed'
  ;*******
  
  ;*******UV plots
  restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/antenna_check/hot_baselines/seti_tiles20_thesis_evenodd_total.sav'
  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  params = GETVAR_SAVEFILE(dir+'metadata/1061316296_params.sav', 'params')
  obs = GETVAR_SAVEFILE(dir+'metadata/1061316296_obs.sav', 'obs')
  uu = params.uu
  vv = params.vv
  tile_a = (*obs.baseline_info).tile_a
  tile_b = (*obs.baseline_info).tile_b
  
  hot = getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/antenna_check/hot_baselines_3000.sav','hot')
  hot_col = hot mod 129 ;tileA
  hot_row = hot / 129 ;tileB
  
  vis_mask = FLTARR(N_elements(tile_a))
  for hot_i=0, N_elements(hot)-1 do begin
    ;print, 'hot = ' + strtrim(hot_i,2)
    where_col = where(tile_a EQ hot_col[hot_i])
    where_row = where(tile_b EQ hot_row[hot_i])
    vis_mask[where_col[where_row]] = 1.
  endfor
  
  uu_masked = uu * vis_mask
  vv_masked = vv * vis_mask
  
  n_bins = 1000.
  uu_range = minmax(uu)
  vv_range = minmax(vv)
  
  pixels = FLTARR(n_bins,n_bins)
  binsize = (uu_range[1]-uu_range[0])/n_bins
  uu_result = histogram(uu_masked, binsize=binsize, min=uu_range[0], reverse_indices=ri_uu)
  ;  n_bins_dec = Ceil((dec_range_sculptor[0]-dec_range_sculptor[1])/binsize)
  
  for bin_i=0,n_bins-2 do begin
    if ri_uu[bin_i] EQ ri_uu[bin_i+1] then continue
    uu_inds = ri_uu[ri_uu[bin_i]:ri_uu[bin_i+1]-1]
    vv_result = histogram(vv_masked[uu_inds], binsize=binsize, min=vv_range_sculptor[0], reverse_indices=ri_vv)
    for bin_j=0,N_elements(vv_result)-2 do begin
      if ri_vv[bin_j] EQ ri_vv[bin_j+1] then continue
      vv_inds = ri_vv[ri_vv[bin_j]:ri_vv[bin_j+1]-1]
      pixels[bin_i,bin_j] = total(vv_result[uu_inds[vv_inds]])
    endfor
  endfor
  stop
  
  ;*******
  
  ;*******1D antenna plots
  longrun_names_match, obs_names=obs_names
  minustwo_obs = where(obs_names[*,2] EQ 'minustwo')
  minusone_obs = where(obs_names[*,2] EQ 'minusone')
  zenith_obs = where(obs_names[*,2] EQ 'zenith')
  plusone_obs = where(obs_names[*,2] EQ 'plusone')
  plustwo_obs = where(obs_names[*,2] EQ 'plustwo')
  tiles20_final = ULONG64(INTARR(129))
  for obs_i=0,1028 do if (where(obs_i EQ minustwo_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then for tile_i=1,128 do tiles20_final[tile_i] = tiles20_final[tile_i] $
    + total(tiles20_total[tile_i,*,obs_i]) + total(tiles20_total[*,tile_i,obs_i])
  cgplot, tiles20_final, ytitle='Counts', xtitle='Tile Number', title='-2 Fall 2013 Residual Even-Odd Visibilities',charsize=1.5
  tiles20_final = ULONG64(INTARR(129))
  for obs_i=0,1028 do if (where(obs_i EQ minusone_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then for tile_i=1,128 do tiles20_final[tile_i] = tiles20_final[tile_i] $
    + total(tiles20_total[tile_i,*,obs_i]) + total(tiles20_total[*,tile_i,obs_i])
  cgoplot, tiles20_final, ytitle='Counts', xtitle='Tile Number', title='-1 Fall 2013 Residual Even-Odd Visibilities',charsize=1.5, color='blue'
  for obs_i=0,1028 do if (where(obs_i EQ zenith_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then for tile_i=1,128 do tiles20_final[tile_i] = tiles20_final[tile_i] $
    + total(tiles20_total[tile_i,*,obs_i]) + total(tiles20_total[*,tile_i,obs_i])
  cgoplot, tiles20_final, ytitle='Counts', xtitle='Tile Number', title='0 Fall 2013 Residual Even-Odd Visibilities',charsize=1.5, color='green'
  for obs_i=0,1028 do if (where(obs_i EQ plusone_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then for tile_i=1,128 do tiles20_final[tile_i] = tiles20_final[tile_i] $
    + total(tiles20_total[tile_i,*,obs_i]) + total(tiles20_total[*,tile_i,obs_i])
  cgoplot, tiles20_final, ytitle='Counts', xtitle='Tile Number', title='+1 Fall 2013 Residual Even-Odd Visibilities',charsize=1.5, color='purple'
  for obs_i=0,1028 do if (where(obs_i EQ plustwo_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then for tile_i=1,128 do tiles20_final[tile_i] = tiles20_final[tile_i] $
    + total(tiles20_total[tile_i,*,obs_i]) + total(tiles20_total[*,tile_i,obs_i])
  cgoplot, tiles20_final, ytitle='Counts', xtitle='Tile Number', title='+2 Fall 2013 Residual Even-Odd Visibilities',charsize=1.5, color='yellow'
  ;*******
  
  
  ;*******1D freq plots
  longrun_names_match, obs_names=obs_names
  minustwo_obs = where(obs_names[*,2] EQ 'minustwo')
  all_3D_20_final = ULONG64(INTARR(384))
  for obs_i=0,1028 do if (where(obs_i EQ minustwo_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then for freq_i=0,383 do all_3D_20_final[freq_i] = all_3D_20_final[freq_i] $
    + total(all_3D_20_total[freq_i,*,obs_i])
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/remove_center_freq_pointing/plots/3D_20_evenodd2_freq_minustwo.png',/quiet,/nomatch
  cgplot, freq_arr,all_3D_20_final, ytitle='Counts', xtitle='Frequency (MHz)',xrange=[freq_arr[0],freq_arr[383]], title='-2 Fall 2013 Residual Even-Odd Visibilities, Center flagged',charsize=1.5
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  minusone_obs = where(obs_names[*,2] EQ 'minusone')
  all_3D_20_final = ULONG64(INTARR(384))
  for obs_i=0,1028 do if (where(obs_i EQ minusone_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then for freq_i=0,383 do all_3D_20_final[freq_i] = all_3D_20_final[freq_i] $
    + total(all_3D_20_total[freq_i,*,obs_i])
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/3D_20_evenodd2_freq_minusone.pdf',/quiet,/nomatch
  cgplot, freq_arr,all_3D_20_final, ytitle='Counts', xtitle='Frequency (MHz)',xrange=[freq_arr[0],freq_arr[383]], title='-1 Fall 2013 Residual Even-Odd Visibilities',charsize=1.5
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  zenith_obs = where(obs_names[*,2] EQ 'zenith')
  all_3D_20_final = ULONG64(INTARR(384))
  for obs_i=0,1028 do if (where(obs_i EQ zenith_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then for freq_i=0,383 do all_3D_20_final[freq_i] = all_3D_20_final[freq_i] $
    + total(all_3D_20_total[freq_i,*,obs_i])
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/3D_20_evenodd2_freq_zenith.pdf',/quiet,/nomatch
  cgplot, freq_arr,all_3D_20_final, ytitle='Counts', xtitle='Frequency (MHz)',xrange=[freq_arr[0],freq_arr[383]], title='0 Fall 2013 Residual Even-Odd Visibilities',charsize=1.5
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  plusone_obs = where(obs_names[*,2] EQ 'plusone')
  all_3D_20_final = ULONG64(INTARR(384))
  for obs_i=0,1028 do if (where(obs_i EQ plusone_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then for freq_i=0,383 do all_3D_20_final[freq_i] = all_3D_20_final[freq_i] $
    + total(all_3D_20_total[freq_i,*,obs_i])
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/3D_20_evenodd2_freq_plusone.pdf',/quiet,/nomatch
  cgplot, freq_arr,all_3D_20_final, ytitle='Counts', xtitle='Frequency (MHz)',xrange=[freq_arr[0],freq_arr[383]], title='+1 Fall 2013 Residual Even-Odd Visibilities',charsize=1.5
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  plustwo_obs = where(obs_names[*,2] EQ 'plustwo')
  all_3D_20_final = ULONG64(INTARR(384))
  for obs_i=0,1028 do if (where(obs_i EQ plustwo_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then for freq_i=0,383 do all_3D_20_final[freq_i] = all_3D_20_final[freq_i] $
    + total(all_3D_20_total[freq_i,*,obs_i])
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/3D_20_evenodd2_freq_plustwo.pdf',/quiet,/nomatch
  cgplot, freq_arr,all_3D_20_final, ytitle='Counts', xtitle='Frequency (MHz)',xrange=[freq_arr[0],freq_arr[383]], title='+2 Fall 2013 Residual Even-Odd Visibilities',charsize=1.5
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  ;*******
  
  ;*******1D freq plots
  restore, '/nfs/eor-00/h1/nbarry/vis_res/Aug23_decon/seti_all_3D_20_thesis_evenodd2.sav'
  longrun_names_match, obs_names=obs_names
  filename2='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23.txt'
  readcol, filename2, obs_ids, format='A', /silent
  match, obs_ids,obs_names[0:63,0],suba,subb

  all_3D_20_final = ULONG64(INTARR(384))
  for obs_i=0,1028 do if (obs_names[obs_i,1] NE 'Oct31') then for freq_i=0,383 do all_3D_20_final[freq_i] = all_3D_20_final[freq_i] $
    + total(all_3D_20_total[freq_i,*,obs_i])
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/minusthree/hot_removed/plots/3D_20_residual_freq.png',/quiet,/nomatch
  cgplot, freq_arr,all_3D_20_final, ytitle='Counts', xtitle='Frequency (MHz)',xrange=[freq_arr[0],freq_arr[383]], title='-3 Fall 2013 Residual Visibilities, Flagging',charsize=1.25
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  ;1D freq leftovers
  
end