pro seti_make_plots, plot_data_file=plot_data_file

  ;Restore data for plotting
  if not keyword_set(plot_data_file) then plot_data_file = '/nfs/mwa-00/h1/nbarry/vis_res/thesis/seti_binned_diff_thesis_evenodd2.sav'
  restore, plot_data_file
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
  
  area=1217.50858531*10000000
  s=39.5904041
  y_arr_dist=area*x_arr_dist*exp(-x_arr_dist^2./(2.*s^2.))/s^2.
  
  
  binned_diff=Long64(binned_diff)
  neg_values=where(binned_diff LT 0, n_count)
  binned_diff[neg_values]=binned_diff[neg_values]+2147483648*2
  x_arr_full=[0,FINDGEN(N_elements(binned_diff))+.5,N_elements(binned_diff)-.5]
  y_arr_full=[binned_diff[0],binned_diff,binned_diff[N_elements(binned_diff)-1]]
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_res_vis_thesis_log_evenodd2_minustwo_fit.pdf',/quiet,/nomatch
  cgplot, x_arr_full,y_arr_full, xrange=[1,2000], yrange=[1,10^10.], xtitle='Residual Visibility Amplitude (Jy)', ytitle='Binned Result', title='-2 Fall 2013 Residual Even-Odd Visibilties', charsize=1.5, /ylog
  cgoplot, x_arr_dist, y_arr_dist, color='blue'
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  stop
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_res_vis_thesis_evenodd2.png',/quiet,/nomatch
  cgplot, x_arr_full,y_arr_full, xrange=[1,200], xtitle='Residual Visibility Amplitude (Jy)', ytitle='Binned Result', title='Fall 2013 Semester Residual Even-Odd Visibilties', charsize=1
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
  
  ;*******Even-Odd timing across pointings
    longrun_names_match,obs_names=obs_names
  restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/seti_all_3D_20_thesis_evenodd2_total.sav'
  all_3D_20 = ULONG64(INTARR(384,56))
  for obs_i=0,1028 do if (where(obs_i EQ minustwo_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then all_3D_20 = all_3D_20 + all_3D_20_total[*,*,obs_i]
  ;quick_image, all_3D_20[*,*,70],title=obs_names[70,0] + ' ' +obs_names[70,1] + ' ' +obs_names[70,2], xtitle = 'frequency index', ytitle = 'time index' ;;;;;test
  quick_image, all_3D_20[*,0:27],freq_arr, INDGEN(28),title='-2 Fall 2013 Residual Even-Odd Visibilities', $
    xtitle = 'Frequency (MHz)', ytitle = 'Time index (Even-Odd 2 second intervals)',charsize = 1.5,$
    savefile = '/nfs/mwa-00/h1/nbarry/vis_res/thesis/vis_res_movie/evenodd2/'+string(obs_i,format='(I4.4)')
  
  
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
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/3D_20_evenodd2_freq_minustwo.pdf',/quiet,/nomatch
  cgplot, freq_arr,all_3D_20_final, ytitle='Counts', xtitle='Frequency (MHz)',xrange=[freq_arr[0],freq_arr[383]], title='-2 Fall 2013 Residual Even-Odd Visibilities',charsize=1.5
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
end