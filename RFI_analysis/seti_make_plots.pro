pro seti_make_plots

  restore,'/nfs/mwa-00/h1/nbarry/vis_res/thesis/seti_binned_diff_thesis_evenodd2.sav'
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

   area=5229.98088072*10000000
  s=38.89738168
    y_arr_dist=area*x_arr_dist*exp(-x_arr_dist^2./(2.*s^2.))/s^2.
  
  
  binned_diff=Long64(binned_diff)
  neg_values=where(binned_diff LT 0, n_count)
  binned_diff[neg_values]=binned_diff[neg_values]+2147483648*2
  x_arr_full=[0,FINDGEN(N_elements(binned_diff))+.5,N_elements(binned_diff)-.5]
  y_arr_full=[binned_diff[0],binned_diff,binned_diff[N_elements(binned_diff)-1]]
  stop
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_res_vis_thesis_log_evenodd2.png',/quiet,/nomatch
  cgplot, x_arr_full,y_arr_full, xrange=[1,2000], xtitle='Residual Visibility Amplitude (Jy)', ytitle='Binned Result', title='Fall 2013 Semester Residual Even-Odd Visibilties', charsize=1, /ylog
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
  
  ;Movie style plots
  longrun_names_match
  restore, '~/vis_res/thesis/seti_all_3D_20_thesis_evenodd2.sav'
  ;quick_image, all_3D_20[*,*,70],title=obs_names[70,0] + ' ' +obs_names[70,1] + ' ' +obs_names[70,2], xtitle = 'frequency index', ytitle = 'time index' ;;;;;test
  for obs_i=0, 1028 do if max(all_3D_20[*,*,obs_i]) GT 0 then quick_image, all_3D_20[*,*,obs_i],title=obs_names[obs_i,0] + ' ' +obs_names[obs_i,1] + ' ' +obs_names[obs_i,2] $
    + ' even-odd thesis cut, start of deviation', xtitle = 'frequency index', ytitle = 'time index',charsize = 1.2, data_range=[0,10], $
    savefile = '/nfs/mwa-00/h1/nbarry/vis_res/vis_res_movie/thesis_300/'+string(obs_i,format='(I4.4)')
  
  
  
end