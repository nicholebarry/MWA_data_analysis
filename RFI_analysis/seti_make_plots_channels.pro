pro seti_make_plots_channels

  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  obs = GETVAR_SAVEFILE(dir+'metadata/1061316296_obs.sav', 'obs')
  freq_arr=(*obs.baseline_info).freq/10.^6.
  
  ;x_arr_dist=FINDGEN(N_elements(binned_diff))+.5
  ;area=5246.81437*100000.
  ;s=52.32255141
  ;y_arr_dist=area*x_arr_dist*exp(-x_arr_dist^2./(2.*s^2.))/s^2.
  
  ;cgplot, [1,2],[1,2],  yrange=[1.,10^6.], xrange=[100,1000],$
  ;  xtitle='Residual Even-Odd Visibility Amplitude (Jy)', ytitle='Binned Result', title='-2 Fall 2013 Even-Odd Residual Visibilties, DC Channels', charsize=1.0,/ylog,/xlog, /NoDATA
  
  rgbcolors=[[69,190,207],[144,23,3],[240,87,249],[45,165,30],[42,25,77],[1,53,8],[167,138,28],[251,193,236],[52,86,193],[246,229,194],[148,33,145],[253,88,89],[239,233,110],$
    [14,128,63],[96,50,14],[176,247,108],[74,126,157],[47,14,36],[212,120,250],[202,29,168],[113,166,89],[137,117,202],[197,120,62],[203,86,131]]
    
  y_array = LONG64(FLTARR(24,4098))
  x_array = LONG64(FLTARR(24,4098))
  
  for i=0,23 do begin
    binned_diff = getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/center_freq_pointing/individual/seti_binned_diff_thesis_evenodd2_'+strtrim(i,2)+'_minustwo.sav','binned_diff')
    binned_diff=Long64(binned_diff)
    neg_values=where(binned_diff LT 0, n_count)
    binned_diff[neg_values]=binned_diff[neg_values]+2147483648*2
    x_arr_full=[0,FINDGEN(N_elements(binned_diff))+.5,N_elements(binned_diff)-.5]
    x_array[i,0:N_elements(x_arr_full)-1] = x_arr_full
    y_arr_full=[binned_diff[0],binned_diff,binned_diff[N_elements(binned_diff)-1]]
       y_array[i,2:N_elements(y_arr_full)-10] = y_arr_full[2:N_elements(y_arr_full)-10]
    
  ;cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/center_freq_pointing/plots/vis_evenodd_res_thesis_loglog_DC_plustwo.png',/quiet,/nomatch
  ;    Device, Decomposed=0
  ;    TVLCT, rgbcolors[0,i], rgbcolors[1,i], rgbcolors[2,i],100
  ;    cgoplot, x_arr_full,y_arr_full,  yrange=[1.,10^6.], xrange=[2,1000],$
  ;      xtitle='Residual Even-Odd Visibility Amplitude (Jy)', ytitle='Binned Result', title='-2 Fall 2013 Even-Odd Residual Visibilties, DC Channel 1', $
  ;      charsize=1.0,/ylog,/xlog,color=100B
  ;    Device, Decomposed=1
  ;, yrange=[1.,10^10.]
  ;cgoplot, x_arr_dist, y_arr_dist, color='blue'
  ;cgoplot, x_arr_full_old,y_arr_full_old, color='blue'
  ;cglegend, title=['Flagged','Unflagged'], color=['black','blue'],location=[.65,.80],charsize=1
  ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  endfor
  quick_image, y_array[*,0:1000],/log,xtitle='DC channel index',ytitle='Jy', title = '-2 Fall 2013 Even-Odd Residual Visibilties, DC Channel', charsize=1, cb_title='Counts'
  stop
end