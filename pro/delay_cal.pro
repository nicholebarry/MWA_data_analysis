pro delay_cal

  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_maxcalsources_nod_zenithpointing_notileflag/calibration/1061316176_cal.sav'
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_maxcalsources_nod_zenithpointing_notileflag/metadata/1061316176_obs.sav'
  freq_arr=(*obs.baseline_info).freq
  ;xaxis = freq_arr*1/30.64e6/80e3
  
  calxx=(*cal.gain[0])-1.
  fft_xx = FLTARR(384,128)
  fft_xx_fold = FLTARR(192,128)
  xaxis=FINDGEN(384)/30.64e6-384./30.64e6/2.
  
  for i=0,127 do fft_xx[*,i]=abs((fft(calxx[*,i],/center)))
  fft_xx[192,76]=0
  for i=0,191 do fft_xx_fold[i,*]=fft_xx[191+i,*]+fft_xx[191-i-1,*]
  
  fft_xx_sum=FLTARR(384)
  for i=0,383 do fft_xx_sum[i]=abs(TOTAL(fft_xx[i,*]))
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/delay_cal.png',/quiet,/nomatch
  
  cgplot, 1e9*xaxis[192:383],abs(fft_xx_fold[*,0]), yrange=[0,.003],charsize=1,/xlog,xrange=[6,3000], xtitle='Delay (ns)',ytitle='Amplitude of gain deviation modes', title='FT of calibration deviation frequency spectrum'
  for i=1, 127 do cgoplot,1e9*xaxis[192:383], fft_xx_fold[*,i], yrange=[0,.001],xrange=[6,3000]
  cgoplot,1e9*xaxis[192:383], mean(fft_xx_fold,dimension=2,/NAN), yrange=[0,.001],xrange=[6,3000],color='red'
  ;cgPS_Open,'/nfs/eor-00/h1/nbarry/delay_cal_sum.png',/quiet,/nomatch
  ;cgplot, 1e9/xaxis,abs(fft_xx_sum),/xlog,xrange=[6,3000],charsize=1, xtitle='Delay (ns)',ytitle='Amplitude of summed gain deviation modes',$
  ;  title='FT of summed calibration deviation frequency spectrum'
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  stop
  
end