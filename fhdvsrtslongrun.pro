pro FHDvsRTSlongrun

  ;  obs_id=['1065448328','1065449184','1065879152','1065880008','1066568464','1066740792','1066741648','1067085448','1067086304','1067257776','1068809584']

  ;obs_id=['1062173424','1062174032','1062174160','1062175008','1062175376','1062175496','1062175624','1062176112','1062176472','1062176600', $
  ;  '1062176840','1062346360','1062518688','1063466496','1062949512','1063121960','1062347336','1062519664','1062950368','1063122696','1061313368', $
  ;  '1064588456','1064590536','1061317640','1063639072']


  obs_id= ['1061658880','1063470400','1065448328','1065449184','1065879152','1065880008','1066568464','1066740792','1066741040','1066741648','1067085448', $
    '1067086304','1067090576','1067257048','1067257776','1067257896','1067258752','1067259000','1067259120','1067259360','1067259848','1067260096','1067260336', $
    '1067260464','1068809584']
    
    
    
  cal_good = getvar_savefile('/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/1061312520_cal.sav','cal')
  
  cal_tot=FLTARR(384,128,N_elements(obs_id))
  
  for obs_i=0,N_elements(obs_id)-1 do begin
    cal_bad = getvar_savefile('/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/'+obs_id[obs_i]+'_cal.sav','cal')
    cal_divide=abs((*cal_bad.gain[0]))/abs((*cal_good.gain[0]))
    for freq_i=0,383 do for tile_i=0,127 do if ~finite(cal_divide[freq_i,tile_i]) then cal_divide[freq_i,tile_i]=0
    cal_divide_norm=cal_divide
    quick_image, cal_divide_norm, data_range=[0,1], title='Percent change of '+obs_id[obs_i]+' calibration amplitude, yy', xtitle='Freq channel index', ytitle='Tile index',window=obs_i, $
      savefile='/nfs/mwa-00/h1/nbarry/Videos/'+string(obs_i,format='(I2.2)')
      
  ;cal_tot[*,*,obs_i]=abs((*cal_bad.gain[0]))
  endfor
  
  stop
  
  
  quick_image, cal_divide_norm, data_range=[0,1], title='Percent change of 1067257048 calibration amplitude', xtitle='Freq channel index', ytitle='Tile index'
end