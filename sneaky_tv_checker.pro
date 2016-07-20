pro sneaky_TV_checker
  ;check for sneaky TV contamination in tiles 16-23 that has a ratio of beginning band/end band in cal amplitudes of greater than 1.3.

  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  obs = GETVAR_SAVEFILE(dir+'metadata/1061316296_obs.sav', 'obs')
  freq_arr=(*obs.baseline_info).freq/10.^6.
  
  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/beardsley_thesis_list.txt'
  readcol, filename, obs_id, format='A', /silent
  
  for obs_i=0,N_elements(obs_id)-1 do begin
    cal = getvar_savefile('/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/'+obs_id[obs_i]+'_cal.sav','cal')
    
    cal_divide_xx=mean(abs((*cal.gain[0])[0:100,*]),dimension=1)/mean(abs((*cal.gain[0])[283:383,*]),dimension=1)
    cal_divide_yy=mean(abs((*cal.gain[1])[0:100,*]),dimension=1)/mean(abs((*cal.gain[1])[283:383,*]),dimension=1)
    ;for tile_i=0,127 do if ~finite(cal_divide_xx[tile_i]) then cal_divide_xx[tile_i]=0
    ;for tile_i=0,127 do if ~finite(cal_divide_yy[tile_i]) then cal_divide_yy[tile_i]=0
    
    flag_xx=0
    flag_yy=0
    
    for tile_i=0,100 do begin
      if ~finite(cal_divide_xx[tile_i]) then cal_divide_xx[tile_i]=0
      if ~finite(cal_divide_yy[tile_i]) then cal_divide_yy[tile_i]=0
      
      ;if (cal_divide_xx[16] GT 1.3) AND (cal_divide_xx[17] GT 1.3) AND (cal_divide_xx[18] GT 1.3) AND (cal_divide_xx[19] GT 1.3) AND (cal_divide_xx[20] GT 1.3) $
      ;  AND (cal_divide_xx[21] GT 1.3) AND (cal_divide_xx[22] GT 1.3) AND (cal_divide_xx[23] GT 1.3) then begin
      
      if (cal_divide_xx[tile_i] GT 1.5) AND (flag_xx EQ 0) then begin
      
        quick_image, abs(*cal.gain[0]), title=obs_id[obs_i]+' calibration amplitude, xx, flagged TV', xtitle='Freq channel index', ytitle='Tile index',window=obs_i, $
          savefile='/nfs/mwa-00/h1/nbarry/TVcheck/all_tiles/xx_'+string(obs_i,format='(I4.4)')
        print, 'xx found: '  +string(obs_i,format='(I4.4)') + ' '  + string(tile_i,format='(I3.3)')
        flag_xx=1
        
      endif
      
      if (cal_divide_yy[tile_i] GT 1.5) AND (flag_yy EQ 0) then begin
      
        quick_image, abs(*cal.gain[1]), title=obs_id[obs_i]+' calibration amplitude, yy, flagged TV', xtitle='Freq channel index', ytitle='Tile index',window=obs_i, $
          savefile='/nfs/mwa-00/h1/nbarry/TVcheck/all_tiles/yy_'+string(obs_i,format='(I4.4)')
        print, 'yy found: '  +string(obs_i,format='(I4.4)') + ' ' + string(tile_i,format='(I3.3)')
        flag_yy=1
        
      endif
      
    endfor
    
  ;if (cal_divide_yy[16] GT 1.3) AND (cal_divide_yy[17] GT 1.3) AND (cal_divide_yy[18] GT 1.3) AND (cal_divide_yy[19] GT 1.3) AND (cal_divide_yy[20] GT 1.3) $
  ;  AND (cal_divide_yy[21] GT 1.3) AND (cal_divide_yy[22] GT 1.3) AND (cal_divide_yy[23] GT 1.3) then begin
    
  ;quick_image, abs(*cal.gain[1]), title=obs_id[obs_i]+' calibration amplitude, yy, flagged TV', xtitle='Freq channel index', ytitle='Tile index',window=obs_i, $
  ;  savefile='/nfs/mwa-00/h1/nbarry/TVcheck/yy_'+string(obs_i,format='(I4.4)')
  ;print, 'yy found: '  +string(obs_i,format='(I4.4)')
    
  ;endif
    
    
  endfor
  
end