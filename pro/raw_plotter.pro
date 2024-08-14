pro raw_plotter

  day=['Aug23']
  day_num=N_elements(day)
  
  pointing_num=[-2,-1,0,1,2,3]
  pointing_name=['-2','-1','0','1','2','3']
  
  parsednames=STRARR(day_num,N_elements(pointing_num))
  for day_i=0,day_num-1 do begin
    parsednames[day_i,*]=[day[day_i]+'_minustwo',day[day_i]+'_minusone',day[day_i]+'_zenith',day[day_i]+'_plusone',day[day_i]+'_plustwo',day[day_i]+'_plusthree']
  endfor
  
  
  obsid_day_pointing=PTRARR((size(pointing_num))[1],day_num,/allocate)
  obsid_count=INTARR((size(pointing_num))[1])
  
  FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
    undefine, obsid
    undefine, beginning_obsid
    
    FOR day_i=0,day_num-1 DO BEGIN
    
      filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/' + parsednames[day_i,j] + '.txt'
      If parsednames[day_i,j] EQ 'Aug23_plusthree' then filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23plusthree.txt' ;TEMP to get Aug23 +3
      
      obs_array='empty'
      readcol, filename, obs_array, format='A', /silent
      
      If obs_array[0] NE 'empty' Then begin
        textfast,obs_array,/read,file_path=filename,string=1
        
        If keyword_set(obsid) then obsid=[obsid,ULONG(obs_array[*])] else obsid=ULONG(obs_array[*])
        If keyword_set(beginning_obsid) then beginning_obsid=[beginning_obsid,ULONG(obs_array[0])] else beginning_obsid=ULONG(obs_array[0])
        *obsid_day_pointing[j, day_i]=(obs_array[*])
        
        
      endif
      
    ENDFOR ; end day for
    
    obsid_count[j]=N_elements(obsid)
    
  endfor ;end pointing for
  
  gain_raw=complex(DBLARR(2,day_num,384,128,(size(pointing_num))[1],20))
  for j=0,(size(pointing_num))[1]-1 do begin
    print, 'Reading in pointing ' + pointing_name[j]
    for day_i=0,day_num-1 do begin
    
      If *obsid_day_pointing[j, day_i] NE !NULL then begin
      
        for obs_i=0, N_elements(*obsid_day_pointing[j, day_i])-1 do begin
        
          for obs_i=0, N_elements(*obsid_day_pointing[j, day_i])-1 do begin
          
           ; filename= '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/'+strtrim((*obsid_day_pointing[j, day_i])[obs_i],2)+'_cal.sav'
            filename= '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/calibration/'+strtrim((*obsid_day_pointing[j, day_i])[obs_i],2)+'_cal.sav'
            restore,filename
            
            
            ;Reinstate the raw gains from the longrun
            filename= '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/'+strtrim((*obsid_day_pointing[j, day_i])[obs_i],2)+'_obs.sav'
            restore,filename
            for pol_i=0,1 do *cal.gain[pol_i]=*cal.gain[pol_i]+*cal.gain_residual[pol_i]
            
            ;Take the bandpass out of the raw gains. Just added the longrun bp to see effects (smoother freq dependence is what I predict)
            cal_bandpass=vis_cal_bandpass(cal,obs,cal_remainder=cal_remainder,cable_bandpass_fit=1, saved_run_bp=1)
            for pol_i=0,1 do *cal.gain[pol_i]=*cal_remainder.gain[pol_i]
            
            gain_raw[0,day_i,*,*,j,obs_i]=*cal.gain[0]  ;pol x day x freq x tile x pointing x obs
            gain_raw[1,day_i,*,*,j,obs_i]=*cal.gain[1]
            
          endfor
          
        endfor
        
      endif
      
    endfor
    
  endfor
  
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/1061316296_obs.sav'
  tile_names=ULONG((*obs.baseline_info).tile_names)
  freq_use=where((*obs.baseline_info).freq_use,nf_use)
  
  day_i=0
  
  ;cgPS_Open,'/nfs/eor-00/h1/nbarry/Aug23_longrunpoly/gain_rave_tile'+strtrim(tile_names[tile_i],2)+'.png',/quiet,/nomatch
  ;cgplot, abs(gain_raw[0,0,freq_use,tile_i,0,0]), yrange=[.8,1.2],xrange=[0,336], title='Scaled longrun raw gains by '+avestep+ ', tile '+strtrim(tile_names[tile_i],2)+' xx -2pointing',$
  ;  xtitle='frequency index', ytitle='scaled gain', charsize=1,/NODATA
  
  pol_names=['xx','yy']
  
  for tile_i=0,127 do begin
  
    for pol_i=0,1 do begin
    
      for j=0,(size(pointing_num))[1]-1 do begin
      
        cgPS_Open,'/nfs/eor-00/h1/nbarry/Aug23_longrunpoly_afterpolfix/gainfit_tile'+strtrim(tile_names[tile_i],2)+'_'+pointing_name[j]+'_'+pol_names[pol_i]+'.png',/quiet,/nomatch
        
        cgplot, abs(gain_raw[pol_i,0,freq_use,tile_i,j,0]), yrange=[.8,2],xrange=[0,336], title='Raw Gains and polyfit, tile '+strtrim(tile_names[tile_i],2)+' '+pointing_name[j]+' '+pol_names[pol_i],$
          xtitle='frequency index', ytitle='gain amplitude', charsize=1,/NODATA
          
        for obs_i=0, N_elements(*obsid_day_pointing[j, day_i])-1 do begin
          cgoplot, abs(gain_raw[pol_i,0,freq_use,tile_i,j,obs_i])
        endfor
        
        restore, '/nfs/eor-00/h1/nbarry/Aug23_longrunpoly_afterpolfix/'+(*obsid_day_pointing[j, day_i])[0]+'_cal.sav'
        cgoplot, abs((*cal.gain[pol_i])[freq_use,tile_i]),color='green'
        
        ;restore, '/nfs/eor-00/h1/nbarry/Aug23_std_test_twopolyquad_polyonly/'+(*obsid_day_pointing[j, day_i])[0]+'_cal.sav'
        ;cgoplot, abs((*cal.gain[pol_i])[freq_use,tile_i]),color='blue'
        
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        ;stop
      endfor
      
    endfor
    
  endfor
  stop
  
end