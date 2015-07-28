pro poly_comparisons, longrun=longrun, lst_compare=lst_compare, raw=raw, by_pointing=by_pointing, quick_check=quick_check, line_plots=line_plots, spread=spread
  ;Program for plotting polyfit by pointings solutions over the longrun to check for stability. Longrun must be set at this time.
  ;lst_compare is not implemented at this time. Raw underplots the raw (minus bandpass by cable by pointing) for each obsid of the
  ;pointing or day and pointing depending on the next keyword.  by_pointing plots each all the pointings over all the days
  ;for a specific tile.  If that is not set, then plots are made for each pointing for each day.
  ;quick_check reads in only the -2 pointing to make plots quicker
  ;line_plots makes the plots described above
  ;spread is a new option to plot spreads of the zeroth order amp params by obs , overplotted with by pointing

  ;days of the long run to use
  day=['Aug23','Aug27','Sep02','Sep04','Sep06','Sep09','Sep11','Sep13','Sep17','Sep19','Sep30','Oct02','Oct04','Oct08','Oct10','Oct15','Oct23','Oct25','Oct29']
  
  ;autos only
  ;day=['Aug23','Aug27']
  
  day_num=N_elements(day)
  rgbcolors=[[69,190,207],[144,23,3],[240,87,249],[45,165,30],[42,25,77],[1,53,8],[167,138,28],[251,193,236],[52,86,193],[246,229,194],[148,33,145],[253,88,89],[239,233,110],$
    [14,128,63],[96,50,14],[176,247,108],[74,126,157],[47,14,36],[212,120,250],[202,29,168]]
    
  ;Directory to print plots -- CAN CHANGE
  outdir='/nfs/eor-00/h1/nbarry/longrun_std_test_twopolyquad/plots/by_tile_raw_underplot_legend/'
  
  parsednames=STRARR(day_num,6)
  for day_i=0,day_num-1 do begin
    If keyword_set(longrun) then parsednames[day_i,*]=[day[day_i]+'_minustwo',day[day_i]+'_minusone',day[day_i]+'_zenith',day[day_i]+'_plusone',day[day_i]+'_plustwo',day[day_i]+'_plusthree']
  endfor
  
  pointing_num=[-2,-1,0,1,2,3]
  pointing_name=['-2','-1','0','1','2','3']
  If keyword_set(quick_check) then pointing_num=[-2]
  If keyword_set(quick_check) then ponting_name=['-2']
  
  ;***************Loop to get the obsids of each day/pointing of the longrun
  obsid_day_pointing=PTRARR((size(pointing_num))[1],day_num,/allocate)
  obsid_count=INTARR((size(pointing_num))[1])
  
  FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
    undefine, obsid
    
    FOR day_i=0,day_num-1 DO BEGIN
    
      If keyword_set(longrun) then filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/' + parsednames[day_i,j] + '.txt'
      
      obs_temp='empty'
      readcol, filename, obs_temp, format='A'
      
      If obs_temp[0] NE 'empty' Then begin
        If keyword_set(obsid) then obsid=[obsid,obs_temp] else obsid=obs_temp
        *obsid_day_pointing[j, day_i]=obs_temp
      endif
      
    ENDFOR ; end day for
    
    obsid_count[j]=N_elements(obsid)
    
  endfor ;end pointing for
  ;***************End of loop to get the obsids of each day/pointing of the longrun
  
  
  ;****************************Read in the polyfit calibration solutions, and make the raw gains minus bp sol if specified
  gain=complex(DBLARR(2,day_num,384,128,(size(pointing_num))[1]))
  If keyword_set(raw) then gain_raw=complex(DBLARR(2,day_num,384,128,(size(pointing_num))[1],20))
  If keyword_set(raw) AND keyword_set(spread) then gain_obs_amp=complex(DBLARR(2,day_num,128,(size(pointing_num))[1],20))
  If keyword_set(spread) then gain_poi_amp=complex(DBLARR(2,day_num,128,(size(pointing_num))[1]))
  
  file_result=0
  If keyword_set(spread) then begin
    file_result1=file_test(outdir+'gain_poi_amp.sav')
    If keyword_set(raw) then file_result2=file_test(outdir+'gain_poi_amp.sav') else file_result2=1
    file_result=file_result1*file_result2
  endif
  If keyword_set(line_plots) then begin
    file_result1=file_test(outdir+'gain.sav')
    If keyword_set(raw) then file_result2=file_test(outdir+'gain_poi_amp.sav') else file_result2=1
    file_result=file_result1*file_result2
  endif
  
  If file_result EQ 0 then begin
  
    for j=0,(size(pointing_num))[1]-1 do begin
      print, 'Reading in pointing ' + pointing_name[j]
      for day_i=0,day_num-1 do begin
      
        If *obsid_day_pointing[j, day_i] NE !NULL then begin
          filename='/nfs/eor-00/h1/nbarry/longrun_std_test_twopolyquad/'+day[day_i]+'/'+strtrim(string((*obsid_day_pointing[j, day_i])[0]),2)+'_cal.sav'
          restore, filename
          
          gain[0,day_i,*,*,j]=*cal.gain[0]
          gain[1,day_i,*,*,j]=*cal.gain[1]
          
          If keyword_set(spread) then begin
            for tile_i=0, 127 do begin
              ;IF (*cal.amp_params[0,tile_i]) NE !NULL then gain_poi_amp[0,day_i,tile_i,j]=(*cal.amp_params[0,tile_i])[0] ;just the amplitude (zeroth order)
              ;IF (*cal.amp_params[1,tile_i]) NE !NULL then gain_poi_amp[1,day_i,tile_i,j]=(*cal.amp_params[1,tile_i])[0] ;just the amplitude (zeroth order)
              gain_poi_amp[0,day_i,tile_i,j]=mean(abs(gain[0,day_i,0:255,tile_i,j]))
              gain_poi_amp[1,day_i,tile_i,j]=mean(abs(gain[1,day_i,0:255,tile_i,j]))
            endfor
          endif
          
          If keyword_set(raw) then begin
            for obs_i=0, N_elements(*obsid_day_pointing[j, day_i])-1 do begin
            
              filename= '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/'+strtrim(string((*obsid_day_pointing[j, day_i])[obs_i]),2)+'_cal.sav'
              restore,filename
              If ~keyword_set(spread) AND keyword_set(line_plots) then begin
                filename= '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/'+strtrim(string((*obsid_day_pointing[j, day_i])[obs_i]),2)+'_obs.sav'
                restore,filename
                for pol_i=0,1 do *cal.gain[pol_i]=*cal.gain[pol_i]+*cal.gain_residual[pol_i]
                
                cal_bandpass=vis_cal_bandpass(cal,obs,cal_remainder=cal_remainder,saved_run_bp=1,cable_bandpass_fit=1)
                for pol_i=0,1 do *cal.gain[pol_i]=*cal_remainder.gain[pol_i]
                
                gain_raw[0,day_i,*,*,j,obs_i]=*cal.gain[0]
                gain_raw[1,day_i,*,*,j,obs_i]=*cal.gain[1]
              endif
              
              If keyword_set(spread) then begin
                gain_arr=complex(DBLARR(2,384,128))
                for tile_i=0,127 do begin
                  for pol_i=0,1 do begin
                    gain_fit=DBLARR(384)
                    phase_fit=DBLARR(384)
                    amp_params=*(cal.amp_params[pol_i,tile_i])
                    phase_params=*(cal.phase_params[pol_i,tile_i])
                    FOR di=0L,2 DO gain_fit+=amp_params[di]*findgen(384)^di
                    FOR di=0L,1 DO phase_fit+=phase_params[di]*findgen(384)^di
                    
                    gain_arr[pol_i,*,tile_i]=gain_fit*Exp(Complex(0,1)*phase_fit)
                  endfor
                  
                  ;IF (*cal.amp_params[0,tile_i]) NE !NULL then gain_obs_amp[0,day_i,tile_i,j,obs_i]=(*cal.amp_params[0,tile_i])[0] ;just the amplitude (zeroth order)
                  ;IF (*cal.amp_params[1,tile_i]) NE !NULL then gain_obs_amp[1,day_i,tile_i,j,obs_i]=(*cal.amp_params[1,tile_i])[0] ;just the amplitude (zeroth order)
                  gain_obs_amp[0,day_i,tile_i,j,obs_i]=mean(abs(gain_arr[0,0:255,tile_i]))
                  gain_obs_amp[1,day_i,tile_i,j,obs_i]=mean(abs(gain_arr[1,0:255,tile_i]))
                endfor
              endif
              
            endfor
          endif
          
        endif else begin
          gain[0,day_i,*,*,j]=1
          gain[1,day_i,*,*,j]=1
        endelse
        
      endfor
    endfor ;end pointing for
    save, gain_obs_amp, filename=outdir+'gain_obs_amp.sav'
    save, gain_poi_amp, filename=outdir+'gain_poi_amp.sav'
    save, gain_raw, filename=outdir+'gain_raw.sav'
  endif else begin
    print, 'Restoring'
    If keyword_set(spread) AND keyword_set(raw) then restore, oudir+'gain_obs_amp.sav'
    If keyword_set(spread) then restore, outdir+'gain_poi_amp.sav'
    If keyword_set(line_plots) AND keyword_set(raw) then restore, outdir+'gain_raw.sav'
    If keyword_set(line_plots) then restore, outdir+'gain.sav'
  endelse
  
  ;****************************End of read in the polyfit calibration solutions, and make the raw gains minus bp sol if specified
  
  for pol_i=0,1 do begin
    if pol_i EQ 0 then pol_name='xx'
    if pol_i EQ 1 then pol_name='yy'
    
    ;line_plots=1
    poi_name=['-2','-1','0','1','2','3']
    If keyword_set(line_plots) then begin
    
      for tile_i=0,127 do begin
      
        If keyword_set(by_pointing) then begin
          savelocation=outdir+strtrim(string(tile_i),2)+'_line_'+pol_name+'.png'
          cgPS_Open,savelocation,/quiet,/nomatch
          
          title='tile ' +strtrim(string(tile_i),2)+ ', longrun polyfit solutions, '+pol_name
          
          cgplot, abs(gain[pol_i,0,*,tile_i,0]),xtitle='Frequency Index',ytitle='Gain', $
            title=title,charsize=1, psym=2, symsize=0.2, yrange=[1,2], xrange=[0,384], /NODATA
            
          If keyword_set(raw) then begin
            for day_i=0,day_num-1 do begin
              for j=0,(size(pointing_num))[1]-1 do begin
                temp=where(abs(gain_raw[pol_i,day_i,0,tile_i,j,*]) NE 0, obs_num)
                If obs_num EQ -1 then obs_num=0
                for obs_i=0,obs_num-1 do cgoplot, abs(gain_raw[pol_i,day_i,*,tile_i,j,obs_i]), color='light grey',psym=2, symsize=0.2
              endfor
            endfor
          endif
          
          Device, Decomposed=0
          for day_i=0,day_num-1 do begin
            for j=0,(size(pointing_num))[1]-1 do begin
            
              TVLCT, rgbcolors[0,j], rgbcolors[1,j], rgbcolors[2,j],100
              cgoplot, abs(gain[pol_i,day_i,*,tile_i,j]), color=100B,psym=2, symsize=0.2
              
            endfor ;end pointing for
          endfor
          Device, Decomposed=1
          
          cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
          
        endif else begin
        
        
          for j=0,((size(pointing_num))[1])-1 do begin
          
            savelocation=outdir+strtrim(string(tile_i),2)+'_'+poi_name[j]+'_line_'+pol_name+'.png'
            cgPS_Open,savelocation,/quiet,/nomatch
            
            title='tile ' +strtrim(string(tile_i),2)+', pointing '+poi_name[j]+', longrun polyfit solutions, '+pol_name
            
            cgplot, abs(gain[pol_i,0,*,tile_i,0]),xtitle='Frequency Index',ytitle='Gain', $
              title=title,charsize=1, psym=2, symsize=0.2, yrange=[1,2], xrange=[0,384], /NODATA
              
            If keyword_set(raw) then begin
              for day_i=0,day_num-1 do begin
                temp=where(abs(gain_raw[pol_i,day_i,0,tile_i,j,*]) NE 0, obs_num)
                If obs_num EQ -1 then obs_num=0
                for obs_i=0,obs_num-1 do cgoplot, abs(gain_raw[pol_i,day_i,*,tile_i,j,obs_i]), color='light grey',psym=2, symsize=0.2
              endfor
            endif
            
            undefine,color_byte
            
            Device, Decomposed=0
            for day_i=0,day_num-1 do begin
              TVLCT, rgbcolors[0,day_i], rgbcolors[1,day_i], rgbcolors[2,day_i],100+day_i
              If keyword_set(color_byte) then color_byte=[color_byte,BYTE(100+day_i)] else color_byte=BYTE(100+day_i)
              cgoplot, abs(gain[pol_i,day_i,*,tile_i,j]), color=color_byte[day_i],psym=2, symsize=0.2
            endfor
            
            cgLegend, Title=day[0:5], $
              Color=color_byte[0:5], Location=[0.2,0.87],charsize=0.6,VSpace=.9,thick=2
            cgLegend, Title=day[6:11], $
              Color=color_byte[6:11], Location=[0.4,0.87],charsize=0.6,VSpace=.9,thick=2
            cgLegend, Title=day[12:18], $
              Color=color_byte[12:18], Location=[0.6,0.87],charsize=0.6,VSpace=.9,thick=2
              
              
            Device, Decomposed=1
            cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
          endfor
          
          
        endelse
        
        
        
      endfor
    endif ;end of line plots
    
    if keyword_set(spread) then begin
    
      for tile_i=0, 127 do begin
        savelocation=outdir+strtrim(string(tile_i),2)+'_spread_'+pol_name+'.png'
        cgPS_Open,savelocation,/quiet,/nomatch
        
        title='tile ' +strtrim(string(tile_i),2)+', longrun amp solutions, '+pol_name
        minmax_array=gain_obs_amp[0,0,tile_i,*,*] ;CAN CHANGE 0 TO day_i
        temp_indices=where(minmax_array NE 0)
        yrange=minmax(minmax_array[temp_indices])
        
        cgplot,[-2,-1,0,1,2,3],[2,2,2,2,2,2], title=title,xtitle='Pointing Index',ytitle='Amplitude (zeroth order of fit)',$
          xrange=[-3,4],yrange=yrange,/NODATA
          
        day_num=1
        
        for day_i=0, day_num-1 do begin
        
          Device, Decomposed=0
          for j=0,(size(pointing_num))[1]-1 do begin
            TVLCT, rgbcolors[0,day_i], rgbcolors[1,day_i], rgbcolors[2,day_i],100
            temp= where(gain_obs_amp[0,day_i,tile_i,j,*] NE 0,obs_count)
            If obs_count EQ 0 then continue
            for obs_i=0, obs_count-1 do cgoplot, j-2,gain_obs_amp[0,day_i,tile_i,j,obs_i],color=100B,psym=2, symsize=0.5
          endfor
          Device, Decomposed=1
          
        endfor
        
        ;to overplot the pointing value
        for day_i=0, day_num-1 do begin
          for j=0,(size(pointing_num))[1]-1 do begin
            cgoplot, j-2,gain_poi_amp[0,day_i,tile_i,j],color='black',psym=7, symsize=0.5
          endfor
        endfor
        
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      endfor
      
    endif
    
    
  endfor
  
  
end