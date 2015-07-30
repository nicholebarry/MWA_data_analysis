pro poly_comparisons_bftemp, longrun=longrun, lst_compare=lst_compare, raw=raw, mean_tile=mean_tile, quick_check=quick_check, line_plots=line_plots, spread=spread
  ;Program for plotting polyfit by pointings solutions over the longrun to check for stability. Longrun must be set at this time.
  ;lst_compare is not implemented at this time. Raw underplots the raw (minus bandpass by cable by pointing) for each obsid of the
  ;pointing or day and pointing depending on the next keyword.  by_pointing plots each all the pointings over all the days
  ;for a specific tile.  If that is not set, then plots are made for each pointing for each day.
  ;quick_check reads in only the -2 pointing to make plots quicker
  ;line_plots makes the plots described above
  ;spread is a new option to plot spreads of the zeroth order amp params by obs , overplotted with by pointing

  ;days of the long run to use
  day=['Aug23','Aug27','Sep02','Sep04','Sep06','Sep09','Sep11','Sep13','Sep17','Sep19','Sep30','Oct02','Oct04','Oct08','Oct10','Oct15','Oct23','Oct25','Oct29']
  
  day_num=N_elements(day)
  rgbcolors=[[69,190,207],[144,23,3],[240,87,249],[45,165,30],[42,25,77],[1,53,8],[167,138,28],[251,193,236],[52,86,193],[246,229,194],[148,33,145],[253,88,89],[239,233,110],$
    [14,128,63],[96,50,14],[176,247,108],[74,126,157],[47,14,36],[212,120,250],[202,29,168]]
    
  ;Directory to print plots -- CAN CHANGE
  outdir='/nfs/eor-00/h1/nbarry/longrun_poly_bftemp/fit_added/'
  
  pointing_num=[-2,-1,0,1,2];,3]
  ;pointing_num=[-5,-4,-3,-2,-1,0,1,2,3,4]
  pointing_name=['-2','-1','0','1','2'];,'3']
  ;pointing_name=['-5','4','-3','-2','-1','0','1','2','3','4']
  
  parsednames=STRARR(day_num,N_elements(pointing_num))
  for day_i=0,day_num-1 do begin
    If keyword_set(longrun) then parsednames[day_i,*]=[day[day_i]+'_minustwo',day[day_i]+'_minusone',day[day_i]+'_zenith',day[day_i]+'_plusone',day[day_i]+'_plustwo'];,day[day_i]+'_plusthree']
  ;parsednames[day_i,*]=[day[day_i]+'_minusfive',day[day_i]+'_minusfour',day[day_i]+'_minusthree',day[day_i]+'_minustwo',$
  ;  day[day_i]+'_minusone',day[day_i]+'_zenith',day[day_i]+'_plusone',day[day_i]+'_plustwo',day[day_i]+'_plusthree',day[day_i]+'_plusfour']
  endfor
  
  If keyword_set(quick_check) then pointing_num=[-2]
  If keyword_set(quick_check) then ponting_name=['-2']
  
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/1061316296_obs.sav'
  tile_names=ULONG((*obs.baseline_info).tile_names)
  
  ;***************Loop to get the obsids of each day/pointing of the longrun
  obsid_day_pointing=PTRARR((size(pointing_num))[1],day_num,/allocate)
  obsid_count=INTARR((size(pointing_num))[1])
  beginning_obsid_count=INTARR((size(pointing_num))[1])
  bftemp_pointing=FLTARR((size(pointing_num))[1],day_num)
  
  ;GET_LUN,lun
  ;OPENW,lun,outdir+'obs_to_query_v2.txt'
  
  FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
    undefine, obsid
    undefine, beginning_obsid
    
    FOR day_i=0,day_num-1 DO BEGIN
    
      If keyword_set(longrun) then filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/' + parsednames[day_i,j] + '.txt'
      
      obs_array='empty'
      readcol, filename, obs_array, format='A'
      
      If obs_array[0] NE 'empty' Then begin
        textfast,obs_array,/read,file_path=filename,string=1
        
        If keyword_set(obsid) then obsid=[obsid,ULONG(obs_array[*])] else obsid=ULONG(obs_array[*])
        If keyword_set(beginning_obsid) then beginning_obsid=[beginning_obsid,ULONG(obs_array[0])] else beginning_obsid=ULONG(obs_array[0])
        *obsid_day_pointing[j, day_i]=(obs_array[*])
      ;PRINTF, lun,(*obsid_day_pointing[j,day_i])[0]
      ;if day_i EQ 0 then stop
      ;If (size(data_array))[1] GT 1 then bftemp_pointing[j,day_i]=mean(FLOAT(data_array[39,*]))
      endif
      
    ENDFOR ; end day for
    
    obsid_count[j]=N_elements(obsid)
    beginning_obsid_count[j]=N_elements(beginning_obsid)
    
  endfor ;end pointing for
  ;***************End of loop to get the obsids of each day/pointing of the longrun
  
  ;FREE_LUN,lun
  
  ;stop
  
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
  
  If file_result EQ 0 then begin
  
    for j=0,(size(pointing_num))[1]-1 do begin
      print, 'Reading in pointing ' + pointing_name[j]
      for day_i=0,day_num-1 do begin
      
        If *obsid_day_pointing[j, day_i] NE !NULL then begin
          filename='/nfs/eor-00/h1/nbarry/longrun_std_test_twopolyquad/'+day[day_i]+'/'+(*obsid_day_pointing[j, day_i])[0]+'_cal.sav'
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
            
              filename= '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/'+(*obsid_day_pointing[j, day_i])[obs_i]+'_cal.sav'
              restore,filename
              If ~keyword_set(spread) AND keyword_set(line_plots) then begin
                filename= '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/'+(*obsid_day_pointing[j, day_i])[obs_i]+'_obs.sav'
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
  endif else begin
    print, 'Restoring'
    If keyword_set(spread) AND keyword_set(raw) then restore, outdir+'gain_obs_amp.sav'
    If keyword_set(spread) then restore, outdir+'gain_poi_amp.sav'
  endelse
  
  ;****************************End of read in the polyfit calibration solutions, and make the raw gains minus bp sol if specified
  
  
  
  filename='/nfs/eor-00/h1/nbarry/longrun_poly_bftemp/obs_queried_v2.txt'
  
  textfast,data_array,/read,file_path=filename,string=1
  
  temperature_array=FLTARR(day_num,128,8)
  ;Want day x tile x pointing
  
  tile_temp=STRARR(8)
  
  for j=0,(size(pointing_num))[1]-1 do begin
    for day_i=0,day_num-1 do begin
    
      If *obsid_day_pointing[j, day_i] NE !NULL then begin
      
        temp_index=where(strmatch(data_array,'*'+(*obsid_day_pointing[j, day_i])[0]+'*') EQ 1)
        split_temp_array=STRARR(N_elements(temp_index),3)
        for i=0,N_elements(temp_index)-1 do begin
          split_temp_array[i,*]=strsplit(data_array[temp_index[i]], '|',/EXTRACT)  ;Split line in data file into receiver, obsid full, and temps per tile
          split_temp_array[i,1]=(strsplit(split_temp_array[i,1], '|',/EXTRACT))[0]  ;Remove extra ms from end of obsid
          split_temp_array[i,0]=STRING(ULONG(split_temp_array[i,0])*10)  ;Add zero to end of receiver for getting tile names easier
          split_temp_array[i,2]=strmid(split_temp_array[i,2],1) ;Remove '{' from beginning of temp data
          str_length=strlen(split_temp_array[i,2]) ;Find length of string for next command
          split_temp_array[i,2]=strmid(split_temp_array[i,2],0,str_length-2) ;Remove '}' from end of temp data
          tile_temp[*]=strsplit(split_temp_array[i,2], ',',/EXTRACT) ;Split up temp data into 8 discrete temperatures per receiver
          
          
          
          for tile_i=0,7 do begin
            tile_index=where(tile_names EQ (ULONG(split_temp_array[i,0])+tile_i+1))
            temperature_array[day_i,tile_index,j]=FLOAT(tile_temp[tile_i]) ;Add per tile data into array as a float
          endfor
          
        endfor
        
      endif
      
    endfor
    
  endfor
  
  
  
  
  for pol_i=0,1 do begin
    if pol_i EQ 0 then pol_name='xx'
    if pol_i EQ 1 then pol_name='yy'
    
    
    poi_name=['-2','-1','0','1','2','3']
    If keyword_set(line_plots) then begin
    
    
      If keyword_set(mean_tile) then begin
        for j=0,((size(pointing_num))[1])-1 do begin
        
          savelocation=outdir+'beamformertemp_'+poi_name[j]+'_line_'+pol_name+'.png'
          cgPS_Open,savelocation,/quiet,/nomatch
          
          title='mean of tiles, pointing '+poi_name[j]+', longrun polyfit solutions, '+pol_name
          
          cgplot, mean(mean(abs(gain[pol_i,0,0:255,*,0]),dimension=4)),moon_dist_pointing[0,0],xtitle='Average Gain',ytitle='Jupiter Dist', $
            title=title,charsize=1, psym=2, symsize=0.2, xrange=[1.2,1.8],yrange=[110,120],/NODATA
            
          If keyword_set(raw) then begin
            for day_i=0,day_num-1 do begin
              temp=where(abs(gain_raw[pol_i,day_i,0,tile_i,j,*]) NE 0, obs_num)
              If obs_num EQ -1 then obs_num=0
              for obs_i=0,obs_num-1 do cgoplot, abs(gain_raw[pol_i,day_i,*,tile_i,j,obs_i]), color='light grey',psym=2, symsize=0.2
            endfor
          endif
          
          Device, Decomposed=0
          for day_i=0,day_num-1 do begin
            TVLCT, rgbcolors[0,day_i], rgbcolors[1,day_i], rgbcolors[2,day_i],100
            cgoplot, mean(mean(abs(gain[pol_i,day_i,0:255,*,j]),dimension=4)),moon_dist_pointing[j,day_i], color=100B,psym=2, symsize=0.6
          endfor
          Device, Decomposed=1
          cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        endfor
        
      endif else begin
      
        ;GET_LUN,lun
        ;OPENW,lun,outdir+'linear_fits.txt'
      
        for tile_i=2,127 do begin
        
          for j=2,((size(pointing_num))[1])-1 do begin
          
            savelocation=outdir+strtrim(string(tile_names[tile_i]),2)+'_'+poi_name[j]+'_beamformertemp_fit_'+pol_name+'.png'
            ;cgPS_Open,savelocation,/quiet,/nomatch
            
            title='Tile ' + strtrim(string(tile_names[tile_i]),2)+', '+poi_name[j]+' pointing, gain temperature correlation, '+pol_name
            
            cgplot, mean(abs(gain[pol_i,0,0:255,tile_i,0])),mean(abs(gain[pol_i,0,0:255,tile_i,0])),xtitle='Average Gain (after bandpass removal)',ytitle='Beamformer Temperature (C)', $
              title=title,charsize=1, psym=2, symsize=0.5,yrange=[10,40],xrange=[1.1,1.8],/NODATA
              
            input=mean(abs(gain[pol_i,*,0:255,tile_i,j]),dimension=3)
            input_index=where(input NE 1)
            result = LINFIT(reform(input[input_index]),reform(temperature_array[input_index,tile_i,j]),chisqr=chisqr)
            ;resistant_mean,sigma,2,res_mean
            
            ;Before bad day calc
            input=mean(abs(gain[pol_i,0:15,0:255,tile_i,j]),dimension=3)
            input_index=where(input NE 1)
            result2 = LINFIT(reform(input[input_index]),reform(temperature_array[input_index,tile_i,j]),chisqr=chisqr2)
            ;After bad day calc
            input=mean(abs(gain[pol_i,16:day_num-1,0:255,tile_i,j]),dimension=3)
            input_index=where(input NE 1)
            result3 = LINFIT(reform(input[input_index]),reform(temperature_array[input_index,tile_i,j]),chisqr=chisqr3)
            
            ;Complicated sorting mechanism
            input=mean(abs(gain[pol_i,*,0:255,tile_i,j]),dimension=3)
            input_index=where(input NE 1)
            result = LINFIT(reform(input[input_index]),reform(temperature_array[input_index,tile_i,j]),chisqr=chisqr)
            
            temperature_temp=reform(temperature_array[input_index,tile_i,j])
            sorted_index=sort(temperature_temp)
            temperature_sorted=temperature_temp[sorted_index]
            input_sorted=reform(input[input_index[sorted_index]])
            
            diff_array=abs(input_sorted[0]-input_sorted[1])
            for diff_i=1, N_elements(input_sorted)-2 do begin
              diff_array=[diff_array, abs(input_sorted[diff_i+1]-input_sorted[diff_i])]
            endfor
            
            grouping_ptr=PTRARR(N_elements(input_sorted), /allocate)
            breaks=where(diff_array GT 1.5*mean(diff_array),break_count)
            If break_count GT 0 then begin
            
              for break_i=0, N_elements(breaks)-1 do begin
                for previous_i=1, breaks[break_i] do begin
                  previous_diff= abs(input_sorted[break_i]-input_sorted[break_i-previous_i])
                  IF previous_diff LT 1.5*mean(diff_array) then begin
                    if keyword_set(grouping) then grouping=[grouping,breaks[break_i]-previous_i] else grouping=breaks[break_i]-previous_i
                  endif
                endfor
                for after_i=1, N_elements(input_sorted)-1-breaks[break_i] do begin
                  after_diff= abs(input_sorted[break_i]-input_sorted[break_i+after_i])
                  IF after_diff LT 1.5*mean(diff_array) then begin
                    if keyword_set(grouping) then grouping=[grouping,breaks[break_i]+after_i] else grouping=breaks[break_i]+after_i
                  endif
                endfor
                *grouping_ptr[breaks[break_i]]=grouping
                undefine, grouping
              endfor
              
            endif else *grouping_ptr[0]=sorted_index
            
            for group_i=0, N_elements(sorted_index)-1 do begin
              If *grouping_ptr[group_i] NE !NULL then begin
                input_index=Uniq(*grouping_ptr[group_i])
                results=LINFIT(reform(input_sorted[input_index]),reform(temperature_sorted[input_index]))
                cgoplot, [1.1,1.4,1.5,1.8],result[0]+result[1]*[1.1,1.4,1.5,1.8]
              endif
              
            endfor
            
            
            ;If (result3[0]-100 GT result2[0]) OR (result3[0]+100 LT result2[0]) then result=result2
            
            
            ;if j NE ((size(pointing_num))[1])-1 then temp_diff=temperature_array[input_index,tile_i,j]-temperature_array[input_index,tile_i,j+1] else $
            ;  temp_diff=temperature_array[input_index,tile_i,j]-temperature_array[input_index,tile_i,j-1]
            ;result_to_print=[0,poi_name[j],tile_i,result,res_mean]
            ;printf,lun, result_to_print
            ;cgoplot, [1.1,1.4,1.5,1.8],result[0]+result[1]*[1.1,1.4,1.5,1.8]
            
            If keyword_set(raw) then begin
              for day_i=0,day_num-1 do begin
                temp=where(abs(gain_raw[pol_i,day_i,0,tile_i,j,*]) NE 0, obs_num)
                If obs_num EQ -1 then obs_num=0
                for obs_i=0,obs_num-1 do cgoplot, abs(gain_raw[pol_i,day_i,*,tile_i,j,obs_i]), color='light grey',psym=2, symsize=0.6
              endfor
            endif
            
            undefine, color_array
            cgLoadCT, 25, clip=[0,190], Ncolors=day_num+1
            Device, Decomposed=0
            for day_i=0,day_num-1 do begin
              ;TVLCT, rgbcolors[0,day_i], rgbcolors[1,day_i], rgbcolors[2,day_i],10+day_i
              cgoplot, mean(abs(gain[pol_i,day_i,0:255,tile_i,j])),temperature_array[day_i,tile_i,j], color=day_i+1,psym=2, symsize=1
              if ~keyword_set(color_array) then color_array=day_i+1 else color_array=[color_array,day_i+1]
              
            endfor
            
            stop
            cgLegend, Title=day, $
              Color=color_array, Location=[0.8,0.87],charsize=0.7,VSpace=.9,thick=2, $
              /center_sym,psyms=2, length=0
              
            Device, Decomposed=1
            ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
          endfor
          
        endfor
      ;FREE_LUN,lun
      endelse
      
      
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