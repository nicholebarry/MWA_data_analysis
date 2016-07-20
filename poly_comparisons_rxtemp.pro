pro poly_comparisons_rxtemp, longrun=longrun, lst_compare=lst_compare, raw=raw, mean_tile=mean_tile, quick_check=quick_check, line_plots=line_plots, spread=spread, $
    make_obs_to_query=make_obs_to_query
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
  outdir='/nfs/eor-00/h1/nbarry/longrun_poly_rxtemp/'
  
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
  
  if keyword_set(make_obs_to_query) then begin
    GET_LUN,lun
    OPENW,lun,outdir+'obs_to_query_v2.txt'
  endif
  
  FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
    undefine, obsid
    undefine, beginning_obsid
    
    FOR day_i=0,day_num-1 DO BEGIN
    
      If keyword_set(longrun) then filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/' + parsednames[day_i,j] + '.txt'
      
      obs_array='empty'
      readcol, filename, obs_array, format='A', /silent
      
      If obs_array[0] NE 'empty' Then begin
        textfast,obs_array,/read,file_path=filename,string=1
        
        If keyword_set(obsid) then obsid=[obsid,ULONG(obs_array[*])] else obsid=ULONG(obs_array[*])
        If keyword_set(beginning_obsid) then beginning_obsid=[beginning_obsid,ULONG(obs_array[0])] else beginning_obsid=ULONG(obs_array[0])
        *obsid_day_pointing[j, day_i]=(obs_array[*])
        if keyword_set(make_obs_to_query) then PRINTF, lun,(*obsid_day_pointing[j,day_i])[0]
        
      endif
      
    ENDFOR ; end day for
    
    obsid_count[j]=N_elements(obsid)
    beginning_obsid_count[j]=N_elements(beginning_obsid)
    
  endfor ;end pointing for
  ;***************End of loop to get the obsids of each day/pointing of the longrun
  
  if keyword_set(make_obs_to_query) then begin
    FREE_LUN,lun
    stop
  endif
  
  ;****************************Read in the polyfit calibration solutions, and make the raw gains minus bp sol if specified
  gain=complex(DBLARR(2,day_num,384,128,(size(pointing_num))[1]))
  If keyword_set(raw) then gain_raw=complex(DBLARR(2,day_num,384,128,(size(pointing_num))[1],20))
  
  for j=0,(size(pointing_num))[1]-1 do begin
    print, 'Reading in pointing ' + pointing_name[j]
    for day_i=0,day_num-1 do begin
    
      If *obsid_day_pointing[j, day_i] NE !NULL then begin
        filename='/nfs/eor-00/h1/nbarry/longrun_std_test_twopolyquad/'+day[day_i]+'/'+(*obsid_day_pointing[j, day_i])[0]+'_cal.sav'
        restore, filename
        
        gain[0,day_i,*,*,j]=*cal.gain[0]
        gain[1,day_i,*,*,j]=*cal.gain[1]
        
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
            
          endfor
        endif
        
      endif else begin
        gain[0,day_i,*,*,j]=1
        gain[1,day_i,*,*,j]=1
      endelse
      
    endfor
  endfor ;end pointing for
  ;save, gain_obs_amp, filename=outdir+'gain_obs_amp.sav'
  ;save, gain_poi_amp, filename=outdir+'gain_poi_amp.sav'
  
  
  ;****************************End of read in the polyfit calibration solutions, and make the raw gains minus bp sol if specified
  
  
  ;****************************Read in and parse temperature data
  
  If ~file_test('/nfs/eor-00/h1/nbarry/rxtemp_chips_ave.sav') then begin
    ;Find the smaller temperature data file with just the right obs to get temperature data.
    ;To get this file, run this code with make_obs_to_query set and then run query.sh
    filename='/nfs/eor-00/h1/nbarry/rx_temps2.txt'
    
    ;Read out temperature data in the form of the string due to the funky format
    textfast,data_array,/read,file_path=filename,string=1, delimiter='|'
    
    temperature_array=FLTARR(day_num,128,4) ;Want day x tile x pointing
    chip_temp_pointing=FLTARR(16,day_num,(size(pointing_num))[1])
    
    tile_temp=STRARR(8)
    
    for j=0,(size(pointing_num))[1]-1 do begin
      for day_i=0,day_num-1 do begin
      
        If *obsid_day_pointing[j, day_i] NE !NULL then begin
        
          chip_temps_per_adfb=FLTARR(16, N_elements(*obsid_day_pointing[j, day_i]),16) ;rec x obs per pointing x slots
          
          for obs_i=0, N_elements(*obsid_day_pointing[j, day_i])-1 do begin
          
            temp_index=where(strmatch(data_array[1,*],(*obsid_day_pointing[j, day_i])[obs_i]) EQ 1,n_count)
            While (n_count EQ 0) do begin
              (*obsid_day_pointing[j, day_i])[obs_i]=strtrim(ULONG((*obsid_day_pointing[j, day_i])[obs_i])+1,2)
              temp_index=where(strmatch(data_array[1,*],(*obsid_day_pointing[j, day_i])[obs_i]) EQ 1,n_count)
            endwhile
            split_adfb1_array=STRARR(N_elements(temp_index),8)
            split_adfb2_array=STRARR(N_elements(temp_index),8)
            
            for i=0,N_elements(temp_index)-1 do begin
              split_adfb1_array[i,*]=strsplit(data_array[2,temp_index[i]], ',',/EXTRACT)  ;Split line in data file into receiver, obsid full, and temps per tile
              split_adfb1_array[i,0]=strmid(split_adfb1_array[i,0],1) ;Remove '{' from beginning of temp data
              str_length=strlen(split_adfb1_array[i,7]) ;Find length of string for next command
              split_adfb1_array[i,7]=strmid(split_adfb1_array[i,7],0,str_length-1) ;Remove '}' from end of temp data
              
              split_adfb2_array[i,*]=strsplit(data_array[3,temp_index[i]], ',',/EXTRACT)  ;Split line in data file into receiver, obsid full, and temps per tile
              split_adfb2_array[i,0]=strmid(split_adfb2_array[i,0],1) ;Remove '{' from beginning of temp data
              str_length=strlen(split_adfb2_array[i,7]) ;Find length of string for next command
              split_adfb2_array[i,7]=strmid(split_adfb2_array[i,7],0,str_length-1) ;Remove '}' from end of temp data
              
              chip_temps_per_adfb[ULONG(data_array[0,temp_index[i]])-1,obs_i,*]=[split_adfb1_array[i,*],split_adfb2_array[i,*]]
            ;chip_temp=mean(FLOAT(chip_temps_per_adfb))
              
              
            ;for tile_i=0,7 do begin
            ;  tile_index=where(tile_names EQ (ULONG(split_temp_array[i,0])+tile_i+1))
            ;  tile_temp_float=Convert_To_Type(tile_temp[tile_i],4)
            ;  temperature_array[day_i,tile_index,j]=tile_temp_float ;Add per tile data into array as a float
            ;endfor
            endfor
            
          endfor
          
          chip_temp_obs=mean(FLOAT(chip_temps_per_adfb),dimension=3)
          chip_temp_pointing[*,day_i,j]=mean(chip_temp_obs,dimension=2)
          
        endif
        
      endfor
      
    endfor
  endif else begin
    restore,'/nfs/eor-00/h1/nbarry/rxtemp_chips_ave.sav'
    print, 'Restored chip temps'
  endelse
  
  
  
  ;****************************End of read in and parse temperature data
  
    ;****************************Read in and parse temperature data
  
  ;Find the smaller temperature data file with just the right obs to get temperature data.
  ;To get this file, run this code with make_obs_to_query set and then run query.sh
  filename='/nfs/eor-00/h1/nbarry/longrun_poly_bftemp/obs_queried_v2.txt'
  
  ;Read out temperature data in the form of the string due to the funky format
  textfast,data_array,/read,file_path=filename,string=1
  
  temperature_array=FLTARR(day_num,128,8) ;Want day x tile x pointing
  
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
          split_temp_array[i,2]=strmid(split_temp_array[i,2],0,str_length-1) ;Remove '}' from end of temp data
          tile_temp[*]=strsplit(split_temp_array[i,2], ',',/EXTRACT) ;Split up temp data into 8 discrete temperatures per receiver
          
          
          
          for tile_i=0,7 do begin
            tile_index=where(tile_names EQ (ULONG(split_temp_array[i,0])+tile_i+1))
            tile_temp_float=Convert_To_Type(tile_temp[tile_i],4)
            temperature_array[day_i,tile_index,j]=tile_temp_float ;Add per tile data into array as a float
          endfor
        endfor
        
      endif
      
    endfor
    
  endfor
  
  ;****************************End of read in and parse temperature data
  
  ;GET_LUN,lun
  ;OPENW,lun,outdir+'linear_fits.txt',width=600
  ;header=['pointing','tile name','pol','Group 1 days','Group 1 y-inter','Group 1 slope','Group 1 chi^2','Group 2 days','Group 2 y-inter','Group 2 slope','Group 2 chi^2']
  ;printf,lun, header
  
  for pol_i=0,1 do begin
    if pol_i EQ 0 then pol_name='xx'
    if pol_i EQ 1 then pol_name='yy'
    
    
    poi_name=['-2','-1','0','1','2','3']
    If keyword_set(line_plots) then begin
    
      for tile_i=0,127 do begin
        ;print, 'Tile '+strtrim(tile_i,2)
        input_str=STRARR(11,((size(pointing_num))[1]))
        
        for j=0,((size(pointing_num))[1])-1 do begin
          undefine, group1results,group2results,group1,group2, breaks
          
          savelocation=outdir+strtrim(string(tile_names[tile_i]),2)+'_'+poi_name[j]+'_bftemp_rxtemp_'+pol_name+'.png'
          cgPS_Open,savelocation,/quiet,/nomatch
          
          title='Tile ' + strtrim(string(tile_names[tile_i]),2)+', '+poi_name[j]+' pointing, gain temperature correlation, '+pol_name
          
          cgplot, mean(abs(gain[pol_i,0,0:255,tile_i,0])),mean(abs(gain[pol_i,0,0:255,tile_i,0])),xtitle='Average Gain (after bandpass removal)',ytitle='Beamformer Temperature (C)', $
            title=title,charsize=1, psym=2, symsize=0.5,yrange=[0,50],xrange=[1,2],/NODATA
            
            
          If keyword_set(linear_fit) then begin
            ;*************Complicated sorting mechanism
            ;Get inputs which are not flagged
            input=mean(abs(gain[pol_i,*,0:255,tile_i,j]),dimension=3)
            input_index=where(input NE 1)
            ;result = LINFIT(reform(input[input_index]),reform(temperature_array[input_index,tile_i,j]),chisqr=chisqr)
            
            ;Sort unflagged inputs in order of coldest to hotest in anticipation that they will be grouped on temperature,
            ;though it is no requirement
            temperature_temp=reform(temperature_array[input_index,tile_i,j])
            sorted_index=sort(temperature_temp)
            temperature_sorted=temperature_temp[sorted_index]
            input_sorted=reform(input[input_index[sorted_index]])
            
            ;Find temperature differences between temp and next hotest temp for all pairs
            If N_elements(input_sorted) GT 1 then begin
              diff_array=abs(input_sorted[0]-input_sorted[1])
              temperature_diff_array=abs(temperature_sorted[0]-temperature_sorted[1])
              beg_slope_array=temperature_diff_array/diff_array
              for diff_i=1, N_elements(input_sorted)-2 do begin
                diff_array=[diff_array, abs(input_sorted[diff_i+1]-input_sorted[diff_i])]
                temperature_diff_array=[temperature_diff_array, abs(temperature_sorted[diff_i+1]-temperature_sorted[diff_i])]
                beg_slope_array=[beg_slope_array,(temperature_sorted[0]-temperature_sorted[diff_i+1])/(input_sorted[0]-input_sorted[diff_i+1])]
              endfor
              
              slope_diff_array=temperature_diff_array/diff_array
              pre_results=LINFIT(reform(input_sorted),reform(temperature_sorted))
              If pre_results[1] GT -50 then pre_results[1]=0
              distance_from_pre=abs((reform(temperature_sorted)-pre_results[0])/pre_results[1]-input_sorted)
              
              
              hypotenuse_array=(sqrt((diff_array*20.)^2.+temperature_diff_array^2.))
              
              ;Initialize grouping array to biggest possible group, with the idea that some entries will be NULL
              ;grouping_ptr=PTRARR(N_elements(input_sorted), /allocate)
              
              ;Find the "breaks", or deviations from the slope, using a arbitrary indicator (1.5* the mean diff)
              diff_logic=INTARR(N_elements(diff_array))
              diff_logic_break=where(diff_array GT .05,break_count)
              IF break_count GT 0 then diff_logic[diff_logic_break]=1
              
              dist_logic=INTARR(N_elements(diff_array))
              dist_logic_break=where(distance_from_pre GT .75,break_count)
              IF break_count GT 0 then dist_logic[dist_logic_break]=1
              if pre_results[1] EQ 0 then dist_logic[dist_logic_break]=0
              
              ;slope_logic=INTARR(N_elements(diff_array))
              ;slope_logic_break=where(slope_diff_array LT .25*mean(slope_diff_array),break_count)
              ;IF break_count GT 0 then slope_logic[slope_logic_break]=1
              
              beg_slope_logic=INTARR(N_elements(diff_array))
              flag=0
              ave_flag=0
              undefine, ave_slope
              for diff_i=0,N_elements(diff_array)-4 do begin
                slope_diff=(beg_slope_array[diff_i+1]-beg_slope_array[diff_i])/beg_slope_array[diff_i]
                slope_diff=[slope_diff,(beg_slope_array[diff_i+2]-beg_slope_array[diff_i+1])/beg_slope_array[diff_i+1]]
                slope_diff=[slope_diff,(beg_slope_array[diff_i+3]-beg_slope_array[diff_i+2])/beg_slope_array[diff_i+2]]
                ;print, slope_diff
                temp=where(abs(slope_diff) LT .15,slope_count)
                If (slope_count EQ 3) and (ave_flag EQ 0) then begin
                  ave_slope=mean(beg_slope_array[diff_i:diff_i+3])
                  ave_flag=1
                endif
              endfor
              If keyword_set(ave_slope) then begin
                for diff_i=4, N_elements(diff_array)-1 do begin
                
                  slope_diff=ave_slope-beg_slope_array[diff_i]
                  
                  If abs(slope_diff) LT abs(.25*ave_slope) then begin
                    ;ave_slope=((ave_slope+beg_slope_array[diff_i])/2+ave_slope)/2
                  
                    IF flag EQ 1 then begin
                      beg_slope_logic[diff_i]=1
                      flag=0
                    endif
                    
                  endif else begin
                  
                    IF (beg_slope_logic[diff_i-1] EQ 0) and (flag EQ 0) then beg_slope_logic[diff_i]=1
                    flag=1
                    
                  endelse
                endfor
              endif
              
              hypotenuse_logic=INTARR(N_elements(hypotenuse_array))
              hypotenuse_logic_break=where(hypotenuse_array LT 1,break_count)
              IF break_count GT 0 then hypotenuse_logic[hypotenuse_logic_break]=-1
              
              break_logic=diff_logic+dist_logic+hypotenuse_logic+beg_slope_logic
              breaks=where(break_logic GT 0, break_count)
              
            endif else break_count=0
            
            ;Begin looping of groups if there is a break as determined above
            If break_count GT 0 then begin
            
              pre_groups=PTRARR(N_elements(breaks)+1,/allocate)
              IF breaks[0] NE 0 then *pre_groups[0]=INDGEN(breaks[0]+1)
              for group_i=0,N_elements(breaks)-2 do begin
                array_size=breaks[group_i+1]-breaks[group_i]
                IF array_size GT 1 then *pre_groups[group_i+1]=INDGEN(array_size)+breaks[group_i]+1 else *pre_groups[group_i+1]=breaks[group_i]+1
              endfor
              last_index=N_elements(breaks)
              IF breaks[last_index-1] NE max(input_index)-1 then *pre_groups[last_index]=INDGEN(max(sorted_index)-breaks[last_index-1])+breaks[last_index-1]+1
              
              for group_i=0, N_elements(breaks)-1, 2 do begin
                if group1 NE !Null then group1=[group1,*pre_groups[group_i]] else group1=*pre_groups[group_i]
                if group2 NE !Null then group2=[group2,*pre_groups[group_i+1]] else group2=*pre_groups[group_i+1]
              endfor
              
            endif else group1=sorted_index
            
            
            If N_elements(group1) GT 1 then begin
              group1results=LINFIT(reform(input_sorted[group1]),reform(temperature_sorted[group1]),chisqr=group1chisqr)
              print, 'Group 1:'
              print, group1
            ;cgoplot, [1.0,1.4,1.5,1.8],group1results[0]+group1results[1]*[1.0,1.4,1.5,1.8]
            endif
            If N_elements(group2) GT 1 then begin
              group2results=LINFIT(reform(input_sorted[group2]),reform(temperature_sorted[group2]),chisqr=group2chisqr)
              print, 'Group 2:'
              print, group2
            ;cgoplot, [1.0,1.4,1.5,1.8],group2results[0]+group2results[1]*[1.0,1.4,1.5,1.8]
            endif
            ;*************End of Complicated sorting mechanism
            If (N_elements(group1) GT 1) AND  (N_elements(group2) GT 1) then begin
              crossing_point=(group2results[0]-group1results[0])/(group1results[1]-group2results[1])
              If (crossing_point GT min(input_sorted)-.1) AND (crossing_point LT max(input_sorted)+.1) then begin
                undefine, group2results,group2
                group1=sorted_index
                group1results=LINFIT(reform(input_sorted[group1]),reform(temperature_sorted[group1]),chisqr=group1chisqr)
                print, 'Group 2 ignored'
              endif
            endif
            
            ;IF group1results NE !Null then cgoplot, [1.0,1.4,1.5,1.8],group1results[0]+group1results[1]*[1.0,1.4,1.5,1.8]
            ;If group2results NE !Null then cgoplot, [1.0,1.4,1.5,1.8],group2results[0]+group2results[1]*[1.0,1.4,1.5,1.8]
            
            
            ;If (result3[0]-100 GT result2[0]) OR (result3[0]+100 LT result2[0]) then result=result2
            
            
            ;if j NE ((size(pointing_num))[1])-1 then temp_diff=temperature_array[input_index,tile_i,j]-temperature_array[input_index,tile_i,j+1] else $
            ;  temp_diff=temperature_array[input_index,tile_i,j]-temperature_array[input_index,tile_i,j-1]
            IF group1results EQ !NULL then begin
              group1results=[0,0]
              group1chisqr=0
            endif
            If group1 EQ !NULL then group1days=0 else begin
              group1days=' '+day[group1[0]]
              for i=1,N_elements(group1)-1 do group1days=group1days+','+day[group1[i]]
            endelse
            
            IF group2results EQ !NULL then begin
              group2results=[0,0]
              group2chisqr=0
            endif
            
            If group2 EQ !NULL then group2days=0 else begin
              group2days=' '+day[group2[0]]
              for i=1,N_elements(group2)-1 do group2days=group2days+','+day[group2[i]]
            endelse
            
            input_str[0,j]=poi_name[j]
            ;input_str[j,1]=Convert_to_type(tile_names[tile_i],7)
            tile_names_str=Convert_to_type(tile_names[tile_i],7)
            group1results_str=Convert_to_type(group1results,7)
            ;group1results_str=Convert_to_type(group1results,7)
            group1chisqr_str=Convert_to_type(group1chisqr,7)
            ;group1chisqr_str=Convert_to_type(group1chisqr,7)
            If group2days EQ 0 then group2days=Convert_to_type(group2days,7)
            group2results_str=Convert_to_type(group2results,7)
            group2chisqr_str=Convert_to_type(group2chisqr,7)
            input_str[*,j]=[poi_name[j],tile_names_str,pol_name,group1days,group1results_str,group1chisqr_str,group2days,group2results_str,group2chisqr_str]
            
          endif ;end linear fitting
          
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
          
          If keyword_set(linear_fit) then begin
            stop
            ;*****************
            If group2 EQ !NULL then group2days=0 else begin
              group2days=' '+day[group2[0]]
              for i=1,N_elements(group2)-1 do group2days=group2days+','+day[group2[i]]
            endelse
            undefine, group1results, group2results
            group1results=LINFIT(reform(input_sorted[group1]),reform(temperature_sorted[group1]),chisqr=group1chisqr)
            If group2 NE !NULL then group2results=LINFIT(reform(input_sorted[group2]),reform(temperature_sorted[group2]),chisqr=group2chisqr)
            IF group1results NE !Null then cgoplot, [1.0,1.4,1.5,1.8],group1results[0]+group1results[1]*[1.0,1.4,1.5,1.8]
            If group2results NE !Null then cgoplot, [1.0,1.4,1.5,1.8],group2results[0]+group2results[1]*[1.0,1.4,1.5,1.8]
            
            If group1 EQ !NULL then group1days=0 else begin
              group1days=' '+day[group1[0]]
              for i=1,N_elements(group1)-1 do group1days=group1days+','+day[group1[i]]
            endelse
            If group2 EQ !NULL then group2days=0 else begin
              group2days=' '+day[group2[0]]
              for i=1,N_elements(group2)-1 do group2days=group2days+','+day[group2[i]]
            endelse
            
            If group1 NE !NULL then begin
              print, 'Group1:   '
              print, group1days
              print, 'y:   ' +strtrim(group1results[0],2)
              print, 'm:   ' +strtrim(group1results[1],2)
              print, 'chi:    ' + strtrim(group1chisqr,2)
              print, group1days+'       '+strtrim(group1results[0],2)+'       '+strtrim(group1results[1],2)+'       '+strtrim(group1chisqr,2)
            endif
            If group2 NE !NULL then begin
              print, 'Group2:   '
              print, group2days
              print, 'y:   ' +strtrim(group2results[0],2)
              print, 'm:   ' +strtrim(group2results[1],2)
              print, 'chi:    ' + strtrim(group2chisqr,2)
              print, group1days+'       '+strtrim(group1results[0],2)+'       '+strtrim(group1results[1],2)+'       '+strtrim(group1chisqr,2)+'       '+$
                group2days+'       '+strtrim(group2results[0],2)+'       '+strtrim(group2results[1],2)+'       '+strtrim(group2chisqr,2)
            endif
            stop
            
            cgPS_Open,savelocation,/quiet,/nomatch
            cgplot, mean(abs(gain[pol_i,0,0:255,tile_i,0])),mean(abs(gain[pol_i,0,0:255,tile_i,0])),xtitle='Average Gain (after bandpass removal)',ytitle='Beamformer Temperature (C)', $
              title=title,charsize=1, psym=2, symsize=0.5,yrange=[0,50],xrange=[1,2],/NODATA
            undefine, color_array
            cgLoadCT, 25, clip=[0,190], Ncolors=day_num+1
            Device, Decomposed=0
            for day_i=0,day_num-1 do begin
              ;TVLCT, rgbcolors[0,day_i], rgbcolors[1,day_i], rgbcolors[2,day_i],10+day_i
              cgoplot, mean(abs(gain[pol_i,day_i,0:255,tile_i,j])),temperature_array[day_i,tile_i,j], color=day_i+1,psym=2, symsize=1
              if ~keyword_set(color_array) then color_array=day_i+1 else color_array=[color_array,day_i+1]
            endfor
            IF group1results NE !Null then cgoplot, [1.0,1.4,1.5,1.8],group1results[0]+group1results[1]*[1.0,1.4,1.5,1.8]
            If group2results NE !Null then cgoplot, [1.0,1.4,1.5,1.8],group2results[0]+group2results[1]*[1.0,1.4,1.5,1.8]
            
            ;*****************
            
            daytitle=day
            for star_i=0, N_elements(day)-1 do begin
              IF strmatch(group2days,'*'+day[star_i]+'*') then daytitle[star_i]=day[star_i]+' (2)'
            endfor
          endif else begin
            ;daytitle=day
            daytitle=STRARR(19)
            for day_i=0, day_num-1 do begin
              If tile_i LE 7 then reciever=0
              If tile_i GE 8 AND tile_i LE 15 then reciever=1
              If tile_i GE 16 AND tile_i LE 23 then reciever=2
              If tile_i GE 24 AND tile_i LE 31 then reciever=3
              If tile_i GE 32 AND tile_i LE 39 then reciever=4
              If tile_i GE 40 AND tile_i LE 47 then reciever=5
              If tile_i GE 48 AND tile_i LE 55 then reciever=6
              If tile_i GE 56 AND tile_i LE 63 then reciever=7
              If tile_i GE 64 AND tile_i LE 71 then reciever=8
              If tile_i GE 72 AND tile_i LE 79 then reciever=9
              If tile_i GE 80 AND tile_i LE 87 then reciever=10
              If tile_i GE 88 AND tile_i LE 95 then reciever=11
              If tile_i GE 96 AND tile_i LE 103 then reciever=12
              If tile_i GE 104 AND tile_i LE 111 then reciever=13
              If tile_i GE 112 AND tile_i LE 119 then reciever=14
              If tile_i GE 120 AND tile_i LE 127 then reciever=15
              daytitle[day_i]=day[day_i]+' '+strtrim(chip_temp_pointing[reciever,day_i,j],2)
            endfor
          endelse
          
          cgLegend, Title=daytitle, $
            Color=color_array, Location=[0.75,0.87],charsize=0.7,VSpace=.9,thick=2, $
            /center_sym,psyms=2, length=0
            
          Device, Decomposed=1
          ;stop
          cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        endfor
      ;printf,lun, input_str
      endfor
      
      
    endif ;end of line plots
    
    
  endfor
;FREE_LUN,lun
  
end
