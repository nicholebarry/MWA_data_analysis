pro poly_comparisons_bftemp_slopeplotter, longrun=longrun, lst_compare=lst_compare, raw=raw, mean_tile=mean_tile, quick_check=quick_check, line_plots=line_plots, spread=spread, $
    make_obs_to_query=make_obs_to_query
  ;Program for plotting polyfit by pointings solutions over the longrun to check for stability. Longrun must be set at this time.
  ;lst_compare is not implemented at this time. Raw underplots the raw (minus bandpass by cable by pointing) for each obsid of the
  ;pointing or day and pointing depending on the next keyword.  by_pointing plots each all the pointings over all the days
  ;for a specific tile.  If that is not set, then plots are made for each pointing for each day.
  ;quick_check reads in only the -2 pointing to make plots quicker
  ;line_plots makes the plots described above
  ;spread is a new option to plot spreads of the zeroth order amp params by obs , overplotted with by pointing
    bypointingdiff=1
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
              stop
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
  
  
  
  stop
  ;****************************End of read in the polyfit calibration solutions, and make the raw gains minus bp sol if specified
  
  
  
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
  ;OPENW,lun,outdir+'linear_fits_days_fixed.txt',width=600
  ;header=['pointing','tile_name','pol','Group_1_days','Group_1_y-inter','Group_1_slope','Group_1_chi^2','Group_2_days','Group_2_y-inter','Group_2_slope','Group_2_chi^2']
  
  ;Using preexisting file to extract information about which tiles have which cable length
  mode_filepath=filepath('mwa_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  
  ;Taking tile information and cross-matching it with the nonflagged tiles array, resulting in nonflagged tile arrays
  ;grouped by cable length
  cable_length_ref=cable_len[Uniq(cable_len,Sort(cable_len))]
  n_cable=N_Elements(cable_length_ref)
  tile_use_arr=Ptrarr(n_cable)
  FOR cable_i=0,n_cable-1 DO tile_use_arr[cable_i]=Ptr_new(where((cable_len EQ cable_length_ref[cable_i])))
  
  slopes_per_cable=PTRARR(n_cable,/allocate)
  slopes_per_cable1=PTRARR(n_cable,/allocate)
  slopes_per_cable2=PTRARR(n_cable,/allocate)
  slope_diff=PTRARR(4,/allocate)
  
  textfast, data_array,read=1,string=1,file_path=outdir+'linear_fits.txt',delimiter=' '
  
  
  for pol_i=0,1 do begin
    if pol_i EQ 0 then pol_name='xx'
    if pol_i EQ 1 then pol_name='yy'
    
    
    poi_name=['-2','-1','0','1','2','3']
    If keyword_set(line_plots) then begin
    
      for tile_i=0,127 do begin
        input_str=STRARR(11,((size(pointing_num))[1]))
        
        savelocation=outdir+strtrim(string(tile_names[tile_i]),2)+'_scatter_'+pol_name+'.png'
        ;cgPS_Open,savelocation,/quiet,/nomatch
        title='Tile ' + strtrim(string(tile_names[tile_i]),2)+', gain temperature slope scatter, '+pol_name
        
        ;cgplot, [1,2,3],[1,2,3],xtitle='Average Gain (after bandpass removal)',ytitle='Beamformer Temperature (C)', $
        ;  title=title,charsize=1, psym=2, symsize=0.5,yrange=[0,50],xrange=[1,2],/NODATA
        
        ;for j=0,((size(pointing_num))[1])-1 do begin
        for j=2,2 do begin
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
          
          day_temp=day[input_index]
          day_sorted=day_temp[sorted_index]
          
          if pol_i EQ 0 then pol_con=0
          if pol_i EQ 1 then pol_con=640
          current_file_index=pol_con+1+5*tile_i+j
          
          file_days1=data_array[3,current_file_index]
          file_days2=data_array[7,current_file_index]
          file_days1_split=strsplit(file_days1,',',/EXTRACT)
          file_days2_split=strsplit(file_days2,',',/EXTRACT)
          
          If keyword_set(bycable) OR keyword_set(bycablebypol) then begin
            For group_i=5,9,4 do begin
              For cable_i=0,n_cable-1 do begin
                temp=where(tile_i EQ *tile_use_arr[cable_i],in_cable_group)
                
                If in_cable_group GT 0 then begin
                
                  IF abs(FLOAT(data_array[group_i,current_file_index])) GT 0 then begin
                  
                    If slopes_per_cable[cable_i] NE !NULL then *slopes_per_cable[cable_i]=[*slopes_per_cable[cable_i],FLOAT(data_array[group_i,current_file_index])] $
                    else *slopes_per_cable[cable_i]=FLOAT(data_array[group_i,current_file_index])
                    
                  endif
                endif
              endfor
            endfor
          endif
          
          If keyword_set(bypointingdiff) then begin
          
            IF abs(FLOAT(data_array[4,current_file_index])) GT 0 then begin
            
              diff_minusone=(FLOAT(data_array[4,current_file_index]))-(FLOAT(data_array[4,current_file_index-1]))
              diff_minustwo=(FLOAT(data_array[4,current_file_index]))-(FLOAT(data_array[4,current_file_index-2]))
              diff_plusone=(FLOAT(data_array[4,current_file_index]))-(FLOAT(data_array[4,current_file_index+1]))
              diff_plustwo=(FLOAT(data_array[4,current_file_index]))-(FLOAT(data_array[4,current_file_index+2]))
              temp_slopes=[diff_minustwo,diff_minusone,diff_plusone,diff_plustwo]
              
              For poi_i=0,3 do begin
                If slope_diff[poi_i] NE !NULL then *slope_diff[poi_i]=[*slope_diff[poi_i],temp_slopes[poi_i]] else *slope_diff[poi_i]=temp_slopes[poi_i]
              endfor
              
            endif
          endif
          
          
          
          ;For group_i=5,9,4 do begin
          For cable_i=0,n_cable-1 do begin
            temp=where(tile_i EQ *tile_use_arr[cable_i],in_cable_group)
            
            If in_cable_group GT 0 then begin
            
              group_i=5
              IF abs(FLOAT(data_array[group_i,current_file_index])) GT 0 then begin
                If slopes_per_cable1[cable_i] NE !NULL then *slopes_per_cable1[cable_i]=[*slopes_per_cable1[cable_i],FLOAT(data_array[group_i,current_file_index])] $
                else *slopes_per_cable1[cable_i]=FLOAT(data_array[group_i,current_file_index])
              endif
              
              group_i=9
              IF abs(FLOAT(data_array[group_i,current_file_index])) GT 0 then begin
                If slopes_per_cable2[cable_i] NE !NULL then *slopes_per_cable2[cable_i]=[*slopes_per_cable2[cable_i],FLOAT(data_array[group_i,current_file_index])] $
                else *slopes_per_cable2[cable_i]=FLOAT(data_array[group_i,current_file_index])
              endif
              
            endif
          endfor
        ;endfor
          
          
          
        ;input_str[0:2,j]=data_array[0:2,current_file_index]
        ;input_str[4:6,j]=data_array[4:6,current_file_index]
        ;input_str[8:10,j]=data_array[8:10,current_file_index]
          
        ;undefine, indexmap1, indexmap2
        ;for day_j=0, N_elements(file_days1_split)-1 do begin
        ;  index_temp=where(strmatch(day,file_days1_split[day_j]) EQ 1, n_count)
        ;  if keyword_set(indexmap1) then indexmap1=[indexmap1,index_temp] else indexmap1=index_temp
        ;
        ;  If file_days2_split[0] NE '0' then begin
        ;    If day_j LT N_elements(file_days2_split) then begin
        ;    index_temp=where(strmatch(day,file_days2_split[day_j]) EQ 1, n_count)
        ;    if keyword_set(indexmap2) then indexmap2=[indexmap2,index_temp] else indexmap2=index_temp
          
        ;    endif
        ;  endif
          
        ;endfor
          
        ;group1days=day_sorted[indexmap1]
        ;group1days_in=strjoin(group1days,',')
        ;If file_days2_split[0] NE '0' then begin
        ;  group2days=day_sorted[indexmap2]
        ;  group2days_in=strjoin(group2days,',')
        ;endif else group2days_in=0
          
        ;input_str[3,j]=group1days_in
        ;input_str[7,j]=group2days_in
          
        endfor
        
      ;printf,lun, input_str
      endfor
      
      
    endif ;end of line plots
    
    
  endfor
  ;FREE_LUN,lun
  
  If keyword_set(bycable) then begin
  
    y_arr_fill_ptr=PTRARR(n_cable,/allocate)
    x_arr_fill_ptr=PTRARR(n_cable,/allocate)
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/slope_hist_bycable_filled.png',/quiet,/nomatch
    color_array=['black','red','blue','purple','brown','orange']
    names=['90m','150m','230m','320m','400m','524m']
    thick_array=[2,2,2,2,2,2]
    cgplot, [1,2],[1,2],xrange=[-150,-30], yrange=[0,80],/NODATA, title='Slopes of bf temps vs. average gains per day, grouped by cable type', $
      xtitle='Fit slope value', ytitle='Density', charsize=1
      
    ;temp=plot( [1,2],[1,2],xrange=[-150,-40], yrange=[0,100],/NODATA)
      
    For cable_ii=0,n_cable-1 do begin
      IF cable_ii EQ 0 then cable_i=5
      IF cable_ii EQ 1 then cable_i=0
      IF cable_ii EQ 2 then cable_i=1
      IF cable_ii EQ 3 then cable_i=2
      IF cable_ii EQ 4 then cable_i=4
      IF cable_ii EQ 5 then cable_i=3
      
      ;Histogram slopes by cable length
      ;
      ;The binsize is arbitrary, but it does give a pretty plot. Could be increased. Must be even.
      binsize=4.
      
      ;Histogram the uv baselines to the binsize provided above, and record that in the reverse indices for applying
      ;to the visibility difference later. Convert the uv from light speed time to wavelengths.
      result=histogram(*slopes_per_cable[cable_i],binsize=binsize,reverse_indices=ri, /NAN, locations=locations, omax=omax)
      
      ;Calculate the visibility sigma (normalized due to real/imag) for each bin using reverse indices
      ;slopes_per_bin=FLTARR(N_elements(result))
      ;for i=0, N_elements(result)-1 do if result[i] GT 0 then slopes_per_bin[i]=(*slopes_per_cable[cable_i])[ri[ri[i]:ri[i+1]-1]]
      
      ;Make the x input and y input pretty for the plot
      y_arr=[result[0], result, result[N_elements(result)-1]]
      x_arr=[locations[0],locations+binsize/2,omax]
      
      undefine, y_arr_fill, x_arr_fill
      for i=0,N_elements(x_arr)-1 do begin
        IF i EQ 0 then begin
          y_arr_fill=0
          x_arr_fill=x_arr[i]-binsize/2
        endif
        y_arr_fill=[y_arr_fill,y_arr[i],y_arr[i]]
        x_arr_fill=[x_arr_fill,x_arr[i]-binsize/2,x_arr[i]+binsize/2]
        IF i EQ N_elements(x_arr)-1 then begin
          y_arr_fill=[y_arr_fill,0]
          x_arr_fill=[x_arr_fill,x_arr[i]+binsize/2]
        endif
      endfor
      *y_arr_fill_ptr[cable_i]=y_arr_fill
      *x_arr_fill_ptr[cable_i]=x_arr_fill
      
      ;Make
      ;cgPS_Open,'/nfs/eor-00/h1/nbarry/freq'+strtrim(freq,2)+'vissigma_inter_time'+strtrim(time_i,2)+'.png',/quiet,/nomatch
      cgoplot, x_arr, y_arr, psym=10, color=color_array[cable_i], thick=thick_array[cable_i]
      
    endfor
    stop
    For cable_ii=0,n_cable-1 do begin
      IF cable_ii EQ 0 then cable_i=5
      IF cable_ii EQ 1 then cable_i=0
      IF cable_ii EQ 2 then cable_i=1
      IF cable_ii EQ 3 then cable_i=2
      IF cable_ii EQ 4 then cable_i=4
      IF cable_ii EQ 5 then cable_i=3
      cgcolorfill, *x_arr_fill_ptr[cable_i], *y_arr_fill_ptr[cable_i], color=color_array[cable_i]
    endfor
    
    al_legend, ['90m','150m','230m','320m','400m','524m'],linestyle=[0,0,0,0,0,0],color=color_array, charsize=1
    
    
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif
  
  
  If keyword_set(bypointingdiff) then begin
  
    y_arr_fill_ptr=PTRARR(5,/allocate)
    x_arr_fill_ptr=PTRARR(5,/allocate)
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/slope_histdiff_bypointing.png',/quiet,/nomatch
    color_array=['black','red','blue','purple'];,'brown','orange']
    names=['Zen-MinusTwo','Zen-MinusOne','Zen-PlusOne','Zen-PlusTwo']
    thick_array=[3,3,3,3]
    cgplot, [1,2],[1,2],xrange=[-50,50], yrange=[0,90],/NODATA, title='Slopes of bf temps vs. ave gains per day, differenced from zenith', $
      xtitle='Fit slope value difference from zenith', ytitle='Density', charsize=1
      
    ;temp=plot( [1,2],[1,2],xrange=[-150,-40], yrange=[0,100],/NODATA)
      
    For poi_i=0,3 do begin
      ;IF cable_ii EQ 0 then cable_i=5
      ;IF cable_ii EQ 1 then cable_i=0
      ;IF cable_ii EQ 2 then cable_i=1
      ;IF cable_ii EQ 3 then cable_i=2
      ;IF cable_ii EQ 4 then cable_i=4
      ;IF cable_ii EQ 5 then cable_i=3
      
      ;Histogram slopes by cable length
      ;
      ;The binsize is arbitrary, but it does give a pretty plot. Could be increased. Must be even.
      binsize=3.
      
      ;Histogram the uv baselines to the binsize provided above, and record that in the reverse indices for applying
      ;to the visibility difference later. Convert the uv from light speed time to wavelengths.
      result=histogram(*slope_diff[poi_i],binsize=binsize,reverse_indices=ri, /NAN, locations=locations, omax=omax)
      
      ;Calculate the visibility sigma (normalized due to real/imag) for each bin using reverse indices
      ;slopes_per_bin=FLTARR(N_elements(result))
      ;for i=0, N_elements(result)-1 do if result[i] GT 0 then slopes_per_bin[i]=(*slopes_per_cable[cable_i])[ri[ri[i]:ri[i+1]-1]]
      
      ;Make the x input and y input pretty for the plot
      y_arr=[0, result, 0]
      x_arr=[locations[0],locations+binsize/2,omax]
      
      undefine, y_arr_fill, x_arr_fill
      for i=0,N_elements(x_arr)-1 do begin
        IF i EQ 0 then begin
          y_arr_fill=0
          x_arr_fill=x_arr[i]-binsize/2
        endif
        y_arr_fill=[y_arr_fill,y_arr[i],y_arr[i]]
        x_arr_fill=[x_arr_fill,x_arr[i]-binsize/2,x_arr[i]+binsize/2]
        IF i EQ N_elements(x_arr)-1 then begin
          y_arr_fill=[y_arr_fill,0]
          x_arr_fill=[x_arr_fill,x_arr[i]+binsize/2]
        endif
      endfor
      *y_arr_fill_ptr[poi_i]=y_arr_fill
      *x_arr_fill_ptr[poi_i]=x_arr_fill
      
      ;Make
      ;cgPS_Open,'/nfs/eor-00/h1/nbarry/freq'+strtrim(freq,2)+'vissigma_inter_time'+strtrim(time_i,2)+'.png',/quiet,/nomatch
      cgoplot, x_arr, y_arr, psym=10, color=color_array[poi_i], thick=thick_array[poi_i]
      
    endfor

    For poi_i=0,3 do begin
      ;IF cable_ii EQ 0 then cable_i=5
      ;IF cable_ii EQ 1 then cable_i=0
      ;IF cable_ii EQ 2 then cable_i=1
      ;IF cable_ii EQ 3 then cable_i=2
      ;IF cable_ii EQ 4 then cable_i=4
      ;IF cable_ii EQ 5 then cable_i=3
      ;cgcolorfill, *x_arr_fill_ptr[poi_i], *y_arr_fill_ptr[poi_i], color=color_array[poi_i]
    endfor
    
    al_legend, names,linestyle=[0,0,0,0],color=color_array, charsize=1
    
    
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif
  
  
  if keyword_set(bypolbycable) then begin
  
    y_arr_fill_ptr=PTRARR(n_cable,/allocate)
    x_arr_fill_ptr=PTRARR(n_cable,/allocate)
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/slope_hist_bycable_'+pol_name+'.png',/quiet,/nomatch
    color_array=['black','red','blue','purple','brown','orange']
    names=['90m','150m','230m','320m','400m','524m']
    thick_array=[2,2,2,2,2,2]
    cgplot, [1,2],[1,2],xrange=[-150,-30], yrange=[0,80],/NODATA, title='Slopes of bf temps vs. gains per day, grouped by cable type, '+pol_name, $
      xtitle='Fit slope value', ytitle='Density', charsize=1
      
    For cable_ii=0,n_cable-1 do begin
      IF cable_ii EQ 0 then cable_i=5
      IF cable_ii EQ 1 then cable_i=0
      IF cable_ii EQ 2 then cable_i=1
      IF cable_ii EQ 3 then cable_i=2
      IF cable_ii EQ 4 then cable_i=4
      IF cable_ii EQ 5 then cable_i=3
      
      ;The binsize is arbitrary, but it does give a pretty plot. Could be increased. Must be even.
      binsize=4.
      
      ;Histogram the uv baselines to the binsize provided above, and record that in the reverse indices for applying
      ;to the visibility difference later. Convert the uv from light speed time to wavelengths.
      result=histogram(*slopes_per_cable[cable_i],binsize=binsize,reverse_indices=ri, /NAN, locations=locations, omax=omax)
      
      ;Make the x input and y input pretty for the plot
      y_arr=[0, result, 0]
      x_arr=[locations[0],locations+binsize/2,omax]
      
      undefine, y_arr_fill, x_arr_fill
      for i=0,N_elements(x_arr)-1 do begin
        IF i EQ 0 then begin
          y_arr_fill=0
          x_arr_fill=x_arr[i]-binsize/2
        endif
        y_arr_fill=[y_arr_fill,y_arr[i],y_arr[i]]
        x_arr_fill=[x_arr_fill,x_arr[i]-binsize/2,x_arr[i]+binsize/2]
        IF i EQ N_elements(x_arr)-1 then begin
          y_arr_fill=[y_arr_fill,0]
          x_arr_fill=[x_arr_fill,x_arr[i]+binsize/2]
        endif
      endfor
      *y_arr_fill_ptr[cable_i]=y_arr_fill
      *x_arr_fill_ptr[cable_i]=x_arr_fill
      
      cgoplot, x_arr, y_arr, psym=10, color=color_array[cable_i], thick=thick_array[cable_i]
      
    endfor
    
    For cable_ii=0,n_cable-1 do begin
      IF cable_ii EQ 0 then cable_i=5
      IF cable_ii EQ 1 then cable_i=0
      IF cable_ii EQ 2 then cable_i=1
      IF cable_ii EQ 3 then cable_i=2
      IF cable_ii EQ 4 then cable_i=4
      IF cable_ii EQ 5 then cable_i=3
    ;cgcolorfill, *x_arr_fill_ptr[cable_i], *y_arr_fill_ptr[cable_i], color=color_array[cable_i]
    endfor
    
    al_legend, ['90m','150m','230m','320m','400m','524m'],linestyle=[0,0,0,0,0,0],color=color_array, charsize=1
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif
  
  
  
  stop
  
  y_arr_fill_ptr=PTRARR(n_cable,/allocate)
  x_arr_fill_ptr=PTRARR(n_cable,/allocate)
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/slope_hist_bycable_bygroup.png',/quiet,/nomatch
  color_array=['black','red','blue','purple','brown','orange']
  names=['90m','150m','230m','320m','400m','524m']
  thick_array=[2,2,2,2,2,2]
  cgplot, [1,2],[1,2],xrange=[-150,-30], yrange=[0,80],/NODATA, title='Slopes of bf temps vs. gains per day, grouped by cable type', $
    xtitle='Fit slope value', ytitle='Density', charsize=1
    
  For group_i=5,9,4 do begin
    If group_i EQ 5 then slopes_per_cable=slopes_per_cable1
    If group_i EQ 9 then slopes_per_cable=slopes_per_cable2
    For cable_ii=0,n_cable-1 do begin
      IF cable_ii EQ 0 then cable_i=5
      IF cable_ii EQ 1 then cable_i=0
      IF cable_ii EQ 2 then cable_i=1
      IF cable_ii EQ 3 then cable_i=2
      IF cable_ii EQ 4 then cable_i=4
      IF cable_ii EQ 5 then cable_i=3
      
      ;The binsize is arbitrary, but it does give a pretty plot. Could be increased. Must be even.
      binsize=4.
      
      ;Histogram the uv baselines to the binsize provided above, and record that in the reverse indices for applying
      ;to the visibility difference later. Convert the uv from light speed time to wavelengths.
      result=histogram(*slopes_per_cable[cable_i],binsize=binsize,reverse_indices=ri, /NAN, locations=locations, omax=omax)
      
      ;Make the x input and y input pretty for the plot
      y_arr=[0, result, 0]
      x_arr=[locations[0],locations+binsize/2,omax]
      
      undefine, y_arr_fill, x_arr_fill
      for i=0,N_elements(x_arr)-1 do begin
        IF i EQ 0 then begin
          y_arr_fill=0
          x_arr_fill=x_arr[i]-binsize/2
        endif
        y_arr_fill=[y_arr_fill,y_arr[i],y_arr[i]]
        x_arr_fill=[x_arr_fill,x_arr[i]-binsize/2,x_arr[i]+binsize/2]
        IF i EQ N_elements(x_arr)-1 then begin
          y_arr_fill=[y_arr_fill,0]
          x_arr_fill=[x_arr_fill,x_arr[i]+binsize/2]
        endif
      endfor
      *y_arr_fill_ptr[cable_i]=y_arr_fill
      *x_arr_fill_ptr[cable_i]=x_arr_fill
      
      cgoplot, x_arr, y_arr, psym=10, color=color_array[cable_i], thick=thick_array[cable_i]
      If group_i EQ 9 then cgcolorfill, *x_arr_fill_ptr[cable_i], *y_arr_fill_ptr[cable_i], color=color_array[cable_i]
      
    endfor
  endfor
  
  For cable_ii=0,n_cable-1 do begin
    IF cable_ii EQ 0 then cable_i=5
    IF cable_ii EQ 1 then cable_i=0
    IF cable_ii EQ 2 then cable_i=1
    IF cable_ii EQ 3 then cable_i=2
    IF cable_ii EQ 4 then cable_i=4
    IF cable_ii EQ 5 then cable_i=3
  ;cgcolorfill, *x_arr_fill_ptr[cable_i], *y_arr_fill_ptr[cable_i], color=color_array[cable_i]
  endfor
  
  al_legend, ['90m','150m','230m','320m','400m','524m'],linestyle=[0,0,0,0,0,0],color=color_array, charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  stop
end
