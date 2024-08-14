pro time_flagger

dir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_int/'

;obs_ids = ['1061312152','1061312272','1061312640','1061317400']
;obs_ids = ['1061317400','1061317520','1061317640','1061317888','1061318128','1061318248','1061318376','1061318496','1061318736','1061318864']

  dir_obsids = '/home/nbarry/MWA/pipeline_scripts/bash_scripts/ozstar/obs_list/'
  day = '1'
  parsednames=[day+'_minustwo',day+'_minusone',day+'_zenith',day+'_plusone',day+'_plustwo']
;  parsednames=[day+'_minustwo']
  n_pointings = N_elements(parsednames)
  obs_ptr=PTRARR(n_pointings,/allocate_heap)
  n_obs_per_pointing = FLTARR(n_pointings)

  ; Get observation ids and put them in a pointing pointer
  FOR poi_i=0,n_pointings-1 DO BEGIN
    readcol, dir_obsids + parsednames[poi_i]+'.txt', obs_temp, format='A', /silent
    *obs_ptr[poi_i]=obs_temp
    n_obs_per_pointing[poi_i] = N_elements(obs_temp)
  ENDFOR

  n_obs_max = max(n_obs_per_pointing)
  n_obs_total = total(n_obs_per_pointing)
  mean_per_time = FLTARR(n_pointings,n_obs_max,56)
  mean_per_time2 = FLTARR(n_pointings,n_obs_max,56)
  mean_per_time_unflagged = FLTARR(n_pointings,n_obs_max,56)

;restore_mean=1
  if ~keyword_set(restore_mean) then begin

  for poi_i=0, n_pointings-1 do begin

    obs_ids = *obs_ptr[poi_i]
    n_obs = N_elements(obs_ids)


    for obs_i=0,n_obs-1 do begin

      vis_cal = getvar_savefile(dir + 'vis_data/'+obs_ids[obs_i]+'_vis_XX.sav','vis_ptr')
      obs = getvar_savefile(dir + 'metadata/'+obs_ids[obs_i]+'_obs.sav','obs')
      params = getvar_savefile(dir + 'metadata/'+obs_ids[obs_i]+'_params.sav','params')
      freq_use = where((*obs.baseline_info).freq_use EQ 1)
      freq_use = where((*obs.baseline_info).freq_use EQ 1)
temp = where((abs((*vis_cal)[freq_use,*])) EQ 0,n_count)
print, n_count
       mean_amp = mean(abs((*vis_cal)[freq_use,*]),dimension=1,/nan)

       zeros = where(mean_amp EQ 0,n_count)
       if n_count GT 0 then mean_amp[zeros] = !Values.F_NAN
       nbaselines = obs.nbaselines
       n_time = obs.n_time

       for time_i=0, n_time-2 do begin

         range = [time_i * nbaselines, (time_i+1) * nbaselines - 1]
         inds_cross = where((params.uu[range[0]:range[1]]^2 + params.vv[range[0]:range[1]]^2) NE 0)
         mean_per_time[poi_i,obs_i,time_i] = mean(mean_amp[range[0]+inds_cross],/nan)

      endfor
stop
    endfor



    mean_per_time_unflagged[poi_i,*,*] = mean_per_time[poi_i,*,*]
    mean_per_time2[poi_i,*,*] = mean_per_time[poi_i,*,*]

  endfor
  
  zeros = where(mean_per_time EQ 0) 
  mean_per_time[zeros] = !Values.F_NAN
  mean_per_time2[zeros] = !Values.F_NAN
stop
  endif else restore, "time_flagger_restore.sav"

  ;Expect thermal noise to be consistent within a pointing
  time_stddev = FLTARR(n_pointings)
  time_median = FLTARR(n_pointings)
  for poi_i=0,n_pointings-1 do begin
    time_stddev[poi_i] = (stddev(mean_per_time[poi_i,*,*],/nan))
    time_median[poi_i] = (median(mean_per_time[poi_i,*,*]))
  endfor

  cutoff = time_median + time_stddev
  cutoff2 = time_median + time_stddev*2
;for obs_i=0, n_obs-1 do begin
;  flag = where(mean_per_time[obs_i,*] GT cutoff[obs_i],n_count)
;  flag2 = where(mean_per_time[obs_i,*] GT cutoff2[obs_i],n_count2)
;  if n_count GT 0 then mean_per_time[obs_i,flag] = !Values.F_NAN
;  if n_count2 GT 0 then mean_per_time2[obs_i,flag2] = !Values.F_NAN
;endfor

flag = PTRARR(n_pointings,/allocate)
flag2 = PTRARR(n_pointings,/allocate)

for poi_i=0,n_pointings-1 do begin
  *flag[poi_i] = where(mean_per_time[poi_i,*,*] GT cutoff[poi_i],n_count)
  *flag2[poi_i] = where(mean_per_time[poi_i,*,*] GT cutoff2[poi_i],n_count2)
  if n_count GT 0 then begin
    temp = reform(mean_per_time[poi_i,*,*])   
    temp[*flag[poi_i]] = !Values.F_NAN
    mean_per_time[poi_i,*,*] = temp
  endif
  if n_count2 GT 0 then begin
    temp = reform(mean_per_time2[poi_i,*,*])   
    temp[*flag2[poi_i]] = !Values.F_NAN
    mean_per_time2[poi_i,*,*] = temp
  endif
endfor


cgps_open,'pointing_time_flag.png',/quiet,/nomatch

position = [0.10, 0.70, 0.45, 0.92]
cgplot, mean_per_time_unflagged[0,*,*], yrange=[20,60], xrange=[0,60],charsize=.9,Position=position, /NoErase,/NoData, title='Pointing -2' 
for obs_i=0,N_elements(*obs_ptr[0])-1 do cgoplot, reform(mean_per_time_unflagged[0,obs_i,*]),XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4
cgoplot, [0,56], [cutoff[0],cutoff[0]],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',thickness=6
cgoplot, [0,56], [cutoff2[0],cutoff2[0]],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',linestyle=2,thickness=6

position = [0.5, 0.70, 0.85, 0.92]
cgplot, mean_per_time_unflagged[1,*,*], yrange=[20,60], xrange=[0,60],charsize=.9,Position=position, /NoErase,/NoData, title='Pointing -1' 
for obs_i=0,N_elements(*obs_ptr[1])-1 do cgoplot, reform(mean_per_time_unflagged[1,obs_i,*]),XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4
cgoplot, [0,56], [cutoff[1],cutoff[1]],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',thickness=6
cgoplot, [0,56], [cutoff2[1],cutoff2[1]],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',linestyle=2,thickness=6

position = [0.10, 0.4, 0.45, 0.620]
cgplot, mean_per_time_unflagged[2,*,*], yrange=[20,60], xrange=[0,60],charsize=.9,Position=position, /NoErase,/NoData, title='Pointing 0' 
for obs_i=0,N_elements(*obs_ptr[2])-1 do cgoplot, reform(mean_per_time_unflagged[2,obs_i,*]),XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4
cgoplot, [0,56], [cutoff[2],cutoff[2]],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',thickness=6
cgoplot, [0,56], [cutoff2[2],cutoff2[2]],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',linestyle=2,thickness=6

position = [0.5, 0.4, 0.85, 0.620]
cgplot, mean_per_time_unflagged[3,*,*], yrange=[20,60], xrange=[0,60],charsize=.9,Position=position, /NoErase,/NoData, title='Pointing +1' 
for obs_i=0,N_elements(*obs_ptr[3])-1 do cgoplot, reform(mean_per_time_unflagged[3,obs_i,*]),XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4
cgoplot, [0,56], [cutoff[3],cutoff[3]],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',thickness=6
cgoplot, [0,56], [cutoff2[3],cutoff2[3]],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',linestyle=2,thickness=6

position = [0.10, 0.1, 0.45, 0.320]
cgplot, mean_per_time_unflagged[4,*,*], yrange=[20,60], xrange=[0,60],charsize=.9,Position=position, /NoErase,/NoData, title='Pointing +2' 
for obs_i=0,N_elements(*obs_ptr[4])-1 do cgoplot, reform(mean_per_time_unflagged[4,obs_i,*]),XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4
cgoplot, [0,56], [cutoff[4],cutoff[4]],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',thickness=6
cgoplot, [0,56], [cutoff2[4],cutoff2[4]],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',linestyle=2,thickness=6

all_cutoff = median(mean_per_time) +stddev(mean_per_time,/nan)
all_cutoff2 = median(mean_per_time) +2*stddev(mean_per_time,/nan)

position = [0.5, 0.1, 0.85, 0.320]
cgplot, mean_per_time_unflagged[0,*,*], yrange=[20,60], xrange=[0,60],charsize=.9,Position=position, /NoErase,/NoData, title='All pointings'
for poi_i=0, n_pointings-1 do for obs_i=0,N_elements(*obs_ptr[poi_i])-1 do cgoplot, reform(mean_per_time_unflagged[poi_i,obs_i,*]),XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4
cgoplot, [0,56], [all_cutoff,all_cutoff],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',thickness=6
cgoplot, [0,56], [all_cutoff2,all_cutoff2],XTICKFORMAT="(A1)", YTICKFORMAT="(A1)",Position=position,/NoErase,xstyle=4,ystyle=4,color='red',linestyle=2,thickness=6

cgtext, .39,.025, 'Time step index',/Normal,charsize=1.2
cgtext, .05,.287, 'Amplitude average across array, Jy',/Normal,charsize=1.2, orientation=90.
cglegend, Title=['Obs','$\sigma$','2$\sigma$'], color=['black','red','red'], linestyle=[0,0,2], location = [.885,.55], charsize=1,Length=0.03


cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage

stop


end
