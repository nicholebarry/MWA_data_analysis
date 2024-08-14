pro gain_res_diff, obsid, type, diff=diff, line_plot_multiplier=line_plot_multiplier, cable_type=cable_type, aug27=aug27
  ;cable_type is a num 0 through 5 (0:90m,1:150m,etc)

  beg_type=type
  
  IF keyword_set(cable_type) then begin
    cable_beg=cable_type
    cable_end=cable_type
  endif else begin
    cable_beg=0
    cable_end=5
  endelse
  
  IF (N_elements(obsid) EQ 1) THEN BEGIN
    IF (obsid EQ 'obs_id_23') THEN BEGIN
      obs_arr=strarr(94)
      filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/obs_id_23.txt'
      OpenR, lun, filename, /Get_Lun
      readf, lun, obs_arr
      Free_lun, lun
      obsid = obs_arr
    ENDIF
  ENDIF
  
  obsid_tot=(size(obsid))[1]
  IF (size(obsid))[0] EQ 0 then obsid_tot=1
  IF ~keyword_set(line_plot_multiplier) then line_plot_multiplier=1.
  
  for obs_i=0,obsid_tot-1 do begin
  
    type=beg_type
    input_obs=strtrim(string(obsid[obs_i]),2)
    
    If ~keyword_set(aug27) then begin
    
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/metadata/' + obsid[obs_i] + '_obs.sav'
      
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/calibration/' + obsid[obs_i] + '_cal.sav'
      no_cable_cal=(*cal.gain_residual[0])
      
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_std/calibration/' + obsid[obs_i] + '_cal.sav'
      cable_cal=(*cal.gain_residual[0])
      
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_iterative/calibration/' + obsid[obs_i] + '_cal.sav'
      iterative=(*cal.gain_residual[0])
      
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_iterative_cablethroughout/calibration/' + obsid[obs_i] + '_cal.sav'
      cable_cal_iterative=(*cal.gain_residual[0])
      
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_pointing_gainrescor/calibration/' + obsid[obs_i] + '_cal.sav'
      saverun=(*cal.gain_residual[0])
      
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_min_cal_5lambda_with_diffuse/calibration/' + obsid[obs_i] + '_cal.sav'
      min_cal_cable=(*cal.gain_residual[0])
      
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_min_cal_5lambda_with_diffuse/calibration/' + obsid[obs_i] + '_cal.sav'
      min_cal_no_cable=(*cal.gain_residual[0])
      
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_predigjump/calibration/' + obsid[obs_i] + '_cal.sav'
      predigjump_no_cable=(*cal.gain_residual[0])
      
    endif else begin
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_Aug27/calibration/' + obsid[obs_i] + '_cal.sav'
      aug27_cable_cal=(*cal.gain_residual[0])
      If obsid[obs_i] EQ '1061656320' then obsid_comp='1061311664'
      If obsid[obs_i] EQ '1061660952' then obsid_comp='1061316296'
      If obsid[obs_i] EQ '1061667664' then obsid_comp='1061323008'
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_Aug27/metadata/' + obsid[obs_i] + '_obs.sav'
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_Aug27/calibration/' + obsid[obs_i] + '_cal.sav'
      no_cable_cal=(*cal.gain_residual[0])
    endelse
    
    case type of
      'cable_cal':cable_cal=cable_cal
      'no_cable_cal':cable_cal=no_cable_cal ;Only makes sense in diff keyword plots
      'iterative':cable_cal=iterative
      'cable_cal_iterative':cable_cal=cable_cal_iterative
      'saverun':cable_cal=saverun
      'min_cal_cable':cable_cal=min_cal_cable
      'min_cal_no_cable':cable_cal=min_cal_no_cable
      'aug27':cable_cal=aug27_cable_cal
      'dig':cable_cal=predigjump_no_cable
      ELSE: print, 'Only cable_cal, iterative, cable_cal_iterative, saverun, min_cal_cable, or min_cal_no_cable offered. Running cable_cal as default'
    endcase
    
    If type EQ 'abs' then begin
      ;cable_cal=abs(cable_cal)
      cable_cal=abs(saverun)
      no_cable_cal=abs(no_cable_cal)
      ;type='abs_cable_cal'
      type='abs_saverun'
      diff='abs_no_cable_cal'
    endif
    
    If keyword_set(diff) then begin
      case diff of
        'cable_cal':no_cable_cal=cable_cal
        'iterative':no_cable_cal=iterative
        'cable_cal_iterative':no_cable_cal=cable_cal_iterative
        'saverun':no_cable_cal=saverun
        'min_cal_cable':no_cable_cal=min_cal_cable
        'min_cal_no_cable':no_cable_cal=min_cal_no_cable
        'none':no_cable_cal[*,*]=0
        ELSE: print, 'Only cable_cal, iterative, cable_cal_iterative, min_cal_cable, min_cal_no_cable, or saverun offered to be differenced from besides no_cable_cal (default). Running no_cable_cal as default'
      endcase
      type=diff+'-'+type
    endif else begin
      if obs_i EQ 0 then type='no_cable_cal-'+type
    endelse
    
    if obs_i EQ 0 then begin
      input_type = '_'+type
      name_type = ' ' +type
    endif
    
    n_cable=6
    
    mode_filepath=filepath(obs[0].instrument+'_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
    textfast,data_array,/read,file_path=mode_filepath,first_line=1
    cable_len=Reform(data_array[2,*])
    
    ;Initialize arrays for excluding flagged tiles
    tile_use_cable=PTRARR(n_cable, /allocate)
    freq_use=where((*obs.baseline_info).freq_use)
    
    ;Group cable lengths in arrays
    *tile_use_cable[0]=where((*obs.baseline_info).tile_use AND cable_len EQ 90)
    *tile_use_cable[1]=where((*obs.baseline_info).tile_use AND cable_len EQ 150)
    *tile_use_cable[2]=where((*obs.baseline_info).tile_use AND cable_len EQ 230)
    *tile_use_cable[3]=where((*obs.baseline_info).tile_use AND cable_len EQ 320)
    *tile_use_cable[4]=where((*obs.baseline_info).tile_use AND cable_len EQ 400)
    *tile_use_cable[5]=where((*obs.baseline_info).tile_use AND cable_len EQ 524)
    
    half_waterfall=1
    
    If keyword_set(half_waterfall) then begin
    
      for cable_i=cable_beg,cable_end do begin
        If cable_i eq 0 then begin
          cable_name='90'
          cable_tiles=*tile_use_cable[0]
        endif
        If cable_i eq 1 then begin
          cable_name='150'
          cable_tiles=*tile_use_cable[1]
        endif
        If cable_i eq 2 then begin
          cable_name='230'
          cable_tiles=*tile_use_cable[2]
        endif
        If cable_i eq 3 then begin
          cable_name='320'
          cable_tiles=*tile_use_cable[3]
        endif
        If cable_i eq 4 then begin
          cable_name='400'
          cable_tiles=*tile_use_cable[4]
        endif
        If cable_i eq 5 then begin
          cable_name='524'
          cable_tiles=*tile_use_cable[5]
        endif
        
        temp_no_cable=FLTARR(N_elements(cable_tiles),384)
        temp_cable=FLTARR(N_elements(cable_tiles),384)
        x_range=INDGEN(((size(cable_tiles))[1]+1))
        y_range=INDGEN(N_elements(freq_use)+1)
        line_array=INTARR(((size(cable_tiles))[1]+1))
        line_array[*]=223
        For freq_i=0, 383 do temp_no_cable[*,freq_i]=no_cable_cal[freq_i,cable_tiles]
        For freq_i=0, 383 do temp_cable[*,freq_i]=cable_cal[freq_i,cable_tiles]
        
        cgPS_Open,'/nfs/eor-00/h1/nbarry/cable_residuals/'+cable_name+'m_waterfall_'+input_obs+input_type+'.png',/quiet,/nomatch
        quick_image, temp_no_cable[*,freq_use]-temp_cable[*,freq_use],x_range,y_range,data_range=[-4.*10^(-2.),4.*10^(-2.)], $
          xtitle='Tile Index (out of cable group)', ytitle='Frequency Index (unflagged only)', title=cable_name+'m gain res diff '+input_obs+name_type,charsize=1
        cgoplot, line_array, Linestyle=2
        cgText, .74,.655,/Normal, 'Digital Gain Jump', color='black', charsize=1
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
      endfor
      
    endif
    
    ;full_waterfall=1
    
    if keyword_set(full_waterfall) then begin
    
      for cable_i=cable_beg,cable_end do begin
        If cable_i eq 0 then begin
          cable_name='90'
          cable_tiles=*tile_use_cable[0]
        endif
        If cable_i eq 1 then begin
          cable_name='150'
          cable_tiles=*tile_use_cable[1]
        endif
        If cable_i eq 2 then begin
          cable_name='230'
          cable_tiles=*tile_use_cable[2]
        endif
        If cable_i eq 3 then begin
          cable_name='320'
          cable_tiles=*tile_use_cable[3]
        endif
        If cable_i eq 4 then begin
          cable_name='400'
          cable_tiles=*tile_use_cable[4]
        endif
        If cable_i eq 5 then begin
          cable_name='524'
          cable_tiles=*tile_use_cable[5]
        endif
        
        temp_no_cable=FLTARR(128,384)
        temp_no_cable[*,*]=-1
        temp_cable=FLTARR(128,384)
        temp_cable[*,*]=0
        x_range=INDGEN((128+1))
        y_range=INDGEN(N_elements(freq_use)+1)
        line_array=INTARR((128+1))
        line_array[*]=223
        For freq_i=0, 383 do temp_no_cable[cable_tiles,freq_i]=no_cable_cal[freq_i,cable_tiles]
        For freq_i=0, 383 do temp_cable[cable_tiles,freq_i]=cable_cal[freq_i,cable_tiles]
        
        cgPS_Open,'/nfs/eor-00/h1/nbarry/cable_residuals/'+cable_name+'m_full_waterfall_'+input_obs+input_type+'.png',/quiet,/nomatch
        quick_image, temp_no_cable[*,freq_use]-temp_cable[*,freq_use],x_range,y_range,data_range=[-4.*10^(-2.),4.*10^(-2.)], $
          xtitle='Tile Index (out of all tiles)', ytitle='Frequency Index (unflagged only)', title=cable_name+'m gain res difference '+input_obs+name_type,charsize=1
        cgoplot, line_array, Linestyle=2
        cgText, .74,.655,/Normal, 'Digital Gain Jump', color='black', charsize=1
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
      endfor
      
      temp_no_cable=FLTARR(128,384)
      temp_no_cable[*,*]=-1
      temp_cable=FLTARR(128,384)
      temp_cable[*,*]=0
      x_range=INDGEN((128+1))
      y_range=INDGEN(N_elements(freq_use)+1)
      line_array=INTARR((128+1))
      line_array[*]=223
      For freq_i=0, 383 do temp_no_cable[*,freq_i]=no_cable_cal[freq_i,*]
      For freq_i=0, 383 do temp_cable[*,freq_i]=cable_cal[freq_i,*]
      
      cgPS_Open,'/nfs/eor-00/h1/nbarry/cable_residuals/all_cables_full_waterfall_'+input_obs+input_type+'.png',/quiet,/nomatch
      quick_image, temp_no_cable[*,freq_use]-temp_cable[*,freq_use],x_range,y_range,data_range=[-4.*10^(-2.),4.*10^(-2.)], $
        xtitle='Tile Index (out of all tiles)', ytitle='Frequency Index (unflagged only)', title='Gain res difference, all tiles '+input_obs+name_type,charsize=1
      cgoplot, line_array, Linestyle=2
      cgText, .74,.655,/Normal, 'Digital Gain Jump', color='black', charsize=1
      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
    endif
    
    line_plot=1
    
    if keyword_set(line_plot) then begin
    
      for cable_i=cable_beg,cable_end do begin
        If cable_i eq 0 then begin
          cable_name='90'
          cable_tiles=*tile_use_cable[0]
          yrange=[-.05,1.85]*line_plot_multiplier
        endif
        If cable_i eq 1 then begin
          cable_name='150'
          cable_tiles=*tile_use_cable[1]
          yrange=[-.05,3.05]*line_plot_multiplier
        endif
        If cable_i eq 2 then begin
          cable_name='230'
          cable_tiles=*tile_use_cable[2]
          yrange=[-.05,2.15]*line_plot_multiplier
        endif
        If cable_i eq 3 then begin
          cable_name='320'
          cable_tiles=*tile_use_cable[3]
          yrange=[-.05,.75]*line_plot_multiplier
        endif
        If cable_i eq 4 then begin
          cable_name='400'
          cable_tiles=*tile_use_cable[4]
          yrange=[-.05,1.65]*line_plot_multiplier
        endif
        If cable_i eq 5 then begin
          cable_name='524'
          cable_tiles=*tile_use_cable[5]
          yrange=[-.05,2.65]*line_plot_multiplier
        endif
        
        line_array_y=(Findgen(50)/49)*(yrange[1]-yrange[0])+yrange[0]
        line_array_x=INTARR(50)
        line_array_x[*]=223
        
        cgPS_Open,'/nfs/eor-00/h1/nbarry/cable_residuals/'+cable_name+'m_line_plot_'+input_obs+input_type+'.png',/quiet,/nomatch
        cgplot, no_cable_cal[freq_use,cable_tiles[0]]-cable_cal[freq_use,cable_tiles[0]],xrange=[0,336],yrange=yrange, $
          xtitle='Frequency Index (unflagged only)', ytitle='Relative Gain Residual Difference', title=cable_name+'m gain res difference '+input_obs+name_type,charsize=1
        for tile_i=1, N_elements(cable_tiles)-1 do cgoplot, no_cable_cal[freq_use,cable_tiles[tile_i]]-cable_cal[freq_use,cable_tiles[tile_i]]+.1*tile_i*line_plot_multiplier
        cgoplot, line_array_x,line_array_y, Linestyle=2
        ;cgText, .6,.05,/Normal, 'DGJ', color='black', charsize=1
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
      endfor
      
    endif
    
    mean_plot=1
    
    if keyword_set(mean_plot) then begin
    
      for cable_i=cable_beg,cable_end do begin
        If cable_i eq 0 then begin
          cable_name='90'
          cable_tiles=*tile_use_cable[0]
        endif
        If cable_i eq 1 then begin
          cable_name='150'
          cable_tiles=*tile_use_cable[1]
        endif
        If cable_i eq 2 then begin
          cable_name='230'
          cable_tiles=*tile_use_cable[2]
        endif
        If cable_i eq 3 then begin
          cable_name='320'
          cable_tiles=*tile_use_cable[3]
        endif
        If cable_i eq 4 then begin
          cable_name='400'
          cable_tiles=*tile_use_cable[4]
        endif
        If cable_i eq 5 then begin
          cable_name='524'
          cable_tiles=*tile_use_cable[5]
        endif
        
        line_array_y=(Findgen(50)/49)*(1)+(-.1)
        line_array_x=INTARR(50)
        line_array_x[*]=223
        
        cgPS_Open,'/nfs/eor-00/h1/nbarry/cable_residuals/'+cable_name+'m_mean_plot_'+input_obs+input_type+'.png',/quiet,/nomatch
        temp_cable=cable_cal[freq_use, *]
        temp_no_cable=no_cable_cal[freq_use,*]
        if beg_type EQ 'abs' then begin
          cgplot, mean(( temp_no_cable[*,cable_tiles] - temp_cable[*,cable_tiles]),dimension=2), title=cable_name+'m gain res diff mean '+input_obs+name_type, xtitle='Frequency Index (only unflagged)', ytitle='Relative Gain Residual Difference',xrange=[0,336],yrange=[-.015,.015],charsize=1
        endif else begin
          cgplot, mean(abs( temp_no_cable[*,cable_tiles] - temp_cable[*,cable_tiles]),dimension=2), title=cable_name+'m gain res diff mean '+input_obs+name_type, xtitle='Frequency Index (only unflagged)', ytitle='Relative Gain Residual Difference',xrange=[0,336],charsize=1
        endelse
        cgoplot, line_array_x,line_array_y, Linestyle=2
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
      endfor
      
    endif
    
    
  endfor ;end obs for
  
end