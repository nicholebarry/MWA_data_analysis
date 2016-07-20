PRO vis_res_seti_3D,save_point=save_point
  ; Script to plot visibility sigmas given various time interleaving or 2 second time visibilities differenced at various time steps

  ;Would like pointing information, LST information, row info

save_point='1250'

  ;Long run directory
  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  
  ;Use this as a proxy to get the obsids that were run
  vis_files=findfile(dir+'vis_data/106*_vis_XX.sav')
  obsids=file_basename(vis_files,'_vis_XX.sav')
  
  ;Cut obsids from long run, 32 hours
  ;filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/beardsley_thesis_list.txt'
  ;readcol, filename, obsids, format='A', /silent
  
  ;Need the frequency array, which is unchanging
  obs = GETVAR_SAVEFILE(dir+'metadata/1061316296_obs.sav', 'obs')
  freq_arr=(*obs.baseline_info).freq/10.^6.
  
  ;Set up histogram arrays
  all_3D_20=LONG(INTARR(384,56,3000))
  all_3D_40=LONG(INTARR(384,56,3000))
  
  ;Autos already flagged, but make sure for safety
  tiles = (*obs.baseline_info).tile_a - (*obs.baseline_info).tile_b
  autos = where(tiles EQ 0)
  
  undefine, tiles
  
  ;Read in histograms from a save point
  If keyword_set(save_point) then begin
  
    all_3D_20=getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/seti_all_3D_sec_20_'+save_point+'.sav','all_3D_20')
    all_3D_40=getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/seti_all_3D_sec_40_'+save_point+'.sav','all_3D_40')
    binned_diff=getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/seti_binned_diff_3D_'+save_point+'.sav','binned_diff')
    obs_i_start=UINT(save_point)+1
    
    IF ~keyword_set(all_3D_20) OR ~keyword_set(all_3D_40) OR ~keyword_set(binned_diff) then message, 'Uh oh'
    
  endif else obs_i_start=0
  
  ;****Restore the necessary information from the standard run to run this script outside of FHD.
  for obs_i=obs_i_start, N_elements(obsids)-1 do begin
    If (obs_i EQ 302)  OR (obs_i EQ 681) OR (obs_i EQ 927) OR (obs_i EQ 1266) OR(obs_i EQ 1393) then continue
    print, obs_i
    
    vis_XX = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_XX.sav', 'vis_ptr') ;restore array of calibrated visibilities
    
    vis_XX_model = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_model_XX.sav', 'vis_model_ptr') ;restore array of model visibilities
    
    vis_XX_res=*vis_XX-*vis_XX_model
    vis_XX_res[autos]=0
    
    undefine, vis_XX, vis_XX_model
    
    binsize=1.
    result=histogram(abs(vis_XX_res),binsize=binsize,locations=locations,omax=omax, /NAN,reverse_indices=ri)
    ;0bs_i = 301, 681, 1393 not recorded for 40sig test
    ;302, 681, 1393 not recorded for 2D test
    
    IF obs_i EQ 0 then begin
      binned_diff=uLong64(result)
    endif else begin
      loc_size=(size(locations))[1]
      binned_size=(size(binned_diff))[1]
      If loc_size LE binned_size then begin
        binned_diff[0:loc_size-1] = binned_diff[0:loc_size-1]+uLong64(result)
      endif else begin
        binned_diff[0:binned_size-1]=binned_diff[0:binned_size-1]+uLong64(result[0:binned_size-1])
        binned_diff=[binned_diff,uLong64(result[binned_size:loc_size-1])]
      endelse
    endelse
    
    ;outlier_min = 650 ; ~20sigma
    ;outlier_min = 1230 ; ~40sigma
    
    IF N_elements(result)-1 GE 650 then begin
      outliers=ri[ri[650]:ri[N_elements(result)]-1]
      s = SIZE(vis_XX_res)
      ncol = s[1]
      col = outliers MOD ncol
      row = outliers / ncol
      
      For loc_i = 0, N_elements(col)-1 do begin
        all_3D_20[col[loc_i],row[loc_i] mod 56,obs_i]=all_3D_20[col[loc_i],row[loc_i] mod 56,obs_i]+uLong64(1)
      endfor
      
    endif
    
    IF N_elements(result)-1 GE 1230 then begin
      outliers=ri[ri[1230]:ri[N_elements(result)]-1]
      s = SIZE(vis_XX_res)
      ncol = s[1]
      col = outliers MOD ncol
      row = outliers / ncol
      
      For loc_i = 0, N_elements(col)-1 do begin
        all_3D_40[col[loc_i],row[loc_i] mod 56,obs_i]=all_3D_40[col[loc_i],row[loc_i] mod 56,obs_i]+uLong64(1)
      endfor
      
    endif
    
    if (obs_i EQ 2200) OR (obs_i EQ 1100) OR (obs_i EQ 2000) OR (obs_i EQ 750) OR (obs_i EQ 1250) OR (obs_i EQ 1500) OR (obs_i EQ 1750) then begin
    
      save, binned_diff, filename='/nfs/eor-00/h1/nbarry/vis_res/seti_binned_diff_3D_'+strtrim(STRING(obs_i),2)+'.sav'
      save, all_3D_20, filename='/nfs/eor-00/h1/nbarry/vis_res/seti_all_3D_sec_20_'+strtrim(STRING(obs_i),2)+'.sav'
      save, all_3D_40, filename='/nfs/eor-00/h1/nbarry/vis_res/seti_all_3D_sec_40_'+strtrim(STRING(obs_i),2)+'.sav'
      
    endif
    
  endfor
  
  save, binned_diff, filename='/nfs/eor-00/h1/nbarry/vis_res/seti_binned_diff_3D_final.sav'
  save, all_3D_20, filename='/nfs/eor-00/h1/nbarry/vis_res/seti_all_3D_sec_20_final.sav'
  save, all_3D_40, filename='/nfs/eor-00/h1/nbarry/vis_res/seti_all_3D_sec_40_final.sav'
  
end