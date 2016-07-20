pro vis_res_seti_timing

  restore,'/nfs/eor-00/h1/nbarry/seti_binned_diff_40sig_1000.sav'
  restore,'/nfs/eor-00/h1/nbarry/seti_all_col_40sig_1000.sav'
  restore,'/nfs/eor-00/h1/nbarry/seti_all_row_40sig_1000.sav'
  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  obs = GETVAR_SAVEFILE(dir+'metadata/1061316296_obs.sav', 'obs')
  
  times=where((*obs.baseline_info).tile_b EQ 1, n_count)
  time_duration=times[1]-1
  time_ordering=ULONG64(INTARR(time_duration+1))
  for time_i=0, time_duration do begin
    tot_slot=ULONG64(0)
    for slot_i=0, N_elements(times)-1 do begin
      tot_slot = tot_slot+all_row[time_i+time_i*slot_i]
    endfor
    time_ordering[time_i]=tot_slot
  endfor
  
  stop
  
end