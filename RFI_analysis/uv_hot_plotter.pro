pro UV_hot_plotter

  restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/antenna_check/hot_baselines/seti_tiles20_thesis_evenodd_total.sav'
  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  params = GETVAR_SAVEFILE(dir+'metadata/1061316296_params.sav', 'params')
  obs = GETVAR_SAVEFILE(dir+'metadata/1061316296_obs.sav', 'obs')
  uu = params.uu
  vv = params.vv
  tile_a = (*obs.baseline_info).tile_a
  tile_b = (*obs.baseline_info).tile_b
  
  hot = getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/antenna_check/hot_baselines_3000.sav','hot')
  hot_col = hot mod 129 ;tileA
  hot_row = hot / 129 ;tileB
  
  longrun_names_match, obs_names=obs_names
  restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/antenna_check/hot_baselines/seti_tiles20_thesis_evenodd_total.sav'
  tiles_20_totalobs = ULONG64(INTARR(129,129))
  for obs_i=0, 1028 do if (obs_names[obs_i,1] NE 'Oct31') then tiles_20_totalobs = tiles_20_totalobs +tiles20_total[*,*,obs_i]
  
  vis_mask = FLTARR(N_elements(tile_a))
  for hot_i=0, N_elements(hot)-1 do begin
    ;print, 'hot = ' + strtrim(hot_i,2)
    where_col = where(tile_a EQ hot_col[hot_i])
    where_row = where(tile_b[where_col] EQ hot_row[hot_i])
    if hot_i EQ 0 then mean_uu_baseline = mean(uu[where_col[where_row]]) else mean_uu_baseline = [mean_uu_baseline,mean(uu[where_col[where_row]])]
    if hot_i EQ 0 then mean_vv_baseline = mean(vv[where_col[where_row]]) else mean_vv_baseline = [mean_vv_baseline,mean(vv[where_col[where_row]])]
    vis_mask[where_col[where_row]] = 1.
  endfor
  
  uu_masked = uu * vis_mask * 299792458.
  vv_masked = vv * vis_mask * 299792458.
  
  n_bins = 100.
  uu_range = minmax(mean_uu_baseline * 299792458.)
  vv_range = minmax(mean_vv_baseline * 299792458.)
  
  pixels = FLTARR(n_bins,n_bins)
  binsize = (uu_range[1]-uu_range[0])/n_bins
  uu_result = histogram(mean_uu_baseline * 299792458., binsize=binsize, min=uu_range[0], reverse_indices=ri_uu, nbins=n_bins)
  ;  n_bins_dec = Ceil((dec_range_sculptor[0]-dec_range_sculptor[1])/binsize)

  for bin_i=0,n_bins-2 do begin
    if ri_uu[bin_i] EQ ri_uu[bin_i+1] then continue
    uu_inds = ri_uu[ri_uu[bin_i]:ri_uu[bin_i+1]-1]
    vv_result = histogram(mean_vv_baseline[uu_inds] * 299792458., binsize=binsize, min=vv_range[0], reverse_indices=ri_vv, nbins=n_bins)
    for bin_j=0,N_elements(vv_result)-2 do begin
      if ri_vv[bin_j] EQ ri_vv[bin_j+1] then continue
      vv_inds = ri_vv[ri_vv[bin_j]:ri_vv[bin_j+1]-1]
      pixels[bin_i,bin_j] = total(tiles_20_totalobs[hot_col[uu_inds[vv_inds]],hot_row[uu_inds[vv_inds]]])
    endfor
  endfor
  
  pixel_centers_uu = FINDGEN(n_bins)*binsize+uu_range[0]+binsize/2.
  pixel_centers_vv = FINDGEN(n_bins)*binsize + vv_range[0] + binsize/2.
  
  stop
im_range=[24*n_bins/50,26*n_bins/50]
  quick_image, congrid(pixels[im_range[0]:im_range[1],im_range[0]:im_range[1]],n_bins*20,n_bins*20), congrid(pixel_centers_uu[im_range[0]:im_range[1]],n_bins*20), congrid(pixel_centers_vv[im_range[0]:im_range[1]],n_bins*20)
end