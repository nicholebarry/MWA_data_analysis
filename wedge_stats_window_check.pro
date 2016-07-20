pro wedge_stats_window_check

  filename='/nfs/mwa-00/h1/nbarry/wedge_stats.txt'
  readcol, filename, obsids, jds, lsts, pointing, window_x, window_y, wedge_res_x,wedge_res_y,gal_wedge_x ,gal_wedge_y ,ptsrc_wedge_x ,ptsrc_wedge_y, $
    Format='(L,I12,I12,I12,F6.6,F6.6,I1,I1,I1,I1,I1,I1)'
  data2=read_ascii(filename);, delimiter=',')
  data = reform(data2.field01[0,*])
  window_ratio = window_y/window_x
  funky_lows = where(window_ratio LT .8)
  funky_extralows = where(window_ratio[funky_lows] GT .7)
  wh_obsids = long64(obsids[funky_lows])
  funky_obs = long64(wh_obsids[funky_extralows])
  stop
  
end