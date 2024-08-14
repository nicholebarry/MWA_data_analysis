pro remake_obslists

days_to_convert = ['10_minustwo','10_minusone','10_zenith','10_plusone','10_plustwo']
n_days = N_elements(days_to_convert)

restore,'../stats_10_11_12_13_14.sav'
flagged_inds = where(pf_full LT .4)
unflagged_inds = where(pf_full GT .4)
obs = obsid_full[unflagged_inds]
;obs = string(obs, format='(I10)')

for day_i=0, n_days-1 do begin
  readcol,days_to_convert[day_i]+'.txt',obs_from_txt,format='(A)'
  
  match,obs,obs_from_txt,suba,subb
  print, N_elements(subb)/float(N_elements(obs_from_txt))
  openw,21,days_to_convert[day_i]+'_ssins.txt'
  printf,21, transpose(obs_from_txt[subb]),format='(A)'
  close, 21
;  writecol,days_to_convert[day_i]+'_ssins.txt',obs_from_txt[subb]

endfor


end
