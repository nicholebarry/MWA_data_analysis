pro parse

filename='/nfs/eor-00/h1/nbarry/pointing_temp.csv'
;OpenR, lun, filename, /Get_Lun 
;readf, lun, obs, pointing
;Free_lun, lun

readcol, filename, obs, pointing, format='LL,I'
obscount=N_elements(obs)


OpenW, lun0, '/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23zenith.txt', /Get_Lun 
OpenW, lun2, '/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23plusone.txt', /Get_Lun 
OpenW, lun4, '/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23plustwo.txt', /Get_Lun 
OpenW, lun6, '/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23plusthree.txt', /Get_Lun 
OpenW, lun1, '/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23minusone.txt', /Get_Lun 
OpenW, lun3, '/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23minustwo.txt', /Get_Lun 

FOR i=0, obscount-1 DO BEGIN
  If pointing[i] eq 0 THEN printf, lun0, obs[i]
  If pointing[i] eq 2 THEN printf, lun2, obs[i]
  If pointing[i] eq 4 THEN printf, lun4, obs[i]
  If pointing[i] eq 6 THEN printf, lun6, obs[i]
  If pointing[i] eq 1 THEN printf, lun1, obs[i]
  If pointing[i] eq 3 THEN printf, lun3, obs[i]
ENDFOR

Free_lun, lun0, lun1, lun2, lun3, lun4, lun6


END