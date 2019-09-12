

textfast,btl,/read,file_path='/home/nbarry/MWA/FHD/obs_list/beardsley_thesis_list.txt',string=1     ;1029
textfast,alltv,/read,file_path='/home/nbarry/MWA/FHD/obs_list/alltv.txt',string=1                   ;311
textfast,occ4,/read,file_path='/home/nbarry/MWA/FHD/obs_list/occupancy4.txt',string=1               ;70
textfast,occ3,/read,file_path='/home/nbarry/MWA/FHD/obs_list/occupancy3.txt',string=1               ;119
textfast,occ2,/read,file_path='/home/nbarry/MWA/FHD/obs_list/occupancy2.txt',string=1               ;182
textfast,allorbcomm,/read,file_path='/home/nbarry/MWA/FHD/obs_list/allorbcomm.txt',string=1         ;270
textfast,medorbcomm,/read,file_path='/home/nbarry/MWA/FHD/obs_list/medorbcomm.txt',string=1         ;132
textfast,highorbcomm,/read,file_path='/home/nbarry/MWA/FHD/obs_list/highorbcomm.txt',string=1       ;58

match, btl, alltv, suba, subb 
btl_zero=btl
btl_zero[suba]=0
inds = where(btl_zero NE 0, n_count)
btl_notv=btl_zero[inds]

test = occ4
match, btl_notv, test, suba, subb
btl_notv_zero = btl_notv
btl_notv_zero[suba]=0
inds = where(btl_notv_zero NE 0, n_count)
print, btl_notv_zero[inds], format='(D11.0)'

;;
btl_notv_noocc4 = btl_notv_zero[inds]
test = highorbcomm
match, btl_notv_noocc4, test, suba, subb
btl_notv_noocc4[suba]=0
inds = where(btl_notv_noocc4 NE 0, n_count)
print, btl_notv_noocc4[inds], format='(D11.0)'
