pro upgrade_map,

restore, 'diffuse_map_Byrne2021.sav'
full_hpx_inds = LINDGEN(3145728)
full_model = FLTARR(3145728) 

;freq=182e6 
;k_boltz = 1.38064852e-23
;vel_c = 299792458.0
;Jystr_to_mK = vel_c^2 / (1e23*2*freq^2*k_boltz) ;Jy/str to mK

;Jy/str doesn't depend on pixel area, so all good
full_model[hpx_inds] = (*model_arr[0])
ud_grade,full_model,interp_model,nside_out=2048,order_in='ring',order_out='ring',bad_data=0
interp_inds = where(interp_model NE 0)
(*model_arr[0]) = interp_model[interp_inds]

full_model[hpx_inds] = (*model_arr[1])
ud_grade,full_model,interp_model,nside_out=2048,order_in='ring',order_out='ring',bad_data=0
(*model_arr[1]) = interp_model[interp_inds] 

full_model[hpx_inds] = (*model_arr[2])
ud_grade,full_model,interp_model,nside_out=2048,order_in='ring',order_out='ring',bad_data=0
(*model_arr[2]) = interp_model[interp_inds]

full_model[hpx_inds] = (*model_arr[3])
ud_grade,full_model,interp_model,nside_out=2048,order_in='ring',order_out='ring',bad_data=0
(*model_arr[3]) = interp_model[interp_inds]
end
