pro fivesource

; A polarized source off-zenith
id =[1,2,3,4,5]

ra = [359.8494,4.8494,354.8494,4.8494,354.8494] ; obsra= 359.8494 degrees

dec = [-26.78,-21.78,-21.78,-31.78,-31.78] ; obsdec = -26.783640 degrees

freq = 182.475; MHz, the center of the band.

q_frac = 0
u_frac = 0
v_frac = 0

flux_i = 2
;flux_q = flux_i * q_frac
;flux_u = flux_i * u_frac
;flux_v = flux_i * v_frac

catalog = source_comp_init(id=source_id, frequency=freq, ra=ra, dec=dec) 

catalog.flux.i = flux_i
;catalog.flux.q = flux_q
;catalog.flux.u = flux_u
;catalog.flux.v = flux_v

save_path=filepath('5source.sav',root=rootdir('FHD'),subdir='catalog_data')

; If the the specified file does not exist then procede with the save. If it does exist then stop and display error message.
if ~file_test(save_path) then begin
save, catalog, filename=save_path
endif else message, "Old file about to be overwritten! Make a new one, please."
end
