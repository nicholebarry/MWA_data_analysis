

pars = DBLARR(12,5,17)
for i =0,11 do pars[i,*,*] = reform((*psf.beam_gaussian_params[0,i]),5,17) 
.compile /home/nbarry/MWA/FHD/mwa_beam_gaussian_decomp.pro
.compile /home/nbarry/MWA/FHD/beam_gaussian_decomp.pro

cen=280
pad=1.3
ce = (cen*pad)/2.
psf_image_dim=2800

for i=0,9 do begin $
psf_base_single = gaussian_decomp(FINDGEN(psf_image_dim),FINDGEN(psf_image_dim),reform(pars[i,*,*],5.*17.),model_npix=cen*pad) & $
quick_image, psf_base_single[1200:1600,1200:1600],/log,data_range=[10^(-9.),1.],savefile='/home/nbarry/MWA/image_0'+strtrim(i,2) & $
endfor
for i=10,11 do begin $
psf_base_single = gaussian_decomp(FINDGEN(psf_image_dim),FINDGEN(psf_image_dim),reform(pars[i,*,*],5.*17.),model_npix=cen*pad) & $
quick_image, psf_base_single[1200:1600,1200:1600],/log,data_range=[10^(-9.),1.],savefile='/home/nbarry/MWA/image_'+strtrim(i,2) & $
endfor


cgplot, psf.freq[*]/1e6, pars[*,2,0]/pars[0,2,0], xtitle='freq', ytitle='relative params val', title='x sigma', yrange=[.8,1.1]
for i=1,11 do cgoplot, psf.freq[*]/1e6, pars[*,2,i]/pars[0,2,i]
cgoplot, psf.freq[*]/1e6, psf.freq[0]/psf.freq[*],color='blue'

cgplot, psf.freq[*]/1e6, pars[*,4,0]/pars[0,4,0], xtitle='freq', ytitle='relative params val', title='y sigma', yrange=[.8,1.1]
for i=1,11 do cgoplot, psf.freq[*]/1e6, pars[*,4,i]/pars[0,4,i]
cgoplot, psf.freq[*]/1e6, psf.freq[0]/psf.freq[*],color='blue'

cgplot, psf.freq[*]/1e6, (ce-pars[*,1,0])/(ce-pars[0,1,0]), xtitle='freq', ytitle='relative params val', title='x offset', yrange=[.8,1.1]
for i=1,11 do cgoplot, psf.freq[*]/1e6, (ce-pars[*,1,i])/(ce-pars[0,1,i])
cgoplot, psf.freq[*]/1e6, psf.freq[0]/psf.freq[*],color='blue'

lobes = [5:17]
lobes = [0,lobes]

cgplot, psf.freq[*]/1e6, pars[*,2,0]/pars[0,2,0], xtitle='freq', ytitle='relative params val', title='x sigma', yrange=[.8,1.1]
for i=1,12 do cgoplot, psf.freq[*]/1e6, pars[*,2,lobes[i]]/pars[0,2,lobes[i]]
cgoplot, psf.freq[*]/1e6, psf.freq[0]/psf.freq[*],color='blue'


quick_image, (*(*psf.image_info).image_power_beam_arr[0,11])[1200:1600,1200:1600],/log,data_range=[10^(-9.),1.],window=1
quick_image, psf_base_single[1200:1600,1200:1600],/log,data_range=[10^(-9.),1.],window=2
mwa_beam_gaussian_decomp, (cen*pad)/2., psf.pix_horizon, parinfo=parinfo, parvalues=p, freq=psf.freq[11]
psf_base_single_theory = gaussian_decomp(FINDGEN(psf_image_dim),FINDGEN(psf_image_dim),p,model_npix=cen*pad)
quick_image, psf_base_single_theory[1200:1600,1200:1600],/log,data_range=[10^(-9.),1.],window=3
