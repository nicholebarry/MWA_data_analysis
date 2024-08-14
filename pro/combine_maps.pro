pro combine_maps

catalog_dir='/fred/oz048/nbarry/catalogs/'

restore, catalog_dir + 'diffuse_map_Byrne2021_interp2048.sav' 

;Default ordering is RING when written from pygdsm 
;;             Global sky model in healpix format, with NSIDE=1024. Output map
;;            is in galactic coordinates, ring format.
fitsfast, gsm_arr, action='read', file_path=catalog_dir+'gsm_182mhz_Jysr_nomono_nogalaxy_2048.fits'
fitsfast, gsm_ra, action='read', file_path=catalog_dir+'gsm_182mhz_Jysr_nomono_nogalaxy_2048_ra.fits'
fitsfast, gsm_dec, action='read', file_path=catalog_dir+'gsm_182mhz_Jysr_nomono_nogalaxy_2048_dec.fits'
gsm_nside=2048
gsm_inds = LINDGEN(nside2npix(gsm_nside))
;gsm_arr = getvar_savefile(catalog_dir+'gsm2016_182MHz_nomono_interp2048.sav','MODEL_ARR')
;gsm_inds = getvar_savefile('gsm2016_182MHz_nomono_interp2048.sav','hpx_inds')
;gsm_nside = getvar_savefile('gsm2016_182MHz_nomono_interp2048.sav','nside')

;pix2vec_ring,gsm_nside,gsm_inds,gsm_coords 
;vec2ang,gsm_coords,gsm_dec,gsm_ra,/astro
;glactc,gsm_ra,gsm_dec,2000.,gsm_ra,gsm_dec,2, /degree ;dec seems to differ between these formulation and python
 
;pix2ang_ring, 1024, gsm_inds, gsm_gtheta, gsm_gphi 
;gsm_gb = 0.5 * !pi - gsm_theta
;glactc,gsm_phi,gsm_dec,2000.,gsm_gphi,gsm_gtheta,2
;gsm_theta = 0.5 * !pi - gsm_dec

;FHD seems to assume ring ordering, so assume Ruby gave me ring.
pix2ang_ring, nside, hpx_inds, byrne_theta, byrne_phi
byrne_dec = (0.5 * !pi - byrne_theta)*(180./!pi)
byrne_ra = byrne_phi *(180./!pi)

step_size = 2.
n_dec = 180 / step_size
n_ra = 360 / step_size
for dec_i=0,n_dec-1 do begin

byrne_d_inds = where((byrne_dec+90 GE dec_i*step_size) AND $
                     (byrne_dec+90 LE dec_i*step_size+step_size), n_bd_count)

for ra_i=0,n_ra-1 do begin

  if n_bd_count GT 0 then begin
    byrne_r_inds = wherE((byrne_ra[byrne_d_inds] GE ra_i*step_size) AND $
                         (byrne_ra[byrne_d_inds] LE ra_i*step_size+step_size), n_br_count)
  endif else n_br_count=0

  if (n_br_count NE 0) then begin

;if loop has reached here, then byrne is defined in that square degree
;remove gsm

gsm_r_inds = wherE((gsm_ra GE ra_i*step_size) AND $
                   (gsm_ra LE ra_i*step_size+step_size) AND $
                   (gsm_dec+90 GE dec_i*step_size) AND $
                   (gsm_dec+90 LE dec_i*step_size+step_size), n_r_count)
;gsm_rd_inds = wherE((gsm_dec[gsm_r_inds]+90 GE dec_i*step_size) AND (gsm_dec[gsm_r_inds]+90 LE dec_i*step_size+step_size), n_rd_count)


;stop
;if N_elements(valid_inds) NE 0 then valid_inds = [valid_inds,gsm_inds[gsm_r_inds[gsm_rd_inds]]] else valid_inds = gsm_inds[gsm_r_inds[gsm_rd_inds]]

if n_r_count GT 0 then gsm_inds[gsm_r_inds]=-1

endif
;result = histogram(byrne_ra,binsize=1,reverse_indices=byrne_ri)
;temp = where(result EQ 0)
;;defines the min and max where byrne is not defined
;ra_range= minmax(temp *1.)
;;defines the min and max where byrne is defined
;dec_range =  minmax(byrne_dec)

;gsm_ra_inds = where((gsm_ra GT ra_range[0]) AND (gsm_ra LT ra_range[1]))

;gsm_dec_inds = where((gsm_dec LT dec_range[0]) OR (gsm_dec GT dec_range[1]))
;valid_inds = [gsm_dec_inds,gsm_ra_inds]

endfor
endfor
stop
;valid_inds = valid_inds[uniq(valid_inds, sort(valid_inds))]

;this is still in galactic projection
hpx_inds = valid_inds
ra_vals = gsm_ra[valid_inds] ;these are celestial projection
dec_vals = gsm_dec[valid_inds] ;these are celestial projection
model_arr = gsm_arr[valid_inds]
nside=2048
units='kelvin'
coord_sys='galactic'
save, model_arr,hpx_inds,nside, units, coord_sys, filename='gsm2016_182MHz_nomono_interp2048_subset.sav'
fitsfast, model_arr, file_path='/fred/oz048/nbarry/catalogs/gsm2016_182MHz_nomono_interp2048_subset_I.fits', action='write'
fitsfast, ra_vals, file_path='/fred/oz048/nbarry/catalogs/gsm2016_182MHz_nomono_interp2048_subset_ra.fits', action='write'
fitsfast, dec_vals, file_path='/fred/oz048/nbarry/catalogs/gsm2016_182MHz_nomono_interp2048_subset_dec.fits', action='write'
fitsfast, hpx_inds, file_path='/fred/oz048/nbarry/catalogs/gsm2016_182MHz_nomono_interp2048_subset_inds.fits', action='write'

ra_vals = byrne_ra ;these are celestial projection
dec_vals = byrne_dec ;these are celestial projection
fitsfast, *model_arr[0], file_path='/fred/oz048/nbarry/catalogs/byrne_interp2048_I.fits', action='write'
fitsfast, *model_arr[1], file_path='/fred/oz048/nbarry/catalogs/byrne_interp2048_Q.fits', action='write'
fitsfast, *model_arr[2], file_path='/fred/oz048/nbarry/catalogs/byrne_interp2048_U.fits', action='write'
fitsfast, *model_arr[3], file_path='/fred/oz048/nbarry/catalogs/byrne_interp2048_V.fits', action='write'
fitsfast, ra_vals, file_path='/fred/oz048/nbarry/catalogs/byrne_interp2048_ra.fits', action='write'
fitsfast, dec_vals, file_path='/fred/oz048/nbarry/catalogs/byrne_interp2048_dec.fits', action='write'

end
