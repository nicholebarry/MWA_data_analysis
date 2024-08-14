pro analytic_degrid;,obs,psf,source_array,jones

filedir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_gaussbeam_1fits_bmt10000_dim34/'
restore, filedir+'beams/1061316296_beams.sav'
restore, filedir+'beams/1061316296_jones.sav'
restore, filedir+'metadata/1061316296_params.sav'
restore, filedir+'calibration/1061316296_cal.sav'
restore, '/fred/oz048/MWA/CODE/FHD/fhd_nb_gaussbeam_1fits_bmt10000_dim34/metadata/1061316296_obs.sav'
source_array = cal.skymodel.source_list  

frequency_array=(*obs.baseline_info).freq
i=Complex(0,1)
dimension=Long(obs.dimension)
elements=Long(obs.elements)

uu=params.uu
vv=params.vv
kbinsize=obs.kpix
kx_arr=temporary(uu)/kbinsize
ky_arr=temporary(vv)/kbinsize
nbaselines=obs.nbaselines
n_samples=obs.n_time
n_freq_use=N_Elements(frequency_array)
n_freq=Long(obs.n_freq)
n_freq_bin=N_Elements(freq_bin_i)
;psf_dim2=2*psf_dim
;group_arr=reform(psf.id[polarization,freq_bin_i,*])
;beam_arr=*psf.beam_ptr
n_sources = N_elements(source_array)

conj_i=where(ky_arr GT 0,n_conj)
IF n_conj GT 0 THEN BEGIN
    kx_arr[conj_i]=-kx_arr[conj_i]
    ky_arr[conj_i]=-ky_arr[conj_i]
ENDIF

n_freq = N_elements(frequency_array)
n_baselines = N_elements(kx_arr)
x_base=reform(Float(frequency_array#kx_arr),n_freq*n_baselines)/kbinsize;/dimension ;in relative pixel units
y_base=reform(Float(frequency_array#ky_arr),n_freq*n_baselines)/kbinsize;/elements ;in relative pixel units

x_vis=FLTARR(size(x_base,/dimension))
y_vis=FLTARR(size(x_base,/dimension))

source_array_use=Stokes_cnv(source_array,jones,/inverse,/no_extend,_Extra=extra)
x_flux=source_array_use.flux.XX
y_flux=source_array_use.flux.YY
x_source=2*!pi*(source_array_use.x - dimension/2.)/dimension ;in relative pixel units
y_source=2*!pi*(source_array_use.y - elements/2.)/elements ;in relative pixel units
x_beam=reform(((*obs.primary_beam_area[0]))#(FLTARR(size(kx_arr,/dimension))+1.),n_freq*n_baselines)
;phases for all sources per baseline
;for source_i=0, n_sources-1 do begin 

t1= Systime(1)
for source_i=0,9 do begin &$
  phase = x_base*x_source[source_i] + y_base*y_source[source_i] &$
  ;x_vis += x_flux[source_i]*x_beam*(Complex(cos(phase),sin(phase))) &$
  source_wave = x_flux[source_i]*(Complex(cos(phase),sin(phase))) &$
  ;for freq_i=0,n_freq-1 do source_wave[n_baselines*freq_i:(n_baselines*(freq_i+1))-1] *= (*obs.primary_beam_area[0])[freq_i] &$
  ;x_vis += source_wave &$
  x_vis += x_beam * source_wave &$
endfor 
print, Systime(1) - t1

x_vis = reform(x_vis,n_freq,n_baselines)
save,x_vis, filename='/fred/oz048/MWA/CODE/FHD/fhd_nb_gaussbeam_1fits_bmt10000_dim34/x_vis_999.sav'

end
