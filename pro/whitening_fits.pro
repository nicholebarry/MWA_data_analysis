pro whitening_fits,dir=dir,fits=fits,obsids=obsids,deep2_inds=deep2_inds,cal_cut_inds=cal_cut_inds,dir2=dir2, overplot=overplot, phase=phase

  dir='/nfs/mwa-03/r1/EoR2013/fhd_nb_whitening/calibration/'
  
  dir2 = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/'
  
  
  ;Things to set per run
  save_path=dir+'mode_params_stdcompare/'
  ;text_title=' tile, M!Io,a,1flag!N'
  text_title=' tile, M!Ilow!N'
  
  ;if keyword_set(overplot) then text_title=' tile, M!Io,std!N and M!Ia,phasesep,extraremove,crosspoly!N'
  if keyword_set(overplot) then text_title=' tile, M!Ihigh!N and M!Ilow!N'
  
  cal_files_text='_cal.sav'
  if keyword_set(dir2) then cal_files2_text='_cal.sav'
  cal_dir_added=1
  
  
  
  file_path='~/MWA/IDL_code/obs_list/whitening_filter_zenith_subset.txt
  textfast,obsids,header,file_path=file_path,/read,/string
  
  obsids_name=obsids
  obsids=ulong(obsids)
  nobs=n_elements(obsids)
  ntile=128
  npol=2
  
  fits = fltarr(npol,nobs,ntile,3) ; last index for mode_i,amp,phase
  
  for obs_i=0,nobs-1 do begin
    print,obs_i
    cal_files=dir+obsids_name+cal_files_text
    cal = getvar_savefile(cal_files[obs_i],'cal')
    for pol_i=0,npol-1 do begin
      for tile_i=0,ntile-1 do begin
        if (cal.mode_params[pol_i,tile_i] ne !null) then begin
          fits[pol_i,obs_i,tile_i,*] = (*cal.mode_params[pol_i,tile_i])
        endif
      endfor
    endfor
  endfor
  file_mkdir, save_path, /NOEXPAND_PATH
  save,filename=save_path+'fits_save.sav',fits
  ;endelse
  
  file_path='~/MWA/IDL_code/obs_list/Aug23zenith.txt
  textfast,obsids_name,header,file_path=file_path,/read,/string
  
  if keyword_set(dir2) then begin
    fits2 = fltarr(npol,nobs,ntile,3) ; last index for mode_i,amp,phase
    cal_files2=dir2+obsids_name+cal_files2_text
    If keyword_set(cal_dir_added) then cal_files2=dir2+'calibration/'+obsids_name+cal_files2_text
    ;if file_test(dir2+'mode_params_plots/fits_save.sav') then begin
    ;  restore,dir2+'mode_params_plots/fits_save.sav'
    ;endif else begin
    for obs_i=0,nobs-1 do begin
      print,obs_i
      cal = getvar_savefile(cal_files2[obs_i],'cal')
      for pol_i=0,npol-1 do begin
        for tile_i=0,ntile-1 do begin
          if (cal.mode_params[pol_i,tile_i] ne !null) then begin
            fits2[pol_i,obs_i,tile_i,*] = (*cal.mode_params[pol_i,tile_i])
          endif
        endfor
      endfor
    endfor
  ;endelse
    
  endif ;end dir2 if
  
  tile_control=19
  tile_change=[20,21,22,23]
  
  ; get tile info
  restore, '/nfs/eor-03/r1/EoR2013/fhd_nb_std_test_twopolyquad/metadata/1061316296_obs.sav'
  tile_names=(*obs.baseline_info).tile_names
  tile_inds=where(fits2[0,0,*,0] ne 0)
  
  
  mode_filepath=filepath(obs.instrument+'_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  tile_cable_inds=where(cable_len EQ 150)
  match, tile_inds, tile_cable_inds, suba
  tile_inds=tile_inds[suba]
  
  binsize=.1
  result_mode=histogram(fits[*,*,tile_change,0],binsize=binsize,locations=locations,omax=omax)
  y_arr_mode=[0, result_mode, 0]
  x_arr_mode=[locations[0]-binsize/2,locations+binsize/2,locations[N_elements(locations)-1]+binsize]
  
  binsize=.1
  result_mode2=histogram(fits[*,*,tile_control,0],binsize=binsize,locations=locations,omax=omax)
  y_arr_mode2=[0, result_mode2, 0]
  x_arr_mode2=[locations[0]-binsize/2,locations+binsize/2,locations[N_elements(locations)-1]+binsize]
  
  binsize=.001
  result_amp=histogram(fits[*,*,tile_change,1],binsize=binsize,locations=locations,omax=omax)
  y_arr_amp=[0, result_amp, 0]
  x_arr_amp=[locations[0]-binsize/2,locations+binsize/2,locations[N_elements(locations)-1]+binsize]
  
  binsize=.001
  result_amp2=histogram(fits[*,*,tile_control,1],binsize=binsize,locations=locations,omax=omax)
  y_arr_amp2=[0, result_amp2, 0]
  x_arr_amp2=[locations[0]-binsize/2,locations+binsize/2,locations[N_elements(locations)-1]+binsize]
  
stop

cgPS_Open,'/nfs/eor-00/h1/nbarry/whitening_amp_zenith.png',/quiet,/nomatch
cgplot, x_arr_amp2, y_arr_amp2, psym=10,title='Distribution of fit amplitudes for zenith pointing',xtitle='Amplitude (fraction of normalized gain)', ytitle='Density', charsize=1,yrange=[0,15], xrange=[.005,.03]
cgoplot,x_arr_amp, y_arr_amp, psym=10,color='blue'
cglegend, title=['20 Oct 2015 - Tile 34','20 Oct 2015 - Tile 35,36,37,38'], color=['black','blue'],length=.05,location=[.155,.85],charsize=1
cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage  
  
cable_vf=.81
c_light=299792458.
obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/1061316176_obs.sav','obs')
freq_arr=(*obs.baseline_info).freq
n_freq=384.
bandwidth=(Max(freq_arr)-Min(freq_arr))*n_freq/(n_freq-1)
x_arr_mode2*=(1/bandwidth)*(c_light*cable_vf/2.)
x_arr_mode*=(1/bandwidth)*(c_light*cable_vf/2.)

cgPS_Open,'/nfs/eor-00/h1/nbarry/whitening_mode_zenith.png',/quiet,/nomatch
cgplot, x_arr_mode2, y_arr_mode2, psym=10,title='Distribution of fit modes for zenith pointing',xtitle='Mode (effective cable length)', ytitle='Density', charsize=1, yrange=[0,25]
cgoplot,x_arr_mode, y_arr_mode, psym=10,color='blue'
cglegend, title=['20 Oct 2015 - Tile 34','20 Oct 2015 - Tile 35,36,37,38'], color=['black','blue'],length=.05,location=[.5,.85],charsize=1
cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage  
  
end