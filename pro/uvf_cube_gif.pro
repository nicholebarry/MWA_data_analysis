pro uvf_cube_gif

  filedir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_gaussdecomp_Fornax_dft/'
  cube_file = '1125953248_even_gridded_uvf.sav'

  restore, filedir + cube_file
 
  data_range=[10^(-1.),10^2.]

  for freq_i=0,191 do begin 
    plotfile=filedir+string(freq_i,format='(I03)')
    data = *model_uv_arr[0,freq_i]
    quick_image, abs(data), data_range=data_range, savefile=plotfile,/log

  endfor
  
end
