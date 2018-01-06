pro vis_noise_plotter

  poi_names=['Aug23_minustwo','Aug23_minusone','Aug23_zenith','Aug23_plusone','Aug23_plustwo']
  fileroot='/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/metadata/'
  poi_color=['red','blue','green','aqua','purple']
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_sigma.png',/quiet,/nomatch
  cgplot, [1,1],[1,1], /nodata, yrange=[0,60],xrange=[0,384], xtitle='Freq index',ytitle='vis noise (Jy)',charsize=1
  
  for poi_i=0, 4 do begin
    textfast, obsids,/read,/string, file_path='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/'+poi_names[poi_i]+'.txt'
    
    for obs_i=0, N_elements(obsids)-1 do begin
      if ~file_test(fileroot+obsids[obs_i]+'_obs.sav') then continue
      restore, fileroot+obsids[obs_i]+'_obs.sav'
      cgoplot, (*obs.vis_noise)[0,*], color=poi_color[poi_i]
    endfor  
    
  endfor
  cglegend, title=['-2','-1','0','1','2'],color=poi_color, location=[.7,.85],charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
end