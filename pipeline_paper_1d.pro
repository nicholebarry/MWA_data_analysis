pro pipeline_paper_1d

  ;limit_plot=1
  if keyword_set(limit_plot) then begin

    cut = 'Aug23_longrunstyle'
    chans='ch0-127_'
    win='fullimg_'
    spec_window='bh_'
    basefile='/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/Combined_obs_'+cut+'_cubeXX__even_odd_joint_'+win+chans
    cubes=['dirty','model','res']
    pols=['yy']
    averemove='averemove_'
    endfile='_'+averemove+'sw'+spec_window+'dencorr_kperplambda10-50_1dkpower.idlsave'

    ;limit_percent=0.9772 ;2 sigma
    limit_percent = 0.97725
    hubble_param=0.71

    n_cubes=1
    color_num = [10,12,14,16]
    rgbcolors = [[137,117,202],[114,166,89],[73, 142, 217],[197,120,62]]
    ;73, 142, 217 little brighter blue, 54,107,163 little darker blue
    ;rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101]];,[198,118,63]]

    ;ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,28,38,39,40,41,42,43,54,53,55]+4

    ;*ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4

    ;for n_i=0, N_elements(color_num)-1 do $
    ;  TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]
    ;
    ;color_array = [10B,12B,14B,16B]
    color_array=['black','royal blue','firebrick','purple','green']

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 5 else thickness = 3



    file=basefile+cubes[2]+'_'+pols+endfile
    restore,file
    ; find the best limit in this file
    n_k=n_elements(k_edges)
    k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    delta=power*(k^3.)/(2.*!pi^2.)

    power_model = getvar_savefile(basefile+cubes[1]+'_'+pols+endfile,'power')
    delta_model=power_model*(k^3.)/(2.*!pi^2.)
    power_dirty = getvar_savefile(basefile+cubes[0]+'_'+pols+endfile,'power')
    delta_dirty=power_dirty*(k^3.)/(2.*!pi^2.)

    dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
    dsigma[0] = !Values.F_INFINITY

    limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf((delta)/dsigma/sqrt(2)))*sqrt(2))+delta
    limits_abs=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(abs(delta)/dsigma/sqrt(2)))*sqrt(2))+abs(delta)
    lim=min(limits,ind)
    header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut+' '+chans+' '+cubes+' '+win+' '+pols+')'
    print,header

    zeros = where(delta EQ 0,n_count)
    delta[zeros] = !Values.F_INFINITY
    delta[0] = !Values.F_INFINITY
    delta_model[0] = !Values.F_INFINITY
    delta_dirty[0] = !Values.F_INFINITY
    ;delta[*ranges[0]] = !Values.F_INFINITY
    dsigma[0] = !Values.F_INFINITY
    ;dsigma[*ranges[0]] = !Values.F_INFINITY
    limits[0] = !Values.F_INFINITY
    ;limits[*ranges[0]] = !Values.F_INFINITY
    k=k/hubble_param

    ;fiducial theory
    restore,'/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave'
    delta_eor=power*(k_centers^3.)/(2.*!pi^2.)
    k_centers=k_centers/hubble_param

    if keyword_set(pdf_plot) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/pipe_1d.pdf',/quiet,/nomatch



    xtitle_use='k (h Mpc$\up-1$)'
    ytitle_use='$\Delta$$\up2$ (mK$\up2$)'


    low_range=10^1.
    low_x_range=.08


    cgplot, k,delta_dirty,/xlog,/ylog,psym=10, xrange=[.008,1.9], yrange=[10^1.,10^9.],ytitle=ytitle_use, $
      xtitle=xtitle_use, charsize =1.1, color=color_array[0],thick=thickness,/noerase, title='1D Diagnostic Power Spectrum'
    cglegend,title=['calibrated data yy','model yy','residual yy','1$\sigma$ thermal noise','EoR signal'],$
      color=[color_array[0],color_array[1],color_array[2],color_array[3],color_array[4]], location=[.17,.5], $
      charsize=1.1, thick=thickness,linestyle=[0,0,0,2,0],length=.05
    ;cgoplot, k,delta,/xlog,/ylog,psym=10,color=color_array[1],thick=thickness
    ;cgoplot, k,delta_model,/xlog,/ylog,psym=10,color=color_array[2],thick=thickness

    error_low = (dsigma*2.); < (delta);*.99999
    lows = where(delta LT 0, n_count)
    if n_count GT 0 then error_low[lows]=abs(delta[lows])
    lows = where((abs(delta) - error_low) LT low_range,n_count)
    if n_count GT 0 then error_low[lows]=abs(delta[lows])-low_range
    error_high = dsigma*2.
    lows = where(k LT low_x_range,n_count)
    if n_count GT 0 then error_high[lows]=!Values.F_INFINITY
    if n_count GT 0 then limits[lows]=!Values.F_INFINITY
    delta_k=k[1]/2.
    for k_i=0, n_k-2 do $
      cgColorFill, [k[k_i]-delta_k, k[k_i]+delta_k, k[k_i]+delta_k,k[k_i]-delta_k], $
      [limits[k_i], limits[k_i], abs(delta[k_i])-error_low[k_i],abs(delta[k_i])-error_low[k_i]], $
      Color='grey',/checkforfinite

    cgoplot, k,delta_dirty,/xlog,/ylog,psym=10,color=color_array[0],thick=thickness
    cgoplot, k,delta,/xlog,/ylog,psym=10,color=color_array[2],thick=thickness
    cgoplot, k,delta_model,/xlog,/ylog,psym=10,color=color_array[1],thick=thickness

    ;cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=color_use,thick=thickness,position=position_use
    cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_array[3],thick=thickness
    ;cgoplot, k,delta,/xlog,/ylog,color='black',psym=10,thick=thickness-1,position=position_use
    cgoplot, k_centers,delta_eor,/xlog,/ylog,color=color_array[4],thick=thickness-1,position=position_use



    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage

    stop
  endif

  filepath_justEoR = '/Users/nabarry/MWA/data/fhd_nb_sim_BH2_BH2_justEoR/Combined_obs_zenith_cubeXX__even_odd_joint_fullimg_res_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave'
  filepath_recover = '/Users/nabarry/MWA/data/fhd_nb_sim_BH2_BH2_signalloss_noforegrounds_over/Combined_obs_zenith_cubeXX__even_odd_joint_fullimg_res_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave'
  filepath_foregrounds = '/Users/nabarry/MWA/data/fhd_nb_sim_BH2_BH2_signalloss_exactforegrounds_6000/Combined_obs_zenith_cubeXX__even_odd_joint_fullimg_res_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave'
  filepath_calerrors = '/Users/nabarry/MWA/data/fhd_nb_sim_BH2_BH2_signalloss_exactforegrounds_6000_over/Combined_obs_zenith_cubeXX__even_odd_joint_fullimg_res_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave'
filepath_foregrounds_coarse = '/Users/nabarry/MWA/data/fhd_nb_sim_BH2_BH2_signalloss_exactforegrounds_6000_coarse/Combined_obs_zenith_cubeXX__even_odd_joint_fullimg_res_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave'

  restore,filepath_justEoR
  ; find the best limit in this file
  n_k=n_elements(k_edges)
  k_justEoR=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
  delta_justEoR=power*(k_justEoR^3.)/(2.*!pi^2.)

  restore,filepath_recover
  ; find the best limit in this file
  n_k=n_elements(k_edges)
  k_recover=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
  delta_recover=power*(k_recover^3.)/(2.*!pi^2.)

  restore,filepath_foregrounds
  ; find the best limit in this file
  n_k=n_elements(k_edges)
  k_foregrounds=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
  delta_foregrounds=power*(k_foregrounds^3.)/(2.*!pi^2.)

  restore,filepath_calerrors
  ; find the best limit in this file
  n_k=n_elements(k_edges)
  k_calerrors=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
  delta_calerrors=power*(k_calerrors^3.)/(2.*!pi^2.)
  
  restore,filepath_foregrounds_coarse
  ; find the best limit in this file
  n_k=n_elements(k_edges)
  k_foregrounds_coarse=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
  delta_foregrounds_coarse=power*(k_foregrounds_coarse^3.)/(2.*!pi^2.)
  
  restore,'/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave'
  delta_eor=power*(k_centers^3.)/(2.*!pi^2.)
      hubble_param=0.71

k_centers=k_centers/hubble_param
k_foregrounds=k_foregrounds/hubble_param
k_calerrors=k_calerrors/hubble_param
k_recover=k_recover/hubble_param
k_justEoR=k_justEoR/hubble_param

  xtitle_use='k (h Mpc$\up-1$)'
  ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
  color_array=['forest green','firebrick','blue','goldenrod','purple']
  pdf_plot=1
  if keyword_set(pdf_plot) then thickness = 6 else thickness = 3

  if keyword_set(pdf_plot) then cgPS_Open,'/Users/nabarry/MWA/data/sim_delta_update.pdf',/quiet,/nomatch

  cgplot, k_foregrounds,delta_foregrounds,/xlog,/ylog,psym=10, xrange=[.06,1.5], yrange=[10^0.,10^7.],ytitle=ytitle_use, $
    xtitle=xtitle_use, charsize =1.1, color=color_array[0],thick=thickness,/noerase, title='In-situ simulation', aspect=.5
    
    ;cgoplot, k_foregrounds_coarse,delta_foregrounds_coarse,/xlog,/ylog,psym=10,color='forest green',thick=thickness, linestyle=3
    cgoplot, k_calerrors,delta_calerrors,/xlog,/ylog,psym=10,color=color_array[1],thick=thickness
    cgoplot, k_recover,delta_recover,/xlog,/ylog,psym=10,color=color_array[2],thick=thickness
    cgoplot, k_justEoR,delta_justEoR,/xlog,/ylog,psym=10,color=color_array[3],thick=thickness, linestyle=3
    cgoplot, k_centers,delta_eor,/xlog,/ylog,psym=10,color=color_array[4],thick=thickness
    
    cglegend,title=['input EoR','EoR in pipeline','recovered EoR','foregrounds, EoR','cal errors, foregrounds, EoR'],$
      color=[color_array[4],color_array[3],color_array[2],color_array[0],color_array[1]], location=[.57,1.05], $
      charsize=1.1, thick=thickness,linestyle=[0,3,0,0,0],length=.05
    
        if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
        stop

end