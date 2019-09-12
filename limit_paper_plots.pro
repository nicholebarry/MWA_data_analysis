pro limit_paper_plots

  ;sim_plots=1
  if keyword_set(sim_plots) then begin
    filename_24 = '/Users/nabarry/MWA/data/fhd_nb_sim_BH2_BH2_signalloss_exactforegrounds_6000/' ;left over foregounds, perfect
    filename_3 = '/Users/nabarry/MWA/data/fhd_nb_sim_BH2_BH2_justEoR/' ;just EoR


    obsname='zenith'
    ;obsname_1='1061319472'
    obsname_3='zenith'

    cube_type='res'
    power_24 = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_fullimg_'+cube_type+'_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_fullimg_'+cube_type+'_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_24 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    hubble_param=.71
    k_24=k_24/hubble_param
    weights = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_fullimg_'+cube_type+'_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave','weights')
    dsigma = (k_24^3.)/(2.*!pi^2.)/sqrt(weights)

    cube_type='res'
    power_3 = getvar_savefile(filename_3 + 'Combined_obs_'+obsname_3+'_cubeXX__even_odd_joint_fullimg_'+cube_type+'_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_3 + 'Combined_obs_'+obsname_3+'_cubeXX__even_odd_joint_fullimg_'+cube_type+'_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_3 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    k_3=k_3/hubble_param

    power_eor = getvar_savefile('/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave','power')
    k_centers_eor = getvar_savefile('/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave','k_centers')
    k_centers_eor=k_centers_eor/hubble_param

    delta=1
    if ~keyword_set(delta) then ytitle = 'P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)' else ytitle = '$\Delta$$\up2$ (mK$\up2$)'
    delta_24=power_24*(k_24^3.)/(2.*!pi^2.)
    delta_3=power_3*(k_3^3.)/(2.*!pi^2.)
    delta_eor=power_eor*(k_centers_eor^3.)/(2.*!pi^2.)

    cgPS_Open,'/Users/nabarry/MWA/data/EoR_recovered_limit_noise_lowchange.pdf',/quiet,/nomatch
    cgplot, k_24, (delta_24),/ylog,/xlog,psym=10,xrange=[.085,1.5],charsize=1.3, xtitle = 'k (!8h!x / Mpc)',yrange=[10^0., 10^6.],$
      ytitle=ytitle,title = 'In-situ simulation', color='forest green',thick=7, aspect=.5;thick=3
    ;cgoplot, k_24, dsigma, psym=10, linestyle =2, color='forest green', thick=7
    cgoplot, k_3, delta_3,/ylog,/xlog,psym=10, color='orange',thick=11, linestyle=0
    cgoplot, k_centers_eor, delta_eor,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='purple',thick=7
    cglegend, title=['input EoR', 'EoR in pipeline','foregrounds, EoR'],$
      color=['purple','orange','forest green'],location=[.45,.7], charsize=1.1,thick=7, linestyle=[0,0,0];,/background,/box
    ;cglegend, title=['input Gaussian EoR','recovered Gaussian EoR','no filter','Tukey filter','Blackman-Harris filter','Blackman-Harris filter, 100$\lambda$ extent'],$
    ;  color=['purple','dark green','navy','firebrick','orange','turquoise'],location=[.42,.85], charsize=1.1


    cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif

  ;cath_limit_compare=1
  if keyword_set(cath_limit_compare) then begin
    ;coarse_harm_width = 3
    ;wedge_angles=120d
    n_cubes=2

    cut = STRARR(n_cubes)
    cut[0]='cath_cubes'
    cut[1] = 'btl_noalltv_noocc4_zenith'

    chans = strarr(n_cubes)
    chans[0]='ch9-125_' & band='low'
    chans[1]='ch9-126_' & band='low'

    win = strarr(n_cubes)
    win[0]=''
    win[1]='fullimg_'

    spec_window='bh_'

    basefile = STRARR(n_cubes)
    ;basefile[0]=['/Users/nabarry/MWA/data/fhd_apb_EoR0_high_sem1_1/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[0]]
    basefile[1]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_lowchange/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans[1]]

    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[1]=[0,1,6,7,8,9,17,18,19,20,21,22,29,30,31,32,33,34,46]
    ;*ranges[1]=[7,8,9,10,11,21,22,23,24,25,26,36,37,38,39,40,51,52,53,54]+6

    *ranges[1]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[1]=[0]
    ;*ranges[3]=[6,7,8,9,10,11,21,22,23,24,25,26,36,37,38,39,40,41,50,51,52]

    cubes=['res']
    pols=['xx','yy']

    endfile = STRARR(n_cubes)
    endfile[1]='_averemove_sw'+spec_window+'dencorr_no_120deg_wedge_kperplambda18-80_kpar0.15-200_1dkpower.idlsave'
    ;endfile[1]='_averemove_sw'+spec_window+'dencorr_no_115deg_wedge_kperp0.025-200_kpar0.15-200_1dkpower.idlsave' ;TEMP

    limit_percent=0.9772 ;2 sigma
    hubble_param=0.71

    color_num = [10,12,14,16,18,20]
    ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101],[198,118,63],[197,120,62]]

    for n_i=0, n_cubes-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B,20B]

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 6 else thickness = 3

    for j=0,N_elements(pols)-1 do begin
      for cube_i=0,n_cubes-1 do begin

        if cube_i EQ 0 then begin
          textfast, cath_data, file_path='/Users/nabarry/MWA/data/cath_data_z7_errors_140obs.txt',first_line=1,/read
          ;textfast, cath_data, file_path='/Users/nabarry/MWA/data/cath_data_z7_lowkperp.txt',first_line=1,/read ;TEMP
          k = cath_data[0,*]
          if j EQ 0 then limits=(cath_data[2,*])^2.
          if j EQ 1 then limits=(cath_data[1,*])^2.
          if j EQ 0 then dsigma=(cath_data[4,*])^2.
          if j EQ 1 then dsigma=(cath_data[3,*])^2.
          ;stop
        endif
        if cube_i EQ 1 then begin
          file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
          restore,file
          n_k=n_elements(k_edges)
          k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
          ; find the best limit in this file
          delta=power*(k^3.)/(2.*!pi^2.)
          dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
          dsigma[0] = !Values.F_INFINITY
          if (cube_i EQ 1) then begin
            dsigma[0] = !Values.F_INFINITY
            dsigma[*ranges[cube_i]]= !Values.F_INFINITY
          endif
          inds = where(dsigma EQ 0,n_count)
          if n_count GT 0 then dsigma[inds]=!Values.F_INFINITY
          ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
          limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
          lim=min(limits,ind)
          header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans[cube_i]+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
          print,header

          k=k/hubble_param
        endif


        if (cube_i EQ 0) and keyword_set(pdf_plot) and (j EQ 0) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_lowchange/limit_Aug23_cath_lowchange_140obs_switchpol.pdf',/quiet,/nomatch

        position1=[0.1, 0.6, 0.5, 0.9]

        position4=[0.5, 0.6, 0.9, 0.9]

        if (cube_i EQ 0) and (j EQ 0) then begin
          position_use=position1
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 0) then position_use=position1

        if (cube_i EQ 0) and (j EQ 1) then begin
          position_use=position4
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use=''
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(A1)"
        endif
        if (cube_i EQ 1) and (j EQ 1) then position_use=position4

        color_use=color_array[cube_i]

        if ((cube_i EQ 0)) AND (j EQ 0) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.12,2], yrange=[10^3.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
          cglegend,title=['FHD/$\epsilon$ppsilon 2$\sigma$ upper limit', 'noise level'],$
            color=[color_array[1],color_array[1]], location=[.52,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1.1
        endif
        if ((cube_i EQ 0)) AND (j EQ 1) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.12,2], yrange=[10^3.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
          ;cglegend,title=['2016 2$\sigma$ upper limit','noise level','updated 2$\sigma upper limit','noise level'],$
          ;  color=[color_array[1],color_array[1],color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2,0,2],length=.03,vspace=1
          cglegend,title=['RTS/CHIPS 2$\sigma$ upper limit', 'noise level'],$
            color=[color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1.1
        endif

        if (cube_i NE 0) then $
          cgoplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.7], yrange=[10^2.,10^8.],color=color_use,thick=thickness,position=position_use
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use

      endfor

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif


  if keyword_set(noise_panel) then begin
    ps_wrapper,'/fred/oz048/MWA/CODE/FHD/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z','btl_noalltv_noocc4',full_image=1,freq_ch_range=[9,126],/pdf,plot_1d_delta=1,/exact_obsname,kperp_range_lambda_1dave=[15,80],kpar_range_1dave=[0.15,200],kperp_lambda_plot_range=[7,100],wedge_angles=103.7,set_data_ranges=0,sigma_range=[3*10^6.,10^9.],nev_range=[3*10^6.,10^9.],noise_range=[3*10^6.,10^9.]

    ps_wrapper,'/fred/oz048/MWA/CODE/FHD/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z','btl_noalltv_noocc4',freq_ch_range=[9,126],/pdf,plot_1d_delta=1,/exact_obsname,kperp_range_lambda_1dave=[18,80],kpar_range_1dave=[0.15,200],wedge_angles=103.7,set_data_ranges=0,sigma_range=[3*10^6.,10^9.],nev_range=[3*10^6.,10^9.],noise_range=[3*10^6.,10^9.],nnr_range=[.5,2.],ps_foldername='ps_true',full_image=1

  endif

  ;limit_plot=1
  if keyword_set(limit_plot) then begin


    cut = 'btl_noalltv_noocc4'
    ;cut = 'beardsley_thesis_list_zenith'
    cut2 = 'sub_cubes'
    chans='ch9-126_'
    win='fullimg_'
    spec_window='bh_'
    ;basefile='/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_fftchange_spectralwindow/Combined_obs_'+cut+'_cubeXX__even_odd_joint_'+win+chans
    ;basefile2='/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_fftchange_spectralwindow/Combined_obs_'+cut2+'_cubeXX__even_odd_joint_'+win+chans
    basefile='/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_lowchange/Combined_obs_'+cut+'_cubeXX__even_odd_joint_'+win+chans
    basefile2='/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_lowchange/Combined_obs_'+cut2+'_cubeXX__even_odd_joint_'+win+chans
    cubes=['res']
    pols=['xx','yy']
    averemove='_averemove_'
    kperp_end='80'
    kperp_start='18'
    dencorr = 'dencorr_'
    horizon='120'
    endfile=averemove+'sw'+spec_window+dencorr + 'no_'+horizon+'deg_wedge_kperplambda'+kperp_start+'-'+kperp_end+'_kpar0.15-200_1dkpower.idlsave'

    ;limit_percent=0.9772 ;2 sigma
    limit_percent = 0.97725
    hubble_param=0.71

    n_cubes=1
    color_num = [10,12,14,16,18]
    rgbcolors = [[137,117,202],[114,166,89],[73, 142, 217],[197,120,62],[165,42,42]]
    ;73, 142, 217 little brighter blue, 54,107,163 little darker blue
    ;rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101]];,[198,118,63]]

    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,28,38,39,40,41,42,43,54,53,55]+4

    *ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4

    for n_i=0, N_elements(color_num)-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B]

    ;pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 5 else thickness = 3

    for j=0,N_elements(pols)-1 do begin

      file=basefile+cubes+'_'+pols[j]+endfile
      ;file2=basefile2+cubes+'_'+pols[j]+endfile
      restore,file
      ;k_edges2 = getvar_savefile(file2, 'k_edges')
      ; find the best limit in this file
      n_k=n_elements(k_edges)
      ;n_k2=n_elements(k_edges2)
      k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
      ;k2=(k_edges2[1:(n_k2-1)]+k_edges2[0:(n_k2-2)])/2.
      ;power2 = getvar_savefile(file2,'power')
      delta=power*(k^3.)/(2.*!pi^2.)
      ;delta2=power2*(k2^3.)/(2.*!pi^2.)

      ;weights2 = getvar_savefile(file2,'weights')
      dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
      dsigma[0] = !Values.F_INFINITY
      ;ab_scaled=1
      if keyword_set(ab_scaled) then begin
        print, "Artificial scaling applied"
        dsigma = dsigma * 3.4
      endif
      ;dsigma2=(k2^3.)/(2.*!pi^2.)/sqrt(weights2)
      ;dsigma2[0] = !Values.F_INFINITY

      limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf((delta)/dsigma/sqrt(2)))*sqrt(2))+delta
      ;limits2=dsigma2*(inverf(limit_percent-(1.-limit_percent)*erf((delta2)/dsigma2/sqrt(2)))*sqrt(2))+delta2
      limits_abs=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(abs(delta)/dsigma/sqrt(2)))*sqrt(2))+abs(delta)
      lim=min(limits,ind)
      header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut+' '+chans+' '+cubes+' '+win+' '+pols[j]+')'
      print,header

      zeros = where(delta EQ 0,n_count)
      ;zeros2 = where(delta2 EQ 0,n_count)
      delta[zeros] = !Values.F_INFINITY
      delta[0] = !Values.F_INFINITY
      ;delta2[zeros2] = !Values.F_INFINITY
      ;delta2[0] = !Values.F_INFINITY
      delta[*ranges[0]] = !Values.F_INFINITY
      ;delta2[*ranges[0]] = !Values.F_INFINITY
      dsigma[0] = !Values.F_INFINITY
      dsigma[*ranges[0]] = !Values.F_INFINITY
      ;dsigma2[0] = !Values.F_INFINITY
      ;dsigma2[*ranges[0]] = !Values.F_INFINITY
      limits[0] = !Values.F_INFINITY
      limits[*ranges[0]] = !Values.F_INFINITY
      ;limits2[0] = !Values.F_INFINITY
      ;limits2[*ranges[0]] = !Values.F_INFINITY
      k=k/hubble_param
      ;k2=k2/hubble_param

      ;fiducial theory
      restore,'/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave'
      delta_eor=power*(k_centers^3.)/(2.*!pi^2.)
      k_centers=k_centers/hubble_param

      brad_update=1
      if keyword_set(brad_update) then begin
        ;GalaxyParam_LF+NF+Tau
        ;Fourier mode (k [Mpc^-1]), 68th, 95th and 99th percentiles on the power spectrum in units of mK^2
        brad_k = [0.02094395,0.03183814,0.0477084,0.06481693,0.0846983,0.1122589,0.1512655,0.2037167,0.2749778,0.3717481,0.5018322,0.6776517,0.9149306,1.235268,1.6676,2.207298,2.79783]
        brad_68 = [2.06836030818,9.33650238046,13.5354885655,27.9094808815,31.8379169811,35.6549002097,$
          37.7766299375,34.7228577191,22.5256589937,19.1373452916,19.0715957566,18.8232315706,18.1639861695,$
          18.6493671092,21.4087292929,28.7580668946,40.4441293602]
        brad_95 = [3.00633765723,12.8642551767,18.6682120582,38.322192403,43.5651253669,49.1075159329,$
          51.9213041607,48.790639303,31.3894142399,26.2910138229,26.8292876974,26.5396841581,26.0301170622,$
          27.1653881665,31.5002616162,42.7424864461,60.1632270707]
        brad_99 = [3.277767,13.62931,19.91,40.66932,46.2304,52.1699,55.14123,52.86583,34.2097,28.67557,29.69817,29.57222,$
          29.46968,30.93632,35.9361,48.86067,68.80418]

        ;k (Mpc^-1) and fiducial
        brad_k_2 = [2.094395e-02,3.183814e-02,4.770840e-02,6.481693e-02,8.469830e-02,1.122589e-01, 1.512655e-01,2.037167e-01,2.749778e-01,3.717481e-01,$
          5.018322e-01,6.776517e-01,9.149306e-01,1.235268e+00,1.667600e+00,2.207298e+00,2.797836e+00]
        brad_fiducial = [3.148450e-01,1.371862e+00,2.112182e+00,4.716543e+00,6.822401e+00,1.036046e+01,1.477729e+01,9.501991e+00,8.793401e+00,7.287165e+00,$
          6.430653e+00,5.709567e+00,5.390083e+00,5.408442e+00,6.280324e+00,8.959497e+00,1.353618e+01]

        ;  n_kb=n_elements(brad_k)
        ;brad_k=(brad_k[1:(n_kb-1)]+brad_k[0:(n_kb-2)])/2.

        brad_k_2 = brad_k_2/hubble_param
        brad_k = brad_k/hubble_param
      endif

      if (j EQ 0) and keyword_set(pdf_plot) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/limit_update_lowcheck_full.pdf',/quiet,/nomatch

      position1=[0.1, 0.6, 0.5, 0.9]

      position2=[0.1, 0.45, 0.5, 0.6]

      position4=[0.5, 0.6, 0.9, 0.9]

      position5=[0.5, 0.45, 0.9, 0.6]

      if (j EQ 0) then begin
        position_use=position1
        ;xtitle_use=''
        xtitle_use='k (h Mpc$\up-1$)'
        ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
        XTICKFORMAT_use="(I0)"
        YTICKFORMAT_use="(I0)"
      endif

      if (j EQ 1) then begin
        position_use=position4
        ;xtitle_use=''
        xtitle_use='k (h Mpc$\up-1$)'
        ytitle_use=''
        XTICKFORMAT_use="(I0)"
        YTICKFORMAT_use="(A1)"
      endif

      color_use=color_array[2]
      low_range=10^0.
      low_x_range=.14

      if (j EQ 0) then begin
        cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^0.,10^8.],ytitle=ytitle_use, $
          xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,/noerase
        XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
        cglegend,title=['fiducial theory','95% confidence'],$
          color=[color_array[3],color_array[3]], location=[.52,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1
      endif
      if (j EQ 1) then begin
        cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^0.,10^8.],ytitle=ytitle_use, $
          xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
        XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
        cglegend,title=['measured power','2$\sigma$ upper limit','noise level'],$
          color=['black',color_array[2],color_array[2]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,0,2],length=.03,vspace=1.1
      endif
      ;cgoplot, k2,abs(limits2),/xlog,/ylog,psym=10,color=color_array[1],thick=4,position=position_use, linestyle=3

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

      cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=color_use,thick=thickness,position=position_use
      cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use
      cgoplot, k,delta,/xlog,/ylog,color='black',psym=10,thick=thickness-1,position=position_use
      if ~keyword_set(brad_update) then cgoplot, k_centers,delta_eor,/xlog,/ylog,color=color_array[3],thick=thickness-1,position=position_use else begin
        cgoplot, brad_k_2,brad_fiducial,/xlog,/ylog,color=color_array[3],thick=thickness-1,position=position_use
        cgoplot, brad_k,brad_95,/xlog,/ylog,color=color_array[3],thick=thickness-1,position=position_use,linestyle=2;,psym=10

      endelse



      if (j EQ 0) then begin

        ;        cgplot, k2,abs(limits2),/xlog,/ylog,psym=10, xrange=[.14,1.9], c,ytitle=ytitle_use, YMINOR=1, $
        ;          xtitle='k (h Mpc$\up-1$)', charsize =.8, color=color_array[1],thick=thickness,position=position2,XTICKFORMAT=XTICKFORMAT_use,/noerase
        ;          cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=color_use,thick=thickness,position=position2
        ;          cgoplot, k2,dsigma2,/xlog,/ylog,psym=10,linestyle=2, color=color_array[1],thick=thickness,position=position2
        ;          cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position2

        ;        cglegend,title=['noise level'],$
        ;          color=[color_array[1]], location=[.52,.58], charsize=.7, thick=thickness,linestyle=[2],length=.03,vspace=1
      endif
      if (j EQ 1) then begin
        ;        cgplot, k2,abs(limits2),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, YMINOR=1, $
        ;          xtitle='k (h Mpc$\up-1$)', charsize =.8, color=color_array[1],thick=thickness,position=position5,XTICKFORMAT=XTICKFORMAT_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
        ;          cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=color_use,thick=thickness,position=position5
        ;          cgoplot, k2,dsigma2,/xlog,/ylog,psym=10,linestyle=2, color=color_array[1],thick=thickness,position=position5
        ;          cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position5

        ;        cglegend,title=['pre-SSINS 2$\sigma$ upper limit'],$
        ;          color=[color_array[1]], location=[.12,.58], charsize=.7, thick=thickness,linestyle=[0],length=.03,vspace=1.1
      endif

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage

    ;full doubles, no LS, 18-80, 0.15-200, 120deg
    ;#Limit: 9373.3024 mK^2, at k = 0.20330634 h Mpc^-1 (btl_noalltv_noocc4 ch9-126_ res fullimg_ xx)
    ;#Limit: 3273.5603 mK^2, at k = 0.20330634 h Mpc^-1 (btl_noalltv_noocc4 ch9-126_ res fullimg_ yy)

    ;with LS
    ;#Limit: 10282.498 mK^2, at k = 0.20330593 h Mpc^-1 (btl_noalltv_noocc4 ch9-126_ res fullimg_ xx)
    ;#Limit: 3890.4814 mK^2, at k = 0.20330593 h Mpc^-1 (btl_noalltv_noocc4 ch9-126_ res fullimg_ yy)

    ;pointing looks
    ;#Limit: 37670.519 mK^2, at k = 0.20330634 h Mpc^-1 (btl_noalltv_noocc4_zenith ch9-126_ res fullimg_ xx)
    ;#Limit: 3319.2405 mK^2, at k = 0.2323501 h Mpc^-1 (btl_noalltv_noocc4_zenith ch9-126_ res fullimg_ yy)
    ;#Limit: 29834.284 mK^2, at k = 0.20330634 h Mpc^-1 (beardsley_thesis_list_zenith ch9-126_ res fullimg_ xx)
    ;#Limit: 14089.952 mK^2, at k = 0.20330634 h Mpc^-1 (beardsley_thesis_list_zenith ch9-126_ res fullimg_ yy)
    ;xx: .79 yy: 4.24
    ;;
    ;#Limit: 34977.434 mK^2, at k = 0.20330593 h Mpc^-1 (btl_noalltv_noocc4_plusone ch9-126_ res fullimg_ xx)
    ;#Limit: 17538.733 mK^2, at k = 0.20330593 h Mpc^-1 (btl_noalltv_noocc4_plusone ch9-126_ res fullimg_ yy)
    ;#Limit: 29759.657 mK^2, at k = 0.20330593 h Mpc^-1 (beardsley_thesis_list_plusone ch9-126_ res fullimg_ xx)
    ;#Limit: 20092.925 mK^2, at k = 0.20330593 h Mpc^-1 (beardsley_thesis_list_plusone ch9-126_ res fullimg_ yy)
    ;
    ;#Limit: 52942.289 mK^2, at k = 0.20330593 h Mpc^-1 (btl_noalltv_noocc4_plustwo ch9-126_ res fullimg_ xx)
    ;#Limit: 19037.562 mK^2, at k = 0.20330593 h Mpc^-1 (btl_noalltv_noocc4_plustwo ch9-126_ res fullimg_ yy)
    ;#Limit: 42706.043 mK^2, at k = 0.17426223 h Mpc^-1 (beardsley_thesis_list_plustwo ch9-126_ res fullimg_ xx)
    ;#Limit: 24459.931 mK^2, at k = 0.20330593 h Mpc^-1 (beardsley_thesis_list_plustwo ch9-126_ res fullimg_ yy)
    ;
    ;#Limit: 27264.747 mK^2, at k = 0.20330593 h Mpc^-1 (beardsley_thesis_list_minusone ch9-126_ res fullimg_ xx)
    ;#Limit: 22136.987 mK^2, at k = 0.20330593 h Mpc^-1 (beardsley_thesis_list_minusone ch9-126_ res fullimg_ yy)
    ;#Limit: 23460.126 mK^2, at k = 0.20330593 h Mpc^-1 (btl_noalltv_noocc4_minusone ch9-126_ res fullimg_ xx)
    ;#Limit: -202607.45 mK^2, at k = 1.4231415 h Mpc^-1 (btl_noalltv_noocc4_minusone ch9-126_ res fullimg_ yy)

    ;#Limit: 27994.08 mK^2, at k = 0.23234963 h Mpc^-1 (beardsley_thesis_list_minustwo ch9-126_ res fullimg_ xx)
    ;#Limit: 19262.499 mK^2, at k = 0.23234963 h Mpc^-1 (beardsley_thesis_list_minustwo ch9-126_ res fullimg_ yy)

    stop
  endif


  ;  v2_limit_plot=1
  if keyword_set(v2_limit_plot) then begin

    cut = 'btl_noalltv_noocc4'
    cut2 = 'sub_cubes'
    chans='ch9-126_'
    win='fullimg_'
    spec_window='bh_'
    basefile='/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_lowchange/Combined_obs_'+cut+'_cubeXX__even_odd_joint_'+win+chans
    basefile2='/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_lowchange/Combined_obs_'+cut2+'_cubeXX__even_odd_joint_'+win+chans
    cubes=['res']
    pols=['xx','yy']
    averemove='averemove_'
    endfile='_'+averemove+'sw'+spec_window+'dencorr_no_110deg_wedge_kperplambda15-80_kpar0.15-200_1dkpower.idlsave'

    ;limit_percent=0.9772 ;2 sigma
    limit_percent = 0.97725
    hubble_param=0.71

    n_cubes=1
    color_num = [10,12,14,16]
    rgbcolors = [[137,117,202],[114,166,89],[54,107,163],[197,120,62]]
    ;rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101]];,[198,118,63]]

    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,28,38,39,40,41,42,43,54,53,55]+4

    *ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4

    for n_i=0, N_elements(color_num)-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B]

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 5 else thickness = 3

    for j=0,N_elements(pols)-1 do begin

      file=basefile+cubes+'_'+pols[j]+endfile
      file2=basefile2+cubes+'_'+pols[j]+endfile
      restore,file
      k_edges2 = getvar_savefile(file2, 'k_edges')
      ; find the best limit in this file
      n_k=n_elements(k_edges)
      n_k2=n_elements(k_edges2)
      k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
      k2=(k_edges2[1:(n_k2-1)]+k_edges2[0:(n_k2-2)])/2.
      power2 = getvar_savefile(file2,'power')
      delta=power*(k^3.)/(2.*!pi^2.)
      delta2=power2*(k2^3.)/(2.*!pi^2.)

      weights2 = getvar_savefile(file2,'weights')
      dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
      dsigma[0] = !Values.F_INFINITY
      dsigma2=(k2^3.)/(2.*!pi^2.)/sqrt(weights2)
      dsigma2[0] = !Values.F_INFINITY

      limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf((delta)/dsigma/sqrt(2)))*sqrt(2))+delta
      limits2=dsigma2*(inverf(limit_percent-(1.-limit_percent)*erf((delta2)/dsigma2/sqrt(2)))*sqrt(2))+delta2
      limits_abs=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(abs(delta)/dsigma/sqrt(2)))*sqrt(2))+abs(delta)
      lim=min(limits,ind)
      header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut+' '+chans+' '+cubes+' '+win+' '+pols[j]+')'
      print,header

      zeros = where(delta EQ 0,n_count)
      zeros2 = where(delta2 EQ 0,n_count)
      delta[zeros] = !Values.F_INFINITY
      delta[0] = !Values.F_INFINITY
      delta2[zeros2] = !Values.F_INFINITY
      delta2[0] = !Values.F_INFINITY
      delta[*ranges[0]] = !Values.F_INFINITY
      delta2[*ranges[0]] = !Values.F_INFINITY
      dsigma[0] = !Values.F_INFINITY
      dsigma[*ranges[0]] = !Values.F_INFINITY
      dsigma2[0] = !Values.F_INFINITY
      dsigma2[*ranges[0]] = !Values.F_INFINITY
      limits[0] = !Values.F_INFINITY
      limits[*ranges[0]] = !Values.F_INFINITY
      limits2[0] = !Values.F_INFINITY
      limits2[*ranges[0]] = !Values.F_INFINITY
      k=k/hubble_param
      k2=k2/hubble_param

      ;fiducial theory
      restore,'/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave'
      delta_eor=power*(k_centers^3.)/(2.*!pi^2.)
      k_centers=k_centers/hubble_param

      if (j EQ 0) and keyword_set(pdf_plot) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/round_fix/limit_update_9-126_blue_comp.pdf',/quiet,/nomatch

      position1=[0.1, 0.6, 0.5, 0.9]

      position2=[0.1, 0.52, 0.5, 0.6]
      position2b=[0.1, 0.44, 0.5, 0.52]

      position4=[0.5, 0.6, 0.9, 0.9]

      position5=[0.5, 0.52, 0.9, 0.6]
      position5b=[0.5, 0.44, 0.9, 0.52]

      if (j EQ 0) then begin
        position_use=position1
        xtitle_use=''
        ;xtitle_use='k (h Mpc$\up-1$)'
        ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
        XTICKFORMAT_use="(I0)"
        YTICKFORMAT_use="(I0)"
      endif

      if (j EQ 1) then begin
        position_use=position4
        xtitle_use=''
        ;xtitle_use='k (h Mpc$\up-1$)'
        ytitle_use=''
        XTICKFORMAT_use="(I0)"
        YTICKFORMAT_use="(A1)"
      endif

      color_use=color_array[2]
      low_range=10^1.
      low_x_range=.14

      if (j EQ 0) then begin
        cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^1.,10^8.],ytitle=ytitle_use, $
          xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT="(A1)",/noerase
        XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
        cglegend,title=['2$\sigma$ upper limit','noise level'],$
          color=[color_array[2],color_array[2]], location=[.52,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1
      endif
      if (j EQ 1) then begin
        cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^1.,10^8.],ytitle=ytitle_use, $
          xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT="(A1)",YTICKFORMAT=YTICKFORMAT_use,/noerase
        XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
        cglegend,title=['measured power','fiducial theory'],$
          color=['black',color_array[3]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,0],length=.03,vspace=1.1
      endif
      ;cgoplot, k2,abs(limits2),/xlog,/ylog,psym=10,color=color_array[1],thick=4,position=position_use, linestyle=3

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

      cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=color_use,thick=thickness,position=position_use
      cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use
      cgoplot, k,delta,/xlog,/ylog,color='black',psym=10,thick=thickness-1,position=position_use
      cgoplot, k_centers,delta_eor,/xlog,/ylog,color=color_array[3],thick=thickness-1,position=position_use



      if (j EQ 0) then begin

        plotlimits=limits-limits2
        zeros=where(plotlimits LT 0,n_count)
        if n_count GT 0 then plotlimits[zeros]=10^1.
        cgplot, k2,plotlimits,/xlog,/ylog,psym=10, yrange=[10^2.,10^6.],xrange=[.14,1.9],ytitle=ytitle_use, YMINOR=1, xstyle=4,$
          xtitle='k (h Mpc$\up-1$)', charsize =.8, color=color_array[1],thick=thickness,position=position2,XTICKFORMAT=XTICKFORMAT_use,/noerase

        ;  cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=color_use,thick=thickness,position=position2
        cgoplot, k2,dsigma2,/xlog,/ylog,psym=10,linestyle=2, color=color_array[1],thick=thickness,position=position2
        ; cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position2

        cglegend,title=['noise level'],$
          color=[color_array[1]], location=[.52,.58], charsize=.7, thick=thickness,linestyle=[2],length=.03,vspace=1

        plotlimits=limits2-limits
        zeros=where(plotlimits LT 0,n_count)
        if n_count GT 0 then plotlimits[zeros]=10^1.
        cgplot, k2,plotlimits,/xlog,/ylog,psym=10, yrange=[10^6.,10^2.],xrange=[.14,1.9],ytitle=ytitle_use, YMINOR=1,xstyle=8,$
          xtitle='k (h Mpc$\up-1$)', charsize =.8, color=color_array[1],thick=thickness,position=position2b,XTICKFORMAT=XTICKFORMAT_use,/noerase
        cgoplot, k2,dsigma2,/xlog,/ylog,psym=10,linestyle=2, color=color_array[1],thick=thickness,position=position2b
      endif
      if (j EQ 1) then begin
        plotlimits=limits-limits2
        zeros=where(plotlimits LT 0,n_count)
        ;if n_count GT 0 then plotlimits[zeros]=!Values.F_INFINITY
        if n_count GT 0 then plotlimits[zeros]=10^1.
        cgplot, k2,plotlimits,/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^6.],ytitle=ytitle_use, YMINOR=1, xstyle=4, $
          xtitle='k (h Mpc$\up-1$)', charsize =.8, color=color_array[1],thick=thickness,position=position5,XTICKFORMAT=XTICKFORMAT_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
        ; cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=color_use,thick=thickness,position=position5
        cgoplot, k2,dsigma2,/xlog,/ylog,psym=10,linestyle=2, color=color_array[1],thick=thickness,position=position5
        ; cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position5

        cglegend,title=['Limit change with SSINS'],$;title=['2$\sigma$ upper limit difference pre- and post-SSINS'],$
          color=[color_array[1]], location=[.12,.58], charsize=.7, thick=thickness,linestyle=[0],length=.03,vspace=1.1

        plotlimits=limits2-limits
        ;                zeros=where(plotlimits LT 0,n_count)
        ;                if n_count GT 0 then plotlimits[zeros]=!Values.F_INFINITY
        zeros=where(plotlimits LT 0,n_count)
        if n_count GT 0 then plotlimits[zeros]=10^1.
        cgplot, k2,plotlimits,/xlog,/ylog,psym=10, yrange=[10^6.,10^2.],xrange=[.14,1.9],ytitle=ytitle_use, YMINOR=1,xstyle=8,$
          xtitle='k (h Mpc$\up-1$)', charsize =.8, color=color_array[1],thick=thickness,position=position5b,YTICKFORMAT="(A1)",/noerase
        cgoplot, k2,dsigma2,/xlog,/ylog,psym=10,linestyle=2, color=color_array[1],thick=thickness,position=position5b
      endif

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage


    stop
  endif


  ;ab_limit_compare_1band=1
  if keyword_set(ab_limit_compare_1band) then begin
    ;coarse_harm_width = 3
    ;wedge_angles=120d
    n_cubes=2

    cut = STRARR(n_cubes)
    cut[0]='wedge_cut_plus_res_cut'
    cut[1] = 'sub_cubes'
    ;cut[0]='wedge_cut_plus_res_cut_p0'
    ;cut[1] = 'beardsley_thesis_list_zenith'

    chans = strarr(n_cubes)
    chans[0]='ch9-126_' & band='low'
    chans[1]='ch9-126_' & band='low'

    win = strarr(n_cubes)
    win[0]=''
    win[1]='fullimg_'
    ;win[1]=''

    spec_window='bh_'

    basefile = STRARR(n_cubes)
    basefile[0]=['/Users/nabarry/MWA/data/fhd_apb_EoR0_high_sem1_1/ps_checkout/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[0]]
    basefile[1]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_lowchange/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans[1]]

    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[1]=[0,1,6,7,8,9,17,18,19,20,21,22,29,30,31,32,33,34,46]
    ;*ranges[1]=[7,8,9,10,11,22,23,24,25,26,37,38,39,40,41,52,53,54]+6 ;cbw3
    *ranges[1]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    *ranges[0]=[51]
    ;*ranges[3]=[6,7,8,9,10,11,21,22,23,24,25,26,36,37,38,39,40,41,50,51,52]

    cubes=['res']
    pols=['xx','yy']

    endfile = STRARR(n_cubes)
    ;endfile[0]='_swbh_dencorr_no_120deg_wedge_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'
    endfile[0]='_bh_dencorr_no_120deg_wedge_cbw3_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'
    endfile[1]='_averemove_sw'+spec_window+'dencorr_no_120deg_wedge_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'

    limit_percent=0.9772 ;2 sigma
    hubble_param=0.71

    color_num = [10,12,14,16,18,20]
    ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101],[198,118,63],[197,120,62]]

    for n_i=0, n_cubes-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B,20B]

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 6 else thickness = 3
    ;scale_factor1 = 1.33069
    scale_factor1 = 1.28213
    scale_factor2 = 1.22765

    for j=0,N_elements(pols)-1 do begin
      for cube_i=0,n_cubes-1 do begin

        file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
        restore,file
        ; find the best limit in this file
        n_k=n_elements(k_edges)
        k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
        delta=power*(k^3.)/(2.*!pi^2.)
        dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
        if cube_i NE 0 then dsigma[0] = !Values.F_INFINITY
        ;if (cube_i EQ 1) then begin
        ;dsigma[0] = !Values.F_INFINITY
        dsigma[*ranges[cube_i]]= !Values.F_INFINITY
        ;endif
        inds = where(dsigma EQ 0,n_count)
        if n_count GT 0 then dsigma[inds]=!Values.F_INFINITY
        ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
        limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta

        if keyword_set(scale_factor1) and cube_i EQ 0 then begin
          if j EQ 0 then scale_factor = scale_factor1 else scale_factor = scale_factor2
          limits = limits / scale_factor
          dsigma = dsigma / scale_factor
          delta= delta / scale_factor
        endif

        lim=min(limits,ind)
        header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans[cube_i]+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
        print,header

        k=k/hubble_param

        if (cube_i EQ 0) and keyword_set(pdf_plot) and (j EQ 0) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_lowchange/limit_AB_compare_9-126_lowchange.pdf',/quiet,/nomatch

        position1=[0.1, 0.6, 0.5, 0.9]

        position4=[0.5, 0.6, 0.9, 0.9]

        if (cube_i EQ 0) and (j EQ 0) then begin
          position_use=position1
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 0) then position_use=position1

        if (cube_i EQ 0) and (j EQ 1) then begin
          position_use=position4
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use=''
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(A1)"
        endif
        if (cube_i EQ 1) and (j EQ 1) then position_use=position4

        color_use=color_array[cube_i]


        if ((cube_i EQ 0)) AND (j EQ 0) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8


          cglegend,title=['2016 2$\sigma$ upper limit','noise level'],$
            color=[color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1.1
        endif
        if ((cube_i EQ 0)) AND (j EQ 1) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
          ;cglegend,title=['2016 2$\sigma$ upper limit','noise level','updated 2$\sigma upper limit','noise level'],$
          ;  color=[color_array[1],color_array[1],color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2,0,2],length=.03,vspace=1
          cglegend,title=['updated 2$\sigma$ upper limit','noise level'],$
            color=[color_array[1],color_array[1]], location=[.52,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1
        endif

        if (cube_i NE 0) then $
          cgoplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.7], yrange=[10^2.,10^8.],color=color_use,thick=thickness,position=position_use
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use

      endfor

      ;#Limit: 22408.321 mK^2, at k = 0.23044845 h Mpc^-1 (wedge_cut_plus_res_cut ch9-127_ res  xx)
      ;#Limit: 9223.0543 mK^2, at k = 0.23044845 h Mpc^-1 (sub_cubes ch9-127_ res fullimg_ xx)
      ;2.42960 better
      ;#Limit: 27702.209 mK^2, at k = 0.23044845 h Mpc^-1 (wedge_cut_plus_res_cut ch9-127_ res  yy)
      ;#Limit: 7904.3834 mK^2, at k = 0.23044845 h Mpc^-1 (sub_cubes ch9-127_ res fullimg_ yy)
      ;3.50466 better

      ;update
      ;#Limit: 24231.534 mK^2, at k = 0.2342833 h Mpc^-1 (wedge_cut_plus_res_cut ch9-125_ res  xx)
      ;#Limit: 9774.8649 mK^2, at k = 0.2342833 h Mpc^-1 (sub_cubes ch9-125_ res fullimg_ xx)
      ;2.47896 better
      ;#Limit: 28403.79 mK^2, at k = 0.2342833 h Mpc^-1 (wedge_cut_plus_res_cut ch9-125_ res  yy)
      ;#Limit: 8086.5888 mK^2, at k = 0.2342833 h Mpc^-1 (sub_cubes ch9-125_ res fullimg_ yy)
      ;3.51246 better

      ;final update
      ;#Limit: 16695.796 mK^2, at k = 0.20330593 h Mpc^-1 (wedge_cut_plus_res_cut ch9-126_ res  xx)
      ;#Limit: 8859.4568 mK^2, at k = 0.20330593 h Mpc^-1 (sub_cubes ch9-126_ res fullimg_ xx)
      ;;1.8845 better
      ;#Limit: 19336.688 mK^2, at k = 0.20330593 h Mpc^-1 (wedge_cut_plus_res_cut ch9-126_ res  yy)
      ;#Limit: 8676.4582 mK^2, at k = 0.26139334 h Mpc^-1 (sub_cubes ch9-126_ res fullimg_ yy)
      ;;2.2286

      ;final final update (with LS, and checking out eppsilon for ab)
      ;#Limit: 20418.655 mK^2, at k = 0.21817447 h Mpc^-1 (wedge_cut_plus_res_cut ch9-126_ res  xx)
      ;#Limit: 9675.958 mK^2, at k = 0.20330634 h Mpc^-1 (sub_cubes ch9-126_ res fullimg_ xx)
      ;2.11025 better
      ;#Limit: 23749.79 mK^2, at k = 0.24721821 h Mpc^-1 (wedge_cut_plus_res_cut ch9-126_ res  yy)
      ;#Limit: 8585.4886 mK^2, at k = 0.26139386 h Mpc^-1 (sub_cubes ch9-126_ res fullimg_ yy)
      ;2.7662711 better

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif


  ;ab_repro=1
  if keyword_set(ab_repro) then begin

    ;scale_factor1 = 1.28213
    ;scale_factor2 = 1.22765

    n_cubes=1

    cut = STRARR(n_cubes)
    cut[0]='wedge_cut_plus_res_cut'

    chans = strarr(n_cubes)
    chans[0]='ch9-126_' & band='low'

    win = strarr(n_cubes)
    win[0]=''

    spec_window='bh_'

    basefile = STRARR(n_cubes)
    basefile[0]=['/Users/nabarry/MWA/data/fhd_apb_EoR0_high_sem1_1/ps_checkout/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[0]]

    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[0]=[0,1,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+8
    *ranges[0]=[51]

    cubes=['res']
    pols=['xx','yy']

    endfile = STRARR(n_cubes)
    endfile[0]='_bh_dencorr_no_120deg_wedge_cbw3_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'

    limit_percent=0.9772 ;2 sigma
    hubble_param=0.71


    color_num = [10,12,14,16,18,20]
    ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[[197,120,62]],[198,118,63],[197,120,62]]

    for n_i=0, 4 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B,20B]

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 6 else thickness = 3

    for j=0,N_elements(pols)-1 do begin
      for cube_i=0,n_cubes-1 do begin

        file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
        restore,file
        ; find the best limit in this file
        n_k=n_elements(k_edges)
        k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
        delta=power*(k^3.)/(2.*!pi^2.)
        dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
        dsigma[0] = !Values.F_INFINITY
        ;if (cube_i EQ 1) then begin
        dsigma[0] = !Values.F_INFINITY
        dsigma[*ranges[cube_i]]= !Values.F_INFINITY

        ;endif
        inds = where(dsigma EQ 0,n_count)
        if n_count GT 0 then dsigma[inds]=!Values.F_INFINITY
        ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
        limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
        delta[*ranges[cube_i]]= !Values.F_INFINITY
        delta[0] = !Values.F_INFINITY

        zeros = where(delta EQ 0,n_count)
        ;zeros2 = where(delta2 EQ 0,n_count)
        delta[zeros] = !Values.F_INFINITY
        delta[0] = !Values.F_INFINITY
        ;delta2[zeros2] = !Values.F_INFINITY
        ;delta2[0] = !Values.F_INFINITY
        delta[*ranges[0]] = !Values.F_INFINITY
        ;delta2[*ranges[0]] = !Values.F_INFINITY
        dsigma[0] = !Values.F_INFINITY
        dsigma[*ranges[0]] = !Values.F_INFINITY
        ;dsigma2[0] = !Values.F_INFINITY
        ;dsigma2[*ranges[0]] = !Values.F_INFINITY
        limits[0] = !Values.F_INFINITY
        limits[*ranges[0]] = !Values.F_INFINITY

        if keyword_set(scale_factor1) and cube_i EQ 0 then begin
          if j EQ 0 then scale_factor = scale_factor1 else scale_factor = scale_factor2
          limits = limits / scale_factor
          dsigma = dsigma / scale_factor
          delta= delta / scale_factor
        endif

        ;fiducial theory
        restore,'/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave'
        delta_eor=power*(k_centers^3.)/(2.*!pi^2.)
        k_centers=k_centers/hubble_param

        brad_update=1
        if keyword_set(brad_update) then begin
          ;GalaxyParam_LF+NF+Tau
          ;Fourier mode (k [Mpc^-1]), 68th, 95th and 99th percentiles on the power spectrum in units of mK^2
          brad_k = [0.02094395,0.03183814,0.0477084,0.06481693,0.0846983,0.1122589,0.1512655,0.2037167,0.2749778,0.3717481,0.5018322,0.6776517,0.9149306,1.235268,1.6676,2.207298,2.79783]
          brad_68 = [2.06836030818,9.33650238046,13.5354885655,27.9094808815,31.8379169811,35.6549002097,$
            37.7766299375,34.7228577191,22.5256589937,19.1373452916,19.0715957566,18.8232315706,18.1639861695,$
            18.6493671092,21.4087292929,28.7580668946,40.4441293602]
          brad_95 = [3.00633765723,12.8642551767,18.6682120582,38.322192403,43.5651253669,49.1075159329,$
            51.9213041607,48.790639303,31.3894142399,26.2910138229,26.8292876974,26.5396841581,26.0301170622,$
            27.1653881665,31.5002616162,42.7424864461,60.1632270707]
          brad_99 = [3.277767,13.62931,19.91,40.66932,46.2304,52.1699,55.14123,52.86583,34.2097,28.67557,29.69817,29.57222,$
            29.46968,30.93632,35.9361,48.86067,68.80418]

          ;k (Mpc^-1) and fiducial
          brad_k_2 = [2.094395e-02,3.183814e-02,4.770840e-02,6.481693e-02,8.469830e-02,1.122589e-01, 1.512655e-01,2.037167e-01,2.749778e-01,3.717481e-01,$
            5.018322e-01,6.776517e-01,9.149306e-01,1.235268e+00,1.667600e+00,2.207298e+00,2.797836e+00]
          brad_fiducial = [3.148450e-01,1.371862e+00,2.112182e+00,4.716543e+00,6.822401e+00,1.036046e+01,1.477729e+01,9.501991e+00,8.793401e+00,7.287165e+00,$
            6.430653e+00,5.709567e+00,5.390083e+00,5.408442e+00,6.280324e+00,8.959497e+00,1.353618e+01]

          ;  n_kb=n_elements(brad_k)
          ;brad_k=(brad_k[1:(n_kb-1)]+brad_k[0:(n_kb-2)])/2.

          brad_k_2 = brad_k_2/hubble_param
          brad_k = brad_k/hubble_param
        endif

        lim=min(limits,ind)
        header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans[cube_i]+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
        print,header

        k=k/hubble_param

        if (cube_i EQ 0) and keyword_set(pdf_plot) and (j EQ 0) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_apb_EoR0_high_sem1_1/limit_AB_reproduction_checkout_notscaled.pdf',/quiet,/nomatch

        ;position1=[0.1, 0.6, 0.5, 0.9]
        position1=[0.1,  0.6, 0.5,   0.9]

        ;position4=[0.5, 0.6, 0.9, 0.9]
        position4=[0.1,   0.3, 0.5,  0.6 ]

        if (cube_i EQ 0) and (j EQ 0) then begin
          position_use=position1
          ;xtitle_use='k (h Mpc$\up-1$)'
          xtitle_use=''
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          ;XTICKFORMAT_use="(I0)"
          XTICKFORMAT_use="(A1)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 0) then position_use=position1

        if (cube_i EQ 0) and (j EQ 1) then begin
          position_use=position4
          xtitle_use='k (h Mpc$\up-1$)'
          ;ytitle_use=''
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          XTICKFORMAT_use="(I0)"
          ;YTICKFORMAT_use="(A1)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 1) then position_use=position4

        color_use=color_array[cube_i]


        if ((cube_i EQ 0)) AND (j EQ 0) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[3.*10^0.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.35, position_use[3]-.027, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
          cglegend,title=['measured power','fiducial theory', '95% confidence'],$
            color=['black',color_array[3],color_array[3]], location=[.12,.58], charsize=.7, thick=thickness,linestyle=[0,0,2],length=.03,vspace=1
        endif
        if ((cube_i EQ 0)) AND (j EQ 1) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[3.*10^0.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.35, position_use[3]-.027, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
          ;cglegend,title=['2016 2$\sigma$ upper limit','noise level','updated 2$\sigma upper limit','noise level'],$
          ;  color=[color_array[1],color_array[1],color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2,0,2],length=.03,vspace=1
          cglegend,title=['2016 2$\sigma$ upper limit','noise level'],$
            color=[color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1.1
        endif




        low_range=3.*10^0.
        low_x_range=.14

        error_low = (dsigma*2.); < (delta);*.99999
        lows = where(delta LT 0, n_count)
        if n_count GT 0 then error_low[lows]=abs(delta[lows])
        lows = where((abs(delta) - error_low) LT low_range,n_count)
        if n_count GT 0 then error_low[lows]=abs(delta[lows])-low_range
        error_high = dsigma*2.
        lows = where(k LT low_x_range,n_count)
        if n_count GT 0 then error_high[lows]=!Values.F_INFINITY
        if n_count GT 0 then limits[lows]=!Values.F_INFINITY
        delta_k=(k[1]-k[0])/2.
        for k_i=0, n_k-2 do $
          cgColorFill, [k[k_i]-delta_k, k[k_i]+delta_k, k[k_i]+delta_k,k[k_i]-delta_k], $
          [limits[k_i], limits[k_i], abs(delta[k_i])-error_low[k_i],abs(delta[k_i])-error_low[k_i]], $
          Color='grey',/checkforfinite

        cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=color_use,thick=thickness,position=position_use
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use
        cgoplot, k,delta,/xlog,/ylog,color='black',psym=10,thick=thickness-1,position=position_use
        ;cgoplot, k_centers,delta_eor,/xlog,/ylog,color=color_array[3],thick=thickness-1,position=position_use
        cgoplot, brad_k_2,brad_fiducial,/xlog,/ylog,color=color_array[3],thick=thickness-1,position=position_use
        cgoplot, brad_k,brad_95,/xlog,/ylog,color=color_array[3],thick=thickness-1,position=position_use,linestyle=2;,psym=10

        cgoplot,k,abs(limits),/xlog,/ylog,psym=10,color=color_use
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use

      endfor

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif


  ;ssins_limit_compare_1band=1
  if keyword_set(ssins_limit_compare_1band) then begin
    ;coarse_harm_width = 3
    ;wedge_angles=120d
    n_cubes=2

    cut = STRARR(n_cubes)
    ;cut[0]='sub_cubes'
    cut[0]='beardsley_thesis_list_zenith'
    cut[1] = 'btl_noalltv_noocc4_zenith'

    chans = strarr(n_cubes)
    chans[0]='ch9-126_' & band='low'
    chans[1]='ch9-126_' & band='low'

    win = strarr(n_cubes)
    win[0]='fullimg_'
    win[1]='fullimg_'
    ;win[1]=''

    spec_window='bh_'

    basefile = STRARR(n_cubes)
    basefile[0]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_lowchange/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[0]]
    basefile[1]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_lowchange/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans[1]]

    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[1]=[0,1,6,7,8,9,17,18,19,20,21,22,29,30,31,32,33,34,46]
    ;*ranges[1]=[7,8,9,10,11,22,23,24,25,26,37,38,39,40,41,52,53,54]+6 ;cbw3
    *ranges[1]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    *ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[1]=[0]
    ;*ranges[3]=[6,7,8,9,10,11,21,22,23,24,25,26,36,37,38,39,40,41,50,51,52]

    cubes=['res']
    pols=['xx','yy']

    endfile = STRARR(n_cubes)
    endfile[0]='_averemove_sw'+spec_window+'dencorr_no_120deg_wedge_kperplambda18-80_kpar0.15-200_1dkpower.idlsave'
    endfile[1]='_averemove_sw'+spec_window+'dencorr_no_120deg_wedge_kperplambda18-80_kpar0.15-200_1dkpower.idlsave'

    limit_percent=0.9772 ;2 sigma
    hubble_param=0.71

    color_num = [10,12,14,16,18,20]
    ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    ;73, 142, 217 lighter blue
    rgbcolors=[[114,166,89],[73, 142, 217],[202,86,139],[91,169,101],[198,118,63],[197,120,62]]

    for n_i=0, n_cubes-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B,20B]

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 6 else thickness = 3
    ;scale_factor1 = 1.33069
    ;scale_factor1 = 1.28213
    ;scale_factor2 = 1.22765

    for j=0,N_elements(pols)-1 do begin
      for cube_i=0,1 do begin

        file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
        restore,file
        ; find the best limit in this file
        n_k=n_elements(k_edges)
        k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
        delta=power*(k^3.)/(2.*!pi^2.)
        dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
        dsigma[0] = !Values.F_INFINITY
        ;if (cube_i EQ 1) then begin
        dsigma[0] = !Values.F_INFINITY
        dsigma[*ranges[cube_i]]= !Values.F_INFINITY
        ;endif
        inds = where(dsigma EQ 0,n_count)
        if n_count GT 0 then dsigma[inds]=!Values.F_INFINITY
        ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
        limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta

        if keyword_set(scale_factor1) and cube_i EQ 0 then begin
          if j EQ 0 then scale_factor = scale_factor1 else scale_factor = scale_factor2
          limits = limits / scale_factor
          dsigma = dsigma / scale_factor
          delta= delta / scale_factor
        endif

        lim=min(limits,ind)
        header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans[cube_i]+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
        print,header

        k=k/hubble_param

        if (cube_i EQ 0) and keyword_set(pdf_plot) and (j EQ 0) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/round_fix/ssins_lowchange_zenith.pdf',/quiet,/nomatch

        ;position1=[0.1, 0.6, 0.5, 0.9]
        position1=[0.1,  0.6, 0.5,   0.9]

        position4=[0.5, 0.6, 0.9, 0.9]
        ;position4=[0.1,   0.3, 0.5,  0.6 ]

        if (cube_i EQ 0) and (j EQ 0) then begin
          position_use=position1
          xtitle_use='k (h Mpc$\up-1$)'
          ;xtitle_use=''
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          ;XTICKFORMAT_use="(I0)"
          XTICKFORMAT_use="(A1)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 0) then position_use=position1

        if (cube_i EQ 0) and (j EQ 1) then begin
          position_use=position4
          xtitle_use='k (h Mpc$\up-1$)'
          ;xtitle_use=''
          ytitle_use=''
          ;ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(A1)"
          ;YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 1) then position_use=position4

        color_use=color_array[cube_i]


        if ((cube_i EQ 0)) AND (j EQ 0) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use,$
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,/noerase
          ;XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
          XYOuts, position_use[0]+0.35, position_use[3]-.28, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
          ;cglegend,title=['post-SSINS 2$\sigma$ upper limit','noise level'],$
          cglegend,title=['RFI-removed 2$\sigma$ upper limit, zenith pointing','noise level'],$
            color=[color_array[1],color_array[1]], location=[.52,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1 ;location=[.12,.58]
        endif
        if ((cube_i EQ 0)) AND (j EQ 1) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,/noerase,YTICKFORMAT=YTICKFORMAT_use;,/noerase
          ;XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
          XYOuts, position_use[0]+0.35, position_use[3]-.28, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
          ;cglegend,title=['2016 2$\sigma$ upper limit','noise level','updated 2$\sigma upper limit','noise level'],$
          ;  color=[color_array[1],color_array[1],color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2,0,2],length=.03,vspace=1
          ;cglegend,title=['pre-SSINS 2$\sigma$ upper limit','noise level'],$
          cglegend,title=['updated 2$\sigma$ upper limit, zenith pointing','noise level'],$
            color=[color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1.1 ;location=[.12,.88]
        endif

        low_range=10^2.
        low_x_range=.14

        ;if cube_i EQ 1 then begin
        ;        error_low = (dsigma*2.); < (delta);*.99999
        ;        lows = where(delta LT 0, n_count)
        ;        if n_count GT 0 then error_low[lows]=abs(delta[lows])
        ;        lows = where((abs(delta) - error_low) LT low_range,n_count)
        ;        if n_count GT 0 then error_low[lows]=abs(delta[lows])-low_range
        ;        error_high = dsigma*2.
        ;        lows = where(k LT low_x_range,n_count)
        ;        if n_count GT 0 then error_high[lows]=!Values.F_INFINITY
        ;        if n_count GT 0 then limits[lows]=!Values.F_INFINITY
        ;        delta_k=k[1]/2.
        ;        for k_i=0, n_k-2 do $
        ;          cgColorFill, [k[k_i]-delta_k, k[k_i]+delta_k, k[k_i]+delta_k,k[k_i]-delta_k], $
        ;          [limits[k_i], limits[k_i], abs(delta[k_i])-error_low[k_i],abs(delta[k_i])-error_low[k_i]], $
        ;          Color='grey',/checkforfinite
        ;endif

        if (cube_i NE 0) then $
          cgoplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.7], yrange=[10^2.,10^8.],color=color_use,thick=thickness,position=position_use
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use

      endfor
    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif

  ;ab_limit_compare=1
  ;NOT USED
  if keyword_set(ab_limit_compare) then begin
    ;coarse_harm_width = 3
    ;wedge_angles=120d
    n_cubes=6

    cut = STRARR(n_cubes)
    cut[0]='wedge_cut_plus_res_cut'
    cut[1] = 'sub_cubes'
    cut[2]='wedge_cut_plus_res_cut'
    cut[3] = 'sub_cubes'
    cut[4]='wedge_cut_plus_res_cut'
    cut[5] = 'sub_cubes'



    chans = strarr(n_cubes)
    chans[0]='ch0-95_' & band='low'
    chans[1]='ch0-95_' & band='low'
    chans[2]='ch48-143_' & band='mid'
    chans[3]='ch48-143_' & band='mid'
    chans[4]='ch96-191_' & band='high'
    chans[5]='ch96-191_' & band='high'

    win = strarr(n_cubes)
    win[0]=''
    win[1]='fullimg_'
    win[2]=''
    win[3]='fullimg_'
    win[4]=''
    win[5]='fullimg_'

    spec_window='bh_'

    basefile = STRARR(n_cubes)
    basefile[0]=['/Users/nabarry/MWA/data/fhd_apb_EoR0_high_sem1_1/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[0]]
    basefile[1]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans[1]]
    basefile[2]=['/Users/nabarry/MWA/data/fhd_apb_EoR0_high_sem1_1/Combined_obs_'+cut[2]+'_cubeXX__even_odd_joint_'+win[2]+chans[2]]
    basefile[3]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/Combined_obs_'+cut[3]+'_cubeXX__even_odd_joint_'+win[3]+chans[3]]
    basefile[4]=['/Users/nabarry/MWA/data/fhd_apb_EoR0_high_sem1_1/Combined_obs_'+cut[4]+'_cubeXX__even_odd_joint_'+win[4]+chans[4]]
    basefile[5]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/Combined_obs_'+cut[5]+'_cubeXX__even_odd_joint_'+win[5]+chans[5]]

    ranges=PTRARR(n_cubes,/allocate)
    *ranges[1]=[0,1,6,7,8,9,17,18,19,20,21,22,29,30,31,32,33,34,46]+4
    *ranges[3]=[0,1,6,7,8,9,17,18,19,20,21,22,29,30,31,32,33,34,46]+4
    *ranges[5]=[0,1,6,7,8,9,17,18,19,20,21,22,29,30,31,32,33,34,46]+4
    ;*ranges[3]=[6,7,8,9,10,11,21,22,23,24,25,26,36,37,38,39,40,41,50,51,52]

    cubes=['res']
    pols=['xx','yy']

    endfile = STRARR(n_cubes)
    endfile[0]='_bh_dencorr_no_120deg_wedge_cbw3_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'
    endfile[1]='_averemove_sw'+spec_window+'dencorr_no_120deg_wedge_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'
    endfile[2]='_bh_dencorr_no_120deg_wedge_cbw3_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'
    endfile[3]='_averemove_sw'+spec_window+'dencorr_no_120deg_wedge_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'
    endfile[4]='_bh_dencorr_no_120deg_wedge_cbw3_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'
    endfile[5]='_averemove_sw'+spec_window+'dencorr_no_120deg_wedge_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'

    limit_percent=0.9772 ;2 sigma
    hubble_param=0.71

    color_num = [10,12,14,16,18,20]
    ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101],[198,118,63],[197,120,62]]

    for n_i=0, n_cubes-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B,20B]

    pdf_plot=0
    if keyword_set(pdf_plot) then thickness = 5 else thickness = 3

    !Y.OMargin = [2, 8]
    !X.OMargin = [2, 6]
    !P.Multi = [0, 2, 3]

    for j=0,N_elements(pols)-1 do begin
      for cube_i=0,n_cubes-1 do begin

        file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
        restore,file
        ; find the best limit in this file
        n_k=n_elements(k_edges)
        k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
        delta=power*(k^3.)/(2.*!pi^2.)
        dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
        if (cube_i EQ 1) OR (cube_i EQ 3) OR (cube_i EQ 5) then begin
          dsigma[0] = !Values.F_INFINITY
          dsigma[*ranges[cube_i]]= !Values.F_INFINITY
        endif
        ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
        limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
        lim=min(limits,ind)
        header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans[cube_i]+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
        print,header

        ;      if (cube_i NE 3) AND (cube_i NE 1) then begin
        ;      ;remove first bins
        ;      ;limits=limits[3:N_elements(limits)-1]
        ;      limits[*ranges[cube_i]] = !Values.F_NAN
        ;      ;k=k[3:N_elements(k)-1]/hubble_param
        k=k/hubble_param
        ;      endif else begin
        ;
        ;        limits=limits[3:N_elements(limits)-1]
        ;      limits[*ranges[cube_i]] = !Values.F_NAN
        ;
        ;      k=k[3:N_elements(k)-1]/hubble_param
        ;        endelse



        if (cube_i EQ 0) and keyword_set(pdf_plot) and (j EQ 0) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/limit_AB_compare.png',/quiet,/nomatch

        position1=[0.1, 0.7, 0.5, 0.9]
        position2=[0.1, 0.5, 0.5, 0.7]
        position3=[0.1, 0.3, 0.5, 0.5]

        position4=[0.5, 0.7, 0.9, 0.9]
        position5=[0.5, 0.5, 0.9, 0.7]
        position6=[0.5, 0.3, 0.9, 0.5]

        if (cube_i EQ 0) and (j EQ 0) then begin
          position_use=position1
          xtitle_use=''
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          XTICKFORMAT_use="(A1)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 0) then position_use=position1
        if (cube_i EQ 2) and (j EQ 0) then begin
          position_use=position2
          xtitle_use=''
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          XTICKFORMAT_use="(A1)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 3) and (j EQ 0) then position_use=position2
        if (cube_i EQ 4) and (j EQ 0) then begin
          position_use=position3
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 5) and (j EQ 0) then position_use=position3

        if (cube_i EQ 0) and (j EQ 1) then begin
          position_use=position4
          xtitle_use=''
          ytitle_use=''
          XTICKFORMAT_use="(A1)"
          YTICKFORMAT_use="(A1)"
        endif
        if (cube_i EQ 1) and (j EQ 1) then position_use=position4
        if (cube_i EQ 2) and (j EQ 1) then begin
          position_use=position5
          xtitle_use=''
          ytitle_use=''
          XTICKFORMAT_use="(A1)"
          YTICKFORMAT_use="(A1)"
        endif
        if (cube_i EQ 3) and (j EQ 1) then position_use=position5
        if (cube_i EQ 4) and (j EQ 1) then begin
          position_use=position6
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use=''
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(A1)"
        endif
        if (cube_i EQ 5) and (j EQ 1) then position_use=position6

        color_use=color_array[cube_i]

        if ((cube_i EQ 0) OR (cube_i EQ 2) OR (cube_i EQ 4)) AND (j EQ 0) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.16,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =1.25, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]-.02, 'This is the Plot Title', /Normal, Alignment=0.5, Charsize=1
        endif
        if ((cube_i EQ 0) OR (cube_i EQ 2) OR (cube_i EQ 4)) AND (j EQ 1) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.16,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =1.25, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]-.02, 'This is the Plot Title', /Normal, Alignment=0.5, Charsize=1
        endif

        if (cube_i NE 0) then $
          cgoplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.15,1.7], yrange=[10^2.,10^8.],color=color_use,thick=thickness,position=position_use
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use

      endfor
      stop
    endfor
    stop
  endif


  ;all_limit_compare_1band=1
  ;NOT USED
  if keyword_set(all_limit_compare_1band) then begin
    ;coarse_harm_width = 3
    ;wedge_angles=120d
    n_cubes=3

    cut = STRARR(n_cubes)
    cut[*]='btl_noalltv_noocc4'

    chans = strarr(n_cubes)
    chans[0]='ch9-126_' & band='low'
    chans[1]='ch9-126_' & band='low'
    chans[2]='ch9-126_' & band='low'
    ;chans[3]='ch9-123_' & band='low'
    ;chans[4]='ch9-122_' & band='low'
    ;chans[5]='ch9-121_' & band='low'
    ;chans[6]='ch8-120_' & band='low'


    win = strarr(n_cubes)
    win[*]='fullimg_'

    spec_window='bh_'

    basefile = STRARR(n_cubes)
    basefile[0]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/even_odd_fix/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[0]]
    basefile[1]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/round_fix/force_dft/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans[1]]
    basefile[2]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/round_fix/Combined_obs_'+cut[2]+'_cubeXX__even_odd_joint_'+win[2]+chans[2]]
    ;basefile[3]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/even_odd_fix/Combined_obs_'+cut[3]+'_cubeXX__even_odd_joint_'+win[3]+chans[3]]
    ;basefile[4]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/even_odd_fix/Combined_obs_'+cut[4]+'_cubeXX__even_odd_joint_'+win[4]+chans[4]]
    ;basefile[5]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/even_odd_fix/Combined_obs_'+cut[5]+'_cubeXX__even_odd_joint_'+win[5]+chans[5]]
    ;basefile[6]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/even_odd_fix/Combined_obs_'+cut[6]+'_cubeXX__even_odd_joint_'+win[6]+chans[6]]


    ranges=PTRARR(n_cubes,/allocate)
    *ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    *ranges[1]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    *ranges[2]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[3]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[4]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[5]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[6]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4


    cubes=['res']
    pols=['xx','yy']

    endfile = STRARR(n_cubes)
    endfile[*]='_averemove_sw'+spec_window+'dencorr_no_115deg_wedge_kperplambda15-80_kpar0.15-200_1dkpower.idlsave'

    limit_percent=0.9772 ;2 sigma
    hubble_param=0.71

    color_num = [10,12,14,16,18,20,22]
    ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101],[198,118,63],[197,120,62],[203,86,131]]

    for n_i=0, n_cubes-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B,20B,22B]

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 6 else thickness = 3

    for j=0,N_elements(pols)-1 do begin
      for cube_i=0,n_cubes-1 do begin

        file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
        restore,file
        ; find the best limit in this file
        n_k=n_elements(k_edges)
        k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
        delta=power*(k^3.)/(2.*!pi^2.)
        dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
        dsigma[0] = !Values.F_INFINITY
        ;if (cube_i EQ 1) then begin
        dsigma[0] = !Values.F_INFINITY
        dsigma[*ranges[cube_i]]= !Values.F_INFINITY
        ;endif
        inds = where(dsigma EQ 0,n_count)
        if n_count GT 0 then dsigma[inds]=!Values.F_INFINITY
        ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
        limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
        lim=min(limits,ind)
        header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans[cube_i]+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
        print,header

        k=k/hubble_param

        if (cube_i EQ 0) and keyword_set(pdf_plot) and (j EQ 0) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/freq_compare_even_odd_fix.pdf',/quiet,/nomatch

        position1=[0.1, 0.6, 0.5, 0.9]

        position4=[0.5, 0.6, 0.9, 0.9]

        if (cube_i EQ 0) and (j EQ 0) then begin
          position_use=position1
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 0) then position_use=position1

        if (cube_i EQ 0) and (j EQ 1) then begin
          position_use=position4
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use=''
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(A1)"
        endif
        if (cube_i EQ 1) and (j EQ 1) then position_use=position4

        color_use=color_array[cube_i]

        if ((cube_i EQ 0)) AND (j EQ 0) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
          cglegend,title=['ch9-126','ch9-125','ch9-124'],$
            color=[color_array[0],color_array[1],color_array[2]], location=[.52,.88], charsize=.7, thick=thickness,linestyle=0,length=.03,vspace=1
        endif
        if ((cube_i EQ 0)) AND (j EQ 1) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
          ;cglegend,title=['2016 2$\sigma$ upper limit','noise level','updated 2$\sigma upper limit','noise level'],$
          ;  color=[color_array[1],color_array[1],color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2,0,2],length=.03,vspace=1
          cglegend,title=['ch9-126','ch9-125','ch9-124'],$
            color=[color_array[0],color_array[1],color_array[2]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=0,length=.03,vspace=1.1
        endif

        if (cube_i NE 0) then $
          cgoplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.7], yrange=[10^2.,10^8.],color=color_use,thick=thickness,position=position_use
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use

      endfor

      ;#Limit: 22408.321 mK^2, at k = 0.23044845 h Mpc^-1 (wedge_cut_plus_res_cut ch9-127_ res  xx)
      ;#Limit: 9223.0543 mK^2, at k = 0.23044845 h Mpc^-1 (sub_cubes ch9-127_ res fullimg_ xx)
      ;2.42960 better
      ;#Limit: 27702.209 mK^2, at k = 0.23044845 h Mpc^-1 (wedge_cut_plus_res_cut ch9-127_ res  yy)
      ;#Limit: 7904.3834 mK^2, at k = 0.23044845 h Mpc^-1 (sub_cubes ch9-127_ res fullimg_ yy)
      ;3.50466 better

      ;update
      ;#Limit: 24231.534 mK^2, at k = 0.2342833 h Mpc^-1 (wedge_cut_plus_res_cut ch9-125_ res  xx)
      ;#Limit: 9774.8649 mK^2, at k = 0.2342833 h Mpc^-1 (sub_cubes ch9-125_ res fullimg_ xx)
      ;2.47896 better
      ;#Limit: 28403.79 mK^2, at k = 0.2342833 h Mpc^-1 (wedge_cut_plus_res_cut ch9-125_ res  yy)
      ;#Limit: 8086.5888 mK^2, at k = 0.2342833 h Mpc^-1 (sub_cubes ch9-125_ res fullimg_ yy)
      ;3.51246 better

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif


  ;limit_compare_dtv=1
  if keyword_set(limit_compare_dtv) then begin
    ;coarse_harm_width = 3
    ;wedge_angles=120d
    n_cubes=2

    cut = STRARR(n_cubes)
    cut[0]='sub_cubes'
    cut[1] = 'btl_noalltv_noocc4'

    chans = strarr(n_cubes)
    chans[0]='ch9-126_' & band='low'
    chans[1]='ch9-126_' & band='low'

    win = strarr(n_cubes)
    win[0]='fullimg_'
    ;win[1]='fullimg_'
    win[1]='fullimg_'

    spec_window='bh_'

    basefile = STRARR(n_cubes)
    basefile[0]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/round_fix/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[0]]
    basefile[1]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/round_fix/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans[1]]

    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[1]=[0,1,6,7,8,9,17,18,19,20,21,22,29,30,31,32,33,34,46]
    ;*ranges[1]=[7,8,9,10,11,22,23,24,25,26,37,38,39,40,41,52,53,54]+6 ;cbw3
    *ranges[1]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    *ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[1]=[0]
    ;*ranges[3]=[6,7,8,9,10,11,21,22,23,24,25,26,36,37,38,39,40,41,50,51,52]

    cubes=['res']
    pols=['xx','yy']

    endfile = STRARR(n_cubes)
    endfile[0]='_averemove_sw'+spec_window+'dencorr_no_115deg_wedge_kperplambda15-80_kpar0.15-200_1dkpower.idlsave'
    endfile[1]='_averemove_sw'+spec_window+'dencorr_no_115deg_wedge_kperplambda15-80_kpar0.15-200_1dkpower.idlsave'

    limit_percent=0.9772 ;2 sigma
    hubble_param=0.71

    color_num = [10,12,14,16,18,20]
    ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101],[198,118,63],[197,120,62]]

    for n_i=0, n_cubes-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B,20B]

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 6 else thickness = 3

    for j=0,N_elements(pols)-1 do begin
      for cube_i=0,n_cubes-1 do begin

        file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
        restore,file
        ; find the best limit in this file
        n_k=n_elements(k_edges)
        k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
        delta=power*(k^3.)/(2.*!pi^2.)
        dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
        dsigma[0] = !Values.F_INFINITY
        ;if (cube_i EQ 1) then begin
        dsigma[0] = !Values.F_INFINITY
        dsigma[*ranges[cube_i]]= !Values.F_INFINITY
        ;endif
        inds = where(dsigma EQ 0,n_count)
        if n_count GT 0 then dsigma[inds]=!Values.F_INFINITY
        ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
        limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
        lim=min(limits,ind)
        header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans[cube_i]+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
        print,header

        k=k/hubble_param

        if (cube_i EQ 0) and keyword_set(pdf_plot) and (j EQ 0) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/limit_compare_dtv_power.pdf',/quiet,/nomatch

        position1=[0.1, 0.6, 0.5, 0.9]

        position4=[0.5, 0.6, 0.9, 0.9]

        if (cube_i EQ 0) and (j EQ 0) then begin
          position_use=position1
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 0) then position_use=position1

        if (cube_i EQ 0) and (j EQ 1) then begin
          position_use=position4
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use=''
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(A1)"
        endif
        if (cube_i EQ 1) and (j EQ 1) then position_use=position4

        color_use=color_array[cube_i]

        if ((cube_i EQ 0)) AND (j EQ 0) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
          cglegend,title=['updated 2$\sigma$ upper limit, DTV removed','noise level'],$
            color=[color_array[1],color_array[1]], location=[.52,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1
        endif
        if ((cube_i EQ 0)) AND (j EQ 1) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
          ;cglegend,title=['2016 2$\sigma$ upper limit','noise level','updated 2$\sigma upper limit','noise level'],$
          ;  color=[color_array[1],color_array[1],color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2,0,2],length=.03,vspace=1
          cglegend,title=['updated 2$\sigma$ upper limit','noise level'],$
            color=[color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1.1
        endif

        if (cube_i NE 0) then $
          cgoplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.7], yrange=[10^2.,10^8.],color=color_use,thick=thickness,position=position_use
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use

      endfor

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif


  ;limit_compare_fullimg=1
  ;NOT USED
  if keyword_set(limit_compare_fullimg) then begin
    ;coarse_harm_width = 3
    ;wedge_angles=120d
    n_cubes=2

    cut = STRARR(n_cubes)
    cut[0]='Aug23_longrunstyle'
    cut[1] = 'Aug23_longrunstyle'

    chans = strarr(n_cubes)
    chans[0]='' & band='low'
    chans[1]='' & band='low'

    win = strarr(n_cubes)
    win[0]='fullimg_'
    ;win[1]='fullimg_'
    win[1]=''

    spec_window='bh_'

    basefile = STRARR(n_cubes)
    basefile[0]=['Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[0]]
    basefile[1]=['Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans[1]]

    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[1]=[0,1,6,7,8,9,17,18,19,20,21,22,29,30,31,32,33,34,46]
    ;*ranges[1]=[7,8,9,10,11,22,23,24,25,26,37,38,39,40,41,52,53,54]+6 ;cbw3
    ;*ranges[1]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[1]=[0]
    ;*ranges[3]=[6,7,8,9,10,11,21,22,23,24,25,26,36,37,38,39,40,41,50,51,52]

    cubes=['res']
    pols=['xx','yy']

    endfile = STRARR(n_cubes)
    endfile[0]='_averemove_sw'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'
    endfile[1]='_averemove_sw'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'

    limit_percent=0.9772 ;2 sigma
    hubble_param=0.71

    color_num = [10,12,14,16,18,20]
    ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101],[198,118,63],[197,120,62]]

    for n_i=0, n_cubes-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B,20B]

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 6 else thickness = 3

    for j=0,N_elements(pols)-1 do begin
      for cube_i=0,n_cubes-1 do begin

        file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
        restore,file
        ; find the best limit in this file
        n_k=n_elements(k_edges)
        k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
        delta=power*(k^3.)/(2.*!pi^2.)
        dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
        dsigma[0] = !Values.F_INFINITY
        ;if (cube_i EQ 1) then begin
        dsigma[0] = !Values.F_INFINITY
        dsigma[*ranges[cube_i]]= !Values.F_INFINITY
        ;endif
        inds = where(dsigma EQ 0,n_count)
        if n_count GT 0 then dsigma[inds]=!Values.F_INFINITY
        ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
        limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
        lim=min(limits,ind)
        header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans[cube_i]+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
        print,header

        k=k/hubble_param

        if (cube_i EQ 0) and keyword_set(pdf_plot) and (j EQ 0) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/limit_compare_fullimg.pdf',/quiet,/nomatch

        position1=[0.1, 0.6, 0.5, 0.9]

        position4=[0.5, 0.6, 0.9, 0.9]

        if (cube_i EQ 0) and (j EQ 0) then begin
          position_use=position1
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 0) then position_use=position1

        if (cube_i EQ 0) and (j EQ 1) then begin
          position_use=position4
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use=''
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(A1)"
        endif
        if (cube_i EQ 1) and (j EQ 1) then position_use=position4

        color_use=color_array[cube_i]

        if ((cube_i EQ 0)) AND (j EQ 0) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
          cglegend,title=['updated 2$\sigma$ upper limit, full img','noise level'],$
            color=[color_array[1],color_array[1]], location=[.52,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1
        endif
        if ((cube_i EQ 0)) AND (j EQ 1) then begin
          cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
          ;cglegend,title=['2016 2$\sigma$ upper limit','noise level','updated 2$\sigma upper limit','noise level'],$
          ;  color=[color_array[1],color_array[1],color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2,0,2],length=.03,vspace=1
          cglegend,title=['updated 2$\sigma$ upper limit, postage','noise level'],$
            color=[color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1.1
        endif

        if (cube_i NE 0) then $
          cgoplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.7], yrange=[10^2.,10^8.],color=color_use,thick=thickness,position=position_use
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use

      endfor

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif


  ;k0_compare=1
  if keyword_set(k0_compare) then begin
    ;coarse_harm_width = 3
    ;wedge_angles=120d
    n_cubes=2

    cut = STRARR(n_cubes)
    cut[0]='wedge_cut_plus_res_cut'
    cut[1] = 'sub_cubes'

    chans = strarr(n_cubes)
    chans[0]='ch9-125_' & band='low'
    chans[1]='ch9-125_' & band='low'

    win = strarr(n_cubes)
    win[0]=''
    ;win[1]='fullimg_'
    win[1]='fullimg_'

    spec_window='bh_'

    basefile = STRARR(n_cubes)
    basefile[0]=['/Users/nabarry/MWA/data/fhd_apb_EoR0_high_sem1_1/Combined_obs_wedge_cut_plus_res_cut_cubeXX__even_odd_joint_ch9-125_dirty_']
    basefile[1]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/Combined_obs_sub_cubes_cubeXX__even_odd_joint_fullimg_ch9-125_dirty_']

    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[1]=[0,1,6,7,8,9,17,18,19,20,21,22,29,30,31,32,33,34,46]
    ;*ranges[1]=[7,8,9,10,11,22,23,24,25,26,37,38,39,40,41,52,53,54]+6 ;cbw3
    ;*ranges[1]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[1]=[0]
    ;*ranges[3]=[6,7,8,9,10,11,21,22,23,24,25,26,36,37,38,39,40,41,50,51,52]

    cubes=['res']
    pols=['xx','yy']

    endfile = STRARR(n_cubes)
    endfile[0]='_swbh_dencorr_k0power.idlsave'
    endfile[1]='_swbh_dencorr_k0power.idlsave'

    limit_percent=0.9772 ;2 sigma
    hubble_param=0.71

    color_num = [10,12,14,16,18,20]
    ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101],[198,118,63],[197,120,62]]

    for n_i=0, n_cubes-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B,20B]

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 6 else thickness = 3

    for j=0,N_elements(pols)-1 do begin
      for cube_i=0,n_cubes-1 do begin

        file=basefile[cube_i]+pols[j]+endfile[cube_i]
        restore,file
        ; find the best limit in this file
        n_k=n_elements(k_edges)
        k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
        delta=power*(k^3.)/(2.*!pi^2.)
        dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
        dsigma[0] = !Values.F_INFINITY
        ;if (cube_i EQ 1) then begin
        dsigma[0] = !Values.F_INFINITY
        dsigma[*ranges[cube_i]]= !Values.F_INFINITY
        ;endif
        inds = where(dsigma EQ 0,n_count)
        if n_count GT 0 then dsigma[inds]=!Values.F_INFINITY
        ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
        limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
        lim=min(limits,ind)
        header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans[cube_i]+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
        print,header

        k=k/hubble_param

        if (cube_i EQ 0) and keyword_set(pdf_plot) and (j EQ 0) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/k0_compare.pdf',/quiet,/nomatch

        position1=[0.1, 0.6, 0.5, 0.9]

        position4=[0.5, 0.6, 0.9, 0.9]

        if (cube_i EQ 0) and (j EQ 0) then begin
          position_use=position1
          xtitle_use='k perp (h Mpc$\up-1$)'
          ytitle_use='P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)'
          ;XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 0) then position_use=position1

        if (cube_i EQ 0) and (j EQ 1) then begin
          position_use=position4
          xtitle_use='k perp (h Mpc$\up-1$)'
          ytitle_use=''
          ;XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(A1)"
        endif
        if (cube_i EQ 1) and (j EQ 1) then position_use=position4

        color_use=color_array[cube_i]

        if ((cube_i EQ 0)) AND (j EQ 0) then begin
          cgplot, k,power,/xlog,/ylog,psym=10, xrange=[0.001,.2], yrange=[10^13.,10^17.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,/noerase;,XTICKFORMAT=XTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
          cglegend,title=['updated k0'],$
            color=[color_array[1]], location=[.52,.88], charsize=.7, thick=thickness,linestyle=[0],length=.03,vspace=1
        endif
        if ((cube_i EQ 0)) AND (j EQ 1) then begin
          cgplot, k,power,/xlog,/ylog,psym=10, xrange=[0.001,.2], yrange=[10^13.,10^17.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
          ;cglegend,title=['2016 2$\sigma$ upper limit','noise level','updated 2$\sigma upper limit','noise level'],$
          ;  color=[color_array[1],color_array[1],color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2,0,2],length=.03,vspace=1
          cglegend,title=['2016 k0'],$
            color=[color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0],length=.03,vspace=1.1
        endif

        if (cube_i NE 0) then $
          cgoplot, k,power,/xlog,/ylog,psym=10, xrange=[0.001,.2], yrange=[10^13.,10^17.],color=color_use,thick=thickness,position=position_use
        ;cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use

      endfor

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif


  ;sim_models_freqband=1
  ;NOT USED
  if keyword_set(sim_models_freqband) then begin
    ;coarse_harm_width = 3
    ;wedge_angles=120d
    n_cubes=3

    cut = STRARR(n_cubes)
    ;cut[0]='zenith'
    ;cut[1] = 'zenith'
    ;cut[2] = 'zenith'
    cut[0]='btl_noalltv_noocc4'
    cut[1] = 'btl_noalltv_noocc4'
    cut[2] = 'btl_noalltv_noocc4'

    chans = strarr(n_cubes)
    ;chans[0]='ch9-125_' & band='low'
    ;chans[1]='ch9-126_' & band='low'
    chans[0]='ch9-126_' & band='low'
    chans[1]='ch9-126_' & band='low'
    chans[2]='ch9-126_' & band='low'

    win = strarr(n_cubes)
    win[0]='fullimg_'
    ;win[1]='fullimg_'
    win[1]='fullimg_'
    win[2]='fullimg_'

    spec_window = strarr(n_cubes)
    spec_window[0]='bh_'
    spec_window[1]='bh_'
    spec_window[2]='bh_'

    horizon = strarr(n_cubes)
    horizon[0] = '115deg_'
    horizon[1] = '115deg_'
    horizon[2] = '115deg_'

    basefile = STRARR(n_cubes)
    ;basefile[0]=['/Users/nabarry/MWA/data/fhd_nb_model_GLEAM_nosidelobes/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[0]]
    ;basefile[1]=['/Users/nabarry/MWA/data/fhd_nb_model_GLEAM_nosidelobes/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans[1]]
    ;basefile[2]=['/Users/nabarry/MWA/data/fhd_nb_model_GLEAM_nosidelobes/Combined_obs_'+cut[2]+'_cubeXX__even_odd_joint_'+win[2]+chans[2]]
    basefile[0]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_dft_tests_zupdate/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[0]]
    basefile[1]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_rotate_then_fold/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[1]]
    basefile[2]=['/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_nols/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[2]]

    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[1]=[0,1,6,7,8,9,17,18,19,20,21,22,29,30,31,32,33,34,46]
    ;*ranges[1]=[7,8,9,10,11,22,23,24,25,26,37,38,39,40,41,52,53,54]+6 ;cbw3
    ;*ranges[1]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[1]=[0]
    ;*ranges[3]=[6,7,8,9,10,11,21,22,23,24,25,26,36,37,38,39,40,41,50,51,52]

    cubes=['res']
    pols=['xx','yy']

    endfile = STRARR(n_cubes)
    ;endfile[0]='_averemove_sw'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'
    ;endfile[1]='_averemove_sw'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'
    ;endfile[2]='_averemove_sw'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-70_kpar0.15-200_1dkpower.idlsave'
    endfile[0]='_averemove_sw'+spec_window[0]+'dencorr_no_'+horizon[0]+'wedge_kperplambda15-80_kpar0.2-200_1dkpower.idlsave'
    endfile[1]='_averemove_sw'+spec_window[1]+'dencorr_no_'+horizon[1]+'wedge_kperplambda15-80_kpar0.15-200_1dkpower.idlsave'
    endfile[2]='_averemove_sw'+spec_window[2]+'dencorr_no_'+horizon[2]+'wedge_kperplambda15-80_kpar0.15-200_1dkpower.idlsave'

    limit_percent=0.9772 ;2 sigma
    hubble_param=0.71

    color_num = [10,12,14,16,18,20]
    ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101],[198,118,63],[197,120,62]]

    for n_i=0, n_cubes-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B,20B]

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 6 else thickness = 3

    for j=0,N_elements(pols)-1 do begin
      for cube_i=0,n_cubes-1 do begin

        file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
        restore,file
        ; find the best limit in this file
        n_k=n_elements(k_edges)
        k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
        delta=power*(k^3.)/(2.*!pi^2.)
        dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
        dsigma[0] = !Values.F_INFINITY
        ;if (cube_i EQ 1) then begin
        dsigma[0] = !Values.F_INFINITY
        dsigma[*ranges[cube_i]]= !Values.F_INFINITY
        ;endif
        inds = where(dsigma EQ 0,n_count)
        if n_count GT 0 then dsigma[inds]=!Values.F_INFINITY
        ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
        limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
        lim=min(limits,ind)
        header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans[cube_i]+' ' +horizon[cube_i]+' '+spec_window[cube_i]+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
        print,header

        k=k/hubble_param

        if (cube_i EQ 0) and keyword_set(pdf_plot) and (j EQ 0) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/res_limit_fftchange.pdf',/quiet,/nomatch

        position1=[0.1, 0.6, 0.5, 0.9]

        position4=[0.5, 0.6, 0.9, 0.9]

        if (cube_i EQ 0) and (j EQ 0) then begin
          position_use=position1
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(I0)"
        endif
        if (cube_i EQ 1) and (j EQ 0) then position_use=position1

        if (cube_i EQ 0) and (j EQ 1) then begin
          position_use=position4
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use=''
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(A1)"
        endif
        if (cube_i EQ 1) and (j EQ 1) then position_use=position4

        color_use=color_array[cube_i]

        if ((cube_i EQ 0)) AND (j EQ 0) then begin
          cgplot, k,delta,/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^7.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W, z=7', /Normal, Alignment=0.5, Charsize=.8
          cglegend,title=['bh2','noise level'],$
            color=[color_array[2],color_array[1]], location=[.52,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1
        endif
        if ((cube_i EQ 0)) AND (j EQ 1) then begin
          cgplot, k,delta,/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^7.],ytitle=ytitle_use, $
            xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT=XTICKFORMAT_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
          XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
          ;cglegend,title=['2016 2$\sigma$ upper limit','noise level','updated 2$\sigma upper limit','noise level'],$
          ;  color=[color_array[1],color_array[1],color_array[0],color_array[0]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2,0,2],length=.03,vspace=1
          cglegend,title=['blm','bh'],$
            color=[color_array[0],color_array[1]], location=[.12,.88], charsize=.7, thick=thickness,linestyle=[0,2],length=.03,vspace=1.1
        endif

        if (cube_i NE 0) then $
          cgoplot, k,delta,/xlog,/ylog,psym=10, xrange=[.14,1.7], yrange=[10^2.,10^8.],color=color_use,thick=thickness,position=position_use
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use

      endfor

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif

  ;noise_check_plot=1
  if keyword_set(noise_check_plot) then begin
    ;coarse_harm_width = 3
    ;wedge_angles=120d
    n_cubes=3

    cut = STRARR(n_cubes)
    cut[0]='1061316296'
    cut[1] = 'zenith'
    cut[2] = 'zenith'
    ;cut[3] = '1061316296'

    chans = strarr(n_cubes)
    chans[0]='ch9-126_' & band='low'
    chans[1]='ch9-126_' & band='low'
    chans[2]='ch9-126_' & band='low'

    win = strarr(n_cubes)
    win[0]='noimgclip_'
    ;win[1]='fullimg_'
    ;win[1]='fullimg_'
    win[2]='fullimg_'
    ;win[3]='fullimg_'

    spec_window='bh_'

    basefile = STRARR(n_cubes)
    basefile[0]=['/Users/nabarry/MWA/data/uvf_noise_check/'+cut[0]+'_gridded_uvf__even_odd_joint_'+win[0]];+chans[0]]
    basefile[1]=['/Users/nabarry/MWA/data/uvf_noise_check/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]];+chans[1]]
    basefile[2]=['/Users/nabarry/MWA/data/noise_check/BH2/weights/Combined_obs_'+cut[2]+'_cubeXX__even_odd_joint_'+win[2]]
    ;basefile[3]=['/Users/nabarry/MWA/data/uvf_noise_check/uvf_tapered_weights/'+cut[3]+'_gridded_uvf__even_odd_joint_'+win[3]];+chans[0]]


    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[1]=[0,1,6,7,8,9,17,18,19,20,21,22,29,30,31,32,33,34,46]
    ;*ranges[1]=[7,8,9,10,11,22,23,24,25,26,37,38,39,40,41,52,53,54]+6 ;cbw3
    ;*ranges[1]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4
    ;*ranges[1]=[0]
    ;*ranges[3]=[6,7,8,9,10,11,21,22,23,24,25,26,36,37,38,39,40,41,50,51,52]

    cubes=['res']
    pols=['xx','yy']

    endfile = STRARR(n_cubes)
    endfile[0]='_averemove_sw'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'
    endfile[1]='_averemove_sw'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'
    endfile[2]='_averemove_sw'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'
    ;endfile[3]='_averemove_sw'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'

    limit_percent=0.9772 ;2 sigma
    hubble_param=0.71

    color_num = [10,12,14,16,18,20]
    ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    rgbcolors=[[160,32,240],[114,166,89],[202,86,139],[198,118,63],[197,120,62]]

    ;    for n_i=0, n_cubes-1 do $
    ;      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    ;    color_array = [10B,12B,14B,16B,18B,20B]
    color_array=['purple','forest green','orange','blue']
    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 7 else thickness = 3

    for j=0,N_elements(pols)-1 do begin
      for cube_i=0,n_cubes-1 do begin

        file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
        restore,file
        ; find the best limit in this file
        n_k=n_elements(k_edges)
        k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
        delta=power*(k^3.)/(2.*!pi^2.)
        dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
        dsigma[0] = !Values.F_INFINITY
        ;if (cube_i EQ 1) then begin
        dsigma[0] = !Values.F_INFINITY
        dsigma[*ranges[cube_i]]= !Values.F_INFINITY
        ;endif
        inds = where(dsigma EQ 0,n_count)
        if n_count GT 0 then dsigma[inds]=!Values.F_INFINITY
        ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
        limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
        lim=min(limits,ind)
        header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans[cube_i]+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
        print,header

        k=k/hubble_param

        if (cube_i EQ 0) and keyword_set(pdf_plot) and (j EQ 0) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/noise_check_weights_update.pdf',/quiet,/nomatch


        if (cube_i EQ 0) and (j EQ 0) then begin
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
          ;XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(I0)"
        endif

        if (cube_i EQ 0) and (j EQ 1) then begin
          xtitle_use='k (h Mpc$\up-1$)'
          ytitle_use=''
          XTICKFORMAT_use="(I0)"
          YTICKFORMAT_use="(A1)"
        endif

        color_use=color_array[cube_i]

        ;limits=delta

        if ((cube_i EQ 0)) AND (j EQ 0) then begin
          cgplot, k,dsigma,/xlog,/ylog,psym=10, xrange=[.05,1.2], yrange=[10^3.,10^9.],ytitle=ytitle_use, aspect=0.5, $
            xtitle=xtitle_use, color=color_use,thick=thickness,linestyle=2,/noerase,title='Noise contamination, off-zenith',charsize=1.3
          ;XYOuts, position_use[0]+0.2, position_use[3]+.02, 'E-W', /Normal, Alignment=0.5, Charsize=.8
          cglegend,title=['observation noise','AB2016 observation noise','updated observation noise'],$
            ;  cglegend,title=['uvf','postage stamp image','BH2 kernel', 'BH2 kernel uvf'],$
            color=[color_array[0],color_array[1],color_array[2]], location=[.16,.7], charsize=1.1, thick=thickness,linestyle=[2,2,2];,length=.03;,vspace=1
        endif


        if (cube_i NE 0) then $
          cgoplot, k,dsigma,/xlog,/ylog,psym=10, xrange=[.14,1.7], yrange=[10^2.,10^9.],color=color_use,thick=thickness,linestyle=2
        ;cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use
        ;stop

      endfor
      if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
      stop

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif

  ;ssins_gps_obs=1
  if keyword_set(ssins_gps_obs) then begin
    readcol, '/Users/nabarry/MWA/FHD/obs_list/beardsley_thesis_list.txt', btl, FORMAT='(A10)'
    readcol, '/Users/nabarry/MWA/FHD/obs_list/beardsley_thesis_list_jd.txt', btl_jd, FORMAT='(D10.10)'
    readcol, '/Users/nabarry/MWA/FHD/obs_list/alltv.txt', alltv, FORMAT='(A10)'
    readcol, '/Users/nabarry/MWA/FHD/obs_list/btl_noalltv.txt', btl_noalltv, FORMAT='(A10)'
    readcol, '/Users/nabarry/MWA/FHD/obs_list/tv50.txt', tv50, FORMAT='(A10)'
    readcol, '/Users/nabarry/MWA/FHD/obs_list/tv100.txt', tv100, FORMAT='(A10)'
    readcol, '/Users/nabarry/MWA/FHD/obs_list/tv200.txt', tv200, FORMAT='(A10)'
    readcol, '/Users/nabarry/MWA/FHD/obs_list/occupancy4.txt', occupancy4, FORMAT='(A10)'
    readcol, '/Users/nabarry/MWA/FHD/obs_list/Long_Run_LST.csv', obs_long_lst, lst_hrs_temp, lst_rad_temp, FORMAT='(A10)',DELIMITER=',',skipline=1
    readcol, '/Users/nabarry/MWA/FHD/obs_list/Long_Run_LST.csv', obs_long_lst_temp, lst_hrs, lst_rad, FORMAT='(F10.15)',DELIMITER=',',skipline=1

    for obs_i=0, N_elements(obs_long_lst)-1 do begin
      obs_in = where(strmatch(btl,obs_long_lst[obs_i]) EQ 1,n_count)
      if n_count GT 0 then begin
        if n_elements(obs_index) GT 0 then obs_index=[obs_index,obs_in] else obs_index=obs_in
      endif
    endfor
    for obs_i=0, N_elements(occupancy4)-1 do begin
      obs_in = where(strmatch(btl[obs_index],occupancy4[obs_i]) EQ 1,n_count)
      if n_count GT 0 then begin
        if n_elements(obs_index_2) GT 0 then obs_index_2=[obs_index_2,obs_in] else obs_index_2=obs_in
      endif
    endfor
    for obs_i=0, N_elements(tv50)-1 do begin
      obs_in = where(strmatch(btl[obs_index],tv50[obs_i]) EQ 1,n_count)
      if n_count GT 0 then begin
        if n_elements(obs_index_tv50) GT 0 then obs_index_tv50=[obs_index_tv50,obs_in] else obs_index_tv50=obs_in
      endif
    endfor
    for obs_i=0, N_elements(tv100)-1 do begin
      obs_in = where(strmatch(btl[obs_index],tv100[obs_i]) EQ 1,n_count)
      if n_count GT 0 then begin
        if n_elements(obs_index_tv100) GT 0 then obs_index_tv100=[obs_index_tv100,obs_in] else obs_index_tv100=obs_in
      endif
    endfor
    for obs_i=0, N_elements(tv200)-1 do begin
      obs_in = where(strmatch(btl[obs_index],tv200[obs_i]) EQ 1,n_count)
      if n_count GT 0 then begin
        if n_elements(obs_index_tv200) GT 0 then obs_index_tv200=[obs_index_tv200,obs_in] else obs_index_tv200=obs_in
      endif
    endfor
    for obs_i=0, N_elements(alltv)-1 do begin
      obs_in = where(strmatch(btl[obs_index],alltv[obs_i]) EQ 1,n_count)
      if n_count GT 0 then begin
        if n_elements(obs_index_alltv) GT 0 then obs_index_alltv=[obs_index_alltv,obs_in] else obs_index_alltv=obs_in
      endif
    endfor


    julian_date_jan = 2451545.0
    ;btl_jd_round = round(btl_jd - 2451545.0)
    btl_jd_round = (btl_jd - 2451545.0)

    ytitle = 'Julian Date (+2451545.0)'
    xtitle = 'LST (hours)'
    title = 'SSINS cut'

    steps = 5
    redVector = REPLICATE(255, steps)
    blueVector = REPLICATE(0, steps)
    scaleFactor = FINDGEN(steps) / (steps - 1)
    beginNum = 255
    endNum = 0
    greenVector = beginNum + (endNum - beginNum) * scaleFactor

    greenVector[0]=223

    color_num = [8,10,12,14,16]

    for n_i=0, 4 do $
      TVLCT, redVector[n_i], greenVector[n_i], blueVector[n_i], color_num[n_i]

    color_array = [8B,10B,12B,14B,16B]

    pdf_plot=1
    if keyword_set(pdf_plot) then cgPS_Open,'/Users/nabarry/MWA/ssins_cut.pdf',/quiet,/nomatch
    cgplot, lst_hrs, btl_jd_round[obs_index], psym = 'Open Circle', linestyle=6, aspect = .7, yrange=[4978,5085], symsize=.8, xtitle=xtitle, ytitle=ytitle, title=title, charsize=1.1
    cgoplot, lst_hrs[obs_index_2], btl_jd_round[obs_index[obs_index_2]], psym = 'Filled Circle', linestyle=6, color= 'black', symsize=.8
    cgoplot, lst_hrs[obs_index_alltv], btl_jd_round[obs_index[obs_index_alltv]], psym = 'Filled Circle', linestyle=6, color = color_array[0], symsize=.8
    cgoplot, lst_hrs[obs_index_tv200], btl_jd_round[obs_index[obs_index_tv200]], psym = 'Filled Circle', linestyle=6, color =  color_array[1], symsize=.8
    cgoplot, lst_hrs[obs_index_tv100], btl_jd_round[obs_index[obs_index_tv100]], psym = 'Filled Circle', linestyle=6, color= color_array[2], symsize=.8
    cgoplot, lst_hrs[obs_index_tv50], btl_jd_round[obs_index[obs_index_tv50]], psym = 'Filled Circle', linestyle=6, color= color_array[3], symsize=.8
    cglegend,title=['>40% occupied','TV top 50','TV top 100','TV top 200','all TV'],$
      color=['black',color_array[0],color_array[1],color_array[2],color_array[3]], location=[.92,.7], charsize=.7, linestyle=[6,6,6,6,6], psyms=[16,16,16,16,16],length=0
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage

    stop
  endif




  limit_compare=1
  if keyword_set(limit_compare) then begin

    n_papers=10
    papers_array = PTRARR(n_papers,/allocate)
    name_array = STRARR(n_papers)

    ;Organized by paper: mK value, center of k bin/range, and redshift
    ali_2015 = transpose([22.4^2.,.325,8.4])
    *papers_array[0] = ali_2015
    name_array[0] = 'Ali, 2015'

    paciga_2013 = transpose([6.15E4,.5,8.6])
    *papers_array[1] = paciga_2013
    name_array[1] = 'Paciga, 2013'

    parsons_2014 = transpose([41^2.,.27,7.7])
    *papers_array[2] = parsons_2014
    name_array[2] = 'Parsons, 2014'

    patil_2017 = FLTARR(15,3)
    patil_2017[*,0] = [131.5^2.,242.1^2.,220.9^2.,337.4^2.,407.7^2.,86.4^2.,144.2^2.,184.7^2.,$
      296.1^2.,342.0^2.,79.6^2.,108.8^2.,148.6^2.,224.0^2.,366.1^2.]
    patil_2017[*,1] = [.053,.067,.083,.103,.128,.053,.067,.083,.103,.128,.053,.067,.083,.103,.128]
    patil_2017[*,2] = [8.3,8.3,8.3,8.3,8.3,9.15,9.15,9.15,9.15,9.15,10.1,10.1,10.1,10.1,10.1]
    patil_2017_short = FLTARR(6,3)
    patil_2017_short[*,0] = [131.5^2.,220.9^2.,86.4^2.,184.7^2.,79.6^2.,148.6^2.]
    patil_2017_short[*,1] = [.053,.083,.053,.083,.053,.083]
    patil_2017_short[*,2] = [8.3,8.3,9.15,9.15,10.1,10.1]
    patil_2017_extra_short = FLTARR(2,3)
    patil_2017_extra_short[*,0] = [79.6^2.,148.6^2.]
    patil_2017_extra_short[*,1] = [.053,.083]
    patil_2017_extra_short[*,2] = [10.1,10.1]
    *papers_array[3] = patil_2017_short
    name_array[3] = 'Patil, 2017'

    jacobs_2015 = FLTARR(4,3)
    jacobs_2015[*,0] = [1.64E4,6.9E3,6.9E3,2.4E3]
    jacobs_2015[*,1] = [.2,.1,.2,.2]
    jacobs_2015[*,2] = [10.29,8.54,7.94,7.55]
    *papers_array[4] = jacobs_2015
    name_array[4] = 'Jacobs, 2015'

    dillon_2014 = FLTARR(11,3)
    dillon_2014[*,0] = [2.6E7,1.16E6,8.64E5,6.7E5,1.3E6,1.28E7,5.26E7,5.67E8,4.58E6,2.93E8,6.92E8]
    dillon_2014[*,1] = [.058,.06,.063,.065,.0678,.0712,.0715,.149,.15,.15,.089]
    dillon_2014[*,2] = [11.68,10.868,10.153,9.518,8.444,7.985,7.57,7.1896,6.84,6.52,6.23]
    *papers_array[5] = dillon_2014
    name_array[5] = 'Dillon, 2014'

    dillon_2015 = FLTARR(3,3)
    dillon_2015[*,0] = [3.8E4,3.69E4,4.67E4]
    dillon_2015[*,1] = [.18,.18,.16]
    dillon_2015[*,2] = [6.4,6.8,7.25]
    *papers_array[6] = dillon_2015
    name_array[6] = 'Dillon, 2015'

    beardsley_2016 = FLTARR(9,3)
    beardsley_2016[*,0] = [3.67E4,2.70E4,3.56E4,3.02E4,4.70E4,3.22E4,3.2E4,2.6E4,2.5E4]
    beardsley_2016[*,1] = [.231,.27,.24,.24,.20,.24,.16,.14,.14]
    beardsley_2016[*,2] = [7.1,7.1,6.8,6.8,6.5,6.5,7.1,6.8,6.5]
    beardsley_2016_short = FLTARR(6,3)
    beardsley_2016_short[*,0] = [2.70E4,3.02E4,3.22E4,3.2E4,2.6E4,2.5E4]
    beardsley_2016_short[*,1] = [.27,.24,.24,.16,.14,.14]
    beardsley_2016_short[*,2] = [7.1,6.8,6.5,7.1,6.8,6.5]
    beardsley_2016_extra_short = FLTARR(4,3)
    beardsley_2016_extra_short[*,0] = [2.70E4,3.02E4,3.2E4,2.5E4]
    beardsley_2016_extra_short[*,1] = [.27,.24,.16,.14]
    beardsley_2016_extra_short[*,2] = [7.1,6.8,7.1,6.5]
    *papers_array[7] = beardsley_2016_short
    name_array[7] = 'Beardsley, 2016'

    ;barry_2019 = FLTARR(1,3)
    ;barry_2019[0,0] = 2.2E3
    ;barry_2019[0,1] = 0.23
    ;barry_2019[0,2]=7.0
    barry_2019 = FLTARR(6,3)
    barry_2019[*,1]=[0.17426223,0.20330593,0.23234963,0.26139334,0.29043704,0.31948075]
    barry_2019[*,0]=[9176.4701,3273.5603,5681.5812,9564.6415,15549.955,20607.678]
    barry_2019[*,2]=[7,7,7,7,7,7]
    *papers_array[8] = barry_2019
    name_array[8] = 'Barry, 2019'


    ;paper_collab = FLTARR(3,3)
    ;paper_collab[*,0] = [41^2.,22.4^2.,2.4E3]
    ;paper_collab[*,1] = [.27,.325,.2]
    ;paper_collab[*,2] = [7.7,8.4,7.55]
    ;*papers_array[4] = paper_collab
    ;name_array[4] = 'Paper collab.'
    paper_collab = FLTARR(6,3)
    paper_collab[*,0] = [8.09e4,1.28e5,3.61e4,1.46e5,3.8e6,2.41e6]
    paper_collab[*,1] = [0.39,0.53,0.31,0.36,0.334,0.33]
    paper_collab[*,2] = [7.49,8.13,8.37,8.68,9.93,10.88]
    *papers_array[4] = paper_collab
    name_array[4] = 'Kolopanis, 2019'

    phaseII = FLTARR(3,3)
    phaseII[*,0] = [9.28E3,5.59E3,2.39E3]
    phaseII[*,1] = [.25,.29,.59]
    phaseII[*,2] = [7.1,6.8,6.5]
    *papers_array[9] = phaseII
    name_array[9] = 'Li, 2019 (in review)'

    ;fiducial theory
    restore,'/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave'
    delta_eor=power*(k_centers^3.)/(2.*!pi^2.)
    hubble_param=0.71
    k_centers=k_centers/hubble_param

    brad_update=1
    if keyword_set(brad_update) then begin
      ;GalaxyParam_LF+NF+Tau
      ;Fourier mode (k [Mpc^-1]), 68th, 95th and 99th percentiles on the power spectrum in units of mK^2
      brad_k = [0.02094395,0.03183814,0.0477084,0.06481693,0.0846983,0.1122589,0.1512655,0.2037167,0.2749778,0.3717481,0.5018322,0.6776517,0.9149306,1.235268,1.6676,2.207298,2.79783]
      brad_68 = [2.06836030818,9.33650238046,13.5354885655,27.9094808815,31.8379169811,35.6549002097,$
        37.7766299375,34.7228577191,22.5256589937,19.1373452916,19.0715957566,18.8232315706,18.1639861695,$
        18.6493671092,21.4087292929,28.7580668946,40.4441293602]
      brad_95 = [3.00633765723,12.8642551767,18.6682120582,38.322192403,43.5651253669,49.1075159329,$
        51.9213041607,48.790639303,31.3894142399,26.2910138229,26.8292876974,26.5396841581,26.0301170622,$
        27.1653881665,31.5002616162,42.7424864461,60.1632270707]
      brad_99 = [3.277767,13.62931,19.91,40.66932,46.2304,52.1699,55.14123,52.86583,34.2097,28.67557,29.69817,29.57222,$
        29.46968,30.93632,35.9361,48.86067,68.80418]

      ;k (Mpc^-1) and fiducial
      brad_k_2 = [2.094395e-02,3.183814e-02,4.770840e-02,6.481693e-02,8.469830e-02,1.122589e-01, 1.512655e-01,2.037167e-01,2.749778e-01,3.717481e-01,$
        5.018322e-01,6.776517e-01,9.149306e-01,1.235268e+00,1.667600e+00,2.207298e+00,2.797836e+00]
      brad_fiducial = [3.148450e-01,1.371862e+00,2.112182e+00,4.716543e+00,6.822401e+00,1.036046e+01,1.477729e+01,9.501991e+00,8.793401e+00,7.287165e+00,$
        6.430653e+00,5.709567e+00,5.390083e+00,5.408442e+00,6.280324e+00,8.959497e+00,1.353618e+01]

      brad_k_2 = brad_k_2/hubble_param
      brad_k = brad_k/hubble_param

      n_k=n_elements(brad_k)
      brad_k=(brad_k[1:(n_k-1)]+brad_k[0:(n_k-2)])/2.
      ;n_k=n_elements(brad_k_2)
      ;brad_k_2=(brad_k_2[1:(n_k-1)]+brad_k_2[0:(n_k-2)])/2.

      k_centers = brad_k_2
      delta_eor=brad_fiducial


    endif

    sym_array = ['Filled Diamond','Filled Bowtie','Filled Up Triangle','Filled Lower Half Circle','Filled Up Triangle','Filled Circle','Filled Laying Bar','Filled Star','Filled Square','Filled Diamond']
    sym_num = [14,24,17,40,17,16,28,46,15,24]

    ;for n_i=0, 4 do $
    ;  TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    ;color_array = [10B,12B,14B,16B,18B,20B]
    color_array = ['purple', 'teal', 'blue', 'purple','orange', 'teal', 'blue', 'firebrick','green','gold']
    xtitle_use='k (h Mpc$\up-1$)'
    ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
    ;plot_2d=1
    if keyword_set(plot_2d) then begin
      nick_version = 1
      if keyword_set(nick_version) then begin
        cgPS_Open,'/Users/nabarry/MWA/data/limits_nickversion.pdf',/quiet,/nomatch
        sym_array = ['Filled Diamond','Filled Bowtie','Filled Up Triangle','Filled Lower Half Circle','Filled Up Triangle','Filled Circle','Filled Laying Bar','Filled Star','Filled Square']
        paper_ind = [-1,0,-1,1,2,3,4,5,6,7,8,9]
        redshift_ranges=[6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5] ;11 redshift ranges
        redshift_colors = ['green','teal','blue','purple','hot pink','firebrick','orange red','orange','gold','yellow']

        cgplot, (*papers_array[0])[0,1],(*papers_array[0])[0,0], /xlog,/ylog, xrange=[.04,1.0], yrange=[1.,1E7], aspect=.5, $
          YGridStyle=1,YTicklen=1.0, xtitle=xtitle_use, ytitle = ytitle_use,/NoData

        for papers_i=0, n_papers-1 do begin
          arr = *papers_array[papers_i]

          if (papers_i EQ 0) OR (papers_i EQ 2) then continue

          for arr_i=0,N_elements(arr[*,0])-1 do begin

            ind = where(round(arr[arr_i,2]/.5)*.5 EQ redshift_ranges,n_count)
            cgoplot, arr[arr_i,1],arr[arr_i,0],psym=sym_array[paper_ind[papers_i]],thick=3, color=redshift_colors[ind], /xlog,/ylog, $
              symsize=2.0,charsize =1.

          endfor
        endfor

        cglegend, Title=['z=6.25-6.75','z=6.75-7.25','z=7.25-7.75','z=7.75-8.25'],$
           color=redshift_colors[0:3], location=[.15,.4], psym=[15,15,15,15],charsize=1,Length=0,symsize=1.5
         cglegend, Title=['z=8.25-8.75','z=8.75-9.25','z=9.25-9.75','z=9.75-10.25'],$
           color=redshift_colors[4:7], location=[.35,.4], psym=[15,15,15,15],charsize=1,Length=0,symsize=1.5
         cglegend, Title=['z=10.25-10.75','z=10.75-11.25'],$
           color=redshift_colors[8:9], location=[.55,.4], psym=[15,15],charsize=1,Length=0,symsize=1.5
        cglegend, Title=name_array[[1,3]], color=['grey','grey'], location=[.07,.12], psym=sym_num[[0,1]], charsize=1,Length=0,symsize=1.5
        cglegend, Title=name_array[[5,6]], color=['grey','grey'], location=[.28,.12], psym=sym_num[[3,4]], charsize=1,Length=0,symsize=1.5
        cglegend, Title=name_array[[7,4]], color=['grey','grey'], location=[.5,.12], psym=sym_num[[5,2]], charsize=1,Length=0,symsize=1.5
        cglegend, Title=name_array[[9,8]], color=['grey','grey'], location=[.7,.12], psym=sym_num[[7,6]], charsize=1,Length=0,symsize=1.5

        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        stop
      endif else begin
        cgPS_Open,'/Users/nabarry/MWA/data/limits_full4.pdf',/quiet,/nomatch
        for papers_i=0, n_papers-1 do begin
          k_arr = *papers_array[papers_i]

          if (papers_i EQ 0) OR (papers_i EQ 2) then continue

          if papers_i EQ 1 then cgplot, k_arr[*,1],k_arr[*,0],psym=sym_array[papers_i],thick=3, color=color_array[papers_i], /xlog,/ylog, xrange=[.04,.65], yrange=[1.,1E7], aspect=.5, $
            YGridStyle=1,symsize=1.5,YTicklen=1.0, xtitle=xtitle_use, ytitle = ytitle_use,charsize =1.
          cgoplot, k_arr[*,1],k_arr[*,0],psym=sym_array[papers_i],thick=3, color=color_array[papers_i],symsize=2.0
          PLOTSYM, 1 ,2, thick=3
          cgoplot, k_arr[*,1],k_arr[*,0],psym=8,thick=3, color=color_array[papers_i],symsize=1.5
          ;cgoplot, k_centers,delta_eor,/xlog,/ylog,color='brown',thick=4
          ;cgoplot, brad_k,brad_95,/xlog,/ylog,color='brown',thick=4,psym=10
          ;        delta_k=brad_k[1]/2.
          ;        error_low=10
          ;       for k_i=1, n_elements(brad_k)-2 do $
          ;         cgColorFill, [brad_k[k_i]-(brad_k[k_i-1]+brad_k[k_i])/2., brad_k[k_i]+(brad_k[k_i+1]+brad_k[k_i])/2., brad_k[k_i]+(brad_k[k_i+1]+brad_k[k_i])/2.,brad_k[k_i]-(brad_k[k_i-1]+brad_k[k_i])/2.], $
          ;         [brad_95[k_i], brad_95[k_i], 1.,1.], $
          ;         Color='grey',/checkforfinite
          cgoplot, brad_k,brad_95,/xlog,/ylog,color='brown',thick=4,psym=10


          cglegend, Title=name_array[[1,3]], color=color_array[[1,3]], location=[.02,.12], psym=sym_num[[1,3]], charsize=1.,Length=0,symsize=1.5
          cglegend, Title=name_array[[5,6]], color=color_array[[5,6]], location=[.20,.12], psym=sym_num[[5,6]], charsize=1.,Length=0,symsize=1.5
          cglegend, Title=name_array[[7,4]], color=color_array[[7,4]], location=[.32,.12], psym=sym_num[[7,4]], charsize=1.,Length=0,symsize=1.5
          cglegend, Title=name_array[[4,8]], color=color_array[[4,8]], location=[.55,.12], psym=sym_num[[4,8]], charsize=1.,Length=0,symsize=1.5

          for n_i=0, 1 do $
            TVLCT, 197, 120, 62, 10
          cgoplot, k_centers,delta_eor,/xlog,/ylog,color=10B,thick=4


        endfor
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        stop
      endelse
    endif else begin

      for rotz_i=0,30 do begin

        cgPS_Open,"/Users/nabarry/MWA/data/rotate_2019/limits_all_fold_c_"+string(strtrim(rotz_i,2),FORMAT='(I03)')+".png",FONT=1, Charsize=3.0;,/quiet,/nomatch

        ; Set the 3D coordinate space with axes.
        ;x = k
        ;y = z
        ;z = mk^2
        cgSurf, DIST(5), /NODATA, /SAVE, XRANGE=[.04,0.8], $
          YRANGE=[6,12], ZRANGE=[100, 1e6], XSTYLE=1, $
          YSTYLE=1, ZSTYLE=1, CHARSIZE=2.0, $
          POSITION=[0.1, 0.1, 0.95, 0.95, 0.1, 0.95], $
          XTICKLEN=1, YTICKLEN=1, XGRIDSTYLE=1, YGRIDSTYLE=1, $
          xtitle=xtitle_use, ytitle='redshift',ztitle=ytitle_use,title=title, $
          /zlog, /xlog, rotz=30. - float(rotz_i)/2.
        cgAXIS, XAXIS=1, /T3D, CHARSIZE=2.0
        cgAXIS, YAXIS=1, /T3D, CHARSIZE=2.0

        for papers_i=0, n_papers-1 do begin
          if (papers_i EQ 0) OR (papers_i EQ 2) then continue
          arr = *papers_array[papers_i]
          z = arr[*,0]
          x = arr[*,1]
          y = arr[*,2]
          FOR j=0,N_elements(z)-1 DO begin
            if z[j] GT 1.0E6 then continue
            cgPlotS, x[j], y[j], z[j], COLOR=color_array[papers_i], PSYM=sym_array[papers_i], SYMSIZE=2.5, /T3D
          ENDFOR
          FOR j=0,N_elements(z)-1 DO begin
            if z[j] GT 1.0E6 then continue
            cgPlotS, [x[j], x[j]], [y[j], y[j]], [100, z[j]], COLOR=color_array[papers_i], /T3D
          ENDFOR
            

        endfor

        cglegend, Title=name_array[[1,3]], color=color_array[[1,3]], location=[.02,.05], psym=sym_num[[1,3]], charsize=1.25,Length=0,symsize=1.5
        cglegend, Title=name_array[[5,6]], color=color_array[[5,6]], location=[.20,.05], psym=sym_num[[5,6]], charsize=1.25,Length=0,symsize=1.5
        cglegend, Title=name_array[[7,4]], color=color_array[[7,4]], location=[.42,.05], psym=sym_num[[7,4]], charsize=1.25,Length=0,symsize=1.5
        cglegend, Title=name_array[[8]], color=color_array[[8]], location=[.69,.05], psym=sym_num[[8]], charsize=1.25,Length=0,symsize=1.5

        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent;,/nomessage
      endfor
      stop
    endelse
  endif

  ;previous limit #Limit: 27845.762 mK^2, at k = 0.20330593 h Mpc^-1 (wedge_cut_plus_res_cut ch9-126_ res  yy) (published 2.7)
  ;1.28 for flux scale drop
  ;        dsigma = dsigma * 3.4 noise level scaling
  ;#Limit: 5561.384 mK^2, at k = 0.23234963 h Mpc^-1 (btl_noalltv_noocc4 ch9-126_ res fullimg_ yy) ;scaled with estimated Adam noise level
  ;#Limit: 2158.5883 mK^2, at k = 0.23234963 h Mpc^-1 (btl_noalltv_noocc4 ch9-126_ res fullimg_ yy)

  ;checking_sim_plots=1
  if keyword_set(checking_sim_plots) then begin
    filename_24 = '/Users/nabarry/MWA/data/fhd_nb_sim_BH2_BH2_justEoR/' ;left over foregounds, perfect
    filename_3 = '/Users/nabarry/MWA/data/fhd_nb_sim_BH2_BH2_justEoR/ps_doublearr_stdft/' ;just EoR


    obsname='zenith'
    ;obsname_1='1061319472'
    obsname_3='zenith'

    cube_type='res'
    power_24 = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_fullimg_'+cube_type+'_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_fullimg_'+cube_type+'_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_24 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    weights = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_fullimg_'+cube_type+'_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave','weights')
    dsigma = (k_24^3.)/(2.*!pi^2.)/sqrt(weights)

    cube_type='res'
    power_3 = getvar_savefile(filename_3 + 'Combined_obs_'+obsname_3+'_cubeXX__even_odd_joint_fullimg_'+cube_type+'_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_3 + 'Combined_obs_'+obsname_3+'_cubeXX__even_odd_joint_fullimg_'+cube_type+'_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_3 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.

    power_eor = getvar_savefile('/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave','power')
    k_centers_eor = getvar_savefile('/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave','k_centers')

    delta=1
    if ~keyword_set(delta) then ytitle = 'P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)' else ytitle = '$\Delta$$\up2$ (mK$\up2$)'
    delta_24=power_24*(k_24^3.)/(2.*!pi^2.)
    delta_3=power_3*(k_3^3.)/(2.*!pi^2.)
    delta_eor=power_eor*(k_centers_eor^3.)/(2.*!pi^2.)

    hubble_param=0.71
    k_3=k_3/hubble_param
    k_24=k_24/hubble_param
    k_centers_eor=k_centers_eor/hubble_param

    cgPS_Open,'/Users/nabarry/MWA/data/EoR_recovered_limit_noise_fold.pdf',/quiet,/nomatch
    cgplot, k_24, (delta_24),/ylog,/xlog,psym=10,xrange=[.05,1.2],charsize=1.3, xtitle = 'k (!8h!x / Mpc)',yrange=[10^0., 10^3.],$
      ytitle=ytitle,title = 'In-situ simulation', color='forest green',thick=7, aspect=.5;thick=3
    cgoplot, k_24, dsigma, psym=10, linestyle =2, color='forest green', thick=7
    cgoplot, k_3, delta_3,/ylog,/xlog,psym=10, color='orange',thick=11, linestyle=0
    cgoplot, k_centers_eor, delta_eor,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='purple',thick=7
    cglegend, title=['input EoR', 'EoR in pipeline','foregrounds, EoR'],$
      color=['purple','orange','forest green'],location=[.45,.7], charsize=1.1,thick=7, linestyle=[0,0,0];,/background,/box
    ;cglegend, title=['input Gaussian EoR','recovered Gaussian EoR','no filter','Tukey filter','Blackman-Harris filter','Blackman-Harris filter, 100$\lambda$ extent'],$
    ;  color=['purple','dark green','navy','firebrick','orange','turquoise'],location=[.42,.85], charsize=1.1


    cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif

  ;brad_limit_plot=1
  if keyword_set(brad_limit_plot) then begin

    cut = 'btl_noalltv_noocc4'
    cut2 = 'sub_cubes'
    chans='ch9-126_'
    win='fullimg_'
    spec_window='bh_'
    ;basefile='/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_fftchange_spectralwindow/Combined_obs_'+cut+'_cubeXX__even_odd_joint_'+win+chans
    ;basefile2='/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_fftchange_spectralwindow/Combined_obs_'+cut2+'_cubeXX__even_odd_joint_'+win+chans
    basefile='/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_doublearr_stdft/Combined_obs_'+cut+'_cubeXX__even_odd_joint_'+win+chans
    basefile2='/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_doublearr_stdft/Combined_obs_'+cut2+'_cubeXX__even_odd_joint_'+win+chans
    cubes=['res']
    pols=['yy','xx']
    averemove='_averemove_'
    kperp_end='80'
    kperp_start='18'
    dencorr = 'dencorr_'
    endfile=averemove+'sw'+spec_window+dencorr + 'no_120deg_wedge_kperplambda'+kperp_start+'-'+kperp_end+'_kpar0.15-200_1dkpower.idlsave'

    ;limit_percent=0.9772 ;2 sigma
    limit_percent = 0.97725
    hubble_param=0.71

    n_cubes=1
    color_num = [10,12,14,16,18]
    rgbcolors = [[137,117,202],[114,166,89],[73, 142, 217],[197,120,62],[165,42,42]]
    ;73, 142, 217 little brighter blue, 54,107,163 little darker blue
    ;rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101]];,[198,118,63]]

    ranges=PTRARR(n_cubes,/allocate)
    ;*ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,28,38,39,40,41,42,43,54,53,55]+4

    *ranges[0]=[0,9,10,11,12,13,23,24,25,26,27,38,39,40,41,42,52,53,54,53,55]+4

    for n_i=0, N_elements(color_num)-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]

    color_array = [10B,12B,14B,16B,18B]

    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 5 else thickness = 3

    for j=0,0 do begin

      file=basefile+cubes+'_'+pols[j]+endfile
      ;file2=basefile2+cubes+'_'+pols[j]+endfile
      restore,file
      ;k_edges2 = getvar_savefile(file2, 'k_edges')
      ; find the best limit in this file
      n_k=n_elements(k_edges)
      ;n_k2=n_elements(k_edges2)
      k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
      ;k2=(k_edges2[1:(n_k2-1)]+k_edges2[0:(n_k2-2)])/2.
      ;power2 = getvar_savefile(file2,'power')
      delta=power*(k^3.)/(2.*!pi^2.)
      ;delta2=power2*(k2^3.)/(2.*!pi^2.)

      ;weights2 = getvar_savefile(file2,'weights')
      dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
      dsigma[0] = !Values.F_INFINITY
      ;ab_scaled=1
      if keyword_set(ab_scaled) then begin
        print, "Artificial scaling applied"
        dsigma = dsigma * 3.4
      endif
      ;dsigma2=(k2^3.)/(2.*!pi^2.)/sqrt(weights2)
      ;dsigma2[0] = !Values.F_INFINITY

      limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf((delta)/dsigma/sqrt(2)))*sqrt(2))+delta
      ;limits2=dsigma2*(inverf(limit_percent-(1.-limit_percent)*erf((delta2)/dsigma2/sqrt(2)))*sqrt(2))+delta2
      limits_abs=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(abs(delta)/dsigma/sqrt(2)))*sqrt(2))+abs(delta)
      lim=min(limits,ind)
      header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut+' '+chans+' '+cubes+' '+win+' '+pols[j]+')'
      print,header

      zeros = where(delta EQ 0,n_count)
      ;zeros2 = where(delta2 EQ 0,n_count)
      delta[zeros] = !Values.F_INFINITY
      delta[0] = !Values.F_INFINITY
      ;delta2[zeros2] = !Values.F_INFINITY
      ;delta2[0] = !Values.F_INFINITY
      delta[*ranges[0]] = !Values.F_INFINITY
      ;delta2[*ranges[0]] = !Values.F_INFINITY
      dsigma[0] = !Values.F_INFINITY
      dsigma[*ranges[0]] = !Values.F_INFINITY
      ;dsigma2[0] = !Values.F_INFINITY
      ;dsigma2[*ranges[0]] = !Values.F_INFINITY
      limits[0] = !Values.F_INFINITY
      limits[*ranges[0]] = !Values.F_INFINITY
      ;limits2[0] = !Values.F_INFINITY
      ;limits2[*ranges[0]] = !Values.F_INFINITY
      k=k/hubble_param
      ;k2=k2/hubble_param

      brad_update=1
      if keyword_set(brad_update) then begin
        ;GalaxyParam_LF+NF+Tau
        ;Fourier mode (k [Mpc^-1]), 68th, 95th and 99th percentiles on the power spectrum in units of mK^2
        brad_k = [0.02094395,0.03183814,0.0477084,0.06481693,0.0846983,0.1122589,0.1512655,0.2037167,0.2749778,0.3717481,0.5018322,0.6776517,0.9149306,1.235268,1.6676,2.207298,2.79783]
        brad_68 = [2.06836030818,9.33650238046,13.5354885655,27.9094808815,31.8379169811,35.6549002097,$
          37.7766299375,34.7228577191,22.5256589937,19.1373452916,19.0715957566,18.8232315706,18.1639861695,$
          18.6493671092,21.4087292929,28.7580668946,40.4441293602]
        brad_95 = [3.00633765723,12.8642551767,18.6682120582,38.322192403,43.5651253669,49.1075159329,$
          51.9213041607,48.790639303,31.3894142399,26.2910138229,26.8292876974,26.5396841581,26.0301170622,$
          27.1653881665,31.5002616162,42.7424864461,60.1632270707]
        brad_99 = [3.277767,13.62931,19.91,40.66932,46.2304,52.1699,55.14123,52.86583,34.2097,28.67557,29.69817,29.57222,$
          29.46968,30.93632,35.9361,48.86067,68.80418]

        ;k (Mpc^-1) and fiducial
        brad_k_2 = [2.094395e-02,3.183814e-02,4.770840e-02,6.481693e-02,8.469830e-02,1.122589e-01, 1.512655e-01,2.037167e-01,2.749778e-01,3.717481e-01,$
          5.018322e-01,6.776517e-01,9.149306e-01,1.235268e+00,1.667600e+00,2.207298e+00,2.797836e+00]
        brad_fiducial = [3.148450e-01,1.371862e+00,2.112182e+00,4.716543e+00,6.822401e+00,1.036046e+01,1.477729e+01,9.501991e+00,8.793401e+00,7.287165e+00,$
          6.430653e+00,5.709567e+00,5.390083e+00,5.408442e+00,6.280324e+00,8.959497e+00,1.353618e+01]

        k_centers = brad_k_2/hubble_param
        brad_k = brad_k/hubble_param
        n_k2=n_elements(brad_k)
        brad_k=(brad_k[1:(n_k2-1)]+brad_k[0:(n_k2-2)])/2.
        ;k_centers = brad_k/hubble_param
        delta_eor=brad_fiducial
      endif

      if (j EQ 0) and keyword_set(pdf_plot) then cgPS_Open,'/Users/nabarry/MWA/data/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/ps_rtf_cupdate/brad_limit_update_fid.pdf',/quiet,/nomatch


      if (j EQ 0) then begin
        ;position_use=position1
        ;xtitle_use=''
        xtitle_use='k (h Mpc$\up-1$)'
        ytitle_use='$\Delta$$\up2$ (mK$\up2$)'
        XTICKFORMAT_use="(I0)"
        YTICKFORMAT_use="(I0)"
      endif

      color_use=color_array[2]
      low_range=10^0.
      low_x_range=.14

      if (j EQ 0) then begin
        cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^0.,10^8.],ytitle=ytitle_use, $
          xtitle=xtitle_use, charsize =.9, color=color_use,thick=thickness,/noerase, title= 'z=7',aspect=0.5
        cglegend,title=['2$\sigma$ upper limit','noise level','fiducial theory','95% confidence'],$
          color=[color_array[2],color_array[2],color_array[3],color_array[4]], location=[.15,.72], charsize=.9, thick=thickness,linestyle=[0,2,0,0],length=.03,vspace=1
      endif
      if (j EQ 1) then begin
        cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^1.,10^8.],ytitle=ytitle_use, $
          xtitle=xtitle_use, charsize =.8, color=color_use,thick=thickness,position=position_use,XTICKFORMAT="(A1)",YTICKFORMAT=YTICKFORMAT_use,/noerase
        XYOuts, position_use[0]+0.2, position_use[3]+.02, 'N-S, z=7', /Normal, Alignment=0.5, Charsize=.8
        cglegend,title=['measured power','fiducial theory','95% confidence'],$
          color=['black','brown',color_array[3]], location=[.1,.9], charsize=.7, thick=thickness,linestyle=[0,0],length=.03,vspace=1.1
      endif
      ;cgoplot, k2,abs(limits2),/xlog,/ylog,psym=10,color=color_array[1],thick=4,position=position_use, linestyle=3

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

      cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=color_use,thick=thickness,position=position_use
      cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position_use
      ;cgoplot, k,delta,/xlog,/ylog,color='black',psym=10,thick=thickness-1,position=position_use
      cgoplot, k_centers,delta_eor,/xlog,/ylog,color=color_array[3],thick=thickness-1,position=position_use
      cgoplot, brad_k,brad_95,/xlog,/ylog,color=color_array[4],thick=thickness-1,psym=10


      if (j EQ 0) then begin

        ;        cgplot, k2,abs(limits2),/xlog,/ylog,psym=10, xrange=[.14,1.9], c,ytitle=ytitle_use, YMINOR=1, $
        ;          xtitle='k (h Mpc$\up-1$)', charsize =.8, color=color_array[1],thick=thickness,position=position2,XTICKFORMAT=XTICKFORMAT_use,/noerase
        ;          cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=color_use,thick=thickness,position=position2
        ;          cgoplot, k2,dsigma2,/xlog,/ylog,psym=10,linestyle=2, color=color_array[1],thick=thickness,position=position2
        ;          cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position2

        ;        cglegend,title=['noise level'],$
        ;          color=[color_array[1]], location=[.52,.58], charsize=.7, thick=thickness,linestyle=[2],length=.03,vspace=1
      endif
      if (j EQ 1) then begin
        ;        cgplot, k2,abs(limits2),/xlog,/ylog,psym=10, xrange=[.14,1.9], yrange=[10^2.,10^8.],ytitle=ytitle_use, YMINOR=1, $
        ;          xtitle='k (h Mpc$\up-1$)', charsize =.8, color=color_array[1],thick=thickness,position=position5,XTICKFORMAT=XTICKFORMAT_use,YTICKFORMAT=YTICKFORMAT_use,/noerase
        ;          cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=color_use,thick=thickness,position=position5
        ;          cgoplot, k2,dsigma2,/xlog,/ylog,psym=10,linestyle=2, color=color_array[1],thick=thickness,position=position5
        ;          cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_use,thick=thickness,position=position5

        ;        cglegend,title=['pre-SSINS 2$\sigma$ upper limit'],$
        ;          color=[color_array[1]], location=[.12,.58], charsize=.7, thick=thickness,linestyle=[0],length=.03,vspace=1.1
      endif

    endfor
    if keyword_set(pdf_plot) then cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage

    stop
  endif

end