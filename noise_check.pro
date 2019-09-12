pro noise_check

  noise_check_plot=1
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

        sigma=2
        dsigma = dsigma * sigma
        title_arr = ['observation noise','AB2016 observation noise','updated observation noise']
        print, title_arr[cube_i]
        print,[transpose(k),transpose(dsigma)]

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

end