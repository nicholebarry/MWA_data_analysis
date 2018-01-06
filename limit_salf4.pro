pro limit_salf4

  ;coarse_harm_width = 3
  ;wedge_angles=120d
  n_cubes=4
  
  cut = STRARR(n_cubes)
  ;cut[2]='wedge_cut_plus_res_cut'
  ;cut[2]='Aug23_longrunstyle'
  ;cut[2] = 'beardsley_thesis_list_8shist'
  cut[3] = 'beardsley_thesis_list'
  cut[2] = 'beardsley_thesis_list'
  ;cut[1] = 'beardsley_thesis_list_firstthird'
  cut[1] = 'beardsley_thesis_list'
  ;cut[0] = 'beardsley_thesis_lis_noautocutt'
  cut[0] = 'beardsley_thesis_list'
  
  chans = strarr(n_cubes)
  chans[0]='ch0-95_' & band='low'
  ;chans='ch0-127_'
  chans[2]='ch96-191_' & band='high'
  chans[1]='ch48-143_' & band='mid'
  chans[3]='ch10-127_'
  ;chans=''
  
  win = strarr(n_cubes)
  win[3]='tk_'
  win[2]='tk_'
  win[1]='tk_'
  win[0]='tk_'
  
  spec_window='bh_'
  
  basefile = STRARR(n_cubes)
  ;basefile[2]=['/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps_jan2016/Combined_obs_'+cut[2]+'_cubeXX__even_odd_joint_'+win[2]+chans]
  ;basefile[2]=['/nfs/mwa-05/r1/EoRuvfits/EoR2013/fhd_nb_2013longrun_savedbp/ps_large/Combined_obs_'+cut[2]+'_cubeXX__even_odd_joint_tk_'+chans+'_']
  basefile[3]=['/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal/ps_redo/Combined_obs_'+cut[3]+'_cubeXX__even_odd_joint_'+win[3]+chans[3]]
  basefile[2]=['/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal/ps/Combined_obs_'+cut[2]+'_cubeXX__even_odd_joint_'+win[2]+chans[2]]
  ;basefile[1]=['/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_Aug2017_globalbp_w_cable_w_digjump/ps/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans]
  ;basefile[1]=['/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal/ps/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+'ch0-127_']
  basefile[1]=['/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/ps/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans[1]]
  basefile[0]=['/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal/ps/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans[0]]
  
  ranges=PTRARR(n_cubes,/allocate)
  *ranges[0]=[5,6,7,8,9,17,18,19,20,21,29,30,31,32,33,41,42]
  *ranges[1]=[5,6,7,8,9,17,18,19,20,21,29,30,31,32,33,41,42,43]
  *ranges[2]=[5,6,7,8,9,17,18,19,20,21,29,30,31,32,33,41,42,43]
  *ranges[3]=[6,7,8,9,10,11,21,22,23,24,25,26,36,37,38,39,40,41,50,51,52]
  
  cubes=['res']
  pols=['yy']
  
  endfile = STRARR(n_cubes)
  ;endfile[2]='_bh_dencorr_no_120deg_wedge_cbw3_kperplambda10-60_kpar0.15-200_1dkpower.idlsave'
  ;endfile[2]='_averemove_bh_dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'
  endfile[3]='_averemove_'+spec_window+'dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave'
  endfile[2]='_averemove_'+spec_window+'dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave'
  endfile[1]='_averemove_'+spec_window+'dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave'
  endfile[0]='_averemove_'+spec_window+'dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave'
  
  limit_percent=0.9772 ;2 sigma
  hubble_param=0.71
  
  ;color_array=['black','blue','green','purple']
  color_num = [10,12,14,16]
  ;rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
  rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101]];,[198,118,63]]
  
  for n_i=0, n_cubes-1 do $
    TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]
    
  color_array = [10B,12B,14B,16B]
  
  for j=0,N_elements(pols)-1 do begin
    for cube_i=0,n_cubes-1 do begin
    
      file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
      restore,file
      ; find the best limit in this file
      n_k=n_elements(k_edges)
      k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
      delta=power*(k^3.)/(2.*!pi^2.)
      dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
      ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
      limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
      lim=min(limits,ind)
      header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
      print,header
      
      ;remove first bins
      limits=limits[3:N_elements(limits)-1]
      limits[*ranges[cube_i]] = !Values.F_NAN
      k=k[3:N_elements(k)-1]/hubble_param
      
      ;if (cube_i EQ 0) then cgPS_Open,'/nfs/mwa-00/h1/nbarry/'+pols[j]+'_limit_comparison_salf_1.png',/quiet,/nomatch
      if cube_i EQ 3 then cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.15,1.7], yrange=[10^2.,10^7.],ytitle='$\Delta$$\up2$ (mK$\up2$)', $
        xtitle='k (h Mpc$\up-1$)', title='UW preliminary EoR limit', charsize =1.25, color=color_array[cube_i],thick=5; else $
      ;    cgoplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.08,1.5],color=color_array[cube_i],thick=3
      if cube_i EQ 3 then  cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_array[cube_i],thick=5
      
      
    endfor
    ;cglegend, title=['z=7.1','z=6.8','z=6.5','z=7'],color=color_array, location=[.2,.85], charsize=1, thick=3
    cglegend, title=['z=7 upper limit','1$\sigma$ noise'],color=[color_array[3],color_array[3]], location=[.18,.85], charsize=1.25, thick=5, linestyle=[0,2]
    ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    
    ;add_others=1
    if keyword_set(add_others) then begin
      n_papers=8
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
      *papers_array[7] = beardsley_2016_short
      name_array[7] = 'Beardsley, 2016'
      
      sym_array = ['Filled Diamond','Filled Bowtie','Filled Up Triangle','Filled Lower Half Circle','Filled Standing Bar','Filled Circle','Filled Laying Bar','Filled Star']
      sym_num = [14,24,17,40,26,16,28,46]
      
      color_num = [10,12,14,16,18,20]
      ;rgbcolors = [[9,14,181],[70,104,206],[52,115,145],[60,158,149],[60,158,102],[52,175,63]]
      ;rgbcolors = [[60,158,190],[60,158,160],[60,158,145],[60,158,120],[60,158,102],[60,158,63]]
      rgbcolors = [[45,31,198],[52,116,173],[40,183,183],[39,158,96],[47,158,39]]
      
      ;for n_i=0, 4 do $
      ;  TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]
        
      ;color_array = [10B,12B,14B,16B,18B,20B]
      color_array = ['purple', 'teal', 'blue', 'firebrick','purple', 'teal', 'blue', 'firebrick']
      
      for papers_i=0, n_papers-1 do begin
        k_arr = *papers_array[papers_i]
        
        ;inds1 = where(k_arr[*,2] GT 6,n_count1)
        inds = where((k_arr[*,2] LE 8.0) AND (k_arr[*,2] GE 6.0),n_count)
        ;inds3 = where((k_arr[*,2] LT 8.0) AND (k_arr[*,2] GE 7.0),n_count3)
        ;inds4 = where((k_arr[*,2] LT 9.0) AND (k_arr[*,2] GE 8.0),n_count4)
        ;inds5 = where((k_arr[*,2] LT 10.0) AND (k_arr[*,2] GE 9.0),n_count5)
        ;inds5 = where(k_arr[*,2] GE 9.0,n_count5)
        
        if n_count GT 0 then cgoplot, k_arr[inds,1],k_arr[inds,0],psym=sym_array[papers_i],thick=3, color=color_array[papers_i]
        ;if n_count1 GT 0 then cgoplot, k_arr[inds1,1],k_arr[inds1,0],psym=sym_array[papers_i],thick=3, color=color_array[0]
        ;if n_count2 GT 0 then cgoplot, k_arr[inds2,1],k_arr[inds2,0],psym=sym_array[papers_i],thick=3, color=color_array[1]
        ;if n_count3 GT 0 then cgoplot, k_arr[inds3,1],k_arr[inds3,0],psym=sym_array[papers_i],thick=3, color=color_array[2]
        ;if n_count4 GT 0 then cgoplot, k_arr[inds4,1],k_arr[inds4,0],psym=sym_array[papers_i],thick=3, color=color_array[3]
        ;if n_count5 GT 0 then cgoplot, k_arr[inds5,1],k_arr[inds5,0],psym=sym_array[papers_i],thick=3, color=color_array[4]
        ;if n_count6 GT 0 then cgoplot, k_arr[inds6,1],k_arr[inds6,0],psym=sym_array[papers_i],thick=3, color=color_array[5]
        
        
      endfor
      cglegend, Title=name_array[[2,4,5,6,7]], color=color_array[[2,4,5,6,7]], location=[.67,.32], psym=sym_num[[2,4,5,6,7]], charsize=1.25,Length=0

      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    endif
    
    stop
  endfor
  stop
  
end