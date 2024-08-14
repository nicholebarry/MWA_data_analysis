pro export_1dkpower

  ;coarse_harm_width = 3
  ;wedge_angles=120d
  n_cubes=3
  
  cut = STRARR(n_cubes)
  cut[2]='wedge_cut_plus_res_cut'
  ;cut[2]='Aug23_longrunstyle'
  ;cut[2] = 'beardsley_thesis_list_8shist'
  ;cut[2] = 'beardsley_thesis_list'
  ;cut[1] = 'beardsley_thesis_list_firstthird'
  cut[1] = 'beardsley_thesis_list'
  ;cut[0] = 'beardsley_thesis_lis_noautocutt'
  cut[0] = 'beardsley_thesis_list_iono3_cut'
  
  chans='ch0-95_' & band='low'
  ;chans='ch0-127_'
  ;chans='ch96-191_' & band='high'
  ;chans='ch48-143_' & band='mid'
  ;chans='ch10-127_'
  ;chans=''
  
  win = strarr(n_cubes)
  win[2]=''
  win[1]=''
  win[0]=''
  
  spec_window='bh_'
  
  basefile = STRARR(n_cubes)
  basefile[2]=['/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps_jan2016/Combined_obs_'+cut[2]+'_cubeXX__even_odd_joint_'+win[2]+chans]
  ;basefile[2]=['/nfs/mwa-05/r1/EoRuvfits/EoR2013/fhd_nb_2013longrun_savedbp/ps_large/Combined_obs_'+cut[2]+'_cubeXX__even_odd_joint_tk_'+chans+'_']
  ;basefile[2]=['/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal/ps/Combined_obs_'+cut[2]+'_cubeXX__even_odd_joint_'+win[2]+chans]
  ;basefile[1]=['/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_Aug2017_globalbp_w_cable_w_digjump/ps/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans]
    ;basefile[1]=['/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal/ps/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+'ch0-127_']
    basefile[1]=['/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/ps/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+win[1]+chans]
        basefile[0]=['/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/ps/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans]
    
  ;basefile[0]=['/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/ps_nonmaster/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans]
  ;basefile[0]=['/nfs/mwa-05/r1/EoRuvfits/EoR2013/fhd_nb_2013longrun_savedbp/ps_large/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+win[0]+chans]
  
  cubes=['res']
  pols=['xx','yy']
  
  endfile = STRARR(n_cubes)
  endfile[2]='_bh_dencorr_no_120deg_wedge_cbw3_kperplambda10-60_kpar0.15-200_1dkpower.idlsave'
  ;endfile[2]='_averemove_bh_dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'
  ;endfile[2]='_averemove_'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'
  endfile[1]='_averemove_'+spec_window+'dencorr_no_120deg_wedge_kperplambda10-60_kpar0.15-200_1dkpower.idlsave'
  endfile[0]='_averemove_'+spec_window+'dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave'
  ;endfile[1]='_averemove_'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'
  ;endfile[0]='_averemove_'+spec_window+'dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'
  
  limit_percent=0.9772 ;2 sigma
  hubble_param=0.71
  
  weights_save = FLTARR(n_cubes,92)
  
  for j=0,1 do begin
    for cube_i=0,n_cubes-1 do begin
    
      file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
      restore,file
      ; find the best limit in this file
      n_k=n_elements(k_edges)
      k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
      delta=power*(k^3.)/(2.*!pi^2.)
      dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
      weights_save[cube_i,0:N_elements(weights)-1]=weights
      ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
      limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
      lim=min(limits,ind)
      header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut[cube_i]+' '+chans+' '+cubes+' '+win[cube_i]+' '+pols[j]+')'
      print,header
      
      ;power=[0.,power] ; first element a placeholder to match size of k_edges
      ;weights=[0.,weights]
      ;foo=transpose([[k_edges],[power],[weights]])
      ;outfile='~/'+band+'_'+cubes+'_'+pols[j]+'.txt'
      ;textfast,foo,[header,'',''],file_path=outfile,/write
      
      ;if (cube_i EQ 0) then cgPS_Open,'/nfs/mwa-00/h1/nbarry/'+pols[j]+'_supercut_limit_comparison'+chans+'.png',/quiet,/nomatch
      if cube_i EQ 0 then begin
        cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.08,1.5], yrange=[10^3.,10^8.],ytitle='$\Delta$$\up2$ (mK$\up2$)', $
          xtitle='k (h Mpc$\up-1$)', title=pols[j] + ' Limit Comparisons ' + chans, charsize =1.5, linestyle=1
        cgoplot, k,limits,/xlog,/ylog,psym=10
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2
      endif
      if cube_i EQ 1 then begin
        cgoplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.08,1.5],color='blue', linestyle=1
        cgoplot, k,(limits),/xlog,/ylog,psym=10, xrange=[.08,1.5],color='blue'
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2,color='blue'
      endif
      if cube_i EQ 2 then begin
        cgoplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.1,1.5],color='green', linestyle=1
        cgoplot, k,(limits),/xlog,/ylog,psym=10, xrange=[.08,1.5],color='green'
        cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2,color='green'
        cglegend, title=['Autocal + Tukey + Sampled UV','Std + Tukey',"Autocal + Reg"],color=['black','blue','green'], location=[.2,.85], charsize=1
      endif
    ;if (cube_i EQ 2) then cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    ;textfast,foo,file_path=outfile,/write
    endfor
    stop
  endfor
  stop
  
end