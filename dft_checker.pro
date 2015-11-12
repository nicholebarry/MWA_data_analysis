pro dft_checker

  gain_arr=complex(FLTARR(384))
  i_comp=Complex(0,1)
  n_freq=380
  
  gain_arr=10*Sin(2*!Pi*5*(((FINDGEN(100))/100))+!Pi);+i_comp;+5 ;frequency of 2, once a second for 380/380 seconds, amplitude of 10
  ;gain_arr=(10*exp(-i_comp*2*!Pi*5*((FINDGEN(100))/100)))
  ;gain_arr[0:1]=0
  ;gain_arr[382:383]=0
  
  mode0=5 ; start with nominal cable length
  dmode=.05 ; pretty fine
  ;dmode=0.1
  nmodes=90 ; range around the central mode to test
  modes=(dindgen(nmodes)-nmodes/2)*dmode+mode0 ; array of modes to try
  
  modes=rebin(modes,nmodes,100)
  
  ;gainr=rebin(transpose(reform(real_part(gain_arr[2:381]))),nmodes,n_freq)
  ;gaini=rebin(transpose(reform(imaginary(gain_arr[2:381]))),nmodes,n_freq) ; and this...
    gainr=rebin(transpose(reform(real_part(gain_arr))),nmodes,100)
  gaini=rebin(transpose(reform(imaginary(gain_arr))),nmodes,100) ; and this...
  ;*************
  gain_temp=gainr+i_comp*gaini ; for some reason I cant rebin complex numbers
  ;freq_mat=rebin(transpose(freq_use),nmodes,nf_use) ; this too...
  ;*************MOD
  ;freq_mat=rebin(transpose(Findgen(n_freq)+2),nmodes,n_freq) ; this too...
  freq_mat=rebin(transpose(Findgen(100)),nmodes,100) ; this too...
  ;**************
  test_fits=2*Total(exp(i_comp*2.*!Pi/100*modes*(freq_mat))*gain_temp,2)
  
  gainr=rebin(transpose(reform(real_part(test_fits))),nmodes,100)
  gaini=rebin(transpose(reform(imaginary(test_fits))),nmodes,100)
  test_fits_temp=gainr+i_comp*gaini
  
  amp_use=max(abs(test_fits),mode_ind)/100
  mode_i=modes[mode_ind,0]
  test_fits1=Total(exp(-i_comp*2.*!Pi/100*abs(mode_i)*(freq_mat))*test_fits_temp/100,1)/1000

  ;test_fits2=Total(exp(i_comp*2.*!Pi*(modes/bandwidth)*(freq_mat/n_freq*bandwidth+min(freq_arr)))*gain_temp,2)
  
  amp_use=max(abs(test_fits),mode_ind)/100
  phase_use=atan(test_fits[mode_ind],/phase)
  mode_i=modes[mode_ind,0]
  ;gain_mode_fit=amp_use*exp(-i_comp*2.*!Pi*(mode_i*(findgen(384))/384)+i_comp*phase_use)
  
  stop
  
  
end