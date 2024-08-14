pro beam_inds, psf=psf, inds_arr=inds_arr

n_freq = 384
nbaselines = 8128
n_pol=2
psf_dim=80
psf_resolution=50

inds_arr=Ptrarr(n_pol,n_freq,nbaselines)
inds_arr_single = PTRARR(psf_resolution+1, psf_resolution+1)

for pol_i=0,n_pol-1 do begin
for freq_i=0,n_freq-1 do begin

for i=0, psf_resolution do begin
for j=0, psf_resolution do begin

beam_single = (*(*(*psf.beam_ptr)[pol_i,freq_i,0])[i,j])

inds_arr_single[i,j] = PTR_new( where( abs(reform(*(*(*psf.beam_ptr)[pol_i,freq_i,0])[i,j],psf_dim,psf_dim)) EQ 0 ) )
(*(*(*psf.beam_ptr)[pol_i,freq_i,0])[i,j]) = beam_single[*inds_arr_single[i,j]]


endfor
endfor

(*psf.beam_ptr)[pol_i,freq_i,*] = (*psf.beam_ptr)[pol_i,freq_i,0]
inds_arr[pol_i,freq_i,*] = Ptr_new(inds_arr_single)

endfor
print, 'One more pol!'
endfor




return

end
