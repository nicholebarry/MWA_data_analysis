pro kz_fft

restore, '/nfs/eor-00/h1/nbarry/1061316296_even_gridded_uvf.sav'

;Get cube slice line
modelsize=size(*model_arr[0,0])
zcube=fltarr(modelsize[1],modelsize[2],192)
for i=0,191 do zcube[*,*,i]=(*model_arr[0,i])
;zcube_fft=fft(zcube1,DIMENSION=3)

for i=0,191 do i_use[i]=where((freq_bin_i2 EQ i) AND (freq_use GT 0),nf_use)
redshifts_cube=(1420.4057517*10^6.-(*obs.baseline_info).freq[i_use])/(*obs.baseline_info).freq[i_use]
cosmology_measures,redshifts_cube,comoving_dist_los=comoving_dist_los

for i=0,191 do yaxis[i]=zcube[0,0,i]
zcube_fft=fft(zcube, DIMENSION=3)


dchi = abs(comoving_dist_los[191]-comoving_dist_los[1])/(192-1)
kz_max = 2*pi / (2*dchi); % factor of 2 for positive and negative

dk=(2*!Pi)/(max(comoving_dist_los[*])-min(comoving_dist_los[*]))

;kz_axis = -kz_max:dkz:kz_max; % kz (k_parallel) values for each voxel

end