pro vis_add_noise

;Clear important variables
undefine, vis_ptr, vis_xx, vis_yy, random_arr, real_noise, imaginary_noise, noisy_vis_xx, noisy_vis_yy

;Restore the obs structure to get the noise standard deviation per visibility
restore, '/data4/MWA/FHD_Aug23/fhd_nb_bubbles_low_noise_retry/metadata/1061316296_obs.sav'
vis_noise=*obs.vis_noise

;Read in xx visibilities
restore, '/data4/MWA/FHD_Aug23/fhd_nb_bubbles_low_noise_retry/vis_data/1061316296_vis_model_XX.sav'
vis_xx = *vis_model_ptr

;Read in yy visibilities
restore, '/data4/MWA/FHD_Aug23/fhd_nb_bubbles_low_noise_retry/vis_data/1061316296_vis_model_YY.sav'
vis_yy = *vis_model_ptr
visibility_num=size(vis_xx)

;Get an array of 2,2,384,largenumber (real/imaginary, pol, freq, visibilities) of random numbers based off the machine time tag
random_arr= randomn(systime(/seconds),2,2,384,visibility_num[2]) ;real/imaginary, pol, freq, visibilities

;Set up an array of pol,freq,visibilities. The real and imaginary parts have the same visibility noise. Each visibility
;apparently has the same noise standard deviation...divide by the sq root of 30000, which is about 1000hrs of observation
vis_noise_large=DBLARR(2,384,visibility_num[2])
for i=0, visibility_num[2]-1 do vis_noise_large[*,*,i]=vis_noise[*,*]/sqrt(30000.)

;Get the real and imaginary noise, which are uncorrelated.
real_noise=vis_noise_large*random_arr[0,*,*,*]
imaginary_noise=vis_noise_large*random_arr[1,*,*,*]

;Add noise to each visibility based on real/imaginary and polarization
noisy_vis_xx=(real_part(vis_xx)+real_noise[0,*,*])+Complex(0,1)*(imaginary(vis_xx)+imaginary_noise[0,*,*])
noisy_vis_yy=(real_part(vis_yy)+real_noise[1,*,*])+Complex(0,1)*(imaginary(vis_yy)+imaginary_noise[1,*,*])
;noisy_vis_xx=(real_noise[0,*,*])+Complex(0,1)*(imaginary_noise[0,*,*])
;noisy_vis_yy=(real_noise[1,*,*])+Complex(0,1)*(imaginary_noise[1,*,*])
undefine, vis_model_ptr

;Save the visibilities in pointers (like the input format) and save in seperate xx and yy savefiles.
vis_model_ptr=ptr_new(/allocate_heap)
*vis_model_ptr=noisy_vis_xx
save, vis_model_ptr, filename='/data4/MWA/FHD_Aug23/fhd_nb_bubbles_low_noise_retry/vis_data/1061316296_vis_model_XX.sav'
*vis_model_ptr=noisy_vis_yy
save, vis_model_ptr, filename='/data4/MWA/FHD_Aug23/fhd_nb_bubbles_low_noise_retry/vis_data/1061316296_vis_model_YY.sav'


end