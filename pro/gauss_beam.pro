


pix_hor=obs.dimension/antenna1.psf_scale
cen=pix_hor/2

p = [$
  1.,cen,pix_hor/14,cen,pix_hor/14,$ ;pr
  0.07,cen-pix_hor/3.5,pix_hor/14,cen,pix_hor/7,$ ;sw1
  0.01,cen-pix_hor*13/14,pix_hor/14,cen,pix_hor/7,$ ;sw2
  0.07,cen+pix_hor/3.5,pix_hor/14,cen,pix_hor/7,$ ;se1
  0.01,cen+pix_hor*13/14,pix_hor/14,cen,pix_hor/7,$;se2
  0.07,cen,pix_hor/7,cen+pix_hor/3.5,pix_hor/14,$ ;sn1
  0.01,cen,pix_hor/7,cen+pix_hor*13/14,pix_hor/14,$ ;sn2
  0.07,cen,pix_hor/7,cen-pix_hor/3.5,pix_hor/14,$ ;ss1
  0.01,cen,pix_hor/7,cen-pix_hor*13/14,pix_hor/14,$ ;ss2
  0.003,cen+pix_hor/3.5,pix_hor/14,cen+pix_hor/3.5,pix_hor/14,$ ;sne1
  0.003,cen-pix_hor/3.5,pix_hor/14,cen+pix_hor/3.5,pix_hor/14,$ ;snw1
  0.003,cen+pix_hor/3.5,pix_hor/14,cen-pix_hor/3.5,pix_hor/14,$ ;sse1
  0.003,cen-pix_hor/3.5,pix_hor/14,cen-pix_hor/3.5,pix_hor/14];ssw1

x = FINDGEN(pix_hor)
y = FINDGEN(pix_hor)
range = [psf_image_dim/2-pix_hor/2+1,psf_image_dim/2+pix_hor/2-1]


nrms = MPFIT2DFUN('gaussian_decomp_beam', x, y, abs(image_power_beam[range[0]:range[1],range[0]:range[1]]), 1 , p, weights=1d)

pix_hor=obs.dimension/antenna1.psf_scale
psf_image_dim = antenna1.psf_image_dim
cen=pix_hor
 
 p = [$
   1.,cen,pix_hor/14.4,cen,pix_hor/14.4,$ ;pr
   -0.04,cen,pix_hor/28,cen,pix_hor/28,$ ;pr
   -0.04,cen,pix_hor/28,cen-pix_hor/7.,pix_hor/28,$ ;pr
   -0.04,cen,pix_hor/28,cen+pix_hor/7.,pix_hor/28,$ ;pr
   -0.04,cen-pix_hor/7.,pix_hor/28,cen,pix_hor/28,$ ;pr
   -0.04,cen+pix_hor/7.,pix_hor/28,cen.,pix_hor/28,$ ;pr
   0.05,cen-pix_hor/3.65,pix_hor/28,cen,pix_hor/14.4,$ ;sw1
   0.007,cen-pix_hor/2.3,pix_hor/28,cen,pix_hor/14.4,$ ;sw2
   0.05,cen+pix_hor/3.65,pix_hor/28,cen,pix_hor/14.4,$ ;se1
   0.007,cen+pix_hor/2.3,pix_hor/28,cen,pix_hor/14.4,$;se2
   0.063,cen,pix_hor/14.4,cen+pix_hor/3.65,pix_hor/28,$ ;sn1
   0.02,cen,pix_hor/14.4,cen+pix_hor/2.3,pix_hor/28,$ ;sn2
   0.063,cen,pix_hor/14.4,cen-pix_hor/3.65,pix_hor/28,$ ;ss1
   0.02,cen,pix_hor/14.4,cen-pix_hor/2.3,pix_hor/28,$ ;ss2
   0.003,cen+pix_hor/3.65,pix_hor/28,cen+pix_hor/3.65,pix_hor/28,$ ;sne1
   0.003,cen-pix_hor/3.65,pix_hor/28,cen+pix_hor/3.65,pix_hor/28,$ ;snw1
   0.003,cen+pix_hor/3.65,pix_hor/28,cen-pix_hor/3.65,pix_hor/28,$ ;sse1
   0.003,cen-pix_hor/3.65,pix_hor/28,cen-pix_hor/3.65,pix_hor/28];ssw1

;   -0.04,cen,pix_hor/28,cen,pix_hor/28,$ ;pr
;   -0.04,cen,pix_hor/28,cen-pix_hor/7.,pix_hor/28,$ ;pr
;   -0.04,cen,pix_hor/28,cen+pix_hor/7.,pix_hor/28,$ ;pr
;   -0.04,cen-pix_hor/7.,pix_hor/28,cen,pix_hor/28,$ ;pr
;   -0.04,cen+pix_hor/7.,pix_hor/28,cen.,pix_hor/28,$ ;pr
 
 x = FINDGEN(cen*2)
 y = FINDGEN(cen*2)
 range = [psf_image_dim/2-cen,psf_image_dim/2+cen-1]
 nrms = MPFIT2DFUN('gaussian_decomp_beam', x, y, abs(image_power_beam[range[0]:range[1],range[0]:range[1]]), 1 , p, weights=1d) 

result=gaussian_decomp_beam(FINDGEN(psf_image_dim), FINDGEN(psf_image_dim), nrms, model_npix=cen*2)
