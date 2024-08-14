pro five_antenna_model


restore, '1061316784_vis_model_XX.sav'

uvfits_read,hdr,params,layout,vis_arr,vis_weights,file_path_vis='/fred/oz048/MWA/data/2013/v5_1/1061316784.uvfits',n_pol=n_pol,silent=silent,error=error,_Extra=extra

defined_tiles = [1,11,21,51,101]

defined_inds = params.antenna1
defined_inds[*] = 0

inds0 = where(params.antenna1 EQ defined_tiles[0],n_count)
inds0auto = where(params.antenna2[inds0] EQ defined_tiles[0],n_count)
inds0a = where(params.antenna2[inds0] EQ defined_tiles[1],n_count)
inds0b = where(params.antenna2[inds0] EQ defined_tiles[2],n_count)
inds0c = where(params.antenna2[inds0] EQ defined_tiles[3],n_count)
inds0d = where(params.antenna2[inds0] EQ defined_tiles[4],n_count)

defined_inds[inds0[inds0auto]] = 1
defined_inds[inds0[inds0a]] = 1
defined_inds[inds0[inds0b]] = 1
defined_inds[inds0[inds0c]] = 1
defined_inds[inds0[inds0d]] = 1

inds1 = where(params.antenna1 EQ defined_tiles[1],n_count)
inds1auto = where(params.antenna2[inds1] EQ defined_tiles[0],n_count)
inds1a = where(params.antenna2[inds1] EQ defined_tiles[1],n_count)
inds1b = where(params.antenna2[inds1] EQ defined_tiles[2],n_count)
inds1c = where(params.antenna2[inds1] EQ defined_tiles[3],n_count)
inds1d = where(params.antenna2[inds1] EQ defined_tiles[4],n_count)

defined_inds[inds1[inds1auto]] = 1
defined_inds[inds1[inds1a]] = 1
defined_inds[inds1[inds1b]] = 1
defined_inds[inds1[inds1c]] = 1
defined_inds[inds1[inds1d]] = 1

inds2 = where(params.antenna1 EQ defined_tiles[2],n_count)
inds2auto = where(params.antenna2[inds2] EQ defined_tiles[0],n_count)
inds2a = where(params.antenna2[inds2] EQ defined_tiles[1],n_count)
inds2b = where(params.antenna2[inds2] EQ defined_tiles[2],n_count)
inds2c = where(params.antenna2[inds2] EQ defined_tiles[3],n_count)
inds2d = where(params.antenna2[inds2] EQ defined_tiles[4],n_count)

defined_inds[inds2[inds2auto]] = 1
defined_inds[inds2[inds2a]] = 1
defined_inds[inds2[inds2b]] = 1
defined_inds[inds2[inds2c]] = 1
defined_inds[inds2[inds2d]] = 1

inds3 = where(params.antenna1 EQ defined_tiles[3],n_count)
inds3auto = where(params.antenna2[inds3] EQ defined_tiles[0],n_count)
inds3a = where(params.antenna2[inds3] EQ defined_tiles[1],n_count)
inds3b = where(params.antenna2[inds3] EQ defined_tiles[2],n_count)
inds3c = where(params.antenna2[inds3] EQ defined_tiles[3],n_count)
inds3d = where(params.antenna2[inds3] EQ defined_tiles[4],n_count)

defined_inds[inds3[inds3auto]] = 1
defined_inds[inds3[inds3a]] = 1
defined_inds[inds3[inds3b]] = 1
defined_inds[inds3[inds3c]] = 1
defined_inds[inds3[inds3d]] = 1

inds4 = where(params.antenna1 EQ defined_tiles[4],n_count)
inds4auto = where(params.antenna2[inds4] EQ defined_tiles[0],n_count)
inds4a = where(params.antenna2[inds4] EQ defined_tiles[1],n_count)
inds4b = where(params.antenna2[inds4] EQ defined_tiles[2],n_count)
inds4c = where(params.antenna2[inds4] EQ defined_tiles[3],n_count)
inds4d = where(params.antenna2[inds4] EQ defined_tiles[4],n_count)

defined_inds[inds4[inds4auto]] = 1
defined_inds[inds4[inds4a]] = 1
defined_inds[inds4[inds4b]] = 1
defined_inds[inds4[inds4c]] = 1
defined_inds[inds4[inds4d]] = 1


inds = where(defined_inds EQ 1, n_count)
inds0 = where(params.antenna1[inds] EQ 128,n_count)
if n_count GT 0 then defined_inds[inds[inds0]] = 0
inds0 = where(params.antenna2[inds] EQ 128,n_count)
if n_count GT 0 then defined_inds[inds[inds0]] = 0
inds = where(defined_inds EQ 1, n_count)

end

