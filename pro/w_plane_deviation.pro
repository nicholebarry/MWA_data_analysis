
; M. Katz 1/26/04
; IDL function to perform a least-squares fit a plane, based on
;   Ax + By + C = z
;
;   ABC = plane_fit(x, y, z, error=error)
;
function plane_fit, x, y, z, error=error, noerror=noerror,noshow=noshow

  tx2 = total(x^2)
  ty2 = total(y^2)
  txy = total(x*y)
  tx = total(x)
  ty = total(y)
  N = n_elements(x)
  
  A = [[tx2, txy, tx], $
    [txy, ty2, ty], $
    [tx,  ty,  N ]]
    
  b = [total(z*x), total(z*y), total(z)]
  
  out = invert(a) # b
  
  if not keyword_set(noshow) then begin
    print, 'Plane Fit: Ax + By + C = z'
    print, 'A = ', out(0)
    print, 'B = ', out(1)
    print, 'C = ', out(2)
  endif
  
  if not keyword_set(noerror) then begin
    error = stdev(out(0)*x + out(1)*y + out(2) - z)
    if not keyword_set(noshow) $
      then print, 's = ', error
  endif
  
  return, out ;--- [A,B,C]
end

function SCATTER3D, x,y,z, PostScript=postscript,xtitle=xtitle,ytitle=ytitle,ztitle=ztitle,title=title

  ; Set the PostScript keyword to send draw the plot in a PostScript file
  ; instead of on the display.
  IF Keyword_Set(postscript) THEN cgPS_Open, FONT=1, Charsize=3.0,'/nfs/eor-00/h1/nbarry/scatter3D.png',/quiet,/nomatch
  
  ; Create the random data. Set the seed so you see what I see.
  ;seed = 1L
  ;x = RANDOMU(seed, 32)
  ;y = RANDOMU(seed, 32)
  ;z = EXP(-3 * ((x - 0.5)^2 + (y - 0.5)^2))
  
  ; Load a color table and create colors for the scatterplot.
  cgLoadCT, 33
  zcolors = BYTSCL(z)
  
  ; Set the 3D coordinate space with axes.
  cgSurf, DIST(5), /NODATA, /SAVE, XRANGE=[min(x),max(x)], $
    YRANGE=[min(y),max(y)], ZRANGE=[min(z), max(z)], XSTYLE=1, $
    YSTYLE=1, ZSTYLE=1, CHARSIZE=2.0, $
    POSITION=[0.1, 0.1, 0.95, 0.95, 0.1, 0.95], $
    XTICKLEN=1, YTICKLEN=1, XGRIDSTYLE=1, YGRIDSTYLE=1,xtitle=xtitle, ytitle=ytitle,ztitle=ztitle,title=title
  cgAXIS, XAXIS=1, /T3D, CHARSIZE=2.0
  cgAXIS, YAXIS=1, /T3D, CHARSIZE=2.0
  
  ; Plot the random points in 3D space with a filled circle shape.
  ;phi = Findgen(32) * (!PI * 2 / 32.)
  ;phi = [ phi, phi(0) ]
  cgPlotS, x, y, z, PSYM=16, COLOR=zcolors, SYMSIZE=2.5, /T3D
  
  ; Connect the data points to the XY plane of the plot.
  FOR j=0,N_elements(z)-1 DO cgPlotS, [x[j], x[j]], [y[j], y[j]], [min(z), z[j]], $
    COLOR=zcolors[j], /T3D
    
  ; Close the PostScript file and clean-up, if required.
  IF Keyword_Set(postscript) THEN cgPS_Close,Density=300,Resize=100.,/allow_transparent,/nomessage;, /PNG
  
END

function plot_plane,data_array,coeffs, deviations_z,PostScript=postscript,xtitle=xtitle,ytitle=ytitle,ztitle=ztitle,title=title

  IF Keyword_Set(postscript) THEN cgPS_Open, FONT=1, Charsize=3.0,'/nfs/eor-00/h1/nbarry/plot_plane.png',/quiet,/nomatch
  
  x_elements=floor(max(data_array[1,*])/10.)-floor(min(data_array[1,*])/10.)
  y_elements=floor(max(data_array[2,*])/10.)-floor(min(data_array[2,*])/10.)
  plane_surface=FLTARR(x_elements,y_elements)
  x_surface = FINDGEN(x_elements)*10.+floor(min(data_array[1,*]))
  y_surface = FINDGEN(y_elements)*10.+floor(min(data_array[2,*]))
  for x_i=0, x_elements-1 do plane_surface[x_i,*] = coeffs[0]*x_surface[x_i]+coeffs[1]*y_surface[*]+coeffs[2]
  
  ;Set up axes and whatnot, no data
  cgSurf, plane_surface,x_surface,y_surface, /Nodata, /SAVE, XRANGE=[min(data_array[1,*]),max(data_array[1,*])], $
    YRANGE=[min(data_array[2,*]),max(data_array[2,*])], ZRANGE=[min(data_array[3,*]), max(data_array[3,*])], XSTYLE=1, $
    YSTYLE=1, ZSTYLE=1, CHARSIZE=2.0, $
    POSITION=[0.1, 0.1, 0.95, 0.95, 0.1, 0.95], $
    XTICKLEN=1, YTICKLEN=1, XGRIDSTYLE=1, YGRIDSTYLE=1,xtitle=xtitle, ytitle=ytitle,ztitle=ztitle,title=title
    
  ; Load a color table and create colors for the tiles
  cgLoadCT, 33
  zcolors = BYTSCL(deviations_z)
  
  ;circles where the tiles are, colored by deviation
  cgPlotS, reform(data_array[1,*]),reform(data_array[2,*]), reform(data_array[3,*]), PSYM=16, COLOR=zcolors, SYMSIZE=2.5, /T3D
  
  ;actually plot the surface
  cgSurf, plane_surface,x_surface,y_surface, /Noerase, /SAVE, XRANGE=[min(data_array[1,*]),max(data_array[1,*])], $
    YRANGE=[min(data_array[2,*]),max(data_array[2,*])], ZRANGE=[min(data_array[3,*]), max(data_array[3,*])], XSTYLE=1, $
    YSTYLE=1, ZSTYLE=1, CHARSIZE=2.0, $
    POSITION=[0.1, 0.1, 0.95, 0.95, 0.1, 0.95], $
    XTICKLEN=1, YTICKLEN=1, XGRIDSTYLE=1, YGRIDSTYLE=1,xtitle=xtitle, ytitle=ytitle,ztitle=ztitle,title=title
    
  plot_above_inds = where(deviations_z LT 0.)
  cgPlotS, reform(data_array[1,plot_above_inds]),reform(data_array[2,plot_above_inds]), reform(data_array[3,plot_above_inds]), PSYM=16, COLOR=zcolors[plot_above_inds], SYMSIZE=2.5, /T3D
  
  IF Keyword_Set(postscript) THEN cgPS_Close,Density=300,Resize=100.,/allow_transparent,/nomessage;, /PNG
  
end

pro w_plane_deviation

  filename = '/nfs/eor-00/h1/nbarry/MWA/projects/flagged_project/antenna_locations_textfast.txt'
  textfast, data_array,/read,file_path=filename
  
  solar_tile_list=[78,79,80,87,88,95,96,97,98,104,105,112,113,122,123,124]
  solar_tile_inds = INTARR(128)+1
  solar_tile_inds[solar_tile_list]=0
  ;data_array = data_array[*,where(solar_tile_inds EQ 1)]
  
  coeffs = plane_fit(reform(data_array[1,*]),reform(data_array[2,*]),reform(data_array[3,*]),error=error)
  
  plane_z = coeffs[0]*reform(data_array[1,*])+coeffs[1]*reform(data_array[2,*])+coeffs[2]
  deviations_z = plane_z - reform(data_array[3,*])
  
  plane = plot_plane(data_array, coeffs, deviations_z, xtitle='East-West (meters)',ytitle='North-South (meters)',$
    ztitle='Height (meters)',title='Coplaner fit for MWA128',postscript=1)
  
  scatter = scatter3D(reform(data_array[1,*]),reform(data_array[2,*]),deviations_z,xtitle='East-West (meters)',$
    ytitle='North-South (meters)', ztitle='Deviations from coplaner (meters)',title='Devations from coplaner for MWA128',postscript=1)
    
  ;Fov =30deg, max and min angle deviation for min freq and max freq in high band  
  max_ang_freq = abs(deviations_z*(cos(!Pi/12.)-1.)/1.5218*2.*!Pi*2.)  
  min_ang_freq = abs(deviations_z*(cos(!Pi/12.)-1.)/1.7845*2.*!Pi*2.)  
  
  print, 'Stddev max: ' + strtrim(stddev(max_ang_freq),2) + ', r: ' + strtrim(1.-.5*stddev(max_ang_freq)^2.,2)
  print, 'Max: ' + strtrim(max(max_ang_freq),2) + ', r: ' + strtrim(cos(max(max_ang_freq)),2)
  print, 'Stddev min: ' + strtrim(stddev(min_ang_freq), 2) + ', r: ' + strtrim(1.-.5*stddev(min_ang_freq)^2.,2)
  print, 'Max: ' + strtrim(max(min_ang_freq),2) + ', r: ' + strtrim(cos(max(min_ang_freq)),2)
  
  stop
  
end

