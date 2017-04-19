FUNCTION conv,f1,f2,method=method,nosort=nosort,cut=cut
;+
; NAME:
; CONV
;
; PURPOSE:
; This function performs the mathematical convolution of two 
; input functions f1 and f2, each one defined as fltarr(2,npts).
; Both functions must have the same abscissas axes and the same
; number of points. (Use CONV2 if the abscissas are not the same)
;
; CATEGORY:
; Mathematics
;
; CALLING SEQUENCE:
;
; result = CONV(f1,f2)
;
; INPUTS:
; f1 and f2: the input functions
;
; OPTIONAL INPUTS:
; 
; KEYWORD PARAMETERS:
; METHOD: 0 (default) calculates convolution using the FFT  
;   1 calculates convolution performing the integrals
;
; OUTPUTS:
;       the convoluted function defined as
;       h(t) = int[ f2(y) f1(t-y) dy ]
;       Be careful with the order, f2 is "fixed" and f1 is "moving"!
;
;
; SIDE EFFECTS:
; the array of the convoluted function is fltarr(2,2*npts), being
; npts the number of points for the input functions
;
; PROCEDURE:
;
; EXAMPLE:
;
;   a = conv(f1,f2,METHOD=1)
;
; MODIFICATION HISTORY:
;   Written by: M. Sanchez del Rio and C. Ferrero
; April, 1993 
; 95-03-15 MSR adds nosort and cut keywords and correct a bug that
;   shifted the result data when the zero is not centered in the 
;   incoming abscissas array.
;       96/01/05 PF+MSR correct abscissas shift.
;-
;
;
; makes the convolution of two sets of data (each set is a fltarr(2,npts) )
;
;
;
on_error,2
;
; test correctness of array dimension and abscissas arrays
;
if (n_elements(f1(0,*)) ne n_elements(f2(0,*))) then begin
  print,' Error from CONV: mismatch of array dimensions. Abort'
  return,0
endif
diff = abs(f1(0,1) - f2(0,1))
if (diff gt 1e-4) then begin
  print,' Error from CONV: mismatch of abscissas arrays. Abort'
  return,0
endif

if not(keyword_set(nosort)) then begin
  xsort = sort(f1(0,*))
  f1(0,*)=f1(0,xsort)
  f1(1,*)=f1(1,xsort)
  xsort = sort(f2(0,*))
  f2(0,*)=f2(0,xsort)
  f2(1,*)=f2(1,xsort)
endif
;
; define abscissas arrays
;
npts = n_elements(f1(0,*))
index=-1
for i=1L,npts-1 do if f1(0,i) EQ 0 then index=i

if index GT 0 then npts2 = (index+(index+1))*2 else begin
  npts2 = 2*npts
  index=npts/2
endelse

;print,'npts: ',npts
;print,'npts2: ',npts2
;print,'index: ',index
step = f1(0,1) - f1(0,0)
xx=fltarr(npts2)
xx(0:npts-1)=f1(0,*)
j=0L
for i=npts,npts2-1 do begin
  xx(i) = xx(npts-1)+step*(j+1)
  j=j+1
endfor
xx = xx + f1(0,0)
;
; calculates convolution of z2 and z3 functions
;
z2 = fltarr(npts2)
z2(0:npts-1) = f1(1,*)
z3 = fltarr(npts2)
z3(0:npts-1) = f2(1,*)
h = fltarr(npts2)

;plot,xx,z2
;plotfile,f1,ki=2,psym=1
;pause
;
; calculate the convolution. Select one method:
;    METHOD=0 Using the convolution theorem
;    METHOD=1 Making the integral
;
if not(keyword_set(method)) then begin
  ;
  ; makes the convolution using the convolution theorem
  ;
  ff1 = fft(z2,-1)
  ff2 = fft(z3,-1)
  ff3 = ff1*ff2         ;*conj(ff2)
  h = fft(ff3,+1)
  h = float(h) *npts2
endif else begin
  ;
  ;
  ; makes the integral
  ;
  for j=0,npts2-1 do begin
    h(j) = 0
    for i=0,j do begin
        n= j - i
        h(j) = h(j) + z2(n)*z3(i)
    endfor
  endfor
endelse
;
; normalization factor = 1/sqrt( integral(f1)*integral(f2) )
;
int1 = total(f1(1,*))  ;*step 
int2 = total(f2(1,*))  ;*step
aa= total(h)*step/(int1*int2)
h = h*aa
;
; prepare the output
;
res = fltarr(2,npts2)
res(0,*) = xx  ;- 0.5*( f1(0,npts-1) - f1(0,0) ) +  58.0796

res(1,*) = h
if keyword_set(cut) then res=res(*,index:index+npts-1)
return,res
;
end