; JDM March 30, 1997
; JDM Feb 24, 1998    Based  on clean.pro
; JDM 		      Zero padding.
; JDM 2005Nov23		FIXED u->-u problem !
;		Validated against binary_disks.pro
; AS LONG IMAGE is East to the LEFT, North UP
; JDM 2007Jan01. increased default zero padding for improved accuracy.

;
; extract_data
;   This routine extracts the fourier amps and phases of a given
;      image.
;
;  extract_data, image1, vis, phases (degs), scale=scale (mas/pixel), u=u, v=v
;
pro extract_data, image1, vis, phases,scale=im_scale, u=u, v=v,cubic=cubic,$
 help=help,nohan=nohan,nopad=nopad,fft_scale=scale,amp=a,theta=t,get=get,$
 noprint=noprint

if (keyword_set(cubic) eq 0) then cubic =0

if (keyword_set(help) eq 1) then begin
 print,'extract_data, image, vis, phases (degs), scale=im_scale (mas/pixel), u=u, v=v,/cubic,/nohan,/nopad
return
endif

n_holes=(size(u))(1)
phases=fltarr(n_holes,n_holes)
vis=fltarr(n_holes,n_holes)
ydim=(size(image1))(2)
xdim=(size(image1))(1)

han=hanning(xdim,ydim)
if (keyword_set(nohan) eq 1) then  han=1.0

;scale=xdim/max(b_lengths)*.8/2.

im_scale_mas=mas2rad(im_scale)
if (keyword_set(nopad) eq 1) then factor=1 else factor=4


scale=1.0/(xdim*factor*im_scale_mas)

if (keyword_set(get) eq 1) then begin
if (keyword_set(noprint) eq 0) then begin
 print,'Fourier Scale:  1 pixel in Fourier Space = ',scale/1e5,' SFUs'
 print,'                Highest Spatial Frequency= ',scale/1e5*xdim/2.,' SFUs'
endif
endif
sfus=sqrt(u^2+v^2)
maxsfu=max(sfus/1e5)
minsfu=min(  (sfus/1e5)(where(sfus gt 0)))
if (keyword_set(get) eq 1) then begin
 if (keyword_set(noprint) eq 0) then $
  print,'DATA : Lowest SFU = ',minsfu,'     Highest SFU = ',maxsfu
endif
image=fltarr(xdim*factor,ydim*factor)
image(0:xdim-1,0:ydim-1)=image1*han
; Attempt to do bilinear and tri-linear tinerpolation.
; Line up based
image_fft=fft(    shift(image,-xdim/2 +1,-ydim/2 +1 ),/inv)

image_fft=shift(image_fft,xdim*factor/2.,ydim*factor/2)
realimage=float(image_fft)
imagimage=imaginary(image_fft)
if (keyword_set(get) eq 1) then begin
  ri2at,realimage,imagimage,a,t
  theta=mod360(t)
  a=shift(a,-xdim*factor/2,-ydim*factor/2)
  t=shift(t,-xdim*factor/2,-ydim*factor/2)
endif
; JDM Added MINUS SIGNS to U to deal with that u points left. 
realdata=interpolate(realimage,  (-u/scale)+xdim*factor/2,(v/scale)+ydim*factor/2,cubic=cubic)
imagdata=interpolate(imagimage,  (-u/scale)+xdim*factor/2,(v/scale)+ydim*factor/2,cubic=cubic)
; JDM Added MINUS SIGN to imagdata to be consistent with closure phase convention
ri2at,realdata,-imagdata,vis,phases
;stop
;q=full_uv(mean_x=u/scale+xdim/2,mean_y=v/scale+ydim/2,val=realdata)
;q1=full_uv(mean_x=u/scale+xdim/2,mean_y=v/scale+ydim/2,val=realdata1)
phases=mod360(phases)
end







    





