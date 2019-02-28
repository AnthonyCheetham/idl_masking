; JDM 2007Jan01
; extract_data_direct
;   This routine extracts the fourier amps and phases of a given
;      image by direct calculation not interpolation
;
;  extract_data_direct, image1, vis, phases (degs), scale=scale (mas/pixel), $
;       u=u, v=v
;
pro extract_data_direct, image1, vis0, phases0,scale=im_scale, u=u1, v=v1

if (keyword_set(cubic) eq 0) then cubic =0

if (keyword_set(help) eq 1) then begin
 print,'extract_data_direct, image, vis, phases (degs), scale=im_scale (mas/pixel), u=u, v=v
return
endif

u=u1(*)
v=v1(*)

dimy=(size(image1))(2)
dimx=(size(image1))(1)

im_scale_mas=mas2rad(im_scale)
x=-1d0*im_scale_mas*(findgen(dimx) -dimx/2) ; positive to left
y= 1d0*im_scale_mas*(findgen(dimy) - dimy/2)

make_2d,x,y,xx,yy
xx=xx(*)
yy=yy(*)
vals = image1(*) ;/total(image1*1.0)
nval=n_elements(vals)

u0=u(*)
v0=v(*)
nuv=n_elements(u0)

vals2d= rebin(reform(vals,1,nval),nuv,nval)
xtransform = complex(cos(-2.0*!pi*u0#xx),sin(-2.0*!pi*u0#xx))
ytransform = complex(cos(-2.0*!pi*v0#yy),sin(-2.0*!pi*v0#yy))

vcomplex = total(vals2d*xtransform*ytransform,2)
phases = atan(vcomplex,/phase)
phases = mod360(phases*180/!pi)
vis = abs(vcomplex)

vis0=u1 ; make sure vis0 and phases0 have some ndim format as input u1 structure
vis0(*)=vis
phases0=u1
phases0(*)=phases


end







    





