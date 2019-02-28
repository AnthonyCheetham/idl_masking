; JDM 2001Dec01		Returns Vis and Phase of Binary Disk for u,v
; JDM 2006Oct22		Modified for uniform ellipse
;           
; INPUTS:
; Model of a uniform ellipse
; Params:
;      a(0)=The Average of Major/Minor Axes in units of milliarcseconds
; If diameter is less than zero than invert result
;      a(1)=Axis Ratio (Major/Minor Axis)
;      a(2)=Angle of major Axis  (Degrees -- East of North
;      a(3)=The zero visibility intercept.  Generally should be 1.0
;
;
; u,v : units of Number of wavelengths (i.e., rad-1)
;
; outputs: 
;   visib  -- visibility
;   phases -- phases in degrees (-180,180)


pro uv_uniform_ellipse,u,v,params,visib,phases

Minor = mas2rad(2*params(0)/(1.0+ params(1)))
Major =(Minor*params(1))
diameter=mas2rad(abs(params(0)))
intercept=params(3)
Ratio=params(1)
Angle=params(2)
angle_nofw=90.0-angle

pi=3.1415926

; First we should rotate our coordinate system so that
; uprime is in direction of Major Axis
; vprime is in direction of Minor Axis
;
uprime=u*cos(angle_nofw*pi/180.)+ v*sin(angle_nofw*pi/180.)
vprime=v*cos(angle_nofw*pi/180.)- u*sin(angle_nofw*pi/180.)

; Now Scale these Appropriately to allow use of
; use J1
uprime=uprime*Major/diameter
vprime=vprime*minor/diameter
krad=sqrt(uprime*uprime+vprime*vprime)

index=where(krad eq 0,count)
if (count gt 0) then krad(index)=1e-15

      f=2*beselj(pi*diameter*krad,1) / (pi*diameter*krad)
negin=where(f lt 0,ctneg)

if (params(0) lt 0) then f=(intercept/abs(f)) else f=abs(f)*intercept
visib=f ; visib can be negative ONLY if you pass a NEGATIVE INTERCEPT.... 
	; best notation?
phases=f
phases(*)=0.
if ctneg ne 0 then phases(negin)=180.00
;stop

end
