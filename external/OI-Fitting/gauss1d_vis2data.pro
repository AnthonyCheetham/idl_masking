;2003Sep28 JDM
;2005Nov20 JDM
;2005Nov23 JDM

function gauss1d_vis2data,params,vis2data=vis2data,file=file,$
 vis2model=model_vis2data
  
; For use with oidata library by acting on vis2data from
; extract_vis2data
; alternatively, can take a oidata file directly!
; INPUTS:
; Model 
; Params:
;params=[point_flux,gauss_flux, gauss FWHM (mas)]


if (keyword_set(vis2data) eq 0 and keyword_set(file) eq 0) then begin
 print,'Must input some data or filename'
return,-1
endif

if (keyword_set(vis2data) eq 0) then begin
  extract_vis2data,file=file,vis2data
endif

;else assume we got vis2data

model_vis2data=vis2data
model_vis2data.vis2err=0.

pi=3.1415926
diameter=abs(mas2rad(params(2)))
intercept=params(1)
x=vis2data.baseline/vis2data.eff_wave ; sfu
index=where(x eq 0,count)
if (count gt 0) then x(index)=1e-15

;      f=exp(-3.57*diameter^2*x^2)
      f=exp(-3.55971*(diameter*x)^2)
if (params(2) gt 0) then f=abs(intercept*f) else f=abs(intercept/f)
f=f+params(0)
if count gt 0 then f(index)=params(0)+params(1)
model_vis2data.vis2data=f*f
;print,params,'**  ',f*f

return,model_vis2data
end




