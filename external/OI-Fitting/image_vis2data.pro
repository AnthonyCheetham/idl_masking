;2003Sep28 JDM
;2005Nov20 JDM
;2005Nov23 JDM

function image_vis2data,im,vis2data=vis2data,file=file, scale=scale, $
  direct = direct 

;if /direct then use DIRECT calculation instead of interp on FFT.

if keyword_set(scale) eq 0 then scale = 1.0 ;mas per pix

; For use with oidata library by acting on vis2data from
; extract_vis2data
; alternatively, can take a oidata file directly!
; INPUTS:
; Model 
; Params:


if (keyword_set(vis2data) eq 0 and keyword_set(file) eq 0) then begin
 print,'Must input some data or filename'
return,-1
endif

if (keyword_set(vis2data) eq 0) then begin
  extract_vis2data,file=file,vis2data
endif

;else assume we got vis2data

model_vis2data=vis2data
if keyword_set(direct) eq 0 then begin
  extract_data,im,visib,phases,scale=scale,u=vis2data.u,v=vis2data.v,/cubic,$
    /nohan
endif else begin
  extract_data_direct,im,visib,phases,scale=scale,u=vis2data.u,v=vis2data.v
endelse


visib_real =      visib*cos(phases*!pi/180.)
visib_imag =      visib*sin(phases*!pi/180.)

model_vis2data.vis2data=visib_real^2 + visib_imag^2
model_vis2data.vis2err=0.
;stop
return,model_vis2data
end




