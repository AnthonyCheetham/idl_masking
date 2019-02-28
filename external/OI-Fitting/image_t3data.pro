;2003Sep28 JDM
;2005Nov20 JDM
;2005Nov23 jdm


function image_t3data,im,t3data=t3data,file=file,scale=scale, direct=direct

; For use with oidata library by acting on t3data from
; extract_t3data
; alternatively, can take a oidata file directly!

; INPUTS:
; Model
; Params:
;params=[point_flux,gauss FWHM major, gauss FWHM MINOR, Major Angle,skew,theta0-deg]

if (keyword_set(t3data) eq 0 and keyword_set(file) eq 0) then begin
 print,'Must input some data or filename'
return,-1
endif

if (keyword_set(t3data) eq 0) then begin
  extract_t3data,file=file,t3data
endif

;else assume we got  t3data

model_t3data=t3data
if keyword_set(direct) eq 0 then begin
  extract_data,im,visib1,phases1,scale=scale,u=t3data.u1,v=t3data.v1,/cubic,$
    /nohan
endif else begin
  extract_data_direct,im,visib1,phases1,scale=scale,u=t3data.u1,v=t3data.v1
endelse


visib1_real =      visib1*cos(phases1*!pi/180.)
visib1_imag =      visib1*sin(phases1*!pi/180.)

if keyword_set(direct) eq 0 then begin
  extract_data,im,visib2,phases2,scale=scale,u=t3data.u2,v=t3data.v2,/cubic,$
    /nohan
endif else begin
  extract_data_direct,im,visib2,phases2,scale=scale,u=t3data.u2,v=t3data.v2
endelse

visib2_real =      visib2*cos(phases2*!pi/180.)
visib2_imag =      visib2*sin(phases2*!pi/180.)

if keyword_set(direct) eq 0 then begin
  extract_data,im,visib3,phases3,scale=scale,u=t3data.u3,v=t3data.v3,/cubic,$
    /nohan
endif else begin
 extract_data_direct,im,visib3,phases3,scale=scale,u=t3data.u3,v=t3data.v3
endelse


visib3_real =      visib3*cos(phases3*!pi/180.)
visib3_imag =      visib3*sin(phases3*!pi/180.)

ri2at,visib1_real,visib1_imag,visib1,phases1
ri2at,visib2_real,visib2_imag,visib2,phases2
ri2at,visib3_real,visib3_imag,visib3,phases3

model_t3data.t3amp=visib1*visib2*visib3
model_t3data.t3amperr=0.
model_t3data.t3phi=mod360(phases1+phases2+phases3)
model_t3data.t3phierr=0.
return,model_t3data
end




