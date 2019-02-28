;2005Oct22	JDM

function uniform_ellipse_t3data,params,t3data=t3data,file=file

; For use with oidata library by acting on t3data from
; extract_t3data
; alternatively, can take a oidata file directly!

; There are four parameters to this fit.
;      a(0)=The Average of Major/Minor Axes in units of milliarcseconds (MAS)
; If diameter is less than zero than invert result
;      a(1)=Axis Ratio (Major/Minor Axis)
;      a(2)=Angle of major Axis  (Degrees -- East of North
;      a(3)=The zero visibility intercept.  Generally should be 1.0
;

if (keyword_set(t3data) eq 0 and keyword_set(file) eq 0) then begin
 print,'Must input some data or filename'
return,-1
endif

if (keyword_set(t3data) eq 0) then begin
  extract_t3data,file=file,t3data
endif

;else assume we got  t3data

model_t3data=t3data

uv_uniform_ellipse,t3data.u1,t3data.v1,params,visib1,phases1
uv_uniform_ellipse,t3data.u2,t3data.v2,params,visib2,phases2
uv_uniform_ellipse,t3data.u3,t3data.v3,params,visib3,phases3


model_t3data.t3amp=visib1*visib2*visib3
model_t3data.t3amperr=0.
model_t3data.t3phi=mod360(phases1+phases2+phases3)
model_t3data.t3phierr=0.

return,model_t3data
end

