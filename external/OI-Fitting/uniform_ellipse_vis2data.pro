;	    JDM 2006Oct22
;          
; %%%%% This is a function to return the visibility of an elliptical uniform
;       disk as a function of baseline.
; 	This will utilize the trick that its just a scaling from J1.
;
function uniform_ellipse_vis2data, params,vis2data=vis2data,file=file

; There are four parameters to this fit.
;      a(0)=The Average of Major/Minor Axes in units of milliarcseconds (MAS)
; If diameter is less than zero than invert result
;      a(1)=Axis Ratio (Major/Minor Axis)
;      a(2)=Angle of major Axis  (Degrees -- East of North
;      a(3)=The zero visibility intercept.  Generally should be 1.0
;

if (keyword_set(vis2data) eq 0 and keyword_set(file) eq 0) then begin
 print,'Must input some data or filename'
return,-1
endif

if (keyword_set(vis2data) eq 0) then begin
  extract_vis2data,file=file,vis2data
endif

;else assume we got vis2data

model_vis2data=vis2data

uv_uniform_ellipse,vis2data.u,vis2data.v,params,visib,phases

model_vis2data.vis2data=visib^2
model_vis2data.vis2err=0.

return,model_vis2data

end
