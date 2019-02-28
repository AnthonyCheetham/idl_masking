;

function uniform_ellipse_mpfunc, params, dp, vis2data0=vis2data0,t3data0=t3data0,$
                           vis2model0=vis2model0,t3model0=t3model0


;print,'Params: ',params
; note vis2data/t3data are from extract_vis2data calls, not the
;  original binary tables.
;vis2model0=vis2data0
vis2yes=0
t3yes=0

if keyword_set(vis2data0) ne 0 then begin
 vis2model0 = uniform_ellipse_vis2data(params,vis2data=vis2data0)
 vis2yes=1
endif

;t3model0=t3data0
if (keyword_set(t3data0) ne 0) then  begin
 t3model0=uniform_ellipse_t3data(params,t3data=t3data0)
 t3yes=1
endif

if (t3yes eq 1 and vis2yes eq 1) then begin
residuals = [ (vis2data0.vis2data - vis2model0.vis2data)/(vis2data0.vis2err),$
              angle_diff(t3data0.t3phi,t3model0.t3phi)/t3data0.t3phierr]
endif

if (t3yes eq 0 and vis2yes eq 1) then begin
residuals = [ (vis2data0.vis2data - vis2model0.vis2data)/(vis2data0.vis2err)]
              
endif

if (t3yes eq 1 and vis2yes eq 0) then begin
residuals = [ angle_diff(t3data0.t3phi,t3model0.t3phi)/t3data0.t3phierr]
endif



return,residuals


end

