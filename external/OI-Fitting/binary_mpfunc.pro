;

function binary_mpfunc, params, dp, vis2data0=vis2data0,t3data0=t3data0,$
                           vis2model0=vis2model0,t3model0=t3model0


;print,'Params: ',params
; note vis2data/t3data are from extract_vis2data calls, not the
;  original binary tables.
;vis2model0=vis2data0
vis2model0 = binary_vis2data(params,vis2data=vis2data0)


;t3model0=t3data0

t3model0=binary_t3data(params,t3data=t3data0)


residuals = [ (vis2data0.vis2data - vis2model0.vis2data)/(vis2data0.vis2err),$
              angle_diff(t3data0.t3phi,t3model0.t3phi)/t3data0.t3phierr]

return,residuals


end

