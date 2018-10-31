;; this program does the actual calibration of differential V2 and vis
;; quantities for gpi pol data. It is called by gpi_pol_double_diff.pro.
;;
;; Setting uncertainties=1 will return the standard deviation of the samples.
;;
;; Setting average=1 will average the v2 after the first difference measurement
;;
;; Author: ACC 2014.

function gpi_pol_cal_v2,v2_h_0,v2_v_0,v2_h_45,v2_v_45,olog,uncertainties=uncertainties,average_together=average_together
;;This does amplitude quantities (like V2).

;;First cal step:
v2_hv=v2_h_0/v2_v_0
v2_vh=v2_h_45/v2_v_45

;;Average
if keyword_set(average_together) then begin
    V2_hv=mean(v2_hv,dim=1)
    v2_vh=mean(v2_vh,dim=1)
endif

;;Second cal step:
v2_hv_cal=sign(v2_hv)*sign(v2_vh)*(abs(v2_hv)/abs(v2_vh))^0.25 ;; ^0.25 converts to visibility from V^4.
;;The sign and abs terms above preserve the sign of the negative visibilities, which are needed to avoid biasing noisy data.

;;If we want uncertainties, calculate the standard deviation
if keyword_set(uncertainties) then begin

    ;;if we arent averaging the uncertainties, then we need to reshape the samples to remove the first dimension
    if not keyword_set(average_together) then begin
        sz=size(v2_hv_cal)
        v2_hv_cal=reform(v2_hv_cal,[long(sz[1])*sz[2],sz[3]])
    endif

    v2_hv_var=variance(v2_hv_cal,dim=1)
    v2_hv_stdev=sqrt(v2_hv_var)
    ;;and then return it as v2_hv_cal
    v2_hv_cal=v2_hv_stdev
endif

return,v2_hv_cal
end
