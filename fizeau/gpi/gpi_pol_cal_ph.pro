;; this program does the actual calibration of differential phase
;; quantities for gpi pol data. It is called by gpi_pol_double_diff.pro.
;;
;; Setting uncertainties=1 will return the standard deviation of the samples.
;;
;; Setting average=1 will average the v2 after the first difference measurement
;;
;; Author: ACC 2014.

function gpi_pol_cal_ph,ph_h_0,ph_v_0,ph_h_45,ph_v_45,uncertainties=uncertainties,average_together=average_together
;;This function does the phase quantities (like clp and phase)

;;First cal step:
ph_hv=ph_h_0-ph_v_0
ph_vh=ph_h_45-ph_v_45

;;Average them
if keyword_set(average_together) then begin
    ph_hv=mean(ph_hv,dim=1)
    ph_vh=mean(ph_vh,dim=1)
endif

;;Second cal step:
ph_hv_cal=(ph_hv-ph_vh)

;;If we want uncertainties, calculate the standard deviation
if keyword_set(uncertainties) then begin

    ;;if we arent averaging the uncertainties, then we need to reshape the samples to remove the first dimension
    if not keyword_set(average_together) then begin
        sz=size(ph_hv_cal)
        ph_hv_cal=reform(ph_hv_cal,[long(sz[1])*sz[2],sz[3]])
    endif

   ph_hv_var=variance(ph_hv_cal,dim=1)
   ph_hv_stdev=sqrt(ph_hv_var)
   ;;and then return it as ph_hv_cal
   ph_hv_cal=ph_hv_stdev
endif

return,ph_hv_cal
end

