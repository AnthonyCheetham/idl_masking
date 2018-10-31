;######## RELEASE CLIP HERE ######## 
;+
; inquire(paramname,infostr,ix=ix)
;
; inquire will return a default setting (for example, a template) by
;   using the observing log info structure.

; Input:  name         - the name of the thing you want
;         ix           - index for addressing infostr(if needed)
; Returns: infostr     - a structure with everything you wanted to know
; 
; created                                                    PGT 22Oct05
;
; This is really example code to show how to build it at
; present
function inquire_simu,paramname,infostr,ix

camera=infostr.instrument[ix]
case paramname of
   'template': begin
      mask_xtra=''
      datasize=fix(infostr.cube_sz[ix,0])
      mf_name='simu/mf_golay9_kp.idlvar'
      return,mf_name
   end
   'clog': return,make_clog(infostr)
   'n_blocks':return,0 ;;??
   'pscale':stop        ;return,14.0;; PIXEL SCALE
   'mask':return,infostr.mask
endcase
return,0

end
 

