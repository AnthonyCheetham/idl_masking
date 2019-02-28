; JDM
;  ri2at.pro
;
;   REAL/Imaginary to Amplitude/Phase
;
; 2011Jun17 JDM changed variables.
; 	seemed conflict with ISAP library....

pro ri2at,jdm_re,jdm_im,jdm_am,jdm_ph,help=help

if (keyword_set(help) eq 1) then begin
  print,' REAL/Imaginary to Amplitude/Phase (phase in Degrees)
 print,'pro ri2at,re,im,am,ph
 return
endif
pi=3.14159265358979323846264338327950288419716939937d
;stop
jdm_ph=jdm_re
jdm_ph(*)=0.0
jdm_am=abs(complex(jdm_re,jdm_im))
in=where(jdm_re ne 0,ct)
if (ct ne 0) then $
  jdm_ph(in)=atan(jdm_im(in)/jdm_re(in))*180.d /pi
in=where(jdm_re lt 0,count)
if (count gt 0) then $
  jdm_ph(in)=jdm_ph(in) + 180.d
jdm_ph=(jdm_ph+360.0000d) mod 360.d

index=where(jdm_re eq 0 and jdm_im lt 0.0,ct) 
if (ct gt 0) then begin
  jdm_ph(index)=270.d
  jdm_am(index)=abs(jdm_im(index))
endif
index=where(jdm_re eq 0 and jdm_im gt 0.0,ct)
if (ct gt 0) then begin
  jdm_ph(index)=90d
  jdm_am(index)=abs(jdm_im(index))
endif


end


