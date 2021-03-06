;$Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $
;$Log: inquire.pro,v $
;Revision 1.7  2010/04/01 04:04:00  gekko
;Added release clip line
;
;Revision 1.6  2009/03/03 07:11:12  snert
;MJI commiting stuff for Conica.
;
;Revision 1.5  2008/05/21 23:11:28  mireland
;Not quite sure what all these changes are - but now they are commmited anyway.
;
;Revision 1.4  2007/07/23 08:01:45  gekko
;add conica PGT
;
;Revision 1.3  2006/12/04 06:03:35  mireland
;Added T-ReCS functionality to fizeau.
;
;Revision 1.2  2005/12/20 21:51:55  mireland
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and $Log: inquire.pro,v $
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and Revision 1.7  2010/04/01 04:04:00  gekko
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and Added release clip line
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and Revision 1.6  2009/03/03 07:11:12  snert
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and MJI commiting stuff for Conica.
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and Revision 1.5  2008/05/21 23:11:28  mireland
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and Not quite sure what all these changes are - but now they are commmited anyway.
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and Revision 1.4  2007/07/23 08:01:45  gekko
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and add conica PGT
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and Revision 1.3  2006/12/04 06:03:35  mireland
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and Added T-ReCS functionality to fizeau.
;Added $Id: inquire.pro,v 1.7 2010/04/01 04:04:00 gekko Exp $ and to important files, and added the g18 mask
;for PHARO's inquire.
;
;######## RELEASE CLIP HERE ######## 
; inquire(paramname,infostr,ix=ix)
;
; inquire will return a default setting (for example, a template) by
;   using the observing log info structure.

; Input:  name         - the name of the thing you want
;         ix           - index for addressing infostr(if needed)
; Returns: infostr     - a structure with everything you wanted to know
; 
; created                                                    PGT 22Oct05

function inquire,name,infostr,ix=ix

if(keyword_set(ix) eq 0) then ix=0

camera=infostr.instrument[ix]
camera=strupcase(camera)

case camera of

  'NIRC2': return, inquire_nirc2(name, infostr, ix)
  'NIRC':  return, inquire_nirc (name, infostr, ix)
  'LWS':   return, inquire_lws  (name, infostr, ix)
  'PHARO': return, inquire_pharo(name, infostr, ix)
  'TRECS': return, inquire_trecs(name, infostr, ix)
  'CONICA': return, inquire_conica(name, infostr, ix)
  'JWST':   return,  inquire_jwst(name,  infostr,  ix)
  'LAMP': return, inquire_lamp(name,infostr,ix)
  'MMT':   return,  inquire_mmt(name,  infostr,  ix)
  'GPI': return,inquire_gpi(name,infostr,ix)
  'SIMU': return,inquire_simu(name,infostr,ix)
  'SPHERE': return,inquire_sphere(name,infostr,ix)
; etc (add more...)

endcase


end

