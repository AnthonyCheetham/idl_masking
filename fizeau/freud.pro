;######## RELEASE CLIP HERE ######## 
;
; freud (head)
;
; freud will read the header passed to it, returning a structure with
; useful information extracted

; Input:  header     - the fits header file to be searched
; Returns: output     - a structure with everything you wanted to know
; 
; created                                                    PGT 01Dec03
; hacked-up to be simpler for LWS                            MJI 11Jan06

function freud, head, olog = olog,  ix = ix

; *** not tested yet ***

s='' &  i=0 &  f=float(0.0) &  d=double(0.0)

headinfo={instrument:s,nax1:i,nax2:i,nax3:i,t_int:f,coadd:i,filter:s,slit:s,optic_cfg:s, $
          lyot:s,grism:s,source_name:s,utc:s,date:s,jd:d,elevation:f,airmass:f,pa:f, $
          ra:d,dec:d,equinox:f, mask:s,  raoff:f, decoff:f,  del_elev:f,  del_pa:f}
;Allow a dummy headinfo to be returned.
if (head[0] eq '') then return,  headinfo

; Firstly figure out what sort of data we have ... 
iname=sxpar(head,'INSTRUME',count=count)
if(count eq 0) then iname=sxpar(head,'CURRINST',count=count)  ; this is for NIRC2

; NOTE need to add ra, dec, equinox to all but nirc2

if(strpos(iname,'NIRC') ne -1 and strpos(iname,'NIRC2') eq -1) then begin
  headinfo.instrument='NIRC'
  headinfo.nax1=sxpar(head,'NAXIS1')
  headinfo.nax2=sxpar(head,'NAXIS2')
  headinfo.nax3=sxpar(head,'NAXIS3')
  headinfo.t_int=sxpar(head,'TINT')
  headinfo.coadd=sxpar(head,'M56COAD0')
  headinfo.filter=strtrim(sxpar(head,'FILTER'),2)
  headinfo.slit=strtrim(sxpar(head,'SLTNAME'),2)
  headinfo.source_name=strtrim(sxpar(head,'TARGNAME'),2)
  headinfo.utc=strtrim(sxpar(head,'UTC'),2)
  headinfo.date=strtrim(sxpar(head,'DATE-OBS'),2)
  headinfo.jd=sxpar(head,'MJD-OBS')
  headinfo.elevation=sxpar(head,'EL')
  headinfo.airmass=sxpar(head,'AIRMASS')
  if(headinfo.jd lt 50500.) then $       ; This selects Jan97 & before where
     headinfo.pa=sxpar(head,'PARANG') $  ; U242N/U379N had different header keyword
     else headinfo.pa=sxpar(head,'PARANTEL')
endif

if(strpos(iname,'LWS') ne -1) then begin
  ;Straight to olog...
  olog.instrument[ix] = 'LWS'
  olog.nax1[ix]=sxpar(head,'NAXIS1')
  olog.nax2[ix]=sxpar(head,'NAXIS2')
  olog.nax3[ix]=sxpar(head,'NAXIS3')
  olog.t_int[ix]=sxpar(head,'FRMTIME')*sxpar(head,'FRMCOADD')
  olog.coadd[ix]=sxpar(head,'FRMCOADD') ;add 'CHPCOADD' ?
  olog.filter[ix]=strtrim(sxpar(head,'FILNAME'),2)
  olog.optic_cfg[ix]=strtrim(sxpar(head,'GRANAME'),2)
  olog.source_name[ix]=strtrim(sxpar(head,'TARGNAME'),2)
  olog.utc[ix]=strtrim(sxpar(head,'UTC'),2)
  olog.date[ix]=strtrim(sxpar(head,'DATE-OBS'),2)
  olog.jd[ix]=sxpar(head,'MJD-OBS')
  olog.elevation[ix]=sxpar(head,'EL')
  olog.airmass[ix]=sxpar(head,'AIRMASS')
  olog.pa[ix]= (sxpar(head,'PARANTEL') + sxpar(head,'ROTPPOSN') - 69.0) mod 360.0
  headinfo.instrument = 'LWS'
endif

if(strpos(iname,'TReCS') ne -1) then begin
  ;Straight to olog...
  ;!!! We need WINDOW for the pixel scale, and both filter wheels...
  headinfo.instrument = 'TRECS'
  olog.instrument[ix] = 'TRECS'
  olog.nax1[ix]=sxpar(head,'NAXIS1')
  olog.nax2[ix]=sxpar(head,'NAXIS2')
  olog.nax3[ix]=sxpar(head,'NAXIS3')
  olog.source_name[ix]=strtrim(sxpar(head,'OBJECT'),2)
  olog.utc[ix]=strtrim(sxpar(head,'TIME-OBS'),2)
  olog.date[ix]=strtrim(sxpar(head,'DATE-OBS'),2)
  olog.filter[ix]=strtrim(sxpar(head,'FILTER1'),2)
  olog.elevation[ix]=sxpar(head,'ELEVATIO')
  olog.airmass[ix]=sxpar(head,'AIRMASS')
  olog.t_int[ix]=sxpar(head,'OBJTIME')
  olog.coadd[ix]=sxpar(head,'FRMCOADD') ;add 'CHPCOADD' ?
  olog.jd[ix]=sxpar(head,'MJD-OBS')
  az =  sxpar(head,'AZIMUTH')
 ; lat =  -(30.+14./60.+16.8/3600.)
 ; altaz2hadec_rot,  olog.elevation[ix],  az,  lat,  ha,  dec,  par
 ; olog.pa[ix]= sxpar(head,'PA') + par ;!!! sign etc not done here !!!
  getrot,head,rot,cdelt
  print,  cdelt
  olog.pa[ix]= -rot
  maskname =  sxpar(head,'LYOT')
  case maskname of  
    'Ciardi  ' :  olog.mask =  't7'
     else     :  begin
     	    	    print,"NO MASK IN DATA"
		    stop
    	         endcase
  endcase
  olog.optic_cfg[ix]=strtrim(sxpar(head,'FILTER2'),2) ;!!! Check this...
endif

if(strpos(iname,'PHARO') ne -1) then begin
  headinfo.instrument  = 'PHARO'
  olog.instrument[ix]  = 'PHARO'
  olog.ra[ix]          = strtrim(sxpar(head,'CRVAL1'),2)
  olog.dec[ix]         = strtrim(sxpar(head,'CRVAL2'),2)
  olog.ha[ix]          = strtrim(sxpar(head,  'HOURANGL'), 2)
  olog.filter[ix]      = strtrim(sxpar(head, 'FILTER'),2)
  olog.slit[ix]        = strtrim(sxpar(head, 'SLIT'),2)
  olog.optic_cfg[ix]   = strtrim(sxpar(head, 'CAROUSEL'),2)
  olog.lyot[ix]        = strtrim(sxpar(head, 'LYOT'),2)
  olog.grism[ix]       = strtrim(sxpar(head, 'GRISM'),2)
  olog.source_name[ix] = strtrim(sxpar(head, 'OBJECT'),2)
  olog.utc[ix]         = strtrim(sxpar(head, 'TIME-OBS'),2)
  olog.date[ix]        = strtrim(sxpar(head, 'DATE-OBS'),2)
  olog.nax1[ix]        = sxpar(head, 'NAXIS1')
  olog.nax2[ix]        = sxpar(head, 'NAXIS2')
  olog.nax3[ix]        = sxpar(head, 'NAXIS3')
  olog.equinox[ix]     = sxpar(head, 'EQUINOX')
  olog.airmass[ix]     = sxpar(head, 'AIR_MASS')
  olog.t_int[ix]       = sxpar(head, 'T_INT') * 1e-3
  ;From Stanmire Metchev's document:
  ; www.astro.caltech.edu/palomar/200inch/palao/Pharo/pharo_plate_scale.pdf
  ; and Mike's checks with Xi Cep from the first run.
;  olog.pa[ix]          = sxpar(head, 'CR_ANGLE') + 25.45
  olog.pa[ix]          = sxpar(head, 'CR_ANGLE') + 25.08 ;New Metchev value
  reads,strmid(olog.date[0],0,4),thisyear               ; Is it after P3K install?
  if(thisyear ge 2012.) then begin 
     olog.pa[ix]          = (sxpar(head, 'CR_ANGLE')  -140.) mod 360 ;New P3K value BUT NOTE MUST ALSO FLIP X !!!
  endif 
  olog.coadd[ix]       = 1
  olog.raoff[ix]       =  sxpar(head, 'RA_OFFS')
  olog.decoff[ix]       =  sxpar(head, 'DEC_OFFS')
;  headinfo.raoff = sxpar(head, 'RAOFF')
;  headinfo.decoff = sxpar(head, 'DECOFF')

  maskname = sxpar(head, 'LYOT')
  case maskname of
      '18H     ' : olog.mask = 'p18'
      '9H - S2 ' : olog.mask = 'p9'
      '9H      ' : olog.mask = 'p9'
      else       : olog.mask = ''
  endcase
endif

if(strpos(iname,'NIRC2') ne -1) then begin
  headinfo.instrument='NIRC2'
  headinfo.nax1=sxpar(head,'NAXIS1')
  headinfo.nax2=sxpar(head,'NAXIS2')
  ; headinfo.nax3=sxpar(head,'NAXIS3')
  headinfo.t_int=sxpar(head,'ITIME')
  headinfo.coadd=sxpar(head,'COADDS')
  headinfo.filter=strtrim(sxpar(head,'FILTER'),2)
  headinfo.slit=strtrim(sxpar(head,'SLITNAME'),2)
  headinfo.lyot=strtrim(sxpar(head,'PMSNAME'),2)
  headinfo.optic_cfg=strtrim(sxpar(head,'CAMNAME'),2)
  headinfo.source_name=strtrim(sxpar(head,'TARGNAME'),2)
  headinfo.utc=strtrim(sxpar(head,'UTC'),2)
  headinfo.date=strtrim(sxpar(head,'DATE-OBS'),2)
  headinfo.jd=sxpar(head,'MJD-OBS')
  headinfo.elevation=sxpar(head,'EL')
  headinfo.airmass=sxpar(head,'AIRMASS')
  headinfo.pa=360.+sxpar(head,'PARANG')+sxpar(head,'ROTPPOSN') $
                  -sxpar(head,'EL')-sxpar(head,'INSTANGL')   
                                      ; THIS formula not checked yet!
  headinfo.ra=sxpar(head,'RA')
  headinfo.dec=sxpar(head,'DEC')
  headinfo.equinox=sxpar(head,'EQUINOX')
  headinfo.raoff = sxpar(head, 'RAOFF')
  headinfo.decoff = sxpar(head, 'DECOFF')
endif

; a bit hacked for now. Need more intelligence to work out filters (won't do for J)
if(strpos(iname,'CONICA') ne -1) then begin
  headinfo.instrument='CONICA'
  headinfo.nax1=sxpar_conica(head,'NAXIS1')
  headinfo.nax2=sxpar_conica(head,'NAXIS2')
  nax=sxpar_conica(head,'NAXIS')
  if(nax gt 2) then headinfo.nax3=sxpar_conica(head,'NAXIS3') 
  headinfo.t_int=sxpar_conica(head,'SODETDIT') ; HIERARCH ESO DET DIT
  ;headinfo.coadd=sxpar_conica(head,'COADDS')
  ; "filter" is tricky. Usually one of two wheels, but I have not yet covered case for J!
  wheel5=strtrim(sxpar_conica(head,'SOPTI5ID'),2) ; HIERARCH ESO INS OPTI6 ID
  wheel6=strtrim(sxpar_conica(head,'SOPTI6ID'),2) ; HIERARCH ESO INS OPTI6 ID
  if( (strpos(wheel5,"empty"))[0] eq -1) then headinfo.filter=wheel5 else headinfo.filter=wheel6
  ;headinfo.slit=strtrim(sxpar_conica(head,'SLITNAME'),2)
  headinfo.mask=strtrim(sxpar_conica(head,'SOPTI3ID'),2) ; HIERARCH ESO INS OPTI3 ID
  headinfo.lyot=strtrim(sxpar_conica(head,'SOPTI7ID'),2) ; ok, *NOT* lyot - this is the camera S13 or L27
  headinfo.optic_cfg=strtrim(sxpar_conica(head,'SOPTI4ID'),2) ; HIERARCH ESO INS OPTI4 ID - eg 'Wollaston_45'
  headinfo.source_name=strtrim(sxpar_conica(head,'OBJECT'),2) 
  if( (strpos(headinfo.source_name,"name not set"))[0] ne -1) then $
          headinfo.source_name=strtrim(sxpar_conica(head,'OOBSNAME'),2)
  headinfo.utc=strtrim(sxpar_conica(head,'UTC'),2)
  headinfo.date=strtrim(sxpar_conica(head,'DATE-OBS'),2)
  headinfo.jd=sxpar_conica(head,'MJD-OBS')
  headinfo.elevation=sxpar_conica(head,'EL')
  headinfo.airmass=sxpar_conica(head,'AIRMASS')
  ; Now work out sky PA of data cube. Not completely straightforward... why no start/end for Alt??
  rotstart=sxpar_conica(head,'ROTSTART')    ; HIERARCH ESO ADA ABSROT START
  rotend  =sxpar_conica(head,'BSROTEND')    ; HIERARCH ESO ADA ABSROT END
  pastart =sxpar_conica(head,'ANGSTART')    ; HIERARCH ESO TEL PARANG START
  paend   =sxpar_conica(head,'ARANGEND')    ; HIERARCH ESO TEL PARANG END
  alt     =sxpar_conica(head,'SOTELALT')    ; HIERARCH ESO TEL ALT
  instrument_offset= -0.55
  headinfo.pa=(rotstart+rotend)/2.+alt-(180.-(pastart+paend)/2.) + instrument_offset
  headinfo.ra=sxpar_conica(head,'RA')
  headinfo.dec=sxpar_conica(head,'DEC')
  headinfo.equinox=sxpar_conica(head,'GEQUINOX') ; HIERARCH ESO TEL TARG EQUINOX
  ;headinfo.raoff = sxpar_conica(head, 'RAOFF')
  ;headinfo.decoff = sxpar_conica(head, 'DECOFF')
endif

if (strpos(iname,'GPI') ne -1) then begin
   headinfo.instrument='GPI'
   olog.instrument[ix]='GPI'
   olog.nax1[ix]=sxpar(head,'NAXIS1')
   olog.nax2[ix]=sxpar(head,'NAXIS2')
   olog.nax3[ix]=sxpar(head,'NAXIS3')
   olog.t_int[ix]=sxpar(head,'ITIME0')*1e-6
   olog.coadd[ix]=sxpar(head,'COADDS')
   filt=strtrim(sxpar(head,'OBSMODE'),2)
   if filt eq 'NRM_Y' then filter='Y'
   if filt eq 'NRM_H' or filt eq 'H_direct' then filter='H'
   if filt eq 'NRM_J' then filter='J'
   if filt eq 'NRM_K1' then filter='K1'
   if filt eq 'NRM_K2' then filter='K2'
   if filt eq 'Y_coron' then filter='Y'

   olog.filter[ix]=filter
   olog.slit[ix]=''                 ;;no slit?
   olog.optic_cfg[ix]=sxpar(head,'')  ;;??what do we need?
   olog.lyot[ix]=sxpar(head,'LYOTMASK')
   olog.grism[ix]=sxpar(head,'')      ;;??
   olog.source_name[ix]=sxpar(head,'OBJECT') 
   olog.utc[ix]=sxpar(head,'UTSTART')
   olog.date[ix]=sxpar(head,'DATE-OBS')
   olog.jd[ix]=sxpar(head,'MJD-OBS')
   olog.elevation[ix]=sxpar(head,'ELEVATIO')  ;;Not important?
   olog.airmass[ix]=sxpar(head,'')    ;;Not important?
   
   ;;After talking to Marshall Perrin, the position angle is complicated. Par_ang is measured from the ifs rotation. i.e. 90-24.5 degrees from up.
   crpa=sxpar(head,'CRPA') ;;cass rotator pos ang. Should be near zero!
   ipa=sxpar(head,'PA')  ;;GPIs pos ang. Should be fixed. 
   iaa=sxpar(head,'IAA');;instrument alignment angle for when it is on backwards or on side port.
   par_ang=sxpar(head,'AVPARANG') ;;parralactic angle
   lenslet_ang=90.-23.5 ;;this was taken from the gpi pipeline in 2014. This gets north to be up.
   gpi_north_offset=-1.00 ;;degrees, +/- 0.03
   ;;this is what I now believe the formula to be, but needs checking (Jul 2014)
   ; pa=par_ang + crpa - lenslet_ang -gpi_north_offset

   ;; Alternatively, getrot does all this for us!
   getrot,head,rot,/silent
   olog.pa[ix]=(720 + rot+gpi_north_offset) mod 360
   ; print,'Hack in freud. Need to use getrot instead of this manual formula.'
   ; olog.pa[ix]=(720 + pa) mod 360
   olog.ra[ix]=sxpar(head,'CRVAL1')
   olog.dec[ix]=sxpar(head,'CRVAL2')
   olog.equinox[ix]=2000.0      
   ;;Now, what mask:
   mask=sxpar(head,'APODIZER')
   if mask eq 'APOD_NRM_G6208' then mask='g10s40' else mask='None'
   olog.mask=mask
   olog.hwp[ix]=sxpar(head,'WPANGLE')
   olog.del_pa[ix]=sxpar(head,'')    ;; This is calculated in qbe_gpi
   olog.del_elev[ix]=sxpar(head,'')  ;; This is calculated in qbe_gpi
   ;; the following aren't used anywhere in the pipeline, so ACC didn't add them
   olog.raoff[ix]=0
   olog.decoff[ix]=0
endif

if (strpos(sxpar(head,'TELESCOP'),'simu') ne -1) then begin
   headinfo.instrument='simu'
   olog.instrument[ix]='simu'
   olog.nax1[ix]=sxpar(head,'NAXIS1')
   olog.nax2[ix]=sxpar(head,'NAXIS2')
   olog.nax3[ix]=sxpar(head,'NAXIS3')
   olog.t_int[ix]=sxpar(head,'TINT')*1e-6
   olog.coadd[ix]=sxpar(head,'COADDS')
   olog.filter[ix]=sxpar(head,'filter')
   olog.source_name[ix]='Fake'

   olog.utc[ix]=sxpar(head,'OTIME')
   olog.date[ix]=sxpar(head,'ODATE')
;   olog.jd[ix]=sxpar(head,'MJD-OBS')

   olog.pa[ix]=sxpar(head,'orient')         ;;Is that right?
   olog.ra[ix]=sxpar(head,'CRVAL1')
   olog.dec[ix]=sxpar(head,'CRVAL2')
   olog.equinox[ix]=2000.0      ;sxpar(head,'')  ;; Should be in header somewhere...
;;Now, what mask:
   mask='golay9';sxpar(head,'APODIZER')
   olog.mask=mask
endif

if (strpos(iname,'SPHERE') ne -1) then begin
  ;;I'm hoping this is the same for IRDIS, IFS and ZIMPOL.
  headinfo.instrument='SPHERE'
  olog.instrument[ix]='SPHERE'
  olog.nax1[ix]=sxpar_conica(head,'NAXIS1')
  olog.nax2[ix]=sxpar_conica(head,'NAXIS2')
  nax=sxpar_conica(head,'NAXIS')
  if(nax gt 2) then olog.nax3[ix]=sxpar_conica(head,'NAXIS3')
  olog.t_int[ix]=sxpar_conica(head,'TSEQ1DIT')

  ;; This is going to get messy. SPHERE has way too many filter and 
  ;;  instrument combinations, and it saves the setup of all 3 cameras in every header.
  ;; These are the locations of the mask, filter1 and filter2 for each camera
  ;   irdis_par=['1OPTI1NAME','1OPTI2NAME','S1FILTNAME']
  ;   ifs_par=['OPTI14NAME','2OPTI2NAME'] ;(only 1 filter)
  ;   zimpol_par=['4OPTI9NAME','3OPTI5NAME','3OPTI6NAME']

  ;; Work out what camera it is
  camera=sxpar_conica(head,'ESODETID',nchar=8) ;; IRDIS and IFS have it here
  if typename(camera) ne typename('string') then camera=sxpar_conica(head,'ETDEV1ID') ;; ZIMPOL has it here
  camera=strcompress(camera,/remove_all)
  olog.optic_cfg[ix]=camera

  if camera eq 'IFS' then begin
    ;; !!!!!!!!!!!
    ;; !!  IFS
    ;; !!!!!!!!!!!
    filter=strsplit(sxpar_conica(head,'2OPTI2NAME',nchar=10),'PRI_',/regex,/extract)
    mask= strcompress(sxpar_conica(head,'OPTI14NAME',nchar=10),/remove_all)

    north_offset=0.;-1.764 ;;(+/- 0.01 ) from SVT, in manual. USE ZERO AND CORRECT LATER
    instrument_offset=-100.46;;deg (+/- 0.13) from manual
    pupil_offset=135.87;;deg (+/- 0.03) from manual
  endif else if camera eq 'IRDIS' then begin
    ;; !!!!!!!!!!!
    ;; !!  IRDIS
    ;; !!!!!!!!!!!
    filt1=strcompress(sxpar_conica(head,'1OPTI2NAME',nchar=10),/remove_all)
    filt2=strcompress(sxpar_conica(head,'S1FILTNAME',nchar=10),/remove_all)
    mask= strcompress(sxpar_conica(head,'INS1OPTI1TYPE',nchar=13),/remove_all)

    if filt1 eq 'P0-90' then filter='DPI_'+filt2 else $
    if filt1 eq 'CLEAR' then filter='CI_'+filt2 else filter=filt1

    north_offset=0.;-1.764 ;;(+/- 0.01 ) from SVT, in manual. USE ZERO AND CORRECT LATER
    instrument_offset=0. ;; from manual (pupil offset takes this into account already)
    pupil_offset=135.87;;deg (+/- 0.03) from manual

  endif else if camera eq 'ZIMPOL' then begin
    ;; !!!!!!!!!!!
    ;; !!  ZIMPOL
    ;; !!!!!!!!!!!
    mask=strcompress(sxpar_conica(head,'4OPTI19NAME',nchar=11),/remove_all)
    filt1=strcompress(sxpar_conica(head,'3OPTI5NAME',nchar=10),/remove_all)
    filt2=strcompress(sxpar_conica(head,'3OPTI6NAME',nchar=10),/remove_all)
    obstype=strcompress(sxpar_conica(head,'DPRTECH',nchar=7),/remove_all)

    if mask eq 'OPEN' then begin
      print,'  WARNING: ZIMPOL says there is no SAM mask...'
      mask='FREE'
    endif

    if strpos(obstype,'IMAGE') ne -1 then begin
      ;; ZIMPOL is in CI mode!
      filter=filt1+'_'+filt2
    endif else if strpos(obstype,'POLARIMETRY') ne -1 then begin
      ;; ZIMPOL is in Pol mode!
      filter=filt1+'_'+filt2;; +'_pol' ;; However the _pol suffix is added in inquire_sphere.pro

    endif else begin
      ;; 
      print,'ZIMPOL is using a weird mode: ', obstype
      stop
    endelse
    
    ;; What are these values? They are not in the manual...
    north_offset=0.
    instrument_offset=0.
    pupil_offset=0.

  endif else begin
    print,'ERROR: Unknown SPHERE camera detected!'
    stop
  endelse
  filter=strcompress(filter,/remove_all)

  olog.filter[ix]=filter
  case mask of
   'MASK_SAM': olog.mask='7Hole'
   'V_SAM7' : olog.mask='7Hole' ;; Unfortunately that's its name in ZIMPOL at the moment
   'ST_SAM': olog.mask='7Hole'
   ELSE: stop ;; the mask is not known
  endcase
  olog.lyot[ix]=strtrim(sxpar_conica(head,'ILT2NAME'),2) ; This is *NOT* the Lyot, but actually the ND filter
  olog.source_name[ix]=strtrim(sxpar_conica(head,'OBJECT'),2) 
  if( (strpos(olog.source_name,"name not set"))[0] ne -1) then $
          olog.source_name=strtrim(sxpar_conica(head,'OOBSNAME'),2)
  date=strsplit(sxpar_conica(head,'DATE-OBS'),'T',/extract)
  olog.utc[ix]=date[1]
  olog.date[ix]=date[0]
  olog.jd[ix]=sxpar_conica(head,'MJD-OBS')
  
  olog.airmass[ix]=(sxpar_conica(head,'LAIRMEND')+sxpar_conica(head,'IRMSTART'))/2

  ;; Half wave plates / Pol keywords
  zimp_hwp=sxpar_conica(head,'3DROT2POSANG',nchar=12) ;; ZIMPOL HWP
  cpi_hwp1=sxpar_conica(head,'4DROT1POSANG',nchar=12) ;; Common path HWP1
  cpi_hwp2=sxpar_conica(head,'4DROT3POSANG',nchar=12) ;; Common path HWP2
  stokes=strcompress(sxpar_conica(head,'STOKES',nchar=6),/remove_all) ;; This should work for IRDIS and ZIMPOL
  if stokes eq 'Qplus' then hwp=0 $
  else if stokes eq 'Qminus' then hwp=45 $
  else if stokes eq 'Uplus' then hwp=22.5 $
  else if stokes eq 'Uminus' then hwp=67.5 $
  else hwp=0
  olog.hwp[ix]=hwp

  olog.elevation[ix]='';;not important?
  olog.grism[ix]='';not important
  olog.slit[ix]='';not important

  ; Now work out sky PA of data cube. Not completely straightforward... why no start/end for Alt??
  rotstart=sxpar_conica(head,'ROTSTART')    ; HIERARCH ESO ADA ABSROT START
  rotend  =sxpar_conica(head,'BSROTEND')    ; HIERARCH ESO ADA ABSROT END
  pastart =sxpar_conica(head,'ANGSTART')    ; HIERARCH ESO TEL PARANG START
  paend   =sxpar_conica(head,'ARANGEND')    ; HIERARCH ESO TEL PARANG END
  alt     =sxpar_conica(head,'SOTELALT')    ; HIERARCH ESO TEL ALT
  ; This is the equation for NACO. Is it the same?
  ;olog.pa[ix]=(rotstart+rotend)/2.+alt-(180.-(pastart+paend)/2.) + pupil_offset + instrument_offset + north_offset
  ;; Here is what the SPHERE manual suggests, which gives much better results on HD 142527.
  ; olog.pa[ix]=(pastart+paend)/2. + pupil_offset + instrument_offset + north_offset

  ;; Actually, we want to use the RA and DEC values from the header to calculate the parang,
  ;; since SPHERE had some problems with the header values.
  ;; This calculation is complicated, so it has its own routine:

  olog.pa[ix] = calc_sphere_parang(head) + pupil_offset + instrument_offset + north_offset
  print,'Using True north offset of:',north_offset,'deg. PA:',olog.pa[ix]

  olog.ra[ix]=sxpar_conica(head,'RA')
  olog.dec[ix]=sxpar_conica(head,'DEC')
  olog.equinox[ix]=sxpar_conica(head,'GEQUINOX') ; HIERARCH ESO TEL TARG EQUINOX
endif

if(headinfo.instrument eq '') then $
  print,'*** ERROR - Freud could not understand the Header'
return,headinfo

end
 

