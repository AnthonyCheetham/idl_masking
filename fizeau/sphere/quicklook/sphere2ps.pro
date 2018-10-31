; Quick script to quicklook at data
; PGT 2010
; edited for SPHERE/IRDIS: ACC 2015
;
; This is a very quick-and-dirty quicklook eavesdrop.
; For the moment, just analyze the SUMMED frame at the end of the cube
; caldir = directory with bad pix and flats
; datadir = directory with data

;; This is just for IRDIS!

function sphere2ps,filename,caldir=caldir,datadir=datadir

; Directory containing bad pix and flats
if not keyword_set(caldir) then caldir='./'

if not keyword_set(datadir) then datadir='./'

print,''
print,'Reading File: ',filename
incube=float(reform(readfits(datadir+filename,head)))
info=size(incube)
dimx=info[1]
dimy=info[2]
nframes=info[3]

;;get the header info so we know what we're dealing with
hinfo=freud(head)

stop

; sum the frames
sumframe=total(incube,3)


; Formula for position angle taken from /import/spiral1/snert/code/masking/fizeau/freud.pro
rotstart = sxpar_conica(head,'ROTSTART')
rotend   = sxpar_conica(head,'BSROTEND')
pastart  = sxpar_conica(head,'ANGSTART')
paend    = sxpar_conica(head,'ARANGEND')
alt      = sxpar_conica(head,'SOTELALT')
instrument_offset = -0.55
pos_ang = (rotstart+rotend)/2.+alt-(180.-(pastart+paend)/2.) + instrument_offset
pos_ang_lessPA = (rotstart+rotend)/2.+alt + instrument_offset

; These are approximate rotations for given mask
if(hinfo.mask eq  '7Holes') then maskrot=104.
if(hinfo.mask eq  '9Holes') then maskrot=47.5
if(hinfo.mask eq  '18Holes') then maskrot=78.0

tolerance=1.5  ; scream if mask rot is more than this out!
rotdiff=(pos_ang_lessPA - maskrot)
if(abs(rotdiff) lt tolerance) then $
   print,'Pupil Rotation =',pos_ang_lessPA,' ... looks OK' $
   else print,'Pupil Rotation =',pos_ang_lessPA,' ##### WARNING - THIS IS NOT STANDARD #####' 



; Now we need to find an appropriate flat file ...

  ssky=fltarr(dimx,dimy)    ; set this to zero in case we can't get a sky
  if(dimx eq 1024) then begin
     restore,sdir+'flat1024.idlvar'
     subsz=512
     ;medskylv=median(incube[4:30,*,*])
     medskylv=median(sumframe[4:30,*])
  endif
  if(dimx eq 512) then begin
     restore,sdir+'flat512.idlvar'
     subsz=256
     ;medskylv=median(incube[4:30,*,*])
     medskylv=median(sumframe[4:30,*])
     if(hinfo.filter eq 'L_prime') then begin
        restore,sdir+'flat512_L.idlvar'
        restore,sdir+'sky512_L.idlvar'    ;restores ssky
      endif
  endif
  if(dimx eq 256) then begin
     restore,sdir+'flat256.idlvar'
     subsz=128
     ;medskylv=median(incube[2:15,*,*])
     medskylv=median(sumframe[2:15,*])
     if(hinfo.filter eq 'L_prime') then begin
        restore,sdir+'flat256_L.idlvar'
        restore,sdir+'sky256_L.idlvar'    ;restores ssky
     endif
     if(hinfo.filter eq 'M_prime') then begin
        restore,sdir+'flat256_L.idlvar'
        restore,sdir+'sky256_M.idlvar'    ;restores ssky
     endif
     if(hinfo.filter eq 'NB_3.74') then begin
        restore,sdir+'flat256_L.idlvar'
        restore,sdir+'sky256_NB_3.74.idlvar'    ;restores ssky
     endif
     if(hinfo.filter eq 'NB_4.05') then begin
        restore,sdir+'flat256_L.idlvar'
        restore,sdir+'sky256_NB_4.05.idlvar'    ;restores ssky
     endif
  endif


  ;median of supersky 
  medssky=median(ssky[2:15,*])
  ;medssky=median(ssky)
  if(abs(medssky) gt 1e-2) then $  ; i.e. we have a real ssky
     ssky = ssky*medskylv/medssky  ; scale the ssky to the actual flux

  pspec=fltarr(subsz,subsz)

  ; subtract a scaled supersky if we have one
  speck=sumframe-ssky
  speck=speck/flat
  
  ; Get rid of bad pixels
  ;speck(bad_pixels)=median(speck)
  speck=fix_bad_pixels(speck,bad=bad_pixels)
  speck = sigma_filter_nirc2(speck,5,n_sigma=n_sigma,/all,/iterate) ;,/mon)


  wset,1
  !p.multi=[0,2,1]
  image_cont,speck,/nocont,/asp, tit=filename

  ; Shift to center
  ; Trim edges and centerlines for centering...
   dummy=speck
   border=10 & minlv=min(speck)
   dummy[0:border,*]=minlv & dummy[dimx-border:*,*]=minlv
   dummy[*,0:border]=minlv & dummy[*,dimy-border:*]=minlv
   dummy[dimx/2-border:dimx/2+border,*]=minlv
   dummy[*,dimy/2-border:dimy/2+border]=minlv
   dummy=smooth(dummy,fix(subsz/25))
   mx=max(dummy,mix)
   x0= mix mod dimx & y0= mix/dimy
   speck=shift(speck,dimx/4-x0,dimy/4-y0)
   speck_bgr=median(speck[dimx*.7:dimx*.8,dimy*.7:dimy*.8])
   speck=speck[0:subsz-1,0:subsz-1]
   speck=speck-speck_bgr

  ; blank out chip periphery for data analysis
   window=fltarr(subsz,subsz)
   window(border:subsz-border,border:subsz-border)=1
   window=smooth(window,border,/edge)
   speck=speck*window

  print,'Peak Pixel Value =',max(speck),'  for file ',filename
  if(max(speck) gt 10000) then print," *** WARNING - NONLINEAR REGIME ***"

  ; Make up Hanning window ...
   d=shift(dist(subsz,subsz),subsz/2,subsz/2)
   hsize=nint(subsz/1.8)  
   h=hanning(hsize)
   h=shift(h,hsize/2)
   chan=fltarr(subsz,subsz)
   chan(where(d le hsize/2)) = h(round(d(where(d le hsize/2))))

  image_cont,speck,/nocont,/asp
  !p.multi=0

   ; Now form the power spectrum ...
   ft=fft(speck*chan)
   ps=(abs(ft))^2
   wset,0
   image_cont,alog10(shift(ps,subsz/2,subsz/2)),/nocont,/asp, tit='log power spec'

outdata=ps 

return,outdata

end
