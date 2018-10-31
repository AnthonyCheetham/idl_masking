;;Program for cubing sphere data.
;;Developed by ACC starting in June 2015.
;;The aim is to make one program for IFS, IRDIS and ZIMPOL, though that may change.
;; Inputs:
;;     -data_dir: directory containing the data
;;     -jobname: an identifier to be added to the cubeinfo files etc
;;     -flatpath: the location of the flat field (actually optional)
;;     -frames: the frame numbers of the data
;;     -skies: the frame numbers of the sky files (not implemented!)
;;     -tsize: simultaneously the angular size of the target in each frame, and a 
;;         way to distinguish between targets. Make a unique number for each target.
;; Options:
;;     -prefix/extn: file names must be prefix+frames[ix]+extn
;;     -cutsize: the size of the output frames (default 256)
;;     -no_centre/no_center: don't centre the data frames
;;     -pol: process as polarized data (this probably doesnt need to be a keyword...)
;;     -discard_sigma: settings for flagbad_conica, which rejects frames based on n-sigma 
;;         outliers in: [x_pos,y_pos,total_counts,peak_pixel,fraction of peak median]. 
;;         Set to a negative number to not reject based on that criteria.
;;     -merge_irdis: Merge the irdis channels into 1 cube so they are treated the same. 
;;         Default behaviour is to keep them separate, in case they are different 
;;         wavelengths or have different systematics.
;;     -centre_once_per_cube: If set, each frame in a cube will be centred by the same amount.
;;
;;
;;Todo:
;;    - ZIMPOL Pol
;;    - ZIMPOL imaging
;;    - Speed everything up by cutting the data early on, then centering it later. This may cause some image artefacts...
;; For pol data, each cube should be saved with _0 = horiz, _1=vertical ? (or vice versa? We should talk to the ZIMPOL team to check!)

pro qbe_sphere,data_dir,jobname,flatpath,frames,skies,tsize,prefix=prefix,extn=extn,centre_once_per_cube=centre_once_per_cube,$
    cutsize=cutsize,no_centre=no_centre,no_center=no_center,pol=pol,discard_sigma=discard_sigma,merge_irdis=merge_irdis

if keyword_set(prefix) eq 0 then prefix=''
if keyword_set(extn) eq 0 then extn='.fits'
if keyword_set(cutsize) eq 0 then cutsize=256
;;these numbers inherited from NACO, but they arent good defaults
if not keyword_set(discard_sigma) then discard_sigma=[-2., -2., 3.5,  3.5, 0.4]
if keyword_set(no_center) then no_centre=1
if keyword_set(pol) then pol=1 else pol=0;;this will be used later
;;some hardcoded parameters:
n_sigma_irdis=6    ;; for sigma_filter_nirc2 (cosmic ray removal)
box_width_irdis=5 ;; for sigma_filter_nirc2 (cosmic ray removal)
n_sigma_ifs=6    ;; for sigma_filter_nirc2 (cosmic ray removal)
box_width_ifs=15 ;; for sigma_filter_nirc2. Need a bigger area since the bad pixels come in blocks

n_blocks=n_elements(frames)
if keyword_set(pol) then n_cubes=2*n_blocks else n_cubes=n_blocks

;;first, populate the olog file and flat field the data
s=replicate('',n_cubes) & s2=strarr(n_cubes,2) &  i=replicate(0,n_cubes) &  f=replicate(float(0.0),n_cubes) &  d=replicate(double(0.0),n_cubes) & fi = intarr(2, n_cubes)
;;make two ologs: an actual one, and one to save quantities for every frame so they can be averaged together.
olog={instrument:s,nax1:i,nax2:i,nax3:i,t_int:f,coadd:i,filter:s,   $
      slit:s,optic_cfg:s,lyot:s,grism:s,source_name:s,utc:s,        $
      date:s,jd:d,elevation:f,del_elev:f,airmass:f,pa:f,del_pa:d,   $
      ra:d,dec:d,equinox:f,  ha:d,raoff:d,  decoff:d,rawdir:'',     $
      tsize:tsize,mask:'',frames:frames,skies:skies,uflip:1.0,      $
      cubedir:'', frames_ix:fi,adate:'',comments:'',proc_hist:'',   $
      cal4src:intarr(n_cubes,n_cubes), quad:i,logflag:0,            $
      dk_fname:s2,cube_fname:s2,cube_tsize:f,         $
      cube_sz:fltarr(n_cubes,3),flat_file:strarr(1),hwp:f,          $
      pas:fltarr(n_blocks),hwps:fltarr(n_blocks)}

;; Define data structure to hold all the frame statistics:
f=fltarr(n_cubes,2)     ;; frame statistic and err
fl= fltarr(n_cubes,10,2);; average and errors. src#, aperture size, 0=mean/1=error
;;the 10 on the previous line comes from qbe_conica, so dont blame me!
framestats={xpk:f,ypk:f,totflx:f,pkflx:f,skybgr:f,phot:fl}

;; Prepare the sky frame before we work on the data

if keyword_set(skies) then begin
  n_skies = n_elements(skies)
  print,'Constructing sky frame from',n_skies,'files'
  supersky = 0
  for sky_ix=0,n_skies-1 do begin
    sky = readfits(skies[sky_ix])
    if (size(sky))[0] eq 3 then sky = mean(sky,dim=3)
    supersky += sky
  endfor

  supersky /= n_skies

  ;; Save the sky for diagnostic purposes
  save,supersky,file='sky.idlvar'

endif else print,'No sky frame set! Skipping sky subtraction'

;; loop through data files
for cube_ix=0,n_blocks-1 do begin
  print,'Reading file',cube_ix+1,' of ',n_blocks
  ;;read the data
  in_name=data_dir+prefix+strtrim(string(frames[cube_ix],format="(I4.4)"))+extn
  data_in=readfits(in_name,head,/silent)

  datasz=size(data_in)
  if datasz[0] eq 2 then nframes=1 else nframes=datasz[3]

  ;;collect statistics about each frame (xpeak,ypeak,tot_flux,pk_flux,sky_bg)
  fstats=fltarr(5,nframes)

  ;;get the header
  if keyword_set(pol) then begin
    ;; since the 2 pol channels have different olog indices we need the entries repeated (easiest to run freud again)
    dummy=freud(head,olog=olog,ix=2*cube_ix)
    dummy=freud(head,olog=olog,ix=2*cube_ix+1)
  endif else dummy=freud(head,olog=olog,ix=cube_ix)


  ;;clean the data: (IRDIS ONLY)
  if keyword_set(skies) then begin
    for frame_ix=0,nframes-1 do data_in[*,*,frame_ix] -= supersky
  endif

  if keyword_set(flatpath) then begin
    print,'  Cleaning cube'
    restore,flatpath
    for frame_ix=0,nframes-1 do begin

      ;flat field
      tempim=data_in[*,*,frame_ix]
      tempim/=flat

      ;remove bad pixels and cosmic rays
      ; IRDIS has a big block of bad pixels, and if this isn't picked up by the bad pixel map,
      ; then performing the steps in this order seems to help
      tempim[bad_pixels]=median(tempim)
      dummy=sigma_filter_nirc2(tempim,box_width_irdis,n_sigma=n_sigma_irdis,/all,/iterate,bad_pixels=wbad)
      tempim=fix_bad_pixels(tempim,bad_pixels=[bad_pixels,wbad],neighborhood=5)
      
      ; print,'    Hack to rely on the input bad pixel map to speed up the cleaning'
      ; tempim=fix_bad_pixels(tempim,bad_pixels=bad_pixels,neighborhood=5)

      data_in[*,*,frame_ix]=tempim
    endfor
  endif

  ;;Put the instrument specific code here:
  camera=olog.optic_cfg[(pol+1)*cube_ix]
  if camera eq 'IRDIS' then begin
    ;; For IRDIS, we need to split the data into the two channels
    print,'  IRDIS data detected. Splitting frame in two'
    frame1=data_in[0:1023,*,*]
    frame2=data_in[1024:*,*,*]
    datasz[1]=datasz[1]/2
    if keyword_set(pol) then merge_irdis=0
    if keyword_set(merge_irdis) then begin
      print,'merging the IRDIS channels!'
      data_in=fltarr(datasz[1],datasz[2],nframes*2,1)
      data_in[*,*,0:nframes-1]=reform(frame1,[datasz[1],datasz[2],nframes])
      data_in[*,*,nframes:*]=reform(frame2,[datasz[1],datasz[2],nframes])
      nchannels=1
    endif else begin
      data_in=fltarr(datasz[1],datasz[2],nframes,2)
      data_in[*,*,*,0]=reform(frame1,[datasz[1],datasz[2],nframes])
      data_in[*,*,*,1]=reform(frame2,[datasz[1],datasz[2],nframes])
      nchannels=2
    endelse
  endif else if camera eq 'IFS' then begin
    print,'  IFS data detected.'
    print,'  Suppressing over/underflow errors in sigma_filter_nirc2'
    nwavs=datasz[4]
    nchannels=1 ;;for centering (we dont want the wavelength slices centred independently...)

    ;;we can put the IFS bad pixel code here
    for frame_ix=0,nframes-1 do begin
      for wav=0,nwavs-1 do begin
        tempim=data_in[*,*,frame_ix,wav]
        ;remove bad pixels and cosmic rays
        ; IFS bad pixels often come in blocks rather than individual pixels.
        !except=0
        dummy=sigma_filter_nirc2(tempim,box_width_ifs,n_sigma=n_sigma_ifs,/all,/iterate,bad_pixels=wbad)
        !except=1
        tempim2=fix_bad_pixels(tempim,bad_pixels=wbad,neighborhood=15)
        data_in[*,*,frame_ix,wav]=tempim
      endfor
    endfor
  endif else if camera eq 'ZIMPOL' then begin
    print,'  ZIMPOL data detected'
    stop;;need to work on this
  endif else begin
    print,'Unknown instrument name:'
    print,olog.instrument[(pol+1)*cube_ix]
    stop
  endelse

  if not keyword_set(no_centre) then begin
    ;;centre data:
    print,'  Centring data'
    xpeaks=[] & ypeaks=[]

    if keyword_set(centre_once_per_cube) then begin
     ;;loop through channels (IRDIS or hwp). IFS wavelength channels dont need to be centred independently
      for channel_ix=0,nchannels-1 do begin
        summed_cube=total(data_in[*,*,*,channel_ix],3);;sum over frames

        if camera eq 'IFS' then summed_cube=median(total(data_in,3),dim=3) ;; median over wavelength channels (don't need to centre them separately)

        bw = 30 ;; border of the detector to be avoided
        myKer = shift(exp(-(dist(11,11)/(2.0))^2), 5,5)
        temp0 = convol(summed_cube, myKer)
        dimx=datasz[1] & dimy=datasz[2]
        mx = max(temp0[bw:dimx-bw, bw:dimy-bw], mxy)
        ind = array_indices(temp0[bw:dimx-bw, bw:dimy-bw], mxy)
        xpeak = ind[0]+bw & ypeak = ind[1]+bw
        
        xpeaks=replicate(xpeak,nframes) ; Recording Peak Position.
        ypeaks=replicate(ypeak,nframes)

        ;;Now centre the frames
        if camera eq 'IFS' then begin
          data_in = shift(data_in,[dimx/2,dimy/2,0,0]-[xpeak,ypeak,0,0])
        endif else begin
          data_in[*,*,*,channel_ix]=shift(data_in[*,*,*,channel_ix],[dimx/2,dimy/2,0]-[xpeak,ypeak,0])
        endelse

      endfor
    endif else begin

      ;;loop through frames:
      for frame_ix=0,nframes-1 do begin
        ;;loop through channels (IRDIS or hwp). IFS wavelength channels dont need to be centred independently
          for channel_ix=0,nchannels-1 do begin
            im=data_in[*,*,frame_ix,channel_ix]

            if camera eq 'IFS' then im=median(data_in[*,*,frame_ix,*],dim=4) ;; median over wavelength channels (don't need to centre them separately)
            
            bw = 30 ;; border of the detector to be avoided
            myKer = shift(exp(-(dist(11,11)/(2.0))^2), 5,5)
            temp0 = convol(im, myKer)
            dimx=datasz[1] & dimy=datasz[2]
            mx = max(temp0[bw:dimx-bw, bw:dimy-bw], mxy)
            ind = array_indices(temp0[bw:dimx-bw, bw:dimy-bw], mxy)
            xpeak = ind[0]+bw & ypeak = ind[1]+bw
            
            xpeaks=[xpeaks,xpeak] ; Recording Peak Position.
            ypeaks=[ypeaks,ypeak]

            im=shift(im,[dimx/2,dimy/2]-[xpeak,ypeak])
            data_in[*,*,frame_ix,channel_ix]=im
          endfor
        ;;also record peak position (just mean the channels rather than display them all)
        fstats[0,frame_ix]=mean(xpeaks)
        fstats[1,frame_ix]=mean(ypeaks)
      endfor
    endelse
  endif else begin
    print, '  Data not centred!'
    fstats[0,frame_ix]=cutsize/2
    fstats[1,frame_ix]=cutsize/2
  endelse

  ;;cut the images
  cube=data_in[datasz[1]/2-cutsize/2:datasz[1]/2+cutsize/2-1,datasz[2]/2-cutsize/2:datasz[2]/2+cutsize/2-1,*,*]

  ;Fill in framestats 
  for frame_ix=0,nframes-1 do begin
    fstats[2,frame_ix]=total(cube[*,*,frame_ix,*])  ;; keep count of total flux,
    fstats[3,frame_ix]=max(cube[*,*,frame_ix,*])    ;; peak flux and
    fstats[4,frame_ix]=median(cube[*,*,frame_ix,*]) ;; background (actually median flux)
  endfor
  
  ;now flag bad data using code stolen from qbe_conica
  good_frames=flagbad_conica(fstats,cube,discard_sigma,tot_bad)
  good_index=where(good_frames ge 0)
  bad_frames=where(good_frames lt 0)
  ;; %%%%% Plot out frame statistics
  identifier = jobname+"_"+string(frames[cube_ix],format="(I4.4)")+"a.ps"
  fnumstr=[string(frames[cube_ix],format="(I4.4)")]

  plota_nirc2,fstats,jobname,olog,strcompress(string(cube_ix),/remove_all),    $
             '','',bad_frames=bad_frames,identifier=identifier

  ;;remove the bad frames:
  cube=cube[*,*,good_index,*]
  fstats=fstats[*,*,good_index]

  ;; Save the files out
  if not keyword_set(pol) then begin
    filestring=[string(cube_ix,format="(I4.4)")]
    filename='cube'+filestring+'.fits'
    olog.cube_fname[cube_ix,0]=filename
    olog.cube_fname[cube_ix,1]=filestring
    writefits,filename,cube

    olog.cube_tsize[cube_ix]=tsize[cube_ix]
    olog.cube_sz[cube_ix,*]= (size(cube))[1:3]

  endif else begin
    ;;these arrays are n_files instead of n_cubes long, so that the numbers from each file are preserved (mainly useful for non-cube mode instruments)
    olog.pas[cube_ix]=olog.pa[2*cube_ix]
    olog.hwps[cube_ix]=olog.hwp[2*cube_ix]

    filestring0=[string(cube_ix,format="(I4.4)")]+'_0'
    filestring1=[string(cube_ix,format="(I4.4)")]+'_1'
    filename0='cube'+filestring0+'.fits'
    filename1='cube'+filestring1+'.fits'
    cube0=cube[*,*,*,0]
    cube1=cube[*,*,*,1]
    olog.cube_fname[2*cube_ix,0]=filename0
    olog.cube_fname[2*cube_ix,1]=filestring0
    olog.cube_fname[2*cube_ix+1,0]=filename1
    olog.cube_fname[2*cube_ix+1,1]=filestring1
    writefits,filename0,cube0
    writefits,filename1,cube1

    olog.cube_tsize[2*cube_ix]=tsize[cube_ix]
    olog.cube_sz[2*cube_ix,*]= (size(cube))[1:3]
    olog.cube_tsize[2*cube_ix+1]=tsize[cube_ix]
    olog.cube_sz[2*cube_ix+1,*]= (size(cube))[1:3]

  endelse

endfor

; Generate a single postscript file from all "a.ps" files for ease of clockthru
; (optional this can be commented out)
spawn, 'psjoin *a.ps > aplots.ps'


;;now handle cal4src and tsize
cal4src =  intarr(n_cubes,  n_cubes) ; for each source. Set this var in flagging later.
targ  =  where(tsize lt 0,  complement = calib)
if(calib[0] ne -1 and targ[0] ne -1) then $
  for i =  0, n_elements(targ)-1 do cal4src[calib,targ[i]] =  1  
olog.tsize=tsize


; %%%%% Save the olog and framestats data structures for later use
save,olog,framestats,file='cubeinfo'+jobname+'.idlvar'

end