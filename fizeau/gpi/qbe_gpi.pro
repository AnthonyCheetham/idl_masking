;;Program for cubing GPI data
;;Developed by ACC over Dec 2013 - 2014
;;Lots of code borrowed or inspired by the nirc2, conica routines
;;
;; show_images: Will display images and their power spectra
;; wav_channel: Process a single channel only (Not guaranteed to work!)
;; prefix, extn: prefix and file extension.
;; cutsize: the size of the output images (we usually want 256)
;; no_centre/no_center: turn of the centring algorithm
;; pol: Process polarised data (default is for spectrally dispersed data)
;; discard_sigma: Used to flag bad data. set=-1 to not flag any data!
;;    discard_sigma=[x_pos,y_pos,total_counts,peak_pixel,fraction of peak median]
;;    default is  [-2., -2., 3.5,  3.5, 0.4]
;;    Negative elements are not actually used.
;; nwav = number of wavelength channels (default = 37)

;; Things to do:
;;              - fix bad pixels
;;              - use actual flat fields?
;;              - combine wavelength channels before centering, then centre all by the same amount (currently does one channel only)
pro qbe_gpi,data_dir,jobname,flatpath,frames,skies,tsize,show_images=show_images,wav_channel=wav_channel,prefix=prefix,extn=extn,cutsize=cutsize,no_centre=no_centre,no_center=no_center,pol=pol,discard_sigma=discard_sigma,nwav=nwav

;;Here is where all of the fixed parameters go
box_width=5 ;; bad pixel identification (if a pixel is n_sigma away from its neighbours in box_width it is bad)
n_sigma=7   ;; bad pixel identification
badpix_neighbourhood=15 ;; bad pixel replacement (the region to draw the replacement from)
noise_amp =1e-4 ;; the relative amplitude of the noise to be added to the nan region (so the bad pixel check doesnt flag them)

if keyword_set(no_center) then no_centre=1

if keyword_set(prefix) eq 0 then prefix='S20120919S'
if keyword_set(extn) eq 0 then extn='.fits'
if keyword_set(show_images) then !p.multi=[0,2,1]

;;these numbers inherited from NACO, but they arent good defaults
if not keyword_set(discard_sigma) then discard_sigma=[-2., -2., 3.5,  3.5, 0.4]

if flatpath ne '' then restore,flatpath else flat =1.

nblocks=n_elements(uniq(tsize))
if keyword_set(pol) then n_cubes=8*nblocks else n_cubes=nblocks
wblock=uniq(tsize)
if keyword_set(nwav) eq 0 then nwav=37

if keyword_set(pol) then nwav=2
if keyword_set(pol) then pol=1 else pol=0 ;;we will use this in a clever way later
if keyword_set(wav_channel) then nwav=n_elements(wav_channel)

n_frames_total=n_elements(tsize)

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
      pas:fltarr(n_frames_total),hwps:fltarr(n_frames_total)}

;; Define data structure to hold all the frame statistics:
f=fltarr(n_cubes,2)     ;; frame statistic and err
fl= fltarr(n_cubes,10,2);; average and errors. src#, aperture size, 0=mean/1=error
;;the 10 on the previous line comes from qbe_conica, so dont blame me!
framestats={xpk:f,ypk:f,totflx:f,pkflx:f,skybgr:f,phot:fl}

if keyword_set(cutsize) then imsz=cutsize else imsz=281 ;; assume images are 281x281 pixels

start_ix_frame=0 ;; counts what frame number we're up to for pol data

for blocknum=0,nblocks-1 do begin
   block_files=where(tsize eq tsize[wblock[blocknum]])

   tsz=tsize[block_files[0]]
   nf=n_elements(block_files)

   ;;get the cube ready
   cube=fltarr(imsz,imsz,nf,nwav)

   ;;make an olog structure to save info for each file, so that
   ;; they can be combined into a single structure for each block
   ;; (or 8 for pol- vert and horiz for each hwp rot)
   cube_olog=make_olog(nf,tsize,frames,skies)
   
   ;;collect statistics about each frame (xpeak,ypeak,tot_flux,pk_flux,sky_bg)
   fstats=fltarr(5,nf)

   for ix=0,nf-1 do begin
      ;;read in the data and the two headers (info is spread over both)
      in_name=data_dir+prefix+strtrim(string(frames[block_files[ix]],format="(I4.4)"))+extn

      head=headfits(in_name)
      data_in=readfits(in_name,head2,exten=1,/silent)

      ;; Make sure the images are the right size
      if (size(data_in))[1] ne 281 then begin
         current_sz = (size(data_in))[1]
         ;; Use the median value of the top right corner as a placeholder for any missing pixels
         med_val = median(data_in[0:20,0:20,*])
         im_temp = fltarr(281,281,(size(data_in))[3])+med_val 

         im_temp[281/2-current_sz/2:281/2+current_sz/2-1,281/2-current_sz/2:281/2+current_sz/2-1,*] = data_in
         data_in = im_temp
      endif
      
      if ix eq 0 then first_header=[head,head2]
      if ix eq nf-1 then last_header=[head,head2]
      if ix eq 0 then hwp=fltarr(nf)

      print,'    Cleaning and centering image',ix+1,' of',nf

      ;;take the wavelength channels of interest and discard the rest
      if keyword_set(wav_channel) then begin         
         data_in=data_in[*,*,wav_channel]
         if (size(flat))[0] eq 3 then flat=flat[*,*,wav_channel]
      endif

      if (size(flat))[0] eq 2 then begin
         for wav=0,nwav-1 do data_in[*,*,wav]/=flat
      endif else if (size(flat))[0] eq 3 then begin
         data_in=data_in/flat
      endif
      
      ;;median subtract and set the nan coordinates to zero, then fix the bad pixels/cosmic rays
      nbad=0
      for wav=0,nwav-1 do begin
         temp=data_in[*,*,wav]
         goodpix=where(finite(data_in[*,*,wav]),complement=nans)

         fstats[2,ix]=total(temp[goodpix])  ;; keep count of total flux,
         fstats[3,ix]=max(temp[goodpix])    ;; peak flux and
         fstats[4,ix]=median(temp[goodpix]) ;; background (actually median flux)
         
         temp-=median(temp[goodpix])
         temp[nans]=0+noise_amp*(max(temp[goodpix])-min(temp[goodpix]))*randomn(seed,n_elements(nans)) ;;add a tiny amount of noise so they arent flagged as bad in the next step

         ;;remove bad pixels and cosmic rays
         dummy=sigma_filter_nirc2(temp,box_width,n_sigma=n_sigma,/all,/iterate,bad_pixels=wbad)
         ;; remove the nan pixels
         tempim=fix_bad_pixels(temp,bad_pixels=wbad,neighborhood=badpix_neighbourhood)
         nbad+=n_elements(wbad)

         data_in[*,*,wav]=tempim
         
      endfor
      print,'  ',nbad,' bad pixels flagged'
      ;;normalise peak flux by channel
      fstats[2,ix]/=nwav

      ;;centre data:
      if not keyword_set(no_centre) then begin
         xpeaks=[] & ypeaks=[]
         for wav=0,nwav-1 do begin
            im=data_in[*,*,wav]
            
            bw = 30 ;; border of the detector to be avoided
            myKer = shift(exp(-(dist(11,11)/(2.0))^2), 5,5)
            temp0 = convol(im, myKer)
            dimx=(size(im))[1] & dimy=(size(im))[2]
            mx = max(temp0[bw:dimx-bw, bw:dimy-bw], mxy)
            ind = array_indices(temp0[bw:dimx-bw, bw:dimy-bw], mxy)
            xpeak = ind[0]+bw & ypeak = ind[1]+bw
            
            xpeaks=[xpeaks,xpeak] ; Recording Peak Position.
            ypeaks=[ypeaks,ypeak]

            

            if keyword_set(pol) and wav eq 1 then begin
               ;; if it is polarized data, we want to shift both frames by
               ;; the same amount so that differential phase works properly.
               xpeak=xpeaks[0]
               ypeak=ypeaks[0]
            endif
            im=shift(im,[dimx/2,dimy/2]-[xpeak,ypeak])
            data_in[*,*,wav]=im
         endfor

         ;;also record peak position
         fstats[0,ix]=median(xpeaks)
         fstats[1,ix]=median(ypeaks)
      endif

      if keyword_set(show_images) then begin
         ;;turn the Nan's into zeroes (just for testing)
         goodpix=where(finite(data_in))
         w=array_indices(data_in,goodpix)
         sz=size(data_in)
         disp_im=fltarr(sz[1],sz[2])
         disp_im[w[0,*],w[1,*]]=data_in[w[0,*],w[1,*]]
         image_cont,disp_im[*,*],/n,/a,tit=in_name
         image_cont,abs(shiftfft(disp_im,281/2)),/n,/a
         wait,0.2
      endif

      ;;cut the image
      if keyword_set(cutsize) then data_in=data_in[281/2-imsz/2:281/2+imsz/2-1,281/2-imsz/2:281/2+imsz/2-1,*]
      
      ;;flip the data so east is CCW of north (as for every other instrument)
      data_in=reverse(data_in)

      ;;add to cube
      cube[*,*,ix,*]=data_in
      
      ;;get info for each file, to be put into a single olog
      ;; structure later (or 8 for pol)
      this=freud([head,head2],ix=ix,olog=cube_olog)
   endfor

   ;;now flag bad data using code stolen from qbe_conica
   good_frames=flagbad_gpi(fstats,cube,discard_sigma,tot_bad)
   good_index=where(good_frames ge 0)
   bad_frames=where(good_frames lt 0)
   ;; %%%%% Plot out frame statistics
   identifier = jobname+"_"+string(frames[block_files[0]],format="(I4.4)")+"a.ps"
   fnumstr=[string(frames[blocknum],format="(I4.4)")]

   plota_nirc2,fstats,jobname,cube_olog,strcompress(string(ix),/remove_all),    $
               '','',bad_frames=bad_frames,identifier=identifier
   if keyword_set(show_images) then !p.multi=[0,2,1]
   
   ;; Remove bad frame data from all arrays
   if ((size(bad_frames))[0] eq 0) then num_bad=0 else num_bad=n_elements(bad_frames)

   if(num_bad gt 0) then begin
      print,"Removing ",num_bad," flagged frames of bad data from this cube"
      cube=cube(*,*,good_index,*)
      fstats=fstats(*,*,good_index)
      print,'Warning: cube_olog was not updated with removed frames in qbe_gpi!!'
      ;stop
   endif    
   ;; %%%%% Reduce Aperture Photometry to averages 
   ;; %%%%% and Populate framestats data structure.
   framestats.xpk[blocknum,*]        =[mean(fstats[0,*]),stdev(fstats[0,*])] 
   framestats.ypk[blocknum,*]        =[mean(fstats[1,*]),stdev(fstats[1,*])] 
   framestats.totflx[blocknum,*]     =[mean(fstats[2,*]),stdev(fstats[2,*])] 
   framestats.pkflx[blocknum,*]      =[mean(fstats[3,*]),stdev(fstats[3,*])] 
   framestats.skybgr[blocknum,*]     =[mean(fstats[4,*]),stdev(fstats[4,*])] 
   
   ;;may be an extra dimension if we are only looking at one wavelength
   cube=reform(cube)

   ;;combine the cube_olog entries into one olog entry
   ;; first get the majority of info from the last file
   ;; this also cleverly handles the pol data as well by
   ;; duplicating the entries for all hwp angles by assuming pol=1 or pol=0
   start_ix=blocknum*(1+7*pol)
   end_ix=start_ix+7*pol
   for pol_ix=start_ix,end_ix do this=freud([head,head2],ix=pol_ix,olog=olog)
   olog.del_pa[start_ix:end_ix]=0.
   ;;save the averaged quantities
   olog.pa[start_ix:end_ix]=mean(cube_olog.pa)
   ;;save the frame-by-frame quantities. For non-pol data, these will be the same as the averaged quantities
   end_ix_frame=start_ix_frame+nf-1
   olog.pas[start_ix_frame:end_ix_frame]=cube_olog.pa
   olog.hwps[start_ix_frame:end_ix_frame]= (cube_olog.hwp mod 360)

   if olog.del_pa[start_ix] gt 180. then begin
      print,'WARNING: the instrument PA passes through zero (which wasnt allowed for in qbe_gpi)' ;; olog.pa, olog.del_pa above are the culprits
      stop
   endif
   olog.elevation[start_ix:end_ix]=mean(cube_olog.elevation)
   olog.del_elev[start_ix:end_ix]=max(cube_olog.elevation)-min(cube_olog.elevation)
   olog.utc[start_ix:end_ix]=cube_olog.utc[nf/2] ;;cant average strings, so take the middle
   olog.cube_tsize[start_ix:end_ix]=tsz
   olog.jd[start_ix:end_ix]=mean(cube_olog.jd)

   if keyword_set(pol) then begin
      ;;collect files with the same hwp rotation
      ;;I'm assuming we're only ever going to want 4: 0,22,45,68
      hwp = (cube_olog.hwp mod 360)
      sort_hwp=hwp[sort(hwp)]
      uniq_hwp=sort_hwp[uniq(sort_hwp)]
      if n_elements(uniq_hwp) ne 4 then begin
         print,'WARNING: number of half wave plate rotations is not 4!'
         stop
      endif
      
      for hwp_rot=0,n_elements(uniq_hwp)-1 do begin
         w0=where(hwp eq uniq_hwp[hwp_rot])
         save_cube1=cube[*,*,w0,0]
         save_cube2=cube[*,*,w0,1]

         olog.hwp[blocknum*8+hwp_rot*2:blocknum*8+hwp_rot*2+1]=uniq_hwp[hwp_rot]
         olog.cube_sz[blocknum*8+hwp_rot*2,*]=(size(save_cube1))[1:3]
         olog.cube_sz[blocknum*8+hwp_rot*2+1,*]=(size(save_cube2))[1:3]

         ;; Save the files out
         filestring=[string(blocknum*4+hwp_rot,format="(I4.4)")] ;;assume 4 hwp rotations
         filename='cube'+filestring+'_0.fits'
         olog.cube_fname[blocknum*8+2*hwp_rot,0]=filename
         olog.cube_fname[blocknum*8+2*hwp_rot,1]=filestring+'_0'
         writefits,filename,save_cube1

         filestring=[string(blocknum*4+hwp_rot,format="(I4.4)")]
         filename='cube'+filestring+'_1.fits'
         olog.cube_fname[blocknum*8+2*hwp_rot+1,0]=filename
         olog.cube_fname[blocknum*8+2*hwp_rot+1,1]=filestring+'_1'
         writefits,filename,save_cube2
      endfor

   endif else begin

      ;; Save the files out
      filestring=[string(blocknum,format="(I4.4)")]
      filename='cube'+filestring+'.fits'
      olog.cube_fname[blocknum,0]=filename
      olog.cube_fname[blocknum,1]=filestring
      writefits,filename,cube
   endelse
   start_ix_frame=end_ix_frame+1

   print,'Finished block',blocknum+1,' of',nblocks
endfor

;; Generate a single postscript file from all "a.ps" files for ease of clockthru
;; (optional this can be commented out)
spawn, 'psjoin *a.ps > aplots.ps'

;;and some file independent stuff:
olog.cube_sz=[281,281,nf]
olog.logflag=0
save,olog,framestats,file='cubeinfo'+jobname+'.idlvar'

if keyword_set(show_images) then !p.multi=0

end
