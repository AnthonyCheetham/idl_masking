;; ACC 2016

;; This program takes a list of files with dithered NIRC2 Lp data, and turns them into a flat field and bad pixel map
;; The idea is that the background level is much higher than the dark current, so all of the flux (except for the star)
;; is basically giving you a flat field anyway.

;; Step 1: Find the star in each frame and mask it out (while flagging bad pixels)
;; Step 2: Average the remaining frames together
;; Step 4: Divide by mean flux to get the flat field

; USAGE:
;;
;;INPUTS
;;  files: filenames for the frames to use
;; cutoff: number of times a pixel has to be flagged to be counted as bad
;;
;;OUTPUTS
;;   flat: the flat field
;;    bad: the bad pixel map
;;
;;OPTIONS
;;sigma : this is the number of standard deviations away from the mean
;;        that a pixel has to be before it is called bad. Obviously
;;        assumes that it is roughly Gaussian. Default=5
;; smooth_width: The width of the window used to smooth the image to find the star. Default=7 pixels
;; full_auto: Set this to 1 to automate the bad pixel rejection

pro dither_flat_nirc2,files,flat,bad_pixels,cutoff,sigma=sigma,smooth_width=smooth_width,$
    star_mask_size=star_mask_size,full_auto=full_auto

if not keyword_set(sigma) then sigma = 5
if not keyword_set(smooth_width) then smooth_width = 7
if not keyword_set(star_mask_size) then star_mask_size = 100


nflats=n_elements(files)

;; Initialise arrays to collect the flat and keep track of how many frames contribute to each pixel
flat=0
pixelcount=0
badpixcount=0

;; Loop through the files
for ix=0,nflats-1 do begin
    ;; Load file
    image=readfits(files[ix],head,/silent)

    if ix eq 0 then begin
        pixelcount=0*image
        badpixcount=0*image
    endif

    ;; Find the star
    dimx=(size(image))[1]
    dimy=(size(image))[2]
    temp = (smooth(image,nint(smooth_width),/edge_truncate))[smooth_width:dimx-smooth_width, smooth_width:dimy-smooth_width]
    star_centre=array_indices(temp,where(temp eq max(temp)))
    star_centre+=[smooth_width+1,smooth_width+1] ;; the actual centre position (since we truncated the image before smoothing)

    ;; mask out the star. Need to work the coordinates out properly
    minx=max([star_centre[0]-star_mask_size,0])
    maxx=min([star_centre[0]+star_mask_size,dimx-1])
    miny=max([star_centre[1]-star_mask_size,0])
    maxy=min([star_centre[1]+star_mask_size,dimy-1])

    ;; Mask out the star
    image[minx:maxx,miny:maxy]=0

    ;; use the plothist code from lampflat2 to work out the bad pixels, making sure to only plot the pixels without the star!
       ;;Anthony's quick way
   hist=histogram(image,locations=xaxis)
   plothist,image,xhist,yhist,/noplot
   med=median(image)
   sd=sqrt(variance(image))
   
   ;;truncate the histogram at a certain number of standard deviations
   ;;from the median
   max_ok=med+sigma*sd
   min_ok=med-sigma*sd
   
  
   ;;Now display a plot so we can check that it worked
   plot,xaxis,hist,title='Click outside area to accept bounds, or inside to switch to manual bounds'
   oplot,[min_ok,min_ok],[0,5*max(hist)],line=1
   oplot,[max_ok,max_ok],[0,5*max(hist)],line=2
   if keyword_set(full_auto) then v1 =0 else cursor,v1,a,/device
   wait,0.2
   
   ;;the edge of the plot should be 60 in device coordinates, so use
   ;;that
   if v1 lt 60 then good=where(image lt max_ok and image gt min_ok ,complement=w) else begin
       try_again:
       plothist,image,tit='Click on left then right Edge of allowed range'
       cursor,v1,a &wait,.6
       print,' Lower Bound ',v1
       cursor,v2,a & wait,.6
       print,' Upper Bound ',v2
       good=where(image gt v1 and image lt v2)
       if good[0] eq -1 then begin
           print,'You selected zero pixels! Try again!'
           goto,try_again
       endif
       plothist,image[good],tit='SECOND CHANCE: Click on left then right Edge of allowed range'
       cursor,v1,a &wait,.6
       print,' Lower Bound 2',v1
       cursor,v2,a & wait,.6
       print,' Upper Bound 2',v2
       good=where(image gt v1 and image lt v2)
       if good[0] eq -1 then begin
           print,'You selected zero pixels! Try again!'
           goto,try_again
       endif
       w=where(image lt v1 or image gt v2)
       plothist,image[good],tit='Final Histogram'
   endelse
   flat=flat+image
   badpixcount(w)=badpixcount(w)+1

    ;; Add it to the big array
    pixelcount+=1 ;; the easy way is to add 1 to each pixel then subtract 1 from the masked pixels
    pixelcount[minx:maxx,miny:maxy]-=1

endfor

bad_pixels=where(badpixcount gt cutoff,complement=goodpix)
flat=flat/pixelcount

;; Now normalise the flat
flat/=median(flat[goodpix])
flat[bad_pixels]=1.

;; Flag extra bad pixels as the ones that are significantly different from the regular flat
flat_sigma=median(abs(flat-1.)) ;; use the median absolute deviation, which is robust to outliers
flat_hotpix=where(flat gt (1+7*flat_sigma))
flat_coldpix=where(flat lt (1-7*flat_sigma))

;; combine the bad pixels
badpixmap=0*flat
badpixmap[bad_pixels]=1
badpixmap[flat_hotpix]=1
badpixmap[flat_coldpix]=1

bad_pixels=where(badpixmap gt 0.5)
print,n_elements(bad_pixels),' bad pix flagged'

flat[bad_pixels]=1.

end