;; This function does the basic processing for raw zimpol images
;; It extracts the two sub-frames from the raw images, removes
;; the stripe from each image and then bins the x direction so 
;; the pixel scales are the same
;;
;; (The following lines assume no extra x or y binning beyond the default values of 2)
;; Input: image: a raw 1156 x 1024 x NDIT image
;;        head:  the header for the image extension
;; Output: return_im: a 537 x 512 x NDIT x 2 image
;;
function split_zimpol_im,image,head

;; We regularly run into trouble with floats and doubles unless we do this
image=double(image)
init_sz=size(image)

;; First use the prescan/overscan (middle) regions to bias subtract the images
n_prescan=sxpar_conica(head,'UT1PRSCX') ;; number of pixels
n_overscan=sxpar_conica(head,'UT1OVSCX') ;; number of pixels
overscan_centre=init_sz[1]/2

;; Loop through the exposures and bias subtract them
for ix=0,init_sz[3]-1 do begin
    ;; Left side
    image[0:overscan_centre-1,*,ix]-=mean(image[overscan_centre-lindgen(n_overscan/2),*,ix])
    ;; Right side
    image[overscan_centre:*,*,ix]-=mean(image[overscan_centre+lindgen(n_overscan/2),*,ix])
endfor

;; Now remove the overscan and prescan areas
;; First, the middle section
x_ix=indgen(init_sz[1]-n_overscan) ;; an array to mark the x indices that we want
x_ix[overscan_centre-n_overscan/2:*]+=n_overscan ;; shift all the values after the start of the overscan region by n_overscan
image=image[x_ix,*,*]
;; Then, the prescan (easy, assuming it starts at zero)
image=image[n_prescan:*,*,*]


sz=size(image)

;; Now split the rows to extract the two images
ix1=2*indgen(sz[2]/2)

im1=image[*,ix1,*]
im2=image[*,ix1+1,*]

if sz[0] eq 3 then begin
    return_im=fltarr(sz[1],sz[2]/2,sz[3],2)
    return_im[*,*,*,0]=im1
    return_im[*,*,*,1]=im2
endif else begin
    return_im=fltarr(sz[1],sz[2]/2,2)
    return_im[*,*,0]=im1
    return_im[*,*,1]=im2
endelse

;; and bin it so the directions have the same pixel scale
if sz[0] eq 3 then begin
    return_im=rebin(return_im,[sz[1]/2,sz[2]/2,sz[3],2])
endif else if sz[0] eq 4 then begin
    stop ;; this shouldnt happen. If it does, I need to rewrite the splitting code above to work for 4D...
    ; return_im=rebin(return_im,[sz[1]/2,sz[2],sz[3],sz[4]])
endif

return,return_im
end