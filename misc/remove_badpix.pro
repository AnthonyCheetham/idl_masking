;;This program will remove bad pixels by ensuring that the power
;;outside of the region of interest for masking is minimized. What it
;;actually does is projects the signal in the support defined by the
;;region with zero power in principle on to the subspace spanned by
;;the bad pixels themselves. The metric tensor on this subspace is
;;inverted, enabling the inner products with the basis vectors to be
;;converted into coefficients in front of the basisvectors.
;;
;;Inputs:
;; in:  The input image.
;; u,v: The (u,v) coordinates of the baselines in Fourier pixels, i.e.
;;      linear_dimension/lambda*radians_per_pixel*size_of_in
;; d  : The mask hole diameter in the same units.
;; badpix : 
;;      A vector of bad pixels (i.e. equal to y_coord*size + x_coord)
;;
;;NB With the Sydney masking code conventions, one has to pass -u and
;;v, not u and v straight from the mf file.

function remove_badpix, in, u,v, d, badpix
;;Copy the input image to a new array.
im=in
;;Do nothing if there are no bad pixels.
if (badpix[0] eq -1) then return, in
im[badpix]=0 ;Start off with something simple
;;Preliminaries: set up our mask.
xsz = (size(im))[1]
ysz = (size(im))[2]
if (xsz ne ysz) then begin
 print, "Error: input array must be square!"
 stop
endif
p = make_pupil(xsz,2.2*d)
mask = make_pupil(xsz,4*d)
for i=0,n_elements(u)-1 do begin
 mask += shift(p,v[i],u[i])
 mask += shift(p,-v[i],-u[i])
endfor
mask = shift(mask eq 0,xsz/2,xsz/2)
;;Now find the projection vectors and metric tensor.
nbad = n_elements(badpix)
proj = dcomplexarr(xsz,ysz,nbad)
for i=0,nbad-1 do begin
 a = dblarr(xsz,ysz)
 a[badpix[i]]=1
 proj[*,*,i] = fft(a)*mask
endfor
metric = dblarr(nbad,nbad)
for i=0,nbad-1 do for j=0,nbad-1 do $
 metric[i,j]=total(proj[*,*,i]*conj(proj[*,*,j]))/total(modsq(proj[*,*,i]))
metric = cov2cor(metric)
;;As we have to invert the metric, we can't have a determinant near zero.
if determ(metric) lt 0.01 then begin
 print, 'Error: misbehaved bad pixels!'
 stop
endif
coeffs = dblarr(nbad)
ftim = fft(im)
;;Find the coefficients by projecting onto the bad pixel basis vectors.
for i=0,nbad-1 do coeffs[i]=total(proj[*,*,i]*conj(ftim))/total(modsq(proj[*,*,i]))
;;
im[badpix]-=invert(metric)#coeffs

return, im
end
