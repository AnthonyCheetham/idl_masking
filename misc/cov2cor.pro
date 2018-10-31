;This function returns a correlation matrix when given a covariance
;matrix
;sig returns the standard deviations of the parameters (the sqrt of
;the diagonal)

function cov2cor, cov, sig=sig
sz=size(cov)
n = sz[1]
if (n ne (size(cov))(2)) then begin
   print, 'Input matrix must be square.'
   stop
endif

if sz[0] eq 3 then begin
   n3=sz[3]
   x = indgen(n)
   cor=0*cov
   sig=reform(cor[*,0,*])
   for k=0,n3-1 do begin
      ks=rebin([k],[n])
      sig[*,k] = sqrt(cov[x,x,ks])
      cor[*,*,k] = cov[*,*,k]
      for i = 0,n-1 do for j = 0,n-1 do cor[i,j,k] = cov[i,j,k]/sig[i,k]/sig[j,k]
   endfor
endif else begin
   x = indgen(n)
   sig = sqrt(cov[x,x])
   cor = cov
   for i = 0,n-1 do for j = 0,n-1 do cor[i,j] = cov[i,j]/sig[i]/sig[j]

endelse
return, cor

end
