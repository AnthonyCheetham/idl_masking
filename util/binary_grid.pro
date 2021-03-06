;Tau6/LkCa4_0247                 42.52   78.73   31.61

;;This procedure does a binary grid search, outputting the detection
;;limits for a binary, as well as finding the best fitting binary.
;;
;;OPTIONS:
;;  printfile:   print results to a file instead of the screen.
;;  usevis:      0 Don't use amplitudes.
;;               1 Use amplitudes 
;;               2 Use amplitudes after fitting for seeing and adding
;;                errors.
;;   annuli:     Set this to an array containing the edges of each annulus
;;               that you want to use for the monte-carlo
;;               simulations. Default=[10,20,40,80,160,240,320]
;;   nfr:      
;;   nosim:      Set to skip the monte-carlo simulations
;;   gpi:        assume there are 37 spectral channels, 17 of which
;;                are independent (i.e. scale errors by sqrt(37/17)
;;   quiet_mode: Don't plot anything (useful for batch jobs), but save
;;                output so you can look at it later. 
;;                Set quiet_mode=save_file_name 
;;                (otherwise it will use binary_grid_output.idlvar)
;;   
;NB: -Common blocks can be reset with .reset_session
; - Need a routine that finds the Hessian matrix inverse...
;TODO:
; - This script doesn't work for moderate contrast-ratio
;   binaries. Needs fixing...?
; - Use visibility to help constrain upper limits for near 1:1
;   binaries at \lambda/2D. We can consider visibility amplitude to
;   give an independent constraint on the contrast ratio. This is the
;   same game: simulate V^2 data and fit to them. The catch is that
;   since we're talking about the non-linear regime, we'll need
;   non-linear fitting...

;;Options:
;; use_cperr_bin: For the Monte Carlo sims, use the closure phase
;;     errors for the binary model instead of the null hypothesis. If
;;     you have a strong detection, the cperr_null will be huuuge, and
;;     mean that your contrast limits will be wrong. If you have a
;;     detection, you should set this flag! This won't work with cp_cov
;;     or proj

;;!!!NB !!! For the Upper Sco paper, the limits at wide separations
;;used the following hacks:
;  rat995_160_240 = m[s[rnum-1]]/0.9  ;*0.9  due to the super-gaussian window
;  rat995_240_320 = m[s[rnum-1]]/0.67 ;*0.67 due to the super-gaussian
;  window
pro binary_grid,  infile,  printfile = printfile, no_bfinfo=no_bfinfo, usevis = usevis_in,  nfr = nfr,  nosim = nosim,  $
 init_crat = init_crat,  fix_crat = fix_crat,  cp_cov =  cp_cov_in, sep_range=sep_range,pa_range=pa_range,crat_range=crat_range,  $
 apriori = apriori_in,  proj = proj_in,  cor = cor,  covar = covar,  bestp = bestp,  perror = perror, scale_err=scale_err,$
 confidence_int = confidence_int, pa_keck9h=pa_keck9h, highc_only=highc_only, force_min_sep=force_min_sep0, $
 nsim=nsim, tom_this_doenst_work_for_nirc2=tom_this_doenst_work_for_nirc2, t3errscale=t3errscale,annuli=annuli,maxsep=maxsep,$
 use_cperr_bin=use_cperr_bin,nosave=nosave,minsep=minsep, extra_error=extra_error,gpi=gpi,quiet_mode=quiet_mode,t3data_in=t3data_in,$
final_perror=final_perror
common t3block,  t3data,  apriori, vis2data, usevis,  cp_cinv,  proj, proj_err, force_min_sep

if (keyword_set(usevis_in) eq 0) then usevis_in = 0
usevis = usevis_in
if (keyword_set(apriori_in) eq 0) then apriori_in = {val:[-1, -1, -1],  err:[0, 0, 0]}
apriori = apriori_in
if (keyword_set(printfile)) then openw,  1,  printfile, width=150 else printfile = ''
;;Convert errors in radians^(-2) to degrees^(-2)
if not keyword_set(cp_cov_in) then cp_cov = [-1.] else cp_cov = cp_cov_in/(!pi/180.)^2
cp_cov_null = cp_cov
cp_cinv = invert(cp_cov)
if not keyword_set(confidence_level) then confidence_level = 0.999
p = printfile
if (p ne '') then print,  '********** Binary fit for file: ',  p,  '**********' $
else print,  '********** Starting Binary fit. No save file **********'
if not keyword_set(nsim) then nsim =long(10000) 
if (keyword_set(init_crat) eq 0) then init_crat =  250
if keyword_set(fix_crat) then init_crat = fix_crat
extract_t3data,  file = infile,  t3data

if keyword_set(gpi) then begin
  print,'Assuming GPI data with 37 wavelength channels (17 independent)'
  t3data.t3phierr*=sqrt(37./17.)
endif

if (keyword_set(t3errscale)) then t3data.t3phierr *= t3errscale
w = where(t3data.flag)
if w[0] ne -1 then t3data[w].t3phierr = 1e4

;;If we have the matrix, project the closure-phases onto the
;;linearly dependent set of closure-phases.
if (keyword_set(proj_in)) then begin
 proj = proj_in
 nproj = (size(proj))[2]
 if (cp_cov[0] eq -1) then stop
 cp_cov = transpose(proj)#cp_cov#proj
 ix =  indgen(nproj)
 proj_err = sqrt(cp_cov[ix, ix])
 cp_cov = [-1]
 cp_cinv = [-1]
endif else begin
  proj = [-1]
  proj_err = 0
endelse
extract_vis2data,  file = infile,  vis2data
read_oidata,infile,oiarray,oitarget
target_name=oitarget.target
target_name=target_name[0]

;;ACC HACK to input a modified t3data array (thats not in an oifits file):
if keyword_set(t3data_in) then begin
   t3data=t3data_in
   print,'Anthony left a hack in binary_grid...'
endif

;Enforce a maximum SNR of 5.0 on visibility data.
vis2data.vis2err = sqrt(vis2data.vis2err^2 + (0.2*vis2data.vis2data)^2)
;;Find the number of rotations on the sky (i.e. independent cubes) by
;;finding how many times the first baseline in the file is repeated,
;;where we define "baseline" as containing uniqu stations indices.
;;WARNING: breaks for a >100 hole mask.
ix = vis2data.sta_index[0] + 100*vis2data.sta_index[1]
nrotations = n_elements(where(ix eq ix[0]))
nv2 = n_elements(vis2data)/nrotations
nt3 = n_elements(t3data)/nrotations
;Number of independent degrees of freedom in the data. The complicated
;expression at the end is (n_holes - 1)
ndf =  nrotations*(nv2 - (3*nt3/nv2+1) )

  if (usevis le 0) then begin
    print, nrotations, (nv2 - (3*nt3/nv2+1) ), ndf, $
      format = '("Degrees of Freedom:",I5," cubes x ",I5," independent clp data = ",I6)'
    if (p ne '') then printf,1,nrotations, (nv2 - (3*nt3/nv2+1) ), ndf, $
      format = '("Degrees of Freedom:",I5," cubes x ",I5," independent clp data = ",I6)'
endif

if (keyword_set(nfr)) then begin
 ;Median single-frame noise in radians
 med_noise =  median(t3data.t3phierr)*sqrt(nfr)/180*!pi
 ;;Weight the number of dependent and independent closure-phases 
 ;;so that the number of degrees of freedom is intermediate between
 ;;these values for and SNR of 1.0. 
 ndf =  (n_elements(t3data)*med_noise + ndf*(1./med_noise))$
       /(                   med_noise +      1./med_noise)
endif

r =  sqrt(vis2data.u^2 + vis2data.v^2)
maxr =  max(r)
if (keyword_set(minsep) eq 0) then minsep =  rad2mas(1./4/maxr)
if (n_elements(force_min_sep0) gt 0) then force_min_sep=force_min_sep0 $
else force_min_sep=minsep
;;The next lines are a little hacked-up, but works for the 9h mask,
;;where the baseline separation really defines the field-of-view...
s =  sort(r) 
if (keyword_set(maxsep) eq 0) then maxsep =  rad2mas(1./2./min([r, sqrt(!pi*r[s[5]]^2/12.)]))
print,                     minsep, maxsep, $
                           format = '("Min and max separations in search (mas): ", 2F7.1)'
if (p ne '') then printf,1,minsep, maxsep, $
                         format = '("Min and max separations in search (mas): ", 2F7.1)'
;Set the minimum angle to be 1 radian of fringe phase at the max
;baseline...
delsep =  rad2mas(1/2./!pi/maxr)      ;1 radian. This is kind-of arbitrarily chosen...
nr   =  round((maxsep-minsep)/delsep) ;+2
delang = (delsep/maxsep)*180./!pi
nang =  round( 360./delang ) 
delang = 360./nang           
errtest=t3data.t3phierr

params =  fltarr(5, nr, nang)
params[3, *, *] = 0.1
params[4, *, *] = 0.1
for i = 0, nr-1 do for j = 0, nang-1 do begin
; if (i le 3) then  params[0, i, j] = minsep+i*delsep/2. $;separation: finer sampling
; else params[0, i, j] = minsep+(i-2)*delsep ;separation
   params[0, i, j] = minsep+i*delsep ;separation
   params[1, i, j] = j*delang        ;position angle
   params[2, i, j] = init_crat       ;brightness ratio
endfor

;;ACC: save the model cps now to save time later
nphi =  n_elements(t3data)
modcp =  fltarr(nphi, nr*nang)

;Now search for the best chi^2 over the grid...
chi2_arr =  fltarr(nr, nang)
chi2_vis = fltarr(nr, nang)
t0=systime(1)
for i =  0, nr-1 do for j = 0, nang-1 do begin
   modelt3 = binary_t3data(params[*, i, j],t3data=t3data)
   modcp[*, i+long(nr)*j] = modelt3.t3phi
   
   cpresid = mod360(modelt3.t3phi - t3data.t3phi)
   if (proj[0] ne -1) then begin
      chi2_arr[i, j] = total((cpresid#proj)^2/proj_err^2)
   endif else if (cp_cinv[0] eq -1) then chi2_arr[i, j] =  total(cpresid^2/t3data.t3phierr^2) else $
      chi2_arr[i, j] = transpose(cpresid)#cp_cinv#cpresid

   if j eq 0 then print,'done',i,'of',nr
   if (usevis gt 0) then begin
      modelv2 = binary_vis2atm(params[*, i, j],vis2data=vis2data)
      chi2_vis[i, j] = total((modelv2.vis2data-vis2data.vis2data)^2/vis2data.vis2err^2)
; Also need to add to the number of degrees of freedom when using vis+clps
      ndf =  nrotations*(nv2 + (nv2 - (3*nt3/nv2+1)))
   endif
endfor

if (usevis gt 0) then begin
   print, nrotations,nv2, (nv2 - (3*nt3/nv2+1)), ndf, $
          format = '("Degrees of Freedom:",I5," cubes x (",I5," vis + ",I5," clp data) = ",I6)'
   if (p ne '') then printf,1,nrotations,nv2, (nv2 - (3*nt3/nv2+1)), ndf, $
                            format = '("Degrees of Freedom:",I5," cubes x (",I5," vis + ",I5," clp data) = ",I6)'
endif
print,'Time taken for chi2 loop:',systime(1)-t0

;;Include visibility in finding the best initial fit.
;;Filter the chi^2 so that we don't actually look within lambda/2D, 
;;as a good fit at lambda/3D will be found from a wider grid search
;;anyway.
w =  where(params[0, *, 0] lt rad2mas(1./2/maxr))
filter_array =  fltarr(nr, nang)
filter_array[*] = 1.0
if (w[0] ne -1) then for i = 0, nang-1 do filter_array[w, i] =  10.0
dummy =  min((chi2_arr+chi2_vis)*filter_array, m)
in = array_indices(chi2_arr, m)

if (apriori.val[0] eq -1) then $
 apriori = {val:double([params[0, in[0], in[1]], params[1, in[0], in[1]], params[2, in[0], in[1]]]),err:double([1e3, 1e3, 1e5])}
if keyword_set(fix_crat) then begin
 apriori.err[2] = 0.01
 apriori.val[2] = init_crat
endif
xi = 0.1*identity(3)
bestp = apriori.val

powell,  bestp,  xi,  1e-4,  bestchi2,  'binary_oifits_chi2'
;;OK - no idea why this is so, but a negative contrast ratio is
;;     equivalent to a positive one.
if (bestp[2] lt 0) then bestp[2] *= -1

;Now we check (crudely?) that it was correct to try the "low contrast
;ratio" branch.
if not keyword_set(highc_only) then begin
 highp =  bestp
 highp[2] = 1.2
 xi = 0.1*identity(3)
 powell,  highp,  xi,  1e-4,  highchi2,  'binary_oifits_chi2'
 if (highp[2] lt 0) then begin
  highp[2] *= -1
  highp[1] = mod360(highp[1] + 360)
 endif
 if (usevis gt 0) then begin
  modelv2 =  binary_vis2atm([bestp, 0.1, 0.1],vis2data = vis2data)
  bestchi2 += total((vis2data.vis2data - modelv2.vis2data)^2/vis2data.vis2err^2)
  modelv2 =  binary_vis2atm([highp, 0.1, 0.1],vis2data = vis2data)
  highchi2 += total((vis2data.vis2data - modelv2.vis2data)^2/vis2data.vis2err^2)
 endif
 if (highchi2 lt bestchi2) then begin
  print,  'Using low contrast ratio solution!'
  lowp = bestp
  bestp = highp
 endif
endif
if not keyword_set(quiet_mode) then begin
   window,0
   image_cont, exp(-(chi2_arr-min(chi2_arr))/2./sqrt(n_elements(t3data)/ndf)/4.), /noc,  tit = 'Lhood at Contrast: '+strtrim(init_crat, 2),$
               xv = params[0, *, 0],  yv = params[1, 0, *],  xtit = 'Sep (mas)',  ytit = 'PA (degrees)'
endif
print,                                           "                       Sep(mas) PA(degs) Contrast"
if (p ne '') then printf,  1, "                       Sep(mas) PA(degs) Contrast"
print,  bestp,  format = '("Initial best solution: ", 3F9.2)'
if (p ne '') then printf,  1, bestp,  format = '("Initial best solution: ", 3F9.2)'

modelt3 = binary_t3data([bestp, 0.1, 0.1],t3data=t3data)
good = where(t3data.t3phierr lt 60)
;Find the closure-phase error to add in quadrature so that chi^2=1
extra_error =  [0, 0.01 * 10.^(alog10(300.)*findgen(29)/20.) ] ;;lots of errors to try
newchi2 =  fltarr(30)
newchi2_null = fltarr(30)
cpresid = double(mod360(t3data.t3phi-modelt3.t3phi))

if (cp_cinv[0] eq -1) then begin
   if (proj[0] eq -1) then begin
      ;;what would the chi2 (for both model and null) be if we used each value of extra_error
      for i =  0, 29 do newchi2[i] =  total(cpresid^2/(t3data.t3phierr^2+extra_error[i]^2)) ;
      for i =  0, 29 do newchi2_null[i] =  total(t3data.t3phi^2/(t3data.t3phierr^2+extra_error[i]^2))
      ;;turn chi squared into a number that would give reduced chi
      ;;squared of one.
      newchi2 = newchi2/n_elements(good)*ndf/(ndf-3.) ;The -3 is for the three parameters that are fit.
      newchi2_null = newchi2_null/n_elements(good)
   endif else begin
      for i =  0, 29 do newchi2[i] =  total((cpresid#proj)^2/(proj_err^2+extra_error[i]^2)) 
      for i =  0, 29 do newchi2_null[i] =  total((t3data.t3phi#proj)^2/(proj_err^2+extra_error[i]^2)) 
      newchi2 = newchi2/(nproj-3)
      newchi2_null = newchi2_null/nproj
   endelse
   extra_error_null = interpol(extra_error,  newchi2_null,  1) 
   extra_error =  interpol(extra_error,newchi2,1)
endif

if (cp_cinv[0] eq -1) then begin
   if (newchi2[0] lt 1.0) then begin
      extra_error =  0.0 
      print,  'New Chi2: ',  newchi2[0]
      if (p ne '') then printf,  1, 'New Chi2: ',  newchi2[0]
   endif
   if (newchi2_null[0] lt 1.0) then begin
      extra_error_null = 0.
      print,  'New Chi2: ',  newchi2_null[0]
      if (p ne '') then printf,  1, 'New Chi2 (null): ',  newchi2_null[0]
   endif
   print,  extra_error,  extra_error_null, $
           format = '("Extra error for binary and single star soln (degs): ", 2F7.2)' 
   if (p ne '') then printf, 1, extra_error,  extra_error_null, $
                             format = '("Extra error for binary and single star soln (degs): ", 2F7.2)'
   ;;Add the extra error in
                                ;extra_error_null = 0.0;!!! For testing.
   cperr_null = sqrt(t3data.t3phierr^2+extra_error_null^2)*sqrt(n_elements(good)/double(ndf))
   cperr_bin = sqrt(t3data.t3phierr^2+extra_error^2)*sqrt(n_elements(good)/double((ndf-3.))) ;;-3 due to extra params
   ;;Now we need to scale the error by sqrt(n_clps/ndf) (i.e. n_holes/3) to take into account the fact that not all clps are independent.
   if (keyword_set(scale_err)) then begin
      t3data.t3phierr = t3data.t3phierr*sqrt(newchi2[0]*n_elements(good)/ndf) 
   endif else t3data.t3phierr = sqrt(t3data.t3phierr^2 + extra_error^2)*sqrt(n_elements(good)/double(ndf))
   proj_err_null = sqrt(proj_err^2+extra_error_null^2)
   proj_err = sqrt(proj_err^2+extra_error^2)
endif else begin
   ;;OK - really what I want to do is add in extra errors to the
   ;;     covariance matrix.
; newchi2 = transpose(cpresid)#cp_cinv#cpresid/(n_elements(good)-3)
; newchi2_null = transpose(t3data.t3phi)#cp_cinv#t3data.t3phi/(n_elements(good)-3)
   print,  'Inverting and calculating extra error...'
   cor = cov2cor(cp_cov,  sig = sig)
   for i =  0, 29 do begin
      temp_cinv =  invert(double(extra_error[i]^2*cor + cp_cov))
      newchi2[i] =  transpose(cpresid)#temp_cinv#cpresid/(n_elements(good)-3)
      newchi2_null[i] = transpose(t3data.t3phi)#temp_cinv#t3data.t3phi/(n_elements(good)-3)
   endfor
   if (newchi2_null[0] lt 1.0) then begin
      print,  'Chi2 for 1-star solution: ',  newchi2_null[0]
      if (p ne '') then printf,  1, 'Chi2 for 1-star solution: ',  newchi2_null[0]
      cp_cov_null = cp_cov
   endif else begin
      extra_error_null = interpol(extra_error,newchi2_null,1) 
      if (newchi2[0] gt 1) then extra_error =  interpol(extra_error,newchi2,1) else extra_error = 0.
      scale = sqrt(newchi2[0])
      print,  extra_error,  extra_error_null, $
              format = '("Extra error for binary and single star soln (degs): ", 2F7.2)' 
      if (p ne '') then printf, 1, extra_error,  extra_error_null, $
                                format = '("Extra error for binary and single star soln (degs): ", 2F7.2)'
      cp_cov += extra_error^2*cor
      cp_cov_null += extra_error_null^2*cor
   endelse
endelse

cp_cinv = invert(double(cp_cov))
cp_cinv_null = invert(double(cp_cov_null))
;!!! Temp line: what is the delta chi^2 with cp_cinv_null?
;cp_cinv = cp_cinv_null
;!!!
;;now that we have a better estimate of the errors, fit again using powell
powell,  bestp,  xi,  1e-2,  fmin,  'binary_oifits_chi2'
chi2_bin =  binary_oifits_chi2(bestp, sig = perror,  corr = cor,  covar = covar)
print,  bestp,                        format = '("Final best solution:   ", 3F9.2)'
if (p ne '') then printf, 1, bestp,   format = '("Final best solution:   ", 3F9.2)'
print,  perror,                       format = '("Errors:                ", 3F9.2)'
if (p ne '') then printf, 1, perror,  format = '("Errors:                ", 3F9.2)'
print,  printfile, 2.5*alog10(bestp[2]), 1.085*perror[2]/bestp[2], bestp[0], perror[0], bestp[1], perror[1],$
  format = '(A15, 3(F7.3, "$\pm$", F5.3, " & "), "\\ %binary_grid fit")'
if (p ne '') then printf,  1,  printfile, 2.5*alog10(bestp[2]), 1.085*perror[2]/bestp[2], bestp[0], perror[0], bestp[1], perror[1],$
  format = '(A15, 3(F7.3, "$\pm$", F5.3, " & "), "\\ %binary_grid fit")'

;;now use the same errors from the binary case to work out the chi
;;squared for the null hypothesis, so we can compare them
temp = bestp
final_perror=perror
temp[2] = 1e5
temp2 = apriori.err[2]
apriori.err[2] = 1e5
chi2_null =  binary_oifits_chi2(temp, sig = perror)
apriori.err[2] = temp2
print,  chi2_bin,                         format = '("Chi-squared (binary): ",F9.2)'
if (p ne '') then printf,  1,  chi2_bin,  format = '("Chi-squared (binary): ",F9.2)'
print,                  chi2_null,        format = '("Chi-squared (null):   ",F9.2)'
if (p ne '') then printf,  1,  chi2_null, format = '("Chi-squared (null):   ",F9.2)' 
print, sqrt(chi2_null-chi2_bin),          format = '("Significance:         ",F9.2,"  = sqrt(chi2_null - chi2_bin)")'
if (p ne '') then printf,  1,  sqrt(chi2_null-chi2_bin), format = '("Significance:        ",F9.2,"  = sqrt(chi2_null - chi2_bin)")'

;;Save a new oifits file of the residuals only.
modelt3 = binary_t3data([bestp, 0.1, 0.1],t3data=t3data)
if not keyword_set(nosave) then begin
   read_oidata, infile, oiarray, oitarget,oiwavelength,oivis,oivis2, oit3
   for i=0,n_elements(oit3)-1 do *(oit3[i].t3phi) = t3data[i].t3phi-modelt3[i].t3phi
   write_oidata, infile+'.resid',oiarray,oitarget,oiwavelength,oivis,oivis2,oit3
endif

;;Now, if this solution is at a very small separation, see how many
;;sigma it really is...
if (bestp[0] lt rad2mas(2/3./maxr)) then begin
   old_apriori = apriori
   apriori.val[*] = bestp
   apriori.err[*] = [0.1,  1e3, 1e4]
   bestp2 = bestp
   dummy2 =  binary_oifits_chi2(bestp2, sig = perror2,  corr = cor2)
   print,  'Flux error at fixed separation: ',  string(perror2[2],  format = '(F6.2)'),  $
           ' and S/N: ', string(bestp2[2]/perror2[2],  format = '(F6.2)')
   apriori.val[0] = minsep
   xi = 0.1*identity(3)
                                ;Now make sure we start at a good place... 
   bestp2[0] =  apriori.val[0]
   powell,  bestp2,  xi,  1e-2,  fmin,  'binary_oifits_chi2'
   dummy2 =  binary_oifits_chi2(bestp2, sig = perror2,  corr = cor2)
   print, bestp2, format =                         '("Solution at minimum separation: ",3F7.2)'
   if (p ne '') then printf, 1, bestp2,   format = '("Solution at minimum separation: ",3F7.2)'
   print,  perror2, format =                       '("Errors:                         ",3F7.2)'
   if (p ne '') then printf, 1, perror2,  format = '("Errors:                         ",3F7.2)' 
   print,  fmin, format =                          '("Chi2:                           ",3F7.2)'
   if (p ne '') then printf, 1, fmin,  format =    '("Chi2:                           ",3F7.2)' 
endif

;;ACC edit: cperr_sim will be the errors used for the Monte-Carlo sims.
;; If you have a binary, it makes sense to use cperr_bin, but if you
;; have a point source, you should use cperr_null. Of course, you find
;; out whether you have a binary after you choose this, making this whole
;; choice dodgy. PENDING BETTER SOLUTION! 
if keyword_set(use_cperr_bin) then cperr_sim=cperr_bin else begin
  if cp_cov[0] eq -1 then cperr_sim=cperr_null
endelse
; print,'ACC hack in binary grid!'
; cperr_sim=errtest*sqrt(n_elements(good)/double(ndf))

;Now we have some good errors, so we can do Monte-Carlo...
;We will take advantage of the fact that for a given sep and PA, the
;fitted contrast ratio is total(modelt3.t3phi*t3data.t3phi/t3phierr^2)/total(modelt3.t3phi^2/t3phierr^2)
if (keyword_set(nosim) eq 0) then begin
   t3sim = t3data
   maxrats = fltarr(2*nr,  nsim)
   max_gtrats =  fltarr(2*nr-1, nsim) ;Maximum ratios for any separations larger than the r values
   
   ;;ACC: this used to calculate the model closure phases, but since we already do it in the first chi squared loop I have changed it to save it out there to save time.
   ;nphi =  n_elements(t3data)
   ;modcp =  fltarr(nphi, nr*nang)
   ;print,  'Calculating model closure-phases'
   ;;calculate the model closure phases for each combination of
   ;;sep and pa
   ;;for i = 0, nr -1 do for j = 0, nang-1 do modcp[*, i+long(nr)*j] = (binary_t3data(params[*, i, j],t3data=t3sim)).t3phi

   norm =  fltarr(nr*long(nang))

   ;;ACC: calculate null hypothesis chi squared for each model (i.e. (cp/error)^2)
   ;;*** Are you MJI? If not you only care about the diag_matrix on the following line (1)*** 
   if (cp_cinv_null[0] eq -1) then begin
      ;;ACC: I changed the next line to something that becomes many orders of magnitude faster when you have lots of clps (e.g. GPI).
      for i = long(0), long(nr)*nang-1 do norm[i] = total((modcp[*, i]/cperr_sim)^2)
   endif else begin
      dm = cp_cinv_null
      for i = long(0), long(nr)*nang-1 do norm[i] = modcp[*, i]#dm#modcp[*, i]
   endelse
   
   if (proj[0] ne -1) then begin
      for i = long(0), long(nr)*nang-1 do norm[i] = (modcp[*,i]#proj)^2#(1./proj_err^2)
   endif

   norm =  reform(norm, nr, nang,  /over)
   ;;Finally, if we have an input covariance matrix, we use this to
   ;;simulate correlated variables...
   if (cp_cinv_null[0] ne -1) then begin
      diag = la_eigenql(cp_cinv_null,  eigenvectors = unitary)
      diag = diag_matrix(1./sqrt(diag))
                                ;unitary = transpose(unitary)
   endif
   
   ;;now start the monte-carlo simulations    
   tp = transpose(proj)
   tbegin =  systime(1)
   for k = 0, nsim-1 do begin
      if (proj[0] ne -1) then begin
         proj_sim = proj_err*randomn(seed, nproj)
         dm = diag_matrix(1./proj_err^2)
         crats = 1./reform(params[2, *, *])*reform(proj_sim#dm#tp#modcp, nr, nang)/norm
      endif else begin 
         ;;generate a random scaling for the null closure phase error
         ;;*** Are you MJI? If not you only care about the following line (3)*** 
         if (cp_cinv_null[0] eq -1) then t3sim.t3phi =  cperr_sim*randomn(seed, nphi) $
         else t3sim.t3phi = unitary#diag#randomn(seed, nphi)
         ;;*** Are you MJI? If not you only care about the following line(4) *** 
         
         ;;ACC edit: Improved speed for non-covariance matrix case
         if  (cp_cinv_null[0] eq -1) then crats = 1./reform(params[2, *, *])*reform((t3sim.t3phi/cperr_sim^2)#modcp, nr, nang)/norm else crats = 1./reform(params[2, *, *])*reform(t3sim.t3phi#dm#modcp, nr, nang)/norm
      endelse
      
      ;;The following formula is a crude correction for using this
      ;;high-contrast limit formula for low contrast and separation.
      ;;... we have to invert:
      ;; Fs_old = Fs*(1-3*Fs+2*Fs^2), where Fs_old = crat_old/(1+crat_old)
      ;; and Fs = crat/(1+crat)
      crats = crats/((1.0 - 3.0*abs(crats))>0.4) ;;Prevent divergence...
      
;;     new_crats =  regrid_fft(crats, 2*nr, 2*nang, /xrev)
      new_crats=rebin(crats,2*nr,2*nang)
      ;;    if (max(crats) gt 1) then stop
      maxrats[*, k] =  max(new_crats,  dim = 2)
      for i = 0, 2*nr-2 do max_gtrats[i, k] =  max(maxrats[i:2*nr-2, k])
      if (k eq 9) then begin
         tend =  tbegin + nsim/10.*(systime(1)-tbegin)
         print,  'First 10 simulations complete. Will finish at: ',  systime(0, tend) 
      endif
   endfor
   rsim = rebin(reform(params[0, *, 0]), 2*nr)
   rsim =  rsim[0:2*nr-2]
   for i = 0, 2*nr-2 do s =  sort(max_gtrats[i, *])
   
   ;;Now for the confidence interval relating to the "detection"
   if (keyword_set(confidence_int)) then begin
      w =  where(rsim gt confidence_int[0] and rsim lt confidence_int[1])
      if (w[0] ne -1) then begin
         m = max(maxrats[w, *], dim = 1)
         s =  sort(m)
         lev = (where(m[s] gt 1./bestp[2]))[0]
         if (lev[0] eq -1) then lev = n_elements(s)
         print,  'Confidence level of detection (%): ',  100.*lev/float(n_elements(s))
      endif
   endif
   ;;Now for the confidence intervals
   ;;specific to the USco paper and
   ;;subsequent papers...
   rnum   = long(confidence_level*nsim)
   rnum99 = long(0.99*nsim)
   if not keyword_set(annuli) then annuli = [10,20,40,80,160,240,320]
   n_annuli = n_elements(annuli)-1
   rat_annuli = fltarr(n_annuli)
   rat99_annuli = fltarr(n_annuli)
   for i=0,n_annuli-1 do begin
      w =  where(rsim ge annuli[i] and rsim lt annuli[i+1])
      if (w[0] ne -1) then begin
         m=max(maxrats[w,*], dim=1)
         s =  sort(m)
         rat_annuli[i]   = m[s[rnum-1]]
         rat99_annuli[i] = m[s[rnum99-1]]
      endif else begin
         rat_annuli[i] = 1.
         rat99_annuli[i]=1.
      endelse
   endfor
   
   ;;Output text... 
   print, 'mrg.oifits file used: ',infile
   ;;set up the formats for printing the output
   f1='(F4.1,"% limits:      MJD    ",'+strcompress(string(n_annuli),/remove_all)+'(I3.0,"-",I3.0," "))'
   f2='(A15, "& ", F7.1, '+strcompress(string(n_annuli),/remove_all)+'("   ", F5.2))'
   f3='("99% only       & ", F7.1, '+strcompress(string(n_annuli),/remove_all)+'("   ", F5.2))'
   ;;set up a temporary array so the limits can be printed easily
   temp_annuli=fltarr(2*(n_annuli))
   for temp_ix=0,n_annuli-1 do temp_annuli[2*temp_ix:2*temp_ix+1]=[annuli[temp_ix],annuli[temp_ix+1]]
   print,confidence_level[*]*100.0,temp_annuli,  format = f1
   print, target_name,t3data[0].mjd,  -2.5*alog10(rat_annuli),    format =f2
   print, t3data[0].mjd,  -2.5*alog10(rat99_annuli),  format =f3
   
   if (p ne '') then begin
      f1='(F4.1,"% limits:      MJD    ",'+strcompress(string(n_annuli),/remove_all)+'(I3.0,"-",I3.0," "))'
      f2='(A15, "& ", F7.1, '+strcompress(string(n_annuli),/remove_all)+'("   ", F5.2),"\\ LIMITS")'
      f3='("99% only       & ", F7.1, '+strcompress(string(n_annuli),/remove_all)+'("   ", F5.2), "\\ L99")'
      printf,1, 'mrg.oifits file used: ',infile
      printf,1, confidence_level*100.0,  format = f1
      printf,1, target_name, t3data[0].mjd, -2.5*alog10(rat_annuli), format =f2
      printf,1, t3data[0].mjd, -2.5*alog10(rat99_annuli), format =f3
   endif  
   annulus_outer = (where(annuli gt bestp[0]))[0]
   if (annulus_outer le 0) then begin
      outstr= 'Unknown Confidence Level'
   endif else begin
      if (bestp[2] gt 1./rat99_annuli[annulus_outer-1]) then outstr = 'Low (<99%) confidence.' $
      else if (bestp[2] gt 1./rat_annuli[annulus_outer-1]) then outstr = string(confidence_level*100, format='("Between 99 and ",F4.1,"% limits.")') $
      else outstr = string(confidence_level*100,'("Detected at ", F4.1, "% limits!")')
   endelse
   print, outstr
   if (p ne '') then printf, 1, outstr
   maglims999=-2.5*alog10(rat_annuli)
   maglims99=-2.5*alog10(rat99_annuli)
   
;;plot estimated detection limits
   rat99_annuli_mags=-2.5*alog10(rat99_annuli)
   annuli_est=(annuli[0:*]+annuli[1:*])/2
   if not keyword_set(quiet_mode) then begin
      window,1
      plot,annuli_est,1./rat99_annuli,yr=[0,max(1./rat99_annuli)],xr=[0,max(annuli+5)],title='Estimated Detection limits',ytitle='Contrast Ratio',xtitle='(Approximate) Separation (mas)'
      oplot,annuli_est,1./rat_annuli,color=100
      ;;and plot the companion
      oplot,[-100,bestp[0]],[0,bestp[2]],psym=5

      ;;plot in magnitudes
      plot,annuli_est,2.5*alog10(1./rat99_annuli),xr=[0,max(annuli+5)],title='Estimated Detection limits',ytitle='Contrast (mags)',xtitle='(Approximate) Separation (mas)'
      oplot,annuli_est,2.5*alog10(1./rat_annuli),color=100
      ;;and plot the companion
      oplot,[-100,bestp[0]],[0,2.5*alog10(bestp[2])],psym=5
   endif
endif else begin
   magslims999=-1
   magslims99=-1
endelse

if (keyword_set(printfile)) then close,  1


; print,'aaa'
if not keyword_set(no_bfinfo) and keyword_set(printfile) and keyword_set(tom_this_doenst_work_for_nirc2) then begin
  pos=strpos(printfile,'.txt',/reverse_search)
  bfinfo=strmid(printfile,0,strlen(printfile)-strlen('.txt'))+'.idlvar'
  wavelength=oiwavelength.eff_wave[0]
  wavelength=*wavelength
  save,target_name,wavelength,extra_error,extra_error_null,bestp,perror,chi2_bin,chi2_null,bestp2,perror2,fmin,infile,maglims999,maglims99,filename=bfinfo
endif
print,'Contrast (mags):',2.5*alog10(bestp[2]),format='(A16,F5.2)'
if keyword_set(quiet_mode) then begin
   modcp=0 ;;to save space!
   if type(quiet_mode) ne type('string') then name='binary_grid_output.idlvar' else name=quiet_mode
   save,/variables,filename=name
endif

end
