;;ACC's edit of binary grid (from binary grid2) that should be able to
;;do a 3d grid search (hence cube) over separation, position angle and
;;contrast ratio. Best fit parameters and uncertainties are calculated
;;from the likelihood distributions marginalized over remaining variables.
;;Keywords:
;;  infile: The input oifits file. We don't input the
;;    closure-phase covariance matrix to this program - only the data
;;    that is part of the oifits format.
;;  init_crat: The contrast ratio at which the grid search is done. NB
;;    in the high contrast regime, you'll get the same best
;;    solution no matter what init_crat is chosen.
;;  boxsize: The size of the grid box in milli-arcsec
;;  like: An output likelihood map.
;;  minsep: The resolution of the grid.
;; Author: ACC 2012
pro binary_cube,  infile,like=like, extra_error=extra_error,crat_range=crat_range,sep_range=sep_range,pa_range=pa_range,nparams=nparams,save_name=save_name

extract_t3data,  file = infile,  t3data
extract_vis2data,  file = infile,  vis2data
;;Deal with "bad" triangles in a simplistic way.
w = where(t3data.flag)

r =  sqrt(vis2data.u^2 + vis2data.v^2)
maxr =  max(r)
minr=min(r)

;;automatic ranges if not set:
if keyword_set(pa_range) eq 0 then pa_range=[0.,360.]
if keyword_set(sep_range) eq 0 then sep_range=[rad2mas(1./4/maxr),rad2mas(1./minr)];;assuming we are looking for something not too resolved
if keyword_set(crat_range) eq 0 then crat_range=[20,400];;arbitrary?

;if (w[0] ne -1) then t3data.t3phierr = 1e3 ;;commenting out is a hack!
if keyword_set(extra_error) then t3data.t3phierr = sqrt(t3data.t3phierr^2 + extra_error^2)
extract_vis2data,  file = infile,  vis2data
read_oidata,infile,oiarray,oitarget
target_name=oitarget.target
target_name=target_name[0]

;;Enforce a maximum SNR of 4.0 on visibility data.
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

r =  sqrt(vis2data.u^2 + vis2data.v^2)
maxr =  max(r)
if (keyword_set(minsep) eq 0) then minsep =  rad2mas(1./8/maxr)

;;get the chi squared array ready
if keyword_set(nparams) eq 0 then nparams=50
if (size(nparams))[1] lt 3 then nparams=[nparams[0],nparams[0],nparams[0]]
chi2_arr=fltarr(nparams[0],nparams[1],nparams[2])

;;what are test params? Determine which one is fixed, then make the
;;other two 2D arrays where one varies in one direction, and the other
;;in the other
min_sep=sep_range[0]
max_sep=sep_range[1]
min_pa=pa_range[0]
max_pa=pa_range[1]
min_crat=crat_range[0]
max_crat=crat_range[1]
seps=min_sep+(max_sep-min_sep)*findgen(nparams[0])/(nparams[0]-1)
pas=min_pa+(max_pa-min_pa)*findgen(nparams[1])/(nparams[1]-1)
crats=min_crat+(max_crat-min_crat)*findgen(nparams[2])/(nparams[2]-1)

;;Now search for the best chi^2 over the cube...
for i=0,nparams[0]-1 do begin
    t1=systime(1)
    for j=0,nparams[1]-1 do begin
        for k=0,nparams[2]-1 do begin
            params = [seps[i], pas[j], crats[k],0.01,0.01]
            modelt3 = binary_t3data(params,t3data=t3data)
            cpresid = mod360(modelt3.t3phi - t3data.t3phi)
            chi2_arr[i, j,k] =  total(cpresid^2/t3data.t3phierr^2) 
        endfor
    endfor
    print,'done',i+1,' of',nparams[0],' in',systime(1)-t1,' secs'
endfor

;;Now get ready for the chi^2 contour plot (Make sure you have
;;cgimage!)
xax='Separation (mas)'
yax='Position Angle (deg)'
xr=[min_sep,max_sep]
yr=[min_pa,max_pa]

;cgimage,chi2_arr,Margin=0.07, /Save,ncolors=2000,/keep_aspect,font=0,xrange=xr,yrange=yr,xtitle=xax,ytitle=yax,position=[0.05,0.05,1.,1.],/axes

;;Make reduced chi^2 equal to 1 by scaling errors, and make sure that
;;our final chi^2 variable has a minimum of ndf (crudely takes into
;;account linear dependence of closure-phases)
ndf =  n_elements(vis2data) -(n_elements(t3data)/float(n_elements(vis2data))*3.+1)
chi2_arr2 = chi2_arr/min(chi2_arr)
chi2_arr2 *= ndf
;;The likelihood map.
like = exp(-(chi2_arr2-ndf)/2)

;for i=0,nparams[2]-1 do begin
;    cgimage,like[*,*,i],Margin=0.07, /Save,ncolors=2000,/keep_aspect,font=0,xrange=xr,yrange=yr,xtitle=xax,ytitle=yax,position=[0.05,0.05,1.,1.],/axes,title=crats[i]
;    wait,1
;endfor

;;calculate the best model closure phases (THIS IS WRONG! We want to use the marginalized likelihood distributions, not the chi squared!)
;best=array_indices(chi2_arr,where(chi2_arr eq min(chi2_arr)))
;bestp=[[seps[best[0]],pas[best[1]],crats[best[2]]],0.01,0.01]
;modelt3 = binary_t3data(bestp,t3data=t3data)

;;plots!
;window,0
;plot,modelt3.t3phi,t3data.t3phi,psym=5,charsize=1.5,xtitle='Model Closure Phase (deg)',ytitle='Measured Closure Phase (deg)'
;ploterr,modelt3.t3phi,t3data.t3phi,t3data.t3phierr,psym=5,charsize=1.5,xtitle='Model Closure Phase (deg)',ytitle='Measured Closure Phase (deg)'
;oplot,[-360,360],[-360,360],color=100


;resids=mod360(modelt3.t3phi - t3data.t3phi)
;window,1
;plot,modelt3.t3phi,resids,psym=5,charsize=1.5,xtitle='Model Closure Phase (deg)',ytitle='Closure Phase residual (deg)'
;oplot,[-360,360],[0,0],color=100
;y=x
;image_cont_deluxe, -like, /noc, /asp, xv=-x, yv=y, ytit='Dec (mas)', xtit='RA (mas)'
;oplot, [0], [0], psym=2 

;;Marginalise over each variable:
sep_pa=total(like,3)
sep_crat=total(like,2)
pa_crat=total(like,1)

sep_like=total(sep_pa,2) ;;it doesnt matter which order the totalling is done
pa_like=total(sep_pa,1)
crat_like=total(sep_crat,1)

;;calculate the best parameters from the maximum of the marginalized likelihood distributions:
bestp=[seps[where(sep_like eq max(sep_like))],pas[where(pa_like eq max(pa_like))],crats[where(crat_like eq max(crat_like))]]

;;plots
;;sep vs pa:
;image_cont,sep_pa,/n,/a,xv=seps,yv=pas,tit='Sep vs Pa',xtit='Separation (mas)',ytit='PA (deg)'
;;sep vs crat:
;image_cont,sep_crat,/n,/a,xv=seps,yv=crats,tit='Sep vs crat',xtit='Separation (mas)',ytit='Contrast Ratio'
;;pa vs crat:
;image_cont,pa_crat,/n,/a,xv=pas,yv=crats,tit='Pa vs crat',xtit='PA (deg)',ytit='Contrast Ratio'

;;marginalised plots over each variable:
;;best PA:
;plot,pas,pa_like,xtitle='Position Angle',ytitle='Likelihood'
;;best crat:
;plot,crats,crat_like,xtitle='Contrast Ratio',ytitle='Likelihood'
;;best sep:
;plot,seps,sep_like,xtitle='Separation',ytitle='Likelihood'

;;contour plots over every pair of variables
;for i=0,149 do begin temp2=reform(like[*,*,i]) & temp2[0,0]=1 & image_cont,temp2,/n,/a,xv=seps,yv=pas,tit=string(i) & wait,0.05
;for i=0,149 do begin temp2=reform(like[*,i,*]) & temp2[0,0]=1 & image_cont,temp2,/n,/a,xv=seps,yv=crats,tit=string(i) & wait,0.05
;for i=0,149 do begin temp2=reform(like[i,*,*]) & temp2[0,0]=1 & image_cont,temp2,/n,/a,xv=pas,yv=crats,tit=string(i) & wait,0.05

;;find best pa and errors:
wpa=where(pa_like eq max(pa_like))
low_pa=interpol(pas[0:wpa],pa_like[0:wpa],[exp(-1)])
high_pa=interpol(pas[wpa:*],pa_like[wpa:*],[exp(-1)])
w=array_indices(chi2_arr,where(chi2_arr eq min(chi2_arr)))

;;find confidence intervals. We can find it from the extent of the 68.3% region of the 3d likelihood function or the 1d likelihood function.
confidence_3d=0;;use the 3d likelihood function (i.e. no marginalization)
if keyword_set(confidence_3d) then begin
   x=findgen(500)/499
   y=0*x
   for i=0,499 do begin w=where(like gt x[i]) &  if w[0] ne -1 then y[i]=total(like[w])
   endfor
   y/=total(like) ;;y[i] is the % of points above likelihood of x[i]
   
   ;;take the 68.3% range (1 sigma for 1D Gaussian)
   h=interpol(x,y,[0.683]) ;;h is the likelihood at which 68.3% lies above.
   w_ind=array_indices(like,where(like gt h[0]))

   ;;And print!
   print,'      Min          Best          Max'
   print,seps(min(w_ind[0,*])),bestp[0],seps(max(w_ind[0,*]))
   print,pas(min(w_ind[1,*])),bestp[1],pas(max(w_ind[1,*]))
   print,crats(min(w_ind[2,*])),bestp[2],crats(max(w_ind[2,*]))
   ;;OR
   print,'      Best           +           -'
   print,bestp[0],seps(max(w_ind[0,*]))-bestp[0],bestp[0]-seps(min(w_ind[0,*]))
   print,bestp[1],pas(max(w_ind[1,*]))-bestp[1],bestp[1]-pas(min(w_ind[1,*]))
   print,bestp[2],crats(max(w_ind[2,*]))-bestp[2],bestp[2]-crats(min(w_ind[2,*]))
   
endif else begin ;;otherwise use the 1d functions (like in the T Cha work)
   lims=fltarr(3,2)
   n=1000
   x=findgen(n)/(n-1)/2
   norm_sep_like=sep_like/total(sep_like)
   norm_pa_like=pa_like/total(pa_like)
   norm_crat_like=crat_like/total(crat_like)
   
   cumulative_crat=0*x
   cumulative_pa=0*x
   cumulative_sep=0*x
   for i=0,n-1 do begin
      wsep=where(norm_sep_like gt x[i])
      wpa=where(norm_pa_like gt x[i])
      wcrat=where(norm_crat_like gt x[i])
      if wsep[0] ne -1 then cumulative_sep[i]=total(norm_sep_like[wsep])
      if wpa[0] ne -1 then cumulative_pa[i]=total(norm_pa_like[wpa])
      if wcrat[0] ne -1 then cumulative_crat[i]=total(norm_crat_like[wcrat])
   endfor
   
   ;;take the 68.3% range (1 sigma for 1D Gaussian)
   hsep=interpol(x,cumulative_sep,[0.683]) ;;h is the likelihood at which 68.3% lies above.
   hpa=interpol(x,cumulative_pa,[0.683]) ;;h is the likelihood at which 68.3% lies above.
   hcrat=interpol(x,cumulative_crat,[0.683]) ;;h is the likelihood at which 68.3% lies above.
    w_ind_sep=array_indices(norm_sep_like,where(norm_sep_like gt hsep[0]))
    w_ind_pa=array_indices(norm_pa_like,where(norm_pa_like gt hpa[0]))
    w_ind_crat=array_indices(norm_crat_like,where(norm_crat_like gt hcrat[0]))

    lims[0,*]=[min(seps(w_ind_sep)),max(seps(w_ind_sep))]
    lims[1,*]=[min(pas(w_ind_pa)),max(pas(w_ind_pa))]
    lims[2,*]=[min(crats(w_ind_crat)),max(crats(w_ind_crat))]
   ;;And print!
   print,'      Min          Best          Max'
   print,lims[0,0],bestp[0],lims[0,1]
   print,lims[1,0],bestp[1],lims[1,1]
   print,lims[2,0],bestp[2],lims[2,1]
   ;;OR
   print,'      Best           +           -'
   print,bestp[0],lims[0,1]-bestp[0],bestp[0]-lims[0,0]
   print,bestp[1],lims[1,1]-bestp[1],bestp[1]-lims[1,0]
   print,bestp[2],lims[2,1]-bestp[2],bestp[2]-lims[2,0]
endelse






if keyword_set(save_name) then save,filename=save_name,/variables

stop
end

;binary_cube,infile,crat_range=[100,400],sep_range=[30,200],pa_range=[70,110],nparams=[150,150,
;save,filename='~/code/tcha/mar10_binarycube.idlvar',/all
;binary_cube,infile,crat_range=[120,450],sep_range=[30,200],pa_range=[70,110],nparams=150,save_name='~/code/tcha/mar10_binarycube.idlvar'
;acc_psopen,'~/plots/cube_plots/tcha_mar10_pa_crat'
