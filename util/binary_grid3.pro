;;ACC's edit of binary grid that is able to plot chi squared space
;;for any pair of separation, position angle and contrast.
;;Keywords:
;;  infile: The input oifits file. We don't input the
;;    closure-phase covariance matrix to this program - only the data
;;    that is part of the oifits format.
;;  fix_() : Fixes parameter () to the value given. You must set one
;;  parameter to be fixed! (What about pa=0?)
;;  ()_range: 2 element array =[min,max] of the grid over that parameter
;;  nparams: The number of points in each direction of the grid (need
;;     to make this 3d so we can have different sizes in each direction)
;;  like: An output likelihood map.
;; polar: Plot it in polar coordinates
pro binary_grid3,  infile,  init_crat=init_crat,like=like, extra_error=extra_error,fix_crat=fix_crat,fix_sep=fix_sep,fix_pa=fix_pa,crat_range=crat_range,sep_range=sep_range,pa_range=pa_range,nparams=nparams,polar=polar

if not keyword_set(fix_crat) and not keyword_set(fix_sep) and not keyword_set(fix_pa) then begin
    print,'You didnt fix any params! Fix either crat,sep or pa to do a grid over the others!'
    stop
endif

extract_t3data,  file = infile,  t3data
extract_vis2data,  file = infile,  vis2data
;;Deal with "bad" triangles in a simplistic way.
w = where(t3data.flag)


r =  sqrt(vis2data.u^2 + vis2data.v^2)
maxr =  max(r)

;;automatic ranges if not set:
if keyword_set(pa_range) eq 0 then pa_range=[0.,360.]
if keyword_set(sep_range) eq 0 then sep_range=[rad2mas(1./4/maxr),rad2mas(1./minr)];;assuming we are looking for something not too resolved
if keyword_set(crat_range) eq 0 then init_crat=250.;;arbitrary?

if (w[0] ne -1) then t3data.t3phierr = 1e3;;???
if keyword_set(extra_error) then t3data.t3phierr = sqrt(t3data.t3phierr^2 + extra_error^2)
extract_vis2data,  file = infile,  vis2data
read_oidata,infile,oiarray,oitarget
target_name=oitarget.target
target_name=target_name[0]

;;Enforce a maximum SNR of 4.0 on visibility data. (Mike's comment?)
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

;;get the chi squared array ready
if keyword_set(nparams) eq 0 then nparams=50
chi2_arr=fltarr(nparams,nparams)

;;what are test params? Determine which one is fixed, then make the
;;other two 2D arrays where one varies in one direction, and the other
;;in the other. Make the fixed parameter
if keyword_set(fix_crat) then begin
    min_sep=sep_range[0]
    max_sep=sep_range[1]
    min_pa=pa_range[0]
    max_pa=pa_range[1]
    x=min_sep+(max_sep-min_sep)*findgen(nparams)/(nparams-1)
    y=min_pa+(max_pa-min_pa)*findgen(nparams)/(nparams-1)
    crats=fix_crat+fltarr(nparams,nparams)
    make_2d,x,y,seps,pas
endif else if keyword_set(fix_sep) then begin
    xax='Contrast Ratio'
    yax='Position Angle (deg)'
    min_pa=pa_range[0]
    max_pa=pa_range[1]
    min_crat=crat_range[0]
    max_crat=crat_range[1]
    crats_in=min_crat+(max_crat-min_crat)*findgen(nparams)/(nparams-1)
    seps=fix_sep+fltarr(nparams,nparams)
    y=min_pa+(max_pa-min_pa)*findgen(nparams)/(nparams-1)
    make_2d,crats_in,y,crats,pas 
endif else if keyword_set(fix_pa) then begin
    xax='Separation (mas)'
    yax='Contrast Ratio'
    min_sep=sep_range[0]
    max_sep=sep_range[1]
    min_crat=crat_range[0]
    max_crat=crat_range[1]
    x=min_sep+(max_sep-min_sep)*findgen(nparams)/(nparams-1)
    crats_in=min_crat+(max_crat-min_crat)*findgen(nparams)/(nparams-1)
    pas=fix_pa+fltarr(nparams,nparams)
    make_2d,x,crats_in,seps,crats 
endif

;;Now search for the best chi^2 over the grid...
t1=systime(1)
for i=0,nparams-1 do begin
    for j=0,nparams-1 do begin
        params = [seps[i,j], pas[i,j], crats[i,j],0.01,0.01]
        modelt3 = binary_t3data(params,t3data=t3data)
        cpresid = mod360(modelt3.t3phi - t3data.t3phi)
        chi2_arr[i, j] =  total(cpresid^2/t3data.t3phierr^2) 
    endfor
endfor
print,'Time taken:',systime(1)-t1,' secs'

;;Now get ready for the chi^2 contour plot (You need cgimage!)
if keyword_set(fix_crat) then begin
    xax='Separation (mas)'
    yax='Position Angle (deg)'
    xr=[min_sep,max_sep]
    yr=[min_pa,max_pa]
endif else if keyword_set(fix_sep) then begin
    xax='Contrast Ratio'
    yax='Position Angle (deg)'  
    xr=[min_crat,max_crat]
    yr=[min_pa,max_pa]
endif else if keyword_set(fix_pa) then begin
    xax='Separation (mas)'
    yax='Contrast Ratio'
    xr=[min_sep,max_sep]
    yr=[min_crat,max_crat]
endif

;;what are the best params and model cps?
w=where(chi2_arr eq min(chi2_arr))
w=array_indices(chi2_arr,w)
bestp=[seps[w[0],w[1]], pas[w[0],w[1]], crats[w[0],w[1]],0.01,0.01]
modelt3=binary_t3data(bestp,t3data=t3data)
model_cps=modelt3.t3phi
model_cperr=modelt3.t3phierr


;;get the measured closure phases and errors:
cps=t3data.t3phi
cperr=t3data.t3phierr

;;Bin the closure phases:
range=[min([model_cps,t3data.t3phi]),max([model_cps,t3data.t3phi])]
minclp=min(model_cps)
maxclp=max(model_cps)
diff=(maxclp-minclp)
bins=40.
y=dblarr(bins)
x=y
yerr=y
mn=minclp
mx=minclp+diff/bins
bad_bins=0
;;loop over bins
for i=0,bins-1 do begin
    restart_bin_loop:
    ;;save x position of bin and find points in bin:
    x[i]=(mn+mx)/2
    w=where((model_cps lt mx) and (model_cps gt mn))

    ;;if there is no data in this bin, then go to the next one. Count
    ;;the number of bad bins so we can reduce the size of x and y by
    ;;that amount
    if w[0] eq -1 then begin
       bad_bins+=1
       mn+=diff/bins
       mx+=diff/bins
       goto,restart_bin_loop
    endif
    
    y[i]=(moment(t3data[w].t3phi,sdev=temp))[0]
    yerr[i]=temp

    ;;min and max range of each bin for next loop:
    mn+=diff/bins
    mx+=diff/bins
    
    ;;this will ensure we exit the loop at the correct time
    if i eq (bins-1-bad_bins) then goto, exit_bin_loop
endfor
exit_bin_loop:
x=x[0:bins-1-bad_bins]
y=y[0:bins-1-bad_bins]


;cgimage,chi2_arr,Margin=0.07, /Save,ncolors=2000,/keep_aspect,font=0,xrange=xr,yrange=yr,xtitle=xax,ytitle=yax,position=[0.05,0.05,1.,1.],/axes

;;Make reduced chi^2 equal to 1 by scaling errors, and make sure that
;;our final chi^2 variable has a minimum of ndf (crudely takes into
;;account linear dependence of closure-phases). (From binary_grid2)
ndf =  n_elements(vis2data) -(n_elements(t3data)/float(n_elements(vis2data))*3.+1)
chi2_arr2 = chi2_arr/min(chi2_arr)
chi2_arr2 *= ndf
;;The likelihood map.
;like = exp(-(chi2_arr-ndf)/2);;from binary_grid2
like=exp(-(chi2_arr2-min(chi2_arr2))/2./sqrt(n_elements(t3data)/ndf)/4.);; from binary_grid
cgimage,like,Margin=0.07, /Save,ncolors=2000,/keep_aspect,font=0,xrange=xr,yrange=yr,xtitle=xax,ytitle=yax,position=[0.05,0.05,1.,1.],/axes,title=infile

;;get ra vs dec plots ready:
theta=(reform(pas[0,*]))/180.*!pi
r=reform(seps[*,0])
ra_dec=reverse(transpose(polar_surface(like,r,theta,/grid,bounds=[-200,-200,200,200],/quintic))) ;;need to rotate 90deg anti-clockwise to get RA on x axis

;;!!!!!!!!!!!!!!!
;;Plots
;;!!!!!!!!!!!!!!!

;;chi squared array:
;cgimage,chi2_arr,Margin=0.07, /Save,ncolors=2000,/keep_aspect,font=0,xrange=xr,yrange=yr,xtitle=xax,ytitle=yax,position=[0.05,0.05,1.,1.],/axes

;;measured closure phases v model
;plot,model_cps,cps,psym=5,xtitle='Model Closure Phase',ytitle='Measured Closure Phase',charsize=2
;oplot,[-50,50],[-50,50],color=100

;;plot the binned cps:
;ploterr,x,y,yerr,xtitle='Model',ytitle='Data',charsize=2,psym=5
;oplot,range,range,color=100

;;RA v DEC plots: (use only if generated plots of pa vs sep!)
;polar_contour,transpose(like),theta-!pi/2,r,/fill
;cgimage,ra_dec,Margin=0.07,/Save,ncolors=2000,/keep_aspect,font=0,xrange=[-max(seps),max(seps)],yrange=[-max(seps),max(seps)],xtitle='RA (mas)',ytitle='Dec (mas)',position=[0.05,0.05,1.,1.],/axes,title='Polar'


;;for saving:
;psopen,'~/plots/binary_grid_plots/t_cha_mar10_cpedit_sep_pa',/encapsulated,/bold,/helvetica,xsize=20,ysize=15,/color
;(insert plot commands here)
;psclose
stop
end

