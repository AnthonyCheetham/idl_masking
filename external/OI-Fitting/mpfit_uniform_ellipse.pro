pro mpfit_uniform_ellipse,file=file,vis2data=vis2data,t3data=t3data,$
 params=params,sigma=sigma,vis2model=vis2model,t3model=t3model,$
  fixed=fixed,bestnorm=bestnorm,$
  minvis2err=minvis2err,minvis2perr=minvis2perr,mincperr=mincperr,$
  helpme=helpme,plotit=plotit

; set FIXED =[1,1,1,0] to fix all but the 4th parameter.

if keyword_set(helpme) then begin
 print,'mpfit_uniform_ellipse,file=file,vis2data=vis2data,t3data=t3data,$'
 print,'   params=params,sigma=sigma,vis2model=vis2model,t3model=t3model,$'
 print,'    fixed=fixed, bestnorm=bestnorm,$ '
 print,'    minvis2err=minvis2err,minvis2perr=minvis2perr,$'
 print,'    mincperr=mincperr,/plot,/help'
return
endif


; Accepts FILEname OR vis2data/t3data.
if (keyword_set(file) eq 1) then begin
 extract_vis2data,file=file,vis2data
 extract_t3data,file=file,t3data
endif



; parameters:
; Model of a uniform ellipse
; Params:
;      a(0)=The Average of Major/Minor Axes in units of milliarcseconds
; If diameter is less than zero than invert result
;      a(1)=Axis Ratio (Major/Minor Axis)
;      a(2)=Angle of major Axis  (Degrees -- East of North
;      a(3)=The zero visibility intercept.  Generally should be 1.0

; start params
if keyword_set(params) eq 0 then $
params=[1.5,1.0,0.0,1.0]

if keyword_set(minvis2perr) eq 0 then minvis2perr=0.05; mult
if keyword_set(minvis2err) eq 0 then minvis2err=0.02 ; additive
if keyword_set(mincperr)   eq 0 then mincperr=.5
if n_elements(fixed) ne 4 then fixed=[0,0,0,0]

vis2data0=vis2data
if n_elements(t3data) ne 0 then t3data0=t3data
if (n_elements(vis2data) ne 0) then begin
  vis2data0.vis2err=vis2data0.vis2err>(minvis2perr*vis2data0.vis2data)
  vis2data0.vis2err=vis2data0.vis2err>minvis2err
endif

if (n_elements(t3data) ne 0 ) then $
 t3data0.t3phierr=t3data.t3phierr>mincperr


nparms=n_elements(params)
parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                                  limits:[0.D,0],step:0d0,relstep:0.d0,tied:''}, nparms)

parinfo.fixed=fixed

parinfo(1).limited=[1,0]
parinfo(1).limits=[1.d0,0]
parinfo(0).limited=[1,0]
parinfo(0).limits=[1e-5,0]

;parinfo(2).step=.001
;parinfo(3).relstep=.01
;parinfo(4).step=.01


vis2model0=vis2data0
if n_elements(t3data0) ne 0 then t3model0=t3data0
nv=n_elements(vis2data0)

if n_elements(t3data0) ne 0 then $
fa={$
    vis2data0:vis2data0,t3data0:t3data0,$ 
    vis2model0:vis2model0,t3model0:t3model0} $
else fa={vis2data0:vis2data0, vis2model0:vis2model0}

thebest=params

; First fit w/o letting ratio and PA change.
startparams=params
parinfo(1:2).fixed=[1,1]
startparams(1:2)= [1.0,0]
if (total(parinfo.fixed) ne 4) then begin

  best_params = mpfit('uniform_ellipse_mpfunc', startparams, functargs=fa,bestnorm=bestnorm,$
           parinfo=parinfo,perror=perror,covar=covar,status=status,/quiet)
startparams=best_params
endif
parinfo(1:2).fixed=fixed(1:2)
startparams(1:2)=params(1:2)

best_params = mpfit('uniform_ellipse_mpfunc', startparams, functargs=fa,bestnorm=bestnorm,$
           parinfo=parinfo,perror=perror,covar=covar,status=status,/quiet)
;stop
resid=uniform_ellipse_mpfunc(best_params,dp, vis2data0=vis2data0,$
 t3data0=t3data0,$
  vis2model0=vis2model0,t3model0=t3model0)

bestnorm=total( resid*resid)
    bestnorm=bestnorm/n_elements(resid)
 print,"Min Chi2/DOF: ", bestnorm

best_params(2)=mod360(best_params(2))
if best_params(2) le 0 then best_params(2)=best_params(2)+180.
; guarantees PA 0->180

print,'Best Params:',best_params


vis2model=vis2model0
if n_elements(t3model0) ne 0 then t3model=t3model0

params=best_params
sigma=perror

if (keyword_set(plotit) eq 1) then begin
  plotsym,0,1.0
  ploterror,vis2data.sfu/1e6,sqrt(vis2data.vis2data),$
      sqrt(vis2data.vis2data)*.5*vis2data.vis2err/vis2data.vis2data,$
      tit=vis2data(0).target, xtit="Spatial Frequency (MegaLambda)",$
      ytit="Visibility",psym=8,$
      xr=[0,max(vis2data.sfu)/1e6]
  plotsym,0,/fill,.5
  oplot,vis2model.sfu/1e6,sqrt(vis2model.vis2data),psym=1
endif




end

