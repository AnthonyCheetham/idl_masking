;;program to do double differential calibration of visibilites, phases and clps.
;; Uncertainties are calculated by bootstrapping
;;
;;outputs the clps and vis in big arrays:
;; clps= n_bs x 2 (value, uncertainty) x 2 (0 deg,22 deg)
;; vis = n_bl x 2 (value, uncertainty) x 2 (0 deg,22 deg)
;;
;; OPTIONS:
;;   average_together:
;;     Due to the large parallactic angle changes frequently encountered
;;     during GPI observations, this program calibrates frames in pairs
;;     taken close together in time, and then preserves the individual sky
;;     rotations. To reverse this behaviour, use the average_together
;;     option.
;;
;; This version works based on tsize rather than individual files, since it it very difficult to track which parallactic angle values belong to which files when they are jumbled up.
;; This program makes lots of assumptions about how the files are arranged in the arrays, and may break if this changes substantially.

pro gpi_pol_double_diff_old,cubeinfo_file,target_tsizes,nbootstrap=nbootstrap,average_together=average_together,root_dir=root_dir,save_name=save_name

defsysv, '!ROOT_DIR', exists=exists
if exists then root_dir=!ROOT_DIR
if not keyword_set(root_dir) then begin
 print, 'ERROR: Must set root_dir keyword or !ROOT_DIR system variable.'
 stop
endif

if not keyword_set(nbootstrap) then nbootstrap=1000

restore,cubeinfo_file
mf_file=root_dir+'templates/'+plog.mf_file
restore,mf_file

;; Loop over the target tsizes
for tsize_ix=0,n_elements(target_tsizes)-1 do begin
    target_tsize=target_tsizes[tsize_ix]

    ;;find the relevant files
    targnum=where(olog.tsize(uniq(olog.tsize)) eq target_tsize)
    files=plog.bs_names[targnum*8:targnum*8+7] ;;qbe puts them in blocks of 8

    print,'Calibrating files for target with tsize = ',target_tsize
    ;;first the 0 degree quantities
    restore,files[0]
    u0_0=u
    v0_0=v
    cp0_0=atan(bs_all,/phase)*180/!dpi
    v20_0=v2_all
    p0_0=atan(cvis_all,/phase)*180/!dpi
    restore,files[1]
    u0_45=u
    v0_45=v
    cp0_45=atan(bs_all,/phase)*180/!dpi
    v20_45=v2_all
    p0_45=atan(cvis_all,/phase)*180/!dpi

    ;;The 22 degree quantities
    restore,files[2]
    u22_0=u
    v22_0=v
    cp22_0=atan(bs_all,/phase)*180/!dpi
    v222_0=v2_all
    p22_0=atan(cvis_all,/phase)*180/!dpi
    restore,files[3]
    u22_45=u
    v22_45=v
    cp22_45=atan(bs_all,/phase)*180/!dpi
    v222_45=v2_all
    p22_45=atan(cvis_all,/phase)*180/!dpi

    ;;The 45 degree quantities
    restore,files[4]
    u45_0=u
    v45_0=v
    cp45_0=atan(bs_all,/phase)*180/!dpi
    v245_0=v2_all
    p45_0=atan(cvis_all,/phase)*180/!dpi
    restore,files[5]
    u45_45=u
    v45_45=v
    cp45_45=atan(bs_all,/phase)*180/!dpi
    v245_45=v2_all
    p45_45=atan(cvis_all,/phase)*180/!dpi

    ;;The 68 degree quantities
    restore,files[6]
    u68_0=u
    v68_0=v
    cp68_0=atan(bs_all,/phase)*180/!dpi
    v268_0=v2_all
    p68_0=atan(cvis_all,/phase)*180/!dpi
    restore,files[7]
    u68_45=u
    v68_45=v
    cp68_45=atan(bs_all,/phase)*180/!dpi
    v268_45=v2_all
    p68_45=atan(cvis_all,/phase)*180/!dpi

    ;;calculate some useful quantities
    n0=(size(v20_0))[1] ;; number of files
    n22=(size(v222_0))[1] ;; number of files
    n45=(size(v245_0))[1] ;; number of files
    n68=(size(v268_0))[1] ;; number of files
    sz0=size(v20_0) ;; for keeping track of the other dimensions
    sz22=size(v222_0) ;; for keeping track of the other dimensions
    sz45=size(v245_0) ;; for keeping track of the other dimensions
    sz68=size(v268_0) ;; for keeping track of the other dimensions
    bl=sqrt(u^2+v^2)   ;; baseline length
    az=atan(v,u)       ;; azimuth

    ;; If we're not using average_together, then this program wont work with different numbers of files for each hwp angle
    if not keyword_set(average_together) and ((n0 ne n45) or (n22 ne n68)) then begin
    	print,'!!! Warning !!!'
    	print,'The number of frames in the orthogonal polarization channels is not equal!'
    	print,'Calculating the differential quantities frame by frame will not work.'
    	print,'You can try setting average_together to average over frames instead.'
    	stop
    endif

    nbl=(size(u))[1]
    ncp=(size(cp0_0))[2]

    ;;now, calibrate them to get the means
    vis_vh_0=gpi_pol_cal_v2(v20_0,v20_45,v245_0,v245_45,olog,average_together=average_together)
    vis_vh_22=gpi_pol_cal_v2(v222_0,v222_45,v268_0,v268_45,olog,average_together=average_together)

    cp_vh_0=gpi_pol_cal_ph(cp0_0,cp0_45,cp45_0,cp45_45,average_together=average_together)
    cp_vh_22=gpi_pol_cal_ph(cp22_0,cp22_45,cp68_0,cp68_45,average_together=average_together)

    ph_vh_0=gpi_pol_cal_ph(p0_0,p0_45,p45_0,p45_45,average_together=average_together)
    ph_vh_22=gpi_pol_cal_ph(p22_0,p22_45,p68_0,p68_45,average_together=average_together)

    ;; now calculate the uncertainties by bootstrapping. Pick nbootstrap samples
    ;; each of length ndata.
    ;; Unfortunately in making this data-dimension independent, it is now much harder to read.
    ;; But sz[-1]/sz[2] is the number of clps in each frame ()
stop
    samples0 =fix(sz0[-1]/sz0[2]*randomu(seed1,n0,nbootstrap))
    samples22=fix(n22*randomu(seed2,n22,nbootstrap))
    samples45=fix(n45*randomu(seed3,n45,nbootstrap))
    samples68=fix(n68*randomu(seed4,n68,nbootstrap))

    boots_v20_0=reform(v20_0[samples0,*],n0,nbootstrap,nbl)
    boots_v20_45=reform(v20_45[samples0,*],n0,nbootstrap,nbl)
    boots_v222_0=reform(v222_0[samples22,*],n22,nbootstrap,nbl)
    boots_v222_45=reform(v222_45[samples22,*],n22,nbootstrap,nbl)
    boots_v245_0=reform(v245_0[samples45,*],n45,nbootstrap,nbl)
    boots_v245_45=reform(v245_45[samples45,*],n45,nbootstrap,nbl)
    boots_v268_0=reform(v268_0[samples68,*],n68,nbootstrap,nbl)
    boots_v268_45=reform(v268_45[samples68,*],n68,nbootstrap,nbl)
    vis_vh_0_err=gpi_pol_cal_v2(boots_v20_0,boots_v20_45,boots_v245_0,$
              boots_v245_45,/uncertainties,average_together=average_together)
    vis_vh_22_err=gpi_pol_cal_v2(boots_v222_0,boots_v222_45,boots_v268_0,$
              boots_v268_45,/uncertainties,average_together=average_together)

    boots_cp0_0=reform(cp0_0[samples0,*],n0,nbootstrap,ncp)
    boots_cp0_45=reform(cp0_45[samples0,*],n0,nbootstrap,ncp)
    boots_cp22_0=reform(cp22_0[samples22,*],n22,nbootstrap,ncp)
    boots_cp22_45=reform(cp22_45[samples22,*],n22,nbootstrap,ncp)
    boots_cp45_0=reform(cp45_0[samples45,*],n45,nbootstrap,ncp)
    boots_cp45_45=reform(cp45_45[samples45,*],n45,nbootstrap,ncp)
    boots_cp68_0=reform(cp68_0[samples68,*],n68,nbootstrap,ncp)
    boots_cp68_45=reform(cp68_45[samples68,*],n68,nbootstrap,ncp)
    cp_vh_0_err=gpi_pol_cal_ph(boots_cp0_0,boots_cp0_45,boots_cp45_0,boots_cp45_45,/uncertainties,average_together=average_together)
    cp_vh_22_err=gpi_pol_cal_ph(boots_cp22_0,boots_cp22_45,boots_cp68_0,boots_cp68_45,/uncertainties,average_together=average_together)

    boots_ph0_0=reform(p0_0[samples0,*],n0,nbootstrap,nbl)
    boots_ph0_45=reform(p0_45[samples0,*],n0,nbootstrap,nbl)
    boots_ph22_0=reform(p22_0[samples22,*],n22,nbootstrap,nbl)
    boots_ph22_45=reform(p22_45[samples22,*],n22,nbootstrap,nbl)
    boots_ph45_0=reform(p45_0[samples45,*],n45,nbootstrap,nbl)
    boots_ph45_45=reform(p45_45[samples45,*],n45,nbootstrap,nbl)
    boots_ph68_0=reform(p68_0[samples68,*],n68,nbootstrap,nbl)
    boots_ph68_45=reform(p68_45[samples68,*],n68,nbootstrap,nbl)

    ph_vh_0_err=gpi_pol_cal_ph(boots_ph0_0,boots_ph0_45,boots_ph45_0,boots_ph45_45,/uncertainties,average_together=average_together)
    ph_vh_22_err=gpi_pol_cal_ph(boots_ph22_0,boots_ph22_45,boots_ph68_0,boots_ph68_45,/uncertainties,average_together=average_together)

    if not keyword_set(average_together) then begin
        vis_vh_0_err=rebin(reform(vis_vh_0_err,[1,nbl]),[n0,nbl])
        vis_vh_22_err=rebin(reform(vis_vh_22_err,[1,nbl]),[n22,nbl])
        cp_vh_0_err=rebin(reform(cp_vh_0_err,[1,ncp]),[n0,ncp])
        cp_vh_22_err=rebin(reform(cp_vh_22_err,[1,ncp]),[n22,ncp])
        ph_vh_0_err=rebin(reform(ph_vh_0_err,[1,nbl]),[n0,nbl])
        ph_vh_22_err=rebin(reform(ph_vh_22_err,[1,nbl]),[n22,nbl])
    endif else begin
        vis_vh_0=reform(vis_vh_0,[1,nbl])
        vis_vh_22=reform(vis_vh_22,[1,nbl])
        cp_vh_0=reform(cp_vh_0,[1,ncp])
        cp_vh_22=reform(cp_vh_22,[1,ncp])
        ph_vh_0=reform(ph_vh_0,[1,nbl])
        ph_vh_22=reform(ph_vh_22,[1,nbl])
        vis_vh_0_err=reform(vis_vh_0_err,[1,nbl])
        vis_vh_22_err=reform(vis_vh_22_err,[1,nbl])
        cp_vh_0_err=reform(cp_vh_0_err,[1,ncp])
        cp_vh_22_err=reform(cp_vh_22_err,[1,ncp])
        ph_vh_0_err=reform(ph_vh_0_err,[1,nbl])
        ph_vh_22_err=reform(ph_vh_22_err,[1,nbl])
    endelse

    ;;and put them together into big arrays
    vis=[[[vis_vh_0]],[[vis_vh_22]]]
    vis_err=[[[vis_vh_0_err]],[[vis_vh_22_err]]]
    cp=[[[cp_vh_0]],[[cp_vh_22]]]
    cp_err=[[[cp_vh_0_err]],[[cp_vh_22_err]]]
    ph=[[[ph_vh_0]],[[ph_vh_22]]]
    ph_err=[[[ph_vh_0_err]],[[ph_vh_22_err]]]

    ;;the old way:
    ;vis=[[[vis_vh_0],[vis_vh_0_err]],[[vis_vh_22],[vis_vh_22_err]]]
    ;cp=[[[cp_vh_0],[cp_vh_0_err]],[[cp_vh_22],[cp_vh_22_err]]]
    ;phase=[[[ph_vh_0],[ph_vh_0_err]],[[ph_vh_22],[ph_vh_22_err]]]

    ;;rotate the uv coords by the parallactic angle
    ;; the un-rotated uv coords are all the same, so take the first one
    u=u0_0
    v=v0_0

    ;;find the relevant elements in the olog arrays
    w_tsize=where(olog.tsize eq target_tsize)
    parangs=olog.pas[w_tsize]
    hwps=olog.hwps[w_tsize]

    if keyword_set(average_together) then begin
        ;; if we averaged frames together, then average the parallactic angles together
        parang=mean(parangs)

        ;;now rotate the uv coords
        unew =  olog.uflip*u*cos(parang*!pi/180.0)  + v*sin(parang*!pi/180.0)
        vnew = -olog.uflip*u*sin(parang*!pi/180.0) +  v*cos(parang*!pi/180.0)

        u=rebin(reform(unew,[1,nbl,1]),[1,nbl,2])
        v=rebin(reform(vnew,[1,nbl,1]),[1,nbl,2])

        bs_u=dblarr(1,n_bispect,3,2)
        bs_v=dblarr(1,n_bispect,3,2)

        print,'  Differential quantities calculated with a single parallactic angle.'
        print,'  Max parallactic angle change between combined frames:',max(parangs)-min(parangs),'deg'

    endif else begin
        ;;otherwise rotate them frame-pair by frame-pair

        sort_hwps=hwps(sort(hwps))
        uniq_hwps=sort_hwps(uniq(sort_hwps))
        w0=where(hwps eq uniq_hwps[0]) ;; this is the index of hwps and parangs corresponding to 0 degrees (or the smallest hwp rotation)
        w22=where(hwps eq uniq_hwps[1])
        w45=where(hwps eq uniq_hwps[2])
        w68=where(hwps eq uniq_hwps[3])

        ;;loop over frames for the 0-45 deg quantities
        delpas=[]
        u0=dblarr(n0,nbl)
        v0=dblarr(n0,nbl)
        for ix=0,n0-1 do begin
            parang=(parangs[w0[ix]]+parangs[w45[ix]])/2
            delpa=abs(parangs[w0[ix]]-parangs[w45[ix]])
            delpas=[delpas,delpa]

            ;;now rotate the uv coords
            u0[ix,*] =  olog.uflip*u*cos(parang*!dpi/180.0)  + v*sin(parang*!dpi/180.0)
            v0[ix,*] = -olog.uflip*u*sin(parang*!dpi/180.0) +  v*cos(parang*!dpi/180.0)
        endfor

        ;;loop over frames for the 22-68 deg quantities
        u22=dblarr(n22,nbl)
        v22=dblarr(n22,nbl)
        for ix=0,n22-1 do begin
            parang=(parangs[w22[ix]]+parangs[w68[ix]])/2
            delpa=abs(parangs[w22[ix]]-parangs[w68[ix]])
            delpas=[delpas,delpa]

            ;;now rotate the uv coords
            u22[ix,*] =  olog.uflip*u*cos(parang*!pi/180.0)  + v*sin(parang*!pi/180.0)
            v22[ix,*] = -olog.uflip*u*sin(parang*!pi/180.0) +  v*cos(parang*!pi/180.0)
        endfor
        print,'  Differential quantities calculated frame-by-frame in pairs.'
        print,'  Max parallactic angle change between combined frames:',max(delpas),'deg'
        print,'  Parallactic angle change during observations:',max(parangs)-min(parangs),'deg'

        ;;Put the coords together
        u=[[[u0]],[[u22]]]
        v=[[[v0]],[[v22]]]

        bs_u=dblarr(n0,n_bispect,3,2)
        bs_v=dblarr(n0,n_bispect,3,2)
    endelse

    ;; and combine the uv coords for the closure phases
    bs_u[*,*,0,*]=u[*,bs2bl_ix[0,*],*]
    bs_u[*,*,1,*]=u[*,bs2bl_ix[1,*],*]
    bs_u[*,*,2,*]=u[*,bs2bl_ix[2,*],*]
    bs_v[*,*,0,*]=v[*,bs2bl_ix[0,*],*]
    bs_v[*,*,1,*]=v[*,bs2bl_ix[1,*],*]
    bs_v[*,*,2,*]=v[*,bs2bl_ix[2,*],*]

    ; Save everything as an idlvar file.
    this_save_name='pol'+string(targnum,format='(I04)')+'.idlvar'
    print,'  Saving as '+this_save_name
    save,vis,vis_err,cp,cp_err,ph,ph_err,u,v,bs_u,bs_v,file=this_save_name

    ;;now the merged quantities
    if tsize_ix eq 0 then begin
        vis_mrg=vis
        vis_err_mrg=vis_err
        cp_mrg=cp
        cp_err_mrg=cp_err
        ph_mrg=ph
        ph_err_mrg=ph_err
        u_mrg=u
        v_mrg=v
        bs_u_mrg=bs_u
        bs_v_mrg=bs_v
    endif else begin
        vis_mrg=[vis_mrg,vis]
        vis_err_mrg=[vis_err_mrg,vis_err]
        cp_mrg=[cp_mrg,cp]
        cp_err_mrg=[cp_err_mrg,cp_err]
        ph_mrg=[ph_mrg,ph]
        ph_err_mrg=[ph_err_mrg,ph_err]
        u_mrg=[u_mrg,u]
        v_mrg=[v_mrg,v]
        bs_u_mrg=[bs_u_mrg,bs_u]
        bs_v_mrg=[bs_v_mrg,bs_v]
    endelse

    print,''
endfor
;;now save the merged quantities
vis=vis_mrg
vis_err=vis_err_mrg
cp=cp_mrg
cp_err=cp_err_mrg
ph=ph_mrg
ph_err=ph_err_mrg
u=u_mrg
v=v_mrg
bs_u=bs_u_mrg
bs_v=bs_v_mrg
targnum=where(olog.tsize(uniq(olog.tsize)) eq target_tsizes[0])
if not keyword_set(save_name) then save_name='pol'+string(targnum,format='(I04)')+'mrg.idlvar'
print,'Saving merged file as '+save_name
save,vis,vis_err,cp,cp_err,ph,ph_err,u,v,bs_u,bs_v,file=save_name
end
