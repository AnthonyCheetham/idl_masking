;;;A program to calibrate ZIMPOL differential imaging (not polarimetry) data
;;
;; INPUTS:
;;    cubeinfo_file: the output cubeinfo file from the IDL pipeline
;;    target_tsizes: 2 x N array of tsize values to calibrate together (i.e. with inverted filters)
;;    wavelength_ix: The index of the channel that we care about. The other one will be treated as a PSF reference.

pro calibrate_zimpol_diff_im,cubeinfo_file,target_tsizes,wavelength_ix=wavelength_ix,nbootstrap=nbootstrap,save_python=save_python,root_dir=root_dir,save_dir=save_dir


defsysv, '!ROOT_DIR', exists=exists
if exists then root_dir=!ROOT_DIR
if not keyword_set(root_dir) then begin
 print, 'ERROR: Must set root_dir keyword or !ROOT_DIR system variable.'
 stop
endif

flag_val = 0 ;; False in IDL. (python uses 70)

if not keyword_set(nbootstrap) then nbootstrap=1000
if not keyword_set(wavelength_ix) then wavelength_ix=0

if not keyword_set(save_dir) then save_dir = './'

restore,cubeinfo_file
mf_file=root_dir+'templates/'+plog.mf_file
restore,mf_file

;; Set up the arrays to save the final values
oif_to_merge = []

;; Loop over the target tsizes
for tsize_ix=0,(n_elements(target_tsizes)/2)-1 do begin

    target0_tsize = target_tsizes[0,tsize_ix]
    target1_tsize = target_tsizes[1,tsize_ix]

    ;;find the relevant files
    targ0_ix = where(olog.tsize eq target0_tsize)
    targ1_ix = where(olog.tsize eq target1_tsize)
    files0 = plog.bs_names[targ0_ix]
    files1 = plog.bs_names[targ1_ix]


    print,'Calibrating files for target with tsize = ',target0_tsize,'with tsize = ',target1_tsize
    print,'Calibration is applied in pairs of files.'

    if n_elements(files0) ne n_elements(files1) then begin
        print,'calibrate_zimpol_diff_im assumes that files are calibrated in pairs! But the number of files in each configuration is not the same!'
        stop
    endif

    ;; Loop through the files
    for file_ix=0,n_elements(files1)-1 do begin
        print,'  Calibrating files: ',files0[file_ix],' ',files1[file_ix]
      
          ;; Load both files into memory
        restore,files0[file_ix]
        u_0=u
        v_0=v
        cp_0=atan(bs_all,/phase)*180/!dpi
        v2_0=v2_all
        p_0=atan(cvis_all,/phase)*180/!dpi
        restore,files1[file_ix]
        u_1=u
        v_1=v
        cp_1=atan(bs_all,/phase)*180/!dpi
        v2_1=v2_all
        p_1=atan(cvis_all,/phase)*180/!dpi


        ;;; Bootstrap to get the values and uncertainties!!!
        sz_cp = size(cp_0)
        sz_cp = sz_cp[1:sz_cp[0]]
        sz_v2 = size(v2_0)
        sz_v2 = sz_v2[1:sz_v2[0]]
        n_frames_0 = (size(cp_0))[1]
        n_frames_1 = (size(cp_1))[1]
        boot_ix_0 = fix(n_frames_0*randomu(seed0,n_frames_0,nbootstrap))
        boot_ix_1 = fix(n_frames_1*randomu(seed1,n_frames_1,nbootstrap))
        cp_0_boot = reform(cp_0[boot_ix_0,*,*],[n_frames_0,nbootstrap,sz_cp[1:*]])
        cp_1_boot = reform(cp_1[boot_ix_1,*,*],[n_frames_1,nbootstrap,sz_cp[1:*]])
        v2_0_boot = reform(v2_0[boot_ix_0,*,*],[n_frames_0,nbootstrap,sz_v2[1:*]])
        v2_1_boot = reform(v2_1[boot_ix_1,*,*],[n_frames_1,nbootstrap,sz_v2[1:*]])


        ;; Now calibrate them with the first differential step:
        ;;(note ~wavelength_ix should be the opposite channel)
        cp_0_cal_boot = cp_0_boot[*,*,*,wavelength_ix] - cp_0_boot[*,*,*,~wavelength_ix]
        cp_1_cal_boot = cp_1_boot[*,*,*,wavelength_ix] - cp_1_boot[*,*,*,~wavelength_ix]
        v2_0_cal_boot = v2_0_boot[*,*,*,wavelength_ix]/v2_0_boot[*,*,*,~wavelength_ix]
        v2_1_cal_boot = v2_1_boot[*,*,*,wavelength_ix]/v2_1_boot[*,*,*,~wavelength_ix]

        ;; Take the mean over frames
        cp_0_cal_boot_mean = mean(cp_0_cal_boot,dim=1)
        cp_1_cal_boot_mean = mean(cp_1_cal_boot,dim=1)
        v2_0_cal_boot_mean = mean(v2_0_cal_boot,dim=1)
        v2_1_cal_boot_mean = mean(v2_1_cal_boot,dim=1)
        

        ;; And apply the second differential step
        ;; This should be + and * if the signal channel is wavelength_ix in both cubes.
        ;; If the channels are swapped, then we need to use - and /
        cp_boot = (cp_0_cal_boot_mean + cp_1_cal_boot_mean)/2
        v2_boot = sqrt(v2_0_cal_boot_mean * v2_1_cal_boot_mean) ;; let's keep them as square visibilities...

        ;; Now use this to get the mean value and uncertainties
        cp_mean_boot = median(cp_boot,dim=1)
        cp_uncert_boot = stddev(cp_boot,dim=1)

        v2_mean_boot = median(v2_boot,dim=1)
        v2_uncert_boot = stddev(v2_boot,dim=1)

        plothist,cp_0_boot[*,*,0],col='red',bin=0.1,peak=1
        plothist,cp_0_cal_boot,/overplot,col='green',bin=0.1,peak=1
        plothist,cp_boot,/overplot,col='black',bin=0.1,peak=1

        print,'       Init. scatter     First calib    Second calib'
        print,'    ',stdev(cp_0_boot[*,*,0]),stdev(cp_0_cal_boot),stdev(cp_boot)
        print,'       RMS clp           Median uncertainty'
        print,'    ',sqrt(mean(cp_mean_boot^2)),median(cp_uncert_boot)

        parang0 = olog.pa[targ0_ix[file_ix]] - 0.5*olog.del_pa[targ0_ix[file_ix]]
        parang1 = olog.pa[targ1_ix[file_ix]] - 0.5*olog.del_pa[targ1_ix[file_ix]]
        print,'    Parallactic angle difference between files:',parang1-parang0

        ; stop

        thepa = (parang0 +parang1) / 2.
        u1 =  olog.uflip*u_0*cos(thepa*!pi/180.0)  + v_0*sin(thepa*!pi/180.0)
        v1 = -olog.uflip*u_0*sin(thepa*!pi/180.0) +  v_0*cos(thepa*!pi/180.0)

        ;;;;;;;;;;;;;;
        ;; Save as oifits
        ;;;;;;;;;;;;;;

        source_str = olog.source_name[targ0_ix[file_ix]]
        p = strpos(source_str, '  ')
        while (p gt 0) do begin
            source_str =  strmid(source_str, 0, p)+strmid(source_str, p+1)
            p = strpos(source_str, '  ')
        endwhile
  	  ;Now, replace spaces with underscores.
        p = strpos(source_str, ' ')
        while (p gt 0) do begin
            source_str =  strmid(source_str, 0, p)+'_'+strmid(source_str, p+1)
            p = strpos(source_str, ' ')
        endwhile
   	    ;Finally, get rid of any / characters
        p = strpos(source_str, '/')
        while (p gt 0) do begin
            source_str =  strmid(source_str, 0, p)+strmid(source_str, p+1)
            p = strpos(source_str, '/')
        endwhile
        if (source_str ne '') then source_str =  source_str+'_'
        oif_name =  source_str + olog.cube_fname[targ0_ix[file_ix], 1] + '_diff.oifits'

        ;; This is all copied from calibrate_v2_cp.pro, and edited slightly

    ;First, the oiarray: not used properly here (needs one row per hole)
        define_oiarray, oiarray_unit
        oiarray = replicate(oiarray_unit, 1)
        oiarray[0].extver = 1
        oiarray[0].arrname =olog.mask
        oiarray[0].frame = "GEOCENTRIC"
        oiarray[0].arrayx =  10.0
        oiarray[0].arrayy =  20.0
        oiarray[0].arrayz =  30.0
        oiarray[0].tel_name = ["Dummy Table"]
        oiarray[0].sta_name = ["Dummy Table"]
        oiarray[0].sta_index = [0]
        oiarray[0].diameter = .45
        oiarray[0].staxyz = [0., 0., 0.]
        
    ;Next the target star - only use 1 target star per save file.
        define_oitarget, oitarget_unit
        oitarget = replicate(oitarget_unit, 1)
        oitarget[0].target_id = 0
        oitarget[0].target = olog.source_name[targ0_ix[file_ix]]
        oitarget[0].raep0 = olog.ra[targ0_ix[file_ix]]
        oitarget[0].decep0 = olog.dec[targ0_ix[file_ix]]
        oitarget[0].equinox = olog.equinox[targ0_ix[file_ix]]
        oitarget[0].ra_err = .05
        oitarget[0].dec_err = .08
        oitarget[0].veltyp= 'UNKNOWN'
        oitarget[0].veldef= 'OPTICAL'

    ;Next the wavelength - only use one per save file.
        define_oiwavelength, oiwavelength_unit, nwave = 1
        oiwavelength = oiwavelength_unit
        oiwavelength.extver = 1
        oiwavelength.insname = olog.instrument[targ0_ix[file_ix]]
        oiwavelength.eff_wave = ptr_new(filter[0])
        oiwavelength.eff_band = ptr_new(filter[1])

	;Next the v^2 information
	    define_oivis2, oivis2_unit, nwave = 1
	    oivis2 = replicate(oivis2_unit, n_elements(v2_mean_boot))
	    oivis2(*).extver = 1
	    oivis2(*).date_obs = olog.date[targ0_ix[file_ix]]
	    oivis2(*).arrname = olog.mask
	    oivis2(*).insname = olog.instrument[targ0_ix[file_ix]]
	    oivis2(*).target_id = 0
	    oivis2(*).time = olog.utc[targ0_ix[file_ix]]
	    oivis2(*).mjd   = olog.jd[targ0_ix[file_ix]]
	    oivis2(*).int_time = olog.t_int[targ0_ix[file_ix]]
	    for i =  0, n_elements(v2_mean_boot)-1 do begin
	       oivis2[i].vis2data = ptr_new(v2_mean_boot[i])
	       oivis2[i].vis2err = ptr_new(v2_uncert_boot[i])
	       oivis2[i].flag = ptr_new(flag_val)
	    endfor
	    oivis2[*].ucoord    = u1[*,wavelength_ix]*filter[0]
	    oivis2[*].vcoord    = v1[*,wavelength_ix]*filter[0]
	    oivis2[*].sta_index  = bl2h_ix

    ;Last, the closure phase information
        define_oit3, oit3_unit, nwave = 1
        oit3 = replicate(oit3_unit, n_bispect)
        oit3(*).extver = 1
        oit3(*).date_obs = olog.date[targ0_ix[file_ix]]
        oit3(*).arrname = olog.mask
        oit3(*).insname = olog.instrument[targ0_ix[file_ix]]
        oit3(*).target_id = 0
        oit3(*).time =  olog.utc[targ0_ix[file_ix]]
        oit3(*).mjd   = olog.jd[targ0_ix[file_ix]]
        oit3(*).int_time = olog.t_int[targ0_ix[file_ix]]
        for i =  0,  n_bispect-1 do begin
           ; oit3[*].t3amp    = ptr_new(aquan.bs_amp[i])     ; ptr_new(0.0)
           ; oit3[*].t3amperr = ptr_new(aquan.bs_amp_err[i]) ;ptr_new(1.0)
           oit3[i].t3phi     = ptr_new(cp_mean_boot[i])
           oit3[i].t3phierr = ptr_new(cp_uncert_boot[i])
           oit3[i].flag     = ptr_new(flag_val)
         endfor
        oit3[*].u1coord    = reform(u1[bs2bl_ix[0, *]]*filter[0])
        oit3[*].v1coord    = reform(v1[bs2bl_ix[0, *]]*filter[0])
        oit3[*].u2coord    = reform(u1[bs2bl_ix[1, *]]*filter[0])
        oit3[*].v2coord    = reform(v1[bs2bl_ix[1, *]]*filter[0])
        
        for i = 0, n_bispect-1 do $
          oit3(i).sta_index  = [bl2h_ix[0, bs2bl_ix[0, i]], $
                                bl2h_ix[1, bs2bl_ix[0, i]],  $
                                bl2h_ix[1, bs2bl_ix[1, i]]]

        ;; Write it out!
        ; write_oidata,oif_name, oiarray, oitarget, oiwavelength, 0, 0, oit3
        write_oidata, oif_name, oiarray, oitarget, oiwavelength, 0, oivis2, oit3

        oif_to_merge = [oif_to_merge,oif_name]
    endfor


endfor

;; Merge the oifits
merge_name = repstr(oif_to_merge[0],'.oifits','_mrg.oifits')
merge_name_python = repstr(oif_to_merge[0],'.oifits','_python_mrg.oifits')

merge_oidata, outfile=merge_name,  infiles=oif_to_merge
;;ACC: fix all the problems caused by merge_oidata
if keyword_set(save_python) then fix_oifits,merge_name,merge_name_python,/python 
fix_oifits,merge_name,merge_name,python=0


end