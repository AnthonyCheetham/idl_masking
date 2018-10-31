;;this just

prefix='SPHER'
calib_dir='~/data/sphere_data/calib/'

window,0
window,1
print,  'Press enter to end...'
done_files=''
!except=0
while (1) do begin
    ;;poll the current directory for files
    spawn,'ls '+prefix+'*.fits*',new_files,errs

    nf=n_elements(new_files)
    for ix=0,nf-1 do begin
        w=where(done_files eq new_files[ix])
        if w eq -1 then begin
            print,'new file found!'
            wait,0.5
            done_files=[done_files,new_files[ix]]

            ;;open the file
            im=readfits(new_files[ix],head,/silent)

            ;;run freud to get header info
            olog=make_olog(1,[0.01],0,0)
            dummy=freud(head,olog=olog)
            camera=olog.optic_cfg
            if camera eq 'IRDIS' then begin
                print,'  IRDIS frame'
                ;irdis2ps,im,caldir

                splitname=strsplit(new_files[ix],'.fit',/extract,/regex)
                extn='.fit'+splitname[1]
                prefix=(strmid(splitname[0].reverse(),4)).reverse() ;;use crazy idl object oriented string methods
                frames=(strmid(splitname[0].reverse(),0,4)).reverse() ;;use crazy idl object oriented string methods
                frames=[fix(frames)]
                skies=0
                tsize=[-0.01]
                show_images=0
                discard_sigma=[-1,-1,-1,-1,-1]
                jobname='quicklook'
                merge_irdis=1

                ;;now pick a flat
                flatpath=calib_dir+'flat_'+olog.filter+'_IRDIS_2048.idlvar'
                wset,1

                qbe_sphere,'',jobname,flatpath,frames,skies,tsize,show_images=show_images,prefix=prefix,extn=extn,cutsize=256,discard_sigma=discard_sigma,merge_irdis=merge_irdis
                wset,0
                calc_bispect,"cubeinfo"+jobname+".idlvar" , root_dir="~/code/masking/"
            endif else if camera eq 'IFS' then begin
                print,'  IFS frame'
            endif else if camera eq 'ZIMPOL' then begin
                print,'  ZIMPOL frame'
            endif else begin
                print,'  camera not recognised:',camera
                print,'  skipping file'
            endelse
            print,'Processing complete. Waiting for new file.'
        endif
    endfor

    if (get_kbrd(0) ne '') then goto,  finish
    wait,0.5
endwhile
finish:

end