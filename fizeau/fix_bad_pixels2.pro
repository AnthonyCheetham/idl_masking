;######## RELEASE CLIP HERE ######## 
function fix_bad_pixels2,image,headinfo,bad_pixels=bad_pixels,root_dir=root_dir
;;New version that uses MJI's magic remove_badpix.pro
;;Since we wont have the olog populated by the time we want to
;;remove bad pixels, we will need to make one ourselves.
;;
; bad_pixels is a list of bad pixels.
if (keyword_set(image) eq 0) then begin
   print,'function fix_bad_pixels2,im,headinfo,bad_pixels=bad_pixels,root_dir=root_dir'
   print,'    bad_pixels is a list of bad pixels.'
   return,-1
endif

headinfo.mask=inquire('mask',  headinfo)
mf_file=inquire('template',headinfo,ix=0)
sz=(size(image))[1] ;;wont work if it's not square...
restore,root_dir+'templates/'+mf_file
d = hole_diam/filter[0]*rad_pixel*sz

newimage=remove_badpix(image,-u,v,d,bad_pixels)

return,newimage
end



