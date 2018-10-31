; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; %             program    bispect.pro 				       %
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;$Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $
;$Log: bispect.pro,v $
;Revision 1.18  2010/04/26 21:48:11  mireland
;Not sure...
;
;Revision 1.17  2010/04/08 21:06:57  dbernat
;oops.  left in a stop
;
;Revision 1.16  2010/04/08 20:50:49  dbernat
;switched ph_all for cvis_all
;
;Revision 1.15  2009/01/16 06:43:37  mireland
;A bunch of changes to again make this code work with Keck/nirc2 data. Some
;of these changes maybe need some more thought in the context of Cornell changes.
;
;Revision 1.14  2008/11/02 21:50:34  mireland
;A couple of bug fixes inspired by Woody.
;
;Revision 1.13  2008/10/24 21:23:59  dbernat
;BS_VAR calc took modsq( double( bs[j] ) ) but that double() takes the real part only and is unintended.
;
;Revision 1.12  2008/10/22 06:22:47  mireland
;calibrate_v2_cp had a conflict and isn't tested. But I think that there are
;basically bugfixes only here.
;
;Revision 1.11  2008/09/17 18:00:48  dbernat
;Baseline Phase information now saved (ph_all)
;
;Revision 1.10  2008/09/09 20:16:29  dbernat
;now saves all bs and v2 and added CP_VAR (var of CPS, not arctan(var of BS), see comments)
;
;Revision 1.9  2008/07/05 06:00:15  mireland
;Minor changes only.
;
;Revision 1.8  2008/05/21 23:12:46  mireland
;Haven't commited in a while - again unsure what the exact changes are...
;
;Revision 1.7  2007/12/06 01:47:39  mireland
;Added a destripe capability to nirc2. Unsure what else...
;
;Revision 1.6  2007/11/22 06:10:17  mireland
;I hope that these changes are all good - haven't commited in a while...
;
;Revision 1.5  2007/06/18 17:19:13  mireland
;Bugfixes...
;
;Revision 1.4  2007/06/15 00:31:28  mireland
;Still working on the covariance matrix stuff...
;
;Revision 1.3  2006/06/16 20:19:18  mireland
;Changed lots of small things. I think that most of the nirc2 changes were Peter's
;from a couple of months ago that he didn't commit. Now there is a n_blocks option
;for calc_bispect.pro (also needed in inquire) that first splits the data into
;n_blocks blocks of data before calculating variances (more realistic errors).
;
;Revision 1.2  2005/12/20 21:51:55  mireland
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and $Log: bispect.pro,v $
;Added $Id: bispect.pro,v 1.17 2010/04/08 21:06:57 dbernat Exp $ and Revision 1.18  2010/04/26 21:48:11  mireland
;Added $Id: bispect.pro,v 1.17 2010/04/08 21:06:57 dbernat Exp $ and Not sure...
;Added $Id: bispect.pro,v 1.17 2010/04/08 21:06:57 dbernat Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.17  2010/04/08 21:06:57  dbernat
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and oops.  left in a stop
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.16  2010/04/08 20:50:49  dbernat
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and switched ph_all for cvis_all
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.15  2009/01/16 06:43:37  mireland
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and A bunch of changes to again make this code work with Keck/nirc2 data. Some
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and of these changes maybe need some more thought in the context of Cornell changes.
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.14  2008/11/02 21:50:34  mireland
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and A couple of bug fixes inspired by Woody.
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.13  2008/10/24 21:23:59  dbernat
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and BS_VAR calc took modsq( double( bs[j] ) ) but that double() takes the real part only and is unintended.
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.12  2008/10/22 06:22:47  mireland
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and calibrate_v2_cp had a conflict and isn't tested. But I think that there are
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and basically bugfixes only here.
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.11  2008/09/17 18:00:48  dbernat
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Baseline Phase information now saved (ph_all)
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.10  2008/09/09 20:16:29  dbernat
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and now saves all bs and v2 and added CP_VAR (var of CPS, not arctan(var of BS), see comments)
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.9  2008/07/05 06:00:15  mireland
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Minor changes only.
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.8  2008/05/21 23:12:46  mireland
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Haven't commited in a while - again unsure what the exact changes are...
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.7  2007/12/06 01:47:39  mireland
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Added a destripe capability to nirc2. Unsure what else...
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.6  2007/11/22 06:10:17  mireland
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and I hope that these changes are all good - haven't commited in a while...
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.5  2007/06/18 17:19:13  mireland
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Bugfixes...
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.4  2007/06/15 00:31:28  mireland
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Still working on the covariance matrix stuff...
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Revision 1.3  2006/06/16 20:19:18  mireland
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and Changed lots of small things. I think that most of the nirc2 changes were Peter's
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and from a couple of months ago that he didn't commit. Now there is a n_blocks option
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and for calc_bispect.pro (also needed in inquire) that first splits the data into
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and n_blocks blocks of data before calculating variances (more realistic errors).
;Added $Id: bispect.pro,v 1.18 2010/04/26 21:48:11 mireland Exp $ and to important files, and added the g18 mask
;for PHARO's inquire.
;

; NB: This program has not had the bispectral bias option properly
; tested yet.
;
;Keywords:
; coherent_bs:      NOT IMPLEMENTED (assumed 1)!
; dark_ps:          Input a dark ps to remove correlated noise (!!! not
;                   tested yet)
; subtract_bs_bias :Set 1 so subtract the bias in the bispectrum
; fluxes:           This will store the fluxes for all frames.
; avar:             The variance in v2 in a single-frame of data due
;                   to atmospheric effects
;
;
;
;phs_v2corr:        A correction to apply to each baseline's V^2 based
;                   on phase slope differences between holes. 
;predictor:         For testing phs_v2corr only

; ACC: this program is edited to calculate V^2 then average, rather than average the powers then divide by flux to calculate V^2.

pro bispect_gpi_edit, ft_arr, mf_pvct, mf_gvct, mf_ix, mf_rmat, mf_imat, v2, v2_cov, bs,bs_var, bs_cov, bs_v2_cov, $
                 bl2h_ix, bs2bl_ix, bl2bs_ix, bscov2bs_ix,cp_var, bs_all, v2_all, cvis_all, closing_tri_pix=closing_tri_pix, $
                 mfc_pvct=mfc_pvct, mfc_gvct=mfc_gvct, fluxes=fluxes,avar=avar, err_avar=err_avar,$
                 coherent_bs=coherent_bs, dark_ps=dark_ps, subtract_bs_bias=subtract_bs_bias, sig2_arr=sig2_arr,$
                 v2_arr=v2_arr,  phs_v2corr = phs_v2corr,  hole_phs = hole_phs,  hole_err_phs = hole_err_phs, $
                 hole_piston = hole_piston,  imsize = imsize,  predictor = predictor,  n_blocks = n_blocks_in,  cp_cov = cp_cov

if (keyword_set(coherent_bs)) then coherent_bs = 1
if keyword_set(closing_tri_pix) then using_pgt_mf=1 else using_pgt_mf=0

smallfloat =  1e-16             ;A very small number that you can still square and not get floating underflows.
;ft_arr is an n_specks array of Fourier transforms...
n_holes = max(bl2h_ix)+1.0
n_baselines = (size(bl2bs_ix))[1]
n_bispect = (size(bs2bl_ix))[2]
n_cov = (size(bscov2bs_ix))[2]
i_ps=size(ft_arr)  
n_ps=i_ps(3)                    ;n_ps = number of speckles.
if (keyword_set(n_blocks_in) eq 0) then n_blocks = n_ps else n_blocks = n_blocks_in
if (n_blocks gt n_ps) then begin
   print,  "WARNING: n_blocks greater than n_ps: setting n_blocks to n_ps."
   n_blocks = n_ps
endif
dim1=i_ps(1)                                  ; # ofpixels in x axix
dim2=i_ps(2)                                  ; # of pixels in yaxis
if i_ps[0] eq 4 then dim3=i_ps(4) else dim3=1 ; # of wavelengths

aveps=fltarr(dim1,dim2,dim3)
avedps=fltarr(dim1,dim2,dim3)   ;Do dark ps subtraction - for noise estimation...
;A one speckle frame option is unweildy and defeats the purpose of
;calculating the covariance matrix etc...
if (i_ps(0) eq 2) then begin
   print, 'This routine (bispect) needs more than 1 speckle frame...'
   stop 
endif

;___________________________________________________________
;INITIALIZE ARRAYS

cvis=dcomplexarr(n_baselines)
rvis = dblarr(n_baselines)
ivis = dblarr(n_baselines)
dark_v2 = dblarr(n_baselines)
v2_arr=dblarr(n_ps,n_baselines,dim3)
ph_arr=dblarr(n_ps,n_baselines,dim3)
phs_arr=dblarr(2,n_ps,n_baselines,dim3)
phserr_arr=dblarr(2,n_ps,n_baselines,dim3)
biasmn=dblarr(n_baselines,dim3)
bs_arr = dcomplexarr(n_ps,n_bispect,dim3)
bs2_arr = dcomplexarr(n_ps,n_bispect,dim3)
cvis_arr = dcomplexarr( n_ps, n_baselines,dim3)
fluxes =  dblarr(n_ps,dim3)



; ###############################
; PART I - Fill up arrays
; ###############################
for wav=0,dim3-1 do begin
   for i=0,n_ps-1 do begin   	
; Build up PowerSpectrum
      ft_frame = ft_arr[*,*,i,wav]
      ps=modsq(ft_frame)
      
; Subtract off dark pspec if option is set:
      if (keyword_set(dark_ps) and (size(dark_ps))[0] eq 3) then dps = dark_ps[*,*,i,wav] $
      else if (keyword_set(dark_ps) and (size(dark_ps))[0] eq 2) then dps = dark_ps $
      else dps = dblarr(dim1,dim2)
      avedps[*,*,wav]+=dps
      aveps[*,*,wav]+=ps
      
      for j=0,n_baselines-1 do begin
         dark_v2[j] = total(mf_gvct[mf_ix[0,j,wav]:mf_ix[1,j,wav]]^2 * dps[mf_pvct[mf_ix[0,j,wav]:mf_ix[1,j,wav]]])
         cvis[j] = total(mf_gvct[mf_ix[0,j,wav]:mf_ix[1,j,wav]] * ft_frame[mf_pvct[mf_ix[0,j,wav]:mf_ix[1,j,wav]]])
         pix=mf_pvct[mf_ix[0,j,wav]:mf_ix[1,j,wav]]
         ftf1=shift(ft_frame,1,0)
         ftf2=shift(ft_frame,-1,0) 
         dummy=total(ft_frame[pix]*conj(ftf1[pix])+conj(ft_frame[pix])*ftf2[pix])
         phs_arr[0,i,j,wav]=atan(dummy, /phase)
         phserr_arr[0,i,j,wav]=1./abs(dummy) ;NB only a relative error
         ftf1=shift(ft_frame,0,1)
         ftf2=shift(ft_frame,0,-1) 
         dummy=total(ft_frame[pix]*conj(ftf1[pix])+conj(ft_frame[pix])*ftf2[pix])
         phs_arr[1,i,j,wav]=atan(dummy, /phase)
         phserr_arr[1,i,j,wav]=1./abs(dummy) ;NB only a relative error
      endfor
      
      ;;Correct for overlapping baselines
      rvis = float(cvis)
      ivis = imaginary(cvis)
      rvis = mf_rmat # rvis
      ivis = mf_imat # ivis
      cvis = dcomplex(rvis,ivis)
      cvis_arr[i,*,wav] = cvis
      
      ;;Save phases
      ph_arr[i,*,wav]=atan(cvis, /phase)
      
      ;;Calculate square visibilities (correct for bias later)
      v2_arr[i,*,wav] = modsq(cvis) - dark_v2
      
      ;;Calculate Bispectrum
      if(using_pgt_mf ne 1) then $
         bs_arr[i,*,wav] = cvis[bs2bl_ix[0,*]]*cvis[bs2bl_ix[1,*]]*conj(cvis[bs2bl_ix[2,*]]) $
      else begin
         ;print,"Anthony hasn't finished this. It doesn't seem like it's working."
         ;stop
         ;; !! This is the April2013 code update which implements  
         ;; the pixel-triangle-loops to explicitly populate the bispectrum
         ;; Firstly calculate up the complex spectrum multiplied by the match filter (mfilter_spec)
         mfilter_spec=fltarr(dim1,dim2)
         mfilter_spec[mfc_pvct[mf_ix[0,0,wav]:mf_ix[1,n_baselines-1,wav]]]=mfc_gvct[mf_ix[0,0,wav]:mf_ix[1,n_baselines-1,wav]]
         mfilter_spec+=shift(reverse(reverse(mfilter_spec),2),1,1)
         mfilter_spec = mfilter_spec * shift(ft_frame,dim1/2,dim2/2)
     
         base_origin=closing_tri_pix[0,0]
         n_closing_tri= (size(closing_tri_pix))[2]
         for this_bs=0,n_bispect-1 do begin
            this_tri=bs2bl_ix[*,this_bs]
            tri_splodge_origin = mfc_pvct[mf_ix[0,this_tri,wav]] ; this gives the origin pixel of the three splodge sampling circles
            splodge_shift= tri_splodge_origin - base_origin  ; this gives the offset between the closing triangle vector and the actual sampling needed
            this_trisampling = closing_tri_pix + transpose(splodge_shift##(replicate(1,n_closing_tri))) ; apply offset to entire big trisampling
            
            ;;ACC fixed a bug here! If using 2 "positive" splodges and one "negative" splodge then the triangle all close. However, if using 3 "positive" splodges, the third one must be flipped to ensure all triangles close. These lines move to the negative splodge rather than try to flip them. Also "conj" was removed from the bs_arr[i,this_bs] line to reflect this change
            spl_offset=array_indices(ft_frame,(splodge_shift+(long(dim1/2)*(dim1+1))))-dim1/2
            this_trisampling[2,*]=this_trisampling[2,*]-2*(spl_offset[0,2]+spl_offset[1,2]*long(dim1))

            ;; OK now we can just compute the bispectrum, triangle by triangle, by using the pre-prepared big-block-of-closing-pixels:
            bs_arr[i,this_bs,wav] = total(mfilter_spec[this_trisampling[0,*]] * mfilter_spec[this_trisampling[1,*]] * mfilter_spec[this_trisampling[2,*]])
         endfor
      endelse
      
      ;;Calculate flux
      fluxes[i,wav] = abs(ft_arr[0,0,i,wav]) - sqrt(dps[0, 0])
   endfor
endfor

ps = aveps/n_ps
dps = avedps/n_ps
;Find the bias. First find the region where there is signal...
signal = dblarr(dim1,dim2,dim3)
for wav=0,dim3-1 do begin
   for i=0,n_baselines-1 do begin
      dummy=fltarr(dim1,dim2)
      pix=mf_pvct[mf_ix[0,i,wav]:mf_ix[1,i,wav]]
      dummy[pix]=1.0
      signal(*,*,wav) += dummy
   endfor
endfor

for wav=0,dim3-1 do signal[*,*,wav] = signal[*,*,wav] + shift(rotate(signal[*,*,wav],2),1,1)
d = dist(dim1)
for wav=0,dim3-1 do signal[where(d lt dim1/16),wav] = 1
;Now, from where there is no signal, find the median
for i = 0,3 do begin
   w = where(signal eq 0)
   bias = median(ps[w])
   dark_bias = median(dps[w])
   signal[where(ps gt 3.0*bias)] = 1.0
endfor
;;w = where(signal ne 0)          ;ie where there is signal.

;Now we know where the signal is, we can find the 'real' bias that
;includes the correlation between neighbouring terms in each
;ft_frame...
autocor_noise = dblarr(dim1, dim2,wav)
for wav=0,dim3-1 do begin
   w = where(signal[*,*,wav] ne 0) ;ie where there is signal.
   for i=0,n_ps-1 do begin
      ft_frame = ft_arr[*,*,i,wav]
      ft_frame[w] = 0.0
      autocor_noise[*,*,wav]+= float(fft(modsq(fft(ft_frame,-1)),1))
   endfor
endfor
autocor_noise = autocor_noise/n_ps
for wav=0,dim3-1 do begin
   for j = 0,n_baselines-1 do begin
      mf = dblarr(dim1,dim2)
      mf[mf_pvct[mf_ix[0,j,wav]:mf_ix[1,j,wav]]] = mf_gvct[mf_ix[0,j,wav]:mf_ix[1,j,wav]]
      autocor_mf = float(fft(modsq(fft(mf,-1)),1))
      biasmn[j,wav] = total(autocor_mf*autocor_noise[*,*,wav])*bias/autocor_noise[0,0,wav]*dim1*dim2      
      v2_arr[*,j,wav] = v2_arr[*,j,wav] - biasmn[j,wav]
   endfor
endfor

;Problem: we've subtracted the dark noise as well here (already
;subtracted as dark_v2 above), so we must add it back. It might look like we've subtracted then
;added a similar thing (which is true) but it's not the same, as the
;structure of dark_ps goes into the different values of dark_v2, but
;dark_bias here is constant for all baselines.
for wav=0,dim3-1 do for j = 0,n_baselines-1 do v2_arr[*,j,wav] = v2_arr[*,j,wav] + total(mf_gvct[mf_ix[0,j,wav]:mf_ix[1,j,wav]]^2)*dark_bias

ft_arr = 0 ;Done with this - free memory. We'll need it for bs_cov! ???
dark_ps=0  ; We have finished with this big data cube.

;; ACC: Now normalise by flux. This ensures that any flux changes don't affect the visibilities or their scatter.
bs_all=complexarr(n_blocks,n_bispect,dim3) ;; this used to say complex...
v2_all=dblarr(n_blocks,n_baselines,dim3);; this used to say double...
cvis_all=complexarr(n_blocks,n_baselines,dim3)
for wav=0,dim3-1 do begin
   for fr=0,n_ps-1 do begin ;;ACC added this loop to test flux
      bs_all[fr,*,wav] = complex( bs_arr[fr,*,wav]/mean( fluxes[fr,wav]^3 ) * n_holes^3 )
      v2_all[fr,*,wav] = v2_arr[fr,*,wav] / mean( fluxes[fr,wav]^2 ) * n_holes^2
      cvis_all[fr,*,wav] = cvis_arr[fr,*,wav] / mean( fluxes[fr,wav] ) * n_holes
   endfor
endfor

;;ACC: this is easier than changing all the following equations
v2_arr=v2_all
cvis_arr=cvis_all
bs_arr=bs_all

; #########################################################
; PART II - Turn Arrays into means and covariance matrices
; #########################################################

v2 = dblarr(n_baselines,dim3)
v2_cov = dblarr(n_baselines,n_baselines,dim3)
bs = dcomplexarr(n_bispect,dim3)
bs_var = dblarr(2,n_bispect,dim3)
bs_cov = dblarr(2,n_cov,dim3)
bs_v2_cov = dblarr(n_baselines,n_holes-2,dim3)

print, 'Calculating mean V^2 and variance...'
v2 = reform(total(v2_arr,1)/n_ps)
v2diff = dblarr(n_blocks, n_baselines,dim3)
for wav=0,dim3-1 do for j = 0,n_baselines-1 do for k = 0, n_blocks-1 do $
   v2diff[k,j,wav] = mean(v2_arr[k*n_ps/n_blocks:(k+1)*n_ps/n_blocks-1,j,wav]) - v2[j,wav]
for wav=0,dim3-1 do for j = 0,n_baselines - 1 do for k = 0,n_baselines-1 do $ 
   v2_cov[j,k,wav] = total(v2diff[*,j,wav] * v2diff[*,k,wav])/(n_blocks - 1.0)/float(n_blocks)
;Now, for the case where we have coherent integration over splodges,
;it is easy to calculate the bias in the variance. (note that the
;variance of the bias is simply the square of it's mean)
x = indgen(n_baselines)
temp_v2_cov=dblarr(n_baselines,dim3)
for wav=0,dim3-1 do temp_v2_cov[*,wav]=v2_cov[x,x,fltarr(n_baselines)+wav]
avar = temp_v2_cov*n_ps - biasmn^2*(1 + 2.0*v2/biasmn)
err_avar = sqrt(2.0/n_ps*double(temp_v2_cov)^2*n_ps^2 + 4.0*temp_v2_cov*biasmn^2) ;Assumes no error in biasmn...
print, 'Calculating mean bispectrum and variance...'
bs = reform(total(bs_arr,1)/n_ps)
;Bispectral bias subtraction. This assumes that the matched filter has been
;correctly normalised...
if (keyword_set(subtract_bs_bias)) then begin
   stop ;;ACC: this probably doesn't work now that we use V^2 and not power...
   bs_bias = v2[bs2bl_ix[0,*],*] + v2[bs2bl_ix[1,*],*] + v2[bs2bl_ix[2,*],*] + mean(fluxes)
   bs = bs - bs_bias
   print, 'Maximum bispectrum bias (as a fraction of bispectral amplitude): ', max(bs_bias/abs(bs))
endif

temp =  dcomplexarr(n_blocks)

for wav=0,dim3-1 do begin
   for j=0,n_bispect-1 do begin
      ;;temp is the complex difference from the mean, shifted so that the real
      ;;axis corresponds to amplitude and the imaginary axis phase.
      temp2 = dcomplex(bs_arr[*,j,wav]-bs[j,wav])*conj(bs[j,wav])
      for k = 0, n_blocks-1 do temp[k] = mean(temp2[k*n_ps/n_blocks:(k+1)*n_ps/n_blocks-1])
      bs_var[0,j,wav] = total(real_part(temp)^2)/double(n_blocks)/(n_blocks - 1.0)/modsq( bs[j,wav] )
      bs_var[1,j,wav] = total(imaginary(temp)^2)/double(n_blocks)/(n_blocks - 1.0)/modsq( bs[j,wav] )
      
   endfor
endfor

print, 'Calculating covariance between power and bispectral amplitude...'
; This complicated thing below calculates the dot product between the
; bispectrum point and its error term ie (x . del_x)/|x| and
; multiplies this by the power error term. Note that this is not the
; same as using absolute value, and that this sum should be zero where
; |bs| is zero within errors.
print, 'Calculating the bispectral covariances...'
for wav=0,dim3-1 do begin
   for j=0,n_baselines-1 do for k = 0,n_holes-3 do $
      bs_v2_cov[j,k,wav] = total(real_part((bs_arr[*,bl2bs_ix[j,k],wav]-bs[bl2bs_ix[j,k],wav])*conj(bs[bl2bs_ix[j,k],wav])) $
                                 *double(v2_arr[*,j,wav] - v2[j,wav]))/abs(bs[bl2bs_ix[j,k],wav])/(n_ps - 1.0)/double(n_ps)
   for j=long(0),n_cov-1 do begin
      temp1 = dcomplex( (bs_arr[*,bscov2bs_ix[0,j],wav] - bs[bscov2bs_ix[0,j],wav])*conj(bs[bscov2bs_ix[0,j],wav]) )
      temp2 = dcomplex( (bs_arr[*,bscov2bs_ix[1,j],wav] - bs[bscov2bs_ix[1,j],wav])*conj(bs[bscov2bs_ix[1,j],wav]) )
      denom = double( abs(bs[bscov2bs_ix[0,j],wav])*abs(bs[bscov2bs_ix[1,j],wav])*(n_ps - 1.0)*double(n_ps) )
      bs_cov[0,j,wav] = total(real_part(temp1)*real_part(temp2))/denom
      bs_cov[1,j,wav] = total(imaginary(temp1)*imaginary(temp2))/denom
   endfor
endfor

if (n_holes le 21 and arg_present(cp_cov)) then begin ;This already takes 2.6 times the memory and time for 21 holes over 18...
   cp_cov = dblarr(n_bispect, n_bispect,dim3)
   for wav=0,dim3-1 do begin
      for i = 0, n_bispect-1 do for j = 0, n_bispect-1 do begin
         temp1 = dcomplex( (bs_arr[*,i,wav] - bs[i,wav])*conj(bs[i,wav]) )
         temp2 = dcomplex( (bs_arr[*,j,wav] - bs[j,wav])*conj(bs[j,wav]) )
         denom = double(modsq(bs[i,wav]))*double(modsq(bs[j,wav]))*(n_ps - 1.0)*double(n_ps) 
         cp_cov[i, j,wav] = total(imaginary(temp1)*imaginary(temp2))/denom
      endfor
   endfor
endif else cp_cov = -1 

;Now normalise all return variables...

;;ACC: This whole section doesn't work if the flux changes slightly between frames...
;
;for wav=0,dim3-1 do begin
;   v2[*,wav] = v2[*,wav]/mean(fluxes[*,wav]^2)*n_holes^2
;   bs[*,wav] = complex(bs[*,wav]/mean(fluxes[*,wav]^3)*n_holes^3)
;   v2_cov[,*,wav] = double(v2_cov[*,*,wav]/mean(fluxes[*,wav]^4))*n_holes^4
;   ix =  indgen(n_baselines)
;   v2_cov[ix, ix,wav] = v2_cov[ix, ix,wav] > smallfloat
;      avar[fr,wav] = avar[fr,wav]/mean(fluxes[fr,wav]^4)*n_holes^4 > smallfloat
;      err_avar[*,wav] = err_avar[fr,wav]/mean(fluxes[fr,wav]^4)*n_holes^4 > smallfloat
;      fluxes[fr,wav] = fluxes[fr,wav]/10000.0 ;Prevent floating overflows...
;      bs_v2_cov[*,*,wav] = real_part(bs_v2_cov[*,*,wav]/mean(fluxes[*,wav]^5)*n_holes^5)/10000.0^5
;      bs_cov[*,*,wav] = float(bs_cov[*,*,wav]/mean(fluxes[*,wav]^6)*n_holes^6)/10000.0^6 ;; this won't work properly if the flux changes between frames...
;      bs_var[*,*,wav] = (float(bs_var[*,*,wav]/mean(fluxes[*,wav]^6)*n_holes^6)/10000.0^6) > smallfloat
;endfor

;Finally, convert baseline variables to hole variables...
 ;1) In the MAPPIT-style, we define the relationship between hole phases 
 ;(or phase slopes) and baseline phases (or phase slopes) by the
 ;use of a matrix, fitmat.
 common phasestuff, fitmat, ph_mn, ph_err
 fitmat = dblarr(n_holes, n_baselines+1)
 for j=0,n_baselines-1 do fitmat[bl2h_ix[0,j],j] = 1.0
 for j=0,n_baselines-1 do fitmat[bl2h_ix[1,j],j] = -1.0
 fitmat[0,n_baselines]=1.0 
 ;Firstly, fit to the phases by doing a weighted least-squares fit
 ;to baseline phasors.
 phasors=exp(complex(0,ph_arr))
 ph_mn = atan(total(phasors,1), /phase)
 ph_err = replicate(1.0,[n_baselines,dim3])
 for wav=0,dim3-1 do for j = 0,n_baselines-1 do ph_err[j,wav] = stdev( ((ph_arr[*,j,wav]-ph_mn[j,wav]+3*!pi) mod (2*!pi)) - !pi)
 ph_err = ph_err/sqrt(n_ps)
 
;;this next step fits the hole pistons. These should be the same for each wavelength?
 ;;this step can be done better:
 hole_piston=fltarr(n_holes)
 hole_piston = amoeba(1e-3, p0=replicate(0.0,n_holes), $
  function_value=fval, function_name='phase_chi2_gpi', scale=1.4)
 if (hole_piston[0] ne -1) then $
  print, 'Phase Chi^2: ', phase_chi2_gpi(hole_piston)/(n_baselines-n_holes+1) $
 else print,  'Error calculating hole pistons...'

;2) fit to the phase slopes using weighted linear regression.
 ;Normalisation:  hole_phs was in radians per Fourier pixel.
 ;Convert to phase slopes in pixels.
 phs_arr = phs_arr/2/!pi*dim1
 hole_phs = dblarr(2,n_ps,n_holes,dim3)
 hole_err_phs=dblarr(2,n_ps,n_holes,dim3)
 for j=0,n_baselines-1 do fitmat[bl2h_ix[1,j],j] = 1.0
 fitmat=fitmat/2.0
 fitmat=fitmat[*,0:n_baselines-1]
 err2=v2_arr ;ie a matrix the same size
 err2_bias=err2
 for wav=0,dim3-1 do begin
    for j = 0,n_ps-1 do begin
       hole_phs[0,j,*,wav]=$
          regress_noc(fitmat, reform(phs_arr[0,j,*,wav]),reform(phserr_arr[0,j,*,wav]), cov,YFIT, MSE, var_yfit)
       dummy=cov2cor(cov,sig=sig)
       hole_err_phs[0,j,*,wav]=sig*sqrt(MSE)
       hole_phs[1,j,*,wav]=$
          regress_noc(fitmat, reform(phs_arr[1,j,*,wav]),reform(phserr_arr[1,j,*,wav]), cov,YFIT, MSE, var_yfit)
       dummy=cov2cor(cov,sig=sig)
       hole_err_phs[1,j,*,wav]=sig*sqrt(MSE)
       err2[j,*,wav] = (hole_phs[0,j,bl2h_ix[0,*],wav]-hole_phs[0,j,bl2h_ix[1,*],wav])^2+$
                       (hole_phs[1,j,bl2h_ix[0,*],wav]-hole_phs[1,j,bl2h_ix[1,*],wav])^2
       err2_bias[j,*,wav] = (hole_err_phs[0,j,bl2h_ix[0,*],wav]-hole_err_phs[0,j,bl2h_ix[1,*],wav])^2+$
                            (hole_err_phs[1,j,bl2h_ix[0,*],wav]-hole_err_phs[1,j,bl2h_ix[1,*],wav])^2
    endfor
 endfor
 
 hole_mnphs=total(hole_phs,2)/n_ps 
 mnerr2=(hole_mnphs[0,bl2h_ix[0,*],*]-hole_mnphs[0,bl2h_ix[1,*],*])^2+$
        (hole_mnphs[1,bl2h_ix[0,*],*]-hole_mnphs[1,bl2h_ix[1,*],*])^2
 phs_v2corr=dblarr(n_baselines,dim3)
 predictor = v2_arr
 for wav=0,dim3-01 do for j=0,n_baselines-1 do predictor[*,j,wav] = err2[*,j,wav] -mean(err2_bias[*,j,wav])
 ;imsize is \lambda/hole_diameter in pixels. A factor of 3.0 was only
 ;roughly correct based on simulations.
 ;2.5 seems to be better based on real data. 
 ;NB there is no window size adjustment here. 
for wav=0,dim3-1 do for j=0,n_baselines-1 do phs_v2corr[j,wav] = mean(exp( -2.5*predictor[*,j,wav]/imsize^2 ))

end
  
