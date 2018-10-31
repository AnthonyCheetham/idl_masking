#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 09:50:43 2018

@author: cheetham
"""

import glob
import numpy as np
import astropy.io.fits as pf

data_dir = '/Users/cheetham/data/sphere_data/NuOct_SAM/Night2//IFS/'
raw_dir = data_dir+'Raw/'
extracted_dir = data_dir+'Targ/'
output_dir = data_dir+'Cubed/'

extn='.fits'
dry_run = False

dc_mode = False # Use the DC filenames instead of the archive ones

#############

# Find the files in the raw directory
files = glob.glob(raw_dir+'SPHER*.fits')
files.sort()

last_base = ''
    
for file_ix,f in enumerate(files):
#    if file_ix < 32:
#        continue
    # work out how many frames there were in the original file to make the array
    # Also get the header
    hdr = pf.getheader(f)
#    nframes = hdr['NAXIS3']
    nframes = hdr['HIERARCH ESO DET NDIT']

    if dc_mode:
        fname = files[file_ix]
        base = fname.split('IFS_')[0]+extn
    else:
        base=files[file_ix]
        
    if base == last_base:
        continue
    
    last_base = base

    #loop through the frames and load them
    for ix in range(nframes):
        
        # Guess their filenames in the extracted directory
        if dc_mode:
#            fname = files[file_ix].replace(extn,'')[:-5]+'{0:05}'.format(ix)+extn
            fname = base.replace(extn,'')+'*_{0:05}'.format(ix)+extn
            fname = fname.replace(raw_dir,extracted_dir)
            fname = glob.glob(fname)[0]
        else:
#             base=strsplit(files[file_ix],data_dir,/regex,/extract)
            fname = base.replace(extn,'')+'_{0:05}'.format(ix)+extn
            fname = fname.replace(raw_dir,extracted_dir)

        #read in the file and add it to the big cube
        frame = pf.getdata(fname)
        # CUBE should be [Nwav,Nframes,291,291]
        if ix == 0:
            cube = np.zeros((frame.shape[0],nframes,frame.shape[1],frame.shape[2]))
        cube[:,ix,:,:] = frame

    outname = base.replace(raw_dir,output_dir)
    
    if not dry_run:
        pf.writeto(outname,cube,header=hdr,output_verify='silentfix',clobber=True)
    print 'File saved as: '+outname
#    raise Exception