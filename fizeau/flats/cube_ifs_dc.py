# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 09:46:56 2016

@author: cheetham
"""

import astropy.io.fits as pf
import numpy as np
import glob,os

# Re-cube the sphere data

wdir = '/Users/cheetham/data/sphere_data/HD142527_SAM/IFS/'
#wdir = '/Users/cheetham/data/sphere_data/HIP_65426_SAM/IFS/'
data_dir = wdir+'Raw/'
save_dir = wdir+'Cubed/'
header_dir = wdir+'header_refs/'

os.chdir(wdir)

folders = glob.glob(wdir+'*/')

data = {}
headers = {}
files = np.sort(glob.glob(data_dir+'SPHER*.fits'))

# Load them and put them together
for f in files:
    loaded_cube,hdr = pf.getdata(f,header=True)
    # Find the right filename from the header
#    loaded_fname = hdr['ORIGFILE'] # This is unfortunately removed by the DC
    loaded_fname = hdr['ARCFILE']
    

    try:
        current_cube = data[loaded_fname]
        data[loaded_fname].append(loaded_cube)
    except:
        data[loaded_fname] = [loaded_cube]
        # Save the header the first time as well
        header_fname = header_dir+f.replace(data_dir,'')
        header_fname = header_fname.rsplit('_',1)[0]+'.fits'
        hdr = pf.getheader(header_fname)
        headers[loaded_fname] = hdr 
    

# Now save them all out
files = np.sort(data.keys())
for ix,f in enumerate(files):
    cube = data[f]
    cube = np.array(cube)
    cube = np.swapaxes(cube,0,1) # make it the right shape
    
    
    hdr = headers[f]
    
    # Fix the CRVAL header keywords
    hdr['CTYPE3'] = 'PIXEL   '
    
    # Make up a fake filename since the DC removes it from the headers
    fname = 'SPHERE_IRDIFS_IFS_OBS999_{0:04d}.fits'
    fname = fname.format(ix) 
    
    # And save out the file
    outname=save_dir+fname
    pf.writeto(outname,cube,header=hdr,clobber=True,
               output_verify='silentfix')