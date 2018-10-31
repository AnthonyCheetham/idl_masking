#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 15:54:30 2017

@author: cheetham
"""

# Combine several SPHERE data cubes, so there are more frames per cube.
# Update: this actually isn't needed because it was the DCs fault (it binned some frames)

import astropy.io.fits as pf
import numpy as np
import glob,os

combine_n = 6

wdir = '/Users/cheetham/data/sphere_data/HIP_65426_SAM/IFS/'
data_dir = wdir+'Cubed_N2/'
save_dir = wdir+'Recubed_N2/'

os.chdir(wdir)

folders = glob.glob(wdir+'*/')

data = {}
headers = {}
files = np.sort(glob.glob(data_dir+'SPHER*.fits'))

# Loop through the datacubes
