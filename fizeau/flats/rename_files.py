#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 15:34:28 2017

@author: cheetham
"""

import os,glob,shutil
import astropy.io.fits as pf
import matplotlib.pyplot as plt

data_dir='/Users/cheetham/data/sphere_data/nuOct_SAM/2018-07-31/Raw/'
#data_dir = '/Users/cheetham/data/sphere_data/ScoCen_SAM_Mike/IRDIS/Raw/'
#data_dir='/Users/cheetham/data/sphere_data/NuOct_SAM/Night2/'

os.chdir(data_dir)


extn = '.fits'
#files = glob.glob('*OBS13*.fits')
files = glob.glob('*'+extn)

# plt.pause(5) # just make sure the file finishes first

files.sort()
nf=len(files)

for ix in range(nf):
    oldname = files[ix]

    head=pf.getheader(oldname)
    if 'ORIGFILE' in head:
        newname = head['ORIGFILE']
#        newname = head['ARCFILE']
        print newname
#        newname=sxpar(head,'ORIGFILE')
#        newname = 'SPHERE_IRDIFS_IRDIS_OBS999_{0:04}.fits'.format(ix+1)
#        raise Exception

#        cmd='mv '+oldname+' '+data_dir+newname
    #    spawn,cmd
        shutil.move(oldname,newname)
    else:
        print 'Couldnt find header keyword in',oldname
#        newname = 'SPHERE_IRDIFS_IFS_OBS999_{0:04}.fits'.format(ix+1)
        
        #        cmd='mv '+oldname+' '+data_dir+newname
    #    spawn,cmd
#        shutil.move(oldname,newname)