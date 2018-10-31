# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 12:12:41 2016

@author: cheetham
"""
import naco_ispy

#wdir='/Users/cheetham/data/sphere_data/HIP_65426_SAM/IFS/Cubed_N2/'
wdir='/Users/cheetham/data/sphere_data/nuOct_SAM/2018-07-11/IRDIS/Raw/'
save_name = wdir+'../header_sci.txt'

header_keys={'tname':'HIERARCH ESO OBS TARG NAME',
      'datetime':'DATE-OBS',
#      'RA':'RA',
#      'Dec':'DEC',
      'RA': 'HIERARCH ESO INS4 DROT2 RA',
      'Dec': 'HIERARCH ESO INS4 DROT2 DEC',
#      'obstype':'HIERARCH ESO DPR TYPE',
      'obstype':'HIERARCH ESO DPR TYPE',
      'filter1':'HIERARCH ESO INS1 FILT NAME',
      'filter2':'HIERARCH ESO INS COMB IFLT',
      'dit':'HIERARCH ESO DET SEQ1 DIT',
      'ndit':'HIERARCH ESO DET NDIT',
      'NAXIS':'NAXIS',
      'ax1':'NAXIS1',
      'ax3':'NAXIS3',
      'nexpo':'HIERARCH ESO SEQ NEXPO',
      'rotstart':'HIERARCH ESO INS4 DROT2 BEGIN',
      'rotend':'HIERARCH ESO INS4 DROT2 END',
      'pastart':'HIERARCH ESO TEL PARANG START',
      'paend':'HIERARCH ESO TEL PARANG END',
      'alt':'ESO TEL ALT',
      'camera':'HIERARCH ESO SEQ ARM',
      'nd':'HIERARCH ESO INS4 FILT2 NAME'
      }

# Quick test
if False:
    import glob
    import astropy.io.fits as pf
    fs = glob.glob(wdir+'*.fits')
    f = fs[0]
    hdr = pf.getheader(f)
    for hk in header_keys.keys():
        name = header_keys[hk]
        if name not in hdr.keys():
            print 'Couldnt find:',name,'in header'
    

all_info=naco_ispy.header.make_header_file(wdir,prefix='SPHER',save_name=save_name,
                                          header_keys=header_keys)
