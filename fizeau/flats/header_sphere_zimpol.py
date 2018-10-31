# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 12:12:41 2016

@author: cheetham
"""
import naco_ispy

#wdir='/Users/cheetham/data/sphere_data/HIP_65426_SAM/IFS/Cubed_N2/'
wdir='/Users/cheetham/data/sphere_data/HD100546_SAM/ZIMPOL/Raw/'

header_keys={'tname':'HIERARCH ESO OBS TARG NAME',
      'datetime':'DATE-OBS',
      'RA':'RA',
      'Dec':'DEC',
#      'RA': 'HIERARCH ESO INS4 DROT2 RA',
#      'Dec': 'HIERARCH ESO INS4 DROT2 DEC',
#      'obstype':'HIERARCH ESO DPR TYPE',
      'obstype':'HIERARCH ESO DPR TYPE',
      'filter1':'HIERARCH ESO INS3 OPTI5 NAME',
      'filter2':'HIERARCH ESO INS3 OPTI6 NAME',
      'dit':'HIERARCH ESO DET DIT1',
      'ndit':'HIERARCH ESO DET NDIT',
      'NAXIS':'NAXIS',
      'ax1':'NAXIS1',
      'ax3':'NAXIS2',
      'nexpo':'HIERARCH ESO TPL NEXP',
      'rotstart':'HIERARCH ESO INS4 DROT2 BEGIN',
      'rotend':'HIERARCH ESO INS4 DROT2 END',
      'pastart':'HIERARCH ESO TEL PARANG START',
      'paend':'HIERARCH ESO TEL PARANG END',
      'alt':'ESO TEL ALT',
      'camera':'HIERARCH ESO SEQ ARM',
      'nd':'HIERARCH ESO INS3 OPTI6 NAME' # actual second filter
      }


all_info=naco_ispy.header.make_header_file(wdir,prefix='SPHER',save_name=wdir+'../header.txt',
                                           header_keys=header_keys)
