#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /data/user/aguilar/icetray/meta-projects/combo/V01-00-02/build

from time import gmtime, strftime
import numpy as np
import time
import random
import astropy
from astropy.time import Time as apTime
import sys, os
import math
from icecube import astro
import pickle
from optparse import OptionParser

import env

#source = 'source /data/ana/BSM/HT_Cascade/FinalAnalysisCode/env.sh'
#dump = '/usr/bin/python -c "import os,pickle;print pickle.dumps(os.environ)"'
#penv = os.popen('%s && %s' %(source,dump))
#env = pickle.loads(penv.read())
#os.environ = env

base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')

from fileList import nfiles, allFiles
from utils import psi_f

#######################
# get and define parameters
#######################

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-a", "--lecut", default='0.15',
                  dest="LECUT", help="Cut on LE BDT")
parser.add_option("-b", "--hecut", default='0.2',
                  dest="HECUT", help="Cut on HE BDT")
(options,args) = parser.parse_args()
        
LECut = float(options.LECUT)
HECut = float(options.HECUT)

Datatype = 'Data'

bins_vars={'energy_rec':np.logspace(1,5,48+1),'psi_rec':np.linspace(0.,np.pi,90+1)}

outfile_path = base_path+'PDFs/'+Datatype+'_scrambledFullSky/'

os.popen('mkdir -p '+outfile_path)
   
outfile = outfile_path+'/PDF_'+Datatype+'_ScrambleFullSky_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_2D'+'.pkl'
outfile_quad = outfile_path+'/PDF_'+Datatype+'_ScrambleFullSky_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_2D'+'_quad.pkl'

if os.path.isfile(outfile):
    print (' ... already exists')
    sys.exit(1)

livetime = {}
livetime['Burnsample'] = 1116125.821572
livetime['Data'] = 28272940. + 30674072. + 31511810.5 + 31150852. + 30059465.  

######## Arrays for all years ######
energy_reco = np.array([])
psi_reco = np.array([])
weight = np.array([])

for year in allFiles[Datatype].keys():
    print ('  year %s' %year)

    fileData = np.load(allFiles[Datatype][year], allow_pickle=True, encoding='latin1').item()

    print (' %i events in sample' %len(fileData['all']['energy_rec']))
    
    BDTScore_LE = np.array(fileData['all']['BDTScore_bb100'])
    BDTScore_HE = np.array(fileData['all']['BDTScore_ww300'])

    BDTmask = (BDTScore_LE>LECut)&(BDTScore_HE>HECut)
    tmp_energy_reco = np.array(fileData['all']['energy_rec'])[BDTmask]
    tmp_psi_reco = np.array(fileData['all']['psi_rec'])[BDTmask]
    tmp_zenith_reco = np.array(fileData['all']['zenith_rec'])[BDTmask]
    tmp_azimuth_reco = np.array(fileData['all']['azimuth_rec'])[BDTmask]       
    tmp_ra_reco = np.array(fileData['all']['RA_rec'])[BDTmask]
    tmp_dec_reco = np.array(fileData['all']['dec_rec'])[BDTmask]       

    print ('  %i events passed BDT cut' %len(tmp_energy_reco))

    RA_GC = 266./180.*np.pi
    dec_GC = -29.*np.pi/180

 #    argument = (np.cos(tmp_psi_reco) - np.cos(np.pi/2.-dec_GC)*np.cos(np.pi/2.-tmp_dec_reco))/\
 #       (np.sin(np.pi/2.-dec_GC)*np.sin(np.pi/2.-tmp_dec_reco))

#    print (tmp_delta_RA, np.arccos( argument ))
    
    #keep_argument_boundary =  np.where(argument<-1.)
    #keep_large_DeltaRA = np.where(tmp_delta_RA>np.pi/2.)
    
    #new_energy_reco = np.append(tmp_energy_reco[keep_argument_boundary], tmp_energy_reco[keep_large_DeltaRA])
    #new_dec_reco = np.append(tmp_dec_reco[keep_argument_boundary] , tmp_dec_reco[keep_large_DeltaRA])
    #new_weight = [2.]*len(new_dec_reco)
    
    #### us this if you don't want to cut in RA:
    new_energy_reco = tmp_energy_reco
    new_dec_reco = tmp_dec_reco
    new_weight = [1.]*len(new_dec_reco)
    
    #####
 
    new_RA_reco = np.random.random_sample((len(new_dec_reco),)) * 2*np.pi
    new_psi_reco = psi_f(new_RA_reco,new_dec_reco)     
    
    energy_reco = np.append(energy_reco, new_energy_reco)
    psi_reco = np.append(psi_reco, new_psi_reco)
    weight = np.append(weight, new_weight)


hist_pdf_signal = np.histogram2d(energy_reco,psi_reco,
                                 bins=(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = weight)
#hist_pdf_signal_quad = np.histogram2d(energy_reco,psi_reco,
#                                 bins=(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = weight**2)    
savefile = open(outfile,'wb')
pickle.dump(hist_pdf_signal, savefile)

#savefile_quad = open(outfile_quad,'wx')
#pickle.dump(hist_pdf_signal_quad, savefile_quad)
print (' ... done')
