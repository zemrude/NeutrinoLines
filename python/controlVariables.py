#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /data/user/sbaur/metaprojects/icerec/build/

from __future__ import division

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

base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')

from fileList import nfiles, allFiles
from IceCube_sim import genie_correction_factor, oversample
from fluxCalculation import phiDM_ann, phiDM_dec, psi_f, osc_Atm_flux_weight, neutrinoOscillator

#######################
# get and define parameters
#######################

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-a", "--lecut", default='0.2',
                  dest="LECUT", help="Cut on LE BDT")
parser.add_option("-b", "--hecut", default='0.2',
                  dest="HECUT", help="Cut on HE BDT")
parser.add_option("-c", "--psicut", default='90',
                  dest="PSICUT", help="Cut on Psi")
(options,args) = parser.parse_args()
        
LECut = float(options.LECUT)
HECut = float(options.HECUT)
PsiCut = float(options.PSICUT)/180*np.pi
    
astro_gamma = -2.89
    
bins_vars={
    'energy_rec':np.logspace(1,5,24+1),
    'psi_rec':np.linspace(0.,np.pi,18+1),
    'azimuth_rec':np.linspace(-2*np.pi,2*np.pi,72+1),
    'zenith_rec':np.linspace(0.,np.pi,18+1),

}

burnsampleTime = 1116125.821572

outfile = base_path+'controlDistributions/Histograms_LEBDT'+options.LECUT+'_HEBDT'+options.HECUT+'_PSI'+options.PSICUT+'degrees.pkl'

if os.path.isfile(outfile):
    print ' ... already exists'
    sys.exit(1)

histograms = {}
unbinned_data = {}

for MCtype in ['corsika','astro','atm']:

    print 'start with',MCtype
    
    weight = np.array([])
    energy_reco = np.array([])
    psi_reco = np.array([])
    azimuth_reco = np.array([])
    zenith_reco = np.array([])
    
    if 'corsika' == MCtype:
        
        fileData = np.load(allFiles['MC']['corsika']).item()

        weight = np.append(weight,fileData['all weight'])
        energy_reco = np.append(energy_reco,fileData['all']['energy_rec'])
        psi_reco = np.append(psi_reco,fileData['all']['psi_rec'])
        
        azimuth_reco = np.append(azimuth_reco,fileData['all']['azimuth_rec'])
        zenith_reco = np.append(zenith_reco,fileData['all']['zenith_rec'])

        BDTScore_LE = np.array(fileData['all']['BDTScore_bb100'])
        BDTScore_HE = np.array(fileData['all']['BDTScore_ww300'])

        weight[BDTScore_LE<LECut] = 0.
        weight[BDTScore_HE<HECut] = 0.

        weight[psi_reco<PsiCut] = 0.

        
    else: 

        for fileType in allFiles['MC']['nominal'].keys():

            fileData = np.load(allFiles['MC']['nominal'][fileType]).item()

            tmp_weight = np.array(fileData['all weight'])[:,0]
            tmp_nu_type = np.array(fileData['all weight'])[:,1]

            tmp_energy_true = np.array(fileData['all']['energy_true'])
            tmp_psi_true = np.array(fileData['all']['psi_true'])

            tmp_energy_reco = np.array(fileData['all']['energy_rec'])
            tmp_psi_reco = np.array(fileData['all']['psi_rec'])

            tmp_zenith_true = np.array(fileData['all']['zenith_true'])

            tmp_zenith_reco = np.array(fileData['all']['zenith_rec'])
            tmp_azimuth_reco = np.array(fileData['all']['azimuth_rec'])
            
            BDTScore_LE = np.array(fileData['all']['BDTScore_bb100'])
            BDTScore_HE = np.array(fileData['all']['BDTScore_ww300'])

            tmp_weight[BDTScore_LE<LECut] = 0.
            tmp_weight[BDTScore_HE<HECut] = 0.
            tmp_weight[tmp_psi_reco<PsiCut] = 0.
            
            if 'nugen' in fileType:
                tmp_weight[tmp_energy_true<190] = 0.
                tmp_weight[tmp_energy_true<195] = tmp_weight[tmp_energy_true<195]*(0.2*tmp_energy_true[tmp_energy_true<195]-38.)
                tmp_weight = tmp_weight/0.5

            elif 'genie' in fileType:
                tmp_weight[tmp_energy_true>195] = 0.
                tmp_weight[tmp_energy_true>190] = tmp_weight[tmp_energy_true>190]*(-0.2*tmp_energy_true[tmp_energy_true>190]+39.)
                tmp_weight = tmp_weight*genie_correction_factor(tmp_nu_type)

            tmp_weight = tmp_weight/nfiles['nominal'][fileType]

            mask = tmp_weight>0.
            
            if 'astro' == MCtype:
                tmp_flux = 6.45 * 1e-18 * pow(tmp_energy_true[mask]/1e5, astro_gamma) / 6.
            else:
                osc = neutrinoOscillator(atmH=20)
                tmp_flux = osc_Atm_flux_weight(tmp_energy_true[mask],tmp_zenith_true[mask],tmp_nu_type[mask],osc)
            
            weight = np.append(weight,tmp_weight[mask]*tmp_flux)
            energy_reco = np.append(energy_reco,tmp_energy_reco[mask])
            psi_reco = np.append(psi_reco,tmp_psi_reco[mask])
            azimuth_reco = np.append(azimuth_reco,tmp_azimuth_reco[mask])
            zenith_reco = np.append(zenith_reco,tmp_zenith_reco[mask])

    histograms[MCtype] = {}
    
    histograms[MCtype]['energy_rec'] = np.histogram(energy_reco,bins=bins_vars['energy_rec'],weights = weight)
    histograms[MCtype]['energy_rec_quad'] = np.histogram(energy_reco,bins=bins_vars['energy_rec'],weights = np.power(weight,2))

    histograms[MCtype]['psi_rec'] = np.histogram(psi_reco,bins=bins_vars['psi_rec'],weights = weight)
    histograms[MCtype]['psi_rec_quad'] = np.histogram(psi_reco,bins=bins_vars['psi_rec'],weights = np.power(weight,2))

    histograms[MCtype]['azimuth_rec'] = np.histogram(azimuth_reco,bins=bins_vars['azimuth_rec'],weights = weight)
    histograms[MCtype]['azimuth_rec_quad'] = np.histogram(azimuth_reco,bins=bins_vars['azimuth_rec'],weights = np.power(weight,2))
    
    histograms[MCtype]['zenith_rec'] = np.histogram(zenith_reco,bins=bins_vars['zenith_rec'],weights = weight)
    histograms[MCtype]['zenith_rec_quad'] = np.histogram(zenith_reco,bins=bins_vars['zenith_rec'],weights = np.power(weight,2))
    
    
    unbinned_data[MCtype] = {}
    unbinned_data[MCtype]['energy_rec'] = energy_reco
    unbinned_data[MCtype]['psi_rec'] = psi_reco
    unbinned_data[MCtype]['azimuth_rec'] = azimuth_reco
    unbinned_data[MCtype]['zenith_rec'] = zenith_reco
    unbinned_data[MCtype]['weight'] = weight
  
    print '... done'
    
    
for Datatype in ['Burnsample']:

    print 'start with',Datatype
    
    energy_reco = np.array([])
    psi_reco = np.array([])
    azimuth_reco = np.array([])
    zenith_reco = np.array([])
    weight = np.array([])

    for year in allFiles[Datatype].keys():
        print '  year',year

        fileData = np.load(allFiles[Datatype][year]).item()
        
        print '   ', len(fileData['all']['energy_rec']), 'events in sample'
        
        tmp_energy_reco = np.array(fileData['all']['energy_rec'])
        tmp_psi_reco = np.array(fileData['all']['psi_rec'])
        tmp_azimuth_reco = np.array(fileData['all']['azimuth_rec'])
        tmp_zenith_reco = np.array(fileData['all']['zenith_rec'])
        
        tmp_BDTScore_LE = np.array(fileData['all']['BDTScore_bb100'])
        tmp_BDTScore_HE = np.array(fileData['all']['BDTScore_ww300'])

        tmp_weight = np.array([1.]*len(tmp_energy_reco))
                
        tmp_weight[tmp_BDTScore_LE<LECut] = 0.
        tmp_weight[tmp_BDTScore_HE<HECut] = 0.
        tmp_weight[tmp_psi_reco<PsiCut] = 0.
                
        print '   ', len(tmp_weight[tmp_weight>0.]), 'events pass cuts'
                              
        energy_reco = np.append(energy_reco,tmp_energy_reco[tmp_weight>0.])
        psi_reco = np.append(psi_reco,tmp_psi_reco[tmp_weight>0.])
        azimuth_reco = np.append(azimuth_reco,tmp_azimuth_reco[tmp_weight>0.])
        zenith_reco = np.append(zenith_reco,tmp_zenith_reco[tmp_weight>0.])
    
    histograms[Datatype] = {}
    histograms[Datatype]['energy_rec'] = np.histogram(energy_reco,bins=bins_vars['energy_rec'])
    histograms[Datatype]['psi_rec'] = np.histogram(psi_reco,bins=bins_vars['psi_rec'])
    histograms[Datatype]['azimuth_rec'] = np.histogram(azimuth_reco,bins=bins_vars['azimuth_rec'])
    histograms[Datatype]['zenith_rec'] = np.histogram(zenith_reco,bins=bins_vars['zenith_rec'])

    unbinned_data[Datatype] = {}
    unbinned_data[Datatype]['energy_rec'] = energy_reco
    unbinned_data[Datatype]['psi_rec'] = psi_reco
    unbinned_data[Datatype]['azimuth_rec'] = azimuth_reco
    unbinned_data[Datatype]['zenith_rec'] = zenith_reco
    
    print '... done'

    
out = {}
out['histograms'] = histograms
out['unbinned_data'] = unbinned_data

savefile = open(outfile,'wx')
pickle.dump(out, savefile)

print 'File',outfile,'saved!'