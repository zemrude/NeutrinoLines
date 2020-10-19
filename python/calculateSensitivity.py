#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /data/user/aguilar/icetray/meta-projects/combo/V01-00-02/build

from time import gmtime, strftime
import numpy as np
import time
import random
import sys, os
import math
from icecube import astro
import pickle
from optparse import OptionParser

from utils import ConfidenceIntervalError, merge_bins

from termcolor import colored, cprint

import importlib
halo_spec = importlib.util.find_spec("halo")
halo_found = halo_spec is not None

if halo_found:
    from halo import Halo

# This is needed to load the environmental variable
import os, subprocess as sp, json
source = 'source /data/ana/BSM/HT_Cascade/FinalAnalysisCode/env.sh'
dump = '/usr/bin/python -c "import os, json;print json.dumps(dict(os.environ))"'
pipe = sp.Popen(['/bin/bash', '-c', '%s && %s' %(source,dump)], stdout=sp.PIPE)
env = json.loads(pipe.stdout.read())
os.environ = env
    
base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')


r"""
This scripts calcualtes the sensitivity 

"""
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-t", "--type",default="annihilation",
                  dest="TYPE", help="Define DM signal type: annihilation or decay")
parser.add_option("-c", "--channel",default="nue",
                  dest="CHANNEL", help="Annihilation channel")
parser.add_option("-x", "--profile",default="NFW",
                  dest="PROFILE", help="DM Profile")
parser.add_option("-s", "--syst",default="nominal",
                  dest="SYST", help="Define systematic variation")
parser.add_option("-a", "--lecut", default='0.2',
                  dest="LECUT", help="Cut on LE BDT")
parser.add_option("-b", "--hecut", default='0.2',
                  dest="HECUT", help="Cut on HE BDT")
parser.add_option("-m", "--mass", default='0',
                  dest="MASS", help="DM mass")
parser.add_option("-o", "--oversampling", default='1',
                  dest="OVERSAMPLING", help="Oversampling factor")
parser.add_option("-e", "--energy-binning", default='2',
                  dest="REBINE", help="rebin factor energy")
parser.add_option("-p", "--psi-binning", default='5',
                  dest="REBINPSI", help="rebin factor psi")
parser.add_option("-l", "--conf-level", default='90',
                  dest="CONFLEVEL", help="confidence level")
parser.add_option("-d", "--llh", default='Effective',
                  dest="LLH", help="LLH type")
parser.add_option("-z", "--signal", default='0', # either sensitivity or xs value
                  dest="SIGNAL", help="Amount of signal contamination")
parser.add_option("-y", "--bkg", default='simulation', # either simulation or data
                  dest="BKG", help="Background")

(options,args) = parser.parse_args()
        
mode = options.TYPE #annihilation, decay, background
channel = options.CHANNEL
profile = options.PROFILE
systematics = options.SYST

LECut = float(options.LECUT)
HECut = float(options.HECUT)
mass = int(options.MASS)
nOversampling = int(options.OVERSAMPLING)

bins_merge_E   = int(options.REBINE)
bins_merge_Psi  = int(options.REBINPSI)

conf_level = int(options.CONFLEVEL)

LLH_type = options.LLH

out_path = os.path.join(base_path, "Sensitivities", mode, channel) 

if not os.path.exists(out_path):
        os.popen('mkdir -p '+out_path)

        
out_file = os.path.join(out_path, LLH_type + '_' + systematics + '_' + mode  + '_' + profile + '_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_'+channel+'_'+str(mass)+'GeV_'+str(bins_merge_E) + '-'+str(bins_merge_Psi)+'_oversampling'+str(nOversampling)+'_'+str(conf_level)+'CL.npy')

print ("Writing output in \r %s"%out_file)
if os.path.isfile(out_file):
    cprint (" ... already exists", 'red')
    sys.exit(0)

if 'Gamma' in systematics:
    systematics_DM = 'nominal'
else: 
    systematics_DM = systematics
    
    
dm_path = os.path.join(base_path, "PDFs", "Signal", mode)
dm_scrambled_path = os.path.join(base_path, "PDFs", "ScrambledSignal", mode)

dm_file = 'PDF_' + systematics_DM + '_DM_' + mode + '_' + profile + 'profile_LEBDT' + str(LECut)+'_HEBDT'+str(HECut)+ '_2D_'+channel+'_'+str(int(mass))+'GeV_oversampling' + str(nOversampling)

dm_scrambled_file = 'PDF_' + systematics_DM + '_DM_FullSkyScrambled_' + mode + '_' + profile + 'profile_LEBDT' + str(LECut)+'_HEBDT'+str(HECut)+ '_2D_'+channel+'_'+str(int(mass))+'GeV_oversampling' + str(nOversampling)

try:
    h_DM = np.load(os.path.join(dm_path, dm_file+'.pkl'), allow_pickle = True, encoding='latin1')
except Exception as e:
    raise Exception("No Dark Matter PDF found")

    h_DM_W2 = np.load(os.path.join(dm_path, dm_file+'_quad.pkl'), allow_pickle = True, encoding='latin1')


h_DM_scrambled = np.load(os.path.join(dm_scrambled_path, dm_scrambled_file+'.pkl'), allow_pickle = True, encoding='latin1')
h_DM_scrambled_W2 = np.load(os.path.join(dm_scrambled_path, dm_scrambled_file+'_quad.pkl'), allow_pickle = True, encoding='latin1')

bkg_type = 'data'

cprint ("Loading the histograms with bkg_type: %s"%bkg_type, "green")

bkg_file = os.path.join(base_path, "PDFs", "Data_scrambledFullSky", 'PDF_Data_ScrambleFullSky_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_2D.pkl')
h_Bkg = np.load(bkg_file,'r', allow_pickle=True, encoding='latin1')
h_Bkg_W2 = [np.array([ConfidenceIntervalError(n)**2 for n in h_Bkg[0]]),h_Bkg[1], h_Bkg[2]]


livetime      = 28272940. + 30674072. + 31511810.5 + 31150852. + 30059465.
hDM     = merge_bins(h_DM, livetime,bins_merge_E, bins_merge_Psi)
hDM_W2    = merge_bins(h_DM_W2, livetime**2,bins_merge_E, bins_merge_Psi) 

hDM_Scram = merge_bins(h_DM_scrambled, livetime,bins_merge_E, bins_merge_Psi)
hDM_Scram_W2 = merge_bins(h_DM_scrambled_W2, livetime**2,bins_merge_E, bins_merge_Psi)

hbkg    = merge_bins(h_Bkg, 1.,bins_merge_E, bins_merge_Psi) 
hbkg_W2  = merge_bins(h_Bkg_W2, 1.,bins_merge_E, bins_merge_Psi)

    
r"""
    use SimpleLLHanalyser    
"""

import LLHAnalyser
                                                 
analysis = LLHAnalyser.Profile_Analyser_Normalised()
analysis.saveMoreOutput()

analysis.loadBackgroundPDF(np.array(hbkg[0]), True)
analysis.loadSignalPDF(np.array(hDM[0]))
analysis.setEmptyBins()

#This is needed even without Signal Subtraction to generate pseudo-trials
analysis.loadSignalScrambledPDF(np.array(hDM_Scram[0]) )   

cprint ("Using %s likelihood method"%LLH_type, 'green')
analysis.setLLHtype(LLH_type)

if LLH_type == 'Effective':
    analysis.loadUncertaintyPDFs(np.array(hbkg_W2[0]),np.array(hDM_W2[0]))
    analysis.allowNegativeSignal()

    
elif LLH_type == 'EffectiveWithSignalSubtraction':
    analysis.loadUncertaintyPDFs(np.array(hbkg_W2[0]),np.array(hDM_W2[0]))
    analysis.loadSignalScrambledUncertaintyPDF(np.array(hDM_Scram_W2[0]))
    analysis.allowNegativeSignal()
    

r"""

    Let's calcualte the sensitivities

"""

#We want a precison of 90 +/- 1 %
results, = analysis.CalculateFrequentistSensitivity(90, 0.5)

print ('Saving to', out_file)

np.save(out_file, results)

