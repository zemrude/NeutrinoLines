#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /data/user/aguilar/icetray/meta-projects/combo/V01-00-02/build

import sys, os
from iminuit import Minuit
import numpy as np; import time; from timeout_decorator import timeout,TimeoutError
import scipy.special as sps
from optparse import OptionParser
from scipy.optimize import fsolve, root
import pickle

try:
    from physt import histogram, binnings, h1, h2, h3
except Exception as e:
    raise Expection("The required `physt` module is missing (https://github.com/janpipek/physt)")
    
# This is needed to load the environmental variable
import os, subprocess as sp, json
source = 'source /data/ana/BSM/HT_Cascade/FinalAnalysisCode/env.sh'
dump = '/usr/bin/python -c "import os, json;print json.dumps(dict(os.environ))"'
pipe = sp.Popen(['/bin/bash', '-c', '%s && %s' %(source,dump)], stdout=sp.PIPE)
env = json.loads(pipe.stdout.read())
os.environ = env

base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')

from termcolor import colored, cprint

def CorrectTypos(dm_data):
    r"""
    Previous unbinned files had a typo in the dictionary names
    This fixes it 
    
    """
    
    old_psi = "psi_rec"
    old_energy = "energy_rec"
    old_weights = "weigths"

    dm_data2 = {}
    for key in dm_data.keys():
        if key == old_energy:
            dm_data2["energy_reco"] = dm_data[old_energy]
        elif key == old_psi:
            dm_data2["psi_reco"] = dm_data[old_psi]
        elif key == old_weights:
            dm_data2["weights"] = dm_data[old_weights]

        else:
            dm_data2[key] = dm_data[key]

    return dm_data2
   
    
def RemoveElement(dm_data, eventid):
    r"""
    Remove elements from a dictionary
    """
    dm_data2 = {}
    for key in dm_data.keys():
        dm_data2[key] = np.delete(dm_data[key], eventid)
    
    return dm_data2
    
    
r"""
    SensitivityCalculation_xi.py
    ----------------------------
    Calculates the Sensitivity with both Liklihood intervals and Frequentist approach
    
    Output
    ------
    numpy file with 4 objects [sens, sensflux, freq_sens, freq_sensflux]
    sens: float, median sensitivity for likelihood interval
    freq_sens: float, median sensitivity for frequentist interval
    sensflux: dict, median sensitivity in terms of sv for annihilation or tau for decay for l.i.
    freq_sensflux: dict, median sensitivity in terms of sv for annihilation or tau for decay for f.i.

"""

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-t", "--type",default="annihilation",
                  dest="TYPE", help="Define type of PDF, default is `annihilation`")
parser.add_option("-c", "--channel",default="nue",
                  dest="CHANNEL", help="Annihilation channel, default is `nue`")
parser.add_option("-x", "--profile",default="NFW",
                  dest="PROFILE", help="DM Profile, default is `NFW`")
parser.add_option("-s", "--syst",default="nominal",
                  dest="SYST", help="Define systematic variation for signal, default is `nominal`")
parser.add_option("-m", "--mass", default='1000',
                  dest="MASS", help="DM mass in GeV, default is `1000` GeV")
parser.add_option("-a", "--lecut", default='-1',
                  dest="LECUT", help="Cut on LE BDT, default is for the HE sample")
parser.add_option("-b", "--hecut", default='0.3',
                  dest="HECUT", help="Cut on HE BDT, default is for the HE sample")
parser.add_option("-l", "--conf-level", default='90',
                  dest="CONFLEVEL", help="confidence level, default is 90\%")
parser.add_option("-d", "--llh", default='Poisson',
                  dest="LLH", help="LLH type, default is `Poisson`")
parser.add_option("-n", "--neg", default=False,
                  dest="NEG", help="Allow fitting negative ns, default is `False`")
parser.add_option("-p", "--prec", default=0.5,
                  dest="PREC", help="Frequentist precission in conf-level, default is `0.5`")

(options,args) = parser.parse_args()
        
mode = options.TYPE #annihilation, decay, background
channel = options.CHANNEL
profile = options.PROFILE
systematics = options.SYST
negative_signal = options.NEG
precision = float(options.PREC)

LECut = float(options.LECUT)
HECut = float(options.HECUT)

LEmasses = [10, 16, 20, 25, 40, 63, 100, 158, 251, 398, 631]
HEmasses = [1000, 1585, 2512, 3981, 6310, 10000, 15850, 25120, 39810]

all_masses = np.append(LEmasses,HEmasses)

mass = int(options.MASS)

if mass not in all_masses:
    print ("There are no PDFs for the mass selected")
    sys.exit(0)
    
#LECuts = [0.15,0.2]
#HECuts = [-1.0,0.3]
if channel in ['nue'] and profile in ['NFW'] and systematics in ['nominal']:
    nOversampling = 200
else:
    nOversampling = -1
#if mass == 631:
#    nOversampling = 100 #22/10/20 the PDF with oversampling 200 is not finished yet for mass 631

conf_level = int(options.CONFLEVEL)

LLH_type = options.LLH

livetime  = 28272940. + 30674072. + 31511810.5 + 31150852. + 30059465.  

out_path = os.path.join(base_path, 'sensitivity', mode, channel, systematics, LLH_type)

if not os.path.exists(out_path):
        os.popen('mkdir -p '+out_path)
        
out_file = os.path.join(out_path, 'Sensitivity_final_' + 
                        '_' + str(LLH_type) + 
                        '_' + systematics + 
                        '_' + mode + 
                        '_' + profile + 
                        '_a' + str(LECut) + 
                        '_b' + str(HECut) + 
                        '_' + channel + 
                        '_m' + str(mass) + 
                        '_p'  + str(precision) + 
                        '_CL' + str(conf_level) + 
                        '_n' + str(negative_signal) + 
                        '.npy')

if os.path.isfile(out_file):
    print ('File {} ... already exists'.format(out_file))
    sys.exit(0)

print("Writing file", colored(os.path.basename(out_file), "green"))
    
##-----------------------------#
#      Loading histograms     #
#-----------------------------#

print ("Loading the Bkg PDFs...")

bkg_file = os.path.join(base_path, "PDFs", "Data_scrambledFullSky",
                        'PDF_Data_ScrambleFullSky_LEBDT' + str(LECut)
                        + '_HEBDT' + str(HECut) + '_2D_physt.pkl')

pdfs = np.load(bkg_file,'r', allow_pickle=True, encoding='latin1')

numpy_hist_Bkg = pdfs[0] #This is the numpy pdf as usual

h_bkg = pdfs[6] #This is the physt object with new binning 


cprint("Done", "green")

if 'Gamma' in systematics:
    systematics_DM = 'nominal'
else: 
    systematics_DM = systematics

dm_file = 'Unbinned_PDF_'+systematics_DM+'_DM_'+mode+'_'+profile+'profile_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+    '_2D_'+channel+'_'+str(int(mass))+'GeV_oversampling'+str(nOversampling)+".pkl"

print ("Opening the DM file %s "%dm_file)
dm_path = os.path.join(base_path, "PDFs", "Unbinned", mode)

dm_data = np.load(os.path.join(dm_path,
                            dm_file),
                            allow_pickle=True, 
                            encoding="latin1")

#We correct the typos in the dictionnary for previous files
dm_data = CorrectTypos(dm_data)

#Check if there are events with unrealistically high weight typically everything greater than 2.
eventid = np.where(dm_data["weights"]*livetime >= 2.)
#Let's keep the normalization though,
sum_dm = np.sum(dm_data["weights"]*livetime)

dm_data = RemoveElement(dm_data, eventid)


print ("Creating the DM histograms")
qbins = h_bkg.bins

h_DM = h2(dm_data['psi_reco'], np.log10(dm_data['energy_reco']), 
             weights = livetime * dm_data['weights'], bins=qbins, 
             axis_names=["$\Psi$", "$\log_{10} E_{rec}$"])

h_DM_scrambled = h2(dm_data['scrambled_psi_reco'], np.log10(dm_data['energy_reco']), 
             weights = livetime * dm_data['weights'], bins=qbins, 
             axis_names=["$\Psi$", "$\log_{10} E_{rec}$"])


cprint("Done", "green")


    
#-----------------------------#
#    use LLHanalyser          #
#-----------------------------#
import LLHAnalyser
                                                 
analysis = LLHAnalyser.Likelihood_Analyser()
analysis.saveMoreOutput()

#Load the necessary histograms. Scrambled signal is needed for generating pseudo-trials
analysis.loadBackgroundPDF(h_bkg.frequencies, True)
analysis.loadSignalPDF(h_DM.frequencies)
analysis.loadSignalScrambledPDF(h_DM_scrambled.frequencies)   

if LLH_type == 'Effective':
    analysis.loadUncertaintyPDFs(h_bkg.errors2, h_DM.errors2)
    
elif LLH_type == 'EffectiveWithSignalSubtraction':
    analysis.loadUncertaintyPDFs(h_bkg.errors2, h_DM.errors2)
    analysis.loadSignalScrambledUncertaintyPDF(h_DM_scrambled.errors2)
    
print ("Using", colored(LLH_type, "green"), "likelihood method")
analysis.setLLHtype(LLH_type)

if (negative_signal):
    analysis.allowNegativeSignal()
    
Ntrials = 10000
sens = analysis.CalculateSensitivity(Ntrials, conf_level)

print('Likelihood interval median sensitivity xi=',sens['median'])
print('Median sensitivity xs=',sens['median']*np.sum(h_bkg)/sum_dm*10**-23)

sensflux = {}
sensflux['mass'] = mass

if mode == 'annihilation':

    sensflux['error_68_low'] = sens['error_68_low'] * np.sum(h_bkg) / sum_dm*10**-23
    sensflux['error_68_high'] = sens['error_68_high'] * np.sum(h_bkg) / sum_dm*10**-23
    sensflux['error_95_low'] = sens['error_95_low'] * np.sum(h_bkg) / sum_dm*10**-23
    sensflux['error_95_high'] = sens['error_95_high'] * np.sum(h_bkg) / sum_dm*10**-23   
    sensflux['median'] = sens['median'] * np.sum(h_bkg) / sum_dm*10**-23

elif mode == 'decay':
    sensflux['error_68_low'] = 1. / sens['error_68_low'] / np.sum(h_bkg) * sum_dm*10**28
    sensflux['error_68_high'] = 1. / sens['error_68_high'] / np.sum(h_bkg) * sum_dm*10**28
    sensflux['error_95_low'] = 1. / sens['error_95_low'] / np.sum(h_bkg) * sum_dm*10**28
    sensflux['error_95_high'] = 1. / sens['error_95_high'] / np.sum(h_bkg) * sum_dm*10**28
    sensflux['median'] = 1. / sens['median'] / np.sum(h_bkg) * sum_dm*10**28
    #sens['5Sigma_DiscoveryPotential'] = 1./analysis.CalculateDiscoveryPotential(5)/np.sum(h_bkg)*sum_dm*10**28
    #sens['Neyman_sensitivity'] = 1./sens['Neyman_sensitivity']/np.sum(h_bkg)*sum_dm*10**28


freq_sens = analysis.CalculateFrequentistSensitivity(conf_level, precision)
freq_sensflux = {}

if mode == 'annihilation':
    freq_sensflux['median'] =  freq_sens[0] * np.sum(h_bkg) / sum_dm*10**-23
elif mode == 'decay':
    freq_sensflux['median'] =  1./ freq_sens[0] / np.sum(h_bkg) * sum_dm*10**28

freq_sensflux['mass'] = mass

print ("Saving results to %s"%out_file)

np.save(out_file, [sens, sensflux, freq_sens, freq_sensflux])
