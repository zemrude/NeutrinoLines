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

from termcolor import colored, cprint


#This is used for the quantile binning
try:
    from physt import histogram, binnings, h1, h2, h3
except Exception as e:
    raise Expection("The required `physt` module is missing (https://github.com/janpipek/physt)")
    
# This is to make the waiting nicer
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

from fileList import nfiles, allFiles
from utils import psi_f

r"""
    
    Create the data pdf from scrambled data. It can use the full sample or the burn sample \
    14/10/2020 To do: Burnsample numpy files do not contain ra_rec, these need to be rebuilt

    Returns:
    --------
    A fine binned numpy histogram 2d
    2 physt binned phyist histogram with quantile binning (same statistical significance per bin)
     with the same old binning, the other with the quantile binning
     
    
"""
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-t", "--type", default="Data",
                  dest="TYPE", help="Type : Data or Burnsample, default is data")
parser.add_option("-a", "--lecut", default='-1.0',
                  dest="LECUT", help="Cut on LE BDT, default is for the HE sample")
parser.add_option("-b", "--hecut", default='0.3',
                  dest="HECUT", help="Cut on HE BDT, default is for the HE sample")
(options,args) = parser.parse_args()
        
LECut = float(options.LECUT)
HECut = float(options.HECUT)
Datatype = options.TYPE

#This are the fine binning for the standard numpy pdf
bins_vars={'energy_rec':np.logspace(1,5,48+1),'psi_rec':np.linspace(0.,np.pi,90+1)}
bins_vars2={'energy_rec':np.logspace(1,5,24+1),'psi_rec':np.linspace(0.,np.pi,10+1)}

outfile_path = os.path.join(base_path, 'PDFs', 'Data_scrambledFullSky')

os.popen('mkdir -p '+outfile_path)
   
outfile = os.path.join(outfile_path, 'PDF_' + Datatype + '_ScrambleFullSky_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_2D_physt'+'.pkl')

if os.path.isfile(outfile):
    print (' ... already exists')
 #   sys.exit(1)

livetime = {}
livetime['Burnsample'] = 1116125.821572
livetime['Data'] = 28272940. + 30674072. + 31511810.5 + 31150852. + 30059465.  

######## Arrays for all years ######

energy_reco = np.array([])
psi_reco = np.array([])
weight = np.array([])

if halo_found:
    spinner = Halo(spinner='dots')

for year in allFiles[Datatype].keys():

    text = " Doing year %s "%year
   
    if halo_found:
        spinner.start(text)
    else:
        print(text)
        
    
    fileData = np.load(allFiles[Datatype][year], allow_pickle=True, encoding='latin1').item()

    text = ' %i events in sample' %len(fileData['all']['energy_rec'])
    if halo_found:
        spinner.stop_and_persist(text=text)
    else:
        print (text)
    
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

    if halo_found:
        spinner.stop()
        spinner.start("Scrambling RA")
    else:
        print ("Scrambling RA")

    
 
    new_energy_reco = tmp_energy_reco
    new_dec_reco = tmp_dec_reco
    new_weight = [1.]*len(new_dec_reco)
    
 
    new_RA_reco = np.random.random_sample((len(new_dec_reco),)) * 2*np.pi
    new_psi_reco = psi_f(new_RA_reco,new_dec_reco)     
    
    energy_reco = np.append(energy_reco, new_energy_reco)
    psi_reco = np.append(psi_reco, new_psi_reco)
    weight = np.append(weight, new_weight)
    
    if halo_found:
        spinner.succeed(" Finished with year!")

#We create the numpy histogram in linear energy and psi
 
print ("Max and min energies")
print (np.max(np.log10(energy_reco)), np.min(np.log10(energy_reco)))
    
    
np_hist_pdf = np.histogram2d(energy_reco,psi_reco, 
                          bins=(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = weight)

#We create the physt histogram in linear psi and log10(energy)! Be aware of that, with final binning but quantile division
#h2 (x , y, "binning_method", (binsx, binsy), [labelsx, labelsy])

hphyst1 = h2(psi_reco, np.log10(energy_reco),"quantile", (10, 24), axis_names=["$\Psi$", "$\log_{10} E_{rec}$"])
hphyst2 = h2(psi_reco, energy_reco, bins=(bins_vars2['psi_rec'], bins_vars2['energy_rec']), axis_names=["$\Psi$", "$\log_{10} E_{rec}$"])
hphyst3 = h2(psi_reco, np.log10(energy_reco),"quantile", (9, 9), axis_names=["$\Psi$", "$\log_{10} E_{rec}$"])
#25/11/2020 New binning, more bins (4x) in energy
hphyst4 = h2(psi_reco, np.log10(energy_reco),"quantile", (9, 36), axis_names=["$\Psi$", "$\log_{10} E_{rec}$"])

#Let's go nuts, and do 8x binings
hphyst5 = h2(psi_reco, np.log10(energy_reco),"quantile", (9, 72), axis_names=["$\Psi$", "$\log_{10} E_{rec}$"])

#Let's go nuts, and do 8x binings
hphyst6 = h2(psi_reco, np.log10(energy_reco),"quantile", (18, 72), axis_names=["$\Psi$", "$\log_{10} E_{rec}$"])


savefile = open(outfile,'wb')

pickle.dump([np_hist_pdf, hphyst1, hphyst2, hphyst3, hphyst4, hphyst5, hphyst6], savefile)

cprint (' ... Done!', 'green')
