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

import importlib
halo_spec = importlib.util.find_spec("halo")
halo_found = halo_spec is not None

# This is to make the waiting nicer

if halo_found:
    from halo import Halo
    
tqdm_spec = importlib.util.find_spec("tqdm")
tqdm_found = tqdm_spec is not None

if tqdm_found:
    import tqdm as tq

# This is needed to load the environmental variable

import os, subprocess as sp, json
source = 'source /data/ana/BSM/HT_Cascade/FinalAnalysisCode/env.sh'
dump = '/usr/bin/python -c "import os, json;print json.dumps(dict(os.environ))"'
pipe = sp.Popen(['/bin/bash', '-c', '%s && %s' %(source,dump)], stdout=sp.PIPE)
env = json.loads(pipe.stdout.read())
os.environ = env
    
base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')
from utils import psi_f, bcolors

r"""
    This script generates an oversampler file from simulation

    Returns:
    --------
    Output numpy file with an oversampling file
    
    
"""

from fileList import nfiles, allFiles
from IceCube_sim import genie_correction_factor


def oversample(tmp_weight, tmp_nu_type,
               tmp_energy_reco, tmp_energy_true,
               tmp_zenith_reco, tmp_zenith_true,
               tmp_azimuth_reco, tmp_azimuth_true,
               nOversampling):

   
    
    r"""
    Function to do the oversampling. For each event it's repeated nsampling times with different times (ie RA)
    
    Parameters
    ----------
    event weight, event_type,
    energy_reco, energy_true,
    
    
    Returns
    -------
    
    """
    #generate random time
    stime = time.mktime(time.strptime("1/1/2010/00/00/00", '%m/%d/%Y/%H/%M/%S'))
    etime = time.mktime(time.strptime("1/1/2015/00/00/00", '%m/%d/%Y/%H/%M/%S'))

    oversampled_RA_reco = np.zeros((nOversampling,len(tmp_weight)))
    oversampled_dec_reco = np.zeros((nOversampling,len(tmp_weight)))
    oversampled_RA_true = np.zeros((nOversampling,len(tmp_weight)))
    oversampled_dec_true = np.zeros((nOversampling,len(tmp_weight)))

    if tqdm_found and host =='local':
        pbar = tq.tqdm(total =  nOversampling)
        
    for i in range(nOversampling):
        #Generate a random time in MJD
        eventTime = stime + random.random() * (etime - stime)
        eventTime = apTime(eventTime, format='unix').mjd
        
        if tqdm_found and host =='local':
            pbar.update(1)
        
        #recalculate RA, dec based on reconstructed zenith,azimuth and new time
        oversampled_RA_reco[i], oversampled_dec_reco[i] = astro.dir_to_equa(tmp_zenith_reco,tmp_azimuth_reco,eventTime)
        oversampled_RA_true[i], oversampled_dec_true[i] = astro.dir_to_equa(tmp_zenith_true,tmp_azimuth_true,eventTime)

    if tqdm_found and host =='local':
        pbar.close()
        
    # These variables remain the same, so we append nOversampling values
    oversampled_weight = np.tile(tmp_weight/nOversampling,nOversampling)
    oversampled_nu_type = np.tile(tmp_nu_type,nOversampling)
    oversampled_energy_reco = np.tile(tmp_energy_reco,nOversampling)
    oversampled_energy_true = np.tile(tmp_energy_true,nOversampling)
        
    return oversampled_weight, oversampled_nu_type, oversampled_energy_reco, oversampled_energy_true, oversampled_RA_reco.flatten(), oversampled_RA_true.flatten(), oversampled_dec_reco.flatten(), oversampled_dec_true.flatten()



#------------------------


usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-a", "--lecut", default='-1',
                  dest="LECUT", help="Cut on LE BDT, default is for the HE sample")
parser.add_option("-b", "--hecut", default='0.3',
                  dest="HECUT", help="Cut on HE BDT, default is for the HE sample")
parser.add_option("-o", "--oversampling", default='200',                  
                  dest="OVERSAMPLING", help="Oversampling factor, default is 200")
parser.add_option("-s", "--syst",default="nominal",
                  dest="SYST", help="Define systematic variation, default is `nominal`")
parser.add_option("-n", "--host",default="local",
                  dest="HOST", help="Define where is it running, `local` or `cluster` default is `local`")

(options,args) = parser.parse_args()


LECut = float(options.LECUT)
HECut = float(options.HECUT)

host = options.HOST

nOversampling = int(options.OVERSAMPLING)
           

outfile_path = os.path.join(base_path, 'PDFs', 'OverSample')

os.popen('mkdir -p ' + outfile_path)
    
systematics = options.SYST
    
filename = 'Oversampled_a' + str(LECut) + '_b' + str(HECut) + '_o' + str(nOversampling)


outfile = os.path.join(outfile_path, filename + '.pkl')

print ("Outfile : ")
cprint(os.path.basename(outfile), "green")

if os.path.isfile(outfile):
    cprint (' ... files already exist!', "red")
    sys.exit(1)
    
#REconstructed arrays
signal_energy_reco = np.array([])
signal_psi_reco = np.array([])
signal_RA_reco = np.array([])
signal_dec_reco = np.array([])
scrambled_signal_psi_reco = np.array([])
#True arrays
signal_weight = np.array([])
signal_nu_type = np.array([])
signal_energy_true = np.array([])
signal_psi_true = np.array([])
signal_RA_true = np.array([])
signal_dec_true = np.array([])

if halo_found and host == "local":
    spinner = Halo(spinner='dots')

for fileType in allFiles['MC'][systematics].keys():
    
    inputfile = allFiles['MC'][systematics][fileType]
    text = "Reading file... %s "%os.path.basename(inputfile)
   
    if halo_found and host == "local":
        spinner.start(text)
    else:
        print(text)
        
    fileData = np.load(inputfile, allow_pickle=True, encoding='latin1').item()

    tmp_weight = np.array(fileData['all weight'])[:,0]
    tmp_nu_type = np.array(fileData['all weight'])[:,1]

    tmp_energy_reco = np.array(fileData['all']['energy_rec'])
    tmp_energy_true = np.array(fileData['all']['energy_true'])

    tmp_zenith_reco = np.array(fileData['all']['zenith_rec'])
    tmp_azimuth_reco = np.array(fileData['all']['azimuth_rec'])

    tmp_zenith_true = np.array(fileData['all']['zenith_true'])
    tmp_azimuth_true = np.array(fileData['all']['azimuth_true'])

    BDTScore_LE = np.array(fileData['all']['BDTScore_bb100'])
    BDTScore_HE = np.array(fileData['all']['BDTScore_ww300'])

    #We skip event not passing the cuts
    tmp_weight[BDTScore_LE<LECut] = 0.
    tmp_weight[BDTScore_HE<HECut] = 0.
    #mask = (BDTScore_LE < LECut) or (BDTScore_HE < HECut)
    
    #We skip event not passing the cuts
    #tmp_weight[mask] = 0.
   
    if 'nugen' in fileType:
        tmp_weight[tmp_energy_true<190] = 0.
        tmp_weight[tmp_energy_true<195] = tmp_weight[tmp_energy_true<195]*(0.2*tmp_energy_true[tmp_energy_true<195]-38.)
        tmp_weight = tmp_weight/0.5

    elif 'genie' in fileType:
        tmp_weight[tmp_energy_true>195] = 0.
        tmp_weight[tmp_energy_true>190] = tmp_weight[tmp_energy_true>190]*(-0.2*tmp_energy_true[tmp_energy_true>190]+39.)
        tmp_weight = tmp_weight*genie_correction_factor(tmp_nu_type)

    tmp_weight = tmp_weight/nfiles[systematics][fileType]
    
    tmp_weight = tmp_weight/nfiles[systematics][fileType]

  
    if len(tmp_weight[tmp_weight>0.]) == 0:
        text = "... nothing to do"
        if halo_found and host == "local":
            spinner.fail(text)
        else:
            cprint (text, "red")    
            continue
      
    text = " File loaded! \nTotal weight before oversampling %f"%np.sum(tmp_weight)
    if halo_found and host == "local":
        spinner.succeed(text)
    else:
        print (text)
        
    text = " Oversampling: %i events with factor %i ..." %(len(tmp_weight[tmp_weight>0.]), nOversampling)
    print (text)
    
    oversampled_weight, oversampled_nu_type, oversampled_energy_reco, oversampled_energy_true, oversampled_RA_reco, oversampled_RA_true, oversampled_dec_reco, oversampled_dec_true = oversample(tmp_weight[tmp_weight>0.], tmp_nu_type[tmp_weight>0.], tmp_energy_reco[tmp_weight>0.], tmp_energy_true[tmp_weight>0.], tmp_zenith_reco[tmp_weight>0.], tmp_zenith_true[tmp_weight>0.], tmp_azimuth_reco[tmp_weight>0.], tmp_azimuth_true[tmp_weight>0.], nOversampling)

    new_weight = oversampled_weight[oversampled_weight>0.]
    new_nu_type = oversampled_nu_type[oversampled_weight>0.]
    new_energy_true = oversampled_energy_true[oversampled_weight>0.]
    new_dec_true = oversampled_dec_true[oversampled_weight>0.]
    new_ra_true = oversampled_RA_true[oversampled_weight>0.]
   
    new_energy_reco = oversampled_energy_reco[oversampled_weight>0.]
    new_dec_reco = oversampled_dec_reco[oversampled_weight>0.]
    new_ra_reco = oversampled_RA_reco[oversampled_weight>0.]
   
    
    
    

    text = " Oversampling done! \n Total weight after oversampling %f"%np.sum(new_weight)
    cprint (text, 'green')

    new_psi_reco = psi_f(oversampled_RA_reco[oversampled_weight>0.],oversampled_dec_reco[oversampled_weight>0.])
    
   
    new_psi_true = psi_f(oversampled_RA_true[oversampled_weight>0.],oversampled_dec_true[oversampled_weight>0.])  
    
    #True variables
    signal_energy_true = np.append(signal_energy_true, new_energy_true)
    signal_weight = np.append(signal_weight,new_weight)     
    signal_nu_type = np.append(signal_nu_type, new_nu_type)
    signal_psi_true = np.append(signal_psi_true, new_psi_true)
    signal_RA_true = np.append(signal_RA_true,new_ra_reco)     
    signal_dec_true = np.append(signal_dec_true,new_ra_reco)     
    
    #Rec variables
    signal_energy_reco = np.append(signal_energy_reco, new_energy_reco)
    signal_psi_reco = np.append(signal_psi_reco, new_psi_reco)
    signal_RA_reco = np.append(signal_RA_reco,new_ra_reco)     
    signal_dec_reco = np.append(signal_dec_reco,new_dec_reco)     
    
    #Now we scrambled, we take a random RA from 0 - 2pi and recalculate a psi_reco
    scrambled_RA_reco = np.random.random_sample((len(new_ra_reco),)) * 2*np.pi
    scrambled_psi_reco = psi_f(scrambled_RA_reco, new_dec_reco)                     

    scrambled_signal_psi_reco = np.append(scrambled_signal_psi_reco, scrambled_psi_reco)
    
#Saving the unbinned data to quickly create PDFs
data = {}
#reconstructed varialbes
data['energy_reco'] = signal_energy_reco
data['psi_reco'] = signal_psi_reco
data['scrambled_psi_reco'] = scrambled_signal_psi_reco
data['RA_reco'] = signal_RA_reco
data['dec_reco'] = signal_dec_reco

#true variables
data['nu_type'] = signal_nu_type
data['weights'] = signal_weight
data['energy_true'] = signal_energy_true
data['psi_true'] = signal_psi_true
data['RA_true'] = signal_RA_true
data['dec_true'] = signal_dec_true

savefile = open(outfile,'wb')
pickle.dump(data, savefile)


print (' ... done')
