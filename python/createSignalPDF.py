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
    This script generates the signal PDFs both for the true signal and it's scrambled version as well as the uncertainties**2

    Returns:
    --------
    2 fine binned numpy histograms 2d with signal PDF and scrambled signal PDF
    2 course binned phyist histogram with quantile binning (same statistical significance per bin)
    with signal PDF and scrambled signal PDF
    
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
parser.add_option("-t", "--type",default="annihilation",
                  dest="TYPE", help="annihilation/ decay, default is annihilation")
parser.add_option("-c", "--channel",default="nue",
                  dest="CHANNEL", help="Annihilation channel, default is nue")
parser.add_option("-p", "--profile",default="NFW",
                  dest="PROFILE", help="DM Profile, default is NFW")
parser.add_option("-a", "--lecut", default='-1',
                  dest="LECUT", help="Cut on LE BDT, default is for the HE sample")
parser.add_option("-b", "--hecut", default='0.3',
                  dest="HECUT", help="Cut on HE BDT, default is for the HE sample")
parser.add_option("-m", "--mass", default='1000',
                  dest="MASS", help="DM mass, default is 1000")
parser.add_option("-o", "--oversampling", default='200',
                  dest="OVERSAMPLING", help="Oversampling factor, default is 200")
parser.add_option("-s", "--syst",default="nominal",
                  dest="SYST", help="Define systematic variation, default is `nominal`")
parser.add_option("-n", "--host",default="local",
                  dest="HOST", help="Define where is it running, `local` or `cluster` default is `local`")

(options,args) = parser.parse_args()
   
    
    
#This takes some time due to the Jfactors
print ("Loading the Jfactors... this takes a little while")
from fluxCalculation import phiDM_ann, phiDM_dec
cprint ("Done!", "green")


mode = options.TYPE
channel = options.CHANNEL
profile = options.PROFILE

LECut = float(options.LECUT)
HECut = float(options.HECUT)

host = options.HOST

mass = int(options.MASS)
nOversampling = int(options.OVERSAMPLING)
           
bins_vars={'energy_rec':np.logspace(1,5,48+1),'psi_rec':np.linspace(0.,np.pi,90+1)}

scrambled_outfile_path = os.path.join(base_path,'PDFs/ScrambledSignal_v2/', mode)
outfile_path = os.path.join(base_path, 'PDFs/Signal_v2/' + mode)
unbinned_path = os.path.join(base_path, 'PDFs/Unbinned/' + mode)

os.popen('mkdir -p ' + scrambled_outfile_path)
os.popen('mkdir -p ' + outfile_path)
os.popen('mkdir -p ' + unbinned_path)
    
systematics = options.SYST
    
template = mode + '_' + profile + 'profile_LEBDT' + str(LECut) + '_HEBDT' + str(HECut) +    '_2D_' + channel + '_' + str(int(mass)) + 'GeV_oversampling' + str(nOversampling)
 
scrambled_filename_template =  'PDF_'+ systematics + '_DM_FullSkyScrambled_' + template
filename_template =  'PDF_' + systematics + '_DM_' + template

outfile = os.path.join(outfile_path,filename_template+'.pkl')
outfile_quad = os.path.join(outfile_path,filename_template+'_quad.pkl')

scrambled_outfile = os.path.join(scrambled_outfile_path,scrambled_filename_template+'.pkl')
scrambled_outfile_quad = os.path.join(scrambled_outfile_path, scrambled_filename_template + '_quad.pkl')

unbinned_outfile = os.path.join(unbinned_path, "Unbinned_" + filename_template +'.pkl')

print ("Signal PDF outfile : ")
cprint(os.path.basename(outfile), "green")
print ("Scrambled PDF outfile : ")
cprint(os.path.basename(scrambled_outfile), "green")
print ("Unbinned PDF outfile : ")
cprint(os.path.basename(unbinned_outfile), "green")

if os.path.isfile(outfile) or os.path.isfile(scrambled_outfile) or os.path.isfile(unbinned_outfile):
    cprint (' ... files already exist!', "red")
    sys.exit(1)



signal_energy_reco = np.array([])
signal_psi_reco = np.array([])
signal_weight = np.array([])
scrambled_signal_psi_reco = np.array([])

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

    tmp_psi_true = np.array(fileData['all']['psi_true'])
    tmp_psi_reco = np.array(fileData['all']['psi_rec'])
    
    BDTScore_LE = np.array(fileData['all']['BDTScore_bb100'])
    BDTScore_HE = np.array(fileData['all']['BDTScore_ww300'])

    #We skip event not passing the cuts
    tmp_weight[BDTScore_LE<LECut] = 0.
    tmp_weight[BDTScore_HE<HECut] = 0.
    
    #We skip events with an energy mach larger than the DM mass to speed up
    tmp_weight[tmp_energy_true>1.1*mass] = 0.  

    if 'nugen' in fileType:
        tmp_weight[tmp_energy_true<190] = 0.
        tmp_weight[tmp_energy_true<195] = tmp_weight[tmp_energy_true<195]*(0.2*tmp_energy_true[tmp_energy_true<195]-38.)
        tmp_weight = tmp_weight/0.5

    elif 'genie' in fileType:
        tmp_weight[tmp_energy_true>195] = 0.
        tmp_weight[tmp_energy_true>190] = tmp_weight[tmp_energy_true>190]*(-0.2*tmp_energy_true[tmp_energy_true>190]+39.)
        tmp_weight = tmp_weight*genie_correction_factor(tmp_nu_type)

    tmp_weight = tmp_weight/nfiles[systematics][fileType]

    print("\nTotal weight before oversampling %f "%np.sum(tmp_weight))
    
    if len(tmp_weight[tmp_weight>0.]) == 0:
        text = "... nothing to do"
        if halo_found and host == "local":
            spinner.fail(text)
        else:
            cprint (text, "red")    
        continue
   
        

    if halo_found and host == "local":
        spinner.succeed(" File loaded!")
    
    
    if (nOversampling == -1):
        cprint ("No oversampling!!")
        oversampled_weight = tmp_weight[tmp_weight>0.]
        oversampled_nu_type = tmp_nu_type[tmp_weight>0.]
        oversampled_energy_true = tmp_energy_true[tmp_weight>0.]
        oversampled_energy_reco = tmp_energy_reco[tmp_weight>0.]
       
        #We assume zenith is declination. This is needed for the scrambled psi 
        oversampled_dec_reco = tmp_zenith_reco[tmp_weight>0.] - np.pi/2. 
        oversampled_psi_true = tmp_psi_true[tmp_weight>0.]
        oversampled_psi_reco = tmp_psi_reco[tmp_weight>0.]
        
    else:
        text = " Oversampling: %i events with factor %i ..." %(len(tmp_weight[tmp_weight>0.]), nOversampling)
        print (text)

    
        oversampled_weight, oversampled_nu_type, oversampled_energy_reco, oversampled_energy_true, oversampled_RA_reco, oversampled_RA_true, oversampled_dec_reco, oversampled_dec_true = oversample(tmp_weight[tmp_weight>0.], tmp_nu_type[tmp_weight>0.], tmp_energy_reco[tmp_weight>0.], tmp_energy_true[tmp_weight>0.], tmp_zenith_reco[tmp_weight>0.], tmp_zenith_true[tmp_weight>0.], tmp_azimuth_reco[tmp_weight>0.], tmp_azimuth_true[tmp_weight>0.], nOversampling)
        
        #Calculate the true and reco angular distance to the GC for the DM flux
    
        oversampled_psi_true = psi_f(oversampled_RA_true[oversampled_weight>0.],oversampled_dec_true[oversampled_weight>0.])               
        oversampled_psi_reco = psi_f(oversampled_RA_reco[oversampled_weight>0.],oversampled_dec_reco[oversampled_weight>0.])
    
        text = " Oversampling done!"
        cprint (text, 'green')
   
    print ("\nTotal weight after oversampling %f "%np.sum(new_weight))
    
    if halo_found and host == "local":
        spinner.start(" Calculating DM fluxes...")
    else:
        print(text)
   
          
                           
    if mode == 'annihilation':
        oversampled_flux = phiDM_ann(mass,channel,profile,
                                     oversampled_nu_type[oversampled_weight>0.],
                                     oversampled_energy_true[oversampled_weight>0.],
                                     oversampled_psi_true)
    elif mode == 'decay':
        oversampled_flux = phiDM_dec(mass,channel,profile,
                                     oversampled_nu_type[oversampled_weight>0.],
                                     oversampled_energy_true[oversampled_weight>0.],
                                     oversampled_psi_true)
    else:
        cprint ('mode %s not implemented! Exit'%mode, "red")
        sys.exit()
    
    #here we calculate weight*flux
    new_weight = oversampled_weight[oversampled_weight>0.]*oversampled_flux
    new_energy_reco = oversampled_energy_reco[oversampled_weight>0.]
    new_dec_reco = oversampled_dec_reco[oversampled_weight>0.]
    #new_ra_reco = oversampled_RA_reco[oversampled_weight>0.]
    
    
    new_psi_reco = oversampled_psi_reco[oversampled_weight>0.]                   
        
        
    
    signal_energy_reco = np.append(signal_energy_reco, new_energy_reco)
    signal_psi_reco = np.append(signal_psi_reco, new_psi_reco)
    signal_weight = np.append(signal_weight,new_weight)     
    
    
    #Now we scrambled, we take a random RA from 0 - 2pi and recalculate a psi_reco
    scrambled_RA_reco = np.random.random_sample((len(new_weight),)) * 2*np.pi
    scrambled_psi_reco = psi_f(scrambled_RA_reco, new_dec_reco)                     
    scrambled_signal_psi_reco = np.append(scrambled_signal_psi_reco, scrambled_psi_reco)
    
        
    
    if halo_found and host == "local":
        spinner.succeed(" Finished with file %s!"%os.path.basename(inputfile))


# Saving the numpy arrays 
        
hist_pdf_signal = np.histogram2d(signal_energy_reco, signal_psi_reco, bins = (bins_vars['energy_rec'], bins_vars['psi_rec']), weights = signal_weight)

hist_pdf_signal_quad = np.histogram2d(signal_energy_reco,signal_psi_reco, bins =(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = np.power(signal_weight,2))

scrambled_hist_pdf_signal = np.histogram2d(signal_energy_reco, scrambled_signal_psi_reco, bins = (bins_vars['energy_rec'], bins_vars['psi_rec']), weights = signal_weight)

scrambled_hist_pdf_signal_quad = np.histogram2d(signal_energy_reco, scrambled_signal_psi_reco, bins =(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = np.power(signal_weight,2))

#Saving the unbinned data to quickly create PDFs
data = {}
data['energy_reco'] = signal_energy_reco
data['psi_reco'] = signal_psi_reco
data['weights'] = signal_weight
data['scrambled_psi_reco'] = scrambled_signal_psi_reco


savefile = open(outfile,'wb')
pickle.dump(hist_pdf_signal, savefile)

savefile_quad = open(outfile_quad,'wb')
pickle.dump(hist_pdf_signal_quad, savefile_quad)

scrambled_savefile = open(scrambled_outfile,'wb')
pickle.dump(scrambled_hist_pdf_signal, scrambled_savefile)

scrambled_savefile_quad = open(scrambled_outfile_quad,'wb')
pickle.dump(scrambled_hist_pdf_signal_quad, scrambled_savefile_quad)

unbinned_savefile = open(unbinned_outfile,'wb')
pickle.dump(data, unbinned_savefile)


print (' ... done')