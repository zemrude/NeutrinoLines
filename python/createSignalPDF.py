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

from termcolor import colored, cprint

import importlib
halo_spec = importlib.util.find_spec("halo")
halo_found = halo_spec is not None

# This is to make the waiting nicer

if halo_found:
    from halo import Halo

    
base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')
from utils import psi_f, bcolors

r"""
    This script generates the signal PDFs both for the true signal and it's scrambled version as well as the uncertainties**2

"""

from fileList import nfiles, allFiles
from IceCube_sim import genie_correction_factor

#This takes some time due to the Jfactors
cprint ("Loading the Jfactors... this takes a little while", "green")
from fluxCalculation import phiDM_ann, phiDM_dec

##################
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

    for i in range(nOversampling):
        #Generate a random time in MJD
        eventTime = stime + random.random() * (etime - stime)
        eventTime = apTime(eventTime, format='unix').mjd
        
        #recalculate RA, dec based on reconstructed zenith,azimuth and new time
        oversampled_RA_reco[i], oversampled_dec_reco[i] = astro.dir_to_equa(tmp_zenith_reco,tmp_azimuth_reco,eventTime)
        oversampled_RA_true[i], oversampled_dec_true[i] = astro.dir_to_equa(tmp_zenith_true,tmp_azimuth_true,eventTime)

    # These variables remain the same, so we append nOversampling values
    oversampled_weight = np.tile(tmp_weight/nOversampling,nOversampling)
    oversampled_nu_type = np.tile(tmp_nu_type,nOversampling)
    oversampled_energy_reco = np.tile(tmp_energy_reco,nOversampling)
    oversampled_energy_true = np.tile(tmp_energy_true,nOversampling)
        
    return oversampled_weight, oversampled_nu_type, oversampled_energy_reco, oversampled_energy_true, oversampled_RA_reco.flatten(), oversampled_RA_true.flatten(), oversampled_dec_reco.flatten(), oversampled_dec_true.flatten()

########################

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-t", "--type",default="annihilation",
                  dest="TYPE", help="annihilation/ decay")
parser.add_option("-c", "--channel",default="nue",
                  dest="CHANNEL", help="Annihilation channel")
parser.add_option("-p", "--profile",default="NFW",
                  dest="PROFILE", help="DM Profile")
parser.add_option("-a", "--lecut", default='0.15',
                  dest="LECUT", help="Cut on LE BDT")
parser.add_option("-b", "--hecut", default='0.2',
                  dest="HECUT", help="Cut on HE BDT")
parser.add_option("-m", "--mass", default='1000',
                  dest="MASS", help="DM mass")
parser.add_option("-o", "--oversampling", default='100',
                  dest="OVERSAMPLING", help="Oversampling factor")
(options,args) = parser.parse_args()
        
mode = options.TYPE
channel = options.CHANNEL
profile = options.PROFILE

LECut = float(options.LECUT)
HECut = float(options.HECUT)

mass = int(options.MASS)
nOversampling = int(options.OVERSAMPLING)
           
bins_vars={'energy_rec':np.logspace(1,5,48+1),'psi_rec':np.linspace(0.,np.pi,90+1)}

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-t", "--type",default="annihilation",
                  dest="TYPE", help="annihilation/ decay")
parser.add_option("-c", "--channel",default="nue",
                  dest="CHANNEL", help="Annihilation channel")
parser.add_option("-p", "--profile",default="NFW",
                  dest="PROFILE", help="DM Profile")
parser.add_option("-a", "--lecut", default='0.15',
                  dest="LECUT", help="Cut on LE BDT")
parser.add_option("-b", "--hecut", default='0.2',
                  dest="HECUT", help="Cut on HE BDT")
parser.add_option("-m", "--mass", default='1000',
                  dest="MASS", help="DM mass")
parser.add_option("-o", "--oversampling", default='100',
                  dest="OVERSAMPLING", help="Oversampling factor")
parser.add_option("-s", "--syst",default="nominal",
                  dest="SYST", help="Define systematic variation")

(options,args) = parser.parse_args()
        
mode = options.TYPE
channel = options.CHANNEL
profile = options.PROFILE

LECut = float(options.LECUT)
HECut = float(options.HECUT)

mass = int(options.MASS)
nOversampling = int(options.OVERSAMPLING)
           
bins_vars={'energy_rec':np.logspace(1,5,48+1),'psi_rec':np.linspace(0.,np.pi,90+1)}

scrambled_outfile_path = os.path.join(base_path,'PDFs/ScrambledSignal/', mode)
outfile_path = os.path.join(base_path, 'PDFs/Signal/' + mode)

if os.path.exists(scrambled_outfile_path) != True:       
    os.popen('mkdir -p ' + scrambled_outfile_path)
if os.path.exists(outfile_path) != True:
    os.popen('mkdir -p ' + outfile_path)
    
systematics = options.SYST
    
template = mode + '_' + profile + 'profile_LEBDT' + str(LECut) + '_HEBDT' + str(HECut) +    '_2D_' + channel + '_' + str(int(mass)) + 'GeV_oversampling' + str(nOversampling)
 
scrambled_filename_template =  'PDF_'+ systematics + '_DM_FullSkyScrambled_' + template
filename_template =  'PDF_' + systematics + '_DM_' + template

outfile = os.path.join(outfile_path,filename_template+'.pkl')
outfile_quad = os.path.join(outfile_path,filename_template+'_quad.pkl')

scrambled_outfile = os.path.join(scrambled_outfile_path,scrambled_filename_template+'.pkl')
scrambled_outfile_quad = os.path.join(scrambled_outfile_path, scrambled_filename_template + '_quad.pkl')

print ("Signal PDF outfile : ")
cprint(os.path.basename(outfile), "green")
print ("Scrambled PDF outfile : ")
cprint(os.path.basename(scrambled_outfile), "green")

if os.path.isfile(outfile) or os.path.isfile(scrambled_outfile):
    cprint (' ... files already exist!', "red")
    sys.exit(1)



signal_energy_reco = np.array([])
signal_psi_reco = np.array([])
signal_weight = np.array([])
scrambled_signal_psi_reco = np.array([])

if halo_found:
    spinner = Halo(spinner='dots')

for fileType in allFiles['MC'][systematics].keys():
    
    inputfile = allFiles['MC'][systematics][fileType]
    text = "Reading file... %s "%os.path.basename(inputfile)
   
    if halo_found:
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

    if len(tmp_weight[tmp_weight>0.]) == 0:
        text = "... nothing to do"
        if halo_found:
            spinner.fail(text)
        else:
            cprint (text, "red")    
        continue
   
    if (nOversampling < 0.):
        if halo_found:
            spinner.succeed(" Finished file")
        else: 
            print("Finished file")
        continue
        
    text = " Oversampling: %i events with factor %i ..." %(len(tmp_weight[tmp_weight>0.]), nOversampling)

    if halo_found:
        spinner.succeed(" File loaded!")
        spinner.start(text)
    else:
        print (text)

    oversampled_weight, oversampled_nu_type, oversampled_energy_reco, oversampled_energy_true, oversampled_RA_reco, oversampled_RA_true, oversampled_dec_reco, oversampled_dec_true = oversample(tmp_weight[tmp_weight>0.], tmp_nu_type[tmp_weight>0.], tmp_energy_reco[tmp_weight>0.], tmp_energy_true[tmp_weight>0.], tmp_zenith_reco[tmp_weight>0.], tmp_zenith_true[tmp_weight>0.], tmp_azimuth_reco[tmp_weight>0.], tmp_azimuth_true[tmp_weight>0.], nOversampling)

    text = " Oversampling done!"
    if halo_found:
        spinner.succeed(text)
        spinner.start(" Calculating DM fluxes...")

    else:
        print(text)
   
    #Calculate the true angular distance to the GC for the DM flux
    
    oversampled_psi_true = psi_f(oversampled_RA_true[oversampled_weight>0.],oversampled_dec_true[oversampled_weight>0.])                     
                           
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
    new_ra_reco = oversampled_RA_reco[oversampled_weight>0.]
    
    oversampled_psi_reco = psi_f(oversampled_RA_reco[oversampled_weight>0.],oversampled_dec_reco[oversampled_weight>0.])
    
    new_psi_reco = oversampled_psi_reco[oversampled_weight>0.]                   

    signal_energy_reco = np.append(signal_energy_reco, new_energy_reco)
    signal_psi_reco = np.append(signal_psi_reco, oversampled_psi_reco)
    signal_weight = np.append(signal_weight,new_weight)     
    
    
    #Now we scrambled, we take a random RA from 0 - 2pi and recalculate a psi_reco
    scrambled_RA_reco = np.random.random_sample((len(new_ra_reco),)) * 2*np.pi
    scrambled_psi_reco = psi_f(scrambled_RA_reco, new_dec_reco)                     
    scrambled_signal_psi_reco = np.append(scrambled_signal_psi_reco, scrambled_psi_reco)
    
        
    
    if halo_found:
        spinner.succeed(" Finished with file!")

hist_pdf_signal = np.histogram2d(signal_energy_reco, signal_psi_reco, bins = (bins_vars['energy_rec'], bins_vars['psi_rec']), weights = signal_weight)

hist_pdf_signal_quad = np.histogram2d(signal_energy_reco,signal_psi_reco, bins =(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = np.power(signal_weight,2))

scrambled_hist_pdf_signal = np.histogram2d(signal_energy_reco, scrambled_signal_psi_reco, bins = (bins_vars['energy_rec'], bins_vars['psi_rec']), weights = signal_weight)

scrambled_hist_pdf_signal_quad = np.histogram2d(signal_energy_reco, scrambled_signal_psi_reco, bins =(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = np.power(signal_weight,2))



savefile = open(outfile,'wb')
pickle.dump(hist_pdf_signal, savefile)

savefile_quad = open(outfile_quad,'wb')
pickle.dump(hist_pdf_signal_quad, savefile_quad)

scrambled_savefile = open(scrambled_outfile,'wb')
pickle.dump(scrambled_hist_pdf_signal, scrambled_savefile)

scrambled_savefile_quad = open(scrambled_outfile_quad,'wb')
pickle.dump(scrambled_hist_pdf_signal_quad, scrambled_savefile_quad)

print (' ... done')