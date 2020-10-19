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

import importlib
halo_spec = importlib.util.find_spec("halo")
halo_found = halo_spec is not None

if halo_found:
    from halo import Halo

    
base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')
from utils import psi_f


print ("loading flux...")
from fileList import nfiles, allFiles
from IceCube_sim import genie_correction_factor

print ("this takes alot of time")

from fluxCalculation import phiDM_ann, phiDM_dec
print ("flux loaded")

#######################
# get and define parameters
#######################

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

outfile_path = base_path+'PDFs/ScrambledSignal/'+mode+'/'

os.popen('mkdir -p '+outfile_path)

filename_template = '/PDF_DM_FullSkyScrambled_'+mode+'_'+profile+'profile_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+    '_2D_'+channel+'_'+str(int(mass))+'GeV_oversampling'+str(nOversampling)
    
outfile = outfile_path+filename_template+'.pkl'
outfile_quad = outfile_path+filename_template+'_quad.pkl'

print ("Start) %.1f GeV %s to %s, %s profile, BDTCuts: (%.2f, %.2f) oversampling: %i"%(mass, mode, channel, profile, LECut, HECut, nOversampling))

print ("outfile : %s"%outfile)
#print 'start', mass,'GeV',mode,'to',channel, profile, 'profile','  BDTCuts:', LECut,HECut,' oversampling:',nOversampling  
if os.path.isfile(outfile):
    print (' ... file already exists')
    sys.exit(1)

    
def oversample(tmp_weight, tmp_nu_type,
               tmp_energy_reco, tmp_energy_true,
               tmp_zenith_reco, tmp_zenith_true,
               tmp_azimuth_reco, tmp_azimuth_true,
               nOversampling):
    
    #generate random time
    stime = time.mktime(time.strptime("1/1/2010/00/00/00", '%m/%d/%Y/%H/%M/%S'))
    etime = time.mktime(time.strptime("1/1/2015/00/00/00", '%m/%d/%Y/%H/%M/%S'))

    oversampled_RA_reco = np.zeros((nOversampling,len(tmp_weight)))
    oversampled_dec_reco = np.zeros((nOversampling,len(tmp_weight)))
    oversampled_RA_true = np.zeros((nOversampling,len(tmp_weight)))
    oversampled_dec_true = np.zeros((nOversampling,len(tmp_weight)))

    for i in range(nOversampling):

        eventTime = stime + random.random() * (etime - stime)
        eventTime = apTime(eventTime,format='unix').mjd

        oversampled_RA_reco[i], oversampled_dec_reco[i] = astro.dir_to_equa(tmp_zenith_reco,tmp_azimuth_reco,eventTime)
        oversampled_RA_true[i], oversampled_dec_true[i] = astro.dir_to_equa(tmp_zenith_true,tmp_azimuth_true,eventTime)

    # append n versions
    oversampled_weight = np.tile(tmp_weight/nOversampling,nOversampling)
    oversampled_nu_type = np.tile(tmp_nu_type,nOversampling)
    oversampled_energy_reco = np.tile(tmp_energy_reco,nOversampling)
    oversampled_energy_true = np.tile(tmp_energy_true,nOversampling)
        
    return oversampled_weight, oversampled_nu_type, oversampled_energy_reco, oversampled_energy_true, oversampled_RA_reco.flatten(), oversampled_RA_true.flatten(), oversampled_dec_reco.flatten(), oversampled_dec_true.flatten()

    
    
######### sensitivity
"""
bins_merge_E = 2
bins_merge_Psi = 5
conf_level = 90
LLH_type = 'effective'
systematics = 'nominal'

sens_path = base_path+'sensitivity/'+mode+'/'+channel+'/'+'nominal'+'/'

sens_oversampling=100
if mass>900:
    sens_oversampling=200

try:
    sens_file = sens_path+'Sensitivity_'+LLH_type+'_'+systematics+'_'+mode+'_'+profile+'profile_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_'+channel+'_'+str(mass)+'GeV_binsmerged_'+str(bins_merge_E)+'-'+str(bins_merge_Psi)+'_oversampling'+str(sens_oversampling)+'_'+str(conf_level)+'CL.npy'

    sensitivity = np.load(sens_file).item()['median']
    
except:
    try:
        sens_oversampling=100
        
        sens_file = sens_path+'Sensitivity_'+LLH_type+'_'+systematics+'_'+mode+'_'+profile+'profile_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_'+channel+'_'+str(mass)+'GeV_binsmerged_'+str(bins_merge_E)+'-'+str(bins_merge_Psi)+'_oversampling'+str(sens_oversampling)+'_'+str(conf_level)+'CL.npy'

        sensitivity = np.load(sens_file).item()['median']
        
    except:
        print 'could not load sensitivity!'
        sys.exit(1)    
    

print 'cross section used for signal injection:',sensitivity
"""
#############
systematics = 'nominal'

signal_energy_reco = np.array([])
signal_psi_reco = np.array([])
signal_weight = np.array([])

if halo_found:
    spinner = Halo(spinner='dots')

                      
for fileType in allFiles['MC'][systematics].keys():
    
    text = "  Opening file %s"%allFiles['MC'][systematics][fileType].split('/')[-1]
    if halo_found:
        spinner.start(text)
    else:
        print(text)
        
    fileData = np.load(allFiles['MC'][systematics][fileType], allow_pickle=True, encoding='latin1').item()

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

    tmp_weight[BDTScore_LE<LECut] = 0.
    tmp_weight[BDTScore_HE<HECut] = 0.
                           
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
        print ('   nothing more to do')
        continue
    
    text = " start oversampling: %i events with factor %i" %(len(tmp_weight[tmp_weight>0.]), nOversampling)

    if halo_found:
        spinner.succeed(" Loaded")
        spinner.start(text)
    else:
        print (text)
    
    oversampled_weight, oversampled_nu_type, oversampled_energy_reco, oversampled_energy_true, oversampled_RA_reco, oversampled_RA_true, oversampled_dec_reco, oversampled_dec_true = oversample(tmp_weight[tmp_weight>0.], tmp_nu_type[tmp_weight>0.], tmp_energy_reco[tmp_weight>0.], tmp_energy_true[tmp_weight>0.], tmp_zenith_reco[tmp_weight>0.], tmp_zenith_true[tmp_weight>0.], tmp_azimuth_reco[tmp_weight>0.], tmp_azimuth_true[tmp_weight>0.], nOversampling)
                           
    if halo_found:
        spinner.succeed(text)
        
    else:
        print(text)

        
    RA_GC = 266./180.*np.pi
                                                      
    tmp_RA_reco_corrected = oversampled_RA_reco-RA_GC
    tmp_RA_reco_corrected[tmp_RA_reco_corrected<0.] += 2*np.pi
    
    # select only events away from the GC
    #oversampled_weight[tmp_RA_reco_corrected<np.pi/2.] = 0.
    #oversampled_weight[tmp_RA_reco_corrected>np.pi*3./2.] = 0.
         
    # preliminary psi                       
    tmp_psi_true = psi_f(oversampled_RA_true[oversampled_weight>0.],oversampled_dec_true[oversampled_weight>0.])                     
                           
    if mode == 'annihilation':
        oversampled_flux = phiDM_ann(mass,channel,profile,
                                     oversampled_nu_type[oversampled_weight>0.],
                                     oversampled_energy_true[oversampled_weight>0.],
                                     tmp_psi_true)
    elif mode == 'decay':
        oversampled_flux = phiDM_dec(mass,channel,profile,
                                     oversampled_nu_type[oversampled_weight>0.],
                                     oversampled_energy_true[oversampled_weight>0.],
                                     tmp_psi_true)
    else:
        print ('mode %s not implemented! Exit'%mode)
        sys.exit()

    ## from here on work with off source events only
    new_weight = oversampled_weight[oversampled_weight>0.]*oversampled_flux
    new_energy_reco = oversampled_energy_reco[oversampled_weight>0.]
    new_dec_reco = oversampled_dec_reco[oversampled_weight>0.]
                       
    new_RA_reco = np.random.random_sample((len(new_dec_reco),)) * 2*np.pi

    new_psi_reco = psi_f(new_RA_reco,new_dec_reco)                     
                           
    signal_energy_reco = np.append(signal_energy_reco, new_energy_reco)
    signal_psi_reco = np.append(signal_psi_reco,new_psi_reco)
    signal_weight = np.append(signal_weight,new_weight)     
    
    if halo_found:
        spinner.succeed(" Finished file")
        
hist_pdf_signal = np.histogram2d(signal_energy_reco,signal_psi_reco,
                                 bins=(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = signal_weight)
hist_pdf_signal_quad = np.histogram2d(signal_energy_reco,signal_psi_reco,
                                  bins=(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = np.power(signal_weight,2))

                           
savefile = open(outfile,'wb')
pickle.dump(hist_pdf_signal, savefile)

savefile_quad = open(outfile_quad,'wb')
pickle.dump(hist_pdf_signal_quad, savefile_quad)

print (' ... done')
