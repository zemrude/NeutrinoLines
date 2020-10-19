import numpy as np
import time
import random
import astropy
from astropy.time import Time as apTime
from icecube import astro
import sys, os

base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')

from utils import psi_f

def genie_correction_factor(nutype_):
    if nutype_ > 0:
        return 1./0.7   
    else:
        return 1./0.3
genie_correction_factor = np.vectorize(genie_correction_factor)


def oversample(tmp_weight, tmp_nu_type,
               tmp_energy_true, tmp_energy_reco,
               tmp_zenith_true, tmp_zenith_reco,
               tmp_azimuth_true, tmp_azimuth_reco,
               nOversampling):
    
    #generate random time
    stime = time.mktime(time.strptime("1/1/2010/00/00/00", '%m/%d/%Y/%H/%M/%S'))
    etime = time.mktime(time.strptime("1/1/2015/00/00/00", '%m/%d/%Y/%H/%M/%S'))

    oversampled_psi_true = np.zeros((nOversampling,len(tmp_weight)))
    oversampled_psi_reco = np.zeros((nOversampling,len(tmp_weight)))

    for i in range(nOversampling):

        eventTime = stime + random.random() * (etime - stime)
        eventTime = apTime(eventTime,format='unix').mjd

        RA_true, decl_true = astro.dir_to_equa(tmp_zenith_true,tmp_azimuth_true,eventTime)
        oversampled_psi_true[i] = psi_f(RA_true,decl_true)

        RA_reco, decl_reco = astro.dir_to_equa(tmp_zenith_reco,tmp_azimuth_reco,eventTime)
        oversampled_psi_reco[i] = psi_f(RA_reco,decl_reco)

    # append n versions
    oversampled_weight = np.tile(tmp_weight/nOversampling,nOversampling)
    oversampled_nu_type = np.tile(tmp_nu_type,nOversampling)
    oversampled_energy_true = np.tile(tmp_energy_true,nOversampling)
    oversampled_energy_reco = np.tile(tmp_energy_reco,nOversampling)

        
    return oversampled_weight, oversampled_nu_type, oversampled_energy_true, oversampled_energy_reco, oversampled_psi_true.flatten(), oversampled_psi_reco.flatten()



def oversample_Corsika(tmp_weight, tmp_energy_reco,
                       tmp_zenith_reco, tmp_azimuth_reco,
                       nOversampling):
    
    #generate random time
    stime = time.mktime(time.strptime("1/1/2010/00/00/00", '%m/%d/%Y/%H/%M/%S'))
    etime = time.mktime(time.strptime("1/1/2015/00/00/00", '%m/%d/%Y/%H/%M/%S'))

    oversampled_psi_reco = np.zeros((nOversampling,len(tmp_weight)))

    for i in range(nOversampling):

        eventTime = stime + random.random() * (etime - stime)
        eventTime = apTime(eventTime,format='unix').mjd

        RA_reco, decl_reco = astro.dir_to_equa(tmp_zenith_reco,tmp_azimuth_reco,eventTime)
        oversampled_psi_reco[i] = psi_f(RA_reco,decl_reco)

    # append n versions
    oversampled_weight = np.tile(tmp_weight/nOversampling,nOversampling)
    oversampled_energy_reco = np.tile(tmp_energy_reco,nOversampling)

    return oversampled_weight, oversampled_energy_reco, oversampled_psi_reco.flatten()
