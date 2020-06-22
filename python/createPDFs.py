#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /data/user/sbaur/metaprojects/icerec/build/

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
parser.add_option("-t", "--type",default="background",
                  dest="TYPE", help="Define type of PDF")
parser.add_option("-c", "--channel",default="nue",
                  dest="CHANNEL", help="Annihilation channel")
parser.add_option("-p", "--profile",default="NFW",
                  dest="PROFILE", help="DM Profile")
parser.add_option("-s", "--syst",default="nominal",
                  dest="SYST", help="Define systematic variation")
parser.add_option("-a", "--lecut", default='0.2',
                  dest="LECUT", help="Cut on LE BDT")
parser.add_option("-b", "--hecut", default='0.2',
                  dest="HECUT", help="Cut on HE BDT")
parser.add_option("-m", "--mass", default='0',
                  dest="MASS", help="DM mass")
parser.add_option("-o", "--oversampling", default='-1',
                  dest="OVERSAMPLING", help="Oversampling factor")
(options,args) = parser.parse_args()
        
mode = options.TYPE #annihilation, decay, background
channel = options.CHANNEL
DMprofile = options.PROFILE
systematicSet = options.SYST

LECut = float(options.LECUT)
HECut = float(options.HECUT)

if mode != 'background':
    mass = int(options.MASS)
    oversampling = int(options.OVERSAMPLING)
    
astro_gamma = {}
astro_gamma['nominal'] = -2.89
astro_gamma['nominalGammaUp'] = -3.09
astro_gamma['nominalGammaDown'] = -2.70
    
bins_vars={'energy_rec':np.logspace(1,5,48+1),'psi_rec':np.linspace(0.,np.pi,90+1)}

outfile_path = base_path+'PDFs/'+mode+'/'

#######################
# define PDF functions
#######################    

def makeDMPDF(mode,channel,profile,systematics,LECut,HECut,mass,nOversampling=-1):

    os.popen('mkdir -p '+outfile_path+channel+'/'+systematicSet+'/')

    filename_template = channel+'/'+systematicSet+'/PDF_'+systematics+'_DM'+mode+'_'+profile+'profile_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+    '_2D_'+channel+'_'+str(int(mass))+'GeV_oversampling'+str(nOversampling)
    
    outfile = outfile_path+filename_template+'.pkl'
    outfile_quad = outfile_path+filename_template+'_quad.pkl'
    
    print 'start',systematics, mass,'GeV',mode,'to',channel, profile, 'profile','  BDTCuts:', LECut,HECut,' oversampling:',nOversampling  
    if os.path.isfile(outfile):
        print ' ... already exists'
        return 1
    
    weight = np.array([])
    nu_type = np.array([])
    energy_true = np.array([])
    energy_reco = np.array([])
    psi_true = np.array([])
    psi_reco = np.array([])

    
    for fileType in allFiles['MC'][systematics].keys():

        fileData = np.load(allFiles['MC'][systematics][fileType]).item()

        tmp_weight = np.array(fileData['all weight'])[:,0]
        tmp_nu_type = np.array(fileData['all weight'])[:,1]

        tmp_energy_true = np.array(fileData['all']['energy_true'])
        tmp_energy_reco = np.array(fileData['all']['energy_rec'])

        tmp_zenith_true = np.array(fileData['all']['zenith_true'])
        tmp_zenith_reco = np.array(fileData['all']['zenith_rec'])

        tmp_azimuth_true = np.array(fileData['all']['azimuth_true'])
        tmp_azimuth_reco = np.array(fileData['all']['azimuth_rec'])

        tmp_psi_true = np.array(fileData['all']['psi_true'])
        tmp_psi_reco = np.array(fileData['all']['psi_rec'])
        
        BDTScore_LE = np.array(fileData['all']['BDTScore_bb100'])
        BDTScore_HE = np.array(fileData['all']['BDTScore_ww300'])

        tmp_weight[BDTScore_LE<LECut] = 0.
        tmp_weight[BDTScore_HE<HECut] = 0.

        tmp_weight[tmp_energy_true>1.5*mass] = 0.
                
        if 'nugen' in fileType:
            tmp_weight[tmp_energy_true<190] = 0.
            tmp_weight[tmp_energy_true<195] = tmp_weight[tmp_energy_true<195]*(0.2*tmp_energy_true[tmp_energy_true<195]-38.)
            tmp_weight = tmp_weight/0.5

        elif 'genie' in fileType:
            tmp_weight[tmp_energy_true>195] = 0.
            tmp_weight[tmp_energy_true>190] = tmp_weight[tmp_energy_true>190]*(-0.2*tmp_energy_true[tmp_energy_true>190]+39.)
            tmp_weight = tmp_weight*genie_correction_factor(tmp_nu_type)

        tmp_weight = tmp_weight/nfiles[systematics][fileType]

        ## skip if no events are left:
        if len(tmp_weight[tmp_weight>0.]) == 0:
            continue
        
        ## now do the oversampling for events with non-zero weight:
        
        if nOversampling==-1:
            oversampled_weight = tmp_weight[tmp_weight>0.]
            oversampled_nu_type = tmp_nu_type[tmp_weight>0.]
            oversampled_energy_true = tmp_energy_true[tmp_weight>0.]
            oversampled_energy_reco = tmp_energy_reco[tmp_weight>0.]
            oversampled_psi_true = tmp_psi_true[tmp_weight>0.]
            oversampled_psi_reco = tmp_psi_reco[tmp_weight>0.]
            
        else:
            oversampled_weight, oversampled_nu_type, oversampled_energy_true, oversampled_energy_reco, oversampled_psi_true, oversampled_psi_reco = oversample(tmp_weight[tmp_weight>0.], tmp_nu_type[tmp_weight>0.],tmp_energy_true[tmp_weight>0.], tmp_energy_reco[tmp_weight>0.], tmp_zenith_true[tmp_weight>0.], tmp_zenith_reco[tmp_weight>0.], tmp_azimuth_true[tmp_weight>0.], tmp_azimuth_reco[tmp_weight>0.], nOversampling)
        
        
        if mode == 'annihilation':
            oversampled_flux = phiDM_ann(mass,channel,profile,oversampled_nu_type,oversampled_energy_true,oversampled_psi_true)
        elif mode == 'decay':
            oversampled_flux = phiDM_dec(mass,channel,profile,oversampled_nu_type,oversampled_energy_true,oversampled_psi_true)
        else:
            print 'mode', mode, 'not implemented! Exit.'
            sys.exit()

        weight = np.append(weight,oversampled_weight*oversampled_flux)
        nu_type = np.append(nu_type,oversampled_nu_type)
        energy_true = np.append(energy_true,oversampled_energy_true)
        energy_reco = np.append(energy_reco,oversampled_energy_reco)
        psi_true = np.append(psi_true,oversampled_psi_true)
        psi_reco = np.append(psi_reco,oversampled_psi_reco)
        

    hist_pdf = np.histogram2d(energy_reco,psi_reco,bins=(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = weight)
    hist_pdf_quad = np.histogram2d(energy_reco,psi_reco,bins=(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = np.power(weight,2))
    
    savefile = open(outfile,'wx')
    pickle.dump(hist_pdf, savefile)
    
    savefile_quad = open(outfile_quad,'wx')
    pickle.dump(hist_pdf_quad, savefile_quad)
    
    print ' ... done'
    
    
def makeBkgPDF(MCtype,systematics,LECut,HECut):
    
    os.popen('mkdir -p '+outfile_path+MCtype+'/'+systematics+'/')
    filename_template = MCtype+'/'+systematics+'/PDF_'+systematics+'_'+mode+'_'+MCtype+'_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_2D'
    
    outfile = outfile_path+filename_template+'.pkl'
    outfile_quad = outfile_path+filename_template+'_quad.pkl'
    
    print 'start',mode,MCtype,'BDTCuts:', LECut,',',HECut

    if os.path.isfile(outfile):
        print ' ... already exists'
        return 1

    
    weight = np.array([])
    nu_type = np.array([])
    energy_true = np.array([])
    energy_reco = np.array([])
    psi_true = np.array([])
    psi_reco = np.array([])

    
    if 'corsika' == MCtype:
        
        fileData = np.load(allFiles['MC']['corsika']).item()

        weight = np.append(weight,fileData['all weight'])
        energy_reco = np.append(energy_reco,fileData['all']['energy_rec'])
        psi_reco = np.append(psi_reco,fileData['all']['psi_rec'])

        BDTScore_LE = np.array(fileData['all']['BDTScore_bb100'])
        BDTScore_HE = np.array(fileData['all']['BDTScore_ww300'])

        weight[BDTScore_LE<LECut] = 0.
        weight[BDTScore_HE<HECut] = 0.

    else: 

        for fileType in allFiles['MC'][systematics].keys():

            fileData = np.load(allFiles['MC'][systematics][fileType]).item()

            tmp_weight = np.array(fileData['all weight'])[:,0]
            tmp_nu_type = np.array(fileData['all weight'])[:,1]

            tmp_energy_true = np.array(fileData['all']['energy_true'])
            tmp_psi_true = np.array(fileData['all']['psi_true'])

            tmp_energy_reco = np.array(fileData['all']['energy_rec'])
            tmp_psi_reco = np.array(fileData['all']['psi_rec'])

            tmp_zenith_true = np.array(fileData['all']['zenith_true'])

            BDTScore_LE = np.array(fileData['all']['BDTScore_bb100'])
            BDTScore_HE = np.array(fileData['all']['BDTScore_ww299'])

            tmp_weight[BDTScore_LE<LECut] = 0.
            tmp_weight[BDTScore_HE<HECut] = 0.
            
            if 'nugen' in fileType:
                tmp_weight[tmp_energy_true<190] = 0.
                tmp_weight[tmp_energy_true<195] = tmp_weight[tmp_energy_true<195]*(0.2*tmp_energy_true[tmp_energy_true<195]-38.)
                tmp_weight = tmp_weight/0.5

            elif 'genie' in fileType:
                tmp_weight[tmp_energy_true>195] = 0.
                tmp_weight[tmp_energy_true>190] = tmp_weight[tmp_energy_true>190]*(-0.2*tmp_energy_true[tmp_energy_true>190]+39.)
                tmp_weight = tmp_weight*genie_correction_factor(tmp_nu_type)

            tmp_weight = tmp_weight/nfiles[systematics][fileType]

            if 'astro' == MCtype:
                if systematics not in astro_gamma.keys():
                    gamma = astro_gamma['nominal']
                else:
                    gamma = astro_gamma[systematics]
                tmp_flux = 6.45 * 1e-18 * pow(tmp_energy_true/1e5, gamma) / 6.
                
            else:
                osc = neutrinoOscillator(atmH=20)
                tmp_flux = osc_Atm_flux_weight(tmp_energy_true,tmp_zenith_true,tmp_nu_type,osc)
            
            weight = np.append(weight,tmp_weight*tmp_flux)
            nu_type = np.append(nu_type,tmp_nu_type)
            energy_true = np.append(energy_true,tmp_energy_true)
            energy_reco = np.append(energy_reco,tmp_energy_reco)
            psi_true = np.append(psi_true,tmp_psi_true)
            psi_reco = np.append(psi_reco,tmp_psi_reco)
            

    hist_pdf = np.histogram2d(energy_reco,psi_reco,bins=(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = weight)
    hist_pdf_quad = np.histogram2d(energy_reco,psi_reco,bins=(bins_vars['energy_rec'], bins_vars['psi_rec']),weights = np.power(weight,2))

    savefile = open(outfile,'wx')
    pickle.dump(hist_pdf, savefile)
    
    savefile_quad = open(outfile_quad,'wx')
    pickle.dump(hist_pdf_quad, savefile_quad)

    print ' ... done'

    
#######################
# produce requested PDFs
#######################        

if systematicSet not in allFiles['MC'].keys():
    if systematicSet not in astro_gamma.keys():
        print 'Systematic variation not defined! Exit.'
        sys.exit()
        
if mode == 'background':
    makeBkgPDF('corsika',systematicSet,LECut,HECut)
    makeBkgPDF('astro',systematicSet,LECut,HECut)
    makeBkgPDF('atm',systematicSet,LECut,HECut)
    

else:
    if 'Gamma' in systematicSet:
        print 'This is not needed! Exit.'
        sys.exit()
        
    makeDMPDF(mode,channel,DMprofile,systematicSet,LECut,HECut,mass,oversampling)
