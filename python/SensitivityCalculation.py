#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /data/user/sbaur/metaprojects/icerec/build/

import sys, os
from iminuit import Minuit
import numpy as np; import time; from timeout_decorator import timeout,TimeoutError
import scipy.special as sps
from optparse import OptionParser
from scipy.optimize import fsolve, root

base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')

#######################
# get and define parameters
#######################

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-t", "--type",default="annihilation",
                  dest="TYPE", help="Define type of PDF")
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
parser.add_option("-l", "--conf-level", default='95',
                  dest="CONFLEVEL", help="confidence level")
parser.add_option("-d", "--llh", default='effective',
                  dest="LLH", help="LLH type")

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

time_tot     = 28272940. + 30674072. + 31511810.5 + 31150852. + 30059465.  

out_path = base_path+'sensitivity/'+mode+'/'+channel+'/'+systematics+'/'
if not os.path.exists(out_path):
    os.popen('mkdir -p '+out_path)
    
out_file = out_path+'Sensitivity_'+LLH_type+'_'+systematics+'_'+mode+'_'+profile+'profile_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_'+channel+'_'+str(mass)+'GeV_binsmerged_'+str(bins_merge_E)+'-'+str(bins_merge_Psi)+'_oversampling'+str(nOversampling)+'_'+str(conf_level)+'CL.npy'

if os.path.isfile(out_file):
    print ' ... already exists'
    sys.exit(0)

##-----------------------------#
#      Loading histograms     #
#-----------------------------#
print "Loading the histograms..."

if 'Gamma' in systematics:
    systematics_DM = 'nominal'
else: 
    systematics_DM = systematics

DM_file_template = '/PDF_'+systematics_DM+'_DM'+mode+'_'+profile+'profile_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+    '_2D_'+channel+'_'+str(int(mass))+'GeV_oversampling'+str(nOversampling)

h_DM = np.load(open(base_path+'PDFs/'+mode+'/'+channel+'/'+systematics_DM+DM_file_template+'.pkl'))
h_DM_quad = np.load(open(base_path+'PDFs/'+mode+'/'+channel+'/'+systematics_DM+DM_file_template+'_quad.pkl'))

h_Bkg = {}
for bkg in ['corsika','atm','astro']:
    Bkg_file = '/PDF_'+systematics+'_background_'+bkg+'_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_2D'
    h_Bkg[bkg] = np.load(open(base_path+'PDFs/background/'+bkg+'/'+systematics+Bkg_file+'.pkl'))
    h_Bkg[bkg+'_quad'] = np.load(open(base_path+'PDFs/background/'+bkg+'/'+systematics+Bkg_file+'_quad.pkl'))
    #h_Bkg[bkg] = np.load(open(base_path+'PDFs/all_old_stages/background/'+bkg+'/'+Bkg_file+'.pkl'))
    #h_Bkg[bkg+'_quad'] = np.load(open(base_path+'PDFs/all_old_stages/background/'+bkg+'/'+Bkg_file+'_quad.pkl'))    

#-------------------------------------#
#-------------------------------------#


def merge_bins(histo , time, nD_ = 2): 
    ''' 
    This function merges the bins in a histogram
    '''
    grid_E   = histo[1] 
    grid_E   = [grid_E[i] for i in range(len(grid_E)) if i%bins_merge_E ==0]
    if nD_ == 2 : 
        grid_psi = histo[2]
        grid_psi = [grid_psi[i] for i in range(len(grid_psi)) if i%bins_merge_Psi ==0]
    lE=len(grid_E);lpsi=len(grid_psi)
    if nD_ == 2:
        tab_events = np.zeros((lE-1,lpsi-1))
        for j in range(bins_merge_E):
                for k in range(lE-1):
                        for j2 in range(bins_merge_Psi):
                                for k2 in range(lpsi-1):
                                        tab_events[k,k2] = tab_events[k,k2] + histo[0][j+k*bins_merge_E,j2+k2*bins_merge_Psi]*time
        return tab_events, grid_E, grid_psi
    else:
        tab_events=np.zeros((lE-1))
        for j in range(bins_merge_E):
                for k in range(lE-1):
                        tab_events[k] = tab_events[k] + histo[0][j+k*bins_merge_Psi]*time
        return tab_events, grid_E
    
    
hatm    = merge_bins(h_Bkg['atm'], time_tot) 
hastro    = merge_bins(h_Bkg['astro'], time_tot) 
hcorsi  = merge_bins(h_Bkg['corsika'], time_tot) 
hDM     = merge_bins(h_DM, time_tot) 

hatm_quad   = merge_bins(h_Bkg['atm_quad'],time_tot**2) 
hastro_quad    = merge_bins(h_Bkg['astro_quad'], time_tot**2) 
hcorsi_quad = merge_bins(h_Bkg['corsika_quad'], time_tot**2)
hDM_quad    = merge_bins(h_DM_quad, time_tot**2) 

tab_events_corsika = np.array(hcorsi[0]).flatten()
tab_events_atm = np.array(hatm[0]).flatten()
tab_events_astro = np.array(hastro[0]).flatten()
tab_events_dm = np.array(hDM[0]).flatten()

tabquad_events_corsika = np.array(hcorsi_quad[0]).flatten()
tabquad_events_atm = np.array(hatm_quad[0]).flatten()
tabquad_events_astro = np.array(hastro_quad[0]).flatten()
tabquad_events_dm = np.array(hDM_quad[0]).flatten()

#-----------------------------#
#        define likelihood    #
#-----------------------------#

def model_array(n1,n2,n3,nsig):
    return n1*tab_events_corsika + n2*tab_events_atm + n3*tab_events_astro + nsig*tab_events_dm

def model_array_quad(n1,n2,n3,nsig):
    return n1**2*tabquad_events_corsika + n2**2*tabquad_events_atm + n3**2*tabquad_events_astro + nsig**2*tabquad_events_dm

def LLH(n1,n2,n3,nsig):
  
    k      = observation
    mu     = model_array(n1,n2,n3,nsig)
    mu2    = np.power(mu,2)
    sigma2 = model_array_quad(n1,n2,n3,nsig)

    bins_to_use = (mu>0.)
    
    if LLH_type == 'effective':
    
        alpha = mu2[bins_to_use]/sigma2[bins_to_use] +1.
        beta  = mu[bins_to_use]/sigma2[bins_to_use]

        # L_eff components
        values = [
            alpha*np.log(beta),
            sps.loggamma(k[bins_to_use]+alpha).real,
            #-sps.loggamma(k+1.0).real,
            -(k[bins_to_use]+alpha)*np.log1p(beta),
            -sps.loggamma(alpha).real,
            ]
        
    elif LLH_type == 'Poisson':
        values = k[bins_to_use]*np.log(mu[bins_to_use])-mu[bins_to_use]
    
    
    return -np.sum(values)


##------------------------------------#
## pseudo x-periments:TS distribution #
##------------------------------------#
nMC=100000
print "Generate",nMC, "background only pseudo-experiments..."

limit={}
limit['lower']=np.zeros(nMC)
limit['upper']=np.zeros(nMC)
fitResult = []
TOut_errors=np.zeros(nMC)

TSdist = np.array([])
allFits = np.array([])
upperLimit = np.array([])

start_time=time.time()

for j in range(nMC):
    #try:
        ## generate nMC poisson distributed pseudo-data
        bkg_exp = 1.*tab_events_corsika + 1.*tab_events_atm + 1.*tab_events_astro
        
        #N_bkg = np.sum(bkg_exp)
        #bkg_exp = bkg_exp/N_bkg
               
        observation=np.zeros(np.shape(bkg_exp)[0])
        for i in range(len(observation)):
            observation[i]=np.random.poisson(bkg_exp[i])
            #observation[i]=np.random.poisson(N_bkg)*bkg_exp[i]
            
        # best fit
        LLHmin_DM=Minuit(LLH,
             nsig=1e-3,n1=1.,n2=1.,n3=1.,
             error_nsig=.1,error_n1=.1,error_n2=.1,error_n3=.1,
             limit_nsig=(-1.,100.),limit_n1=(0.,10.),limit_n2=(0.,10.),limit_n3=(0.,10.),
             errordef=.5,print_level=0)  
        LLHmin_DM.migrad()
        
        DM_fitarg_2=LLHmin_DM.fitarg
        LLHmin_DM_2=Minuit(LLH, errordef=.5, print_level=0,pedantic=True, **DM_fitarg_2)
        LLHmin_DM_2.migrad()

        bestFit = {}
        bestFit['n1']=LLHmin_DM_2.fitarg['n1']
        bestFit['n2']=LLHmin_DM_2.fitarg['n2']
        bestFit['n3']=LLHmin_DM_2.fitarg['n3']
        bestFit['nsig']=LLHmin_DM_2.fitarg['nsig']
        bestFit['LLH']=LLH(bestFit['n1'],bestFit['n2'],bestFit['n3'],bestFit['nsig'])
        
        LLHmin_ref=Minuit(LLH,
                 nsig=0.,fix_nsig = True,
                 n1=1.,n2=1.,n3=1.,
                 error_nsig=.1,error_n1=.1,error_n2=.1,error_n3=.1,
                 limit_nsig=(-10.,100.),limit_n1=(0.,10.),limit_n2=(0.,10.),limit_n3=(0.,10.),
                 errordef=.5,print_level=0)  
        LLHmin_ref.migrad()
            
        TS = 0.
        if bestFit['nsig'] > 0.:
            TS = 2*(LLH(LLHmin_ref.fitarg['n1'],LLHmin_ref.fitarg['n2'],LLHmin_ref.fitarg['n3'],LLHmin_ref.fitarg['nsig'])-bestFit['LLH'])
            
            
        # upper limit calulation
        nIterations =0
        eps_TS=0.005
        eps_param=0.05

        deltaTS = 2.71
        if conf_level==90:
            deltaTS = 1.64
        elif conf_level==95:
            deltaTS = 2.71
        else:
            print("Chosen CL is not defined!")
            sys.exit(1)
            
        param_low=bestFit['nsig']
        param_up=bestFit['nsig']
        param_mean=bestFit['nsig']
        
        dTS=0
        cc=1
        while((dTS<deltaTS) and (nIterations<100)):
            nIterations += 1 

            param_up=param_up+3.*np.abs(param_up)

            LLHmin_fix=Minuit(LLH,
                 nsig=param_up,fix_nsig = True,
                 n1=1.,n2=1.,n3=1.,
                 error_nsig=.1,error_n1=.1,error_n2=.1,error_n3=.1,
                 limit_nsig=(-10.,100.),limit_n1=(0.,10.),limit_n2=(0.,10.),limit_n3=(0.,10.),
                 errordef=.5,print_level=0)  
            LLHmin_fix.migrad()

            if param_up <0.:
                TS_fix = 0.
            else:
                TS_fix = 2*(LLH(LLHmin_ref.fitarg['n1'],LLHmin_ref.fitarg['n2'],LLHmin_ref.fitarg['n3'],LLHmin_ref.fitarg['nsig']) - LLH(LLHmin_fix.fitarg['n1'],LLHmin_fix.fitarg['n2'],LLHmin_fix.fitarg['n3'],param_up) )

            dTS = TS - TS_fix
            
        if nIterations >90:
            continue

        param_low=param_up/4.
        while(cc>0.):
            param_mean=(param_low+param_up)/2.
            LLHmin_fix=Minuit(LLH,
                 nsig=param_mean,fix_nsig = True,
                 n1=1.,n2=1.,n3=1.,
                 error_nsig=.1,error_n1=.1,error_n2=.1,error_n3=.1,
                 limit_nsig=(-10.,100.),limit_n1=(0.,10.),limit_n2=(0.,10.),limit_n3=(0.,10.),
                 errordef=.5,print_level=0)  
            LLHmin_fix.migrad()

            if param_mean <0.:
                TS_fix = 0.
            else:
                TS_fix = 2*(LLH(LLHmin_ref.fitarg['n1'],LLHmin_ref.fitarg['n2'],LLHmin_ref.fitarg['n3'],LLHmin_ref.fitarg['nsig']) - LLH(LLHmin_fix.fitarg['n1'],LLHmin_fix.fitarg['n2'],LLHmin_fix.fitarg['n3'],param_mean) )
            
            dTS = TS - TS_fix
            
            if(dTS<deltaTS):
                param_low=param_mean
                delta_param=(param_up-param_low)/(param_up)
                
                if((dTS>deltaTS-eps_TS) and (delta_param < eps_param)):
                    cc = 0
                    
            if(dTS>deltaTS):
                param_up=param_mean
                delta_param=(param_up-param_low)/(param_up)
                
                if((TS<deltaTS+eps_TS) and (delta_param < eps_param)):
                    cc=0

        ul = param_up
        
        TSdist = np.append(TSdist,TS)
        allFits = np.append(allFits,bestFit)
        upperLimit = np.append(upperLimit,ul)
        
        
p_median = np.percentile(upperLimit, 50)
p_95_low = np.percentile(upperLimit, 2.5)
p_95_high = np.percentile(upperLimit, 97.5)
p_68_low = np.percentile(upperLimit, 16.)
p_68_high = np.percentile(upperLimit, 84.)

dic_brazilian = {}
dic_brazilian['mass'] = mass
dic_brazilian['fits'] = allFits
dic_brazilian['TSdist'] = TSdist
dic_brazilian['nsig_upperLimits'] = upperLimit

if mode == 'annihilation':
    error_68_low = [p_68_low*10**-23]
    error_68_high = [p_68_high*10**-23]
    error_95_low = [p_95_low*10**-23]
    error_95_high = [p_95_high*10**-23]

    dic_brazilian['error_68_low'] = error_68_low
    dic_brazilian['error_68_high'] = error_68_high
    dic_brazilian['error_95_low'] = error_95_low
    dic_brazilian['error_95_high'] = error_95_high   
    dic_brazilian['median'] = p_median*10**-23


elif mode == 'decay':
    error_68_low = [1./p_68_low*10**28]
    error_68_high = [1./p_68_high*10**28]
    error_95_low = [1./p_95_low*10**28]
    error_95_high = [1./p_95_high*10**28]

    dic_brazilian['error_68_low'] = error_68_high
    dic_brazilian['error_68_high'] = error_68_low
    dic_brazilian['error_95_low'] = error_95_high
    dic_brazilian['error_95_high'] = error_95_low  
    dic_brazilian['median'] = [1./p_median*10**28]

print 'median sesitivity: ',p_median,dic_brazilian['median']

print 'save sensitivity to', out_file

np.save(out_file,dic_brazilian)

