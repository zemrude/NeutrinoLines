import numpy as np
from time import time
from scipy import interpolate, integrate
import imp

import sys, os

base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')

#####################################
## Astrophysical stuff
#####################################
def psi_f(RA,decl):
    return np.arccos(np.cos(np.pi/2.-(-29.*np.pi/180))*np.cos(np.pi/2.-decl)\
                      +np.sin(np.pi/2.-(-29.*np.pi/180))*np.sin(np.pi/2.-decl)*\
                       np.cos(RA-266.*np.pi/180))

convf = 3.08567758*10**21 # from kpc to cm

def NFW(r):
    rho_0 = 1.40 
    r_s = 16.1
    return rho_0 / ( r/r_s * (1+r/r_s)**2 ) * 1e7 * 1e-9 * 37.96
NFW = np.vectorize(NFW)

def Burkert(r):
    rho_0 = 4.13 
    r_s = 9.26
    return rho_0 / ( (1+r/r_s) * (1+(r/r_s)**2 ) ) * 1e7 * 1e-9 * 37.96
Burkert = np.vectorize(Burkert)

def Jfactor(profile, psi_in, exp):
    r_sol = 8.5
    s = np.linspace(0,50,10000)
    r = np.sqrt(r_sol**2+s**2-2*r_sol*s*np.cos(psi_in))
    if profile == 'NFW':
        rho = NFW(r)
    elif profile == 'Burkert':
        rho = Burkert(r)
    return integrate.trapz(np.power(rho,exp),s) * convf
Jfactor = np.vectorize(Jfactor)

psi_sample = np.linspace(0,np.pi,1000) 

J = {}
J['NFW'] = {}
J['NFW']['dec']=interpolate.UnivariateSpline(psi_sample,Jfactor('NFW', psi_sample, 1),k=1,s=0)
J['NFW']['ann']=interpolate.UnivariateSpline(psi_sample,Jfactor('NFW', psi_sample, 2),k=1,s=0)
J['Burkert'] = {}
J['Burkert']['dec']=interpolate.UnivariateSpline(psi_sample,Jfactor('Burkert', psi_sample, 1),k=1,s=0)
J['Burkert']['ann']=interpolate.UnivariateSpline(psi_sample,Jfactor('Burkert', psi_sample, 2),k=1,s=0)

#####################################
## DM spectra and flux
#####################################
from PPPC_spectra import dNdlog10x_dict
dNdlog10x=dNdlog10x_dict()

def dNdx(f,key,state,x,mDM): # Assuming ffbar production (no 1/2 for decay taken into account here!)
    return .5*(np.sign(1.-x)+1.)*dNdlog10x[key][state](f*mDM,np.log10(x))[0,:]/(np.log(10)*x)

def dNdE(f,key,state,E_true,mDM):
    return .5*(np.sign(f*mDM-E_true)+1)*dNdx(f,key,state,E_true/(f*mDM),mDM)/(f*mDM)

def phi_ann(profile,key,state,mDM,E_true,psi,sv=10**-23): # sv : <sigmav>
    return 1./(4*np.pi*mDM**2)*(sv/2.)*J[profile]['ann'](psi)*dNdE(1., key,state,E_true,mDM)
phi_ann=np.vectorize(phi_ann)

def phi_dec(profile,key,state,mDM,E_true,psi,tauDM=10**28):
    return 1./(4*np.pi*mDM*tauDM)*J[profile]['dec'](psi)*dNdE(0.5, key,state,E_true,mDM)
phi_dec=np.vectorize(phi_dec)

def phiDM_ann(mass,prtcl,profile,nutype_,nrg,psi):   # Flux of DM at the top of the atmosphere
    return 1/3.*(phi_ann(profile,prtcl,'nue',mass,nrg,psi) + phi_ann(profile,prtcl,'numu',mass,nrg,psi)+ phi_ann(profile,prtcl,'nutau',mass,nrg,psi))
phiDM_ann   = np.vectorize(phiDM_ann)

def phiDM_dec(mass,prtcl,profile,nutype_,nrg,psi):   # Flux of DM at the top of the atmosphere
    return 1/3.*(phi_dec(profile,prtcl,'nue',mass,nrg,psi) + phi_dec(profile,prtcl,'numu',mass,nrg,psi)+ phi_dec(profile,prtcl,'nutau',mass,nrg,psi))
phiDM_dec   = np.vectorize(phiDM_dec)


#####################################
## atmospheric flux
#####################################
from icecube import NewNuFlux
newflux = NewNuFlux.makeFlux('honda2006').getFlux

from icecube.nuCraft_icetray.nuCraft import NuCraft

def newflux2(t, e, cz):
    if e <= 10.:
        return 0.
    else:
        return newflux(t,e,cz) 
newflux2 = np.vectorize(newflux2)
    
def neutrinoOscillator(atmH=None, depth=None):
    # parameters from arxiv1611.01514
    theta23 = np.arcsin(np.sqrt(0.441))/np.pi*180.
    theta13 = np.arcsin(np.sqrt(0.02166))/np.pi*180.
    theta12 = np.arcsin(np.sqrt(0.306))/np.pi*180.;
    DM21   = 7.50e-5
    DM31   = 2.524e-3 + DM21
    
    if (atmH==None) and (depth==None):
        nC=NuCraft((1., DM21, DM31), [(1,2,theta12),(1,3,theta13,0),(2,3,theta23)])
    
    elif (atmH==None) and (depth!=None):
        nC=NuCraft((1., DM21, DM31), [(1,2,theta12),(1,3,theta13,0),(2,3,theta23)],detectorDepth=depth)#, atmHeight=0.)
    
    elif (atmH!=None) and (depth==None):
        nC=NuCraft((1., DM21, DM31), [(1,2,theta12),(1,3,theta13,0),(2,3,theta23)],atmHeight=atmH)#, atmHeight=0.)    
    
    else:
        nC=NuCraft((1., DM21, DM31), [(1,2,theta12),(1,3,theta13,0),(2,3,theta23)],detectorDepth=depth,atmHeight=atmH)    
    
    return nC
    
    
def osc_Atm_flux_weight(energy,zenith,nuType,nC,atm_mode=0):

    numPrec = 6e-4
   
    # assume nu_e at production
    prob_e   = nC.CalcWeights((nuType/np.abs(nuType)*12, energy, zenith), numPrec=numPrec,atmMode=atm_mode)
    flux_e   = newflux2(nuType/np.abs(nuType)*12, energy, np.cos(zenith))

    # assume nu_mu at production
    prob_mu  = nC.CalcWeights((nuType/np.abs(nuType)*14, energy, zenith), numPrec=numPrec,atmMode=atm_mode)
    flux_mu  = newflux2(nuType/np.abs(nuType)*14, energy, np.cos(zenith))

    # assume nu_tau at production
    prob_tau = nC.CalcWeights((nuType/np.abs(nuType)*16, energy, zenith), numPrec=numPrec,atmMode=atm_mode)
    flux_tau   = newflux2(nuType/np.abs(nuType)*16, energy, np.cos(zenith))
    
    # case 1: nu_e sample:
    if np.abs(nuType[0])==12:
        prob_e_e   = np.array(prob_e)[:,0]   # nu_e at production, reconstructed as nu_e
        prob_mu_e  = np.array(prob_mu)[:,0]  # nu_mu at production, reconstructed as nu_e
        prob_tau_e = np.array(prob_tau)[:,0] # nu_tau at production, reconstructed as nu_e

        weight = flux_e*prob_e_e + flux_mu*prob_mu_e + flux_tau*prob_tau_e

        return weight
    
    
    # case 2: nu_mu sample:
    if np.abs(nuType[0])==14:
        prob_e_mu   = np.array(prob_e)[:,1]   # nu_e at production, reconstructed as nu_mu
        prob_mu_mu  = np.array(prob_mu)[:,1]  # nu_mu at production, reconstructed as nu_mu
        prob_tau_mu = np.array(prob_tau)[:,1] # nu_tau at production, reconstructed as nu_mu

        weight = flux_e*prob_e_mu + flux_mu*prob_mu_mu + flux_tau*prob_tau_mu

        return weight
    
    
    # case 3: nu_tau sample:
    if np.abs(nuType[0])==16:
        prob_e_tau   = np.array(prob_e)[:,2]   # nu_e at production, reconstructed as nu_tau
        prob_mu_tau  = np.array(prob_mu)[:,2]  # nu_mu at production, reconstructed as nu_tau
        prob_tau_tau = np.array(prob_tau)[:,2] # nu_tau at production, reconstructed as nu_tau

        weight = flux_e*prob_e_tau + flux_mu*prob_mu_tau + flux_tau*prob_tau_tau

        return weight
