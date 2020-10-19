import numpy as np
from time import time
from scipy import interpolate, integrate
import sys, os
  
from HondaFlux import Honda_Flux

try:
    base_path = os.environ['ANALYSIS_BASE_PATH']
except Exception as e:
    raise Exception("No enviromental variables loaded, try doing `import env` before")
    

sys.path.append(base_path+'python')

#####################################
## Astrophysical stuff
#####################################
from utils import psi_f

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

