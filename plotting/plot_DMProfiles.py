import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pickle
import sys, os

from scipy import integrate

base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'plotting')

import style
style.SetFigureStyle()
style.increaseAxisText(6)
style.increaseLegendText(6)

plt.rcParams.update({'figure.dpi': 300.})


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

convf = 3.08567758*10**21 # from kpc to cm

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



####################################

xx = np.linspace(0.01,50,1000)


fig,ax = plt.subplots(1,1,figsize=(10,8))
ax.plot(xx,NFW(xx),'r-',lw=3, label='NFW')
ax.plot(xx,Burkert(xx),'b-',lw=3, label='Burkert')

ax.set_yscale('log')
ax.set_xlabel('r (kpc)')
ax.set_ylabel(r'$\rho(r)$ (GeV/cm$^3$)')

ax. set_title('Profile comparisons')
leg = ax.legend(frameon = 1, fancybox=False, loc='upper right',labelspacing=0.2)
fig.savefig('plots/DMProfileComparison.png')
fig.savefig('plots/DMProfileComparison.pdf')


####################################


psi = np.linspace(0,np.pi,100)

fig,(ax,ax2) = plt.subplots(1,2,figsize=(20,8))

ax.plot(psi,Jfactor('NFW', psi, 2),'r-',lw=3, label='NFW')
ax.plot(psi,Jfactor('Burkert', psi, 2),'b-',lw=3, label='Burkert')
ax.set_yscale('log')
ax.set_xlabel('$\Psi$')
ax.set_ylabel(r'$J_{\mathrm{\Psi}}^{\mathrm{ann}}$ (GeV$^2$ cm$^{-5}$ sr$^{-1}$)')
leg = ax.legend(frameon = 1, fancybox=False, loc='upper right',labelspacing=0.2)

ax2.plot(psi,Jfactor('NFW', psi, 1),'r-',lw=3, label='NFW')
ax2.plot(psi,Jfactor('Burkert', psi, 1),'b-',lw=3, label='Burkert')
ax2.set_yscale('log')
ax2.set_xlabel('$\Psi$')
ax2.set_ylabel(r'$J_{\mathrm{\Psi}}^{\mathrm{dec}}$ (GeV cm$^{-2}$ sr$^{-1}$)')

leg = ax2.legend(frameon = 1, fancybox=False, loc='upper right',labelspacing=0.2)
fig.savefig('plots/DMJfactorComparison.png')
fig.savefig('plots/DMJfactorComparison.pdf')

####################################
