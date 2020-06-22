from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pickle
import sys, os

if sys.version_info >= (3, 5):
    from math import gcd
else:
    from fractions import gcd
    
from iminuit import Minuit
import scipy.special as sps

from scipy.stats import chisquare, distributions

from optparse import OptionParser


###############################

def model_array(n1,n2,n3):
    return n1*tab_events_corsika + n2*tab_events_atm + n3*tab_events_astro

def model_array_quad(n1,n2,n3):
    return n1**2*tabquad_events_corsika + n2**2*tabquad_events_atm + n3**2*tabquad_events_astro

def LLH(n1,n2,n3):
         
    k      = observation
    mu     = model_array(n1,n2,n3)
    mu2    = np.power(mu,2)
    sigma2 = model_array_quad(n1,n2,n3)
  
    bins_to_use = (k>0.)&(sigma2>0.)
    
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
    
    return -np.sum(values)
    
    
def _compute_prob_inside_method(m, n, g, h):
    """
    Count the proportion of paths that stay strictly inside two diagonal lines.
    Parameters
    ----------
    m : integer
        m > 0
    n : integer
        n > 0
    g : integer
        g is greatest common divisor of m and n
    h : integer
        0 <= h <= lcm(m,n)
    Returns
    -------
    p : float
        The proportion of paths that stay inside the two lines.
    Count the integer lattice paths from (0, 0) to (m, n) which satisfy
    |x/m - y/n| < h / lcm(m, n).
    The paths make steps of size +1 in either positive x or positive y directions.
    We generally follow Hodges' treatment of Drion/Gnedenko/Korolyuk.
    Hodges, J.L. Jr.,
    "The Significance Probability of the Smirnov Two-Sample Test,"
    Arkiv fiur Matematik, 3, No. 43 (1958), 469-86.
    
    """
    # Probability is symmetrical in m, n.  Computation below uses m >= n.
    if m < n:
        m, n = n, m
    mg = m // g
    ng = n // g

    # Count the integer lattice paths from (0, 0) to (m, n) which satisfy
    # |nx/g - my/g| < h.
    # Compute matrix A such that:
    #  A(x, 0) = A(0, y) = 1
    #  A(x, y) = A(x, y-1) + A(x-1, y), for x,y>=1, except that
    #  A(x, y) = 0 if |x/m - y/n|>= h
    # Probability is A(m, n)/binom(m+n, n)
    # Optimizations exist for m==n, m==n*p.
    # Only need to preserve a single column of A, and only a sliding window of it.
    # minj keeps track of the slide.
    minj, maxj = 0, min(int(np.ceil(h / mg)), n + 1)
    curlen = maxj - minj
    # Make a vector long enough to hold maximum window needed.
    lenA = min(2 * maxj + 2, n + 1)
    # This is an integer calculation, but the entries are essentially
    # binomial coefficients, hence grow quickly.
    # Scaling after each column is computed avoids dividing by a
    # large binomial coefficent at the end. Instead it is incorporated
    # one factor at a time during the computation.
    dtype = np.float64
    A = np.zeros(lenA, dtype=dtype)
    # Initialize the first column
    A[minj:maxj] = 1
    for i in range(1, m + 1):
        # Generate the next column.
        # First calculate the sliding window
        lastminj, lastmaxj, lastlen = minj, maxj, curlen
        minj = max(int(np.floor((ng * i - h) / mg)) + 1, 0)
        minj = min(minj, n)
        maxj = min(int(np.ceil((ng * i + h) / mg)), n + 1)
        if maxj <= minj:
            return 0
        # Now fill in the values
        A[0:maxj - minj] = np.cumsum(A[minj - lastminj:maxj - lastminj])
        curlen = maxj - minj
        if lastlen > curlen:
            # Set some carried-over elements to 0
            A[maxj - minj:maxj - minj + (lastlen - curlen)] = 0
        # Peel off one term from each of top and bottom of the binomial coefficient.
        scaling_factor = i * 1.0 / (n + i)
        A *= scaling_factor
    return A[maxj - minj - 1]

def ks_w2(data1, data2, wei1, wei2):
    ix1 = np.argsort(data1)
    ix2 = np.argsort(data2)
    data1 = data1[ix1]
    data2 = data2[ix2]
    wei1 = wei1[ix1]
    wei2 = wei2[ix2]
    data = np.concatenate([data1, data2])
    cwei1 = np.hstack([0, np.cumsum(wei1)/sum(wei1)])
    cwei2 = np.hstack([0, np.cumsum(wei2)/sum(wei2)])
    cdf1we = cwei1[[np.searchsorted(data1, data, side='right')]]
    cdf2we = cwei2[[np.searchsorted(data2, data, side='right')]]
    
    d = np.max(np.abs(cdf1we - cdf2we))
    
    n1 = int(np.sum(wei1))
    n2 = int(np.sum(wei2))
    
    
    g = gcd(n1, n2)
    n1g = n1 // g
    n2g = n2 // g
    prob = -np.inf

    lcm = (n1 // g) * n2
    h = int(np.round(d * lcm))
    d = h * 1.0 / lcm
    if h == 0:
        prob = 1.0
    else:
        prob = 1 - _compute_prob_inside_method(n1, n2, g, h)
    
    return d , prob



######################################

base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'plotting')

import style
style.SetFigureStyle()
latex = style.latex
style.increaseAxisText(6)

plt.rcParams.update({'figure.dpi': 200.})

burnsampleTime = 1116125.821572

for variable in ['energy_rec','psi_rec']:
    
    
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(20,8))
    
    cutValueString = {'LE': 'LEBDT0.15_HEBDT0.2' , 'HE': 'LEBDT-1.0_HEBDT0.3'}
    for s in ['LE','HE']:
       
        infile = base_path+'controlDistributions/Histograms_'+cutValueString[s]+'_PSI90degrees.pkl' 
        data = pickle.load(open(infile,'r'))

        # load histograms
        tab_events_atm    = data['histograms']['atm'][variable][0]*burnsampleTime
        tab_events_astro  = data['histograms']['astro'][variable][0]*burnsampleTime
        tab_events_corsika  = data['histograms']['corsika'][variable][0]*burnsampleTime
        tabquad_events_atm   = data['histograms']['corsika'][variable+'_quad'][0]*burnsampleTime*burnsampleTime
        tabquad_events_astro = data['histograms']['corsika'][variable+'_quad'][0]*burnsampleTime*burnsampleTime
        tabquad_events_corsika = data['histograms']['corsika'][variable+'_quad'][0]*burnsampleTime*burnsampleTime

        observation = data['histograms']['Burnsample'][variable][0]

        # fit background only
        LLHmin=Minuit(LLH,
                 n1=1.,n2=1.,n3=1.,
                 error_n1=.1,error_n2=.1,error_n3=.1,
                 limit_n1=(0.,2.),limit_n2=(0.,2.),limit_n3=(0.,2.),
                 errordef=.5,print_level=0)  
        LLHmin.migrad()

        bestFit = {}
        bestFit['n1']=LLHmin.fitarg['n1']
        bestFit['n2']=LLHmin.fitarg['n2']
        bestFit['n3']=LLHmin.fitarg['n3']
        bestFit['LLH']=LLH(bestFit['n1'],bestFit['n2'],bestFit['n3'])
        
        #load unbinned data
        d_burnsample = data['unbinned_data']['Burnsample']['energy_rec']
        weight_burnsample = np.array([1]*len(data['unbinned_data']['Burnsample']['energy_rec']))

        d_corskia = data['unbinned_data']['corsika']['energy_rec']
        d_atm = data['unbinned_data']['atm']['energy_rec']
        d_astro = data['unbinned_data']['astro']['energy_rec']

        weight_corskia = data['unbinned_data']['corsika']['weight']
        weight_atm = data['unbinned_data']['atm']['weight']
        weight_astro = data['unbinned_data']['astro']['weight']

        bkg = np.append(d_corskia, d_atm)
        bkg = np.append(bkg, d_astro)
        bkg_weight = np.append(bestFit['n1']*weight_corskia*burnsampleTime, bestFit['n2']*weight_atm*burnsampleTime)
        bkg_weight = np.append(bkg_weight, bestFit['n3']*weight_astro*burnsampleTime)

        # calculate p-values
        # 1) weighted two-sample KS
        ks_test = ks_w2(d_burnsample, bkg, weight_burnsample, bkg_weight)

        # 2) chi^2
        expectation = tab_events_corsika+tab_events_atm+tab_events_astro
        chi2_test = chisquare(observation[expectation>0.], expectation[expectation>0.])
        
        # plot
        h_burnsample = data['histograms']['Burnsample'][variable][0]
        bins = data['histograms']['Burnsample'][variable][1]
        binCenters = bins[:-1]+np.diff(bins)/2.
        if variable == 'psi_rec':
            binCenters = binCenters *180/np.pi
            
        h_atm = data['histograms']['atm'][variable][0]
        h_corsika = data['histograms']['corsika'][variable][0]
        h_astro = data['histograms']['astro'][variable][0]
        h_bkg = bestFit['n2']*h_atm + bestFit['n1']*h_corsika + bestFit['n3']*h_astro

        if s == 'LE':
            ax = ax1
        else:
            ax = ax2
            
        ax.errorbar(binCenters, h_burnsample/burnsampleTime, np.sqrt(h_burnsample)/burnsampleTime, xerr=None, color='k', marker='o',markersize=2, lw=2, zorder=10, ls=None, label='Burnsample ({} events)'.format(len(d_burnsample)))

        ax.step(binCenters,bestFit['n1']*h_corsika,where='mid',color = 'orange', label='Corsika', alpha = 0.5, lw=2)
        ax.step(binCenters,bestFit['n2']*h_atm,where='mid',color = 'g', label=r'Atmosperic $\nu$', alpha = 0.5, lw=2)
        ax.step(binCenters,bestFit['n3']*h_astro,where='mid',color = 'b', label=r'Astrophysical $\nu$', alpha = 0.5, lw=2)
        ax.step(binCenters,h_bkg,where='mid',color = 'r', label='Combined background', lw=3)

        ax.plot(0,0,lw=0,label='KS p-value: {:.3f}'.format(ks_test[1]),zorder=100)
        ax.plot(0,0,lw=0,label=r'Pearson $\chi^2$ p-value: {:.3f}'.format(chi2_test[1]),zorder=200)
    
        handles, labels = ax.get_legend_handles_labels()
        
        if variable == 'energy_rec':
            ax.set_xlabel(r"E$_{\mathrm{reco}}$ (GeV)")
            ax.set_xscale('log')
            ax.set_xlim(1e1,1e5)
            leg = ax.legend(handles[::-1], labels[::-1],frameon = 1, fancybox=False, loc='upper right')

        if variable == 'psi_rec':
            ax.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (degrees)")
            ax.set_xlim(0,180)
            leg = ax.legend(handles[::-1], labels[::-1],frameon = 1, fancybox=False, loc='lower left')
            
        ax.set_ylim(1e-8,1e-4)
        ax.set_ylabel(r"rate (Hz)")
        ax.set_yscale('log')
        ax.set_title(s+r' Sample ($\Psi_{\mathrm{reco}}$>90$^\circ$)')
    

    fig.savefig('plots/PreUnblinding_ControlDistribution_'+variable+'.png',bbox='tight')
    fig.savefig('plots/PreUnblinding_ControlDistribution_'+variable+'.pdf',bbox='tight')
