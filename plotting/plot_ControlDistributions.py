import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches

from matplotlib import gridspec

import pickle
import sys, os

import math

from scipy.stats import chisquare, distributions, chi2
from scipy.interpolate import interp1d

from iminuit import Minuit
import scipy.special as sps

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--psi", default='90',
                  dest="PSICUT", help="Cut on Psi")
(options,args) = parser.parse_args()
        

psiCut = float(options.PSICUT)/180.*np.pi

###############################

samples = ['LE','HE']

base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'plotting')

import style
style.SetFigureStyle()
latex = style.latex
style.increaseAxisText(6)

plt.rcParams.update({'figure.dpi': 200.})

livetime = {}
livetime['Burnsample'] = 1116125.821572
livetime['Data'] = 28272940. + 30674072. + 31511810.5 + 31150852. + 30059465.  

bins_vars={
    'energy_rec':np.logspace(1,5,24+1),
    'psi_rec':np.linspace(0.,np.pi,18+1),
}

cutValueString = {'LE': 'LEBDT0.15_HEBDT0.2' , 'HE': 'LEBDT-1.0_HEBDT0.3'}

systematics = ['nominal','DomEffUp','DomEffDown','Ice_HoleIce100_ScatAbs-7','Ice_HoleIce100_Scat+10','Ice_HoleIce100_Abs+10',
              'Ice_HoleIce30_ScatAbs-7','Ice_HoleIce30_Scat+10','Ice_HoleIce30_Abs+10','nominalGammaUp','nominalGammaDown']

###############################

def LLH(n1,n2,n3,syst):
    
    s = systematics[int(syst)]
    
    k      = np.append(tab_events_data_energy,tab_events_data_psi)
    mu     = np.append(
        n1*tab_events_corsika_energy + n2*tab_events_atm_energy[s] + n3*tab_events_astro_energy[s],
        n1*tab_events_corsika_psi + n2*tab_events_atm_psi[s] + n3*tab_events_astro_psi[s],
    )
    mu2    = np.power(mu,2)
    sigma2 = np.append(
        n1**2*tabquad_events_corsika_energy + n2**2*tabquad_events_atm_energy[s] + n3**2*tabquad_events_astro_energy[s],
        n1**2*tabquad_events_corsika_psi + n2**2*tabquad_events_atm_psi[s] + n3**2*tabquad_events_astro_psi[s]
    )
  
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


def cdf_interpolation(data, weight):
    
    ix = np.argsort(data)
    data = data[ix]
    wei = np.array(weight[ix],dtype=float)
    
    cwei = np.cumsum(wei)/np.sum(wei)
   
    cdf = interp1d(data, cwei, kind='linear', bounds_error=False, fill_value=(0,1))
    
    return cdf


def ks_w2(data1, data2, wei1, wei2):

    cdf1 = cdf_interpolation(data1, wei1)
    cdf2 = cdf_interpolation(data2, wei2)
    
    xx = np.linspace(np.min(data1), np.max(data1))
    
    cdf1we = cdf1(xx)
    cdf2we = cdf2(xx)
    cddiffs = cdf1we - cdf2we
    
    minS = np.clip(-np.min(cddiffs), 0, 1)  # Ensure sign of minS is not negative.
    maxS = np.max(cddiffs)
    d = max(minS, maxS)
    
    n1 = int(np.sum(wei1)**2/np.sum(wei1**2))
    n2 = int(np.sum(wei2)**2/np.sum(wei2**2))
    
    n_eff = (n1 * n2 / (n1 + n2))
        
    return d , sps.kolmogorov(np.sqrt(n_eff)*d)
       

for Datatype in ['Data']:#,'Burnsample']:

    for sample in samples:
        infile = base_path+'/controlDistributions/ControlDistributions_unbinned_syst_'+cutValueString[sample]+'.pkl' 

        with open(infile, 'rb') as f:
            data = pickle.load(f)#, encoding="latin1") 

            d_corsika_energy = data['unbinned_data']['corsika']['energy_rec']
            d_corsika_psi = data['unbinned_data']['corsika']['psi_rec']

            d_corsika_energy = d_corsika_energy[d_corsika_psi>psiCut]
            d_corsika_weight_energy = data['unbinned_data']['corsika']['weight'][d_corsika_psi>psiCut]
            d_corsika_weight_psi = data['unbinned_data']['corsika']['weight']

            d_atm_energy = {}
            d_atm_psi = {}
            d_atm_weight = {}

            d_astro_energy = {}
            d_astro_psi = {}
            d_astro_weight = {}

            for s in systematics:

                d_atm_energy[s] = data['unbinned_data']['atm'][s]['energy_rec']
                d_atm_psi[s] = data['unbinned_data']['atm'][s]['psi_rec']
                d_atm_weight[s] = data['unbinned_data']['atm'][s]['weight']

                d_astro_energy[s] = data['unbinned_data']['astro'][s]['energy_rec']
                d_astro_psi[s] = data['unbinned_data']['astro'][s]['psi_rec']
                d_astro_weight[s] = data['unbinned_data']['astro'][s]['weight']

            d_data_energy = data['unbinned_data'][Datatype]['energy_rec']
            d_data_psi = data['unbinned_data'][Datatype]['psi_rec']

            d_data_energy = d_data_energy[d_data_psi>psiCut]
            d_data_psi = d_data_psi[d_data_psi>psiCut]
            
            d_data_weight = np.array([[1]*len(d_data_psi)])


            ##### do histograms

            tab_events_corsika_energy = np.histogram(d_corsika_energy.flatten(),bins=bins_vars['energy_rec'],weights = d_corsika_weight_energy*livetime[Datatype])[0]
            tab_events_corsika_psi = np.histogram(d_corsika_psi[d_corsika_psi>psiCut],bins=bins_vars['psi_rec'],weights = d_corsika_weight_psi[d_corsika_psi.flatten()>psiCut]*livetime[Datatype])[0]


            tabquad_events_corsika_energy = np.histogram(d_corsika_energy.flatten(),bins=bins_vars['energy_rec'],weights = (d_corsika_weight_energy*livetime[Datatype])**2)[0]
            tabquad_events_corsika_psi = np.histogram(d_corsika_psi[d_corsika_psi>psiCut],bins=bins_vars['psi_rec'],weights = (d_corsika_weight_psi[d_corsika_psi.flatten()>psiCut]*livetime[Datatype])**2)[0]

            tab_events_atm_energy = {}
            tab_events_atm_psi = {}
            tab_events_astro_energy = {}
            tab_events_astro_psi = {}
            tabquad_events_atm_energy = {}
            tabquad_events_atm_psi = {}
            tabquad_events_astro_energy = {}
            tabquad_events_astro_psi = {}

            for s in systematics:

                tab_events_atm_energy[s] = np.histogram(d_atm_energy[s][(d_atm_psi[s])>psiCut],bins=bins_vars['energy_rec'],weights = d_atm_weight[s][(d_atm_psi[s])>psiCut]*livetime[Datatype])[0]
                tab_events_atm_psi[s] = np.histogram(d_atm_psi[s][(d_atm_psi[s])>psiCut],bins=bins_vars['psi_rec'],weights = d_atm_weight[s][(d_atm_psi[s])>psiCut]*livetime[Datatype])[0]

                tab_events_astro_energy[s] = np.histogram(d_astro_energy[s][d_astro_psi[s]>psiCut],bins=bins_vars['energy_rec'],weights = d_astro_weight[s][d_astro_psi[s]>psiCut]*livetime[Datatype])[0]
                tab_events_astro_psi[s] = np.histogram(d_astro_psi[s][d_astro_psi[s]>psiCut],bins=bins_vars['psi_rec'],weights = d_astro_weight[s][d_astro_psi[s]>psiCut]*livetime[Datatype])[0]

                tabquad_events_atm_energy[s] = np.histogram(d_atm_energy[s][d_atm_psi[s]>psiCut],bins=bins_vars['energy_rec'],weights = np.power(d_atm_weight[s][d_atm_psi[s]>psiCut]*livetime[Datatype],2))[0]
                tabquad_events_atm_psi[s] = np.histogram(d_atm_psi[s][d_atm_psi[s]>psiCut],bins=bins_vars['psi_rec'],weights = d_atm_weight[s][d_atm_psi[s]>psiCut]*livetime[Datatype]*d_atm_weight[s][d_atm_psi[s]>psiCut]*livetime[Datatype])[0]

                tabquad_events_astro_energy[s] = np.histogram(d_astro_energy[s][d_astro_psi[s]>psiCut],bins=bins_vars['energy_rec'],weights = d_astro_weight[s][d_astro_psi[s]>psiCut]*livetime[Datatype]*d_astro_weight[s][d_astro_psi[s]>psiCut]*livetime[Datatype])[0]
                tabquad_events_astro_psi[s] = np.histogram(d_astro_psi[s][d_astro_psi[s]>psiCut],bins=bins_vars['psi_rec'],weights = d_astro_weight[s][d_astro_psi[s]>psiCut]*livetime[Datatype]*d_astro_weight[s][d_astro_psi[s]>psiCut]*livetime[Datatype])[0]


            tab_events_data_energy, bins_energy = np.histogram(d_data_energy[d_data_psi>psiCut],bins=bins_vars['energy_rec'])
            tab_events_data_psi, bins_psi = np.histogram(d_data_psi[d_data_psi>psiCut],bins=bins_vars['psi_rec'])

            binCenters_energy = bins_energy[:-1]+np.diff(bins_energy)/2.
            binCenters_psi = bins_psi[:-1]+np.diff(bins_psi)/2.
            binCenters_psi = binCenters_psi*180/np.pi


            # fit background only

            bestFit = {}

            for iSyst in range(len(systematics)):

                bestFit[systematics[iSyst]] = {}

                LLHmin=Minuit(LLH,
                                     n1=1.,n2=1.,n3=1.,
                                     error_n1=.1,error_n2=.1,error_n3=.1,
                                     syst = iSyst, fix_syst = True,
                                     limit_n1=(0.,2.),limit_n2=(0.,2.),limit_n3=(0.,2.),
                                     errordef=.5,print_level=0)
                LLHmin.migrad()

                bestFit[systematics[iSyst]]['n1']=LLHmin.fitarg['n1']
                bestFit[systematics[iSyst]]['n2']=LLHmin.fitarg['n2']
                bestFit[systematics[iSyst]]['n3']=LLHmin.fitarg['n3']
                bestFit[systematics[iSyst]]['LLH']=LLH(bestFit[systematics[iSyst]]['n1'],bestFit[systematics[iSyst]]['n2'],bestFit[systematics[iSyst]]['n3'],iSyst)


            ###### 
            d_bkg_energy = {}
            d_bkg_weight_energy = {}
            d_bkg_psi = {}
            d_bkg_weight_psi = {}

            for s in systematics:
                d_bkg_energy[s] = np.append(d_corsika_energy, d_atm_energy[s][d_atm_psi[s]>psiCut])
                d_bkg_weight_energy[s] = np.append(bestFit[s]['n1']*d_corsika_weight_energy*livetime[Datatype], bestFit[s]['n2']*d_atm_weight[s][d_atm_psi[s]>psiCut]*livetime[Datatype]+bestFit[s]['n3']*d_astro_weight[s][d_astro_psi[s]>psiCut]*livetime[Datatype])

                d_bkg_psi[s] = np.append(d_corsika_psi[d_corsika_psi>psiCut], d_atm_psi[s][d_atm_psi[s]>psiCut])
                d_bkg_weight_psi[s] = np.append(bestFit[s]['n1']*d_corsika_weight_psi[d_corsika_psi.flatten()>psiCut]*livetime[Datatype], bestFit[s]['n2']*d_atm_weight[s][d_atm_psi[s]>psiCut]*livetime[Datatype]+bestFit[s]['n3']*d_astro_weight[s][d_astro_psi[s]>psiCut]*livetime[Datatype])

            #################
            h_bkg_energy = {}
            h_bkg_psi = {}
            for s in systematics:
                h_bkg_energy[s] = bestFit[s]['n2']*tab_events_atm_energy[s]/livetime[Datatype] + bestFit[s]['n1']*tab_events_corsika_energy/livetime[Datatype] + bestFit[s]['n3']*tab_events_astro_energy[s]/livetime[Datatype]

                h_bkg_psi[s] = bestFit[s]['n2']*tab_events_atm_psi[s]/livetime[Datatype] + bestFit[s]['n1']*tab_events_corsika_psi/livetime[Datatype] + bestFit[s]['n3']*tab_events_astro_psi[s]/livetime[Datatype]

                
            h_bkg_energy_stat = np.sqrt(
                bestFit['nominal']['n1']**2*tabquad_events_corsika_energy
                + bestFit['nominal']['n2']**2*tabquad_events_atm_energy['nominal']
                + bestFit['nominal']['n3']**2*tabquad_events_astro_energy['nominal']
            )/livetime[Datatype]
                
                
            h_bkg_energy_syst_min = np.array([])    
            h_bkg_energy_syst_max = np.array([]) 
            for i in range(len(h_bkg_energy['nominal'])):
                h_bkg_energy_syst_min = np.append(h_bkg_energy_syst_min, 
                                  h_bkg_energy['nominal'][i] - np.sqrt(
                                      np.abs(np.min([h_bkg_energy[s][i] for s in systematics])-h_bkg_energy['nominal'][i])**2
                                  + h_bkg_energy_stat[i]**2
                                  ))
                h_bkg_energy_syst_max = np.append(h_bkg_energy_syst_max,
                                  h_bkg_energy['nominal'][i] + np.sqrt(
                                      np.abs(np.max([h_bkg_energy[s][i] for s in systematics])-h_bkg_energy['nominal'][i])**2
                                  + h_bkg_energy_stat[i]**2 
                                  ))

            h_bkg_psi_stat = np.sqrt(
                bestFit['nominal']['n1']**2*tabquad_events_corsika_psi
                + bestFit['nominal']['n2']**2*tabquad_events_atm_psi['nominal']
                + bestFit['nominal']['n3']**2*tabquad_events_astro_psi['nominal']
            )/livetime[Datatype]
                                                  
            h_bkg_psi_syst_min = np.array([])    
            h_bkg_psi_syst_max = np.array([]) 
            for i in range(len(h_bkg_psi['nominal'])):
                h_bkg_psi_syst_min = np.append(h_bkg_psi_syst_min,
                                   h_bkg_psi['nominal'][i] - np.sqrt(
                                       np.abs(np.min([h_bkg_psi[s][i] for s in systematics])-h_bkg_psi['nominal'][i])**2
                                        + h_bkg_psi_stat[i]**2 
                                           ))
                h_bkg_psi_syst_max = np.append(h_bkg_psi_syst_max,
                                    h_bkg_psi['nominal'][i] + np.sqrt(
                                           np.abs(np.max([h_bkg_psi[s][i] for s in systematics])-h_bkg_psi['nominal'][i])**2
                                               + h_bkg_psi_stat[i]**2 
                                           ))

  
            # calculate p-values
            # 1) weighted two-sample KS

            KS_pValues_energy = []
            KS_pValues_psi = []
            
            for s in systematics:                
                KS_pValues_energy.append( ks_w2(np.log10(d_data_energy), np.log10(d_bkg_energy[s]), d_data_weight.flatten(), d_bkg_weight_energy[s])[1])
                KS_pValues_psi.append(ks_w2(d_data_psi, d_bkg_psi[s], d_data_weight.flatten(), d_bkg_weight_psi[s])[1])
            
            # 2) chi^2
            expectation_energy = bestFit['nominal']['n1']*tab_events_corsika_energy+bestFit['nominal']['n2']*tab_events_atm_energy['nominal']+bestFit['nominal']['n3']*tab_events_astro_energy['nominal']
            expectation_energy = np.rint(expectation_energy)

            expectation_psi = bestFit['nominal']['n1']*tab_events_corsika_psi+bestFit['nominal']['n2']*tab_events_atm_psi['nominal']+bestFit['nominal']['n3']*tab_events_astro_psi['nominal']
            expectation_psi = np.rint(expectation_psi)

            chi2_syst_energy = np.sum(

                       ( (tab_events_data_energy[expectation_energy>0.] - expectation_energy[expectation_energy>0.]) **2 )

                        / 

                       (
                          (((h_bkg_energy_syst_max[expectation_energy>0.] - h_bkg_energy_syst_min[expectation_energy>0.])
                           *livetime[Datatype] ) **2  )
                       )
                    )
            p_syst_energy = distributions.chi2.sf(chi2_syst_energy, len(expectation_energy[expectation_energy>0.]) - 1)

            chi2_syst_psi = np.sum(
                        ( (tab_events_data_psi[expectation_psi>0.] - expectation_psi[expectation_psi>0.]) **2)
                    /

                    (
                       (((h_bkg_psi_syst_max[expectation_psi>0.] - h_bkg_psi_syst_min[expectation_psi>0.])
                         *livetime[Datatype] ) **2 )
                    )
                 )
            p_syst_psi = distributions.chi2.sf(chi2_syst_psi, len(expectation_psi[expectation_psi>0.]) - 1)



            #################
            #################
            #################

            #fig,(ax1,ax2) = plt.subplots(1,2,figsize=(20,8))

            fig = plt.figure()
            
            fig = plt.figure(figsize=(20,12)) 
            gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1]) 

            ax1_ratio = plt.subplot(gs[2])
            ax2_ratio = plt.subplot(gs[3])

            ax1 = plt.subplot(gs[0], sharex = ax1_ratio)
            ax2 = plt.subplot(gs[1], sharex = ax2_ratio)
            
            fig.subplots_adjust(hspace = 0.1)

            
            #################
            ax1.errorbar(binCenters_energy, tab_events_data_energy/livetime[Datatype], np.sqrt(tab_events_data_energy)/livetime[Datatype], xerr=None, color='k', marker='o',markersize=4, lw=2, zorder=10, ls='none', label=Datatype+' ({} events)'.format(len(d_data_energy)))

            ax1.step(binCenters_energy,bestFit['nominal']['n3']*tab_events_astro_energy['nominal']/livetime[Datatype],where='mid',color = 'b', label=r'Astrophysical $\nu$', alpha = 0.5, lw=2)
            ax1.step(binCenters_energy,bestFit['nominal']['n2']*tab_events_atm_energy['nominal']/livetime[Datatype],where='mid',color = 'g', label=r'Atmosperic $\nu$', alpha = 0.5, lw=2)
            ax1.step(binCenters_energy,bestFit['nominal']['n1']*tab_events_corsika_energy/livetime[Datatype],where='mid',color = 'orange', label='Corsika', alpha = 0.5, lw=2)

            ax1.fill_between(binCenters_energy, h_bkg_energy_syst_min, h_bkg_energy_syst_max, step='mid',color='r',alpha=0.5,lw=0)
            p2 = patches.Patch(color='r', alpha=0.5, linewidth=0)

            p1,= ax1.step(binCenters_energy,h_bkg_energy['nominal'],where='mid',color = 'r', lw=3)

            ax1.plot(0,0,lw=0,label='KS p-value: {:.3f}'.format(max(KS_pValues_energy)),zorder=100)
            ax1.plot(0,0,lw=0,label=r'Pearson $\chi^2$ p-value: {:.3f}'.format(p_syst_energy),zorder=200)

            handles, labels = ax1.get_legend_handles_labels()
            handles.append((p2,p1))
            labels.append(r'Combined background')
            
            order = [5,4,3,6,2,1,0]
            leg = ax1.legend([handles[idx] for idx in order], [labels[idx] for idx in order], frameon = 1, fancybox=False, loc='upper right')

            #ax1.set_xlabel(r"E$_{\mathrm{reco}}$ (GeV)")
            #ax1.set_xscale('log')
            #ax1.set_xlim(1e1,1e5)
            ax1.set_ylim(1e-9,1e-4)
            ax1.set_ylabel(r"rate (Hz)")
            ax1.set_yscale('log')
            plt.setp(ax1.get_xticklabels(), visible=False)

            
            ax1_ratio.errorbar(binCenters_energy, tab_events_data_energy/h_bkg_energy['nominal']/livetime[Datatype], np.sqrt(tab_events_data_energy)/h_bkg_energy['nominal']/livetime[Datatype], xerr=None, color='k', marker='o',markersize=4, lw=2, zorder=10, ls='none', label=Datatype+' ({} events)'.format(len(d_data_energy)))
            
            ax1_ratio.fill_between(binCenters_energy, h_bkg_energy_syst_min/h_bkg_energy['nominal'], h_bkg_energy_syst_max/h_bkg_energy['nominal'], step='mid',color='r',alpha=0.5,lw=0)

            ax1_ratio.step(binCenters_energy,h_bkg_energy['nominal']/h_bkg_energy['nominal'],where='mid',color = 'r', lw=3)
            
            
            ax1_ratio.set_xlabel(r"E$_{\mathrm{reco}}$ (GeV)")
            ax1_ratio.set_xscale('log')
            ax1_ratio.set_xlim(1e1,1e5)
            ax1_ratio.set_ylabel(r"ratio")
            ax1_ratio.set_ylim(0.5,1.5)
            #ax1.set_yscale('log')
            #################
                        
            ax2.errorbar(binCenters_psi, tab_events_data_psi/livetime[Datatype], np.sqrt(tab_events_data_psi)/livetime[Datatype], xerr=None, color='k', marker='o',markersize=4, lw=2, zorder=10, ls='none', label=Datatype+' ({} events)'.format(len(d_data_psi)))

            #ax.step(binCenters_psi,bestFit_energy['n1']*h_corsika,label='Original Corsika',where='mid',color='magenta',lw=1)
            ax2.step(binCenters_psi,bestFit['nominal']['n3']*tab_events_astro_psi['nominal']/livetime[Datatype],where='mid',color = 'b', label=r'Astrophysical $\nu$', alpha = 0.5, lw=2)
            ax2.step(binCenters_psi,bestFit['nominal']['n2']*tab_events_atm_psi['nominal']/livetime[Datatype],where='mid',color = 'g', label=r'Atmosperic $\nu$', alpha = 0.5, lw=2)
            ax2.step(binCenters_psi,bestFit['nominal']['n1']*tab_events_corsika_psi/livetime[Datatype],where='mid',color = 'orange', label='Corsika', alpha = 0.5, lw=2)

            ax2.fill_between(binCenters_psi, h_bkg_psi_syst_min, h_bkg_psi_syst_max, step='mid',color='r',alpha=0.5,lw=0)
            p3 = patches.Patch(color='r', alpha=0.5, linewidth=0)
            p4, = ax2.step(binCenters_psi,h_bkg_psi['nominal'],where='mid',color = 'r', lw=3)

            ax2.plot(0,0,lw=0,label='KS p-value: {:.3f}'.format(max(KS_pValues_psi)),zorder=100)
            ax2.plot(0,0,lw=0,label=r'Pearson $\chi^2$ p-value: {:.3f}'.format(p_syst_psi),zorder=200)
            
            handles, labels = ax2.get_legend_handles_labels()
            handles.append((p4,p3))
            labels.append(r'Combined background')
                                 
            order = [5,4,3,6,2,1,0]
                     
            leg = ax2.legend([handles[idx] for idx in order], [labels[idx] for idx in order],frameon = 1, fancybox=False, loc='upper left')

            #ax2.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (degrees)")
            #ax2.set_xlim(0,180)
            ax2.set_ylim(1e-7,1e-4)
            ax2.set_ylabel(r"rate (Hz)")
            ax2.set_yscale('log')
            plt.setp(ax2.get_xticklabels(), visible=False)

            
            ax2_ratio.errorbar(binCenters_psi, tab_events_data_psi/h_bkg_psi['nominal']/livetime[Datatype], np.sqrt(tab_events_data_psi)/h_bkg_psi['nominal']/livetime[Datatype], xerr=None, color='k', marker='o',markersize=4, lw=2, zorder=10, ls='none', label=Datatype+' ({} events)'.format(len(d_data_psi)))
            
            ax2_ratio.fill_between(binCenters_psi, h_bkg_psi_syst_min/h_bkg_psi['nominal'], h_bkg_psi_syst_max/h_bkg_psi['nominal'], step='mid',color='r',alpha=0.5,lw=0)

            ax2_ratio.step(binCenters_psi,h_bkg_psi['nominal']/h_bkg_psi['nominal'],where='mid',color = 'r', lw=3)
            
            
            ax2_ratio.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (degrees)")
            ax2_ratio.set_xlim(0,180)
            ax2_ratio.set_ylabel(r"ratio")
            ax2_ratio.set_ylim(0.8,1.2)
            #ax1.set_yscale('log')
            
            
            
            
            fig.suptitle(sample+r' Sample ($\Psi_{\mathrm{reco}}$>'+str(int(psiCut/np.pi*180))+'$^\circ$); '+
                         r'$n_{\mu}$='+'{:.2f}, '.format(bestFit['nominal']['n1'])+
                         r'$n_{\nu}^{atm}$='+'{:.2f}, '.format(bestFit['nominal']['n2'])+
                         r'$n_{\nu}^{ast}$='+'{:.2f}, '.format(bestFit['nominal']['n3'])
                         ,fontsize=24)

            fig.savefig('plots/PreUnblinding_ControlDistribution_'+Datatype+'_'+sample+'-Sample_'+options.PSICUT+'.png',bbox='tight')
            fig.savefig('plots/PreUnblinding_ControlDistribution_'+Datatype+'_'+sample+'-Sample_'+options.PSICUT+'.pdf',bbox='tight')


            #################


            fig2,(ax_cdf1,ax_cdf2) = plt.subplots(1,2,figsize=(20,8))

            xx_e = np.linspace(0,5,10000)
            xx_psi = np.linspace(0,np.pi,10000)
            
            
            cdf_energy = {}
            cdf_psi = {}
            for s in systematics:
                cdf_energy[s] = cdf_interpolation(np.log10(d_bkg_energy[s]), d_bkg_weight_energy[s])(xx_e)
                cdf_psi[s] = cdf_interpolation(d_bkg_psi[s], d_bkg_weight_psi[s])(xx_psi)
            
            cdf_energy_syst_min = np.array([])    
            cdf_energy_syst_max = np.array([]) 
            for i in range(len(xx_e)):
                cdf_energy_syst_min = np.append(cdf_energy_syst_min, np.min([cdf_energy[s][i] for s in systematics]))
                cdf_energy_syst_max = np.append(cdf_energy_syst_max, np.max([cdf_energy[s][i] for s in systematics]))

            cdf_psi_syst_min = np.array([])    
            cdf_psi_syst_max = np.array([]) 
            for i in range(len(xx_psi)):
                cdf_psi_syst_min = np.append(cdf_psi_syst_min, np.min([cdf_psi[s][i] for s in systematics]))
                cdf_psi_syst_max = np.append(cdf_psi_syst_max, np.max([cdf_psi[s][i] for s in systematics]))

                
                
            p1, = ax_cdf1.plot(10**xx_e,cdf_energy['nominal'],color='r',lw=2)        
            ax_cdf1.fill_between(10**xx_e, cdf_energy_syst_min, cdf_energy_syst_max, step='mid',color='r',alpha=0.5,lw=0)
            p2 = patches.Patch(color='r', alpha=0.5, linewidth=0)
    
            ax_cdf2.plot(xx_psi*180/np.pi,cdf_psi['nominal'],color='r',lw=2)
            ax_cdf2.fill_between(xx_psi*180/np.pi, cdf_psi_syst_min, cdf_psi_syst_max, step='mid',color='r',alpha=0.5,lw=0)

            cdf_data_energy = cdf_interpolation(np.log10(d_data_energy), d_data_weight.flatten() )(xx_e)
            cdf_data_psi = cdf_interpolation(d_data_psi, d_data_weight.flatten() )(xx_psi)

            ax_cdf1.plot(10**xx_e,cdf_data_energy,color='k',lw=2,label='Data')
            ax_cdf2.plot(xx_psi*180/np.pi,cdf_data_psi,color='k',lw=2,label='Data')

            ax_cdf1.set_xlabel(r"E$_{\mathrm{reco}}$ (GeV)")
            ax_cdf1.set_xscale('log')
            ax_cdf1.set_xlim(1e0,1e5)
            ax_cdf1.set_ylim(0,1)
            ax_cdf1.set_ylabel(r"CDF")

            ax_cdf2.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (degrees)")
            ax_cdf2.set_xlim(0,180)
            ax_cdf2.set_ylim(0,1)
            ax_cdf2.set_ylabel(r"CDF")

            handles, labels = ax_cdf2.get_legend_handles_labels()
            handles.append((p2,p1))
            labels.append(r'Combined background')
            
            leg1 = ax_cdf1.legend(handles, labels,frameon = 1, fancybox=False, loc='upper left')
            leg2 = ax_cdf2.legend(handles, labels,frameon = 1, fancybox=False, loc='upper left')
            
            fig2.suptitle(sample+r' Sample ($\Psi_{\mathrm{reco}}$>'+str(int(psiCut/np.pi*180))+'$^\circ$)' ,fontsize=24)

            fig2.savefig('plots/PreUnblinding_ControlDistribution_'+Datatype+'_'+sample+'-Sample_CDF_'+options.PSICUT+'.pdf',bbox='tight')
            fig2.savefig('plots/PreUnblinding_ControlDistribution_'+Datatype+'_'+sample+'-Sample_CDF_'+options.PSICUT+'.png',bbox='tight')