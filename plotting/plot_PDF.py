from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pickle
import sys, os

from optparse import OptionParser

base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'plotting')

import style
style.SetFigureStyle()
latex = style.latex
style.increaseAxisText(6)

plt.rcParams.update({'figure.dpi': 200.})

def merge_bins(histo, bins_merge, bins_merge2):
    grid_E   = histo[1]
    grid_E   = [grid_E[i] for i in range(len(grid_E)) if i%bins_merge ==0]

    grid_psi = histo[2]
    grid_psi = [grid_psi[i] for i in range(len(grid_psi)) if i%bins_merge2 ==0]

    lE=len(grid_E)
    lpsi=len(grid_psi)

    tab_events = np.zeros((lE-1,lpsi-1))
    for j in range(bins_merge):
        for k in range(lE-1):
            for j2 in range(bins_merge2):
                for k2 in range(lpsi-1):
                    tab_events[k,k2] = tab_events[k,k2] + histo[0][j+k*bins_merge,j2+k2*bins_merge2]
    return tab_events, grid_E, grid_psi

############################################

parser = OptionParser()
parser.add_option("-t", "--type",default="background",
                  dest="TYPE", help="Define type of PDF (background, annihilation, decay)")
parser.add_option("-c", "--channel",default="nu",
                  dest="CHANNEL", help="Annihilation channel")
parser.add_option("-m", "--mass", default='0',
                  dest="MASS", help="DM mass")
parser.add_option("-p", "--profile", default='NFW',
                  dest="PROFILE", help="DM halo profile")
parser.add_option("-s", "--sample", default='LE',
                  dest="SAMPLE", help="data sample (LE or HE)")
parser.add_option("-o", "--oversampling", default='-1',
                  dest="OVERSAMPLING", help="oversampling factor")

(options,args) = parser.parse_args()    


syst = 'nominal'

############################################
# case 1 : background

if options.TYPE == 'background':

	BDTstring = ''
	if options.SAMPLE == 'LE':
		BDTstring = 'LEBDT0.15_HEBDT0.2'
	elif options.SAMPLE == 'HE':
		BDTstring = 'LEBDT-1.0_HEBDT0.3'

	h1 = pickle.load(open(base_path+'PDFs/background/corsika/'+syst+'/PDF_nominal_background_corsika_'+BDTstring+'_2D.pkl','r'))

	fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(30,8))
	plt.rcParams.update({'figure.titlesize': 32})
	plt.suptitle('Atmospheric muon background')

	pcm = ax1.pcolormesh(merge_bins(h1,2,5)[2],merge_bins(h1,2,5)[1],merge_bins(h1,2,5)[0],norm=colors.LogNorm(vmin=1e-12, vmax=1e-6),cmap='Oranges')
	fig.colorbar(pcm, ax=ax1, label='rate (Hz)')
	ax1.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax1.set_ylabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax1.set_yscale('log')

	ax2.step(np.array(merge_bins(h1,2,90)[1][:-1]).flatten(),merge_bins(h1,2,90)[0].flatten(),where='post', lw=3, color = 'k')
	ax2.set_xlabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax2.set_ylabel(r"rate (Hz)")
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.set_ylim(1e-8,1e-4)

	ax3.step(np.array(merge_bins(h1,48,5)[2][:-1]).flatten(),merge_bins(h1,48,5)[0].flatten(),where='post', lw=3, color = 'k')
	ax3.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax3.set_ylabel(r"rate (Hz)")
	ax3.set_yscale('log')
	ax3.set_ylim(1e-8,1e-4)

	fig.subplots_adjust(wspace=0.3)

	fig.savefig('plots/PDF_Corsika.png')
	fig.savefig('plots/PDF_Corsika.pdf')


	############################################

	syst = 'nominal'

	h1 = pickle.load(open(base_path+'PDFs/background/atm/'+syst+'/PDF_'+syst+'_background_atm_'+BDTstring+'_2D.pkl','r'))

	fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(30,8))
	plt.rcParams.update({'figure.titlesize': 32})
	plt.suptitle('Atmospheric neutrino background (Honda 2006)')

	pcm = ax1.pcolormesh(merge_bins(h1,2,5)[2],merge_bins(h1,2,5)[1],merge_bins(h1,2,5)[0],norm=colors.LogNorm(vmin=1e-12, vmax=1e-6),cmap='Oranges')
	fig.colorbar(pcm, ax=ax1, label='rate (Hz)')
	ax1.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax1.set_ylabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax1.set_yscale('log')

	ax2.step(np.array(merge_bins(h1,2,90)[1][:-1]).flatten(),merge_bins(h1,2,90)[0].flatten(),where='post', lw=3, color = 'k')
	ax2.set_xlabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax2.set_ylabel(r"rate (Hz)")
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.set_ylim(1e-10,1e-4)

	ax3.step(np.array(merge_bins(h1,48,5)[2][:-1]).flatten(),merge_bins(h1,48,5)[0].flatten(),where='post', lw=3, color = 'k')
	ax3.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax3.set_ylabel(r"rate (Hz)")
	ax3.set_yscale('log')
	ax3.set_ylim(1e-10,1e-4)

	fig.subplots_adjust(wspace=0.3)


	fig.savefig('plots/PDF_Atm.png')
	fig.savefig('plots/PDF_Atm.pdf')

	############################################

	h1 = pickle.load(open(base_path+'PDFs/background/astro/'+syst+'/PDF_nominal_background_astro_'+BDTstring+'_2D.pkl','r'))

	fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(30,8))
	plt.rcParams.update({'figure.titlesize': 32})
	plt.suptitle('Astrophysical neutrino background (7.5 year HESE)')

	pcm = ax1.pcolormesh(merge_bins(h1,2,5)[2],merge_bins(h1,2,5)[1],merge_bins(h1,2,5)[0],norm=colors.LogNorm(vmin=1e-12, vmax=1e-6),cmap='Oranges')
	fig.colorbar(pcm, ax=ax1, label='rate (Hz)')
	ax1.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax1.set_ylabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax1.set_yscale('log')

	ax2.step(np.array(merge_bins(h1,2,90)[1][:-1]).flatten(),merge_bins(h1,2,90)[0].flatten(),where='post', lw=3, color = 'k')
	ax2.set_xlabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax2.set_ylabel(r"rate (Hz)")
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.set_ylim(1e-10,1e-4)

	ax3.step(np.array(merge_bins(h1,48,5)[2][:-1]).flatten(),merge_bins(h1,48,5)[0].flatten(),where='post', lw=3, color = 'k')
	ax3.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax3.set_ylabel(r"rate (Hz)")
	ax3.set_yscale('log')
	ax3.set_ylim(1e-10,1e-4)

	fig.subplots_adjust(wspace=0.3)

	fig.savefig('plots/PDF_Astro.png')
	fig.savefig('plots/PDF_Astro.pdf')
    
    
    ############################################

	h1 = pickle.load(open(base_path+'PDFs/corsika_KDE/corsika_KDE/nominal/PDF_nominal_background_corsika_KDE_'+BDTstring+'_2D.pkl','r'))

	fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(30,8))
	plt.rcParams.update({'figure.titlesize': 32})
	plt.suptitle('New Corsika: oversample + KDE')

	pcm = ax1.pcolormesh(merge_bins(h1,2,5)[2],merge_bins(h1,2,5)[1],merge_bins(h1,2,5)[0],norm=colors.LogNorm(vmin=1e-12, vmax=1e-6),cmap='Oranges')
	fig.colorbar(pcm, ax=ax1, label='rate (Hz)')
	ax1.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax1.set_ylabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax1.set_yscale('log')

	ax2.step(np.array(merge_bins(h1,2,90)[1][:-1]).flatten(),merge_bins(h1,2,90)[0].flatten(),where='post', lw=3, color = 'k')
	ax2.set_xlabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax2.set_ylabel(r"rate (Hz)")
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.set_ylim(1e-8,1e-4)

	ax3.step(np.array(merge_bins(h1,48,5)[2][:-1]).flatten(),merge_bins(h1,48,5)[0].flatten(),where='post', lw=3, color = 'k')
	ax3.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax3.set_ylabel(r"rate (Hz)")
	ax3.set_yscale('log')
	ax3.set_ylim(1e-8,1e-4)

	fig.subplots_adjust(wspace=0.3)

	fig.savefig('plots/PDF_Corsika_KDE.png')
	fig.savefig('plots/PDF_Corsika_KDE.pdf')
    
    
    
    ############################################

	h1 = pickle.load(open(base_path+'PDFs/corsika_KDE/corsika/nominal/PDF_nominal_background_corsika_'+BDTstring+'_2D.pkl','r'))

	fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(30,8))
	plt.rcParams.update({'figure.titlesize': 32})
	plt.suptitle('New Corsika: oversample')

	pcm = ax1.pcolormesh(merge_bins(h1,2,5)[2],merge_bins(h1,2,5)[1],merge_bins(h1,2,5)[0],norm=colors.LogNorm(vmin=1e-12, vmax=1e-6),cmap='Oranges')
	fig.colorbar(pcm, ax=ax1, label='rate (Hz)')
	ax1.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax1.set_ylabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax1.set_yscale('log')

	ax2.step(np.array(merge_bins(h1,2,90)[1][:-1]).flatten(),merge_bins(h1,2,90)[0].flatten(),where='post', lw=3, color = 'k')
	ax2.set_xlabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax2.set_ylabel(r"rate (Hz)")
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.set_ylim(1e-8,1e-4)

	ax3.step(np.array(merge_bins(h1,48,5)[2][:-1]).flatten(),merge_bins(h1,48,5)[0].flatten(),where='post', lw=3, color = 'k')
	ax3.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax3.set_ylabel(r"rate (Hz)")
	ax3.set_yscale('log')
	ax3.set_ylim(1e-8,1e-4)

	fig.subplots_adjust(wspace=0.3)

	fig.savefig('plots/PDF_Corsika_Oversampling.png')
	fig.savefig('plots/PDF_Corsika_Oversampling.pdf')    
    

############################################
#case 2: signal

if options.TYPE in ['annihilation', 'decay']:
	mass = int(options.MASS)
	channel = options.CHANNEL
	profile = options.PROFILE
	oversampling = options.OVERSAMPLING

	BDTstring = ''
	if options.SAMPLE == 'LE':
		BDTstring = 'LEBDT0.15_HEBDT0.2'
	elif options.SAMPLE == 'HE':
		BDTstring = 'LEBDT-1.0_HEBDT0.3'

	h1 = pickle.load(open(base_path+'PDFs/'+options.TYPE+'/'+channel+'/nominal/PDF_nominal_DM'+options.TYPE+'_'+profile+'profile_'+BDTstring+'_2D_'+channel+'_'+str(mass)+'GeV_oversampling'+oversampling+'.pkl','r'))

	fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(30,8))
	plt.rcParams.update({'figure.titlesize': 32})
	plt.suptitle(r'Signal m$_\mathrm{\chi}$='+str(mass/1000)+' TeV, annihilation, '+profile+' profile')

	pcm = ax1.pcolormesh(merge_bins(h1,2,5)[2],merge_bins(h1,2,5)[1],merge_bins(h1,2,5)[0],norm=colors.LogNorm(vmin=1e-12, vmax=1e-6),cmap='Oranges')
	fig.colorbar(pcm, ax=ax1, label='rate (Hz)')
	ax1.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax1.set_ylabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax1.set_yscale('log')

	ax2.step(np.array(merge_bins(h1,2,90)[1][:-1]).flatten(),merge_bins(h1,2,90)[0].flatten(),where='post', lw=3, color = 'k')
	ax2.set_xlabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax2.set_ylabel(r"rate (Hz)")
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.set_ylim(1e-10,1e-4)

	ax3.step(np.array(merge_bins(h1,48,5)[2][:-1]).flatten(),merge_bins(h1,48,5)[0].flatten(),where='post', lw=3, color = 'k')
	ax3.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax3.set_ylabel(r"rate (Hz)")
	ax3.set_yscale('log')
	ax3.set_ylim(1e-10,1e-4)

	fig.subplots_adjust(wspace=0.3)

	fig.savefig('plots/PDF_Signal_'+str(mass/1000)+'TeV_'+options.TYPE+'_'+profile+'.png')
	fig.savefig('plots/PDF_Signal_'+str(mass/1000)+'TeV_'+options.TYPE+'_'+profile+'.pdf')


	############################################
	## public plot

	fig,ax1 = plt.subplots(1,1,figsize=(10,8))
	plt.rcParams.update({'figure.titlesize': 22})
	plt.suptitle(r'm$_\mathrm{\chi}$='+str(mass/1000)+r' TeV, $\chi\chi\rightarrow\nu\bar{\nu}$, '+profile)

	pcm = ax1.pcolormesh(h1[2],h1[1],h1[0]/np.sum(h1[0]),norm=colors.LogNorm(vmin=1e-6, vmax=1e-2),cmap='Oranges')
	fig.colorbar(pcm, ax=ax1, label='probability')
	ax1.set_xlabel(r"$\Psi_{\mathrm{reco}}$ (rad)")
	ax1.set_ylabel(r"E$_{\mathrm{reco}}$ (GeV)")
	ax1.set_yscale('log')

	plt.text(0.1, 5e4, 'IceCube work in progress', {'color': 'maroon', 'fontsize': 20})

	#fig.tight_layout()
	fig.savefig('plots/PDF_public_'+str(mass/1000)+'TeV_'+options.TYPE+'_'+profile+'.png',bbox_inches = 'tight') 
        fig.savefig('plots/PDF_public_'+str(mass/1000)+'TeV_'+options.TYPE+'_'+profile+'.pdf',bbox_inches = 'tight')
