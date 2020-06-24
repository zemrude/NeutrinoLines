from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches

import pickle
import sys, os

base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'plotting')

import style
style.SetFigureStyle()
latex = style.latex

plt.rcParams.update({'figure.dpi': 300.})

style.increaseAxisText(10)
style.increaseLegendText(1)


def extractMedianLimit(LECut,HECut,mode, channel, syst, profile, m_min, m_max, binning1, binning2, oversampling=100, CL='_90CL'):
    all_masses = [40, 63, 100, 158, 251, 398, 631, 1000, 1585, 2512, 3981, 6310, 10000, 15850, 25120, 39810]

    x = []
    y = []
    for m in all_masses:
        if m < m_min or m>m_max:
            continue
            
        
        sensFile_name = 'Sensitivity_effective_'+syst+'_'+mode+'_'+profile+'profile_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_'+channel+'_'+str(int(m))+'GeV_binsmerged_'+str(binning1)+'-'+str(binning2)+'_oversampling'+str(oversampling)+CL+'.npy'
        sensFile = base_path+'sensitivity/'+mode+'/'+channel+'/'+syst+'/'+sensFile_name
        try: 
            results = np.load(sensFile).item()
            x.append(results['mass'])
            y.append(results['median'])
        except:
            if (LECut>0.) and (m>1000):
                continue # LE Cut and HE mass
            if (HECut>0.2) and (m<1000):
                continue # HE Cut and LE mass
            if channel=='W' and m<180:
                continue # smaller W mass
            else:
                print 'file',sensFile_name,'not found'
                
    return x,y    


def extractBrazilLimit(LECut,HECut,mode, channel,syst, profile, m_min, m_max, binning1, binning2, oversampling=100, CL='_90CL'):
    all_masses = [10, 16, 25, 40, 63, 100, 158, 251, 398, 631, 1000, 1585, 2512, 3981, 6310, 10000, 15850, 25120, 39810]

    x = []
    y = {}
    for key in ['median','error_95_low','error_95_high','error_68_low','error_68_high']: 
        y[key] = []

    for m in all_masses:
        if m < m_min or m>m_max:
            continue
          
        sensFile_name = 'Sensitivity_effective_'+syst+'_'+mode+'_'+profile+'profile_LEBDT'+str(LECut)+'_HEBDT'+str(HECut)+'_'+channel+'_'+str(int(m))+'GeV_binsmerged_'+str(binning1)+'-'+str(binning2)+'_oversampling'+str(oversampling)+CL+'.npy'
        sensFile = base_path+'sensitivity/'+mode+'/'+channel+'/'+syst+'/'+sensFile_name
        try: 
            results = np.load(sensFile).item()
            x.append(results['mass'])
            
            for key in ['median','error_95_low','error_95_high','error_68_low','error_68_high']:
                y[key].append(results[key])

        except:
            if (LECut>0.) and (m>1000):
                continue # LE Cut and HE mass
            if (HECut>0.2) and (m<1000):
                continue # HE Cut and LE mass
            if channel=='W' and m<180:
                continue # smaller W mass
            else:
                print 'file',sensFile_name,'not found'

    return x,y    

def sens_fig_ann():
    fig,ax = plt.subplots(1,1,figsize=(10,8))

    ax.set_ylabel(r"$\langle\sigma v \rangle$ (cm$^{3}$s$^{-1}$)")
    ax.set_xlabel(r"m$_\chi$ (GeV)")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(5,100000)
    ax.set_ylim(1e-25,1e-22)
    ax.tick_params(axis='x', pad=10)
    ax.tick_params(axis='y', pad=10)
    ax.yaxis.labelpad = 20

    return fig, ax


def sens_fig_dec():
    fig,ax = plt.subplots(1,1,figsize=(10,8))

    ax.set_ylabel(r"$\tau_\chi$ (s)")
    ax.set_xlabel(r"m$_\chi$ (GeV)")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(5,100000)
    ax.tick_params(axis='x', pad=10)
    ax.tick_params(axis='y', pad=10)
    ax.yaxis.labelpad = 20

    return fig, ax


def plot_results(mode, channel, profile):

    fig,ax = plt.subplots(1,1,figsize=(10,8))
    
    color = {}
    color['ANTARES']='b'
    color['IC86']='r'
    color['IC22']='g'
    color['IC86_dec']='r'
    color['ICDC']='g'
    color['HESS']='orange'
    color['Fermi']='cyan'
    color['Fermi_dec']='cyan'
    color['HAWK']='orange'
    color['SuperK']='magenta'
    
    ls={}
    ls['ANTARES']='-.'
    ls['IC86']='--'
    ls['IC22']='--'
    ls['IC86_dec']='--'
    ls['ICDC']='--'
    ls['HESS']='-.'
    ls['Fermi']='-.'
    ls['Fermi_dec']='-.'
    ls['HAWK']='-.'
    ls['SuperK']=':'

    label={}
    label['ANTARES']='ANTARES 90% CL [PLB 769 (2017), PLB 796 (2019)]'
    label['IC86']='IceCube Tracks 90% CL [EPJC 77 (2017)]'
    label['IC22']='IceCube 90% CL (Einasto) [PRD 84 (2011)]'
    label['IC86_dec']='IceCube 90% CL (Burkert) [EPJC 78 (2018)]'
    label['ICDC']='IceCube Cascades 90% CL [EPJC 76 (2016)]'
    label['HESS']='HESS 95% CL (Einasto) [PRL 117 (2016)]'
    label['Fermi'] = 'Fermi-LAT + MAGIC 95% CL (dSph) [JCAP 1602 (2016)]'
    label['Fermi_dec'] = 'Fermi-LAT 95% CL [APJ 761 (2012)]'
    label['HAWK'] = 'HAWK 95% CL (dSph) [APJ 853 (2018)]'
    label['SuperK'] = 'SuperK 90% CL [arXiv:2005.05106]'

    if mode == 'annihilation':
        
        for ex in ['HESS']:
            filename = base_path+'resources/paperLimits/annihilation/'+ex+'/'+ex+'_'+channel+'_'+'Einasto'+'.csv'
            try:
                inFile = open(filename, 'r')
                data = np.loadtxt(inFile,delimiter=',')
                x = data[:,0]
                y = data[:,1]
                x,y = zip(*sorted(zip(x, y)))
                ax.plot(x,y, lw=3,ls=ls[ex],color=color[ex],label=label[ex])
            except:
                continue
                
        for ex in ['SuperK','Fermi','ANTARES','ICDC','IC86']:
            filename = base_path+'resources/paperLimits/annihilation/'+ex+'/'+ex+'_'+channel+'_'+profile+'.csv'
            try:
                inFile = open(filename, 'r')
                data = np.loadtxt(inFile,delimiter=',')
                x = data[:,0]
                y = data[:,1]
                x,y = zip(*sorted(zip(x, y)))
                ax.plot(x,y, lw=3,ls=ls[ex],color=color[ex],label=label[ex])
            except:
                continue
                
    if mode == 'decay':

        
        for ex in ['Fermi_dec']:
            try:
                filename = base_path+'resources/paperLimits/decay/Fermi_'+channel+'_NFW.csv'
                inFile = open(filename, 'r')
                data = np.loadtxt(inFile,delimiter=',')
                x = data[:,0]
                y = data[:,1]
                x,y = zip(*sorted(zip(x, y)))
                ax.plot(x,y, lw=2,ls=ls[ex],color=color[ex],label=label[ex])
            except:
                continue
                
        for ex in ['HAWK']:
            try:
                filename = base_path+'resources/paperLimits/decay/HAWK_'+channel+'_dSph.csv'
                inFile = open(filename, 'r')
                data = np.loadtxt(inFile,delimiter=',')
                x = data[:,0]
                y = data[:,1]
                x,y = zip(*sorted(zip(x, y)))
                ax.plot(x,y, lw=2,ls=ls[ex],color=color[ex],label=label[ex])
            except:
                continue
                    
        for ex in ['IC22']:
            try:
                filename = base_path+'resources/paperLimits/decay/'+'IC22'+'_'+channel+'_'+'Einasto'+'.csv'
                inFile = open(filename, 'r')
                data = np.loadtxt(inFile,delimiter=',')
                x = data[:,0]
                y = data[:,1]
                x,y = zip(*sorted(zip(x, y)))
                ax.plot(x,y, lw=3,ls=ls[ex],color=color[ex],label=label[ex])
            except:
                continue
                
                
        for ex in ['IC86_dec']:
            try:
                filename = base_path+'resources/paperLimits/decay/'+'IC86'+'_'+channel+'_'+'Burkert'+'.csv'
                inFile = open(filename, 'r')
                data = np.loadtxt(inFile,delimiter=',')
                x = data[:,0]
                y = data[:,1]
                x,y = zip(*sorted(zip(x, y)))
                ax.plot(x,y, lw=3,ls=ls[ex],color=color[ex],label=label[ex])
            except:
                continue
      
    
    
    # This work
    
    oversampling_LE = -1
    oversampling_HE = -1
    oversampling_HE_Burkert = -1

        
    if channel == 'nue':
        oversampling_LE = 100
        oversampling_HE = 200
        oversampling_HE_Burkert = 100

        if mode == 'decay':
            oversampling_HE = 100
            
        if profile == 'Burkert':
            oversampling_HE = 100
        
    transition = 900
    sens_combined_1 = extractBrazilLimit(0.15,0.2,mode,channel,'nominal',profile,30,transition, 2, 5, oversampling_LE, '_90CL')
    sens_combined_2 = extractBrazilLimit(-1.0,0.3,mode,channel,'nominal',profile,transition,50000, 2, 5, oversampling_HE, '_90CL')    

    sens_combined = []
    sens_combined.append(np.append(np.array(sens_combined_1[0]),np.array(sens_combined_2[0])))
    sens_combined.append({})

    for key in sens_combined_1[1].keys():
        sens_combined[1][key] = np.append(sens_combined_1[1][key],sens_combined_2[1][key])

    p1, = ax.plot(sens_combined[0], sens_combined[1]['median'], 'k-',lw=2,markersize=8)
    ax.fill_between(sens_combined[0], sens_combined[1]['error_95_low'], sens_combined[1]['error_95_high'], color='yellow',alpha=0.8) 
    ax.fill_between(sens_combined[0], sens_combined[1]['error_68_low'], sens_combined[1]['error_68_high'], color='chartreuse',alpha=0.8)
    p2 = patches.Patch(color='yellow', alpha=1, linewidth=0)
    p3 = patches.Patch(color='chartreuse', alpha=1, linewidth=0)
    
    if mode in ['decay','annihilation']:
        sens_combined_2_1 = extractBrazilLimit(0.15,0.2,mode,channel,'nominal','Burkert',30,transition, 2, 5, oversampling_LE, '_90CL')
        sens_combined_2_2 = extractBrazilLimit(-1.0,0.3,mode,channel,'nominal','Burkert',transition,50000, 2, 5, oversampling_HE_Burkert, '_90CL')    

        sens_combined_2 = []
        sens_combined_2.append(np.append(np.array(sens_combined_2_1[0]),np.array(sens_combined_2_2[0])))
        sens_combined_2.append({})

        for key in sens_combined_2_1[1].keys():
            sens_combined_2[1][key] = np.append(sens_combined_2_1[1][key],sens_combined_2_2[1][key])

        p4, = ax.plot(sens_combined_2[0], sens_combined_2[1]['median'], 'k--',lw=2,markersize=8, label='5 year median sensitivity 90% CL (this work, Burkert)')

    handles, labels = ax.get_legend_handles_labels()
    handles.append((p2,p1))
    labels.append(r'5 year median sensitivity 90% CL $\pm 2\,\sigma$ (this work)')
    handles.append((p3,p1))
    labels.append(r'5 year median sensitivity 90% CL $\pm 1\,\sigma$ (this work)')

    if mode == 'annihilation':
        ax.plot([0.01,1e5],[3e-26,3e-26], lw=6,ls='-',color='lightgray',label='')
        plt.text(110, 1.3e-26,  r'relic density', {'color': 'gray', 'fontsize': 18})

        ax.set_ylabel(r"$\langle\sigma \upsilon \rangle$ (cm$^{3}$s$^{-1}$)")
        ax.set_ylim(5e-27,1e-21)
        #plt.text(3000, 1e-21/5, 'IceCube\n'+'work in progress', {'color': 'maroon', 'fontsize': 20})
        plt.text(12, 1e-21/3,  r'$\chi\chi\rightarrow$'+latex[channel]+'\n'+profile+' profile',fontsize=17, color = 'black')
  
    if mode == 'decay':
        ax.set_ylabel(r"$\tau_\chi$ (s)")
        ax.set_ylim(1e24,1e29)
        #plt.text(3000, 1e29/5, 'IceCube\n'+'work in progress', {'color': 'maroon', 'fontsize': 20})
        plt.text(12, 1e29/3,  r'$\chi\rightarrow$'+latex[channel]+'\n'+profile+' profile',fontsize=17, color = 'black')
        
    ax.set_xlabel(r"m$_\chi$ (GeV)")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(8,50000)
    ax.tick_params(axis='x', pad=10)
    ax.tick_params(axis='y', pad=10)
    ax.yaxis.labelpad = 20

    leg = ax.legend(handles[::-1], labels[::-1] ,frameon = 1, fancybox=False, loc='lower left', bbox_to_anchor=(-0.0175, 1., 1.035, 0),labelspacing=0.2, mode="expand")

    fig.tight_layout()
    fig.savefig('plots/Sensitivity_'+mode+'_'+channel+'_'+profile+'.pdf',bbox_inches = 'tight') 
    fig.savefig('plots/Sensitivity_'+mode+'_'+channel+'_'+profile+'.png',bbox_inches = 'tight') 
    
    return fig


###############################################
###############################################


fig,ax = sens_fig_ann()

sens_combined_1 = extractBrazilLimit(0.15,0.2,'annihilation','nue','nominal','NFW',30,900, 2, 5, 100)
sens_combined_2 = extractBrazilLimit(-1.0,0.3,'annihilation','nue','nominal','NFW',900,50000, 2, 5, 200)    

sens_combined = []
sens_combined.append(np.append(np.array(sens_combined_1[0]),np.array(sens_combined_2[0])))
sens_combined.append({})
                             
for key in sens_combined_1[1].keys():
    sens_combined[1][key] = np.append(sens_combined_1[1][key],sens_combined_2[1][key])

sens_LE = extractMedianLimit(0.15,0.2,'annihilation','nue','nominal','NFW',1,100000, 2, 5, 100)
sens_HE = extractMedianLimit(-1.0, 0.3,'annihilation','nue','nominal','NFW',300,100000, 2, 5, 200)

ax.plot(sens_combined[0], sens_combined[1]['median'], 'ko-',lw=3,markersize=8, label='Expected median sensitivity 90% CL')
ax.fill_between(sens_combined[0], sens_combined[1]['error_95_low'], sens_combined[1]['error_95_high'], color='yellow',alpha=1,label='95% containment') 
ax.fill_between(sens_combined[0], sens_combined[1]['error_68_low'], sens_combined[1]['error_68_high'], color='chartreuse',alpha=1,label='68% containment')

ax.plot(sens_LE[0],sens_LE[1],'b--',lw=2,label='LE Sample only')
ax.plot(sens_HE[0],sens_HE[1],'r--',lw=2,label='HE Sample only')

ax.set_ylim(1e-25,1e-23)

leg = ax.legend(frameon = 1, fancybox=False, loc='lower left', bbox_to_anchor=(-0.0175, 1., 1.035, 0),labelspacing=0.2, mode="expand", ncol=1)
fig.tight_layout()
fig.savefig('plots/Sensitivity_Samples_annihilation_nue_NFW.pdf',bbox_inches = 'tight') 
fig.savefig('plots/Sensitivity_Samples_annihilation_nue_NFW.png',bbox_inches = 'tight') 

###############################################
# -------------------- decay

fig,ax = sens_fig_dec()

sens_combined_1 = extractBrazilLimit(0.15,0.2,'decay','nue','nominal','NFW',30,900, 2, 5, 100)
sens_combined_2 = extractBrazilLimit(-1.0,0.3,'decay','nue','nominal','NFW',900,50000, 2, 5, 100)    

sens_combined = []
sens_combined.append(np.append(np.array(sens_combined_1[0]),np.array(sens_combined_2[0])))
sens_combined.append({})
                             
for key in sens_combined_1[1].keys():
    sens_combined[1][key] = np.append(sens_combined_1[1][key],sens_combined_2[1][key])

sens_LE = extractMedianLimit(0.15,0.2,'decay','nue','nominal','NFW',1,100000, 2, 5, 100)
sens_HE = extractMedianLimit(-1.0, 0.3,'decay','nue','nominal','NFW',300,100000, 2, 5, 100)

ax.plot(sens_combined[0], sens_combined[1]['median'], 'ko-',lw=3,markersize=8, label='Expected median sensitivity 90% CL')
ax.fill_between(sens_combined[0], sens_combined[1]['error_95_low'], sens_combined[1]['error_95_high'], color='yellow',alpha=1,label='95% containment') 
ax.fill_between(sens_combined[0], sens_combined[1]['error_68_low'], sens_combined[1]['error_68_high'], color='chartreuse',alpha=1,label='68% containment')

ax.plot(sens_LE[0],sens_LE[1],'b--',lw=2,label='LE Sample only')
ax.plot(sens_HE[0],sens_HE[1],'r--',lw=2,label='HE Sample only')

leg = ax.legend(frameon = 1, fancybox=False, loc='lower left', bbox_to_anchor=(-0.0175, 1., 1.035, 0),labelspacing=0.2, mode="expand", ncol=1)
fig.tight_layout()
fig.savefig('plots/Sensitivity_Samples_decay_nue_NFW.pdf',bbox_inches = 'tight') 
fig.savefig('plots/Sensitivity_Samples_decay_nue_NFW.png',bbox_inches = 'tight') 


###############################################

for profile in ['NFW', 'Burkert']:
    
    oversampling_nue_LE = 100
    oversampling_nue_HE = 100
    if profile == 'NFW':
        oversampling_nue_HE = 200
    

    fig,ax = sens_fig_ann()

    channels = ['tau','mu','b','W']

    sens_1 = {}
    sens_2 = {}
    sens_combined = {}

    c='nue'
    sens_1[c] = extractMedianLimit(0.15,0.2,'annihilation',c,'nominal',profile,30,900, 2, 5, oversampling_nue_LE, '_90CL')
    sens_2[c] = extractMedianLimit(-1.0,0.3,'annihilation',c,'nominal',profile,900,50000, 2, 5, oversampling_nue_HE, '_90CL')    
    sens_combined[c] = []
    sens_combined[c].append(np.append(np.array(sens_1[c][0]),np.array(sens_2[c][0])))
    sens_combined[c].append(np.append(np.array(sens_1[c][1]),np.array(sens_2[c][1])))

    ax.plot(sens_combined[c][0],sens_combined[c][1],lw=3,label=latex[c])


    for c in channels:
        sens_1[c] = extractMedianLimit(0.15,0.2,'annihilation',c,'nominal',profile,30,900, 2, 5, -1, '_90CL')
        sens_2[c] = extractMedianLimit(-1.0,0.3,'annihilation',c,'nominal',profile,900,50000, 2, 5, -1, '_90CL')    
        sens_combined[c] = []
        sens_combined[c].append(np.append(np.array(sens_1[c][0]),np.array(sens_2[c][0])))
        sens_combined[c].append(np.append(np.array(sens_1[c][1]),np.array(sens_2[c][1])))

        ax.plot(sens_combined[c][0],sens_combined[c][1],lw=3,label=latex[c])


    ax.set_ylim(1e-25,1e-21)
    plt.text(12*1e3, 1e-21/2, profile+' profile',fontsize=17, color = 'black')

    leg = ax.legend(frameon = 1, fancybox=False, loc='lower left', bbox_to_anchor=(-0.0175, 1., 1.035, 0),labelspacing=0.2, mode="expand", ncol=1)
    fig.tight_layout()
    fig.savefig('plots/Sensitivity_annihilation_channels_'+profile+'.pdf',bbox_inches = 'tight') 
    fig.savefig('plots/Sensitivity_annihilation_channels_'+profile+'.png',bbox_inches = 'tight') 


    ###############################################
    #--------------------- decay

    fig,ax = sens_fig_dec()

    channels = ['tau','mu','b','W']

    sens_1 = {}
    sens_2 = {}
    sens_combined = {}

    c='nue'
    sens_1[c] = extractMedianLimit(0.15,0.2,'decay',c,'nominal',profile,30,900, 2, 5, oversampling_nue_LE, '_90CL')
    sens_2[c] = extractMedianLimit(-1.0,0.3,'decay',c,'nominal',profile,900,50000, 2, 5, 100, '_90CL')    
    sens_combined[c] = []
    sens_combined[c].append(np.append(np.array(sens_1[c][0]),np.array(sens_2[c][0])))
    sens_combined[c].append(np.append(np.array(sens_1[c][1]),np.array(sens_2[c][1])))

    ax.plot(sens_combined[c][0],sens_combined[c][1],lw=3,label=latex[c])


    for c in channels:
        sens_1[c] = extractMedianLimit(0.15,0.2,'decay',c,'nominal',profile,30,900, 2, 5, -1, '_90CL')
        sens_2[c] = extractMedianLimit(-1.0,0.3,'decay',c,'nominal',profile,900,50000, 2, 5, -1, '_90CL')    
        sens_combined[c] = []
        sens_combined[c].append(np.append(np.array(sens_1[c][0]),np.array(sens_2[c][0])))
        sens_combined[c].append(np.append(np.array(sens_1[c][1]),np.array(sens_2[c][1])))

        ax.plot(sens_combined[c][0],sens_combined[c][1],lw=3,label=latex[c])


    ax.set_ylim(1e22,1e29)
    plt.text(12*1e3, 1e29/2, profile+' profile',fontsize=17, color = 'black')

    leg = ax.legend(frameon = 1, fancybox=False, loc='lower left', bbox_to_anchor=(-0.0175, 1., 1.035, 0),labelspacing=0.2, mode="expand", ncol=1)
    fig.tight_layout()
    fig.savefig('plots/Sensitivity_decay_channels_'+profile+'.pdf',bbox_inches = 'tight') 
    fig.savefig('plots/Sensitivity_decay_channels_'+profile+'.png',bbox_inches = 'tight') 


###############################################


fig = plot_results('annihilation', 'nue', 'NFW')
fig = plot_results('annihilation', 'tau', 'NFW')
fig = plot_results('decay', 'nue', 'NFW')
fig = plot_results('decay', 'tau', 'NFW')
