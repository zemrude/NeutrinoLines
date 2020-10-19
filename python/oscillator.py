try:
    import nuflux
except Exception as e :
    raise Exception("The required `nuflux` project is missing (https://github.com/IceCubeOpenSource/nuflux/tree/master/docs)")

   
newflux = nuflux.makeFlux('honda2006').getFlux

try:
    from icecube import prob3
except Exception as e: 
    raise Exception("The required `prob3` project is missing (https://code.icecube.wisc.edu/projects/icecube/browser/IceCube/sandbox/olivas/prob3)")
#from icecube.nuCraft_icetray.nuCraft import NuCraft

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

def neutrinoOscillator_newParameters(atmH=None, depth=None):
    # parameters from arxiv1611.01514
    theta23 = np.arcsin(np.sqrt(0.582))/np.pi*180.
    theta13 = np.arcsin(np.sqrt(0.02240))/np.pi*180.
    theta12 = np.arcsin(np.sqrt(0.310))/np.pi*180.;
    DM21   = 7.39e-5
    DM31   = 2.525e-3
    
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

    numPrec = 1e-5
   
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

    
    
print ('load Honda flux')

Honda = Honda_Flux(base_path+'/resources/Honda/')
Honda.load()

def getHondaWeight(E,cosZenith,pdg):
    cosZenithBins = np.linspace(1, -1, 21)
    iEnergyBin = np.searchsorted(Honda.getFlux('nuE',0)[0], E)
    iCosZenithBin = np.searchsorted(-cosZenithBins, -cosZenith)

    flux_weight = 0.
    if pdg>0.:
        if np.abs(pdg) == 12:
            flux_weight = Honda.getFlux('nuE',iCosZenithBin-1)[1][iEnergyBin-1]
        elif np.abs(pdg) == 14:
            flux_weight = Honda.getFlux('nuMu',iCosZenithBin-1)[1][iEnergyBin-1]
        else:
            flux_weight = 0.
    if pdg<0.:
        if np.abs(pdg) == 12:
            flux_weight = Honda.getFlux('nuEbar',iCosZenithBin-1)[1][iEnergyBin-1]
        elif np.abs(pdg) == 14:
            flux_weight = Honda.getFlux('nuMubar',iCosZenithBin-1)[1][iEnergyBin-1]
        else:
            flux_weight = 0.

    return flux_weight

getHondaWeight = np.vectorize(getHondaWeight)


def osc_Atm_flux_weight_honda2015(energy,zenith,nuType,nC,atm_mode=0):
    
    numPrec = 1e-5
   
    # assume nu_e at production
    prob_e   = nC.CalcWeights((nuType/np.abs(nuType)*12, energy, zenith), numPrec=numPrec,atmMode=atm_mode)
    flux_e   = getHondaWeight(energy, np.cos(zenith),nuType)/1e4

    # assume nu_mu at production
    prob_mu  = nC.CalcWeights((nuType/np.abs(nuType)*14, energy, zenith), numPrec=numPrec,atmMode=atm_mode)
    flux_mu  = getHondaWeight(energy, np.cos(zenith), nuType)/1e4

    # assume nu_tau at production
    prob_tau = nC.CalcWeights((nuType/np.abs(nuType)*16, energy, zenith), numPrec=numPrec,atmMode=atm_mode)
    flux_tau   = getHondaWeight(energy, np.cos(zenith),nuType)/1e4
    
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
   
