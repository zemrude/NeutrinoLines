import numpy as np
import scipy.special as sps
from iminuit import Minuit

class Profile_Analyser:

    def __init__(self):
        
        self.LLHtype = None

        self.ready = False
        self.signalPDF = None
        self.backgroundPDF = None


        self.signalPDF_uncert2 = None
        self.backgroundPDF_uncert2 = None

        self.signalContaminationPDF = None
        self.subtractSignalContamination = False
        
        self.nBackgroundEvents = 0.
        self.nSignalEvents = 0.

        self.nbins = 0
        
        self.livetime = -1.
        self.observation = None

        self.computedBestFit = False
        self.bestFit = None
        self.TS = None
        
        self.samplingMethod = 'default'
        self.moreOutput = False
        
    def setLivetime(self,lt):
        self.livetime = lt

    def setLLHtype(self,type):
        availableTypes = ['Poisson', 'Effective']
        if type not in availableTypes:
            raise ValueError('LLH type not implemented yet. Choose amongst: '+str(availableTypes))
        else:
            self.LLHtype = type
            
    def saveMoreOutput(self):
        self.moreOutput = True


    def loadBackgroundPDF(self,pdf, verbose = False):
        if self.livetime < 0:
            raise ValueError('Livetime of the analysis is not defined yet. Please do this first!')
            
        self.backgroundPDF = pdf.flatten()*self.livetime
        self.nBackgroundEvents = np.sum(self.backgroundPDF)
        if verbose:
            print('total background events:', self.nBackgroundEvents)
        
        self.nbins = len(pdf)

    def loadSignalPDF(self,pdf):
        if self.livetime < 0:
            raise ValueError('Livetime of the analysis is not defined yet. Please do this first!')
        self.signalPDF = pdf.flatten()*self.livetime
        self.nSignalEvents = np.sum(self.signalPDF)
        print('total signal events:', np.sum(self.signalPDF))
        if self.nbins == len(pdf):
            self.ready = True
        else:
            raise ValueError('Shape of signal pdf does not match the background pdf! Did you initialize the background pdf first?')
            
    def loadUncertaintyPDFs(self,bkg_pdf,sig_pdf):
        self.backgroundPDF_uncert2 = bkg_pdf.flatten()*self.livetime*self.livetime        
        self.signalPDF_uncert2 = sig_pdf.flatten()*self.livetime*self.livetime

        if self.nbins != len(bkg_pdf):
            raise ValueError('Shape of background uncertainty pdf does not match the background pdf!')
        if self.nbins != len(sig_pdf):
            raise ValueError('Shape of background uncertainty pdf does not match the background pdf!')
            
            
    def loadSignalContaminationPDF(self,sig_cont_pdf, sig_cont_unc_pdf):
        
        self.signalContaminationPDF = sig_cont_pdf.flatten()*self.livetime
        self.signalContaminationPDF_uncert2 = sig_cont_unc_pdf.flatten()*self.livetime*self.livetime
        
        if self.nbins != len(sig_cont_pdf):
            raise ValueError('Shape of background uncertainty pdf does not match the background pdf!')
            
        self.subtractSignalContamination = True
        
 
    def sampleObservation(self,n1,nsig):

        if not self.ready:
            raise ValueError('Not all pdfs are correctly loaded!')

        observationPDF = n1*self.backgroundPDF + nsig*self.signalPDF

        self.observation=np.zeros(np.shape(self.backgroundPDF))
        for i in range(len(self.observation)):
            self.observation[i]=np.random.poisson(observationPDF[i])
        
        self.computedBestFit = False

        
    def evaluateLLH(self, n1, nsig):
        
        if self.subtractSignalContamination:
            modelPDF = n1*self.backgroundPDF + nsig*self.signalPDF - nsig*self.signalContaminationPDF
        else:
            modelPDF = n1*self.backgroundPDF + nsig*self.signalPDF
            
        if np.isnan(modelPDF).any():
            print('nan in model array with n1,ns=',n1,nsig, self.computedBestFit)

        if self.LLHtype == 'Poisson':
            bins_to_use = (modelPDF>0.)
            values = self.observation[bins_to_use]*np.log(modelPDF[bins_to_use])-modelPDF[bins_to_use]

        elif self.LLHtype == 'Effective':
        
            if self.subtractSignalContamination:
                modelPDF_uncert2 = n1*n1*self.backgroundPDF_uncert2 + nsig*nsig*self.signalPDF_uncert2 - nsig*nsig*self.signalContaminationPDF_uncert2
            else:
                modelPDF_uncert2 = n1*n1*self.backgroundPDF_uncert2 + nsig*nsig*self.signalPDF_uncert2
            
            bins_to_use = (modelPDF>0.)&(modelPDF_uncert2>0.) 

            alpha = modelPDF[bins_to_use]**2/modelPDF_uncert2[bins_to_use] +1.
            beta  = modelPDF[bins_to_use]/modelPDF_uncert2[bins_to_use]

            values = [
              alpha*np.log(beta),
              sps.loggamma(self.observation[bins_to_use]+alpha).real,
              -(self.observation[bins_to_use]+alpha)*np.log1p(beta),
              -sps.loggamma(alpha).real,
            ]

        else:
            raise ValueError('No valid LLH type defined!')

        return -np.sum(values)

    
    def ComputeBestFit(self):
        LLHmin_DM=Minuit(self.evaluateLLH,
             nsig=1.,n1=1.,
             error_nsig=.01,error_n1=.01,
             limit_nsig=(-10.,100.),limit_n1=(0.,10.),
             errordef=.5,print_level=0)  
        LLHmin_DM.migrad()
        
        self.bestFit = {}
        self.bestFit['n1']=LLHmin_DM.fitarg['n1']
        self.bestFit['nsig']=LLHmin_DM.fitarg['nsig']
        self.bestFit['LLH']=self.evaluateLLH(self.bestFit['n1'],self.bestFit['nsig'])
                
        self.computedBestFit = True
        
        
    def ComputeTestStatistics(self):
        
        if not self.computedBestFit:
            self.ComputeBestFit()

        LLHmin_ref=Minuit(self.evaluateLLH,
             nsig=0.,n1=1.,
             fix_nsig = True,
             error_nsig=.1,error_n1=.1,
             limit_nsig=(-10.,100.),limit_n1=(0.,10.),
             errordef=.5,print_level=0)  
        LLHmin_ref.migrad()
            
        self.TS = 0.

        self.bestFit['LLH_ref'] = self.evaluateLLH(LLHmin_ref.fitarg['n1'],LLHmin_ref.fitarg['nsig'])
        if self.bestFit['nsig'] > 0.:
            self.TS = 2*(self.bestFit['LLH_ref']-self.bestFit['LLH'])
        

    def CalculateUpperLimit(self,conf_level):

        nIterations = 0
        eps_TS=0.005
        eps_param=0.0005

        deltaTS = 2.71
        if conf_level==90:
            deltaTS = 1.64
        elif conf_level==95:
            deltaTS = 2.71
            
        param_low=self.bestFit['nsig']
        param_up=self.bestFit['nsig']
        param_mean=self.bestFit['nsig']
        
        dTS=0
        cc=1
        while((dTS<deltaTS) and (nIterations<100)):
            nIterations += 1 

            param_up=param_up+3.*np.abs(param_up)

            LLHmin_fix=Minuit(self.evaluateLLH,
                 nsig=param_up,fix_nsig = True,
                 n1=1.,
                 error_nsig=.1,error_n1=.1,
                 limit_nsig=(-10.,100.),limit_n1=(0.,10.),
                 errordef=.5,print_level=0)  
            LLHmin_fix.migrad()

            if param_up <0.:
                TS_fix = 0.
            else:
                TS_fix = 2*(self.bestFit['LLH_ref']-self.evaluateLLH(LLHmin_fix.fitarg['n1'],param_up))

            dTS = self.TS - TS_fix

        nIterations = 0
        param_low=param_up/4.
        while((cc>0.)  and (nIterations<100)):
            
            nIterations += 1

            param_mean=(param_low+param_up)/2.
            LLHmin_fix=Minuit(self.evaluateLLH,
                 nsig=param_mean,fix_nsig = True,
                 n1=1., error_nsig=.1,error_n1=.1,
                 limit_nsig=(-10.,100.),limit_n1=(0.,10.),
                 errordef=.5,print_level=0)  
            LLHmin_fix.migrad()

            if param_mean <0.:
                TS_fix = 0.
            else:
                TS_fix = 2*(self.bestFit['LLH_ref']-self.evaluateLLH(LLHmin_fix.fitarg['n1'],param_mean))
                
            dTS = self.TS - TS_fix
            
            if(dTS<deltaTS):

                param_low=param_mean
                delta_param=(param_up-param_low)/(param_up)
                
                if((dTS>deltaTS-eps_TS) and (delta_param < eps_param)):
                    cc = 0
                    
            if(dTS>deltaTS):
                param_up=param_mean
                delta_param=(param_up-param_low)/(param_up)
                
                if((dTS<deltaTS+eps_TS) and (delta_param < eps_param)):
                    cc=0
                    
        return param_up
       
    
    def CalculateSensitivity(self,nTrials, conf_level):

        if self.LLHtype == None:
            raise ValueError('LLH type not defined yet!')

        TS = []
        upperlimits = []
        if self.moreOutput:
            fits = []
        
        for i in range(nTrials):
            self.sampleObservation(1.,0.)
            self.ComputeTestStatistics()
            TS.append(self.TS)

            ul = self.CalculateUpperLimit(conf_level)
            if np.isnan(ul):
                print("Warning: NaN upper limit at trial {i}.\nRepeating trial.".format(i=i))
                i-=1
                continue
            upperlimits.append(ul)
            
            if self.moreOutput:
                fits.append(self.bestFit)
            
        p_median = np.percentile(upperlimits, 50)
        p_95_low = np.percentile(upperlimits, 2.5)
        p_95_high = np.percentile(upperlimits, 97.5)
        p_68_low = np.percentile(upperlimits, 16.)
        p_68_high = np.percentile(upperlimits, 84.)

        dic_brazilian = {}
        dic_brazilian['TS_dist'] = TS
        dic_brazilian['error_68_low'] = p_68_low
        dic_brazilian['error_68_high'] = p_68_high
        dic_brazilian['error_95_low'] = p_95_low
        dic_brazilian['error_95_high'] = p_95_high   
        dic_brazilian['median'] = p_median
        if self.moreOutput:
            dic_brazilian['upperlimits'] = upperlimits
            dic_brazilian['bestFits'] = fits

        return dic_brazilian
    
    
    
    
    
    
class Profile_Analyser_Normalised:

    def __init__(self):
        
        self.LLHtype = None

        self.ready = False
        self.signalPDF = None
        self.backgroundPDF = None
        
        #This is our assumed baseline background PDF
        #This is only filled once, and does not get updated when generating trials
        self.baseline_backgroundPDF = None 
        self.baseline_init = False

        self.signalPDF_uncert2 = None
        self.backgroundPDF_uncert2 = None

        self.signalContaminationPDF = None
        self.signalContaminationPDF_uncert2 = None
        
        self.nTotalEvents = 0.
        self.nSignalEvents = 0.
        
        self.nbins = 0
        
        self.observation = None

        self.computedBestFit = False
        self.bestFit = None
        self.TS = None
        
        self.moreOutput = False
        self.AllowNegativeSignal = False

    def setLLHtype(self,type):
        availableTypes = ['Poisson', 'Effective', 'PoissonWithSignalSubtraction', 'EffectiveWithSignalSubtraction']
        if type not in availableTypes:
            raise ValueError('LLH type not implemented yet. Choose amongst: '+str(availableTypes))
        else:
            self.LLHtype = type
            
    def saveMoreOutput(self):
        self.moreOutput = True

    def allowNegativeSignal(self):
        self.AllowNegativeSignal = True        

    def loadBackgroundPDF(self,pdf, verbose = False):
        self.backgroundPDF = pdf.flatten()/np.sum(pdf)
        self.nTotalEvents = np.sum(pdf)
        if verbose:
            print('Total number of expected background events:', self.nTotalEvents)
        if not self.baseline_init:
            print("Initalizing the baseline pdf")
            self.baseline_backgroundPDF = self.backgroundPDF
            self.baseline_init = True
            
        self.nbins = len(pdf.flatten())

    def loadSignalPDF(self,pdf):
        self.signalPDF = pdf.flatten()/np.sum(pdf)
        self.nSignalEvents = np.sum(pdf)
        if self.nbins == len(pdf.flatten()):
            self.ready = True
        else:
            raise ValueError('Shape of signal pdf does not match the background pdf! Did you initialize the background pdf first?')
    
    def loadUncertaintyPDFs(self,bkg_pdf,sig_pdf):
        self.backgroundPDF_uncert2 = bkg_pdf.flatten()/self.nTotalEvents
        self.signalPDF_uncert2 = sig_pdf.flatten()/self.nSignalEvents

        if self.nbins != len(bkg_pdf.flatten()):
            raise ValueError('Shape of background uncertainty pdf does not match the background pdf!')
        if self.nbins != len(sig_pdf.flatten()):
            raise ValueError('Shape of signal uncertainty pdf does not match the background pdf!')
            

    def loadSignalScrambledPDF(self,pdf):
        #self.signalContaminationPDF = pdf.flatten()/self.nSignalEvents
        #This is a pdf so we should normalize it as such
        self.signalContaminationPDF = pdf.flatten()/np.sum(pdf)
        if self.nbins != len(pdf.flatten()):
            raise ValueError('Shape of signal contamination pdf does not match the background pdf!')
    
    def loadSignalScrambledUncertaintyPDF(self,pdf):
        #TO be check.
        self.signalContaminationPDF_uncert2 = pdf.flatten()/self.nSignalEvents
        if self.nbins != len(pdf.flatten()):
            raise ValueError('Shape of signal contamination uncertainty pdf does not match the background pdf!')
        
    def sampleObservation(self,xi):

        if not self.ready:
            raise ValueError('Not all pdfs are correctly loaded!')
        #We use always the baseline PDF, this does not get updated
        observationPDF = self.nTotalEvents* ((1-xi)*self.baseline_backgroundPDF + xi*self.signalPDF)

        self.observation=np.zeros(np.shape(self.backgroundPDF))
        for i in range(len(self.observation)):
            self.observation[i]=np.random.poisson(observationPDF[i])
        
        self.computedBestFit = False
        
        #So once we have a new dataset, our background estimate has to change! as the scrambled background is taken from data
        #if (self.LLHtype == 'PoissonWithSignalSubtraction' or self.LLHtype == 'EffectiveWithSignalSubtraction'):
        #So this is a trick, in reality we don't use the exact pseudosample, that will require to get individual events and re-scrambled their RA. Instead we are going to construct a PDF using the scrambled signal.
        scrambledPDF = self.nTotalEvents* ((1-xi)*self.baseline_backgroundPDF + xi*self.signalContaminationPDF)
        self.loadBackgroundPDF(scrambledPDF)
        
        
    def evaluateLLH(self, xi):
        
        if self.LLHtype == 'Poisson':

            modelPDF = self.nTotalEvents*((1-xi)*self.backgroundPDF + xi*self.signalPDF)
            
            if np.isnan(modelPDF).any():
                print('nan in model array with xi=',xi, self.computedBestFit)

            bins_to_use = (modelPDF>0.)
            values = self.observation[bins_to_use]*np.log(modelPDF[bins_to_use])-modelPDF[bins_to_use]
            
            
        elif self.LLHtype == 'Effective':
        
            modelPDF = self.nTotalEvents*((1-xi)*self.backgroundPDF + xi*self.signalPDF)
            modelPDF_uncert2 = self.nTotalEvents*((1-xi)*(1-xi)*self.backgroundPDF_uncert2 + xi*xi*self.signalPDF_uncert2)

            if np.isnan(modelPDF).any():
                print('nan in model array with xi=',xi, self.computedBestFit)
                
            bins_to_use = (modelPDF>0.)&(modelPDF_uncert2>0.) 

            alpha = modelPDF[bins_to_use]**2/modelPDF_uncert2[bins_to_use] +1.
            beta  = modelPDF[bins_to_use]/modelPDF_uncert2[bins_to_use]

            values = [
              alpha*np.log(beta),
              sps.loggamma(self.observation[bins_to_use]+alpha).real,
              -(self.observation[bins_to_use]+alpha)*np.log1p(beta),
              -sps.loggamma(alpha).real,
            ]
            

        elif self.LLHtype == 'PoissonWithSignalSubtraction':

            modelPDF = self.nTotalEvents*(self.backgroundPDF + xi*self.signalPDF - xi*self.signalContaminationPDF)
            
            if np.isnan(modelPDF).any():
                print('nan in model array with xi=',xi, self.computedBestFit)

            bins_to_use = (modelPDF>0.)
            values = self.observation[bins_to_use]*np.log(modelPDF[bins_to_use])-modelPDF[bins_to_use]
           
        elif self.LLHtype == 'EffectiveWithSignalSubtraction':
        
            modelPDF = self.nTotalEvents*(self.backgroundPDF + xi*self.signalPDF- xi*self.signalContaminationPDF)
            modelPDF_uncert2 = self.nTotalEvents*(self.backgroundPDF_uncert2 + xi*xi*self.signalPDF_uncert2 - xi*xi*self.signalContaminationPDF_uncert2)

            if np.isnan(modelPDF).any():
                print('nan in model array with xi=',xi, self.computedBestFit)
                
            bins_to_use = (modelPDF>0.)&(modelPDF_uncert2>0.) 

            alpha = modelPDF[bins_to_use]**2/modelPDF_uncert2[bins_to_use] +1.
            beta  = modelPDF[bins_to_use]/modelPDF_uncert2[bins_to_use]

            values = [
              alpha*np.log(beta),
              sps.loggamma(self.observation[bins_to_use]+alpha).real,
              -(self.observation[bins_to_use]+alpha)*np.log1p(beta),
              -sps.loggamma(alpha).real,
            ]
        else:
            raise ValueError('No valid LLH type defined!')
        
        return -np.sum(values)

    
    def ComputeBestFit(self):
        
        if self.AllowNegativeSignal:
            lower_bound = -1.
        else:
            lower_bound = 0.
        
        LLHmin_DM=Minuit(self.evaluateLLH,
             xi=0.1, error_xi=.01,
             limit_xi=(lower_bound,2.),
             errordef=.5,print_level=0)  
        LLHmin_DM.migrad()
        
        self.bestFit = {}
        self.bestFit['xi']=LLHmin_DM.fitarg['xi']
        self.bestFit['LLH']=self.evaluateLLH(self.bestFit['xi'])
                
        self.computedBestFit = True
        
        
    def ComputeTestStatistics(self):
        
        if not self.computedBestFit:
            self.ComputeBestFit()

        self.TS = 0.
        self.bestFit['LLH_ref'] = self.evaluateLLH(0.)
        
        if self.bestFit['xi'] > 0.:
            self.TS = 2*(self.bestFit['LLH_ref']-self.bestFit['LLH'])
            
        #This means that somehow the minimizer didn't reach the best value (which should be xi = 0)
        if self.TS < 0:
            self.TS = 0
        #    self.bestFit['xi'] = 0
    
   # def CalculateFrequentisLimit(self, TS, conf_level): 
        
        
    def CalculateUpperLimit(self,conf_level):

        nIterations = 0
        eps_TS=0.005
        eps_param=0.0005

        deltaTS = 2.71
        if conf_level==90:
            deltaTS = 1.64
        elif conf_level==95:
            deltaTS = 2.71
            
        param_low=self.bestFit['xi']
        param_up=self.bestFit['xi']
        param_mean=self.bestFit['xi']
        
        dTS=0
        cc=1
        while((dTS<deltaTS) and (nIterations<100)):
            nIterations += 1 

            param_up=param_up+3.*np.abs(param_up)

            if param_up <0.:
                TS_fix = 0.
            else:
                TS_fix = 2*(self.bestFit['LLH_ref']-self.evaluateLLH(param_up))
                            
            dTS = self.TS - TS_fix

        nIterations = 0
        param_low=param_up/4.
        while((cc>0.)  and (nIterations<100)):
            
            nIterations += 1

            param_mean=(param_low+param_up)/2.
            if param_mean <0.:
                TS_fix = 0.
            else:
                TS_fix = 2*(self.bestFit['LLH_ref']-self.evaluateLLH(param_mean))
                
            dTS = self.TS - TS_fix
            
            if(dTS<deltaTS):

                param_low=param_mean
                delta_param=(param_up-param_low)/(param_up)
                   
                
                if((dTS>deltaTS-eps_TS) and (delta_param < eps_param)):
                    cc = 0
                    
            if(dTS>deltaTS):
                param_up=param_mean
                delta_param=(param_up-param_low)/(param_up)
                
                if((dTS<deltaTS+eps_TS) and (delta_param < eps_param)):
                    cc=0
                    
        return param_up
       
        
        
    def DoScan(self, ni_min, ni_max, nstep, ts, nTrials, conf_level):
        print(" Scanning injected fraction of events range [%.6f, %.6f]" % (ni_min, ni_max))
        results = []
        step = 0
        frac = 0.
        tolerance = 0.1
        ni_mean = 0
        while (np.abs(frac*100 - conf_level) > tolerance):
            ni_mean = (ni_max + ni_min)/2
          
            TS = []
            for i in range(nTrials):
                self.sampleObservation(ni_mean)
                self.ComputeTestStatistics()
                TS.append(self.TS)

            TS = np.array(TS)
            n = TS[np.where(TS > ts)].size
            ntot = TS.size
            frac = n/ntot
            #frac_error = BinomialError(ntot, n)/ntot

            print(" [%2d] ni %5.6f, n %4d, ntot %4d, C.L %.2f +/- %.2f" %
                (step, ni_mean, n, ntot, frac, 0))

            if(frac * 100 < conf_level):
                ni_min = ni_mean
            else:
                ni_max = ni_mean
           
            results.append([ni_mean, n, ntot])
            
        

        return ni_mean, results    

    def CalculateNeymanSensitivity(self, nTrials, conf_level):
        #Let's first guess using the likelihood intervales
        sens = self.CalculateSensitivity(nTrials, conf_level)
        median_ts = np.percentile(sens['TS_dist'], 50)
        first_guess = sens['median']
        p = conf_level/100
        factor = 1.
        ni_min = first_guess / factor * (1 - p)
        ni_max = first_guess * factor / p
        
        return self.DoScan(ni_min, ni_max, 50, median_ts, nTrials, conf_level)
        
    
    def CalculateSensitivity(self,nTrials, conf_level, doNeyman=False, nTrialsNeyman=100):

        if self.LLHtype == None:
            raise ValueError('LLH type not defined yet!')

        TS = []
        upperlimits = []
        if self.moreOutput:
            fits = []
        
        for i in range(nTrials):
            self.sampleObservation(0.)
            self.ComputeTestStatistics()
            TS.append(self.TS)

            ul = self.CalculateUpperLimit(conf_level)
            if np.isnan(ul):
                print("Warning: NaN upper limit at trial {i}.\nRepeating trial.".format(i=i))
                i-=1
                continue
            upperlimits.append(ul)
            
            if self.moreOutput:
                fits.append(self.bestFit)
            
        p_median = np.percentile(upperlimits, 50)
        p_95_low = np.percentile(upperlimits, 2.5)
        p_95_high = np.percentile(upperlimits, 97.5)
        p_68_low = np.percentile(upperlimits, 16.)
        p_68_high = np.percentile(upperlimits, 84.)

        dic_brazilian = {}
        dic_brazilian['TS_dist'] = TS
        dic_brazilian['error_68_low'] = p_68_low
        dic_brazilian['error_68_high'] = p_68_high
        dic_brazilian['error_95_low'] = p_95_low
        dic_brazilian['error_95_high'] = p_95_high   
        dic_brazilian['median'] = p_median
        if self.moreOutput:
            dic_brazilian['upperlimits'] = upperlimits
            dic_brazilian['bestFits'] = fits
            
            
        if doNeyman:
            print('Adding Neyman sensitivity now')
            first_guess = p_median
                        
            median_TS = np.median(TS)
            if median_TS<0.:
                median_TS = 0.
                
            xi_scan_min = first_guess / 2. * (1 - (conf_level/100.))
            xi_scan_max = first_guess * 2. / (conf_level/100.)
                               
            nIterations = 0
            eps_CL = 0.01
            delta_CL = 1.
            while delta_CL>eps_CL:
                nIterations += 1
                xi_scan_mean=(xi_scan_min+xi_scan_max)/2.
                
                tmp_TS = []
                for i in range(nTrialsNeyman):
                    self.sampleObservation(xi_scan_mean)
                    self.ComputeTestStatistics()
                    tmp_TS.append(self.TS)
                tmp_TS = np.array(tmp_TS)

                n_aboveThreshold = tmp_TS[np.where(tmp_TS > median_TS)].size
                n_tot = tmp_TS.size
                fraction_aboveThreshold = float(n_aboveThreshold)/float(n_tot)
                #frac_aboveThreshold_error = BinomialError(n_tot, n_aboveThreshold)/n_tot
                
                delta_CL = np.abs(conf_level/100.-fraction_aboveThreshold)
                
                if(fraction_aboveThreshold<conf_level/100.):
                    xi_scan_min=xi_scan_mean
                else:
                    xi_scan_max=xi_scan_mean
                                
            dic_brazilian['Neyman_sensitivity'] = xi_scan_mean
            
        return dic_brazilian   
    
    
    def CalculateDiscoveryPotential(self, significance):

        signalStrength = 0.
        TS = 0.
       
        while TS<significance**2:
            tmp_TS = []
            signalStrength += 0.1*self.bestFit['xi']
            for i in range(100):
                self.sampleObservation(signalStrength)
                self.ComputeTestStatistics()
                tmp_TS.append(self.TS)
            TS = np.median(TS)
            
        return signalStrength
        
        
        
        