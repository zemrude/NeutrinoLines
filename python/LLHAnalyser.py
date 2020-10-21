import numpy as np
import scipy.special as sps
from iminuit import Minuit
from utils import ConfidenceIntervalError

from sensitivity_utils import BinomialError, inv_BinomialError

class sensitivity:
    def __init__(self):
        self.xi = 0 # 
        self.sv = 0
        self.error68_low = 0
        self.error68_high = 0
        self.error95_low = 0
        self.error95_high = 0
        self.mass = 0
        
     
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
      
        if removeEmptyBins:
            print("We are removing empty bins")
            
            
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
       
    
    def CalculateSensitivity(self, nTrials, conf_level):

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
    
    
    
    
    
    
class Likelihood_Analyser:

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
        print(" Allowing for negative signal") 

    def loadBackgroundPDF(self, pdf, verbose = False):
        self.backgroundPDF = pdf.flatten()    
        self.nTotalEvents = np.sum(pdf)
        self.backgroundPDF = pdf.flatten()/np.sum(pdf)

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
        #This are sigma**2 hence, the normalization also needs to be squared
        self.backgroundPDF_uncert2 = bkg_pdf.flatten()/self.nTotalEvents**2
        self.signalPDF_uncert2 = sig_pdf.flatten()/self.nSignalEvents**2

        if self.nbins != len(bkg_pdf.flatten()):
            raise ValueError('Shape of background uncertainty pdf does not match the background pdf!')
        if self.nbins != len(sig_pdf.flatten()):
            raise ValueError('Shape of signal uncertainty pdf does not match the background pdf!')
            

    def loadSignalScrambledPDF(self,pdf):
        #We assume that the scrambled PDF has the same number of events as 
        #if (np.sum(pdf) != self.nSignalEvents):
        #    raise ValueError(' Normalization of scrambled signal does not match the signal pdf!')
        self.signalContaminationPDF = pdf.flatten()/self.nSignalEvents
        #This is a pdf so we should normalize it as such
        #self.signalContaminationPDF = pdf.flatten()/np.sum(pdf)
        if self.nbins != len(pdf.flatten()):
            raise ValueError('Shape of signal contamination pdf does not match the background pdf!')
    
    def loadSignalScrambledUncertaintyPDF(self,pdf):
        #TO be check.
        self.signalContaminationPDF_uncert2 = pdf.flatten()/self.nSignalEvents**2
        
        if self.nbins != len(pdf.flatten()):
            raise ValueError('Shape of signal contamination uncertainty pdf does not match the background pdf!')
        
    def sampleObservation(self,xi):

        if not self.ready:
            raise ValueError('Not all pdfs are correctly loaded!')
        #We use always the baseline PDF, this does not get updated
        observationPDF = self.nTotalEvents* ((1-xi)*self.baseline_backgroundPDF + xi*self.signalPDF)

        self.observation=np.zeros(np.shape(self.backgroundPDF))
        for i in range(len(self.observation)):
            if observationPDF[i] == 0:
                self.observation[i]=np.random.poisson(1.02) # Poisson mean with 90% below 2.3
            else:
                self.observation[i]=np.random.poisson(observationPDF[i])
        
        self.computedBestFit = False
        
        r"""
        Once we have a new dataset, our background estimate has to change! as the scrambled background is taken from data
        The following is an approximation: In reality we should the invidual events of the pseudosample and scrambled their RA. Instead we are going to construct a PDF using the scrambled signal.
        """
        
        scrambledPDF = self.nTotalEvents* ((1-xi)*self.baseline_backgroundPDF + xi*self.signalContaminationPDF)
        self.loadBackgroundPDF(scrambledPDF)
        
        
    def evaluateLLH(self, xi):
        
        if self.LLHtype == 'Poisson':

            modelPDF = self.nTotalEvents*((1-xi)*self.backgroundPDF + xi*self.signalPDF)
            
            if np.isnan(modelPDF).any():
                print('nan in model array with xi=',xi, self.computedBestFit)

            bins_to_use = (modelPDF>0.)
            r"""
            Poisson likelihood
            
            logL = sum(k * log lamnda - lambda)
            
            Factorials are kept out as we usually care for DeltaLLH
            
            """
            values = self.observation[bins_to_use]*np.log(modelPDF[bins_to_use])-modelPDF[bins_to_use]
            
            
        elif self.LLHtype == 'Effective':
        
            modelPDF = self.nTotalEvents*((1-xi)*self.backgroundPDF + xi*self.signalPDF)
            #modelPDF_uncert2 = self.nTotalEvents*((1-xi)*(1-xi)*self.backgroundPDF_uncert2 + xi*xi*self.signalPDF_uncert2)
            #I think this is wrong, the uncert2 needs to have ntot**2!
            #But then uncert2 need to be have normalized as n**2
            
            modelPDF_uncert2 = self.nTotalEvents**2*((1-xi)*(1-xi)*self.backgroundPDF_uncert2 + xi*xi*self.signalPDF_uncert2)

            if np.isnan(modelPDF).any():
                print('nan in model array with xi=',xi, self.computedBestFit)
                
            bins_to_use = (modelPDF>0.)&(modelPDF_uncert2>0.) 

            alpha = modelPDF[bins_to_use]**2/modelPDF_uncert2[bins_to_use] + 1.
            beta  = modelPDF[bins_to_use]/modelPDF_uncert2[bins_to_use]

            r"""
            Effective likelihood from:
            ArXiv ePrint:1901.04645
            Formula 3.15
            
            P = b**a * Gamma(k+a)/ (k! (1+b)^(k+a) Gamma(a))
            
            Factorials are kept out as we usually care for DeltaLLH
             
            log L = sum(a * log(b) + logGamma(k + a) - 
            (k + a)*log(1+b) - logGamma(a))
            """
            
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
            
            modelPDF_uncert2 = self.nTotalEvents**2*(self.backgroundPDF_uncert2 + xi*xi*self.signalPDF_uncert2 - xi*xi*self.signalContaminationPDF_uncert2)

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
        
        #As defined in the asimov paper, otherwise is 0
        if self.bestFit['xi'] > 0.:
            self.TS = 2*(self.bestFit['LLH_ref']-self.bestFit['LLH'])
            
        #This means that somehow the minimizer didn't reach the best value (ie is lower than xi = 0, in that case we take xi = 0)
        if self.TS < 0:
            self.TS = 0
    
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
       
        
        
    def DoScan(self, xi_min, xi_max, nstep, ts, conf_level, precision):
        print(" Scanning injected fraction of events range [%.6f, %.6f]" % (xi_min, xi_max))
        results = []
        step = 0
        frac = 0.
        tolerance = 0.1
        xi_mean = 0
        frac_error = 0 
        
        nTrials = inv_BinomialError(precision, conf_level)
        
        print (" Doing %i trials for %i +/- %.1f C.L." %(nTrials, conf_level, precision))
        
        p = conf_level / 100
        while (np.abs(frac - p) > frac_error) and (step < nstep):
            xi_mean = (xi_max + xi_min)/2
            dic_results = {}
            
            TS = []
            for i in range(nTrials):
                self.sampleObservation(xi_mean)
                self.ComputeTestStatistics()
                TS.append(self.TS)

            TS = np.array(TS)
            n = TS[np.where(TS > ts)].size
            ntot = TS.size
            frac = n/ntot
            frac_error = BinomialError(ntot, n)/ntot

            print(" [%2d] xi %5.6f, n %4d, ntot %4d, C.L %.2f +/- %.2f" %
                (step, xi_mean, n, ntot, frac * 100, frac_error * 100))

            if(frac < p):
                xi_min = xi_mean
            else:
                xi_max = xi_mean

            dic_results['xi'] = xi_mean
            dic_results['n'] = n
            dic_results['ntrials'] = ntot
            dic_results['TS_dis'] = TS
            step += 1    
            results.append(dic_results)
                    

        return xi_mean, results    

    def CalculateFrequentistSensitivity(self,  conf_level, precision):
        #Let's first guess using the likelihood intervales
        #To calculate the median a 100 trials seem enough
        sens = self.CalculateSensitivity(100, conf_level)
        median_ts = np.percentile(sens['TS_dist'], 50)
        first_guess = sens['median']
        p = conf_level/100
        factor = 2.
        
        ni_min = first_guess / factor * (1 - p)
        ni_max = first_guess * factor / p
        
        print (" Median TS: %.2f" %median_ts)
        
        return self.DoScan(ni_min, ni_max, 50, median_ts,  conf_level, precision)
        
    
    def CalculateSensitivity(self, nTrials, conf_level):

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
                print(" Warning: NaN upper limit at trial {i}.\nRepeating trial.".format(i=i))
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
        
        
        
        