LEmasses = [40, 63, 100, 158, 251, 398, 631, 1000, 1585]
HEmasses = [1000, 1585, 2512, 3981, 6310, 10000, 15850, 25120, 39810]

DM_mode = ['annihilation']#, 'decay']
DM_profiles = ['NFW']#, 'Burkert']
DM_channel = ['nue' ] #,'b', 'W', 'tau', 'mu']      

systematics = ['nominal']
# 'DomEffUp', 'DomEffDown', 'Ice_HoleIce100_ScatAbs-7' 'Ice_HoleIce100_Scat+10', 'Ice_HoleIce100_Abs+10', 'Ice_HoleIce30_ScatAbs-7', 'Ice_HoleIce30_Scat+10', 'Ice_HoleIce30_Abs+10', 'nominalGammaUp' 'nominalGammaDown']     

likelihood = ['Poisson', 'Effective', 'PoissonWithSignalSubtraction', 'EffectiveWithSignalSubtraction']

jobid = ""

counter = 0
for sys in systematics:
    for channel in DM_channel:
        for profile in DM_profiles:
            for mode in DM_mode:
                for llh in likelhood:
                    for mass in LEmasses:
                        if channel == 'nue':
                            oversampling = 100
                        else:
                            oversampling = -1

                        jobid = "%.2i-sens-%s-%s-%s-%s-LE"%(counter,str(channel), str(mass), str(profile), str(mode))  
                        print("JOB " + jobid + " PDF_SIGNAL.submit")
                        print("VARS " + jobid + " JOBNAME=\"%s\" TYPE=\"%s\" CHANNEL=\"%s\" PROFILE=\"%s\" SYST=\"%s\" LECUT=\"0.15\" HECUT=\"0.2\" MASS=\"%i\" OVERSAMPLING=\"%i\" LLH=\"%s\""%(jobid, mode, channel, profile,sys, mass, oversampling, llh))
                    for mass in HEmasses:
                        if channel == 'nue':
                            oversampling = 200
                        else:
                            oversampling = -1

                        jobid = "%.2i-sens-%s-%s-%s-%s-HE"%(counter,str(channel), str(mass), str(profile), str(mode))  
                        print("JOB " + jobid + " PDF_SIGNAL.submit")
                        print("VARS " + jobid + " JOBNAME=\"%s\" TYPE=\"%s\" CHANNEL=\"%s\" PROFILE=\"%s\" SYST=\"%s\" LECUT=\"-1.0\" HECUT=\"0.3\" MASS=\"%i\" OVERSAMPLING=\"%i\" LLH=\"%s\""%(jobid, mode, channel, profile,sys, mass, oversampling, LLH))
    counter += counter