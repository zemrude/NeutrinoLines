masses = [40, 63, 100, 158, 251, 398, 631, 1000, 1585, 2512, 3981, 6310, 10000, 15850, 25120, 39810]

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
                for llh in likelihood:
                    for mass in masses:
                        
                        jobid = "%.3i-sens-%s-%s-%s-%s-%s"%(counter, str(channel), str(mass), str(profile), str(mode), str(llh))  
                        print("JOB " + jobid + " SENSITIVITY.submit")
                        print("VARS " + jobid + " JOBNAME=\"%s\" TYPE=\"%s\" CHANNEL=\"%s\" PROFILE=\"%s\" SYST=\"%s\" MASS=\"%i\" LLH=\"%s\""%(jobid, mode, channel, profile, sys, mass, llh))
                        counter += 1