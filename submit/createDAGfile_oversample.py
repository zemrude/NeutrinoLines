systematics = ['nominal']

jobid = ""

counter = 0
for sys in systematics:
    oversampling = 200
    jobid = "%.3i-oversampling-LE"%(counter)  
    print("JOB " + jobid + " OVERSAMPLER.submit")
    print("VARS " + jobid + " JOBNAME=\"%s\" SYST=\"%s\" LECUT=\"0.15\" HECUT=\"0.2\"  OVERSAMPLING=\"%i\""%(jobid, sys, oversampling))
    jobid = "%.3i-oversampling-HE"%(counter)  
    print("JOB " + jobid + " OVERSAMPLER.submit")
    print("VARS " + jobid + " JOBNAME=\"%s\" SYST=\"%s\" LECUT=\"-1.0\" HECUT=\"0.3\"  OVERSAMPLING=\"%i\""%(jobid, sys, oversampling))