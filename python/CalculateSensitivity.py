#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /data/user/aguilar/icetray/meta-projects/combo/V01-00-02/build

import sys, os
from iminuit import Minuit
import numpy as np; import time; from timeout_decorator import timeout,TimeoutError
import scipy.special as sps
from optparse import OptionParser
from scipy.optimize import fsolve, root
import pickle

import env


base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')

r"""

Get and define parameters

"""
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-t", "--type",default="annihilation",
                  dest="TYPE", help="Define type of PDF")
parser.add_option("-c", "--channel",default="nue",
                  dest="CHANNEL", help="Annihilation channel")
parser.add_option("-x", "--profile",default="NFW",
                  dest="PROFILE", help="DM Profile")
parser.add_option("-s", "--syst",default="nominal",
                  dest="SYST", help="Define systematic variation")
parser.add_option("-a", "--lecut", default='0.2',
                  dest="LECUT", help="Cut on LE BDT")
parser.add_option("-b", "--hecut", default='0.2',
                  dest="HECUT", help="Cut on HE BDT")
parser.add_option("-m", "--mass", default='0',
                  dest="MASS", help="DM mass")
parser.add_option("-o", "--oversampling", default='1',
                  dest="OVERSAMPLING", help="Oversampling factor")
parser.add_option("-e", "--energy-binning", default='2',
                  dest="REBINE", help="rebin factor energy")
parser.add_option("-p", "--psi-binning", default='5',
                  dest="REBINPSI", help="rebin factor psi")
parser.add_option("-l", "--conf-level", default='90',
                  dest="CONFLEVEL", help="confidence level")
parser.add_option("-d", "--llh", default='Poisson',
                  dest="LLH", help="LLH type")
parser.add_option("-z", "--signal", default='0', # either sensitivity or xs value
                  dest="SIGNAL", help="Amount of signal contamination")
parser.add_option("-y", "--bkg", default='simulation', # either simulation or data
                  dest="BKG", help="Background")

(options,args) = parser.parse_args()
        
mode = options.TYPE #annihilation, decay, background
channel = options.CHANNEL
profile = options.PROFILE
systematics = options.SYST

LECut = float(options.LECUT)
HECut = float(options.HECUT)
mass = int(options.MASS)
nOversampling = int(options.OVERSAMPLING)

bins_merge_E   = int(options.REBINE)
bins_merge_Psi  = int(options.REBINPSI)

conf_level = int(options.CONFLEVEL)

LLH_type = options.LLH
