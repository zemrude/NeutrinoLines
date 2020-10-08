#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /data/user/aguilar/icetray/meta-projects/combo/V01-00-02/build

from time import gmtime, strftime
import numpy as np
import time
import random
import astropy
from astropy.time import Time as apTime
import sys, os
import math
from icecube import astro
import pickle
from optparse import OptionParser
import env

import importlib
halo_spec = importlib.util.find_spec("halo")
halo_found = halo_spec is not None

# This is to make the waiting nicer

if halo_found:
    from halo import Halo

    
base_path = os.environ['ANALYSIS_BASE_PATH']
sys.path.append(base_path+'python')
from utils import psi_f

r"""

  This script generates the signal PDFs both for the true signal and it's scrambled version as well as the uncertainties**2

"""

