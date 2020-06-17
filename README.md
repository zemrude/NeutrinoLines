# Analysis-IC86-BSM-DarkMatter-NeutrinoLine
Analysis code for the search of decaying and annihilating dark matter with focus on neutrino lines

## Code/File structure:

- README.md : this file
- env.sh : file to define environment variables needed by the scripts

data					folder containing numpy files with the processed data and simulations
PDFs					folder containing all PDF in pickle format
sensitivity				folder containing calculated sensitivities in pickle format
resources				folder containing limits by other experiments and the PPPC spectra tables

python		 			various python code for the analysis
	fileList.py			provides filenames and paths for the data files
	IceCube_sim.py			functions related to IceCube simulations
	PPPC_spectra.py			functions to load and interpolate PPPC spectra 
	fluxCalculation.py		helper functions for DM and atmospheric flux calculations
	createPDFs.py			code for the calculation of PDFs
					parameters:
						-t type (background/annihilation/decay)
						--lecut cut on LE BDT score
						--hecut cut on HE BDT score
						-s systematics
						-o oversampling factor
						-c annihilation/decay channel (signal only)
						-m DM mass (signal only)
						-p DM profile NFW/Burkert (signal only)
	SensitivityCalucalation.py	code for sensitivity calculation
	                                parameters:
						-t type (annihilation/decay)
						--lecut cut on LE BDT score
						--hecut cut on HE BDT score
						-s systematics
						-o oversampling factor
						-c annihilation/decay channel
						-m DM mass
						--profile DM profile NFW/Burkert
						-l confidence level
						-d likelihood method (poisson or effective)
						-e energy rebinning
						-p psi rebinning
plotting				python plotting scripts
submit					sample scripts for condor/dagman submission


## Prerequisites

- icerec build, here from /data/user/sbaur/metaprojects/icerec/build/
- python 2, from /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1
- python modules: iminuit, scipy, astropy, OptionParser

## Workflow / examples

1) PDF creation for signal and background, specific for each channel and mass
   example for 1 TeV annihilation to nuE, HE sample:
		./createPDFs.py -t annihilation -p NFW -c nue -m 1000 --lecut -1.0 --hecut 0.3 -s nominal -o 100
		./createPDFs.py -t background --lecut -1.0 -lhecut 0.3 -s nominal

2) Calculate sensitivity for each channel and mass
   example for 1 TeV annihilation to nuE, HE sample:
		./SensitivityCalculation.py -t annihilation --lecut -1.0 --hecut 0.3 -s nominal -o 100 -c nue -m 1000 --profile NFW -l 90 -d effective -e 2 -p 5

3) Unblinded data:
   tbd...

4) Plotting:
   tbd...

