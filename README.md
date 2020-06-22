# Analysis-IC86-BSM-DarkMatter-NeutrinoLine
Analysis code for the search of decaying and annihilating dark matter with focus on neutrino lines

Location: `/data/ana/BSM/HT_Cascade/FinalAnalysisCode/`

## Code/File structure:
- README.md : this file
- env.sh : file to define environment variables needed by the scripts
- data : folder containing numpy files with the processed data and simulations
- PDFs : folder containing all PDF in pickle format
- sensitivity : folder containing calculated sensitivities in pickle format
- resources : folder containing limits by other experiments and the PPPC spectra tables
- python : various python code for the analysis
  - fileList.py : provides filenames and paths for the data files
  - IceCube_sim.py : functions related to IceCube simulations
  - PPPC_spectra.py : functions to load and interpolate PPPC spectra 
  - fluxCalculation.py : helper functions for DM and atmospheric flux calculations
  - createPDFs.py : code for the calculation of PDFs with parameters
    - `-t type` (background/annihilation/decay)
    - `--lecut value` cut value on LE BDT score
    - `--hecut value` cut value on HE BDT score
    - `-s systematics`
    - `-o oversampling_factor`
    - `-c channel` annihilation/decay channel (signal only)
    - `-m mass` (signal only)
    - `-p profile` DM profile NFW/Burkert (signal only)
    
  - SensitivityCalucalation.py : code for sensitivity calculation with parameters:
    - `-t type` (annihilation/decay)
    - `--lecut value` cut value on LE BDT score
    - `--hecut value` cut value on HE BDT score
    - `-s systematics`
    - `-o oversampling_factor`
    - `-c channel` annihilation/decay channel
    - `-m mass`
    - `--profile profile` DM profile NFW/Burkert
    - `-l confidence_level`
    - `-d likelihood_method` (poisson or effective)
    - `-e energy_rebinning`
    - `-p psi_rebinning`
  - controlVariables.py : code to produce histograms of E_rec and psi_rec for pre-Unblinding cross-checks for data consistency with parameters:
    - `--lecut value` cut value on LE BDT score
    - `--hecut value` cut value on HE BDT score
    - `--psicut value` cut value on psi (to select off-source region)

- plotting : python plotting scripts
  - style.py : general style definitions
  - plot_DMProfiles.py : figures of DM halo profiles and J factors
  - plot_PDF.py : plot signal and background PDFs with options:
    - `-t type` (background, annihilation, decay)
    - `-s sample` (LE or HE)
    - `-c channel` (for signal only)
    - `-m mass` (for signal only)
    - `-o oversampling` (for signal only)
  - plot_Sensitivities.py : plot sensitivities
    - tbd
  - plot_ControlDistributions.py : plot pre-Unblinding control distribution and perform goddness of fit tests

- submit : sample scripts for condor/dagman submission


## Prerequisites
- icerec build, here from /data/user/sbaur/metaprojects/icerec/build/
- python 2, from /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1
- python modules: iminuit, scipy, astropy, OptionParser

## Workflow / examples

1. PDF creation for signal and background, specific for each channel and mass; example for 1 TeV annihilation to nuE, HE sample:
```
./createPDFs.py -t annihilation -p NFW -c nue -m 1000 --lecut -1.0 --hecut 0.3 -s nominal -o 100
./createPDFs.py -t background --lecut -1.0 -lhecut 0.3 -s nominal
```

2. Calculate sensitivity for each channel and mass; example for 1 TeV annihilation to nuE, HE sample:
```
./SensitivityCalculation.py -t annihilation --lecut -1.0 --hecut 0.3 -s nominal -o 100 -c nue -m 1000 --profile NFW -l 90 -d effective -e 2 -p 5
```

3. Unblinded data: tbd...

4. Plotting:
```
python plot_PDF.py -t background -s HE
python plot_PDF.py -t annihilation -s HE -c nue -m 1000
```
