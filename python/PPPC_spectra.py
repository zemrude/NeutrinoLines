import numpy as np
from scipy import interpolate
from scipy.integrate import quad

import sys, os

base_path = os.environ['ANALYSIS_BASE_PATH']
spectra_table_path = base_path + 'resources/PPPC4_spectra/'

def dNdlog10x_dict():

    final_state=["numu","positrons","gammas","nutau","antiprotons","antideuterons","nue"]

    dNdlog10x_tab={}
    dNdlog10x_tab["numu"]=np.loadtxt(spectra_table_path+'AtProduction_neutrinos_mu.dat',skiprows=1)
    dNdlog10x_tab["positrons"]=np.loadtxt(spectra_table_path+'AtProduction_positrons.dat',skiprows=1)
    dNdlog10x_tab["gammas"]=np.loadtxt(spectra_table_path+'AtProduction_gammas.dat',skiprows=1)
    dNdlog10x_tab["nutau"]=np.loadtxt(spectra_table_path+'AtProduction_neutrinos_tau.dat',skiprows=1)
    dNdlog10x_tab["antiprotons"]=np.loadtxt(spectra_table_path+'AtProduction_antiprotons.dat',skiprows=1)
    dNdlog10x_tab["antideuterons"]=np.loadtxt(spectra_table_path+'AtProduction_antideuterons.dat',skiprows=1)
    dNdlog10x_tab["nue"]=np.loadtxt(spectra_table_path+'AtProduction_neutrinos_e.dat',skiprows=1)

    # Create dictionaries, easier access to the functions later on.

    dNdlog10x={}
    dNdlog10x["eL"]={};dNdlog10x["eR"]={};dNdlog10x["e"]={}
    dNdlog10x["muL"]={};dNdlog10x["muR"]={};dNdlog10x["mu"]={}
    dNdlog10x["tauL"]={};dNdlog10x["tauR"]={};dNdlog10x["tau"]={}
    dNdlog10x["q"]={};dNdlog10x["c"]={};dNdlog10x["b"]={};dNdlog10x["t"]={}
    dNdlog10x["WL"]={};dNdlog10x["WT"]={};dNdlog10x["W"]={}
    dNdlog10x["ZL"]={};dNdlog10x["ZT"]={};dNdlog10x["Z"]={}
    dNdlog10x["g"]={};dNdlog10x["gamma"]={};dNdlog10x["h"]={}
    dNdlog10x["nue"]={};dNdlog10x["numu"]={};dNdlog10x["nutau"]={}

    key_col={}
    key_col["eL"]=2;key_col["eR"]=3;key_col["e"]=4
    key_col["muL"]=5;key_col["muR"]=6;key_col["mu"]=7
    key_col["tauL"]=8;key_col["tauR"]=9;key_col["tau"]=10
    key_col["q"]=11;key_col["c"]=12;key_col["b"]=13;key_col["t"]=14
    key_col["WL"]=15;key_col["WT"]=16;key_col["W"]=17
    key_col["ZL"]=18;key_col["ZT"]=19;key_col["Z"]=20
    key_col["g"]=21;key_col["gamma"]=22;key_col["h"]=23
    key_col["nue"]=24;key_col["numu"]=25;key_col["nutau"]=26

    #print shape(dNdlog10x_tab["numu"])[0]/179. #62 different masses, 179 values of x

    Z_t=np.zeros((62,179))
    x_mass=np.zeros((62,179))
    y_=np.zeros((62,179))
    for i in range(62):
        x_mass[i,:]=dNdlog10x_tab["numu"][0+i*179,0]
    #print x_mass, shape(x_mass)
    for i in range(62):
        y_[i,:]=dNdlog10x_tab["numu"][0:179,1]
    x_t=y_[0,:]
    y_t=x_mass[:,0]

    # Create interpolating functions

    for key in dNdlog10x:
        #print key, key_col[key]
        #key="eL"
        Z={}
        Z["numu"]=np.zeros((62,179))
        Z["positrons"]=np.zeros((62,179))
        Z["gammas"]=np.zeros((62,179))
        Z["nutau"]=np.zeros((62,179))
        Z["antiprotons"]=np.zeros((62,179))
        Z["antideuterons"]=np.zeros((62,179))
        Z["nue"]=np.zeros((62,179))
    
        for state in final_state:
            for i in range(62):
                Z[state][i,:]=dNdlog10x_tab[state][i*179:(i+1)*179,key_col[key]]
            dNdlog10x[key][state]=interpolate.RectBivariateSpline(y_t, x_t,Z[state],kx=1, ky=1, s=0)
    return dNdlog10x
