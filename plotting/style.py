import matplotlib.pyplot as plt
import seaborn as sns

latex = {
    'nuMu':r"$\nu_\mathrm{\mu}$",
    'nuE':r"$\nu_\mathrm{e}$",
    'nuTau':r"$\nu_\mathrm{\tau}$",
    'nuMubar':r"$\bar{\nu}_\mathrm{\mu}$",
    'nuEbar':r"$\bar{\nu}_\mathrm{e}$",
    'nuTaubar':r"$\bar{\nu}_\mathrm{\tau}$",
    'muon':r"$\mu$",
    'W':r"$W^+W^-$",
    'mu':r"$\mu^+\mu^-$",
    'tau':r"$\tau^+\tau^-$",
    'b':r"$b\bar{b}$",
    'numu':r"$\nu_\mathrm{\mu}\bar{\nu}_\mathrm{\mu}$",
    'nue':r"$\nu\bar{\nu}$"
}

def SetFigureStyle():
    plt.style.use('seaborn-whitegrid')
    plt.rcParams.update({'axes.titlesize': 18})
    plt.rcParams.update({'axes.labelsize': 16})
    plt.rcParams.update({'xtick.labelsize' : 14})
    plt.rcParams.update({'ytick.labelsize' : 14})
    plt.rcParams.update({'legend.fontsize' : 16})
    plt.rcParams.update({'figure.dpi': 600.})
    plt.rcParams.update({'lines.linewidth': 3})

def increaseAxisText(x):
    plt.rcParams.update({'axes.titlesize': 18+x})
    plt.rcParams.update({'axes.labelsize': 16+x})
    plt.rcParams.update({'xtick.labelsize' : 14+x})
    plt.rcParams.update({'ytick.labelsize' : 14+x})

def increaseLegendText(x):
    plt.rcParams.update({'legend.fontsize' : 16+x})
