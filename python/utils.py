import numpy as np

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
    
def ConfidenceIntervalError(n):
    r"""
        Calculates the central 68% CI error instead of using sqrt(n)
        It is more appropiate for low n. For n = 0 returns an  90% upper limit
        
        Parameter
        ---------
        n : interger, np.ndarray
            bin content
        
        Returns
        -------
        error : float, np.ndarray
                error defined as (upper - lower)/2 68% for n > 0 or 90% upper_limit for n = 0
        
    """
    from scipy.stats import chi2
    #  1 - a - b = 0.68
    alpha = 0.159
    beta = 0.159
    
    if isinstance(n, np.ndarray):
        mask = (n == 0)
        s = np.empty(n.size, dtype=float)
        s[mask] = 0.5 * chi2(2).ppf(1 - 0.10)
        s[~mask] = (0.5 * chi2(2*(n[~mask]+1)).ppf(1 - beta) - 0.5 * chi2(2*n[~mask]).ppf(alpha))/2.
        return s
    else:
        if n == 0:
            beta = 0.10
            return 0.5 * chi2(2).ppf(1 - beta) # This is the famous 2.3        
        else:
            error_low = 0.5 * chi2(2*n).ppf(alpha)
            error_high = 0.5 * chi2(2*(n+1)).ppf(1 - beta)
            return (error_high - error_low) / 2.

 
def psi_f(RA,decl):
    return np.arccos(np.cos(np.pi/2.-(-29.*np.pi/180))*np.cos(np.pi/2.-decl)\
                      +np.sin(np.pi/2.-(-29.*np.pi/180))*np.sin(np.pi/2.-decl)*\
                       np.cos(RA-266.*np.pi/180))

def merge_bins(histo , time,  bins_merge_E, bins_merge_Psi,nD_ = 2): 
    ''' 
    This function merges the bins in a histogram
    '''
    grid_E   = histo[1] 
    grid_E   = [grid_E[i] for i in range(len(grid_E)) if i%bins_merge_E ==0]
    if nD_ == 2 : 
        grid_psi = histo[2]
        grid_psi = [grid_psi[i] for i in range(len(grid_psi)) if i%bins_merge_Psi ==0]
    lE=len(grid_E);lpsi=len(grid_psi)
    if nD_ == 2:
        tab_events = np.zeros((lE-1,lpsi-1))
        for j in range(bins_merge_E):
                for k in range(lE-1):
                        for j2 in range(bins_merge_Psi):
                                for k2 in range(lpsi-1):
                                        tab_events[k,k2] = tab_events[k,k2] + histo[0][j+k*bins_merge_E,j2+k2*bins_merge_Psi]*time
        return tab_events, grid_E, grid_psi
    else:
        tab_events=np.zeros((lE-1))
        for j in range(bins_merge_E):
                for k in range(lE-1):
                        tab_events[k] = tab_events[k] + histo[0][j+k*bins_merge_Psi]*time
        return tab_events, grid_E