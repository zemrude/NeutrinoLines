from __future__ import print_function, division

import time
import numpy as np
import scipy.stats
import scipy.special
import scipy.optimize

def inv_BinomialError(error, conf_level):
    r"""Estimates the number of trials for a given precision and a confidence_level
    
    Paramters
    ---------
    error : float
        Desired number of precision in %, error in n/ntot * 100
    conf_level : float
        The confidence level in %, n/ntot * 100
    
    Returns
    -------
    ntot : int
        number of trials
    """
    p = conf_level / 100
    error = error / 100
    
    return int(np.round((1 - p) / error**2))

def BinomialError(ntot, n):
    r"""Computes binomial error (sqrt of variance) for measuring
    n out of ntot trials above a given threshold.

    Parameters
    ----------
    ntot : int, array
        Total number of trials
    n : int, array
        Number of trials exceeding a certain threshold

    Returns
    -------
    en : float, array
        Error on n
    """

    if isinstance(n, np.ndarray):
        return np.sqrt(n * (1 - n.astype(float) / ntot))
    return np.sqrt(n * (1 - float(n) / ntot))

def p2sigma(p):
    r"""Convert 1-sided p-value (probability for observing >= x) to
    significance defined by equivalent integral on a normal distribution

    Parameters
    ----------
    p : float, np.ndarray
        1-sided p-value for observing >= a given measurement

    Returns
    -------
    1-sided significance from a Gaussian curve
    """

    if isinstance(p, np.ndarray):
        s = np.empty(p.size, dtype=float)

        mask = (p >= 1)
        s[mask] = -1000
        s[~mask] = scipy.special.erfcinv(2 * p[~mask]) * np.sqrt(2)
        return s
    # END if (array)

    if p >= 1:
        return -1000
    return scipy.special.erfcinv(2 * p) * np.sqrt(2)


def sigma2p(s):
    r"""Convert 1-sided significance defined by integral of normal
    distribution from s to infinity to 1-sided p-value.

    Parameters
    ----------
    s : float, np.ndarray
        1-sided significance from Gaussian curve

    Returns
    -------
    1-sided p-value for observing >= a given measurement
    """
    return 0.5 * scipy.special.erfc(s / np.sqrt(2))


def ts2p(ts, ndf=1.0, eta=0.5, scale=1.0):
    r"""Convert TS to 1-sided p-value (probability for observing ts >= x)

    Parameters
    ----------
    ts : float
        Test statistic value
    ndf : float
        Effective degrees of freedom, should match number of free parameters
        in the likelihood in the limit of large statistics
    eta : float
        Fraction of TS > 0
    scale : float
        Scaling applied to TS values

    Returns
    -------
    1-sided p-value(s) for observing >= ts
    """

    if any(np.atleast_1d(ts) > 1400):
        raise ValueError("scipy.stats.chi2 cannot handle TS > 1400")

    if isinstance(ts, np.ndarray):
        p = np.ones(ts.size, dtype=float)

        mask = (ts < 0)
        p[mask] = (1 - eta) * \
            scipy.stats.chi2.cdf(-ts[mask] / scale, ndf) + eta
        p[~mask] = eta * scipy.stats.chi2.sf(ts[~mask] / scale, ndf)
        return p
    # END if (array)

    if (ts < 0):
        return (1 - eta) * scipy.stats.chi2.cdf(-ts / scale, ndf) + eta
    return eta * scipy.stats.chi2.sf(ts / scale, ndf)


def p2ts(p, ndf=1.0, eta=0.5, scale=1.0):
    r"""Convert 1-sided p-value (probability for observing >= ts) to
    corresponding ts threshold

    Parameters
    ----------
    p : float
        1-sided p-value
    ndf : float
        Effective degrees of freedom, should match number of free parameters
        in the likelihood in the limit of large statistics
    eta : float
        Fraction of TS > 0
    scale : float
        Scaling applied to TS values

    Returns
    -------
    ts threshold(s) for corresponding p-value(s)
    """
    return scale * scipy.stats.chi2.isf(p / eta, ndf)


def ts2sigma(ts, ndf=1.0, eta=0.5, scale=1.0):
    r"""Convert TS to 1-sided significance from a Gaussian curve.

    Parameters
    ----------
    ts : float
        Test statistic (TS) value
    ndf : float
        Effective degrees of freedom, should match number of free parameters
        in the likelihood in the limit of large statistics
    eta : float
        Fraction of TS > 0
    scale : float
        Scaling applied to TS values

    Returns
    -------
    1-sided significance from a Gaussian curve
    """
    p = ts2p(ts, ndf, eta, scale)
    return p2sigma(p)


def sigma2ts(s, ndf=1.0, eta=0.5, scale=1.0):
    r"""Convert 1-sided significance from a Gaussian curve to
    corresponding TS value.

    Parameters
    ----------
    s : float, np.ndarray
        Significance value(s)
    ndf : float
        Effective degrees of freedom, should match number of free parameters
        in the likelihood in the limit of large statistics
    eta : float
        Fraction of TS > 0
    scale : float
        Scaling applied to TS values

    Returns
    -------
    Test statistic (TS) value corresponding to 1-sided significance s
    """
    p = sigma2p(s)
    return p2ts(p, ndf, eta, scale)



def sensitivity_flux(ts, p, llh, inj, fguess, nstep, nsample=1000,
                     path=None, name="", factor=1, npar_fit=2,
                     par_fit=None, trial_sets=None, **kwargs):
    r"""Returns flux needed to exceed test statistic ts with probability p.

    Parameters
    ----------
    ts : float
        Test statistic (TS) threshold
        Median TS is usually ~0 for sensitivity.
        TS ~25 is around 5 sigma for discovery.
    p : float
        Desired probability for exceeding threshold
        (standard setup is 0.9 for sensitivity, 0.5 for discovery)
    llh : class derived from LLHAnalyzer
        fguess : float
        Estimate on flux needed to reach threshold
    nstep : int
        Number of flux steps for fit
    nsample : int
        Number of signal + background scrambles to use during each step
    path : string
        Path for image output. Default None skips plotting altogether
    name : string
        Name of image output at path
    factor : float
        Flux range of injected signal events [fguess/factor, fguess*factor]
    npar_fit : int
        Number of parameters to use in efficiency curve fit (2 or 3)
    par_fit : int
        Parameters seeds to use in fit. Default None uses hard-coded guesses
        that work in most cases.
    trial_sets : dict
        Dictionary of the form {n_inj : trials}, where `trials' have `n_inj'
        signal injections on average. Default None will perform these trials when
        calculating the sensitivity flux.
        NOTE: using pre-run trials omits the need for `nstep' and `nsample'.
    \*\*kwargs
        Arguments passed to do_trials()
        (ex. {'src_ra' : 3.14, 'src_dec' : 0.00} for point source llh)

    Returns
    -------
    flux : float
        Flux needed to reach threshold
    flux_err : float
        Error on flux derived from fit
    """

   
    results = []
    
    if trial_sets is not None:
        # Use pre-run trials
        for ni in trial_sets.keys():
            trials = trial_sets[ni]

            n = trials['TS'][trials['TS'] > ts].size
            ntot = trials['TS'].size

            results.append([ni, n, ntot])
    else:
        fmin = fguess / factor * (1 - p)
        fmax = fguess * factor / p

        if inj.flux2mu(fmin) < 15:
            fmin = 0

        ni_min = fmin
        ni_max = fmax
        
        print(" Scanning injected fraction of events range [%.2f, %.2f]" % (ni_min, ni_max))

        for step in range(nstep):

            t0 = time.time()

            # number of injected signal events for this step
            ni = ni_min + float(ni_max - ni_min) * step / (nstep - 1)

            # run trials with ns injected signal events
            trials = llh.do_trials(nsample, injector=inj, mean_signal=ni, **kwargs)

            n = trials['TS'][trials['TS'] > ts].size
            ntot = trials['TS'].size

            print(" [%2d] ni %5.2f, n %4d, ntot %4d (%.2f sec)" %
                  (step, ni, n, ntot, time.time() - t0))

            results.append([ni, n, ntot])

    # END for (nstep)

    results = np.transpose(results)

    ni, ni_err, images = fit(results[0],  # injected number of signal events (ni)
                             results[1],  # observed TS > ts
                             results[2],  # total samples taken at each step
                             p, ylabel="fraction TS > %.2f" % ts,
                             npar=npar_fit, par=par_fit,
                             image_base=image_base)

    flux = inj.mu2flux(ni)
    flux_err = flux * ni_err / ni

    print("  ni ------------  %.2f +/- %.2f (p %.2f)" % (ni, ni_err, p))
    print("  flux ----------  %.2e +/- %.2e GeV^-1cm^-2s^-1\n" %
          (flux, flux_err))

    for i in images:
        print("Saved %s" % i)

    return (flux, flux_err)

# END sensitivity_flux()

class EfficiencyCurve:
    r"""Class for efficiency curve based on chi-squared cdf

    Attributes
    ----------
    npar : int
        Number of parameters to fit (2 or 3).
        Setting to 2 fits (ndf, x0).
        Setting to 3 fits (ndf, x0, scale).
    ndf : float
        Effective number of degrees of freedom for chi-squared function
    err_ndf : float
        Error on ndf parameter. Set to 0 until user performs
        fit to determine ndf.
    x0 : float
        X-axis offset
    err_x0 : float
        Error on x0 parameter. Set to 0 until user performs
        fit to determine x0.
    scale : float
        X-axis scale
    err_scale : float
        Error on scale parameter. Set to 0 until user performs
        fit to determine scale.
    random : np.random.RandomState
        Random number generator
    """

    def __init__(self, ndf=10, x0=0, scale=1, seed=1, npar=2):
        r"""Constructor

        Parameters
        ----------
        ndf : float
            Effective degrees of freedom
        x0 : float
            X-axis offset
        scale : float
            X-axis scale
        seed : int
            Seed for random number generator used to create
            random samples from the efficiency curve
        npar : int
            Number of parameters to fit
            (3 == x0, ndf, scale  and  2 == x0, ndf)

        NOTE: Don't read too much into x0, ndf parameters - they're all
              unphysical. We're just using a chi-squared cdf because it
              has the limiting behavior of an exponential approach to 0
              for small x and 1 for large x. A gamma function or
              fdistribution would also work but we chose chi-squared
              because we're used to working with that distribution.
        """

        self.npar = npar  # number of fit parameters

        self.ndf = ndf
        self.x0 = x0
        self.scale = scale

        self.err_ndf = 0.
        self.err_x0 = 0.
        self.err_scale = 0.

        self.random = np.random.RandomState(seed)

    @staticmethod
    def eval_3par(x, ndf, x0, scale):
        r"""Evaluate 3 parameter chi2

        Parameters
        ----------
        x : float, np.ndarray
            X-value(s) to evaluate
        ndf : float
            Effective degrees of freedom
        x0 : float
            X-axis offset
        scale : float
            X-axis scale

        Returns
        -------
        Chi-squared cdf evaluate at x given ndf, x0, scale
        """
        return scipy.stats.chi2.cdf(x, ndf, x0, scale)

    def eval_2par(self, x, ndf, x0):
        r"""Evaluate 2 parameter chi2 with fixed scale

        Parameters
        ----------
        x : float, np.ndarray
            X-value(s) to evaluate
        ndf : float
            Effective degrees of freedom
        x0 : float
            X-axis offset

        Returns
        -------
        Chi-squared cdf evaluate at x given ndf, x0 with fixed scale
        """
        return self.eval_3par(x, ndf, x0, self.scale)

    def eval(self, x):
        r"""Evaluate chi2 using member data for ndf, x0, scale

        Parameters
        ----------
        x : float, np.ndarray
            X-value(s) to evaluate

        Returns
        -------
        Chi-squared cdf evaluate at x given ndf, x0, scale
        """
        return self.eval_3par(x, self.ndf, self.x0, self.scale)

    def inv(self, y):
        r"""Invert chi2 cdf using member data for ndf, x0, scale

        Parameters
        ----------
        y : float, np.ndarray
            Y-value(s) to evaluate

        Returns
        -------
        Inverse of chi-squared cdf
        """
        return scipy.stats.chi2.ppf(y, self.ndf, self.x0, self.scale)

    def fit(self, x, y, ey, par=None, verbose=False):
        r"""Fit x,y data to chi-squared function using chi-squared
        minimization of errors.

        Parameters
        ----------
        x : list, np.ndarray
            X-values of dataset
        y : list, np.ndarray
            Y-values of dataset
        ey : list, np.ndarray
            Error on y values of dataset
        par : list
            Parameter seeds, defaults are hard-coded guesses that work
            for most cases.
        verbose : bool
            Print details when verbose is true

        Returns
        -------
        chi2 : float
            Squared deviation / error from the fit
        dof : int
            Degrees of freedom from the fit
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        ey = np.asarray(ey, dtype=float)

        mask = (ey > 0)

        if par is None:
            # guess ndf seed from median
            ndf = x[np.absolute(y - 0.5).argmin()]
            # force ndf seed >= 10, x0 seed 0
            par = [max([10., ndf]), 0.]

        if self.npar == 3:
            if len(par) < 3:
                # append scale seed 1 when doing 3 par fit
                par.append(1)
            func = self.eval_3par
        else:
            func = self.eval_2par

        popt, pcov = scipy.optimize.curve_fit(
            func, x[mask], y[mask], par, ey[mask])

        self.ndf, self.x0 = popt[:2]
        self.err_ndf, self.err_x0 = np.sqrt(np.diag(pcov))[:2]

        if self.npar == 3:
            self.scale = popt[2]
            self.err_scale = np.sqrt(np.diag(pcov))[2]

        chi2 = ((y[mask] - self.eval(x[mask]))**2 / ey[mask]**2).sum()
        dof = x[mask].size - self.npar

        if dof == 0:
            raise ValueError(("Fit has %d degrees of freedom." % dof) +
                             " Try increasing nstep.")

        if verbose:
            print(("\nBest fit:" +
                   "\n  ndf ----------- %9.2e +/- %9.2e" +
                   "\n  x0 ------------ %9.2e +/- %9.2e") %
                  (self.ndf, self.err_ndf, self.x0, self.err_x0))
            if self.npar == 3:
                print("  scale --------- %9.2e +/- %9.2e" %
                      (self.scale, self.err_scale))
            print("  reduced chi2 --  %.2f (dof %d)" % (chi2 / dof, dof))

        return (chi2, dof)

    def sample(self, x, ntot):
        r"""Sample random data points from efficiency curve.

        Parameters
        ----------
        x : array
            X-axis values of data points
        ntot : array
            Number of trials used at each x point

        Returns
        -------
        y : array
            Measured efficiency for each x
        ey : array
            Error on measured efficiency
        """

        p = self.eval(x)                   # probability at each x
        n = self.random.binomial(ntot, p)  # random count given p and ntot
        en = BinomialError(ntot, n)        # error on random count

        y = n.astype(float) / ntot    # measured efficiency at x
        ey = en.astype(float) / ntot  # error on measured efficieny

        return (y, ey)

# END EfficiencyCurve Class

