import numpy as np
import scipy as sp
#from lmfit                import Model
from scipy.optimize       import curve_fit
from scipy.stats          import chi2
import pickle
from glob import glob
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
import seaborn as sns
from francis.time_integrated_scripts import steady_sensitivity_fits
import healpy as hp
import scipy.stats as st
py_version = int(sys.version[0])
 
palette = ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f']
skymap_files = sorted(glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run*.fits.gz'))
l_ind = skymap_files[0].find("Run")
r_ind = skymap_files[0].find("_nside")

def erfunc(x, a, b):
    return 0.5 + 0.5*sp.special.erf(a*x + b)

def chi2cdf(x,df1,loc,scale):
    func = chi2.cdf(x,df1,loc,scale)
    return func

def fsigmoid(x, a, b):
    return 1.0 / (1.0 + np.exp(-a*(x-b)))

def incomplete_gamma(x, a, scale):
    return sp.special.gammaincc( scale*x, a)

def poissoncdf(x, mu, loc):
    func = sp.stats.poisson.cdf(x, mu, loc)
    return func

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def n_to_flux(N, index, delta_t, smear=True):
    """Convert number of events to a flux
    
    Args:
        N (float): signal strength in number of injected events 
        index (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        float: Flux in default units returned by skylab injector
    """
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/sensitivity/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    if py_version == 3:
        with open(fs[0], 'rb') as f:
            signal_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            signal_trials = pickle.load(f)
    fl_per_one = np.mean(np.array(signal_trials['flux']) / np.array(signal_trials['mean_ninj']))
    return fl_per_one * N

def dec_of_map(index, smear=True):
    """Find the location of the best-fit for an alert

    Args:
        index (int): Alert event index 
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        (float, float): RA, Dec of best-fit location from millipede scan
    """
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/sensitivity/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, 1000.0))
    if py_version == 3:
        with open(fs[0], 'rb') as f:
            signal_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            signal_trials = pickle.load(f)
    ra, dec = np.median(signal_trials['ra']), np.median(signal_trials['dec'])
    return ra, dec

def background_distribution(index, delta_t, smear=True):
    """Obtain background TS distribution for an alert event
    
    Args:
        index (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        array: list of TS values
    """
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    if py_version == 3:
        with open(fs[0], 'rb') as f:
            bg_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            bg_trials = pickle.load(f)
    return bg_trials['ts_prior']

def signal_distribution(index, delta_t, ns, smear=True):
    """Obtain signal TS distribution for an alert event
    
    Args:
        index (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        ns (float): Injected signal strength in number of events
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        dict: Signal trials
    """
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/sensitivity/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    if py_version == 3:
        with open(fs[0], 'rb') as f:
            signal_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            signal_trials = pickle.load(f)
    ret = {}
    msk = np.array(signal_trials['mean_ninj']) == ns
    for k, v in signal_trials.iteritems():
        ret[k] = np.array(v)[msk]
    return ret

def pass_vs_inj(index, delta_t, threshold = 0.5, in_ns = True, with_err = True, trim=-1, smear=True):
    """Calculate the efficiency curve for fraction of TS greater
    than a threshold TS as a function of signal strength
    
    Args:
        index (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        threshold (float, default=0.5): Value of CDF of background to compare
            against (0.5 means compare against median)
        in_ns (bool, default=True): Return in ns or in flux
        with_err (bool, default=True): Include an estimation of the error
            on the passing fraction
        trim (int, default=-1): Trim off the final few points from the curve,
            this sometimes improves the fits
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        (array, array, array): Arrays containing the flux, passing-fractions,
            and errors
    """
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    if py_version == 3:
        with open(fs[0], 'rb') as f:
            bg_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            bg_trials = pickle.load(f)
    bg_trials = bg_trials['ts_prior']
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/sensitivity/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    if py_version == 3:
        with open(fs[0], 'rb') as f:
            signal_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            signal_trials = pickle.load(f)
    bg_thresh = np.percentile(bg_trials, threshold * 100.)
    #print(bg_thresh)
    signal_fluxes, signal_indices = np.unique(signal_trials['mean_ninj'], return_index=True)
    signal_indices = np.append(signal_indices, len(signal_trials['ts_prior']))
    #print signal_indices, signal_fluxes
    if trim != -1 and trim < 0:
        signal_indices = signal_indices[:trim]
        signal_fluxes = signal_fluxes[:trim]
    elif trim > 0:
        signal_indices = signal_indices[:trim + 1]
        signal_fluxes = signal_fluxes[:trim]
    #print signal_trials['ts_prior'], signal_trials['mean_ninj']
    passing = np.array([np.count_nonzero(signal_trials['ts_prior'][li:ri] > bg_thresh) / float(ri - li) for li, ri in zip(signal_indices[:-1], signal_indices[1:])])
    if not with_err:
        return signal_fluxes, passing
    else:
        errs = np.array([np.sqrt(p*(1.-p) / float(ri - li)) for p, li, ri in zip(passing, signal_indices[:-1], signal_indices[1:])])
        ngen = np.array([float(ri - li) for li, ri in zip(signal_indices[:-1], signal_indices[1:])])
        ntrig = passing * ngen
        bound_case_pass = (ntrig + (1./3.)) / (ngen + (2./3.))
        bound_case_sigma = np.sqrt(bound_case_pass*(1. - bound_case_pass) / (ngen + 2))
        errs = np.maximum(errs, bound_case_sigma)
        return signal_fluxes, passing, errs
    
def sensitivity_curve(index, delta_t, threshold = 0.5, in_ns = True, with_err = True, trim=-1, ax = None, 
                    p0 = None, fontsize = 16, conf_lev = 0.9, smear=True, legend=True, text=True):
    """Calculate the sensitivity and plot it for an alert event
    
    Args:
        index (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        threshold (float, default=0.5): Value of CDF of background to compare
            against (0.5 means compare against median)
        in_ns (bool, default=True): Return in ns or in flux
        with_err (bool, default=True): Include an estimation of the error
            on the passing fraction
        trim (int, default=-1): Trim off the final few points from the curve,
            this sometimes improves the fits
        ax (axes, default=None): Use already made axes instance
        p0 (array-like, default=None): initial params for sensitivity fit
        conf_lev (float, default=0.9): Confidence level for the sensitivity
        smear (bool, default=True): Correct for systematics in skymap treatment
    """
    signal_fluxes, passing, errs = pass_vs_inj(index, delta_t, threshold=threshold, in_ns=in_ns, with_err=with_err, trim=trim, smear=smear)
    fits, plist = [], []
    for ffunc in [erfunc, incomplete_gamma, fsigmoid]:
            try:
                fits.append(sensitivity_fit(signal_fluxes, 
                    passing, errs, ffunc, p0=p0, conf_lev=conf_lev))
                plist.append(fits[-1]['pval'])
            except:
                pass
        #print("at least one fit failed")
    #Find best fit of the three, make it look different in plot
    plist = np.array(plist)
    if len(plist) > 0:
        best_fit_ind= np.argmax(plist)
        if fits[best_fit_ind]['chi2'] / fits[best_fit_ind]['dof'] < 5:
            fits[best_fit_ind]['ls'] = '-'
    
    if ax==None:
        fig, ax = plt.subplots()
    
    for fit_dict in fits:
        ax.plot(fit_dict['xfit'], fit_dict['yfit'], 
                 label = r'{}: $\chi^2$ = {:.2f}, d.o.f. = {}'.format(fit_dict['name'], fit_dict['chi2'], fit_dict['dof']),
                ls = fit_dict['ls'])
        if fit_dict['ls'] == '-':
            ax.axhline(conf_lev, color = palette[-1], linewidth = 0.3, linestyle = '-.')
            ax.axvline(fit_dict['sens'], color = palette[-1], linewidth = 0.3, linestyle = '-.')
            if text:
                ax.text(6, 0.5, r'Sens. = {:.2f}'.format(fit_dict['sens']))
    if fits[best_fit_ind]['chi2'] / fits[best_fit_ind]['dof'] > 5:
        inter = np.interp(conf_lev, passing, signal_fluxes)
        ax.axhline(conf_lev, color = palette[-1], linewidth = 0.3, linestyle = '-.')
        ax.axvline(inter, color = palette[-1], linewidth = 0.3, linestyle = '-.')
    ax.errorbar(signal_fluxes, passing, yerr=errs, capsize = 3, linestyle='', marker = 's', markersize = 2)
    if legend:
        ax.legend(loc=4, fontsize = fontsize)
    
def calc_sensitivity(index, delta_t, threshold = 0.5, in_ns = True, with_err = True, trim=-1, 
                    conf_lev = 0.9, p0=None, smear=True):
    """Calculate the sensitivity for an alert event
    
    Args:
        index (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        threshold (float, default=0.5): Value of CDF of background to compare
            against (0.5 means compare against median)
        in_ns (bool, default=True): Return in ns or in flux
        with_err (bool, default=True): Include an estimation of the error
            on the passing fraction
        trim (int, default=-1): Trim off the final few points from the curve,
            this sometimes improves the fits
        p0 (array-like, default=None): initial params for sensitivity fit
        conf_lev (float, default=0.9): Confidence level for the sensitivity
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        dict: Sensitivity dictionary
    """
    signal_fluxes, passing, errs = pass_vs_inj(index, delta_t, threshold=threshold, in_ns=in_ns, with_err=with_err, trim=trim, smear=smear)
    fits, plist = [], []
    for ffunc in [erfunc, incomplete_gamma, fsigmoid]:
            try:
                fits.append(sensitivity_fit(signal_fluxes, 
                    passing, errs, ffunc, p0=p0, conf_lev=conf_lev))
                plist.append(fits[-1]['pval'])
            except:
                pass
    #Find best fit of the three, make it look different in plot
    plist = np.array(plist)
    if len(plist) > 0:
        best_fit_ind= np.argmax(plist)
        if fits[best_fit_ind]['chi2'] / fits[best_fit_ind]['dof'] < 5:
            return fits[best_fit_ind]
    inter = np.interp(conf_lev, passing, signal_fluxes)
    return {'sens': inter, 'name': 'linear_interpolation'}
    
def sensitivity_fit(signal_fluxes, passing, errs, fit_func, p0 = None, conf_lev = 0.9):
    """Fit passing fraction with an analytic function
    
    Args:
        signal_fluxes (array): signal strengths
        passing (array): Passing fractions
        errs (array): Errors on passing fractions
        fit_func (function): Function to fit to data
        p0 (array-like, default=None): initial params for sensitivity fit
        conf_lev (float, default=0.9): Confidence level for the sensitivity

    Returns:
        dict: Sensitivity dictionary
    """
    try:
        name = fit_func.__name__
        name = name.replace("_", " ")
    except:
        name = 'fit'
    signal_scale_fac = np.max(signal_fluxes)
    signal_fls = signal_fluxes / signal_scale_fac
    popt, pcov = curve_fit(fit_func, signal_fls, passing, sigma = errs, p0 = p0, maxfev=4000)
    #print popt
    fit_points = fit_func(signal_fls, *popt)
    chi2 = np.sum((fit_points - passing)**2. / errs**2.)
    dof = len(fit_points) - len(popt)
    xfit = np.linspace(np.min(signal_fls) - 0.5/signal_scale_fac, np.max(signal_fls), 5000)
    yfit = fit_func(xfit, *popt)
    pval = sp.stats.chi2.sf(chi2, dof)
    sens = xfit[find_nearest_idx(yfit, conf_lev)]*signal_scale_fac
    return {'popt': popt, 'pcov': pcov, 'chi2': chi2, 
            'dof': dof, 'xfit': xfit*signal_scale_fac, 'yfit': yfit, 
            'name': name, 'pval':pval, 'ls':'--', 'sens': sens}

def pvals_for_signal(index, delta_t, ns = 1, sigma_units = False, smear=True):
    """Calculate pre trial p-values for a certain injected signal strength
        for a certain alert event
    
    Args:
        index (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        ns (int, default=1): Injected signal strength in number of events
        sigma_units (bool, default=False): Return number of sigma significance
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        array: List of p-values or significances
    """
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    if py_version == 3:
        with open(fs[0], 'rb') as f:
            bg_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            bg_trials = pickle.load(f)
    bg_trials = bg_trials['ts_prior']
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/sensitivity/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    if py_version == 3:
        with open(fs[0], 'rb') as f:
            signal_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            signal_trials = pickle.load(f)
    pvals = [100. - sp.stats.percentileofscore(bg_trials, ts, kind='strict') for ts in signal_trials['ts_prior']]
    pvals = np.array(pvals)*0.01
    pvals = np.where(pvals==0, 1e-6, pvals)
    if not sigma_units:
        return pvals
    else:
        return sp.stats.norm.ppf(1. - (pvals / 2.))

def find_all_sens(delta_t, smear=True, with_disc=True, disc_conf=0.5, 
                    disc_thresh=1.-0.0013, verbose=True):
    """Find the sensitivity for all alerts for a certain analysis
    time-window
    
    Args:
        delta_t (float): Time window (1000. or 864000.)
        smear (bool, default=True): Correct for systematics in skymap treatment
        with_disc (bool, default=True): Also calculate discovery potential
        disc_conf (bool, default=0.5): Confidence level for discovery potential
        disc_thresh (float): p-value for discovery potential (default is 3 sigma)

    Returns:
        array: List of sensitivities
    """
    num_alerts = 276
    sensitivities = np.zeros(num_alerts)
    if with_disc:
        discoveries = np.zeros(num_alerts)
    for ind in range(num_alerts):
        if verbose:
            print(ind, end=' ')
        try:
            sens = n_to_flux(calc_sensitivity(ind, delta_t, smear=smear)['sens'], 
                                ind, delta_t, smear=smear)
            sensitivities[ind] = sens
            if with_disc:
                disc = n_to_flux(calc_sensitivity(ind, delta_t, threshold=disc_thresh, 
                                    conf_lev=disc_conf, smear=smear)['sens'], ind, delta_t, smear=smear)
                discoveries[ind] = disc
            if sens*delta_t*1e6 < 1e-1:
                if verbose:
                    print("Problem calculating sensitivity for alert index {}".format(ind))
        except (IOError, ValueError, IndexError) as err:
            print(err)
    if with_disc:
        return sensitivities, discoveries
    else:
        return sensitivities

def ns_fits_contours(index, delta_t, smear=True, levs = [5., 25., 50., 75., 95.]):
    """Calculate the ns_bias plot contours
    
    Args:
        index (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        smear (bool, default=True): Correct for systematics in skymap treatment
        levs (arr): percentiles used in the ns contours

    Returns:
        (array, array): List of strengths and corresponding bias contours
    """
    smear_str = 'smeared/' if smear else 'norm_prob/'
    levs = np.array(levs)
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/fits/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    if py_version == 3:
        with open(fs[0], 'rb') as f:
            signal_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            signal_trials = pickle.load(f)
    true_inj = np.array(signal_trials['true_ns'])
    ns_fit = np.array(signal_trials['ns_prior'])
    ninjs = np.unique(true_inj)
    if max(ninjs) < 10:
        print('Index {} has max {}'.format(index, max(ninjs)))
    contours = np.stack([np.percentile(ns_fit[true_inj==ninj], levs) for ninj in ninjs])
    return ninjs, contours.T

def ns_fits_contours_plot(index, delta_t, smear=True, levs=[5., 25., 50., 75., 95.],
                      show=False, col='navy green', custom_label = 'Median', ax=None,
                      xlabel=True, ylabel=True, legend=True):
    """Calculate the ns_bias plot contours and plot them
    
    Args:
        index (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        smear (bool, default=True): Correct for systematics in skymap treatment
        levs (arr): percentiles used in the ns contours
    """
    if ax is None:
        fig, ax = plt.subplots()
    ninj, fits = ns_fits_contours(index, delta_t, smear=smear, levs=levs)
    ax.plot(ninj, fits[2], label = custom_label, color = sns.xkcd_rgb[col])
    ax.fill_between(ninj, fits[0], fits[-1], alpha=0.3,
                    label='Central 90\%', color = sns.xkcd_rgb[col], lw=0)
    ax.fill_between(ninj, fits[1], fits[-2], alpha=0.5,
                    label='Central 50\%', color = sns.xkcd_rgb[col], lw=0)
    expectation = ninj
    exp_col = 'dark grey'
    ax.plot(ninj, expectation, ls = '--', color = sns.xkcd_rgb[exp_col])
    if legend:
        ax.legend(loc=4)
    if xlabel:
        ax.set_xlabel(r'$n_{\mathrm{inj}}$')
    ax.set_xlim(0., max(ninj))
    ax.set_ylim(0., 80)
    if ylabel:
        ax.set_ylabel(r'$\hat{n}_{s}$')
    if show:
        plt.show()

def fitting_bias_summary(delta_t, sigs=[2., 5., 10.], smear=True, containment=50.):
    """Calculate the ns_bias plot contours for all alert events
    
    Args:
        delta_t (float): Time window (1000. or 864000.)
        sigs (arr): Injected signal strength values for comparison
        smear (bool, default=True): Correct for systematics in skymap treatment
        containment (float): Central percentage for the fit contours

    Returns:
        (array, array): List of bias and variance of contours
    """
    bias = {sig: [] for sig in sigs}; spread = {sig: [] for sig in sigs};
    levs = [50.-containment / 2., 50., 50.+containment / 2.]
    for ind in range(276):
        try:
            ninjs, contours = ns_fits_contours(ind, delta_t, smear=smear, levs=levs)
        except:
            for sig in sigs:
                bias[sig].append(0.0)
                spread[sig].append(0.0)
            continue
        for sig in sigs:
            try:
                n_ind = np.argwhere(ninjs == sig)[0][0]
                bias[sig].append(contours[1][n_ind])
                spread[sig].append(contours[-1][n_ind] - contours[0][n_ind])
            except:
                bias[sig].append(0.0)
                spread[sig].append(0.0)
    return bias, spread

def background(index, delta_t, smear=True):
    """Obtain the background distributions for an event
    
    Args:
        index (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        smear (bool, default=True): Correct for systematics in skymap treatment
    
    Returns:
        dict: background trials
    """
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    if py_version == 3:
        with open(fs[0], 'rb') as f:
            bg_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            bg_trials = pickle.load(f)
    return bg_trials

def plot_zoom_from_map(skymap, ind, reso=1., cmap=None, draw_contour=True, ax=None, 
                      col_label= r'$\log_{10}$(prob.)'):
    """Plot skymap of an alert event
    
    Args:
        skymap (arr): healpy array
        index (int): Alert event index 
    """
    s, header = hp.read_map(skymap_files[ind], h=True, verbose=False)
    header = {name: val for name, val in header}
    nside = hp.get_nside(s)
    area = np.count_nonzero(s < 64.2) * hp.nside2pixarea(nside) * 180.**2. / (np.pi**2.)
    reso *= int(np.sqrt(area))
    reso = np.max([reso, 1.])
    original_LLH = s
    ra = np.radians(header['RA'])
    dec = np.radians(header['DEC'])
    title = skymap_files[ind][l_ind:r_ind].replace('Run', 'Run ').replace('_', ', Event ')
    if cmap is None:
        pdf_palette = sns.color_palette("Blues", 500)
        cmap = mpl.colors.ListedColormap(pdf_palette)
    if np.count_nonzero(skymap > 0.0) > 1:
        max_color = np.max(skymap)
        min_color = 0.
    else:
        max_color =  -1.8 #5 #max(skymap)
        min_color = -5.  #0.
    
    hp.gnomview(skymap, rot=(np.degrees(ra), np.degrees(dec), 0),
                    cmap=cmap,
                    max=max_color,
                    min=min_color,
                    reso=reso,
                    title=title,
                    notext=True,
                    cbar=False
                    #unit=r""
                    )

    plt.plot(4.95/3.*reso*np.radians([-1, 1, 1, -1, -1]), 4.95/3.*reso*np.radians([1, 1, -1, -1, 1]), color="k", ls="-", lw=3)
    hp.graticule(verbose=False)
    steady_sensitivity_fits.plot_labels(dec, ra, reso)
    con_nside = 256 if area < 5. else 128
    if draw_contour:
        contours = steady_sensitivity_fits.plot_contours(None, original_LLH, levels=[22.2, 64.2], nside = con_nside)
        for contour in np.array(contours).T:
            hp.projplot(contour[0],contour[1],linewidth=1.5,c='k')
    steady_sensitivity_fits.plot_color_bar(cmap = cmap, labels = [min_color, max_color], col_label = col_label)


def background_hotspot_map(ind, delta_t, smear=True):
    """Show where background trials recover hotspots
    
    Args:
        ind (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        smear (bool, default=True): Correct for systematics in skymap treatment
    """
    bg = background(ind, delta_t, smear=smear)
    msk = np.array(bg['ts_prior']) != 0.
    ra, dec = np.array(bg['ra'])[msk], np.array(bg['dec'])[msk]
    theta = np.pi/2. - dec
    inds = hp.ang2pix(256, theta, ra)
    ind_counts = np.unique(inds, return_counts=True)
    reco_hist = np.zeros(hp.nside2npix(256))
    reco_hist[ind_counts[0]] = ind_counts[1]
    plot_zoom_from_map(reco_hist, ind, draw_contour=False, col_label='Counts')

def get_true_fit(ind, delta_t, smear=True):
    """Get unblinded results for an alert followup
    
    Args:
        ind (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        smear (bool, default=True): Correct for systematics in skymap treatment
    """
    name = skymap_files[ind][skymap_files[ind].find('Run')+3:skymap_files[ind].find('_nside')]
    run = name[:name.find('_')]
    event = name[name.find('_') + 1:]
    smear_str = 'smeared/' if smear else 'norm_prob/'
    res_f = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/results/{}*{}*{}_time_window_{:.2e}_results.pickle'.format(smear_str, run, event, delta_t))
    res = np.load(res_f[0])
    return res

##############################################################################
############# REPEAT THIS IN THE STEADY FITS FILE
##############################################################################
def get_true_pval(ind, delta_t, smear=True):
    """Get unblinded p-value for an alert followup
    
    Args:
        ind (int): Alert event index 
        delta_t (float): Time window (1000. or 864000.)
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        float: pre trial p-value for alert followup
    """
    result = get_true_fit(ind, delta_t, smear=smear)
    bg_trials = background(ind, delta_t, smear=smear)
    bg_ts = bg_trials['ts_prior']
    p_val = float(np.count_nonzero(bg_ts >= result['ts'])) / float(len(bg_ts))
    return p_val

def get_true_pval_list(delta_t, smear=True):
    """Get unblinded p-value for all alert followups of a certain time window
    
    Args:
        delta_t (float): Time window (1000. or 864000.)
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        arr: pre trial p-values for alert followup
    """
    problem_inds = [198, 95, 92] if delta_t == 1000. else [198]
    pval_list = []
    for ind in range(len(skymap_files)):
        if ind in problem_inds:
            pval_list.append(1.0)
        else:
            pval = get_true_pval(ind, delta_t, smear=smear)
            pval_list.append(pval)
    return np.array(pval_list)

def get_binomial_p_value_truth(delta_t, smear=True):
    """Calculate the unblinded pre-trial binomial p-value
    
    Args:
        delta_t (float): Time window (1000. or 864000.)
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        float: pre trial binomial p-value for alert followup
    """
    print("CAUTION: ONLY RUN THIS IF YOU HAVE PERMISSION TO LOOK AT REAL DATA")
    obs_p = 1.
    plist = get_true_pval_list(delta_t, smear=smear)
    plist = sorted(plist)
    for i, p in enumerate(plist):
        tmp = st.binom_test(i+1, len(plist), p, alternative='greater')
        if tmp < obs_p and tmp != 0.0:
            if tmp == 0.0:
                print("WHY DOES THE BINOMIAL VALUE EQUAL ZERO")
            obs_p = tmp
    return obs_p
