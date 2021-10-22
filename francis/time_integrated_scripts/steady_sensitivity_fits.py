import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import seaborn as sns
from glob import glob
import healpy as hp
import scipy as sp
#from lmfit                import Model
from scipy.optimize       import curve_fit
from scipy.stats          import chi2
import pickle
import meander
import sys
py_ver = int(sys.version[0])

from francis import utils
f_path = utils.get_francis_path()
utils.initialize_mpl_style()

alerts = pd.read_csv(
    f_path + 'icecube_misc/alert_dataframe.csv'
    )

# Commented path points to where trials were originally stored
# base_path = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/'
base_path = '/data/ana/analyses/NuSources/2021_v2_alert_stacking_FRA/analysis_trials/'
skymap_files = sorted(glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run*.fits.gz'))


palette = ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f']
l_ind = skymap_files[0].find("Run")
r_ind = skymap_files[0].find("_nside")

def erfunc(x, a, b):
    return 0.5 + 0.5*sp.special.erf(a*x + b)

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

def load_bg(ind, smear=True):
    """Obtain background TS distribution for an alert event
    
    Args:
        ind (int): Alert event index 
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        dict: Background trial dictionary
    """
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob(base_path + 'bg/{}index_{}_run_*_steady.pkl'.format(smear_str, ind))
    if len(fs) == 0:
        return None
    if py_ver == 3:
        with open(fs[0], 'rb') as f:
            bg_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            bg_trials = pickle.load(f)
    for k in bg_trials.keys():
        bg_trials[k] = np.array(bg_trials[k])
    return bg_trials 

def load_sig(ind, gamma=2.0, fits=True, smear=True):
    """Obtain signal TS distribution for an alert event
    
    Args:
        ind (int): Alert event index 
        gamma (float, default=2.0): Injected spectral index
        fits (bool, default=True): If false, load poisson fluctuated trials
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        dict: Signal trials
    """
    smear_str = 'smeared/' if smear else 'norm_prob/'
    f_dir = 'fits/{}'.format(smear_str) if fits else 'sensitivity/{}'.format(smear_str)
    fs = glob(base_path + f_dir + 'index_{}_run_*_steady_gamma_{:.1f}.pkl'.format(ind, gamma))
    if len(fs) == 0:
        return None
    if py_ver == 3:
        with open(fs[0], 'rb') as f:
            sig_trials = pickle.load(f, encoding='latin1')
    else:
        with open(fs[0], 'r') as f:
            sig_trials = pickle.load(f)
    for k in sig_trials.keys():
        sig_trials[k] = np.array(sig_trials[k])
    return sig_trials 

def bg_ts_plot(ind, ax=None, smear=True):
    """Plot background TS distribution for an alert event
    
    Args:
        ind (int): Alert event index 
        smear (bool, default=True): Correct for systematics in skymap treatment
    """
    if ax is None:
        fig, ax = plt.subplots()
    bg = load_bg(ind, smear=smear)
    ax.hist(bg['TS'], bins=np.linspace(0., 15., 21), histtype='stepfilled')
    ax.set_yscale('log')
    ax.set_ylim(8e-1, 1e3)
    ax.set_xlabel('TS')
    ax.set_ylabel('Counts')

def bg_ns_gamma_plot(ind, ax=None, smear=True):
    """Plot background distribution of what values of ns and gamma
    get fit for background only trials
    
    Args:
        ind (int): Alert event index 
        smear (bool, default=True): Correct for systematics in skymap treatment
    """
    if ax is None:
        fig, ax = plt.subplots()
    bg = load_bg(ind, smear=smear)
    img = ax.hist2d(bg['nsignal'], bg['gamma'], bins=[np.linspace(0., 40., 16), 
                        np.linspace(1., 5., 16)], norm=LogNorm(vmin=1, vmax=100), cmin=1)
    ax.set_xlabel(r'$\hat{n}_s$')
    ax.set_ylabel(r'$\hat{\gamma}$')

def pass_vs_inj(ind, gamma = 2.0, thresh = 0.5, smear=True, with_err = True):
    """Calculate the efficiency curve for fraction of TS greater
    than a threshold TS as a function of signal strength
    
    Args:
        ind (int): Alert event index 
        gamma (float): injected spectral index
        thresh (float, default=0.5): Value of CDF of background to compare
            against (0.5 means compare against median)
        with_err (bool, default=True): Include an estimation of the error
            on the passing fraction
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        (array, array, array): Arrays containing the flux, passing-fractions,
            and errors
    """
    bg_trials = load_bg(ind, smear=smear)
    signal_trials = load_sig(ind, gamma=gamma, fits=False, smear=smear)
    bg_thresh = np.percentile(bg_trials['TS'], thresh * 100.)
    fluxes = np.unique(signal_trials['flux'])
    flux_inds = [signal_trials['flux'] == flux for flux in fluxes]
    fls = []; passing = []; ntrials = []
    for ii, msk in enumerate(flux_inds):
        if np.count_nonzero(msk) < 15:
            continue
        else:
            fls.append(fluxes[ii])
            ntrials.append(float(len(signal_trials['TS'][msk])))
            passing.append(np.count_nonzero(signal_trials['TS'][msk] > bg_thresh) / ntrials[-1])            
    passing = np.array(passing); ntrials=np.array(ntrials); fls=np.array(fls)
    if not with_err:
        return fls, passing
    else:
        errs = np.array([np.sqrt(p*(1.-p) / num) for p, num in zip(passing, ntrials)])
        ntrig = passing * ntrials
        bound_case_pass = (ntrig + (1./3.)) / (ntrials + (2./3.))
        bound_case_sigma = np.sqrt(bound_case_pass*(1. - bound_case_pass) / (ntrials + 2))
        errs = np.maximum(errs, bound_case_sigma)
        return fls, passing, errs

def steady_sensitivity(ind, gamma=2.0, thresh = 0.5, with_err = True, 
                        smear=True, conf_lev = 0.9, p0=None):
    """Calculate the sensitivity for an alert event
    
    Args:
        ind (int): Alert event index 
        thresh (float, default=0.5): Value of CDF of background to compare
            against (0.5 means compare against median)
        gamma (float): injected spectral index
        with_err (bool, default=True): Include an estimation of the error
            on the passing fraction
        p0 (array-like, default=None): initial params for sensitivity fit
        conf_lev (float, default=0.9): Confidence level for the sensitivity
        smear (bool, default=True): Correct for systematics in skymap treatment

    Returns:
        dict: Sensitivity dictionary
    """
    signal_fluxes, passing, errs = pass_vs_inj(ind, gamma=gamma, thresh=thresh, 
                                                with_err=with_err, smear=smear)
    fits, plist = [], []
    try:
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, fsigmoid, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
    except:
        print("Sigmoid fit failed")
    try:
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, erfunc, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
    except:
        print("Error function fit failed")
    #Find best fit of the three, make it look different in plot
    plist = np.array(plist)
    best_fit_ind= np.argmax(plist)
    return fits[best_fit_ind]

def steady_sensitivity_curve(ind, gamma=2.0, thresh = 0.5, with_err = True, smear=True, ax = None, 
                    p0 = None, fontsize = 16, conf_lev = 0.9):
    """Calculate the sensitivity and plot it for an alert event
    
    Args:
        index (int): Alert event index 
        gamma (float): injected spectral index
        thresh (float, default=0.5): Value of CDF of background to compare
            against (0.5 means compare against median)
        with_err (bool, default=True): Include an estimation of the error
            on the passing fraction
        ax (axes, default=None): Use already made axes instance
        p0 (array-like, default=None): initial params for sensitivity fit
        conf_lev (float, default=0.9): Confidence level for the sensitivity
        smear (bool, default=True): Correct for systematics in skymap treatment
    """
    signal_fluxes, passing, errs = pass_vs_inj(ind, gamma=gamma, thresh=thresh, 
                                                smear=smear, with_err=with_err)
    fits, plist = [], []
    try:
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, fsigmoid, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
    except:
        print("Sigmoid fit failed")
    try:
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, erfunc, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
    except:
        print("Error function fit failed")
    #Find best fit of the three, make it look different in plot
    if len(plist) > 0:
        plist = np.array(plist)
        best_fit_ind= np.argmax(plist)
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
                ax.text(np.percentile(signal_fluxes, 75), 0.6, r'Sens. = {:.1e}'.format(fit_dict['sens']))
    ax.set_xlabel(r'$\phi_{\mathrm{inj}}$')
    ax.set_ylabel(r'Fraction TS $>$ threshold')
    ax.errorbar(signal_fluxes, passing, yerr=errs, capsize = 3, linestyle='', marker = 's', markersize = 2)
    ax.legend(loc=4, fontsize = fontsize)

def fits_contours(ind, gamma=2.0, smear=True, ns=True, levs = [5., 25., 50., 75., 95.]):
    """Calculate the ns_bias plot contours
    
    Args:
        ind (int): Alert event index 
        gamma (float, default=2.0): injected spectral index
        ns (bool, default=True): ns or gamma fits
        smear (bool, default=True): Correct for systematics in skymap treatment
        levs (arr): percentiles used in the ns contours

    Returns:
        (array, array): List of strengths and corresponding bias contours
    """
    sig_trials = load_sig(ind, gamma=gamma, smear=smear)
    key = 'nsignal' if ns else 'gamma'
    msks = [sig_trials['inj_nsignal'] == ni for ni in np.linspace(1., 60., 60)]
    ninjs = []; conf_levs = []
    for ii, msk in enumerate(msks):
        if np.count_nonzero(msk) < 3:
            continue
        else:
            if not ns:
                msk *= sig_trials['nsignal'] != 0.
            ninjs.append(ii)
            conf_levs.append(np.percentile(sig_trials[key][msk], levs))
    return ninjs, np.array(conf_levs).T

def gamma_fits_plots(ind, smear=True, cols=['navy green', 'navy blue', 'orangey red'], ax=None):
    """Calculate the recovered gamma bias plot contours and plot them
    
    Args:
        ind (int): Alert event index 
        smear (bool, default=True): Correct for systematics in skymap treatment
    """
    if ax is None:
        fig, ax = plt.subplots()
    for gamma, col in zip([2.5], ['navy blue']): #zip([2.0, 2.5, 3.0], cols):
        #show = True if gamma == 3.0 else False
        show=False
        fits_contours_plot(ind, gamma=gamma, ns=False, show=show, col=col, smear=smear,
                          custom_label=r'$\gamma_{\mathrm{inj}} = $' + '{:.1f}'.format(gamma), ax=ax)

def fits_contours_plot(ind, gamma=2.0, ns=True, show=False, col='navy green', smear=True,
                      custom_label = 'Median', ax=None):
    """Calculate the ns_bias plot contours and plot them
    
    Args:
        ind (int): Alert event index 
        gamma (float, default=2.0): injected spectral index
        ns (bool, default=True): ns or gamma fits
        smear (bool, default=True): Correct for systematics in skymap treatment
    """
    if ax is None:
        fig, ax = plt.subplots()
    ninj, fits = fits_contours(ind, ns = ns, gamma=gamma, smear=smear)
    ax.plot(ninj, fits[2], label = custom_label, color = sns.xkcd_rgb[col])
    if ns:
        ax.fill_between(ninj, fits[0], fits[-1], alpha=0.3, 
                         label='Central 90\%', color = sns.xkcd_rgb[col], lw=0)
    ax.fill_between(ninj, fits[1], fits[-2], alpha=0.5, 
                     label='Central 50\%' if ns else None,
                     color = sns.xkcd_rgb[col], lw=0)
    expectation = ninj if ns else [gamma]*len(ninj)
    exp_col = 'dark grey' if ns else col
    ax.plot(ninj, expectation, ls = '--', color = sns.xkcd_rgb[exp_col])
    ax.legend(loc=4 if ns else 2)
    ax.set_xlabel(r'$n_{\mathrm{inj}}$')
    ax.set_xlim(0., 60)
    if ns:
        ax.set_ylim(0., 80)
        ax.set_ylabel(r'$\hat{n}_{s}$')
    else:
        ax.set_ylim(1.5, 4.5)
        ax.set_ylabel(r'$\hat{\gamma}$')
    if show:
        plt.show()

def load_map(ind, probs = False):
    """Get the skymap for a certain alert event

    Args:
        ind (int): Alert event index 
        probs (bool, default=False): Return normalized probability (not smeared)
    """
    skymap, header = hp.read_map(skymap_files[ind], 
        h=True, verbose=False, dtype=None)
    header = {name: val for name, val in header}
    if probs:
        skymap = np.exp(-1. * skymap / 2.)
        skymap = np.where(skymap > 1e-70, skymap, 0.0)
        skymap = skymap / np.sum(skymap)
    return skymap, header

def area_of_map(ind):
    """Get the area enclosed by the systematic 90% contour

    Args:
        ind (int): Alert event index 
    
    Returns:
        float: area in square degrees of 90% contour
    """
    skymap, header = load_map(ind, probs=False,)
    nside = hp.get_nside(skymap)
    area = np.count_nonzero(skymap < 64.2) * hp.nside2pixarea(nside) * 180.**2. / (np.pi**2.)
    return area

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
        name = name.replace("incomplete", "inc.")
    except:
        name = 'fit'
    signal_scale_fac = np.max(signal_fluxes)
    signal_fls = signal_fluxes / signal_scale_fac
    popt, pcov = curve_fit(fit_func, signal_fls, passing, sigma = errs, p0 = p0, maxfev=1000)
    fit_points = fit_func(signal_fls, *popt)
    chi2 = np.sum((fit_points - passing)**2. / errs**2.)
    dof = len(fit_points) - len(popt)
    xfit = np.linspace(np.min(signal_fls), np.max(signal_fls), 100)
    yfit = fit_func(xfit, *popt)
    pval = sp.stats.chi2.sf(chi2, dof)
    sens = xfit[find_nearest_idx(yfit, conf_lev)]*signal_scale_fac
    return {'popt': popt, 'pcov': pcov, 'chi2': chi2, 
            'dof': dof, 'xfit': xfit*signal_scale_fac, 'yfit': yfit, 
            'name': name, 'pval':pval, 'ls':'--', 'sens': sens}

def plot_zoom(ind, LLH=False, reso=1., cmap=None, draw_contour=True, ax=None):
    """Plot skymap of an alert event
    
    Args:
        ind (int): Alert event index
        LLH (bool, default=False): plot LLH vs. scaled probs.
    """
    skymap, header = hp.read_map(skymap_files[ind], h=True, verbose=False,
        dtype=None)
    header = {name: val for name, val in header}
    nside = hp.get_nside(skymap)
    area = np.count_nonzero(skymap < 64.2) * hp.nside2pixarea(nside) * 180.**2. / (np.pi**2.)
    reso *= int(np.sqrt(area))
    reso = np.max([reso, 1.])
    original_LLH = None
    if not LLH:
        original_LLH = np.copy(skymap)
        skymap = np.exp(-1. * skymap / 2.)
        skymap = np.where(skymap > 1e-20, skymap, 0.0)
        skymap = skymap / np.sum(skymap) 
    ra = np.radians(header['RA'])
    dec = np.radians(header['DEC'])
    title = skymap_files[ind][l_ind:r_ind].replace('Run', 'Run ').replace('_', ', Event ')
    if cmap is None:
        pdf_palette = sns.color_palette("Blues", 500)
        cmap = mpl.colors.ListedColormap(pdf_palette)
    max_color = 100 if LLH else max(skymap)
    hp.gnomview(skymap, rot=(np.degrees(ra), np.degrees(dec), 0),
                    cmap=cmap,
                    max=max_color,
                    min=0,
                    reso=reso,
                    title=title,
                    notext=True,
                    cbar=False
                    #unit=r""
                    )

    plt.plot(4.95/3.*reso*np.radians([-1, 1, 1, -1, -1]), 4.95/3.*reso*np.radians([1, 1, -1, -1, 1]), color="k", ls="-", lw=3)
    hp.graticule(verbose=False)
    plot_labels(dec, ra, reso)
    original_LLH = skymap if original_LLH is None else original_LLH
    con_nside = 256 if area < 5. else 128
    if draw_contour:
        sample_points = np.array(hp.pix2ang(nside,np.arange(len(original_LLH)))).T
        msk = original_LLH < 110.
        contours = plot_contours(None, original_LLH[msk], 
            sample_points[msk], levels=[22.2, 64.2]) #, nside = con_nside)
        for contour in np.array(contours).T:
            hp.projplot(contour[0],contour[1],linewidth=1.5,c='k')
    col_label=r'$-2\Delta \mathrm{LLH}$' if LLH else "Prob."
    plot_color_bar(cmap = cmap, labels = [0, max_color], col_label = col_label)
    
def plot_labels(src_dec, src_ra, reso, with_axis_labels=True):
    """Add labels to healpy zoom"""
    fontsize = 20
    if np.degrees(src_dec) > 65.:
        ras = [
            np.degrees(src_ra) - reso*3, 
            np.degrees(src_ra), 
            np.degrees(src_ra) + reso*3
            ]
    else:
        ras = [
            np.degrees(src_ra) - reso, 
            np.degrees(src_ra), 
            np.degrees(src_ra) + reso
            ]
    decs = [
        np.degrees(src_dec) - reso, 
        np.degrees(src_dec), 
        np.degrees(src_dec) + reso
        ]
    for ra in ras:
        dec_text = np.pi/2. - src_dec + (np.radians(reso) * 5.5/3.)
        dec_offset = np.abs(np.radians(ra - np.degrees(src_ra)))*np.sin(src_dec)*reso*0.01
        if ra < 360. and ra > 0:
            ra_text = "{:.2f}$^\circ$".format(ra) 
        elif ra > 360.:
            ra_text = "{:.2f}$^\circ$".format(ra - 360.)
        else:
            ra_text = "{:.2f}$^\circ$".format(ra + 360.)
        hp.projtext(dec_text + dec_offset, np.radians(ra),
                    ra_text, lonlat=False,
                    fontsize=20, ha='center')
    for dec in decs:
        scale = np.degrees(hp.rotator.angdist(
            [np.pi/2. - np.radians(dec), src_ra], 
            [np.pi/2. - np.radians(dec), src_ra + np.radians(1.)]
            ))
        ra_text = np.radians(np.degrees(src_ra) + 1.7*reso/scale)
        hp.projtext(np.pi/2. - np.radians(dec), ra_text,
                    "{:.2f}$^\circ$".format(dec), lonlat=False,
                    fontsize=20, ha='right')
    if with_axis_labels:
        plt.text(-1*np.radians(2.4*reso), np.radians(0), r"declination",
                    ha='center', va='center', rotation=90, fontsize=fontsize)
        plt.text(np.radians(0), np.radians(-2*reso), r"right ascension",
                    ha='center', va='center', fontsize=fontsize)
    
    
def plot_contours(proportions, samples, sample_points, levels = None):# , nside = 128):
    r''' Plot containment contour around desired level.
    E.g 90% containment of a PDF on a healpix map

    Parameters:
    -----------
    proportions: list
        list of containment level to make contours for.
        E.g [0.68,0.9]
    samples: array
        array of values read in from healpix map
        E.g samples = hp.read_map(file)
    Returns:
    --------
    theta_list: list
        List of arrays containing theta values for desired contours
    phi_list: list
        List of arrays containing phi values for desired contours
    '''
    if levels is None:
        levels = []
        sorted_samples = list(reversed(list(sorted(samples))))
        nside = hp.pixelfunc.get_nside(samples)
        #sample_points = np.array(hp.pix2ang(nside,np.arange(len(samples)))).T
        for proportion in proportions:
            level_index = (np.cumsum(sorted_samples) > proportion).tolist().index(True)
            level = (sorted_samples[level_index] +
                    (sorted_samples[level_index+1] if level_index+1<len(samples) else 0))/2.0
            levels.append(level)
    # else:
        # samples = hp.pixelfunc.ud_grade(samples, nside)
        # nside = hp.pixelfunc.get_nside(samples)
        # sample_points = np.array(hp.pix2ang(nside,np.arange(len(samples)))).T
    contours_by_level = meander.spherical_contours(sample_points, samples, levels)

    theta_list = []; phi_list=[]
    for contours in contours_by_level:
        for contour in contours:
            theta, phi = contour.T
            phi[phi<0] += 2.0*np.pi
            theta_list.append(theta)
            phi_list.append(phi)

    return theta_list, phi_list

def plot_color_bar(labels=[0., 5e3], col_label=r'$-2 \Delta \mathrm{LLH}$', 
                   range=[0,6], cmap=None, offset=-30):
    fig = plt.gcf()
    ax = fig.add_axes([0.95, 0.2, 0.03, 0.6])
    labels = ['0' if lab == 0. else '{:.1e}'.format(lab) for lab in labels]
    cb = mpl.colorbar.ColorbarBase(ax, orientation="vertical")
    cb.set_label(col_label, labelpad=offset)
    cb.set_ticks([0., 1.])
    cb.set_ticklabels(labels)
    cb.update_ticks()
    
def plot_skymap(ind, LLH=True):
    """Plot skymap of an alert event
    
    Args:
        ind (int): Alert event index
        LLH (bool, default=False): plot LLH vs. scaled probs.
    """
    skymap, header = hp.read_map(skymap_files[ind], h=True, verbose=False,
        dtype=None)
    header = {name: val for name, val in header}
    nside = hp.get_nside(skymap)
    original_LLH = None
    if not LLH:
        original_LLH = np.copy(skymap)
        skymap = np.exp(-1. * skymap / 2.)
        skymap = np.where(skymap > 1e-20, skymap, 0.0)
        skymap = skymap / np.sum(skymap)
    if LLH:
        inf_msk = np.isinf(skymap)
        other_max = np.max(skymap[~inf_msk])
        skymap[inf_msk] = other_max
    pdf_palette = sns.color_palette("Blues", 500)
    cmap = mpl.colors.ListedColormap(pdf_palette)
    hp.mollview(skymap, title="", cbar=True, 
                notext=True, rot=180, hold=False, 
                cmap=cmap, 
                unit = r'$-2 \Delta \mathrm{LLH}$')
    hp.graticule()
    # add some labels
    plt.text(2.0,0., r"$0^\circ$", ha="left", va="center")
    plt.text(1.9,0.45, r"$30^\circ$", ha="left", va="center")
    plt.text(1.4,0.8, r"$60^\circ$", ha="left", va="center")
    plt.text(1.9,-0.45, r"$-30^\circ$", ha="left", va="center")
    plt.text(1.4,-0.8, r"$-60^\circ$", ha="left", va="center")
    plt.text(2.0, -0.15, r"$0\,\mathrm{h}$", ha="center", va="center")
    plt.text(1.333, -0.15, r"$4\,\mathrm{h}$", ha="center", va="center")
    plt.text(.666, -0.15, r"$8\,\mathrm{h}$", ha="center", va="center")
    plt.text(0.0, -0.15, r"$12\,\mathrm{h}$", ha="center", va="center")
    plt.text(-.666, -0.15, r"$16\,\mathrm{h}$", ha="center", va="center")
    plt.text(-1.333, -0.15, r"$20\,\mathrm{h}$", ha="center", va="center")
    plt.text(-2.0, -0.15, r"$0\,\mathrm{h}$", ha="center", va="center")
    original_LLH = skymap if original_LLH is None else original_LLH
    try:
        sample_points = np.array(hp.pix2ang(nside,np.arange(len(original_LLH)))).T
        msk = original_LLH < 110.
        contours = plot_contours(None, original_LLH[msk], sample_points[msk],
            levels=[22.2, 64.2])
        for contour in np.array(contours).T:
            hp.projplot(contour[0],contour[1],linewidth=1.5,c='k')
    except:
        pass
    #plt.draw()

def load_skymap(ind, zoom=True, ax = None):
    if ax is None:
        fig, ax = plt.subplots()
    title = f"Run_{alerts['Run ID'][ind]}_Event_{alerts['Event ID'][ind]}"
    # title = skymap_files[ind][l_ind:r_ind].replace('_', '_Event_').replace('Run', 'Run_')
    zoom_str = 'zoom_LLH' if zoom else 'allsky'
    image = plt.imread('/data/user/apizzuto/fast_response_skylab/alert_event_followup/alert_skymaps/{}_{}.png'.format(title, zoom_str))
    ax.imshow(image); ax.axis('off')


def info_panel_plot(ind, smear=True):
    """Make a whole bunch of plots for one alert event
    
    Args:
        ind (int): Alert event index
        smear (bool, default=True): Correct for systematics in skymap treatment
    """
    fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(28, 8), dpi=250)
    plt.subplots_adjust(wspace=0.2, hspace=0.25)
    load_skymap(ind, ax=axs[0,0])
    bg_ts_plot(ind, smear=smear, ax = axs[0,1])
    bg_ns_gamma_plot(ind, smear=smear, ax = axs[0,2])
    gamma_fits_plots(ind, smear=smear, ax = axs[0,3])
    fits_contours_plot(ind, gamma=2.0, smear=smear, ax = axs[0,4])
    fits_contours_plot(ind, gamma=2.5, smear=smear, col='navy blue', ax = axs[1,0])
    fits_contours_plot(ind, gamma=3.0, smear=smear, col='orangey red', ax = axs[1,1])
    sensitivity_curve(ind, gamma=2.0, smear=smear, ax = axs[1,2])
    sensitivity_curve(ind, gamma=2.5, smear=smear, ax = axs[1,3])
    sensitivity_curve(ind, gamma=3.0, smear=smear, ax = axs[1,4])
    plt.show()


def panel_plot_with_text(ind, smear=True, return_info=False):
    """Make a whole bunch of plots for one alert event, including
    some text information
    
    Args:
        ind (int): Alert event index
        smear (bool, default=True): Correct for systematics in skymap treatment
    """
    fig = plt.figure(constrained_layout=True, figsize=(20,12))
    heights = [1.5, 1, 1]
    gs = gridspec.GridSpec(ncols=3, nrows=3, figure=fig, height_ratios=heights)
    ax_text = fig.add_subplot(gs[0, 2]); ax_text.axis('off')
    skymap, header = load_map(ind)
    txt_size = 30
    ax_text.text(-0.1, 1,"Run, Event: {}, {}".format(header['RUNID'], header['EVENTID']), fontsize=txt_size)
    ax_text.text(-0.1, 0.8, "Time (UTC): {}".format(header['START'][:header['START'].index('.') + 2]), fontsize=txt_size)
    ax_text.text(-0.1, 0.6, "Stream: {}".format(header['I3TYPE']), fontsize=txt_size)
    ax_text.text(-0.1, 0.4, 'Energy, signalness: {} TeV, {}'.format(header['ENERGY'], header['SIGNAL']), fontsize=txt_size)
    ax_text.text(-0.1, 0.2, r'$\delta$: ' + str(header['DEC']) + r'$^{\circ}\; +$' + str(header['DEC_ERR_PLUS']) 
                + r'$^{\circ}\; -$' + str(header['DEC_ERR_MINUS']) + r'$^{\circ}$', fontsize=txt_size)
    ax_text.text(-0.1, 0.0, r'$\alpha$: ' + str(header['RA']) + r'$^{\circ}\; +$' + str(header['RA_ERR_PLUS']) 
                + r'$^{\circ}\; -$' + str(header['RA_ERR_MINUS']) + r'$^{\circ}$', fontsize=txt_size)
    ax1 = fig.add_subplot(gs[0, 0]); load_skymap(ind, ax=ax1)
    ax2 = fig.add_subplot(gs[0, 1]); load_skymap(ind, ax=ax2, zoom=False)
    ax3 = fig.add_subplot(gs[1, 0]); bg_ts_plot(ind, smear=smear, ax=ax3)
    ax4 = fig.add_subplot(gs[1, 1]); bg_ns_gamma_plot(ind, smear=smear, ax=ax4)
    ax5 = fig.add_subplot(gs[1, 2]); gamma_fits_plots(ind, smear=smear, ax = ax5)
    ax7 = fig.add_subplot(gs[2, 0]); fits_contours_plot(ind, smear=smear, gamma=2.5, col='navy blue', ax = ax7)
    ax10 = fig.add_subplot(gs[2, 1]); steady_sensitivity_curve(ind, smear=smear, gamma=2.5, ax = ax10)
    ax13 = fig.add_subplot(gs[2, 2]); steady_sensitivity_curve(ind, smear=smear, gamma=2.5, ax = ax13, conf_lev=0.5, thresh=0.99865)
    if return_info:
        return header
