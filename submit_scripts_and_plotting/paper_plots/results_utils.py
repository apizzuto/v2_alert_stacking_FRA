import numpy as np
import healpy as hp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import scipy.stats as st
import seaborn as sns
import matplotlib.image as mpimg
from matplotlib.lines import Line2D
from matplotlib.path    import Path
from matplotlib.colors import LogNorm
from glob import glob

import warnings, os
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
warnings.filterwarnings("ignore", category=RuntimeWarning)

from francis.time_integrated_scripts.steady_sensitivity_fits import *

class Alert():
    def __init__(self, ind, delta_t=None):
        self.ind = ind
        self.base_trials_path = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/'
        self.alert_df = pd.read_csv(
            '/home/apizzuto/wg-nu-sources/2021_v2_alert_stacking_FRA/francis/icecube_misc/alert_dataframe.csv')
        self.smear_str = 'smeared/'
        self.delta_t = delta_t
        self.results = None
        self.bg_ts = None
        self.pval = None
        self.get_alert_info()
        self.get_background_trials()
        self.get_pvalue()
        
    def get_alert_info(self):
        run_id = self.alert_df['Run ID'][self.ind]
        ev_id = self.alert_df['Event ID'][self.ind]
        mjd = self.alert_df['MJD'][self.ind]
        iso = Time(mjd, format='mjd').iso
        self.run_id = run_id
        self.ev_id = ev_id
        self.mjd = mjd
        self.iso = iso

    def steady_results_str(self):
        return f'index_{self.ind}_run_{self.run_id}_event_{self.ev_id}_steady_seed_1.pkl'

    def transient_results_str(self):
        start_iso = Time(self.mjd - (self.delta_t/86400./2.), format='mjd').iso
        iso_str = start_iso[:10].replace("-" , "_")
        return f'{iso_str}_RUN_{self.run_id}_EVENT_{self.ev_id}_time_window_{self.delta_t:.2e}_results.pickle'

    def get_results_path(self):
        if self.delta_t is None:
            return self.base_trials_path + 'results/' + self.smear_str + self.steady_results_str()
        else:
            return self.base_trials_path + 'results/' + self.smear_str + self.transient_results_str()

    def get_background_path(self):
        delta_t_str = '_steady' if self.delta_t is None else f'_time_{self.delta_t:.1f}'
        dir_str = self.base_trials_path + 'bg/' + self.smear_str
        return dir_str + f'index_{self.ind}_run_{self.run_id}_event_{self.ev_id}{delta_t_str}.pkl'
    
    def get_background_trials(self):
        bg = np.load(self.get_background_path(), allow_pickle=True)
        if self.delta_t is None:
            self.bg_ts = np.asarray(bg['TS'])
        else:
            self.bg_ts =  np.asarray(bg['ts_prior'])

    def get_unblinded_ts(self):
        if self.delta_t is None:
            if self.results is None:
                self.get_steady_unblinded_res()
            self.ts = self.results['TS']
            self.ns = self.results['nsignal']
            self.gamma = self.results['gamma']
            self.ra = self.results['ra']
            self.dec = self.results['dec']
            return self.results['TS']
        else:
            if self.results is None:
                self.get_transient_unblinded_res()
            self.ts = self.results['ts']
            self.ns = self.results['ns']
            return self.results['ts']
        
    def get_steady_unblinded_res(self):
        res = np.load(self.get_results_path(), allow_pickle=True)
        assert res['TS'][0] == res['TS'][1], "Different seeds give different results"
        res = {key: res[key][0] for key in res.keys()}
        self.results = res
        
    def get_transient_unblinded_res(self):
        res = np.load(self.get_results_path(), allow_pickle=True)
        self.results = res
        
    def get_pvalue(self):
        ts = self.get_unblinded_ts()
        time_key = 'steady' if self.delta_t is None else self.delta_t
        res_ts_key = 'TS' if self.delta_t is None else 'ts'
        self.pval = np.count_nonzero(self.bg_ts >= ts) / float(len(self.bg_ts))
        
def get_unblinded_vals(delta_t = None):
    tss = []
    nss = []
    pvals = []
    for ind in range(275):
        if delta_t == 1000.:
            problem_inds = [198, 95, 92]
        elif delta_t == 172800.:
            problem_inds = [198]
        elif delta_t is None:
            problem_inds = [73,  76, 142, 147, 157, 198, 249]
        else:
            print("why why why we here?")
        if ind in problem_inds:
            ts, ns, p = np.nan, np.nan, np.nan
        else:
            al = Alert(ind, delta_t)
            p = al.pval
            ts = al.ts
            ns = al.ns
        tss.append(ts)
        nss.append(ns)
        pvals.append(p)
    return tss, nss, pvals
        
def get_true_binomial_p(delta_t=None, pvals=None, with_inds=False):
    if delta_t == 1000.:
        problem_inds = [198, 95, 92]
    elif delta_t == 172800.:
        problem_inds = [198]
    elif delta_t is None:
        problem_inds = [73,  76, 142, 147, 157, 198, 249]
    else:
        print("why why why we here?")
    if pvals is None:
        pvals, inds = [], []
        for ind in range(275):
            if ind in problem_inds:
                continue
            al = Alert(ind, delta_t)
            pvals.append(al.pval)
            inds.append(ind)
        pvals = np.array(pvals)
        inds = np.asarray(inds)
    else:
        inds = np.r_[:275]
        msk = np.asarray([True if ind not in problem_inds else False for ind in range(275)])
        pvals = np.asarray(pvals)
        inds = inds[msk]
        pvals = pvals[msk]
    idx = np.argsort(pvals)
    inds = inds[idx]
    pvals = np.sort(pvals)
    obs_p = 1.
    tmps = []
    for i, p in enumerate(pvals):
        tmp = st.binom_test(i+1, len(pvals), p, alternative='greater')
        tmps.append(tmp)
        if tmp < obs_p and tmp != 0.0:
            if tmp == 0.0:
                print("WHY DOES THE BINOMIAL VALUE EQUAL ZERO")
            obs_p = tmp
    if with_inds:
        return obs_p, pvals, tmps, inds
    else:
        return obs_p, pvals, tmps

def name_from_mjd(mjd):
    iso = Time(mjd, format='mjd').iso
    yr, month, day = iso.split('-')
    name = f"{yr[2:4]}{month[:2]}{day[:2]}"
    return name

def scale_2d_gauss(arr, sigma_arr, new_sigma):
    tmp = arr**(sigma_arr**2. / new_sigma**2.)/(np.sqrt(2.*np.pi)*new_sigma)* \
                    np.power(np.sqrt(2.*np.pi)*sigma_arr, (sigma_arr**2. / new_sigma**2.))
    return tmp / np.sum(tmp)

def get_spatial_prior_penalty_map(ind, smear=True, nside=256):
    skymap_fits, header = hp.read_map(skymap_files[ind], h=True, verbose=False,
        dtype=None)
    skymap_llh = skymap_fits.copy()
    skymap_fits = np.exp(-1. * skymap_fits / 2.) #Convert from 2LLH to unnormalized probability
    skymap_fits = np.where(skymap_fits > 1e-30, skymap_fits, 0.0)
    skymap_fits = skymap_fits / np.sum(skymap_fits)
    if smear:
        # 64.2 is the 2deltaLLH threshold used for 90% containment. 
        ninety_msk = skymap_llh < 64.2
        init_nside = hp.get_nside(skymap_llh)
        cdf = np.cumsum(np.sort(skymap_fits[ninety_msk][::-1]))
        pixs_above_ninety = np.count_nonzero(cdf> 0.1)
        original_ninety_area = hp.nside2pixarea(init_nside) * pixs_above_ninety
        new_ninety_area = hp.nside2pixarea(init_nside) * np.count_nonzero(skymap_fits[ninety_msk])
        original_ninety_radius = np.sqrt(original_ninety_area / np.pi)
        new_ninety_radius = np.sqrt(new_ninety_area / np.pi)
        scaled_probs = scale_2d_gauss(skymap_fits, original_ninety_radius, new_ninety_radius)
        skymap_fits = scaled_probs
    if hp.pixelfunc.get_nside(skymap_fits)!=nside:
        skymap_fits = hp.pixelfunc.ud_grade(skymap_fits,nside)
        skymap_fits = skymap_fits/skymap_fits.sum()
    ts_penalty = 2. * np.log(skymap_fits / np.max(skymap_fits))
    return ts_penalty

def plot_color_bar(cmap, labels=[0., 5e3], col_label=r'$-2 \Delta \mathrm{LLH}$', 
                   range=[0,6], offset=-13):
    fig = plt.gcf()
    ax = fig.add_axes([0.95, 0.2, 0.03, 0.6])
    labels = ['0' if lab == 0. else '{:.1f}'.format(lab) for lab in labels]
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, orientation="vertical")
    cb.set_label(col_label, labelpad=offset)
    cb.set_ticks([0., 1.])
    cb.set_ticklabels(labels)
    cb.update_ticks()

def get_alert_title(ind):
    run_id = alerts['Run ID'][ind]
    ev_id = alerts['Event ID'][ind]
    return f"Run {run_id}, Event {ev_id}"

def get_unique_id(ind):
    run_id = alerts['Run ID'][ind]
    ev_id = alerts['Event ID'][ind]
    return f"run{run_id:08d}.evt{ev_id:012d}"
    
def plot_penalty_map(ind, reso=None, cmap=None, cmap_range = None,
                     draw_contour=True, ax=None, no_prior=False,
                     title=None, axis_labels=True):
    """Plot skymap of an alert event
    
    Args:
        ind (int): Alert event index
        LLH (bool, default=False): plot LLH vs. scaled probs.
    """
    ts_map = get_spatial_prior_penalty_map(ind)
    contour_90 = np.asarray(np.load(get_contour_file(ind), allow_pickle=True)[-1][0]).T
    extent_dec = np.degrees(np.max(contour_90[0]) - np.min(contour_90[0])) / 2.
    extent_ra = np.degrees(np.max(contour_90[1]) - np.min(contour_90[1])) / 2.
    if reso is None:
        reso = np.max([1, round(extent_ra), round(extent_dec)])
    ra = np.radians(alerts['RA'][ind])
    dec = np.radians(alerts['Dec'][ind])
    if title is None:
        title = get_alert_title(ind)
    if cmap is None:
        cmap = sns.color_palette("mako_r", as_cmap=True)
    if cmap_range is None:
        max_color = max(ts_map)
        max_color = max_color if max_color != 0. else 1.
        min_color = 0.
    else:
        max_color = cmap_range[1]
        min_color = cmap_range[0]
    ts_map[np.isinf(ts_map)] = min_color
    hp.gnomview(ts_map, rot=(np.degrees(ra), np.degrees(dec), 0),
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
    # hp.graticule(verbose=False)
    # plot_labels(dec, ra, reso, with_axis_labels=axis_labels)
    draw_axes(dec, ra, int(reso), axis_labels=axis_labels)

    if draw_contour:
        draw_contour_from_files(ind)
        
    col_label= r'$2\ln(w)$'
    plot_color_bar(cmap, labels = [min_color, max_color], col_label = col_label,
                  )
    
def get_contour_file(ind):
    base_contour_path = '/data/ana/realtime/alert_catalog_v2/contours/'
    name = get_unique_id(ind)
    f_name = base_contour_path + name
    contour_file_name = f_name + f'.contour.pkl'
    if os.path.isfile(contour_file_name):
        return contour_file_name
    return None

def draw_contour_from_files(ind, color='k', lw=2.):
    contour_file = get_contour_file(ind)
    if contour_file is not None:
        cont = np.load(contour_file, allow_pickle=True)
        cont_ls = ['solid', 'dashed']
        for i, contour_level in enumerate(cont):
            for contour in contour_level:
                cont_ra = np.asarray(contour).T[0]
                cont_dec = np.asarray(contour).T[1]
                hp.projplot(np.pi/2. - cont_dec, cont_ra, linewidth=lw, 
                    color=color, linestyle=cont_ls[i], coord='C')
    
def plot_results_zoom(ind, cmap=None, cmap_range = None,
                      draw_contour=True, ax=None, no_prior=False,
                     with_text_box=False, reso=None, crosshair=False,
                     title=None, axis_labels=True):
    """Plot skymap of an alert event
    
    Args:
        ind (int): Alert event index
        LLH (bool, default=False): plot LLH vs. scaled probs.
    """
    ts_map = get_results_map(ind, no_prior=no_prior)
    contour_90 = np.asarray(np.load(get_contour_file(ind), allow_pickle=True)[-1][0]).T
    extent_dec = np.degrees(np.max(contour_90[0]) - np.min(contour_90[0])) / 2.
    extent_ra = np.degrees(np.max(contour_90[1]) - np.min(contour_90[1])) / 2.
    if reso is None:
        reso = np.max([1, round(extent_ra), round(extent_dec)])
    ra = np.radians(alerts['RA'][ind])
    dec = np.radians(alerts['Dec'][ind])
    if title is None:
        title = get_alert_title(ind)
    if cmap is None:
        cmap = sns.color_palette("mako_r", as_cmap=True)
    if cmap_range is None:
        max_color = max(ts_map)
        max_color = max_color if max_color != 0. else 1.
        min_color = 0.
    else:
        max_color = cmap_range[1]
        min_color = cmap_range[0]
    hp.gnomview(ts_map, rot=(np.degrees(ra), np.degrees(dec), 0),
                    cmap=cmap,
                    max=max_color,
                    min=min_color,
                    reso=reso,
                    title=title,
                    notext=True,
                    cbar=False,
                    )
    draw_axes(dec, ra, reso, axis_labels=axis_labels)
    if draw_contour:
        draw_contour_from_files(ind)
    col_label= r'$\max \big(TS_{\mathrm{PS}} + 2\ln(w), 0\big)$' if not no_prior else 'TS$_{\mathrm{PS}}$'
    plot_color_bar(cmap, labels = [min_color, max_color], col_label = col_label)
    if with_text_box:
        run_id = alerts['Run ID'][ind]
        ev_id = alerts['Event ID'][ind]
        al = Alert(ind, None)
        pval = al.pval
        ts = al.ts
        ra_fit = np.degrees(al.ra)
        dec_fit = np.degrees(al.dec)
        ns = al.ns
        gamma = al.gamma
        fontsize=20
        left_text = '\n'.join((
            r'$p =$%.3f' % (pval),
            r'$\Lambda =%.1f$' % (ts),
            r'$\hat{n}_{s} =$%.1f' % (ns),  
            ))
        props = dict(boxstyle='round', facecolor='white', alpha=0.7)
        plt.text(-28., 1.09, left_text, fontsize=fontsize,
            verticalalignment='top', bbox=props)
        right_text = '\n'.join((
            r'$\hat{\gamma} =$%.2f' % (gamma),
            r'$\alpha, \delta =%.1f, %.1f$' % (ra_fit, dec_fit)
            ))
        plt.text(-14., 1.09, right_text, fontsize=fontsize,
            verticalalignment='top', bbox=props)
    if crosshair:
        al = Alert(ind, None)
        ra_fit = al.ra
        dec_fit = al.dec
        markers = [
            reticle(inner=0.35),
            reticle(which='lt'),
            reticle(which='lt', angle=0),
            reticle(inner=0.35),
            reticle(which='rb'),
            reticle(which='rb', angle=0)
            ]
        for marker in markers:
            hp.projplot(
                np.pi/2. - dec_fit, 
                ra_fit, 
                color='w', 
                markersize=35, 
                markeredgewidth=3.0, 
                markeredgecolor='w',
                marker=marker
                )
        
    
def get_results_map(ind, no_prior=False):
    results_base = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/results/smeared/'
    results_path = results_base + 'index_{}_run_*_steady_ts_map_seed_1.pkl'
    results_file = glob(results_path.format(ind))[0]
    results_map = np.load(results_file, allow_pickle=True)
    ts_map = np.zeros(hp.nside2npix(256))
    key = 'TS_spatial_prior_0' if not no_prior else 'TS'
    ts_map[results_map['npix']] = results_map[key]
    return ts_map

def all_sky_contour_map(colorscale='energy', with_problem_inds=False):
    tmp_df = pd.read_csv(
            '/home/apizzuto/wg-nu-sources/2021_v2_alert_stacking_FRA/francis/icecube_misc/alert_dataframe.csv')
    problem_inds = [73,  76, 142, 147, 157, 198, 249]
    zeroes = np.zeros(hp.nside2npix(256))*np.nan
    hp.mollview(zeroes, min=1, cbar=False, badcolor='white', title='', xsize=4000, rot=180)
    hp.graticule(color='grey', lw=0.7, alpha=0.8)
    palette_map = sns.color_palette("rocket", as_cmap=True)
    for ind in range(275):
        if colorscale == 'energy':
            skymap, header = hp.read_map(skymap_files[ind], h=True, verbose=False,
            dtype=None)
            nside = hp.get_nside(skymap)
            header = {name: val for name, val in header}
            energy  
            energy = header['ENERGY']
            min_cb_val = 1e2
            max_cb_val = 5e3
            if energy > max_cb_val:
                energy = max_cb_val
            elif energy < min_cb_val:
                energy = min_cb_val
            scaled_val = (np.log10(energy) - np.log10(min_cb_val)) / (np.log10(max_cb_val) - np.log10(min_cb_val))
            color = palette_map(scaled_val)
        elif colorscale == 'signalness':
            scaled_val = tmp_df.iloc[ind]['Signalness']
            color = palette_map(scaled_val)
        if ind in problem_inds:
            if with_problem_inds:
                color = 'grey'
            else:
                continue
        draw_contour_from_files(ind, color=color, lw=1.)
    if colorscale == 'energy':
        cbar_label = r'$E_{\nu}$ (TeV)'
        norm = LogNorm(vmin=min_cb_val, vmax=max_cb_val)
    elif colorscale == 'signalness':
        cbar_label = 'Signalness'
        norm = mpl.colors.Normalize(vmin=0., vmax=1.)
    legend_els = [Line2D([0], [0], ls = '-', color='k', label=r'50\% cont.'),
            Line2D([0], [0], ls = '--', color='k', label='90\% cont.')]
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
    fig = plt.gcf()
    leg_ax = fig.add_axes([0.5, 0.0, 0.35, 0.06]) 
    leg_ax.legend(loc=(0.2, -0.18), handles=legend_els, ncol = 2, 
                frameon=False)
    leg_ax.axis('off')
    cbaxes = fig.add_axes([0.05, 0.0, 0.45, 0.06]) 
    cb1 = mpl.colorbar.ColorbarBase(cbaxes, cmap=palette_map,
                                    norm=norm,
                                    orientation='horizontal')
    cb1.set_label(cbar_label)
    cb1.ax.tick_params(direction='out', which='both')
    gplane = SkyCoord(frame='galactic', b = np.zeros(5000)*u.degree, l = np.linspace(0.0, 360., 5000)*u.degree)
    gplane_icrs = gplane.icrs
    gcent = SkyCoord(frame='galactic', b = [0.0]*u.degree, l = [0.0]*u.degree)
    gcent_icrs = gcent.icrs
    hp.projplot(np.pi/2. - gplane_icrs.dec.radian, gplane_icrs.ra.wrap_at('360d').radian, c='k', zorder=2, lw=1.)
    hp.projscatter([np.pi/2. - gcent_icrs.dec.radian], [gcent_icrs.ra.wrap_at('360d').radian],
                  marker='o', c='k', s=12)

def reticle(inner=0.5, outer=1.0, angle=0.0, which='lrtb'):
    """Create a reticle or crosshairs marker.

    Parameters
    ----------
    inner : float
        Distance from the origin to the inside of the crosshairs.
    outer : float
        Distance from the origin to the outside of the crosshairs.
    angle : float
        Rotation in degrees; 0 for a '+' orientation and 45 for 'x'.

    Returns
    -------
    path : `matplotlib.path.Path`
        The new marker path, suitable for passing to Matplotlib functions
        (e.g., `plt.plot(..., marker=reticle())`)

    Examples
    --------
    .. plot::
       :context: reset
       :include-source:
       :align: center

        from matplotlib import pyplot as plt
        from ligo.skymap.plot.marker import reticle

        markers = [reticle(inner=0),
                   reticle(which='lt'),
                   reticle(which='lt', angle=45)]

        fig, ax = plt.subplots(figsize=(6, 2))
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylim(-0.5, 0.5)
        for x, marker in enumerate(markers):
            ax.plot(x, 0, markersize=20, markeredgewidth=2, marker=marker)

    """
    angle = np.deg2rad(angle)
    x = np.cos(angle)
    y = np.sin(angle)
    rotation = [[x, y], [-y, x]]
    vertdict = {'l': [-1, 0], 'r': [1, 0], 'b': [0, -1], 't': [0, 1]}
    verts = [vertdict[direction] for direction in which]
    codes = [Path.MOVETO, Path.LINETO] * len(verts)
    verts = np.dot(verts, rotation)
    verts = np.swapaxes([inner * verts, outer * verts], 0, 1).reshape(-1, 2)
    return Path(verts, codes)

def draw_axes(src_dec, src_ra, reso, axis_labels=True):
    plt.plot(
        4.95/3.*reso*np.radians([-1, 1, 1, -1, -1]), 
        4.95/3.*reso*np.radians([1, 1, -1, -1, 1]), 
        color="k", ls="-", lw=3
        )
    ra_scale_factor = 3 if np.degrees(np.abs(src_dec)) < 30. else 10
    num_ra_lines = ra_scale_factor*reso
    num_dec_lines = 3*reso
    ra_axes = np.linspace(
        np.degrees(src_ra)-360.,
        np.degrees(src_ra)+360.,
        721
        )[360%reso::reso]
    ra_axes = ra_axes[
        (ra_axes > (np.degrees(src_ra) - num_ra_lines)) & \
        (ra_axes < (np.degrees(src_ra) + num_ra_lines))
        ]
    ra_axes = np.where(ra_axes > 360., ra_axes - 360., ra_axes)
    dec_axes = np.linspace(
        np.degrees(src_dec) - 180.,
        np.degrees(src_dec) + 180.,
        361.
        )[180%reso::reso]
    dec_axes = dec_axes[
        (dec_axes > (np.degrees(src_dec) - num_dec_lines)) & \
        (dec_axes < (np.degrees(src_dec) + num_dec_lines))
        ]
    dec_axes = dec_axes[(dec_axes > -90.) & (dec_axes < 90.)]
    for tmp_ra in ra_axes:
        tmp_line_ra = np.radians(np.ones(500)*tmp_ra)
        tmp_line_dec = np.radians(np.linspace(dec_axes[0], dec_axes[-1], 500))
        hp.projplot(np.pi/2. - tmp_line_dec, tmp_line_ra, linewidth=0.5, 
                    color='k', linestyle='dotted', coord='C')
    for tmp_dec in dec_axes:
        tmp_line_dec = np.radians(np.ones(500)*tmp_dec)
        tmp_line_ra = np.radians(np.linspace(ra_axes[0], ra_axes[-1], 500))
        hp.projplot(np.pi/2. - tmp_line_dec, tmp_line_ra, linewidth=0.5, 
                    color='k', linestyle='dotted', coord='C')
    plot_labels(src_dec, src_ra, reso, with_axis_labels=axis_labels)
