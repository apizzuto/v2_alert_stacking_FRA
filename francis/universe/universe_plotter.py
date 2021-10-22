import numpy as np
from glob import glob
import pandas as pd
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import healpy as hp
import scipy.stats as st
import scipy as sp
import pickle
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import sys
py_ver = int(sys.version[0])
from francis import utils
utils.initialize_mpl_style()
f_path = utils.get_francis_path()


skymap_files = sorted(glob(
    '/data/ana/realtime/alert_catalog_v2/fits_files/Run1*.fits.gz'
    ))
energy_density = {'transient': {'HB2006SFR': 4.8038e51, 
                        'MD2014SFR': 7.099e51, 
                        #6.196e51, commented bc switch to cascades
                        'NoEvolution': 1.8364e52},
                'steady': {'HB2006SFR': 7.6735e43, 
                        'MD2014SFR': 1.134e44, 
                        #9.897e+43, commented bc switch to cascades
                        'NoEvolution': 2.93335e44}}

energy_density_uncertainty = {'transient': {
                            'MD2014SFR': {'plus_one_sigma': 8.92e51,
                                        'plus_two_sigma': 1.11e52,
                                        'minus_one_sigma': 5.50e51,
                                        'minus_two_sigma': 4.15e51},
                            'HB2006SFR': {},
                            'NoEvolution': {}
                            },
                'steady': {'MD2014SFR': {'plus_one_sigma': 1.425e44,
                                        'plus_two_sigma': 1.776e44,
                                        'minus_one_sigma': 8.785e43,
                                        'minus_two_sigma': 6.627e43},
                          'HB2006SFR': {},
                          'NoEvolution': {}
                          }}

firesong_trials_path = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/ts_distributions/'
followup_trials_path = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/'

class UniversePlotter():
    r'''
    Tool to make some helpful plots

    Parameters:
    -----------
        - delta_t (float): Analysis timescale
        - data_years (float): Years of alert events
        - lumi (str): Luminosity function (standard candle of lognormal)
        - evol (str): Evolution model (ie Hopkins and Beacom, Madau and Dickinson)
    '''
    def __init__(self, delta_t, data_years, lumi, evol, **kwargs):
        self.delta_t = delta_t
        self.transient = True if self.delta_t is not None else False
        self.data_years = data_years
        self.lumi = lumi
        self.evol = evol
        self.background_median_ts = None
        self.background_median_p = None
        self.smeared = kwargs.pop('smeared', True)
        if self.lumi == 'LG':
            self.sigma = kwargs.pop('sigma', 0.4)
        else:
            self.sigma = None
        self.ts_fills = [self.background_median_ts, 50.] if self.transient else [self.background_median_ts, 10.]
        self.lumi_str = {'SC': 'Standard Candle', 'LG': 'Log Normal'}[self.lumi]
        self.evol_str = {'HB2006SFR': 'Hopkins and Beacom 2006 SFR',
                            'MD2014SFR': 'Madau and Dickinson 2014 CSFH'}[self.evol]
        self.ts_path = firesong_trials_path
        self.steady_str = '_delta_t_{:.2e}'.format(self.delta_t) if self.transient else '_steady'
        self.time_window_per_year = (365. * 86400.) / (self.delta_t) if self.transient else 1.
        key = 'transient' if self.transient else 'steady'
        self.energy_density = energy_density[key][evol]
        self.no_evol_energy_density = energy_density[key]['NoEvolution']
        self.energy_density_uncertainty = energy_density_uncertainty[key][evol]
        if self.lumi == 'SC':
            self.evol_lumi_str = 'evol_{}_lumi_{}'.format(self.evol, self.lumi)
        else:
            self.evol_lumi_str = 'evol_{}_lumi_{}_sigma_{:.2f}'.format(self.evol, self.lumi, self.sigma)
        self.densities = np.logspace(-11., -6., 21)
        self.luminosities = np.logspace(49, 62, 53) if self.transient else np.logspace(49., 56., 29)
        if self.transient:
            low_energy_msk = self.luminosities * self.delta_t * self.densities[0] < self.energy_density * 100.
            high_energy_msk = self.luminosities * self.delta_t * self.densities[-1] > self.energy_density * 1e-4
            self.luminosities = self.luminosities[low_energy_msk * high_energy_msk]
        self.seconds_per_year = 365.*86400.
        self.med_TS = None
        self.med_p = None
        self.get_labels()
        self.sigs = None
        self.seed = kwargs.pop('seed', 1)
        self.rng = np.random.RandomState(self.seed)
        self.save_figs = kwargs.pop('save', False)
        if self.save_figs:
            self.savepath = kwargs.pop('savepath', './')
        self.verbose = kwargs.pop('verbose', False)

    def two_dim_sensitivity_plot_ts(self, compare=False, log_ts=False, in_ts=True,
                        ts_vs_p=False, discovery=False):
        r'''Two dimensional contour plot to show the sensitivity of the analysis
        in the luminosity-density plane, highlighting the energy requirements
        for the diffuse flux

        Parameters
        ----------
            - compare (bool): include constraints from previous analyses
            - log_ts (bool): contour colors from log or linear ts values
            - in_ts (bool): TS value or binomial-p value
            - ts_vs_p (bool): compare TS and binomial-p value construction
        '''
        if in_ts or ts_vs_p:
            if self.background_median_ts is None:
                self.get_overall_background_ts()
            if self.med_TS is None:
                self.get_med_TS()
        if (not in_ts) or ts_vs_p:
            if self.background_median_p is None:
                self.get_overall_background_p()
            if self.med_p is None:
                self.get_med_p()
        fig, ax = plt.subplots(figsize=(8,5), dpi=200)
        fig.set_facecolor('w')
        X, Y = np.meshgrid(np.log10(self.densities), np.log10(self.plot_lumis))
        if in_ts or ts_vs_p:
            plot_vals = self.med_TS if not log_ts else np.log10(self.med_TS) 
        else:
            plot_vals = self.med_p if not log_ts else np.log10(self.med_p) 
        if in_ts:
            if log_ts:
                levels = np.linspace(np.log10(self.background_median_ts), np.log10(self.background_median_ts) + 2., 11)
            else:
                levels = np.linspace(self.background_median_ts, self.background_median_ts * 10., 11)
        else:
            if log_ts:
                levels = np.linspace(-8., np.log10(self.background_median_p), 11)
            else:
                levels = np.linspace(0.0, self.background_median_p, 11)
        cmap = self.cmap if in_ts else ListedColormap(self.cmap.colors[::-1])
        extend = 'max' if in_ts else 'min'
        cs = ax.contour(X, Y, plot_vals, cmap=cmap, levels=levels, 
                        #vmin=-0.5, 
                        extend=extend)
        csf = ax.contourf(X, Y, plot_vals, cmap=cmap, levels=levels, 
                        #vmin=-0.5, 
                        extend=extend)
        cbar = plt.colorbar(csf) 
        if in_ts or ts_vs_p:
            cbar_lab = r'Median Stacked TS' if not log_ts else r'$\log_{10}($Median Stacked TS$)$'
        else:
            cbar_lab = r'Median binom. p' if not log_ts else r'$\log_{10}($Median binom. p$)$'
        cbar.set_label(cbar_lab, fontsize = 18)
        cbar.ax.tick_params(axis='y', direction='out')
        if in_ts or ts_vs_p:
            if discovery:
                sens_disc = self.med_TS - self.background_three_sigma_ts
            else:
                sens_disc = self.lower_10 - self.background_median_ts
            cs_ts = ax.contour(X, Y, sens_disc, colors=['k'], 
                            levels=[0.0], linewidths=2., zorder=10)
        if (not in_ts) or ts_vs_p:
            linestyle = 'dashed' if ts_vs_p else 'solid'
            if discovery:
                sens_disc = self.background_three_sigma_p - self.med_p
            else:
                sens_disc = self.background_median_p - self.lower_10_p
            cs_ts = ax.contour(X, Y, sens_disc, colors=['k'], 
                            levels=[0.0], linewidths=2., linestyles=linestyle)
        xs = np.logspace(-11., -6., 1000)

        ys_median = self.energy_density / xs / self.seconds_per_year if self.transient else self.energy_density / xs
        plt.plot(np.log10(xs), np.log10(ys_median), color = 'm', lw=1.5)
        for sig in ['one_sigma', 'two_sigma']:
            upper_factor = self.energy_density_uncertainty['plus_'+sig] / self.energy_density
            lower_factor = self.energy_density_uncertainty['minus_'+sig] / self.energy_density
            alpha = 0.45 if sig is 'two_sigma' else 0.75
            plt.fill_between(np.log10(xs), np.log10(ys_median * lower_factor), np.log10(ys_median * upper_factor), 
                    color = 'm', alpha = alpha, lw=0.0, zorder=5)
        plt.text(-9, 53.6, 'Diffuse', color = 'm', rotation=-28, fontsize=18)
        if compare:
            comp_rho, comp_en, comp_str = self.compare_other_analyses()
            plt.plot(comp_rho, comp_en, color = 'gray', lw=2., zorder=5)
        plt.grid(lw=0.0)
        plt.xlim(-11., -6.)
        if self.transient:
            plt.ylim(50., 55.5)
        else:
            plt.ylim(50., 55.)
        if self.transient:
            if self.delta_t == 1e3:
                time_window_str = r'$\pm 500$ s, '
            else:
                time_window_str = r'$\pm 1$ day, '
        else:
            time_window_str = 'Time integrated, '
        custom_labs = [Line2D([0], [0], color = 'k', lw=2., 
            label='This analysis (' + time_window_str + '{:.1f} yr.)'.format(
                self.data_years
                ))]
        if compare:
            custom_labs.append(Line2D([0], [0], color='grey', lw=2., 
                label=comp_str))
        plt.legend(handles=custom_labs, loc=1).set_zorder(20)
        plt.ylabel(self.lumi_label, fontsize = 22)
        plt.xlabel(self.density_label, fontsize = 22)
        title = self.lumi_str + ', ' + self.evol_str
        plt.title(title)
        if self.save_figs:
            for ftype in ['.png', '.pdf']:
                file_path = self.savepath + 'two_dim_sens' + self.steady_str \
                    + '_' + self.evol \
                    + '_{}_years'.format(int(self.data_years)) \
                    + '_' + self.lumi
                plt.savefig(
                    file_path + ftype, 
                    bbox_inches='tight'
                    )
            if self.verbose:
                print("\t - saved plot to {}".format(file_path))

    def rotated_sensitivity_plot_ts(self, log_ts=False, in_ts=True, ts_vs_p=False, compare=False,
                                    discovery=False):
        r'''Two dimensional contour plot to show the sensitivity of the analysis
        in the rotated luminosity-density plane, highlighting the energy requirements
        for the diffuse flux

        Parameters
        ----------
            - compare (bool): include constraints from previous analyses
            - log_ts (bool): contour colors from log or linear ts values
            - in_ts (bool): TS value or binomial-p value
            - ts_vs_p (bool): compare TS and binomial-p value construction
        '''
        if in_ts or ts_vs_p:
            if self.background_median_ts is None:
                self.get_overall_background_ts()
            if self.med_TS is None:
                self.get_med_TS()
        if (not in_ts) or ts_vs_p:
            if self.background_median_p is None:
                self.get_overall_background_p()
            if self.med_p is None:
                self.get_med_p()
        fig, ax = plt.subplots(figsize=(8,5), dpi=200)
        fig.set_facecolor('w')
        X, Y = np.meshgrid(self.densities, self.plot_lumis)
        Y *= X #Scale by the densities
        X = np.log10(X); Y = np.log10(Y)
        if in_ts or ts_vs_p:
            plot_vals = self.med_TS if not log_ts else np.log10(self.med_TS) 
        else:
            plot_vals = self.med_p if not log_ts else np.log10(self.med_p) 
        if in_ts:
            if log_ts:
                levels = np.linspace(np.log10(self.background_median_ts), np.log10(self.background_median_ts) + 2., 11)
            else:
                levels = np.linspace(self.background_median_ts, self.background_median_ts * 10., 11)
        else:
            if log_ts:
                levels = np.linspace(-10., np.log10(self.background_median_p), 11)
            else:
                levels = np.linspace(0.0, self.background_median_p, 11)
        cmap = self.cmap if in_ts else ListedColormap(self.cmap.colors[::-1])
        extend = 'max' if in_ts else 'min'
        cs = ax.contour(X, Y, plot_vals, cmap=cmap, levels=levels, 
                        #vmin=-0.5, 
                        extend=extend)
        csf = ax.contourf(X, Y, plot_vals, cmap=cmap, levels=levels, 
                        #vmin=-0.5, 
                        extend=extend)
        cbar = plt.colorbar(csf) 
        if in_ts or ts_vs_p:
            cbar_lab = r'Median Stacked TS' if not log_ts else r'$\log_{10}($Median Stacked TS$)$'
        else:
            cbar_lab = r'Median binom. p' if not log_ts else r'$\log_{10}($Median binom. p$)$'
        cbar.set_label(cbar_lab, fontsize = 18)
        cbar.ax.tick_params(axis='y', direction='out')
        if in_ts or ts_vs_p:
            if discovery:
                sens_disc = self.med_TS - self.background_three_sigma_ts
            else:
                sens_disc = self.lower_10 - self.background_median_ts
            cs_ts = ax.contour(X, Y, sens_disc, colors=['k'], 
                            levels=[0.0], linewidths=2.)
        if (not in_ts) or ts_vs_p:
            if discovery:
                sens_disc = self.background_three_sigma_p - self.med_p
            else:
                sens_disc = self.background_median_p - self.lower_10_p
            linestyle = 'dashed' if ts_vs_p else 'solid'
            cs_ts = ax.contour(X, Y, sens_disc, colors=['k'], 
                            levels=[0.0], linewidths=2., linestyles=linestyle)
        xs = np.logspace(-11., -6., 1000)
        ys_median = self.energy_density / xs / self.seconds_per_year if self.transient else self.energy_density / xs
        
        plt.plot(np.log10(xs), np.log10(ys_median*xs), color = 'm', lw=1.5)
        for sig in ['one_sigma', 'two_sigma']:
            upper_factor = self.energy_density_uncertainty['plus_'+sig] / self.energy_density
            lower_factor = self.energy_density_uncertainty['minus_'+sig] / self.energy_density
            alpha = 0.45 if sig is 'two_sigma' else 0.75
            plt.fill_between(
                np.log10(xs), np.log10(ys_median * lower_factor * xs), 
                np.log10(ys_median * upper_factor * xs), 
                color = 'm', alpha = alpha, lw=0.0, zorder=4)
        plt.text(-10, np.log10(np.max(ys_median*xs*upper_factor)*1.1), 'Diffuse', 
                        color = 'm', rotation=0, fontsize=18)
        if compare:
            comp_rho, comp_en, comp_str = self.compare_other_analyses()
            plt.plot(comp_rho, comp_rho+comp_en, color = 'gray', lw=2.) #damn look at that log property
        #plt.grid(lw=0.0)
        plt.xlim(-11., -6.)
        plt.ylim(np.log10(np.min(ys_median * lower_factor * xs)*3e-2), 
            np.log10(np.max(ys_median * upper_factor * xs)*2))
        plt.ylabel(self.scaled_lumi_label, fontsize = 22)
        plt.xlabel(self.density_label, fontsize = 22)
        if self.transient:
            if self.delta_t == 1e3:
                time_window_str = '$\pm 500$ s, '
            else:
                time_window_str = '$\pm 1$ day, '
        else:
            time_window_str = 'Time integrated, '
        custom_labs = [Line2D([0], [0], color = 'k', lw=2., label='This analysis (' + time_window_str + '{:.1f} yr.)'.format(self.data_years))]
        if compare:
            custom_labs.append(Line2D([0], [0], color='grey', lw=2., label=comp_str))
        plt.legend(handles=custom_labs, loc=4).set_zorder(20)
        title = self.lumi_str + ', ' + self.evol_str
        plt.title(title)
        if self.save_figs:
            for ftype in ['.png', '.pdf']:
                file_path = self.savepath + 'rotated_two_dim_sens' \
                    + self.steady_str + '_' + self.evol \
                    + '_{}_years'.format(int(self.data_years)) \
                    + '_' + self.lumi
                plt.savefig(
                    file_path + ftype, 
                    bbox_inches='tight'
                )
            if self.verbose:
                print("\t - saved plot to {}".format(file_path))
    
    def get_labels(self):
        r'''Run during initialization to get the correct units 
        for various plots'''
        self.lumi_label = r'$\log_{10}\Big( \frac{\mathcal{E}}{\mathrm{erg}} \Big)$' if self.transient \
                else r'$\log_{10}\Big( \frac{\mathcal{L}}{\mathrm{erg}\;\mathrm{yr}^{-1}} \Big)$'
        self.density_label = r'$\log_{10}\Big( \frac{\dot{\rho}}{ \mathrm{Mpc}^{-3}\,\mathrm{yr}^{-1}} \Big)$' if self.transient \
                else r'$\log_{10}\Big( \frac{\rho}{ \mathrm{Mpc}^{-3}} \Big)$'
        self.cmap = ListedColormap(sns.light_palette((210, 90, 60), input="husl", n_colors=12))
        self.plot_lumis = self.luminosities / self.time_window_per_year
        self.scaled_lumi_label = r'$\log_{10}\Big( \frac{\mathcal{E}\dot{\rho} }{\mathrm{erg}\;\mathrm{Mpc}^{-3}\,\mathrm{yr}^{-1}} \Big)$' if self.transient \
                else r'$\log_{10}\Big( \frac{\mathcal{L} \rho}{\mathrm{Mpc}^{-3} \mathrm{erg}\;\mathrm{yr}^{-1}} \Big)$'
        self.dens_with_units = r'$\rho$ (Mpc$^{-3}$ yr$^{-1}$)' if self.transient else r'$\rho$ (Mpc$^{-3}$)'
        self.dens_units = r'Mpc$^{-3}$ yr$^{-1}$' if self.transient else r'Mpc$^{-3}$'

    def get_med_TS(self):
        r'''Iterate over the parameter space,
        and extract relevant TS values from distributions
        after trials simulating the parameter space have been run'''
        fmt_path = 'ts_dists_{}year_density_{:.2e}_' + self.evol_lumi_str + \
                        '_manual_lumi_{:.1e}' + self.steady_str + '*.npy'
        shape = (self.luminosities.size, self.densities.size)             
        med_TS = np.zeros(shape); lower_10 = np.zeros(shape)
        for ii, lumi in enumerate(self.luminosities):
            for jj, dens in enumerate(self.densities):
                test_en = lumi*dens*self.delta_t if self.transient else lumi*dens
                if test_en > self.energy_density*5.:
                    lower_10[ii, jj] = self.ts_fills[1]
                    med_TS[ii, jj] = self.ts_fills[1]
                elif test_en < self.energy_density*1e-4:
                    lower_10[ii, jj] = self.background_lower_10_ts
                    med_TS[ii, jj] = self.background_median_ts
                else:
                    try:
                        if py_ver == 3:
                            trials_fs = glob(self.ts_path \
                                + fmt_path.format(self.data_years, dens, lumi)
                                )
                            if len(trials_fs) == 0:
                                raise IOError
                            trials = None
                            for f in trials_fs:
                                if trials is None:
                                    trials = np.load(
                                        f, 
                                        allow_pickle=True, 
                                        encoding='latin1'
                                        )
                                else:
                                    tmp = np.load(
                                        f,
                                        allow_pickle=True,
                                        encoding='latin1'
                                        )
                                    trials = np.hstack([trials, tmp])
                        else:
                            trials_fs = glob(self.ts_path \
                                + fmt_path.format(self.data_years, dens, lumi)
                                )
                            if len(trials_fs) == 0:
                                raise IOError
                            trials = None
                            for f in trials_fs:
                                if trials is None:
                                    trials = np.load(f)
                                else:
                                    tmp = np.load(f)
                                    trials = np.hstack([trials, tmp])
                        lower_10[ii, jj] = np.percentile(trials[0], 10.)
                        med_TS[ii, jj] = np.median(trials[0])
                    except IOError as e:
                        lower_10[ii, jj] = np.nan
                        med_TS[ii, jj] = np.nan
        med_TS = np.where(np.isnan(med_TS), self.background_median_ts, med_TS)
        lower_10 = np.where(np.isnan(lower_10), self.background_lower_10_ts, lower_10)
        self.med_TS = med_TS
        self.lower_10 = lower_10

    def get_med_p(self):
        r'''Iterate over the parameter space,
        and extract relevant Binomial-p values values from distributions
        after trials simulating the parameter space have been run'''
        fmt_path = 'ts_dists_{}year_density_{:.2e}_' + self.evol_lumi_str + \
                        '_manual_lumi_{:.1e}' + self.steady_str + '*.npy'
        shape = (self.luminosities.size, self.densities.size)             
        med_p = np.zeros(shape); lower_10_p = np.zeros(shape)
        for ii, lumi in enumerate(self.luminosities):
            for jj, dens in enumerate(self.densities):
                test_en = lumi*dens*self.delta_t if self.transient else lumi*dens
                comp_test_en = self.energy_density*10.2 if self.transient else self.energy_density*100.
                if test_en > comp_test_en:
                    lower_10_p[ii, jj] = 1e-100
                    med_p[ii, jj] = 1e-100
                elif test_en < self.energy_density*3e-3:
                    lower_10_p[ii, jj] = self.background_lower_10_p
                    med_p[ii, jj] = self.background_median_p
                else:
                    try:
                        if py_ver == 3:
                            trials_fs = glob(self.ts_path \
                                + fmt_path.format(self.data_years, dens, lumi)
                                )
                            if len(trials_fs) == 0:
                                raise IOError
                            trials = None
                            for f in trials_fs:
                                if trials is None:
                                    trials = np.load(
                                        f, 
                                        allow_pickle=True, 
                                        encoding='latin1'
                                        )
                                else:
                                    tmp = np.load(
                                        f,
                                        allow_pickle=True,
                                        encoding='latin1'
                                        )
                                    trials = np.hstack([trials, tmp])
                        else:
                            trials_fs = glob(self.ts_path \
                                + fmt_path.format(self.data_years, dens, lumi)
                                )
                            if len(trials_fs) == 0:
                                raise IOError
                            trials = None
                            for f in trials_fs:
                                if trials is None:
                                    trials = np.load(f)
                                else:
                                    tmp = np.load(f)
                                    trials = np.hstack([trials, tmp])
                        lower_10_p[ii, jj] = np.percentile(trials[2], 90.)
                        med_p[ii, jj] = np.median(trials[2])
                    except IOError as e:
                        print(lumi, dens)
                        lower_10_p[ii, jj] = np.nan
                        med_p[ii, jj] = np.nan
        med_p = np.where(np.isnan(med_p), self.background_median_p, med_p)
        lower_10_p = np.where(np.isnan(lower_10_p), self.background_lower_10_p, lower_10_p)
        self.med_p = med_p
        self.lower_10_p = lower_10_p

    def one_dim_ts_distributions(self, only_gold=False, in_ts = True, log_ts=True):
        r'''Assuming that the diffuse flux is saturated,
        show band plot that scans over density
        and shows the TS distributions'''
        ts_or_p = 0 if in_ts else 1
        ts_inds = (0, 2) if not only_gold else (1, 3)
        levels = []; dens = []
        for density in self.densities:
            try:
                if py_ver == 3:
                    ts_fs = glob(
                        self.ts_path + \
                        'ts_dists_{}year_density_{:.2e}_'.format(
                            self.data_years, density) + \
                        self.evol_lumi_str + self.steady_str + '*.npy'
                        )
                    if len(ts_fs) == 0:
                        raise Exception
                    ts = None
                    for f in ts_fs:
                        if ts is None:
                            ts = np.load(
                                f, 
                                allow_pickle=True, 
                                encoding='latin1'
                                )
                        else:
                            tmp = np.load(
                                f,
                                allow_pickle=True,
                                encoding='latin1'
                                )
                            ts = np.hstack([ts, tmp])
                else:
                    ts_fs = glob(
                        self.ts_path + \
                        'ts_dists_{}year_density_{:.2e}_'.format(
                            self.data_years, density) + \
                        self.evol_lumi_str + self.steady_str + '*.npy'
                        )
                    if len(ts_fs) == 0:
                        raise Exception
                    ts = None
                    for f in ts_fs:
                        if ts is None:
                            ts = np.load(f)
                        else:
                            tmp = np.load(f)
                            ts = np.hstack([ts, tmp])
            except Exception as e:
                print(e)
                #continue
            dens.append(density)
            levels.append(np.percentile(ts[ts_inds[ts_or_p]], [5, 25, 50, 75, 95]))
        levels = np.array(levels).T
        fig = plt.figure(dpi=150, figsize=(8,5))
        fig.set_facecolor('w')
        plt.fill_between(dens, levels[0], levels[-1], alpha = 0.5, 
                            color = sns.xkcd_rgb['light navy blue'], 
                            linewidth = 0.0, label = 'Central 90\%')
        plt.fill_between(dens, levels[1], levels[-2], alpha = 0.75, 
                            color = sns.xkcd_rgb['light navy blue'], 
                            linewidth = 0.0, label = 'Central 50\%')
        plt.plot(dens, levels[2], color = sns.xkcd_rgb['light navy blue'])
        plt.title("{}, {}".format(self.lumi_str, self.evol_str))
        plt.xlabel(self.dens_with_units)
        ylab = 'TS' if in_ts else 'Binomial p'
        plt.ylabel(ylab)
        plt.xscale('log')
        loc = 1 if in_ts else 4
        plt.legend(loc=loc, frameon=False)
        if log_ts:
            plt.yscale('log')
        if self.save_figs:
            for ftype in ['.png', '.pdf']:
                file_path = self.savepath + 'one_dim_ts_dist_' \
                    + self.steady_str + '_' + self.evol \
                    + '_{}_years'.format(int(self.data_years)) \
                    + '_' + self.lumi
                plt.savefig(
                    file_path + ftype, 
                    bbox_inches='tight'
                    )
            if self.verbose:
                print("\t - saved plot to {}".format(file_path))

    def ts_and_ps_plot(self, only_gold=False, log_ts=True):
        r'''Make TS distributions for density, luminosity 
        pairs that saturate the diffuse flux'''
        ts_inds = (0, 2) if not only_gold else (1, 3)
        fig, axs = plt.subplots(ncols=2, nrows=1, dpi=200, sharey=True, figsize=(10,4))
        plt.subplots_adjust(wspace=0.08)
        for density in self.densities[::4]:
            try:
                if py_ver == 3:
                    ts_fs = glob(
                        self.ts_path + \
                        'ts_dists_{}year_density_{:.2e}_'.format(
                            self.data_years, density) + \
                        self.evol_lumi_str + self.steady_str + '*.npy'
                        )
                    if len(ts_fs) == 0:
                        raise Exception
                    ts = None
                    for f in ts_fs:
                        if ts is None:
                            ts = np.load(
                                f, 
                                allow_pickle=True, 
                                encoding='latin1'
                                )
                        else:
                            tmp = np.load(
                                f,
                                allow_pickle=True,
                                encoding='latin1'
                                )
                            ts = np.hstack([ts, tmp])
                else:
                    ts_fs = glob(
                        self.ts_path + \
                        'ts_dists_{}year_density_{:.2e}_'.format(
                            self.data_years, density) + \
                        self.evol_lumi_str + self.steady_str + '*.npy'
                        )
                    if len(ts_fs) == 0:
                        raise Exception
                    ts = None
                    for f in ts_fs:
                        if ts is None:
                            ts = np.load(f)
                        else:
                            tmp = np.load(f)
                            ts = np.hstack([ts, tmp])
            except IOError as e:
                continue
            ts_bins = np.logspace(-1., 2., 31) if log_ts else np.linspace(0., 15., 31)
            axs[0].hist(ts[ts_inds[0]], bins = ts_bins, label = r'$\rho = $' 
                    + '{:.1e}'.format(density) + ' ' + self.dens_units, 
                    histtype='step', lw=2.5, weights = [1./len(ts[0])]*len(ts[0]))
            axs[1].hist(ts[ts_inds[1]], bins = np.logspace(-20., 0, 31), label = r'$\rho = $' 
                    + '{:.1e}'.format(density) + ' ' + self.dens_units, 
                    histtype='step', lw=2.5, weights = [1./len(ts[2])]*len(ts[2]))
        fig.suptitle(self.lumi_str + '\n'+ self.evol_str, y=1.03)
        axs[0].set_ylabel('Probability')
        axs[0].set_xlabel('TS')
        axs[1].set_xlabel(r'Binomial $p$')
        if log_ts:
            axs[0].set_xscale('log')
        axs[1].set_xscale('log')
        axs[0].set_yscale('log'); axs[1].set_yscale('log')
        axs[0].set_ylim(8e-3, 7e-1); axs[1].set_ylim(8e-3, 7e-1)
        axs[1].legend(loc=(1.01, 0.1))
        #plt.show()

    def get_overall_background_ts(self, n_trials=5000):
        r'''Sample alert event background distributions
        to get the overall stacked background ts distribution'''
        if self.background_median_ts is not None:
            return self.background_median_ts
        if self.sigs is None:
            self.sigs = self.load_signalness_array()
        bg_trials = followup_trials_path + 'bg/'
        TSs = []
        for ind in range(len(skymap_files)):
            if self.transient and self.delta_t == 1000.:
                problem_inds = [198, 95, 92]
            elif self.transient:
                problem_inds = [198]
            else:
                problem_inds = [73,  76, 142, 147, 157, 198, 249]
            if ind in problem_inds:
                continue
            else:
                smeared_str = 'smeared/' if self.smeared else 'norm_prob/'
                if self.transient:
                    trials_file = glob(bg_trials + smeared_str 
                                + 'index_{}_*_time_{:.1f}.pkl'.format(ind, self.delta_t))[0]
                    if py_ver == 3:
                        trials = np.load(trials_file, encoding='latin1', 
                            allow_pickle=True)
                    else:
                        trials = np.load(trials_file)
                    ts = self.rng.choice(trials['ts_prior'], size=n_trials)
                else:
                    try:
                        trials_files = glob(bg_trials + smeared_str 
                                    + 'index_{}_steady.pkl'.format(ind))
                        if py_ver == 3:
                            trials = np.load(trials_files[0], encoding='latin1',
                                allow_pickle=True)['TS']
                        else:
                            trials = np.load(trials_files[0])['TS']
                        ts = self.rng.choice(np.array(trials), size=n_trials)
                    except:
                        print("NO STEADY TRIALS FOR INDEX {}".format(ind))
                TSs.append(ts)
        TSs = np.array(TSs)
        stacked_ts = np.multiply(TSs, self.sigs[:, np.newaxis])
        stacked_ts = np.sum(stacked_ts, axis=0) / (stacked_ts.shape[0] - 1.) #-1 because we skip one of the maps
        self.background_median_ts = np.median(stacked_ts)
        self.background_lower_10_ts = np.percentile(stacked_ts, 10.)
        self.background_three_sigma_ts = np.percentile(stacked_ts, 99.87)
        self.stacked_ts = stacked_ts
        return self.background_median_ts

    def get_overall_background_p(self, with_tmps=False):
        r'''Sample alert event background distributions
        to get the overall stacked background binomial-p value distribution'''
        bg_trials = followup_trials_path + 'bg/'
        pvals = []
        n_trials = None
        for ind in range(len(skymap_files)):
            if self.transient and self.delta_t == 1000.:
                problem_inds = [198, 95, 92]
            elif self.transient:
                problem_inds = [198]
            else:
                problem_inds = [73,  76, 142, 147, 157, 198, 249]
            if ind in problem_inds:
                continue
            else:
                smeared_str = 'smeared/' if self.smeared else 'norm_prob/'
                if self.transient:
                    trials_file = glob(bg_trials + smeared_str 
                                + 'index_{}_*_time_{:.1f}.pkl'.format(ind, self.delta_t))[0]
                    if py_ver == 3:
                        trials = np.load(trials_file, 
                            encoding='latin1', allow_pickle=True)
                    else:
                        trials = np.load(trials_file)
                    sorted_trial_inds = np.argsort(np.array(trials['ts_prior']))
                    cdf = np.cumsum(np.sort(np.array(trials['ts_prior'])))
                else:
                    trials_files = glob(bg_trials + smeared_str 
                                + 'index_{}_*_steady.pkl'.format(ind))
                    if py_ver == 3:
                        trials = np.load(trials_files[0], encoding='latin1',
                            allow_pickle=True)['TS']
                    else:
                        trials = np.load(trials_files[0])['TS']
                    sorted_trial_inds = np.argsort(np.array(trials))
                    cdf = np.cumsum(np.sort(np.array(trials)))
                inds = np.linspace(0., 1., len(cdf))
                inds = np.where(inds==0., self.rng.uniform(
                    low=0.0, high=np.min(inds[inds != 0.]), size=1)[0], inds)[::-1]
                ps = np.where(cdf != 0.0, inds, 1.0)
                reverse_arg_sort = np.argsort(sorted_trial_inds)
                ps = ps[reverse_arg_sort]
                pvals.append(ps)
        pvals = np.array(pvals)
        background_binomial = []; counter = 0;
        if with_tmps:
            tmps = []
        for realization in pvals.T:
            if with_tmps:
                tmps.append([])
            counter += 1
            realization = np.sort(realization)
            obs_p = 1.
            for i, p in enumerate(realization):
                tmp = st.binom_test(i+1, len(realization), p, alternative='greater')
                if tmp < obs_p and tmp != 0.0:
                    if tmp == 0.0:
                        print("WHY DOES THE BINOMIAL VALUE EQUAL ZERO")
                    obs_p = tmp
                if with_tmps:
                    tmps[-1].append(tmp)
            background_binomial.append(obs_p)
        background_binomial = np.array(background_binomial)
        binomial_median = np.percentile(background_binomial, 50.)
        binomial_lower_10 = np.percentile(background_binomial, 90.)
        self.background_median_p = np.percentile(background_binomial, 50.)
        self.background_lower_10_p = np.percentile(background_binomial, 90.)
        self.background_three_sigma_p = np.percentile(background_binomial, 0.13)
        self.stacked_p = background_binomial
        if with_tmps:
            return self.background_median_p, tmps
        return self.background_median_p

    def inject_and_fit_dens_lumi_plot(self, dens, lumi, in_ts=True, upper_limit=False):
        r'''Assume a certain density and luminosity, 
        inject it, and see what confidence intervals we 
        can construct
        
        Parameters:
        -----------
            - dens (float): Density of sources
            - lumi (float): Luminosity of sources
        '''
        fmt_path = 'ts_dists_{}year_density_{:.2e}_' + self.evol_lumi_str + \
                        '_manual_lumi_{:.1e}' + self.steady_str + '*.npy'
        if py_ver == 3:
            trials_fs = glob(self.ts_path \
                + fmt_path.format(self.data_years, dens, lumi)
                )
            trials = None
            for f in trials_fs:
                if trials is None:
                    trials = np.load(
                        f, 
                        allow_pickle=True, 
                        encoding='latin1'
                        )
                else:
                    tmp = np.load(
                        f,
                        allow_pickle=True,
                        encoding='latin1'
                        )
                    trials = np.hstack([trials, tmp])
        else:
            trials_fs = glob(self.ts_path \
                + fmt_path.format(self.data_years, dens, lumi)
                )
            if len(trials_fs) == 0:
                raise IOError
            trials = None
            for f in trials_fs:
                if trials is None:
                    trials = np.load(f)
                else:
                    tmp = np.load(f)
                    trials = np.hstack([trials, tmp])
        trials = trials[0] if in_ts else trials[2]
        unblinded_val = self.rng.choice(trials)
        self.inject_and_fit_TS_plot(unblinded_val, in_ts=in_ts, show=False, title=False, upper_limit=upper_limit)
        plt.scatter(np.log10(dens), np.log10(lumi / self.time_window_per_year), 
                        marker='*', color = 'k', s=100)
   
    def inject_and_fit_TS_plot(self, unblinded_val, in_ts=True, show=True, 
                    title=True, upper_limit=False):
        r'''Assume a certain unblinded TS value 
        or binomial p-value and see what confidence intervals we 
        can construct
        
        Parameters:
        -----------
            - unblinded_val (float): Unblinded TS or binomial p-value
            - in_ts (float): Use stacked TS construction instead of binomial p-value
        '''
        fig, ax = plt.subplots(dpi=200)
        X, Y = np.meshgrid(np.log10(self.densities), np.log10(self.plot_lumis))
        cont = self.TS_constraints(unblinded_val, in_ts=in_ts, upper_limit=upper_limit)
        levels = [90., 100.] if upper_limit else [0., 50., 90.]
        csf = ax.contourf(X, Y, cont, cmap=self.cmap, levels = levels)
        xs = np.logspace(-11., -6., 1000)
        ys_median = self.energy_density / xs / self.seconds_per_year if self.transient else self.energy_density / xs
        plt.plot(np.log10(xs), np.log10(ys_median), color = sns.xkcd_rgb['dodger blue'], lw=1.5)
        for sig in ['one_sigma', 'two_sigma']:
            upper_factor = self.energy_density_uncertainty['plus_'+sig] / self.energy_density
            lower_factor = self.energy_density_uncertainty['minus_'+sig] / self.energy_density
            alpha = 0.45 if sig is 'two_sigma' else 0.75
            plt.fill_between(np.log10(xs), np.log10(ys_median * lower_factor), np.log10(ys_median * upper_factor), 
                    color = sns.xkcd_rgb['dodger blue'], alpha = alpha, lw=0.0, zorder=10)
        if not upper_limit:
            legend_elements = [Patch(facecolor=csf.cmap.colors[0], label='50\%'),
                        Patch(facecolor=csf.cmap.colors[-2], label='90\%')]
            ax.legend(handles=legend_elements, loc=3)
        ax.set_ylim(np.log10(ys_median.min()*0.5), np.log10(ys_median.max()*2.))
        ax.grid(lw=0.0)
        ax.set_ylabel(self.lumi_label, fontsize = 22)
        ax.set_xlabel(self.density_label, fontsize = 22)
        if in_ts:
            title_str = 'Observed TS={:.1e}'.format(unblinded_val)
        else:
            title_str = 'Observed binom. p={:.1e}'.format(unblinded_val)
        if title:
            plt.title(title_str)
        if self.save_figs:
            for ftype in ['.png', '.pdf']:
                file_path = self.savepath + 'inject_and_fit' \
                    + self.steady_str + '_' + self.evol \
                    + '_{}_years'.format(int(self.data_years)) \
                    + '_' + self.lumi
                plt.savefig(
                    file_path + ftype, 
                    bbox_inches='tight'
                    )
            if self.verbose:
                print("\t - saved plot to {}".format(file_path))
        if show:
            plt.show()

    def TS_constraints(self, obs_val, in_ts = True, upper_limit=False):
        r'''Based on the observed value, get the 
        frequentist confidence intervals

        Parameters:
        -----------
            - obs_val (float): Observed TS of binomial p-value
            - in_ts (bool): Stacked TS or binomial p-value construction
            - upper_limit (bool): If true, return as value compatible with upper limit
        '''
        fmt_path = 'ts_dists_{}year_density_{:.2e}_' + self.evol_lumi_str + \
                        '_manual_lumi_{:.1e}' + self.steady_str + '*.npy'
        ts_p_ind = 0 if in_ts else 2
        containment = np.zeros((self.luminosities.size, self.densities.size)); 
        for ii, lumi in enumerate(self.luminosities):
            for jj, dens in enumerate(self.densities):
                test_en = lumi*dens*self.delta_t if self.transient else lumi*dens
                if test_en > self.energy_density*5.:
                    containment[ii, jj] = 100.
                elif test_en < self.energy_density*1e-4:
                    containment[ii, jj] = 0.
                else:
                    try:
                        if py_ver == 3:
                            trials_fs = glob(self.ts_path \
                                + fmt_path.format(self.data_years, dens, lumi)
                                )
                            trials = None
                            for f in trials_fs:
                                if trials is None:
                                    trials = np.load(
                                        f, 
                                        allow_pickle=True, 
                                        encoding='latin1'
                                        )
                                else:
                                    tmp = np.load(
                                        f,
                                        allow_pickle=True,
                                        encoding='latin1'
                                        )
                                    trials = np.hstack([trials, tmp])
                        else:
                            trials_fs = glob(self.ts_path \
                                + fmt_path.format(self.data_years, dens, lumi)
                                )
                            if len(trials_fs) == 0:
                                raise IOError
                            trials = None
                            for f in trials_fs:
                                if trials is None:
                                    trials = np.load(f)
                                else:
                                    tmp = np.load(f)
                                    trials = np.hstack([trials, tmp])
                        containment[ii, jj] = sp.stats.percentileofscore(trials[ts_p_ind], obs_val)
                    except IOError as e:
                        containment[ii, jj] = 0.
        if upper_limit and in_ts:
            containment = (50. - containment)*2.
        elif upper_limit and not in_ts:
            containment = (50. + containment)*2.
        else:
            containment = np.abs(50. - containment)*2.
        return containment

    def compare_other_analyses(self):
        r'''Get the sensitivities / upper limits
        from previous IceCube analyses'''
        if self.transient:
            nora_comparison = {}
            for key in ['GRB_lims', 'GRB_diffuse', 'CCSN_lims', 'CCSN_diffuse']:
                nora_tmp = np.genfromtxt(f_path + 'icecube_misc/effective_areas_alerts/Nora_{}.csv'.format(key),
                                        delimiter=', ')
                nora_comparison[key] = nora_tmp
            for key in ['GRB_lims']:
                tmp = np.log10(np.array(list(zip(*nora_comparison[key]))))
                #plt.plot(tmp[0], tmp[1], color = 'grey')
            return tmp[0], tmp[1], 'Multiplets (100 s, 5 yr.)'
        else:
            ps_pap = np.genfromtxt(f_path + 'icecube_misc/effective_areas_alerts/point_source_paper_lims.csv',
                            delimiter=', ')
            tmp = np.array(list(zip(*ps_pap)))
            lums = tmp[0] * self.no_evol_energy_density/self.energy_density #testing diff. spectrum, this is an approximation rn
            return np.log10(tmp[1]), np.log10(lums), '8 yr. point source'

    def load_signalness_array(self):
        r''' Load the signalness by index list, apply
        appropriate masks depending on which analysis is being run'''
        sigs_all = np.load(f_path + 'icecube_misc/effective_areas_alerts/sigs_by_ind.npy')
        if self.transient and self.delta_t == 1000.:
            msk_inds = np.array([60, 79, 228])
        elif self.transient:
            msk_inds = np.array([60])
        else:
            msk_inds = np.array([13, 32, 60, 83, 143, 147])
        sigs = np.delete(sigs_all, msk_inds)
        return sigs

    def brazil_bands(self, in_ts=False, rotated=False, compare=False, 
        with_result=False, result=None):
        r'''Two dimensional contour plot to show the sensitivity of the analysis
        in the luminosity-density plane, with full brazil bands, not just the 
        median upper limit

        Parameters
        ----------
            - in_ts (bool): TS value or binomial-p value
            - rotated (bool): Show typical Kowalski plot (false) or scale y-axis
            - compare (bool): compare to other icecube analyses
        '''
        if in_ts:
            if self.background_median_ts is None:
                self.get_overall_background_ts()
            if self.med_TS is None:
                self.get_med_TS()
        else:
            if self.background_median_p is None:
                self.get_overall_background_p()
            if self.med_p is None:
                self.get_med_p()
        fig, ax = plt.subplots(figsize=(8,5), dpi=200)
        fig.set_facecolor('w')
        X, Y = np.meshgrid(self.densities, self.plot_lumis)
        if rotated:
            Y *= X
        X = np.log10(X); Y = np.log10(Y)

        levels = np.array([16., 50., 84.])
        lower_bound_comp = self.lower_10 if in_ts else self.lower_10_p
        background_dist = self.stacked_ts if in_ts else self.stacked_p
        comp_factor = 1. if in_ts else -1.

        #Plot median UL (sensitivity)
        level = levels[len(levels)//2]
        linestyle='solid'
        reference_val = np.percentile(background_dist, level)
        if in_ts:
            sens_disc = comp_factor * (lower_bound_comp - reference_val)
        else:
            sens_disc = comp_factor * (np.log10(lower_bound_comp) - np.log10(reference_val))
        sens_alpha = 1. if not with_result else 0.5
        cs_ts = ax.contour(X, Y, sens_disc, colors=['k'], 
                        levels=[0.0], 
                        linewidths=2., linestyles=linestyle,
                        alpha=sens_alpha, zorder=6)
        text_loc = [(-9.9, 43.2)] if rotated else [(-9.9, 53.5)]
        ax.clabel(cs_ts, cs_ts.levels, inline=True, fmt='Sensitivity', fontsize=14,
            manual=text_loc)

        if with_result:
            if in_ts:
                upper_lim = comp_factor * (lower_bound_comp - result)
            else:
                upper_lim = comp_factor * (np.log10(lower_bound_comp) - np.log10(result))
            cs_ts_res = ax.contour(X, Y, upper_lim, colors=['k'], 
                        levels=[0.0], 
                        linewidths=2., linestyles=linestyle, zorder=6)
            text_loc = [(-10.4, 43.2)] if rotated else [(-10.5, 53.5)]
            ax.clabel(cs_ts_res, cs_ts.levels, inline=True, 
                fmt=r'Upper limit', fontsize=14, manual=text_loc,
                rightside_up=False, zorder=5)

        # central 68% containment band
        for band_num in range(len(levels)//2):
            lower_level = levels[band_num]
            upper_level = levels[-1*(band_num+1)]

            reference_val_lower = np.percentile(background_dist, lower_level)
            reference_val_upper = np.percentile(background_dist, upper_level)

            if in_ts:
                sens_disc_lower = comp_factor * (lower_bound_comp - reference_val_lower)
                sens_disc_upper = comp_factor * (lower_bound_comp - reference_val_upper)
            else:
                sens_disc_lower = comp_factor * (np.log10(lower_bound_comp) - np.log10(reference_val_lower))
                sens_disc_upper = comp_factor * (np.log10(lower_bound_comp) - np.log10(reference_val_upper))

            sens_disc = sens_disc_lower * sens_disc_upper
          
            cs_ts = ax.contour(X, Y, sens_disc, colors=['k'], 
                                  levels=[0.0], linewidths=1.5, linestyles='dashed',
                                  alpha=sens_alpha, zorder=5)
            csf = ax.contourf(X, Y, sens_disc, cmap=self.cmap, 
                                levels=[np.min(sens_disc), 0.0],
                                alpha=sens_alpha, zorder=5)
            text_loc = [(-9.0, 43.), (-9.0, 43.5)] if rotated else [(-9., 52.), (-9., 53.)]
            ax.clabel(cs_ts, cs_ts.levels, inline=True, fmt=r'$\pm 1 \sigma$', fontsize=14,
                manual=text_loc)
            
        xs = np.logspace(-11., -6., 1000)
        ys_median = self.energy_density / xs / self.seconds_per_year if self.transient else self.energy_density / xs
        
        if not rotated:
            plt.plot(np.log10(xs), np.log10(ys_median), color = sns.xkcd_rgb['dodger blue'], lw=1.5, zorder=1)
            for sig in ['one_sigma', 'two_sigma']:
                upper_factor = self.energy_density_uncertainty['plus_'+sig] / self.energy_density
                lower_factor = self.energy_density_uncertainty['minus_'+sig] / self.energy_density
                alpha = 0.45 if sig is 'two_sigma' else 0.75
                plt.fill_between(np.log10(xs), np.log10(ys_median * lower_factor), np.log10(ys_median * upper_factor), 
                        color = sns.xkcd_rgb['dodger blue'], alpha = alpha, lw=0.0, zorder=1)
        else:
            plt.plot(np.log10(xs), np.log10(ys_median*xs), color = sns.xkcd_rgb['dodger blue'], lw=1.5, zorder=1)
            for sig in ['one_sigma', 'two_sigma']:
                upper_factor = self.energy_density_uncertainty['plus_'+sig] / self.energy_density
                lower_factor = self.energy_density_uncertainty['minus_'+sig] / self.energy_density
                alpha = 0.45 if sig is 'two_sigma' else 0.75
                plt.fill_between(np.log10(xs), np.log10(ys_median * lower_factor * xs), np.log10(ys_median * upper_factor * xs), 
                        color = sns.xkcd_rgb['dodger blue'], alpha = alpha, lw=0.0, zorder=1)
        if rotated:
            plt.text(-10, np.log10(np.max(ys_median*xs*upper_factor)*1.1), 'Diffuse', 
                            color = sns.xkcd_rgb['dodger blue'], rotation=0, fontsize=18)
        else:
            diff_yloc = 53.0 if not self.transient else 53.5
            plt.text(-9, diff_yloc, 'Diffuse', color = sns.xkcd_rgb['dodger blue'], rotation=-28, fontsize=18)
        plt.xlim(-11., -6.)
        if rotated:
            plt.ylim(np.log10(np.min(ys_median*xs)*1e-2), np.log10(np.max(ys_median*xs)*5))
            plt.ylabel(self.scaled_lumi_label, fontsize = 22)
        else:
            if self.transient:
                plt.ylim(50., 55.5)
            else:
                plt.ylim(50., 55.)
            plt.ylabel(self.lumi_label, fontsize = 22)
        if self.transient:
            if self.delta_t == 1e3:
                time_window_str = r'$\pm 500$ s, '
            else:
                time_window_str = r'$\pm 1$ day, '
        else:
            time_window_str = 'Time integrated, '
        if compare:
            comp_rho, comp_en, comp_str = self.compare_other_analyses()
            if not rotated:
                plt.plot(comp_rho, comp_en, color = 'gray', lw=2.)
            else:
                plt.plot(comp_rho, comp_rho+comp_en, color = 'gray', lw=2.) #damn look at that log property
        custom_labs = [Line2D([0], [0], color = 'k', lw=2., label='This analysis (' + time_window_str + '{:.1f} yr.)'.format(self.data_years))]
        if compare:
            custom_labs.append(Line2D([0], [0], color='grey', lw=2., label=comp_str))
        leg_loc = 4 if rotated else 1
        plt.legend(handles=custom_labs, loc=leg_loc).set_zorder(20)
        plt.xlabel(self.density_label, fontsize = 22)
        title = self.lumi_str + ', ' + self.evol_str
        # plt.title(title)
        rot_str = '' if not rotated else 'rotated_'
        if self.save_figs:
            for ftype in ['.png', '.pdf']:
                file_path = self.savepath + rot_str + 'brazil_bands' \
                    + self.steady_str + '_' + self.evol \
                    + '_{}_years'.format(int(self.data_years)) \
                    + '_' + self.lumi
                plt.savefig(
                    file_path + ftype, 
                    bbox_inches='tight'
                    )
            if self.verbose:
                print("\t - saved plot to {}".format(file_path))

    def plot_background_binomial_p(self):
        fig, ax = plt.subplots(dpi=200)
        p_bins = np.logspace(-5., 0., 31)
        plt.hist(self.stacked_p, bins=p_bins)
        plt.xlabel('Binomial p-value')
        plt.ylabel('$N$')
        if self.transient:
            if self.delta_t == 1e3:
                time_window_str = '$\pm 500$ s, '
            else:
                time_window_str = '$\pm 1$ day, '
            plt.title(
                r'$\Delta T = $' + time_window_str \
                + '{} yrs alerts'.format(self.data_years)
                )
        else:
            plt.title(
                r'Time integrated' 
                + ', {} yrs alerts'.format(self.data_years)
                )
        plt.loglog()
        if self.save_figs:
            for ftype in ['.png', '.pdf']:
                file_path = self.savepath + 'binom_p_distribution' \
                    + self.steady_str \
                    + '_{}_years'.format(int(self.data_years))
                plt.savefig(
                    file_path + ftype, 
                    bbox_inches='tight'
                )
            if self.verbose:
                print("\t - saved plot to {}".format(file_path))

    def upper_limit_plot(self):
        pass

    def fit_coverage_plot(self, dens, lumi):
        pass
