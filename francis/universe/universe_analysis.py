import numpy as np
from glob import glob
import pandas as pd
import scipy.stats as st
import pickle
import csv
import sys
from francis.universe.transient_universe import TransientUniverse, SteadyUniverse
from francis.universe.transient_universe import *
from francis import utils
f_path = utils.get_francis_path()

eff_area_path = f_path + 'icecube_misc/effective_areas_alerts/'

# Commented paths point to original file locations
bg_trials = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/'
signal_trials = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/fits/'
# bg_trials = '/data/ana/analyses/NuSources/2021_v2_alert_stacking_FRA/analysis_trials/bg/'
# signal_trials = '/data/ana/analyses/NuSources/2021_v2_alert_stacking_FRA/analysis_trials/fits/'

class UniverseAnalysis():
    r'''Given cosmological parameters, calculate the expected TS distribution
        from triggering short timescale analyses on alert events'''

    def __init__(self, lumi, evol, density, diffuse_flux_norm, diffuse_flux_ind,
                    **kwargs):
        self.lumi = lumi
        self.evol = evol
        self.density = density
        self.diffuse_flux_norm = diffuse_flux_norm
        self.diffuse_flux_ind = diffuse_flux_ind
        self.deltaT = kwargs.pop('deltaT', None)
        self.sigma = kwargs.pop('sigma', 1.0)
        self.transient = True if self.deltaT is not None else False
        if self.deltaT is not None:
            kwargs['timescale'] = self.deltaT
        self.seed = kwargs.pop('seed', 1234)
        if self.transient:
            self.universe = TransientUniverse(self.lumi, self.evol, self.density,
                self.diffuse_flux_norm, self.diffuse_flux_ind, seed=self.seed, sigma=self.sigma,
                **kwargs)
        else:
            self.universe = SteadyUniverse(self.lumi, self.evol, self.density,
                self.diffuse_flux_norm, self.diffuse_flux_ind, seed=self.seed, sigma=self.sigma,
                **kwargs)
        self.smear = kwargs.pop('smeared', True)
        self.smear_str = 'smeared/' if self.smear else 'norm_prob/'
        self.verbose = kwargs.pop('verbose', False)
        self.rng = np.random.RandomState(self.seed)
        self.initialize_universe()

    def print_analysis_info(self):
        r'''Print a message with info about the source once
        the analysis is running'''
        analysis_name = 'Alert event interpretation'
        int_str = '*'*80
        int_str += '\n*' + ' '*78 + '*\n'
        int_str += '*' + ' '*((78-len(analysis_name))//2) + analysis_name + ' '*((78-len(analysis_name))//2 + len(analysis_name)%2) + '*'
        int_str += '\n*' + ' '*78 + '*\n'
        int_str += '*'*80 + '\n'
        int_str += '  '*5 + 'Density: {:.1e}'.format(self.density)
        int_str += '  '*7 + 'Luminosity: {}'.format(self.lumi) + '\n'
        int_str += '  '*5 + 'Evolution: {}'.format(self.evol) 
        time_str = 'Steady' if not self.transient else '{:.1e} s'.format(self.deltaT)
        int_str += '  '*5 + 'Timescale: {}'.format(time_str) + '\n'
        int_str += '  '*5 + 'Diffuse gamma: {:.1f}'.format(self.diffuse_flux_ind)
        int_str += '  '*6 + 'Smearing: {}'.format(self.smear)
        int_str += '\n\n'
        print(int_str)

    #@profile
    def initialize_universe(self):
        """Simulate sources with the given cosmological parameters,
        also find the alert events as well as the additional injected
        events
        """
        if self.verbose:
            print("Simulating universe with specified cosmological parameters")
        self.universe.create_universe()
        self.universe.find_alerts()
        self.universe.find_alert_skymaps()
        self.universe.additional_signal_events()

    #@profile
    def make_alerts_dataframe(self):
        """
        Reformat the results from the simulation into a dataframe
        """
        alerts = {'signalness': [], 'declination': [], 'background': [], 
          'skymap_ind': [], 'stream': [], 'skymap_dec': [],
         'extra_evs': []}
        for k in self.universe.bg_alerts.keys():
            if self.universe.bg_alerts[k][0] > 0:
                alerts['signalness'].extend(self.universe.bg_alerts[k][1])
                alerts['declination'].extend(self.universe.bg_alerts[k][2])
                alerts['background'].extend([True]*self.universe.bg_alerts[k][0])
                alerts['skymap_ind'].extend(self.universe.bg_alerts[k][4]) 
                alerts['skymap_dec'].extend(self.universe.bg_alerts[k][3]) 
                alerts['stream'].extend([k]*self.universe.bg_alerts[k][0])
                alerts['extra_evs'].extend([0]*self.universe.bg_alerts[k][0])
        for k in self.universe.sig_alerts.keys():
            for jj in range(len(self.universe.sig_alerts[k])):
                if self.universe.sig_alerts[k][jj][0] == 0:
                    continue
                else:
                    alerts['signalness'].append(self.universe.sig_alerts[k][jj][1][0])
                    alerts['declination'].append(np.radians(self.universe.sources['dec'][jj]))
                    alerts['background'].append(False)
                    alerts['skymap_ind'].append(self.universe.skymaps[k][jj][1])
                    alerts['skymap_dec'].append(self.universe.skymaps[k][jj][0])
                    alerts['stream'].append(k)
                    alerts['extra_evs'].append(self.universe.extra_events[k][jj])
        alerts = pd.DataFrame(alerts)
        self.alert_df = alerts

    #@profile
    def reinitialize_universe(self):
        """Change the seed and reinitialize everything"""
        if self.verbose:
            print("Recreating universe for more trials, updating seed")
        self.seed = 1 if self.seed is None else self.seed + 1
        self.rng = np.random.RandomState(self.seed)
        self.universe.seed = self.seed
        self.universe.create_universe()
        self.universe.find_alerts()
        self.universe.find_alert_skymaps()
        self.universe.additional_signal_events()

    #@profile
    def calculate_ts(self, only_gold = False, calc_p=True):
        """
        Based off of the additional injected events, sample
        trials and calculate the final test statistics
        """
        ts, sigs, ps = [], [], []
        self.alert_df['TS'] = [None] * len(self.alert_df['background'])
        self.alert_df['pval'] = [None] * len(self.alert_df['background'])
        for index, alert in self.alert_df.iterrows():
            if alert['background']:
                if calc_p:
                    t, p = self.background_alert_trials(alert['skymap_ind'], calc_p=calc_p)
                    ts.append(t); ps.append(p)
                else:
                    ts.append(self.background_alert_trials(alert['skymap_ind'], calc_p=calc_p))
                sigs.append(alert['signalness'])
                self.alert_df.loc[self.alert_df.index == index, 'TS'] = ts[-1]
                if calc_p:
                    self.alert_df.loc[self.alert_df.index == index, 'pval'] = ps[-1]
            else:
                if calc_p:
                    t, p = self.signal_alert_trials(alert['skymap_ind'], alert['extra_evs'], calc_p=calc_p)
                    ts.append(t); ps.append(p)
                else:
                    ts.append(self.signal_alert_trials(alert['skymap_ind'], alert['extra_evs'], calc_p=calc_p))
                sigs.append(alert['signalness'])
                self.alert_df.loc[self.alert_df.index == index, 'TS'] = ts[-1]
                if calc_p:
                    self.alert_df.loc[self.alert_df.index == index, 'pval'] = ps[-1]
        ts, sigs = np.array(ts), np.array(sigs)
        if only_gold:
            gold = []
            for index, alert in self.alert_df.iterrows():
                if 'gold' in alert['stream']:
                    gold.append(True)
                else:
                    gold.append(False)
            gold = np.array(gold)
            ts, sigs = ts[gold], sigs[gold]
        TS = np.sum(sigs * ts) / sigs.size
        self.TS = TS
        return TS

    #@profile        
    def background_alert_trials(self, ind, calc_p=True):
        """If an alert is a background alert, sample from the background
        trials
        """
        if self.transient:
            trials_file = glob(bg_trials + self.smear_str + 'index_{}_*_time_{:.1f}.pkl'.format(ind, self.deltaT))[0]
            if sys.version[0] == '3':
                trials = np.load(trials_file, allow_pickle=True, encoding='latin1')
            else:
                trials = np.load(trials_file)
            ts = self.rng.choice(trials['ts_prior'])
            if calc_p:
                if ts == 0:
                    pval = 1.0
                else:
                    pval = float(np.count_nonzero(np.array(trials['ts_prior']) >= ts)) / np.array(trials['ts_prior']).size
                    if pval == 0.:
                        pval = 1./np.array(trials['ts_prior']).size
        else:
            fs = glob(bg_trials + self.smear_str + 'index_{}_*_steady.pkl'.format(ind))
            trials = np.load(fs[0], allow_pickle=True, encoding='latin1')
            ts = self.rng.choice(trials['TS'])
            if calc_p:
                if ts == 0:
                    pval = 1.0
                else:
                    pval = float(np.count_nonzero(trials['TS'] >= ts)) / trials['TS'].size
                    if pval == 0.:
                        pval = 1./np.array(trials['ts_prior']).size
        del trials
        if calc_p:
            return ts, pval
        else:
            return ts

    #@profile
    def signal_alert_trials(self, ind, N, calc_p = True):
        """If alerts are signal alerts and have additional injected
        events, sample the relevant signal trials"""
        if N == 0:
            ts = self.background_alert_trials(ind, calc_p = False)
        else:
            if self.transient:
                trials_file = glob(signal_trials + self.smear_str + 'index_{}_*_time_{:.1f}.pkl'.format(ind, self.deltaT))[0]
                if sys.version[0] == '3':
                    trials = np.load(trials_file, allow_pickle=True, encoding='latin1')
                else:
                    trials = np.load(trials_file)
            else:
                fs = glob(signal_trials + self.smear_str + 'index_{}_*_steady_gamma_2.5.pkl'.format(ind))
                t_file = fs[0]
                if sys.version[0] == '3':
                    trials = np.load(t_file, allow_pickle=True, encoding='latin1')
                else:
                    trials = np.load(t_file)
                #trials = np.load(signal_trials + 'index_{}_steady.pkl'.format(ind))
            ns_key = 'true_ns' if self.transient else 'inj_nsignal'
            if N <= 10:
                inds = np.argwhere(np.array(trials[ns_key]) == N).flatten()
                if len(inds) == 0:
                    inds = np.argwhere(np.abs(np.array(trials[ns_key]) - N) < 4).flatten()
                    if len(inds) == 0:
                        if self.verbose:
                            print("No trials near {}".format(N))
                        inds = np.argmin(np.abs(np.array(trials[ns_key]) - N)).flatten()    
            else:
                inds = np.argwhere(np.abs(np.array(trials[ns_key]) - N) < 10).flatten()
                if len(inds) == 0:
                    if self.verbose:
                        print("NO TRIALS WITH {} INJECTED EVENTS".format(N))
                    inds = np.argwhere(np.array(trials[ns_key]) == np.max(trials[ns_key])).flatten()
            ts = np.array(trials['ts'])[inds] if self.transient else np.array(trials['TS'])[inds]
            ts = self.rng.choice(ts)
            del trials
        if calc_p:
            pval = self.calculate_trial_pvalue(ind, ts)
            return ts, pval
        else:
            return ts

    #@profile
    def calculate_trial_pvalue(self, ind, TS):
        """Find a p-value from TS value for a specific alert"""
        if TS == 0:
            return 1.
        if self.transient:
            trials_file = glob(bg_trials + self.smear_str + 'index_{}_*_time_{:.1f}.pkl'.format(ind, self.deltaT))[0]
            if sys.version[0] == '3':
                trials = np.load(trials_file, allow_pickle=True, encoding='latin1')
            else:
                trials = np.load(trials_file)
            pval = float(np.count_nonzero(np.array(trials['ts_prior']) >= TS)) / np.array(trials['ts_prior']).size
            if pval == 0.:
                pval = 1./np.array(trials['ts_prior']).size
        else:
            fs = glob(bg_trials + self.smear_str + 'index_{}_*_steady.pkl'.format(ind))
            kwargs = {} if not sys.version[0] == '3' else {'encoding': 'latin1', 'allow_pickle': True}
            trials = np.load(fs[0], **kwargs)
            pval = float(np.count_nonzero(trials['TS'] >= TS)) / trials['TS'].size
            if pval == 0.:
                pval = 1./trials['TS'].size
            del trials
        return pval

    #@profile
    def calculate_binomial_pvalue(self, only_gold=False):
        """With a list of p-values, do a scan over the possible
        number of sources to find the binomial p-value"""
        if self.TS is None:
            self.calculate_ts(only_gold = only_gold, calc_p=True) 
        plist = self.alert_df['pval']
        if only_gold:
            stream_msk = self.alert_df['stream']
            stream_msk = ['gold' in al for al in stream_msk] 
            plist = plist[stream_msk]
        obs_p = 1.
        plist = sorted(plist)
        for i, p in enumerate(plist):
            
            tmp = st.binom_test(i+1, len(plist), p, alternative='greater')
            if tmp < obs_p and tmp != 0.0:
                if tmp == 0.0:
                    print("WHY DOES THE BINOMIAL VALUE EQUAL ZERO")
                obs_p = tmp
        self.binom_p = obs_p
        return obs_p 
