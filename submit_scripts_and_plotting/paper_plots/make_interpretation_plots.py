'''Script to make all interpretation based plots
for the alert followup fast response paper'''

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as st
from francis import utils
from results_utils import *
from francis.universe.universe_plotter import UniversePlotter

mpl.style.use(utils.initialize_mpl_style())
base_output = os.path.abspath(utils.get_francis_path() + '../figures/paper_plots/') + '/'
palette = sns.light_palette((210, 90, 60), input="husl", n_colors=12)

def get_universe_dict():
    print("Initializing all of the results from three different time scales")
    unis = dict()
    binom_scans = dict()
    for delta_t in [1000., 172800., None]:
        uni = UniversePlotter(delta_t, 9.6, 'SC', 'MD2014SFR', save=False)
        pp, tmps = uni.get_overall_background_p(with_tmps=True)
        truth, pvals, arr = get_true_binomial_p(delta_t=delta_t)

        ppost = np.count_nonzero(uni.stacked_p <= truth) / float(len(uni.stacked_p))
        
        unis[delta_t] = uni
        binom_scans[delta_t] = {
            'truth': truth, 
            'pvals': pvals, 
            'arr': arr, 
            'tmps': tmps, 
            'p_post': ppost
            }
        if delta_t is None:
            print('\t- Done with time integrated initialization')
        else:
            print(f'\t- Done with time window {delta_t:.1e}')
    
    return unis, binom_scans

def binomial_scan_plots(unis, binom_scans):
    print('Making binomial scan plots . . . ')
    for delta_t in [1000., 172800., None]:
        uni = unis[delta_t]
        scan = binom_scans[delta_t]

        fig, ax = plt.subplots(dpi=200)

        plt.plot(
            np.r_[:len(scan['arr'])] + 1, 
            scan['arr'], 
            color=palette[-1], 
            lw=2., 
            zorder=20
            )
        thresh_ps = st.norm.sf([0., 1., 2., 3.])
        percents = np.percentile(uni.stacked_p, thresh_ps * 100.)
        for perc, lab in zip(percents, ['Median', r'$+1\sigma$', r'$+2\sigma$', r'$+3\sigma$']):
            plt.axhline(perc, ls='-.', color = sns.xkcd_rgb['dodger blue'], lw=0.85, alpha=0.8)
            plt.text(160, perc*1.1, lab, color = sns.xkcd_rgb['dodger blue'], alpha=0.8)
        plt.yscale('log')
        plt.xlabel('Number of sources')
        plt.ylabel('Binomial p-value')
        plt.ylim(1.1e-5, 1.2e0)
        delta_t_str = 'steady' if delta_t is None else f'delta_t_{delta_t:.0e}'.replace('+', '')
        for ftype in ['.png', '.pdf']:
            plt.savefig(
                base_output + f'binomial_scan_{delta_t_str}{ftype}',
                bbox_inches='tight'
            )
        plt.close()
    print('\t- Done with binomial scan plots')
        
def binomial_distribution_plots(unis, binom_scans):
    print('Making binomial distribution plots')
    for delta_t in [1000., 172800., None]:
        uni = unis[delta_t]
        scan = binom_scans[delta_t]

        uni.plot_background_binomial_p()
        plt.axvline(scan['truth'], color='k', zorder=20, lw=2.)
        plt.text(3e-4, 1e2, f"p-val: {scan['p_post']:.3f}")
        plt.xlim(plt.xlim()[1], plt.xlim()[0])

        delta_t_str = 'steady' if delta_t is None else f'delta_t_{delta_t:.0e}'.replace('+', '')
        for ftype in ['.png', '.pdf']:
            plt.savefig(
                base_output + f'binomial_distribution_{delta_t_str}{ftype}',
                bbox_inches='tight'
            )
        plt.close()
    print('\t- Done with binomial distribution plots')

def sensitivity_and_upper_limits(unis, binom_scans):
    print('Calculating upper limits and making plots')
    for delta_t in [1000., 172800., None]:
        uni = unis[delta_t]
        scan = binom_scans[delta_t]
        delta_t_str = 'steady' if delta_t is None else f'delta_t_{delta_t:.0e}'.replace('+', '')

        uni.brazil_bands(rotated=True)
        for ftype in ['.png', '.pdf']:
            plt.savefig(
                base_output + f'sensitivity_rotated_{delta_t_str}{ftype}',
                bbox_inches='tight'
            )
        plt.close()

        for rotated in [True, False]:
            uni.brazil_bands(
                rotated=rotated, 
                with_result=True, 
                result=scan['truth']
                )
            rotated_str = '_rotated' if rotated else ''
            for ftype in ['.png', '.pdf']:
                plt.savefig(
                    base_output + f'upper_limit{rotated_str}_{delta_t_str}{ftype}',
                    bbox_inches='tight'
                )
            plt.close()
    print('\t- Done with upper limit plots')

if __name__ == "__main__":
    unis, binom_scans = get_universe_dict()
    binomial_scan_plots(unis, binom_scans)
    binomial_distribution_plots(unis, binom_scans)
    sensitivity_and_upper_limits(unis, binom_scans)