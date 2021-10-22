import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os, sys
import matplotlib as mpl
import healpy as hp
from matplotlib.lines import Line2D
import pickle

import francis.sensitivity_fit_functions as sff
from francis import utils
utils.initialize_mpl_style()
f_path = utils.get_francis_path

savepath = f_path + '../figures/followup_plots/transient/'

###############################################################################
##########              Load plotting parameters              #################
###############################################################################
import seaborn as sns
palette = ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f']

utils.initialize_mpl_style()
mpl.rcParams['axes.linewidth'] = 2.
mpl.rcParams['figure.dpi'] = 120

sens_cols = [sns.xkcd_rgb['navy blue'], sns.xkcd_rgb['navy green'], sns.xkcd_rgb['deep magenta']]
disc_cols = [sns.xkcd_rgb['light navy blue'], sns.xkcd_rgb['lighter green'], sns.xkcd_rgb['light magenta']]

###############################################################################
##########              Sensitivity vs. declination           #################
###############################################################################

delta_ts = np.array([1000., 2*86400.])
sensitivities = np.zeros((2, 276))
discoveries = np.zeros((2, 276))

unsmeared_sensitivities = np.zeros((2, 276))
unsmeared_discoveries = np.zeros((2, 276))
for jj, delta_t in enumerate(delta_ts):
    sensitivities[jj], discoveries[jj] = sff.find_all_sens(delta_t, smear=True, disc_conf=0.5,
                                                        disc_thresh=1.-0.0013, verbose=False)

decs_by_ind = np.load(f_path + 'icecube_misc/effective_areas_alerts/decs_by_ind.npy')[1]
time_labels = [r'$10^3$ s', '2 days']
time_integrated = {'sens': [], 'disc': [], 'sinDec': []}
sinDecs = np.linspace(-0.95, 0.95, 39)
decs = np.arcsin(sinDecs) * 180. / np.pi
for dec in decs[:]:
    try:
        with open(f_path + 'icecube_misc/combined_tracks_2.5/results/EstSens_PointSourceTracks_v003p00_10yrs_dec_{:.2f}_degrees.pickle'.format(dec), 'rb') as f:
            data = pickle.load(f, encoding='latin1')
    except:
        try:
            with open(f_path + 'icecube_misc/combined_tracks_2.5/results/EstSens_PointSourceTracks_v003p00_10yrs_dec_{:.1f}_degrees.pickle'.format(dec), 'rb') as f:
                data = pickle.load(f, encoding='latin1')
        except Exception as e:
            continue
    time_integrated['sens'].append(data['sensitivity_flux'])
    time_integrated['disc'].append(data['discovery_flux'])
    time_integrated['sinDec'].append(np.sin(dec * np.pi / 180.))
    #except:
    #    continue
time_integrated['sens'] = np.array(time_integrated['sens'])
time_integrated['disc'] = np.array(time_integrated['disc'])

for ii, delta_t in enumerate(delta_ts[:]):
    try:
        fig = plt.figure(dpi=200)
        fig.set_facecolor('white')
        with open(f_path + 'icecube_misc/ideal_ps_sens/ideal_ps_sensitivity_deltaT_{:.2e}_50CL.pkl'.format(delta_t / 86400.), 'rb') as f:
            ideal = pickle.load(f, encoding='latin1')
        msk = sensitivities[ii]*delta_t*1e6 < 1e2
        if delta_t == 1000.:
            bad_ind_msk = np.array([x not in [92, 95, 198] for x in np.linspace(0, 275, 276)])
        else:
            bad_ind_msk = np.array([x not in [198] for x in np.linspace(0, 275, 276)])
        msk *= bad_ind_msk
        plt.scatter(np.sin(decs_by_ind)[msk], sensitivities[ii][msk]*delta_t*1e6, color=sens_cols[ii], 
                    marker='^', label = 'Alert sensitivity', s=10)
        plt.scatter(np.sin(decs_by_ind)[msk], discoveries[ii][msk]*delta_t*1e6, color=sens_cols[ii], 
                    marker='2', linestyle='--', label = 'Alert discovery', s=10)
        plt.plot(ideal['sinDec'], np.array(ideal['sensitivity'])*delta_t*1e6, lw=3, ls='-', 
                 color=sens_cols[ii], label = 'P.S. sensitivity', alpha=0.7)
        plt.plot(ideal['sinDec'], np.array(ideal['discovery'])*delta_t*1e6, lw=3, ls='--', 
                 color=sens_cols[ii], label = r'P.S. $3\sigma$ discovery' +'\n\t(50\% CL)', alpha=0.7)
        plt.plot(time_integrated['sinDec'], time_integrated['sens'] * 86400. * 10. * 365. * 1e6, color = sns.xkcd_rgb['light grey'],
                alpha = 0.7, ls='-')
        plt.plot(time_integrated['sinDec'], time_integrated['disc'] * 86400. * 10. *365.* 1e6, 
                 color = sns.xkcd_rgb['light grey'],
                alpha = 0.7, ls='--')

        plt.grid(which='both', alpha=0.2, zorder=1)
        plt.yscale('log')
        plt.legend(loc=1, ncol=1, frameon=True, fontsize=14)
        plt.xlabel(r'$\sin \delta$')
        plt.title(r'$\Delta T = $' + time_labels[ii])
        plt.ylabel(r'$E^2 \frac{dN_{\nu+\bar{\nu}}}{dEdA} \Big|_{\mathrm{1 TeV}}$ (GeV cm$^{-2}$)')
        plt.ylim(3e-2, 3e2)
        plt.text(0.03, 1e0, '10 yr. time-integrated', color = sns.xkcd_rgb['light grey'])
        
    except Exception as e:
        print(e)
        pass
    plt.savefig(savepath + 'alert_event_sens_vs_dec_delta_t_{:.1e}_smeared.png'.format(delta_t), 
               bbox_inches='tight', dpi=200)
    plt.close()

###############################################################################
##########              Background trial plots                #################
###############################################################################
bins = np.linspace(0., 20., 21)

bg = np.array(sff.background_distribution(5, 2.*86400., smear=True))
plt.hist(bg, bins=bins, histtype='step', lw=2., label='With systematics')
plt.yscale('log')
plt.xlabel('TS')
plt.ylabel('$N$')
plt.title(r'$\Delta T = 2$ days')
plt.savefig(savepath + 'FRA_alert_event_background_trials_sample.png', 
               bbox_inches='tight', dpi=200)
plt.close()

bins = np.linspace(0., 20., 21)
problem_inds = {1000.: [60., 79., 228], 2.*86400: [60]}

for delta_t in np.array([1000., 2.*86400.]):
    fig, aaxs = plt.subplots(ncols=12, nrows=23, figsize=(17,35), 
                            sharey=True, sharex=True)
    used_axs = []
    axs = aaxs.ravel()
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    for ind in range(276):
        if ind % 25 == 0:
            print(ind, end=' ')
        if ind in problem_inds[delta_t]:
            continue
        try:
            bg = np.array(sff.background_distribution(ind, delta_t, smear=True))
            axs[ind].hist(bg, bins=bins, histtype='step', lw=2.)
            axs[ind].set_yscale('log')
            axs[ind].set_ylim(0.8e0, 2e4)
            if ind / 12 == 22:
                axs[ind].set_xlabel('TS')
            if ind % 12 == 0:
                axs[ind].set_ylabel('$N$')
            used_axs.append(axs[ind])
        except (IOError, ValueError, IndexError) as err:
            print(err)
            # pass
    for ax in axs:
        if ax not in used_axs:
            ax.set_visible(False)
    fig.suptitle(r'$\Delta T = $' + '{:.1e} s'.format(delta_t), y=0.89)
    plt.savefig(savepath + 'alert_event_background_trials_delta_t_{:.1e}_smeared.png'.format(delta_t), 
               bbox_inches='tight', dpi=200)
    plt.close()

###############################################################################
##########              Fit validation plots                  #################
###############################################################################

sff.ns_fits_contours_plot(272, 1000., smear=True, show=False, legend=True)
plt.ylim(0., 40)
plt.savefig(savepath + 'FRA_alert_event_nsignal_fitting_sample.png', 
               bbox_inches='tight', dpi=200)

for delta_t in np.array([1000., 2.*86400.]):
    fig, aaxs = plt.subplots(ncols=12, nrows=23, figsize=(17,35), 
                            sharey=True, sharex=True)
    used_axs = []
    axs = aaxs.ravel()
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    for ind in range(276):
        xlab = ind / 12 == 22
        ylab = ind % 12 == 0
        if ind % 25 == 0:
            print(ind, end=' ')
        try:
            sff.ns_fits_contours_plot(ind, delta_t, smear=True, show=False,
                                 ax=axs[ind], xlabel=xlab, ylabel=ylab, legend=False)
            used_axs.append(axs[ind])
        except (IOError, ValueError, IndexError) as err:
            print(err)
    for ax in axs:
        if ax not in used_axs:
            ax.set_visible(False)
    fig.suptitle(r'$\Delta T = $' + '{:.1e} s'.format(delta_t), y=0.89)
    plt.savefig(savepath + 'alert_event_fit_panels_delta_t_{:.1e}_smeared.png'.format(delta_t), 
               bbox_inches='tight', dpi=200)
    plt.close()

legend_elements = [Line2D([0], [0], color='gray', ls='--', label='Injected'),
                   Line2D([0], [0], marker='x', color='gray', label='Fit (50\% cont.)', ls='')]

for delta_t, delta_t_str in zip([1000., 172800.], 
                                [r'$\Delta T = 10^3$ s', r'$\Delta T = 1.7\times 10^5$ s']):
    bias, spread = sff.fitting_bias_summary(delta_t)
    for ii, sig in enumerate([2.0, 5.0, 10.]):
        msk = np.array(bias[sig]) != 0.
        plt.axhline(sig, ls='--', color=palette[ii])
        plt.errorbar(np.r_[0:276:1][msk], np.array(bias[sig])[msk], yerr=np.array(spread[sig])[msk],
                    marker = 'x', ls='', color=palette[ii], lw=1.)
    plt.ylim(0.0, 14)
    plt.title(delta_t_str)
    plt.xlabel('v2 alert number')
    plt.ylabel(r'$n_s$')
    plt.legend(handles=legend_elements, ncol=2, loc=2)
    plt.savefig(savepath + 'alert_event_fit_summary_delta_t_{:.1e}_smeared.png'.format(delta_t), 
               bbox_inches='tight', dpi=200)
    plt.close()

###############################################################################
##########              Sensitivity curve plots                  ##############
###############################################################################

for delta_t in [1e3, 2.*86400.]:
    fig, aaxs = plt.subplots(ncols=12, nrows=23, figsize=(17,35), 
                            sharey=True)
    used_axs = []
    axs = aaxs.ravel()
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    for ind in range(276):
        try:
            sff.sensitivity_curve(ind, delta_t, ax = axs[ind], in_ns=False, 
                              smear=True, text=False, legend=False)
            used_axs.append(axs[ind])
        except (IOError, ValueError, IndexError) as err:
            pass
    for ax in axs:
        if ax not in used_axs:
            ax.set_visible(False)
    fig.suptitle(r'$\Delta T = $' + '{:.1e} s'.format(delta_t), y=0.89)
    plt.savefig(savepath + 'alert_event_senscurve_summary_delta_t_{:.1e}_smeared.png'.format(delta_t), 
               bbox_inches='tight', dpi=200)
    plt.close()

for delta_t in [1e3, 2.*86400.]:
    fig, aaxs = plt.subplots(ncols=12, nrows=23, figsize=(17,35), 
                            sharey=True)
    used_axs = []
    axs = aaxs.ravel()
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    for ind in range(276):
        try:
            sff.sensitivity_curve(ind, delta_t, ax = axs[ind], in_ns=False, 
                              smear=True, text=False, legend=False,
                                 threshold=1.0-0.0013, conf_lev=0.5)
            used_axs.append(axs[ind])
        except (IOError, ValueError, IndexError) as err:
            pass
    for ax in axs:
        if ax not in used_axs:
            ax.set_visible(False)
    fig.suptitle(r'$\Delta T = $' + '{:.1e} s'.format(delta_t), y=0.89)
    plt.savefig(savepath + 'alert_event_disccurve_summary_delta_t_{:.1e}_smeared.png'.format(delta_t), 
               bbox_inches='tight', dpi=200)
    plt.close()

