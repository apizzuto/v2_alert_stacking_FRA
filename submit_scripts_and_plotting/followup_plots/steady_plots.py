from francis.time_integrated_scripts import steady_sensitivity_fits as ssf
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import seaborn as sns
import pickle

from francis import utils
utils.initialize_mpl_style()
f_path = utils.get_francis_path()
mpl.rcParams['text.usetex'] = True

alerts = pd.read_csv(
    f_path + 'icecube_misc/alert_dataframe.csv'
    )

##############################################################################
#############  For each alert, make a few summary plots          #############
##############################################################################
output_dir = f_path + '../figures/followup_plots/steady/time_integrated_summary_plots/'

# Make all of the summary plots. This takes a looooooong time (O(hours))
just_make_one = True
inds = [0] if just_make_one else list(range(len(alerts)))
for ind in inds:
    if ind % 50 == 0:
        print(ind, end=' ')
    try:
        event_header = ssf.panel_plot_with_text(ind, smear=True, return_info=True)
        plt.savefig(output_dir \
            + f'steady_panel_plot_with_text_{ind}_run' \
            + f'_{event_header["RUNID"]}_event_{event_header["EVENTID"]}.png',
            bbox_inches='tight', dpi=120)
        plt.close()
    except Exception as e:
        print(ind, e)
        plt.close()

##############################################################################
#############  Plot out sensitivity vs. declination              #############
##############################################################################
time_integrated = {'sens': [], 'disc': [], 'sinDec': []}
ideal = {'sens': [], 'disc': [], 'sinDec': []}
sinDecs = np.linspace(-0.95, 0.95, 39)
decs = np.arcsin(sinDecs) * 180. / np.pi
for dec in decs[:]:
    try:
        with open(f_path + 'icecube_misc/combined_tracks_2.5/results/EstSens_PointSourceTracks_v003p00_10yrs_dec_{:.2f}_degrees.pickle'.format(dec), 'rb') as f:
            data = pickle.load(f, encoding='latin1')
    except Exception as e:
        try:
            with open(f_path + 'icecube_misc/combined_tracks_2.5/results/EstSens_PointSourceTracks_v003p00_10yrs_dec_{:.1f}_degrees.pickle'.format(dec), 'rb') as f:
                data = pickle.load(f, encoding='latin1')
        except Exception as ee:
            continue
    time_integrated['sens'].append(data['sensitivity_flux'])
    time_integrated['disc'].append(data['discovery_flux'])
    time_integrated['sinDec'].append(np.sin(dec * np.pi / 180.))
    try:
        with open(f_path + 'icecube_misc/combined_tracks_2.5/results/EstSens_GFUOnline_v001p02_10yrs_dec_{:.2f}_degrees.pickle'.format(dec), 'rb') as f:
            data = pickle.load(f, encoding='latin1')
    except Exception as e:
        try:
            with open(f_path + 'icecube_misc/combined_tracks_2.5/results/EstSens_GFUOnline_v001p02_10yrs_dec_{:.1f}_degrees.pickle'.format(dec), 'rb') as f:
                data = pickle.load(f, encoding='latin1')
        except:
            continue
    ideal['sens'].append(data['sensitivity_flux'])
    ideal['disc'].append(data['discovery_flux'])
    ideal['sinDec'].append(np.sin(dec * np.pi / 180.))
    
time_integrated['sens'] = np.array(time_integrated['sens'])
time_integrated['disc'] = np.array(time_integrated['disc'])

ideal['sens'] = np.array(ideal['sens'])
ideal['disc'] = np.array(ideal['disc'])

alert_sens = {'sin_dec': [], 'sens': [], 'disc': []}
for ii in range(276):
    if ii % 25 == 0:
        print(ii, end=' ')
    try:
        sens = ssf.steady_sensitivity(ii, gamma=2.5, smear=True)['sens']
        disc = ssf.steady_sensitivity(ii, gamma=2.5, smear=True, conf_lev=0.5, thresh=0.99865)['sens']
        skymap, header = ssf.load_map(ii)
        dec = header['DEC']*np.pi / 180.
    except Exception as e:
        print(e)
        print("SENSITIVITY FIT FAILED FOR INDEX :{}".format(ii))
        sens = 0.0
        disc = 0.0
        dec = 0.0
    alert_sens['sens'].append(sens)
    alert_sens['disc'].append(disc)
    alert_sens['sin_dec'].append(np.sin(dec))
for k in alert_sens.keys():
    alert_sens[k] = np.array(alert_sens[k])

utils.initialize_mpl_style()
sens_cols = [sns.xkcd_rgb['navy blue'], sns.xkcd_rgb['navy green'], sns.xkcd_rgb['deep magenta']]
ii=2
try:
    fig = plt.figure(dpi=200)
    fig.set_facecolor('white')

    plt.plot(time_integrated['sinDec'], time_integrated['sens'] * 86400. * 10. * 365. * 1e6, 
             color = sns.xkcd_rgb['light grey'],
            alpha = 0.7, ls='-', zorder=4)
    plt.plot(time_integrated['sinDec'], time_integrated['disc'] * 86400. * 10. *365.* 1e6, 
             color = sns.xkcd_rgb['light grey'],
            alpha = 0.7, ls='--', zorder=4)
    
    plt.plot(ideal['sinDec'], ideal['sens'] * 86400. * 8.607 * 365. * 1e6, 
             color = sens_cols[ii], lw=3, label = 'P.S. sensitivity',
            alpha = 0.7, ls='-', zorder=4)
    plt.plot(ideal['sinDec'], ideal['disc'] * 86400. * 8.607 *365.* 1e6, 
             color = sens_cols[ii], lw=3, label = r'P.S. $3\sigma$ discovery' +'\n\t(50\% CL)',
            alpha = 0.7, ls='--', zorder=4)
    
    plt.scatter(alert_sens['sin_dec'], alert_sens['sens']*86400*8.6*365.*1e6, color=sens_cols[ii], 
                    marker='^', label = 'Alert sensitivity', s=10, zorder=5)
    plt.scatter(alert_sens['sin_dec'], alert_sens['disc']*86400*8.6*365.*1e6, color=sens_cols[ii], 
                    marker='2', linestyle='--', label = 'Alert discovery', s=10, zorder=5)
    
    plt.grid(which='both', alpha=0.2, zorder=1)
    plt.yscale('log')
    plt.legend(loc=1, ncol=1, frameon=True, fontsize=14)
    #plt.loglog()
    plt.xlabel(r'$\sin \delta$')
    plt.ylabel(r'$E^2 \frac{dN_{\nu+\bar{\nu}}}{dEdA} \Big|_{\mathrm{1 TeV}}$ (GeV cm$^{-2}$)')
    #plt.ylim(2e-1, 4e2)
    plt.ylim(3e-2, 3e2)
    plt.text(0.03, 3e-1, '10 yr. time-integrated', color = sns.xkcd_rgb['light grey'])
    plt.title('9.6 years, time-integrated followup')
    #plt.xlim(5e0, 1.5e7)
except Exception as e:
    print(e)
    pass
plt.savefig(f_path + '../figures/followup_plots/steady/alert_event_sens_vs_dec_steady_smeared.png', 
               bbox_inches='tight', dpi=200)
