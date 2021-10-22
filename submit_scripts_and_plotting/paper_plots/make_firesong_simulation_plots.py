import os, copy, pickle
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm
from francis.universe.transient_universe import SteadyUniverse
from francis.universe.transient_universe import *
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from francis import utils
import matplotlib as mpl

base_output = os.path.abspath(utils.get_francis_path() + '../figures/paper_plots/') + '/'
mpl.style.use(utils.initialize_mpl_style())
paired = sns.color_palette('Paired')

def firesong_plots():
    print("Simulating populations of sources")
    pop_params, uni_1, all_sources_1, uni_2, all_sources_2 = get_universe_simulations()
    print("Beginning two dimensional plot")
    two_dimensional_plot(uni_1, all_sources_1, uni_2, all_sources_2)
    for ftype in ['.png', '.pdf']:
        plt.savefig(
            base_output + f'analysis_motivation_firesong_schematic{ftype}',
            bbox_inches='tight'
        )
    plt.close()
    print('\t- Done with two dimensional firesong plot')
    print("Beginning one dimensional projections")
    one_dimensional_projections(uni_1, all_sources_1, uni_2, all_sources_2)
    for ftype in ['.png', '.pdf']:
        plt.savefig(
            base_output + f'firesong_simulation_flux_distribution{ftype}',
            bbox_inches='tight'
        )
    plt.close()
    print("\t- Done with one dimensional projections")

def get_universe_simulations():
    pop_params = {'dens': [1e-8, 1e-6], 'sigma': [0.3, 0.15]}
    uni_1 = SteadyUniverse(
        'LG', 'MD2014SFR', pop_params['dens'][0], 1.5e-8, 2.50,
        sigma=pop_params['sigma'][0], data_years=9.6, seed=5
        )
    uni_1.create_universe()
    all_sources_1 = copy.deepcopy(uni_1.sources)
    uni_1.find_alerts()

    uni_2 = SteadyUniverse(
        'LG', 'MD2014SFR', pop_params['dens'][1], 1.5e-8, 2.50,
        data_years=9.6, sigma=pop_params['sigma'][1], seed=5
        )
    uni_2.create_universe()
    all_sources_2 = copy.deepcopy(uni_2.sources)
    uni_2.find_alerts()
    return pop_params, uni_1, all_sources_1, uni_2, all_sources_2

def two_dimensional_plot(uni_1, all_sources_1, uni_2, all_sources_2):
    gold_norm = uni_1.n_per_dec['HESE_gold'][1] + uni_1.n_per_dec['GFU_gold'][1]
    bronze_norm = uni_1.n_per_dec['HESE_bronze'][1] + uni_1.n_per_dec['GFU_bronze'][1]
    total_norm = gold_norm + bronze_norm

    f_path = utils.get_francis_path()
    eff_area_path = f_path + 'icecube_misc/effective_areas_alerts/'

    src_dec = 0.0

    with open(eff_area_path + 'gfu_online_effective_area_spline.npy', 'rb') as f:
        gfu_eff_spline = pickle.load(f, encoding='latin1')
    en_bins = np.logspace(2., 9., 501)
    ens = centers(en_bins)
    mu_extra = np.sum(gfu_eff_spline(np.log10(ens), np.sin(src_dec*np.pi / 180)) * \
                spectrum(ens, gamma = -2.5, flux_norm=1.0)*np.diff(en_bins)*1e4)
    extra_evs_norm = 1./mu_extra

    z_bins = np.linspace(0., 5.5, 100)
    fl_bins = np.logspace(-17., -10., 100)

    fig = plt.figure(dpi=200, figsize=(7,6))
    fig.set_facecolor('w')
    plt.subplots_adjust(wspace=0.07)

    gs = GridSpec(4,4)
    ax_joint = fig.add_subplot(gs[1:4,0:3])
    ax_marg_y = fig.add_subplot(gs[1:4,3])

    gold_norm = uni_1.n_per_dec['HESE_gold'][1] + uni_1.n_per_dec['GFU_gold'][1]
    bronze_norm = uni_1.n_per_dec['HESE_bronze'][1] + uni_1.n_per_dec['GFU_bronze'][1]
    total_norm = gold_norm + bronze_norm
    paired = sns.color_palette('Paired')
    alphas = [0.1, 0.1]
    custom_lines = [
        Line2D([0], [0], lw=0, marker='s', color=paired[0], 
            label=r'$\rho = 10^{-8}$ Mpc$^{-3}$' \
            + '\n' + r'$\mathrm{LF}=\mathrm{LN}(\sigma=0.3)$'),
        Line2D([0], [0], lw=0, marker='s', color=paired[2], 
            label=r'$\rho = 10^{-6}$ Mpc$^{-3}$' \
            + '\n' + r'$\mathrm{LF}=\mathrm{LN}(\sigma=0.15)$'),
        Line2D([0], [0], lw=0, marker='*', color='k', markersize=8,
            label='Alert')
        ]

    delta_t = 9.6 * 86400. * 365. # fluxes are currently dN/dE (at 1 GeV) * DeltaT for 9.6 years of data
    gev_per_tev = 1000. # Quote flux in per TeV
    norm_at_1_tev = 10.**(-2.5 * 3) # Quote flux at 1 TeV instead of 1 GeV (gamma=2.5)
    fl_conversion = norm_at_1_tev * gev_per_tev / delta_t

    lower_lim = 1e-4

    blues = sns.color_palette("Blues")
    greens = sns.color_palette("Greens")
    whites = [(1., 1., 1.), (1., 1., 1.)]
    both_colors = list(reversed(greens)) + list(whites) + list(blues)
    cmap = LinearSegmentedColormap.from_list('greenblue', both_colors, N=100)

    weights_1 = np.ones_like(all_sources_1['z'])
    weights_2 = -1.*np.ones_like(all_sources_2['z'])
    all_z = np.append(all_sources_1['z'], all_sources_2['z'])
    all_fl = np.append(all_sources_1['flux'], all_sources_2['flux'])
    all_weights = np.append(weights_1, weights_2)

    max_counts_per_bin = 1e5

    ax_joint.hist2d(all_z, all_fl * fl_conversion, 
            bins=[z_bins, fl_bins],
            weights=all_weights, cmap=cmap, # vmin=-max_counts_per_bin, vmax=max_counts_per_bin,
            norm=SymLogNorm(linthresh=1.0, linscale=1.0,
                            vmin=-max_counts_per_bin, vmax=max_counts_per_bin, base=10))

    cbax_1 = fig.add_axes([0.13, 0.7, 0.27, 0.03]) 
    cmap_1 = LinearSegmentedColormap.from_list(
        'myblue',
        list(whites)[:1] + list(blues),
        N=50)

    cbar_1 = mpl.colorbar.ColorbarBase(
        cbax_1, 
        cmap=cmap_1,
        norm=SymLogNorm(
            linthresh=1.0, 
            linscale=1.0,
            vmin=0., 
            vmax=max_counts_per_bin,
            base=10),
        orientation='horizontal',
        ticklocation='top',
        ticks=[1e0, 1e2, 1e4],
        )
    cbar_1.ax.tick_params(direction='out', labelsize=12)
    cbar_1.set_label(r'$N_{\mathrm{sources}}$', size=14)

    cbax_2 = fig.add_axes([0.43, 0.7, 0.27, 0.03]) 
    cmap_2 = LinearSegmentedColormap.from_list(
        'mygreen',
        list(whites)[:1] + list(greens),
        N=50)

    cbar_2 = mpl.colorbar.ColorbarBase(
        cbax_2, 
        cmap=cmap_2,
        norm=SymLogNorm(
            linthresh=1.0, 
            linscale=1.0,
            vmin=0., 
            vmax=max_counts_per_bin,
            base=10),
        orientation='horizontal',
        ticklocation='top',
        ticks=[1e0, 1e2, 1e4],
        )
    cbar_2.ax.tick_params(direction='out', labelsize=12)
    cbar_2.set_label(r'$N_{\mathrm{sources}}$', size=14)

    alert_colors = [sns.xkcd_rgb['navy blue'], sns.xkcd_rgb['navy green']]

    for i, (uni, srcs) in enumerate(zip([uni_1, uni_2], [all_sources_1, all_sources_2])):
        plottable = srcs['flux'] > lower_lim 
        
        dec_inds = [np.argwhere(srcs['dec'] == dec)[0,0] for dec in uni.sources['dec']]
        fl_inds = [np.argwhere(srcs['flux'] == fl)[0,0] for fl in uni.sources['flux']]
        assert dec_inds == fl_inds
        alert_inds = np.asarray(dec_inds)
        alert_msk = np.zeros_like(srcs['dec'], dtype=bool)
        alert_msk[alert_inds] = True
        
        ax_joint.plot(srcs['z'][alert_msk * plottable],
                srcs['flux'][alert_msk * plottable] * fl_conversion,
                linestyle='', marker='*', markersize=5,
                color=alert_colors[i],
                alpha=1.0, zorder=10,
                markeredgewidth=0.3
                )
        
        h, b = np.histogram(
            srcs['flux'][plottable] * fl_conversion, 
            bins=np.logspace(-18, -10, 50)
            )
        cdf = np.cumsum(np.flip(h))
        bin_centers = 10.**(np.log10(b[:-1]) + np.diff(np.log10(b))/2.)
        ax_marg_y.plot(cdf, np.flip(bin_centers), lw=2.5, color=paired[2*i])
        
    ax_joint.axhline(1./total_norm * fl_conversion, ls='dashed', color='grey')
    ax_joint.text(0.8, 0.35/total_norm * fl_conversion, 
            r'$\langle\mathcal{N}^{\mathrm{alert}}_{\delta=0^{\circ}}\rangle$ = 1',
            color='grey', fontsize=14)
    ax_joint.axhline(1./extra_evs_norm * fl_conversion, ls='dotted', color='grey')
    ax_joint.text(2.8, 0.35/extra_evs_norm * fl_conversion, 
            r'$\langle\mathcal{N}^{\mathrm{GFU}}_{\delta=0^{\circ}}\rangle$ = 1',
            color='grey', fontsize=14)
    ax_joint.set_yscale('log')
    ax_joint.set_ylim(lower_lim * fl_conversion, ax_joint.set_ylim()[1])
    ax_joint.set_xlim(-0.08, 5.1)
    ax_joint.set_xlabel(r'$z$')
    ax_joint.set_ylabel(r'$\frac{dN}{dE}\big|_{1 \mathrm{TeV}}$ (TeV$^{-1}$ cm$^{-2}$ s$^{-1}$)')
    ll = ax_joint.legend(
        handles=custom_lines, 
        loc=1, 
        fontsize=12, 
        framealpha=0.95)
    ll.set_zorder(11)
    ax_marg_y.set_yscale('log')
    ax_marg_y.set_xscale('log')
    ax_marg_y.set_ylim(ax_joint.set_ylim())
    ax_marg_y.set_xlabel(
        r'$N_{\mathrm{sources}} > \frac{dN}{dE}$', 
        fontsize=12
        )
    ax_marg_y.set_xticks([1e0, 1e3, 1e6])
    # Turn off tick labels on marginals
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)


def one_dimensional_projections(uni_1, all_sources_1, uni_2, all_sources_2):
    fig, aaxs = plt.subplots(ncols=2, nrows=1, figsize=(8,3), dpi=200, sharey=True)
    fig.subplots_adjust(wspace=0.04)

    delta_t = 9.6 * 86400. * 365. # fluxes are currently dN/dE (at 1 GeV) * DeltaT for 9.6 years of data
    gev_per_tev = 1000. # Quote flux in per TeV
    norm_at_1_tev = 10.**(-2.5 * 3) # Quote flux at 1 TeV instead of 1 GeV (gamma=2.5)
    fl_conversion = norm_at_1_tev * gev_per_tev / delta_t

    axs = aaxs.ravel()
    z_bins = np.linspace(0., 10., 26)
    fl_bins = np.logspace(-18., -9., 26)

    custom_lines = [
        Line2D([0], [0], lw=2., color=paired[0], 
            label='Population A'),
        Line2D([0], [0], lw=2., color=paired[1], 
            label='Population A, Alerts'),
        Line2D([0], [0], lw=2, color=paired[2], 
            label='Population B'),
        Line2D([0], [0], lw=2, color=paired[3], 
            label='Population B, Alerts')
        ]

    for i, (uni, srcs) in enumerate(zip([uni_1, uni_2], [all_sources_1, all_sources_2])):
        dec_inds = [np.argwhere(srcs['dec'] == dec)[0,0] for dec in uni.sources['dec']]
        alert_inds = np.asarray(dec_inds)
        alert_msk = np.zeros_like(srcs['dec'], dtype=bool)
        alert_msk[alert_inds] = True
        
        axs[0].hist(
            srcs['z'], 
            bins=z_bins, 
            histtype='step', 
            lw=2., 
            color=paired[2*i]
        )
        axs[1].hist(
            srcs['flux'] * fl_conversion, 
            histtype='step',
            bins=fl_bins, 
            lw=2., 
            color=paired[2*i]
        )
        
        axs[0].hist(
            srcs['z'][alert_msk],
            bins=z_bins, 
            histtype='step', 
            lw=2., 
            color=paired[2*i + 1]
        )
        axs[1].hist(
            srcs['flux'][alert_msk] * fl_conversion, 
            histtype='step',
            bins=fl_bins, 
            lw=2., 
            color=paired[2*i + 1])
        
        
    axs[0].set_xlim(axs[0].set_xlim()[0], 10.)    
    axs[0].set_xlabel(r'$z$')
    axs[0].set_ylabel(r'$N_{\mathrm{sources}}$ per bin')
    axs[0].set_yscale('log')
    axs[1].set_xlabel(r'$\frac{dN}{dE}\big|_{1 \mathrm{TeV}}$ (TeV$^{-1}$ cm$^{-2}$ s$^{-1}$)')
    axs[1].set_xscale('log')
    axs[1].set_xlim(3e-18, axs[1].set_xlim()[1])
    axs[0].set_ylim(axs[0].set_ylim()[0], 1e7)
    axs[1].legend(
        handles=custom_lines, 
        loc=1, 
        fontsize=10, 
        frameon=True)

if __name__ == '__main__':
    firesong_plots()