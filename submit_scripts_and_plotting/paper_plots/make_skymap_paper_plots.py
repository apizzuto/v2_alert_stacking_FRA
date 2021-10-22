'''Script to make all plots for the fast response
alert followup paper. There is a lot of file I/O 
here, so don't be surprised that many of the plots
take a few minutes to make'''

import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from francis import utils
from francis.universe.transient_universe import SteadyUniverse
from results_utils import *

mpl.style.use(utils.initialize_mpl_style())

base_output = os.path.abspath(utils.get_francis_path() + '../figures/paper_plots/') + '/'
results_df = pd.read_csv(utils.get_francis_path() + 'icecube_misc/alert_dataframe_with_results.csv')

def all_sky_maps():
    '''Make all sky map with all of the contours
    First, we make a version with only the contours that 
    we used in the followup analysis. Next, we make a second
    copy with all alerts, but making those we don't look at gray'''
    print('Making all sky map with alert contours . . . ')
    for with_all_skymaps in [True, False]:
        all_sky_contour_map(
            colorscale='signalness', 
            with_problem_inds=with_all_skymaps)
        problem_ind_str = '_with_gray_extras' if with_all_skymaps else ''
        for ftype in ['.png', '.pdf']:
            plt.savefig(
                base_output + f'all_sky_contours{problem_ind_str}{ftype}',
                bbox_inches='tight'
            )
        plt.close()
    print('\t- Done with all sky map with contours')

def zoomed_in_maps():
    '''Now, make the zoomed in maps of the most significant alert followups'''
    print("Making zoomed in TS maps of the most significant followups")
    binom, pvals, tmps, inds = get_true_binomial_p(
        delta_t=None, 
        pvals=results_df['p_steady'], 
        with_inds=True
        )
    # These are just parameters for how zoomed in to make the plot
    resos = [4, 5, 3]
    for ind, reso in zip(inds[:3], resos):
        unique_id = get_unique_id(ind).replace('.', '_')
        print(f'\t- Making plot for alert {unique_id}')
        plot_results_zoom(
            ind, 
            with_text_box=True,
            cmap_range=[0., 15.],
            reso=reso
            )
        for ftype in ['.png', '.pdf']:
            plt.savefig(
                base_output + f'zoom_results_skymap_{unique_id}{ftype}',
                bbox_inches='tight'
            )
        plt.close()
    print('\t- Done with all zoomed in TS maps')

def schematic_paper_plot():
    '''Make a plot showing how we do the LLH analysis'''
    mpl.rcParams['font.size'] = 30.
    print('\nMaking plot that outlines how we do the LLH analysis')

    plot_results_zoom(20, no_prior=False, cmap_range=[0.0, 10.0], crosshair=True, title='Scan results',
                    axis_labels=False)
    plt.savefig('./tmp_2.png', bbox_inches='tight', dpi=400)
    plt.close()

    plot_results_zoom(20, no_prior=True, cmap_range=[0.0, 10.0], title='Point source information',
                    axis_labels=False)
    plt.savefig('./tmp_0.png', bbox_inches='tight', dpi=400)
    plt.close()

    plot_penalty_map(20, cmap = sns.color_palette("mako", as_cmap=True), cmap_range=[-10., 0.], 
                    title='Alert skymap information',
                    axis_labels=False)
    plt.savefig('./tmp_1.png', bbox_inches='tight', dpi=400)
    plt.close()

    fig, aaxs = plt.subplots(ncols=3, figsize=(17,7), dpi=400)
    plt.subplots_adjust(wspace=0.15)
    axs = aaxs.ravel()

    img = plt.imread('./tmp_0.png')
    axs[0].imshow(img)

    axs[0].text(
        axs[0].set_xlim()[1]*1.01, 
        axs[0].set_ylim()[0]*0.55, 
        r'$+$', 
        fontsize=60)

    img = plt.imread('./tmp_1.png')
    axs[1].imshow(img)

    axs[1].text(
        axs[1].set_xlim()[1]*1.01, 
        axs[1].set_ylim()[0]*0.55, 
        r'$=$', 
        fontsize=60)

    img = plt.imread('./tmp_2.png')
    axs[2].imshow(img)

    for ax in axs:
        ax.axis('off')
        
    for ftype in ['.png', '.pdf']:
            plt.savefig(
                base_output + f'schematic_skymap_figure{ftype}',
                bbox_inches='tight'
            )
    os.remove('./tmp_0.png')
    os.remove('./tmp_1.png')
    os.remove('./tmp_2.png')
    print('\t- Done with schematic figure')

if __name__ == "__main__":
    all_sky_maps()
    zoomed_in_maps()
    schematic_paper_plot()

    


