import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')

from francis.universe.universe_plotter import UniversePlotter
from francis import utils
utils.initialize_mpl_style()
import matplotlib as mpl

def make_analysis_plots(uni):
    """
    Call all relevant class methods to create
    plots needed for a fixed time window analysis

    :type uni: UniversePlotter instance
    :param uni: analysis for which we are making plots
    """
    print("\tMaking plots now . . . ")
    try:
        compare = True if uni.delta_t == 1e3 else False
        uni.two_dim_sensitivity_plot_ts(compare=compare, in_ts=False)
        uni.rotated_sensitivity_plot_ts(in_ts=False, compare=compare)
        uni.plot_background_binomial_p()
        uni.brazil_bands(compare=compare)
        uni.brazil_bands(rotated=True, compare=compare)
        uni.one_dim_ts_distributions(in_ts = False)
    except:
        print("\t\tUNABLE TO MAKE ONE OR MORE PLOTS")
        print("\t\t(likely one dim transient dist)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plot creation for stacked alert followup analysis'
        )
    parser.add_argument(
        '--lumi', type=str, default='SC', 
        help='Luminosity function, default standard candle'
        )
    parser.add_argument(
        '--evol', type=str, default='MD2014SFR', 
        help="Evolution, default Madau Dickinson"
        )
    parser.add_argument(
        '--data_years', type=float, default=8.6,
        help="Years of alert events (default=8.6)"
        )
    args = parser.parse_args()

    for time_window in ['steady', 2.*86400., 1000.]:
        if type(time_window) == str:
            analysis_type = 'Time integrated'
        else:
            analysis_type = 'Transient'
        print("\nInitializing plotter for {} analysis:".format(analysis_type))
        if analysis_type == 'Transient':
            print('\tTime window: {:.2e} seconds'.format(time_window))
        delta_t = None if time_window == 'steady' else time_window
        fig_subdir = 'steady/' if time_window == 'steady' else 'transient/'
        savepath = '/data/user/apizzuto/fast_response_skylab/' + \
            'alert_event_followup/FRANCIS/figures/' + fig_subdir
        uni = UniversePlotter(
            delta_t, args.data_years, args.lumi, args.evol, save=True, 
            savepath=savepath,
            verbose=True
            )
        make_analysis_plots(uni)
    print("\nFinished making all plots!!!")
    print("Check output in {}".format(savepath))