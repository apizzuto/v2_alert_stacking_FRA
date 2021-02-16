#!/usr/bin/env python

'''
Script to generate ntrials number of background only test statistics and save 
information about each trial to outfile

@options: 
    --i index: alert event index
    --ntrials: number of trials to perform
'''
import time, pickle, argparse
import numpy as np 
import pandas as pd

from francis.time_integrated_scripts.config_steady import config

def unblind_steady_map(
    index, seed, smear=True, local_skymap=False,
    verbose=True
    ):
    """DOCSTRING"""

    smear_str = 'smeared/' if smear else 'norm_prob/'
    alert_df = pd.DataFrame.from_csv('/data/user/apizzuto/fast_response_skylab/alert_event_followup/FRANCIS/francis/icecube_misc/alert_dataframe.csv')
    event_id = alert_df.iloc[index]['Event ID']
    run_id = alert_df.iloc[index]['Run ID']
    if not local_skymap:
        outfile = '/data/user/apizzuto/fast_response_skylab/' \
            + 'alert_event_followup/analysis_trials/results/' \
            + '{}index_{}_run_{}_event_{}_steady_seed_{}.pkl'.format(smear_str, index, run_id, event_id, seed)
    else:
        outfile = '/data/user/apizzuto/fast_response_skylab/' \
        + 'alert_event_followup/analysis_trials/results/' \
        + '{}index_{}_run_{}_event_{}_steady_ts_map_seed_{}.pkl'.format(
            smear_str, index, run_id, event_id, seed
            )

    t0 = time.time()
    nside = 2**7
    multillh, spatial_prior = config(
        index, gamma = 2.0, seed = seed, scramble = False, nside=nside, 
        ncpu = 1, injector = False, verbose=verbose, smear=smear,
        remove = True
        )

    t1 = time.time()
    if verbose:
        print('{:.2f} seconds to Initialize Likelihoods'.format(t1 - t0))
        print ("\nRunning fit on real data ...")

    allspots = None
    ii = 1
    n_iter = 2 if not local_skymap else 1
    for results, hotspots in multillh.do_allsky_trials(
        n_iter= n_iter, injector=None, nside=nside, rng_seed = 123*seed + ii,
        spatial_prior=spatial_prior, follow_up_factor = 1,
        scramble = False
        ):
        if verbose:
            print('Trial Number: {}'.format(ii))
        ii += 1
        if not local_skymap:
            if allspots is None:
                allspots = {}
                for k, v in hotspots['spatial_prior_0']['best'].items():
                    allspots[k] = [v]
                if 'pix' not in allspots.keys():
                    allspots['pix'] = [0]
                if 'nside' not in allspots.keys():
                    allspots['nside'] = [0]
            else:
                for k, v in hotspots['spatial_prior_0']['best'].items():
                    allspots[k].append(v)
                if 'pix' not in hotspots['spatial_prior_0']['best'].keys():
                    allspots['pix'].append(0)
                if 'nside' not in hotspots['spatial_prior_0']['best'].keys():
                    allspots['nside'].append(0)
        else:
            allspots = results

    if local_skymap:
        allspots = allspots[allspots['TS'] != 0.]

    dt1 = t1 - t0
    dt = time.time() - t0
    if verbose:
        print("Finished script in {} seconds".format(dt))
        print("Initialization: {} seconds\ntrials: {} seconds".format(
            dt1, (dt-dt1)
            ))

    with open(outfile, 'w') as f:
        pickle.dump(allspots, f, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Alert event followup steady unblinding'
        )
    parser.add_argument(
        '--i', type=int, required=True, help='Alert event index'
        )
    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help="Assorted print statements flag"
        )
    parser.add_argument(
        '--rng', type=int, default=1, help="Random number seed"
        )
    parser.add_argument(
        '--smear', default=False, action='store_true',
        help='Include systematics by smearing norm. prob.'
        )
    parser.add_argument(
        '--local_skymap', default=False, action='store_true',
        help='Instead of returning best fit, return whole map'
        )
    args = parser.parse_args()

    unblind_steady_map(
        args.i, args.rng, smear=args.smear, local_skymap=args.local_skymap,
        verbose=args.verbose
        )
