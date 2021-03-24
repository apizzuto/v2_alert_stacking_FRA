#!/usr/bin/env python
import argparse
import time
import pickle
import numpy as np
import pandas as pd

from francis.time_integrated_scripts.config_steady import config

def run_signal_trials(
    gamma, index, ntrials, seed, verbose=True, 
    fit=True, smear=True
    ):
    """Run the signal trials for a time-integrated alert followup.
    
    Args:
        gamma (float): Injected spectral index
        index (int): Alert event index
        ntrials (int): number of trials
        seed (int): random number seed
        verbose (bool, default=True): print output
        fit (bool, default=True): if True, do not poisson fluctuate fluxes
        smear (bool, default=True): Account for sytematics in the skymap
            millipede localization
    """

    poisson=~fit
    smear_str = 'smeared/' if smear else 'norm_prob/'
    outdir = 'fits' if fit else 'sensitivity'
    alert_df = pd.read_csv('/data/user/apizzuto/fast_response_skylab/alert_event_followup/FRANCIS/francis/icecube_misc/alert_dataframe.csv')
    event_id = alert_df.iloc[index]['Event ID']
    run_id = alert_df.iloc[index]['Run ID']
    outfile = '/data/user/apizzuto/fast_response_skylab/' \
        + 'alert_event_followup/analysis_trials/' \
        + '{}/{}index_{}_run_{}_event_{}_steady_seed_{}_gamma_{}.pkl'.format(
            outdir, smear_str, index, run_id, event_id, seed, gamma
            )

    t0 = time.time()
    nside = 2**7
    multillh, spatial_prior, inj = config(
        index, gamma = gamma, seed = seed, scramble = True, nside=nside,
        ncpu = 1, injector = True, verbose=verbose, smear=smear
        )

    t1 = time.time()
    if verbose:
        print('{:.2f} seconds to Initialize Likelihoods'.format(t1 - t0))
        print ("\nRunning signal injection trials ...")

    allspots = None
    ii = 1
    ninj = np.append(
        np.linspace(1, 10, 10), 
        np.array([15, 20, 25, 30, 40, 50, 60, 70])
        )

    scale_arr = []
    if gamma > 2.0:
        step_size = 10
    else:
        step_size = 3
    for i in range(1,int(20*step_size)+1, step_size):
        scale_arr.append([])
        for j in range(5):
            scale_arr[-1].append(inj.sample(i, poisson=False)[0][0])
    scale_arr = np.median(scale_arr, axis=1)
    try:
        scale_factor = np.min(np.argwhere(scale_arr > 0))*step_size + 1.
    except:
        print("Scale factor thing for prior injector didn't work")
        scale_factor = 1.

    for ni in ninj:
        for results, hotspots in multillh.do_allsky_trials(
            n_iter= ntrials,
            injector=inj,
            mean_signal=ni*scale_factor,
            poisson=poisson, 
            nside=nside, rng_seed = 123*seed + ii,
            spatial_prior=spatial_prior,
            follow_up_factor = 1, return_position=True
            ):
            if verbose:
                print('Trial Number: {}'.format(ii))
            ii += 1
            if allspots is None:
                allspots = {}
                for k, v in hotspots['spatial_prior_0']['best'].items():
                    allspots[k] = [v]
                if 'pix' not in allspots.keys():
                    allspots['pix'] = [0]
                if 'nside' not in allspots.keys():
                    allspots['nside'] = [0]
                if 'inj' in hotspots['spatial_prior_0'].keys():
                    for k, v in hotspots['spatial_prior_0']['inj'].items():
                        allspots['inj_' + k] = [v]
                else:
                    allspots['inj_nsignal'] = [0]
                    allspots['inj_dec'], allspots['inj_ra'] = [0.], [0.]
                allspots['flux'] = [inj.mu2flux(ni*scale_factor)]
            else:
                for k, v in hotspots['spatial_prior_0']['best'].items():
                    allspots[k].append(v)
                if 'pix' not in hotspots['spatial_prior_0']['best'].keys():
                    allspots['pix'].append(0)
                if 'nside' not in hotspots['spatial_prior_0']['best'].keys():
                    allspots['nside'].append(0)
                if 'inj' in hotspots['spatial_prior_0'].keys():
                    for k, v in hotspots['spatial_prior_0']['inj'].items():
                        allspots['inj_' + k].append(v)
                else:
                    allspots['inj_nsignal'].append(0)
                    allspots['inj_dec'].append(0.0)
                    allspots['inj_ra'].append(0.0)
                allspots['flux'].append(inj.mu2flux(ni*scale_factor))

    dt1 = t1 - t0
    dt     = time.time() - t0
    if verbose:
        print("Finished trials in {} seconds".format(dt))
        print("Initialization: {} seconds\ntrials: {} seconds".format(
            dt1, (dt-dt1)
            ))

    with open(outfile, 'wb') as f:
        pickle.dump(allspots, f, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'time integrated alert followup tests'
        )
    parser.add_argument(
        '--i', type=int, required=True, help='Alert event index'
        )
    parser.add_argument(
        '--g', type=float, required=True, help='spectral index'
        )
    parser.add_argument(
        "--ntrials", default=10, type=int,
        help="Number of trials per signal strength (default=10)"
        )
    parser.add_argument(
        '--rng', type=int, default=1, help="Random number seed"
        )
    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help="Assorted print statements flag"
        )
    parser.add_argument(
        '--fit', action='store_true', default=False,
        help="Include poisson fluctuations by default or raise this flag"
        )                
    parser.add_argument(
        '--smear', default=False, action='store_true',
        help='Include systematics by smearing norm. prob., The unblinded version will use this flag'
        )
    args = parser.parse_args()

    run_signal_trials(
        args.g, args.i, args.ntrials, args.rng, verbose=args.verbose, 
        fit=args.fit, smear=args.smear
        )
