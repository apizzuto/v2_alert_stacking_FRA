import pycondor, argparse, sys, os, pwd
from glob import glob
import numpy as np
import pandas as pd

import francis.utils as utils
f_path = utils.get_francis_path()

username = pwd.getpwuid(os.getuid())[0]
if not os.path.exists(f'/scratch/{username}/'):
    os.mkdir(f'/scratch/{username}/')
if not os.path.exists(f'/scratch/{username}/fast_response/'):
    os.mkdir(f'/scratch/{username}/fast_response/')
if not os.path.exists(f'/scratch/{username}/fast_response/condor/'):
    os.mkdir(f'/scratch/{username}/fast_response/condor')

error = f'/scratch/{username}/fast_response/condor/error'
output = f'/scratch/{username}/fast_response/condor/output'
log = f'/scratch/{username}/fast_response/condor/log'
submit = f'/scratch/{username}/fast_response/condor/submit'

low_mem_job = pycondor.Job('low_mem_calc_ts_dists',
    f_path + 'universe/transient_full_2d_ts_sampling.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2, 
    request_memory=4000,
    request_cpus=1,
    extra_lines= ['should_transfer_files = YES', 
        'when_to_transfer_output = ON_EXIT', 
        'Requirements =  (Machine != "node128.icecube.wisc.edu")']
    )

high_mem_job = pycondor.Job('high_mem_calc_ts_dists',
    f_path + 'universe/transient_full_2d_ts_sampling.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2, 
    request_memory=8000,
    request_cpus=1,
    extra_lines= ['should_transfer_files = YES', 
        'when_to_transfer_output = ON_EXIT', 
        'Requirements =  (Machine != "node128.icecube.wisc.edu")']
    )

base_dir = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/ts_distributions/'
for delta_t in [1000., 2.*86400.]:
    # Luminosities in units of erg/yr, assuming constant (non-transient) emission
    for manual_lumi in np.logspace(49, 60, 45):
        # Rate density in Mpc^-3 yr^-1
        for density in np.logspace(-11., -6., 21)[:]:
            for lumi, sigma in [('SC', 1.), ('LG', 0.4)]:
                for evol in ['MD2014SFR']:
                    if density * manual_lumi * delta_t > 8e54:
                        continue
                    elif density * manual_lumi * delta_t < 3e46:
                        continue
                    else:
                        if lumi == 'SC':  
                            fname = base_dir + f'ts_dists_9.6year_density_{density:.2e}_evol_{evol}_lumi_{lumi}_manual_lumi_{manual_lumi:.1e}_delta_t_{delta_t:.2e}.npy'
                        else: 
                            fname = base_dir + f'ts_dists_9.6year_density_{density:.2e}_evol_{evol}_lumi_{lumi}_sigma_{sigma:.2f}_manual_lumi_{manual_lumi:.1e}_delta_t_{delta_t:.2e}.npy'
                        if os.path.exists(fname):
                            continue
                        N = 300 if density < 1e-7 else 120
                        if density > 3e-8:
                            high_mem_job.add_arg('--delta_t={} --density={} --LF={} --n={} --evol={} --manual_lumi={} --sigma={}'.format(delta_t, density, lumi, N, evol, manual_lumi, sigma))
                        else:
                            low_mem_job.add_arg('--delta_t={} --density={} --LF={} --n={} --evol={} --manual_lumi={} --sigma={}'.format(delta_t, density, lumi, N, evol, manual_lumi, sigma))

dagman = pycondor.Dagman('ts_distributions', submit=submit, verbose=2)

dagman.add_job(low_mem_job)
dagman.add_job(high_mem_job)

dagman.build_submit()
