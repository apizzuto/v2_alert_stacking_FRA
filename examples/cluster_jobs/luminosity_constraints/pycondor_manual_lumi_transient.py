import pycondor, argparse, sys, os
from glob import glob
import numpy as np
import pandas as pd

error = '/scratch/apizzuto/fast_response/condor/error'
output = '/scratch/apizzuto/fast_response/condor/output'
log = '/scratch/apizzuto/fast_response/condor/log'
submit = '/scratch/apizzuto/fast_response/condor/submit'

low_mem_job = pycondor.Job('low_mem_calc_ts_dists',
    '/data/user/apizzuto/fast_response_skylab/alert_event_followup/FRANCIS/francis/universe/transient_full_2d_ts_sampling.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2, 
    request_memory=4000,
    request_cpus=1,
    extra_lines= ['should_transfer_files = YES', 'when_to_transfer_output = ON_EXIT', 'Requirements =  (Machine != "node128.icecube.wisc.edu")']
    )

high_mem_job = pycondor.Job('high_mem_calc_ts_dists',
    '/data/user/apizzuto/fast_response_skylab/alert_event_followup/FRANCIS/francis/universe/transient_full_2d_ts_sampling.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2, 
    request_memory=8000,
    request_cpus=1,
    extra_lines= ['should_transfer_files = YES', 'when_to_transfer_output = ON_EXIT', 'Requirements =  (Machine != "node128.icecube.wisc.edu")']
    )

for delta_t in [1000., 2.*86400.]:
    for manual_lumi in np.logspace(49, 60, 45):
    #for manual_lumi in np.logspace(60, 64, 17):
        for density in np.logspace(-11., -6., 21)[:]:
            for lumi in ['SC']: # , 'LG']:
                for evol in ['MD2014SFR']:
                    if density * manual_lumi * delta_t > 3e54:
                        continue
                    elif density * manual_lumi * delta_t < 1e47:
                        continue
                    else:
                        N = 500 if density < 1e-7 else 40
                        if density > 3e-8:
                            high_mem_job.add_arg('--delta_t={} --density={} --LF={} --n={} --evol={} --manual_lumi={}'.format(delta_t, density, lumi, N, evol, manual_lumi))
                        else:
                            low_mem_job.add_arg('--delta_t={} --density={} --LF={} --n={} --evol={} --manual_lumi={}'.format(delta_t, density, lumi, N, evol, manual_lumi))

#for delta_t in [1000., 2.*86400.]:
#    for density in np.logspace(-11., -6., 21)[:]:
#        for lumi in ['SC', 'LG']:
#            for evol in ['MD2014SFR']:
#                N = 500 if density < 1e-7 else 40
#                job.add_arg('--delta_t={} --density={} --LF={} --n={} --evol={}'.format(delta_t, density, lumi, N, evol))                        

dagman = pycondor.Dagman('ts_distributions', submit=submit, verbose=2)

dagman.add_job(low_mem_job)
dagman.add_job(high_mem_job)

dagman.build_submit()
