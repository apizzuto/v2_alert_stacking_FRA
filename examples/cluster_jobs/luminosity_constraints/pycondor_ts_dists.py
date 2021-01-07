import pycondor, argparse, sys, os
from glob import glob
import numpy as np
import pandas as pd

error = '/scratch/apizzuto/fast_response/condor/error'
output = '/scratch/apizzuto/fast_response/condor/output'
log = '/scratch/apizzuto/fast_response/condor/log'
submit = '/scratch/apizzuto/fast_response/condor/submit'

job = pycondor.Job('calc_ts_dists','/data/user/apizzuto/fast_response_skylab/fast-response/trunk/time_integrated_scripts/steady_2d_ts_sampling.py',
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
for density in np.logspace(-11., -6., 21)[:]:
    for lumi in ['SC']: #, 'LG']:
        for evol in ['MD2014SFR']:  #['HB2006SFR', 'MD2014SFR']:
            N = 500 if density < 1e-7 else 40
            job.add_arg('--density={} --LF={} --n={} --evol={}'.format(density, lumi, N, evol))

dagman = pycondor.Dagman('ts_distributions', submit=submit, verbose=2)
dagman.add_job(job)
dagman.build_submit()
