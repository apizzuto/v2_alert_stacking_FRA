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

job = pycondor.Job('calc_ts_dists',
	f_path + 'universe/steady_2d_ts_sampling.py',
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

# Density in Mpc^-3
for density in np.logspace(-11., -6., 21)[:]:
    for lumi, sigma in [('SC', 1), ('LG', 0.4)]:
        for evol in ['MD2014SFR']: 
            N = 500 if density < 1e-7 else 100
            job.add_arg('--density={} --LF={} --n={} --evol={} --sigma={}'.format(density, lumi, N, evol, sigma))

dagman = pycondor.Dagman('ts_distributions', submit=submit, verbose=2)
dagman.add_job(job)
dagman.build_submit()
