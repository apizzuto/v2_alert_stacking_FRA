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

job = pycondor.Job('signal_alerts_steady',
    f_path + 'time_integrated_scripts/signal_trials.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=1, 
    request_memory=8000,
    #request_cpus=5,
    extra_lines= ['should_transfer_files = YES', 
        'when_to_transfer_output = ON_EXIT']
    )

sky_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run1*.fits.gz')

# Could also pass ' --no-smear' as the smear. This converts LLh to a 
# normalized probability. We leave the option as it was tested
# and compared against our smearing procedure, but not used in the 
# final version of the analysis
smear = ''

for rng in range(10):
    for index in range(len(sky_files)):
        # Could include other spectral indices. Fixed to 2.5 because
        # this is already a lot of jobs
        for gamma in [2.5]:
            for fit in [True, False]:
                add_str = ' --no-poisson' if fit else ''
                job.add_arg('--rng={} --i={} --ntrials=10 --g={}{}'.format(rng, index, gamma, smear) + add_str)

dagman = pycondor.Dagman('alert_event_fra_signal_steady', submit=submit, verbose=1)
dagman.add_job(job)
dagman.build_submit()
