import pycondor, argparse, sys, os
from glob import glob
import numpy as np
import pandas as pd

error = '/scratch/apizzuto/fast_response/condor/error'
output = '/scratch/apizzuto/fast_response/condor/output'
log = '/scratch/apizzuto/fast_response/condor/log'
submit = '/scratch/apizzuto/fast_response/condor/submit'

job = pycondor.Job('signal_alerts_steady','/data/user/apizzuto/fast_response_skylab/alert_event_followup/FRANCIS/francis/time_integrated_scripts/signal_trials.py',
			error=error,
			output=output,
			log=log,
			submit=submit,
            getenv=True,
            universe='vanilla',
			verbose=1, 
			request_memory=8000,
            #request_cpus=5,
			extra_lines= ['should_transfer_files = YES', 'when_to_transfer_output = ON_EXIT']
			)

sky_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run1*.fits.gz')

for rng in range(10):
    for index in range(245, len(sky_files)): #ONLY DO A FEW BEFORE WE GET THE FULL SAMPLE
        for gamma in [2.5]:
            for fit in [True, False]:
                for smear in [' --smear']: #, '']:
                    add_str = ' --fit' if fit else ''
                    job.add_arg('--rng={} --i={} --ntrials=10 --g={}{}'.format(rng, index, gamma, smear) + add_str)

dagman = pycondor.Dagman('alert_event_fra_signal_steady', submit=submit, verbose=1)
dagman.add_job(job)
dagman.build_submit()
