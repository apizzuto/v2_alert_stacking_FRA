import pycondor, argparse, sys, os
from glob import glob
import numpy as np
import pandas as pd

error = '/scratch/apizzuto/fast_response/condor/error'
output = '/scratch/apizzuto/fast_response/condor/output'
log = '/scratch/apizzuto/fast_response/condor/log'
submit = '/scratch/apizzuto/fast_response/condor/submit'

job = pycondor.Job('sensitivity_fastresponse_alerts','./alert_event_fits.py',
			error=error,
			output=output,
			log=log,
			submit=submit,
            getenv=True,
            universe='vanilla',
			verbose=2, 
			request_memory=8000,
            request_cpus=5,
			extra_lines= ['should_transfer_files = YES', 'when_to_transfer_output = ON_EXIT', 'Requirements =  (Machine != "node128.icecube.wisc.edu")']
			)

skymap_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run1*.fits.gz')
resubmit_inds = [6, 42, 58, 113, 161, 172, 190, 209, 229]

for index in resubmit_inds: #range(len(skymap_files)):
    for deltaT in np.array([1000.0, 2.*86400.]):
        for smear_str in [' --smear', '']:
            add_str = ' --ntrials=100' if deltaT > 100000. else ''
            job.add_arg('--deltaT={} --index={}'.format(deltaT, index) + add_str + smear_str)

dagman = pycondor.Dagman('alert_event_fra_fits', submit=submit, verbose=2)
dagman.add_job(job)
dagman.build_submit()
