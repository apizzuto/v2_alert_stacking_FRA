import pycondor, argparse, sys, os
from glob import glob
import numpy as np
import pandas as pd
import francis.utils as utils
f_path = utils.get_francis_path()

error = '/scratch/apizzuto/fast_response/condor/error'
output = '/scratch/apizzuto/fast_response/condor/output'
log = '/scratch/apizzuto/fast_response/condor/log'
submit = '/scratch/apizzuto/fast_response/condor/submit'

job = pycondor.Job('background_fastresponse_alerts',
    f_path + 'transient_scripts/alert_event_background.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2, 
    request_memory=8000,
    request_cpus=5,
    extra_lines= ['should_transfer_files = YES', 
        'when_to_transfer_output = ON_EXIT', 
        'Requirements =  (Machine != "node128.icecube.wisc.edu")']
    )

skymap_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run1*.fits.gz')

for index in range(len(skymap_files)):
    for deltaT in np.array([1000., 2.*86400.]):
        for smear_str in [' --smear']:
            if deltaT > 100000.:
                ntr = 2500
            elif deltaT > 1000.:
                ntr = 5000
            else:
                ntr = 10000
            job.add_arg('--deltaT={} --index={} --ntrials={} {}'.format(deltaT, index, ntr, smear_str))

dagman = pycondor.Dagman('alert_event_fra_background', submit=submit, verbose=2)
dagman.add_job(job)
dagman.build_submit()
