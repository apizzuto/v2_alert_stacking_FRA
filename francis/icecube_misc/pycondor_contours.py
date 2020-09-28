import pycondor, argparse, sys, os.path
from glob import glob
import numpy as np
import pandas as pd

error = '/scratch/apizzuto/fast_response/condor/error'
output = '/scratch/apizzuto/fast_response/condor/output'
log = '/scratch/apizzuto/fast_response/condor/log'
submit = '/scratch/apizzuto/fast_response/condor/submit'

job = pycondor.Job('Contours','./contours_from_maps.py',
			error=error,
			output=output,
			log=log,
			submit=submit,
            getenv=True,
            universe='vanilla',
			verbose=2, 
			request_memory=4000,
            #request_cpus=10,
			extra_lines= ['should_transfer_files = YES', 'when_to_transfer_output = ON_EXIT', 'Requirements =  (Machine != "node128.icecube.wisc.edu")']
			)

for ind in range(110):
    for nside_pow in [5,6,7,8,9,10]:
        job.add_arg('--ind={} --nsp={}'.format(ind, nside_pow))

dagman = pycondor.Dagman('Contours_from_skymaps', submit=submit, verbose=2)
dagman.add_job(job)
dagman.build_submit()
