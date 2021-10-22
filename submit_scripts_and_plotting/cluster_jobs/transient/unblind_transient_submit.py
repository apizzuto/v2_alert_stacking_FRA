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

job = pycondor.Job('sensitivity_fastresponse_alerts',
	f_path + 'transient_scripts/alert_event_unblind.py',
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

# Could also pass ' --no-smear' as the smear_str. This converts LLh to a 
# normalized probability. We leave the option as it was tested
# and compared against our smearing procedure, but not used in the 
# final version of the analysis
smear_str = ''

for index in range(244, len(skymap_files)):
	# Time windows in seconds
    for deltaT in np.array([1000.0, 2.*86400.]):
        job.add_arg('--deltaT={} --index={}'.format(deltaT, index) + smear_str)

dagman = pycondor.Dagman('alert_event_fra_unblind', submit=submit, verbose=2)
dagman.add_job(job)
dagman.build_submit()
