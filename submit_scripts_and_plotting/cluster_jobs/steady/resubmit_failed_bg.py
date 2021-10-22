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

job = pycondor.Job('background_fastresponse_alerts_steady',
    f_path + 'time_integrated_scripts/steady_background_trials.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2, 
    request_memory=8000,
    #request_cpus=5,
    extra_lines= ['+AccountingGroup = "1_week.apizzuto"', 
        'should_transfer_files = YES', 
        'when_to_transfer_output = ON_EXIT']
    )

sky_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run1*.fits.gz')

# Could also pass ' --no-smear' as the smear. This converts LLh to a 
# normalized probability. We leave the option as it was tested
# and compared against our smearing procedure, but not used in the 
# final version of the analysis
smear = ''

resubmit_key_pairs = [(15, 16), (32, 7), (37, 6), (43, 8), (43, 15), (65, 2), (72, 15), (78, 3), (100, 3), (133, 5), (133, 6), (133, 15), (162, 16), (172, 0), (172, 1), (172, 2), (172, 5), (172, 7), (172, 12), (172, 13), (172, 17), (172, 18), (172, 19), (174, 18), (195, 3), (203, 19), (213, 2), (213, 7), (213, 8), (213, 10), (213, 11), (222, 11), (224, 2), (240, 9), (245, 5), (258, 18)]

for index, rng in resubmit_key_pairs:
    job.add_arg('--rng={} --i={} --ntrials=200{}'.format(rng, index, smear))

dagman = pycondor.Dagman('alert_event_fra_background_steady', submit=submit, verbose=2)
dagman.add_job(job)
dagman.build_submit()
