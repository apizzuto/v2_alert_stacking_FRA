#!/usr/bin/env python

import numpy as np
import healpy as hp
import os, sys, argparse, pickle
from astropy.time import Time
from astropy.time import TimeDelta
from numpy.lib.recfunctions import append_fields
from astropy.io import fits
from glob import glob

from fast_response.FastResponseAnalysis import FastResponseAnalysis

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--index', type=int,default=None,
                    help='skymap index')
parser.add_argument('--deltaT', type=float, default=None,
                    help='Time Window in seconds')
parser.add_argument('--no-smear', default=False, action='store_true',
                    help='Do not include systematics, instead directly convert LLH to norm. prob.')
args = parser.parse_args()

# Commented path is the original trials location
# output_paths = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/results/'
output_base = os.path.join(os.path.expandvars("$PWD"), "analysis_trials/")
if not os.path.exists(output_base):
    os.mkdir(output_base)
output_paths = output_base + 'results/'
if not os.path.exists(output_paths):
    os.mkdir(output_paths)

skymap_files = sorted(glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run1*.fits.gz'))
skymap_fits, skymap_header = hp.read_map(skymap_files[args.index], h=True, verbose=False)

skymap_header = {name: val for name, val in skymap_header}
ev_mjd = skymap_header['EVENTMJD']
ev_run, ev_id = skymap_header['RUNID'], skymap_header['EVENTID']
source = {"Skipped Events": [(ev_run, ev_id)]}
source['Name'] = "RUN {} EVENT {} time window {:.2e}".format(str(skymap_header['RUNID']), str(skymap_header['EVENTID']), args.deltaT)
source['alert_type'] = 'track'

deltaT = args.deltaT / 86400.
event_mjd = ev_mjd
start_mjd = event_mjd - (deltaT / 2.)
stop_mjd = event_mjd + (deltaT / 2.)
start = Time(start_mjd, format='mjd').iso
stop = Time(stop_mjd, format='mjd').iso

f = FastResponseAnalysis(skymap_files[args.index], start, stop, save=True, 
                            alert_event=True, smear=not args.no_smear, **source)
inj = f.initialize_injector(gamma=2.5) #just put this here to initialize f.spatial_prior
ts = f.unblind_TS()
smear_str = 'smeared/' if not args.no_smear else 'norm_prob/'
if not os.path.exists(output_paths + smear_str):
    os.mkdir(output_paths + smear_str)
#res = f.save_results(alt_path = output_paths + smear_str)
f.plot_ontime()
f.calc_pvalue()
f.make_dNdE()
f.plot_tsd()
results = f.save_results()
f.generate_report()
