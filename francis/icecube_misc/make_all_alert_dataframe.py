r"""Construct a dataframe of all of the relevant
information for alert events from reading in the
details of the .fits files"""

import numpy as np
from glob import glob
import pandas as pd
from astropy.time import Time
import healpy as hp

alert_fs = sorted(glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run*.fits.gz'))
df_dict = {'Event ID': [], 'Run ID': [], 'Dec': [], 'RA': [], 'Signalness': [],
          'MJD': [], 'Stream': [], 'Cut': [], 'Dec err': [], 'RA err': []}

for i in range(len(alert_fs)):
    skymap_fits, skymap_header = hp.read_map(alert_fs[i], h=True, verbose=False)
    skymap_header = {name: val for name, val in skymap_header}
    df_dict['Event ID'].append(skymap_header['EVENTID'])
    df_dict['Run ID'].append(skymap_header['RUNID'])
    df_dict['Dec'].append(skymap_header['DEC'])
    df_dict['RA'].append(skymap_header['RA'])
    df_dict['Signalness'].append(skymap_header['SIGNAL'])
    df_dict['MJD'].append(skymap_header['EVENTMJD'])
    if 'gfu' in skymap_header['I3TYPE'].lower():
        stream = 'GFU'
    elif 'hese' in skymap_header['I3TYPE'].lower():
        stream = 'HESE'
    elif 'ehe' in skymap_header['I3TYPE'].lower():
        stream = 'EHE'
    else:
        stream = None
        print("What stream is this??")
    df_dict['Stream'].append(stream)
    if 'gold' in skymap_header['I3TYPE'].lower():
        cut = 'GOLD'
    elif 'bronze' in skymap_header['I3TYPE'].lower():
        cut = 'BRONZE'
    else:
        cut = None
        print("REALLY WTF")
    df_dict['Cut'].append(cut)
    dec_err = skymap_header['DEC_ERR_PLUS'] + skymap_header['DEC_ERR_MINUS']
    ra_err = skymap_header['RA_ERR_PLUS'] + skymap_header['RA_ERR_MINUS']
    df_dict['Dec err'].append(dec_err)
    df_dict['RA err'].append(ra_err)

alert_df = pd.DataFrame.from_dict(df_dict)
from francis import utils
f_path = utils.get_francis_path()
alert_df.to_csv(f_path + 'icecube_misc/alert_dataframe.csv')
