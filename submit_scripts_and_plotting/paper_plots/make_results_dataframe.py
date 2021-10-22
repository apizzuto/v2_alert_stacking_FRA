import pandas as pd
import numpy as np
from francis import utils
from results_utils import *

base_output = os.path.abspath(utils.get_francis_path() + 'icecube_misc/') 
alert_df = pd.read_csv('/home/apizzuto/wg-nu-sources/2021_v2_alert_stacking_FRA/francis/icecube_misc/alert_dataframe.csv')

for delta_t in [1000., 172800.]:
    tss, nss, pvals = get_unblinded_vals(delta_t=delta_t)
    alert_df[f'TS_{delta_t:.1e}'] = tss
    alert_df[f'ns_{delta_t:.1e}'] = nss
    alert_df[f'p_{delta_t:.1e}'] = pvals

tss, nss, gammas, pvals, ras, decs = [], [], [], [], [], [] 
for ind in range(275):
    problem_inds = [73,  76, 142, 147, 157, 198, 249]
    if ind in problem_inds:
        for val_list in [tss, nss, gammas, pvals, ras, decs]:
            val_list.append(np.nan)
    else:
        try:
            al = Alert(ind, None)
            for val, val_list in [(al.ts, tss), (al.ns, nss), (al.gamma, gammas),
                                 (al.pval, pvals), (al.ra, ras), (al.dec, decs)]:
                val_list.append(val)
        except:
            for val_list in [tss, nss, gammas, pvals, ras, decs]:
                val_list.append(np.nan)
    
for name, val_list in [('TS', tss), ('ns', nss), ('p', pvals),
                      ('ra', ras), ('dec', decs), ('gamma', gammas)]:
    alert_df[f'{name}_steady'] = np.asarray(val_list)


alert_df.to_csv(base_output + '/alert_dataframe_with_results.csv')