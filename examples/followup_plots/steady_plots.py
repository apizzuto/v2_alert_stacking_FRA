from francis.time_integrated_scripts import steady_sensitivity_fits as ssf
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('/home/apizzuto/Nova/scripts/novae_plots.mplstyle')
import pandas as pd

alerts = pd.read_csv(
    '/data/user/apizzuto/fast_response_skylab/alert_event_followup/FRANCIS/' \
    + 'francis/icecube_misc/alert_dataframe.csv'
    )

output_dir = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/alert_skymaps/'
for ind in range(68, 276):
    if ind % 25 == 0:
        print(ind, end=' ')
    run_id = str(alerts['Run ID'][ind])
    event_id = str(alerts['Event ID'][ind])
    
    if run_id == '123986':
        print(str(ind) + ' is a bad ind')
        continue
    elif run_id == '129434':
        print(str(ind) + ' is a bad ind')
        continue
    elif run_id == '133781':
        print(str(ind) + ' is a bad ind')
        continue
    base_str = f'{output_dir}Run_{run_id}_Event_{event_id}'
    
    ssf.plot_skymap(ind, LLH=True)
    plt.savefig(base_str + '_allsky.png', dpi=120, bbox_inches='tight')
    plt.close()
    
    ssf.plot_zoom(ind, LLH=True)
    plt.savefig(base_str + '_zoom_LLH.png', dpi=120, bbox_inches='tight')
    plt.close()
    
    ssf.plot_zoom(ind, LLH=False)
    plt.savefig(base_str + '_probs.png', dpi=120, bbox_inches='tight')
    plt.close()

output_dir = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/' \
    + 'time_integrated_summary_plots/'

# Make all of the summary plots. This takes a looooooong time (O(hours))
for ind in range(len(alerts)):
    if ind % 50 == 0:
        print(ind, end=' ')
    try:
        event_header = ssf.panel_plot_with_text(ind, smear=True, return_info=True)
        plt.savefig(output_dir \
            + f'steady_panel_plot_with_text_{ind}_run' \
            + f'_{event_header["RUNID"]}_event_{event_header["EVENTID"]}.png',
            bbox_inches='tight', dpi=120)
        plt.close()
    except Exception as e:
        print(ind, e)

# Make the sensitivity vs. declination plot. This takes ~ 10 minutes
