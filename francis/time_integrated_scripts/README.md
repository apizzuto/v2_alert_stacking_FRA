# Time-integrated followup scripts

This directory contains the scripts to run trials for the time-integrated followup for each alert event included in the v2 alert catalog. 

Scripts:
`steady_background_trials.py`: Runs background trials. Args:
* `--ntrials`: Number of trials (default 1000)
* `--i`: Alert event index, identifier for alert events ordered by date (from 0 to 276)
* `--rng`: Random number seed (default 1)
* `--verbose`: Flag to raise if you want assorted printed output
* `--smear`: Raise this flag if you want to account for systematics when converting the millipede LLH map to a spatial prior

`signal_trials.py`: Runs background trials. Args:
* `--ntrials`: Number of trials (default 1000)
* `--g`: Injected spectral index 
* `--i`: Alert event index, identifier for alert events ordered by date (from 0 to 276)
* `--rng`: Random number seed (default 1)
* `--verbose`: Flag to raise if you want assorted printed output
* `--smear`: Raise this flag if you want to account for systematics when converting the millipede LLH map to a spatial prior
* `--fit`: If this flag is raised, then don't poisson fluctuate fluxes. I use this flag when performing round-trip tests (ns-bias tests) and I do not raise this flag when I want to calculate sensitivity

`steady_alert_unblind.py`: Runs background trials. Args:
* `--i`: Alert event index, identifier for alert events ordered by date (from 0 to 276)
* `--rng`: Random number seed (default 1)
* `--verbose`: Flag to raise if you want assorted printed output
* `--smear`: Raise this flag if you want to account for systematics when converting the millipede LLH map to a spatial prior
* `--local_skymap`: If you raise this flag, then a local scan of TS will be saved in an area near the alert. If you don't, then it will just output the information on the hotspot. Note: if you don't raise this flag, we actually run 2 trials with different seeds to make sure the fits are somewhat stable

The other three `.py` files in this repository are not meant to be run, they are there to provide helper functions to the scripts in this directory as well as to some of the plotting scripts in other directories. 



